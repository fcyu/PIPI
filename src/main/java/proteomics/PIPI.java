package proteomics;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Output.WritePepXml;
import proteomics.PTM.InferPTM;
import ProteomicsLibrary.Binomial;
import ProteomicsLibrary.PrepareSpectrum;
import proteomics.Types.*;
import proteomics.Index.BuildIndex;
import proteomics.Parameter.Parameter;
import proteomics.Spectrum.PreSpectra;
import ProteomicsLibrary.MassTool;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;

import java.io.*;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.sql.*;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.locks.ReentrantLock;

public class PIPI {

    private static final Logger logger = LoggerFactory.getLogger(PIPI.class);
    public static final String versionStr = "1.4.6";
    static final boolean useXcorr = true;

    public static final int[] debugScanNumArray = new int[]{};

    public static void main(String[] args) {
        long startTime = System.nanoTime();

        // Process inputs
        if (args.length != 2) {
            help();
        }

        // Set parameters
        String parameterPath = args[0].trim();
        String spectraPath = args[1].trim();

        logger.info("Running PIPI version {}.", versionStr);

        String dbName = null;
        String hostName = "unknown-host";
        try {
            hostName = InetAddress.getLocalHost().getHostName();
            logger.info("Computer: {}.", hostName);
        } catch (UnknownHostException ex) {
            logger.warn("Cannot get the computer's name.");
        }

        try {
            logger.info("Spectra: {}, parameter: {}.", spectraPath, parameterPath);

            dbName = String.format(Locale.US, "PIPI.%s.%s.temp.db", hostName, new SimpleDateFormat("yyyy-MM-dd-HH-mm-ss").format(Calendar.getInstance().getTime()));
            new PIPI(parameterPath, spectraPath, dbName);
        } catch (Exception ex) {
            ex.printStackTrace();
            logger.error(ex.toString());
        } finally {
            if (dbName != null) {
                (new File(dbName)).delete();
                (new File(dbName + "-wal")).delete();
                (new File(dbName + "-shm")).delete();
            }
        }

        double totalHour = (double) (System.nanoTime() - startTime) * 1e-9 / 3600;
        logger.info("Running time: {} hours.", totalHour);
        logger.info("Done!");
    }

    private PIPI(String parameterPath, String spectraPath, String dbName) throws Exception {
        // Get the parameter map
        Parameter parameter = new Parameter(parameterPath);
        Map<String, String> parameterMap = parameter.returnParameterMap();
        double ms2Tolerance = Double.valueOf(parameterMap.get("ms2_tolerance"));
        double ms1Tolerance = Double.valueOf(parameterMap.get("ms1_tolerance"));
        double leftInverseMs1Tolerance = 1 / (1 + ms1Tolerance * 1e-6);
        double rightInverseMs1Tolerance = 1 / (1 - ms1Tolerance * 1e-6);
        int ms1ToleranceUnit = Integer.valueOf(parameterMap.get("ms1_tolerance_unit"));
        double minClear = Double.valueOf(parameterMap.get("min_clear_mz"));
        double maxClear = Double.valueOf(parameterMap.get("max_clear_mz"));
        String percolatorPath = parameterMap.get("percolator_path");
        boolean outputPercolatorInput = (Integer.valueOf(parameterMap.get("output_percolator_input")) == 1);

        // print all the parameters
        logger.info("Parameters:");
        for (String k : parameterMap.keySet()) {
            logger.info("{} = {}", k, parameterMap.get(k));
        }

        // Check if Percolator can be executed.
        if (!(new File(percolatorPath)).exists()) {
            throw new NullPointerException(String.format(Locale.US, "Cannot find Percolator from %s.", percolatorPath));
        }

        if (!(new File(percolatorPath)).canExecute()) {
            throw new Exception(String.format(Locale.US, "Percolator (%s) exits but cannot be executed.", percolatorPath));
        }

        String[] tempArray = parameterMap.get("ms_level").split(",");
        Set<Integer> msLevelSet = new HashSet<>(tempArray.length + 1, 1);
        for (String temp : tempArray) {
            msLevelSet.add(Integer.valueOf(temp));
        }

        String labelling = "N14";
        if (parameterMap.get("15N").trim().contentEquals("1")) {
            logger.info("N15 mode is on...");
            labelling = "N15";
        }

        if (parameterMap.get("add_decoy").contentEquals("0")) {
            logger.warn("add_decoy = 0. Won't search the decoy sequences and estimate FDR.");
        }

        if (parameterMap.get("add_contaminant").contentEquals("0")) {
            logger.warn("add_contaminant = 0. Won't search the build-in contaminant proteins.");
        }

        logger.info("Indexing protein database...");
        BuildIndex buildIndex = new BuildIndex(parameterMap, labelling, true, parameterMap.get("add_decoy").contentEquals("1"), parameterMap.get("add_contaminant").contentEquals("1"));
        MassTool massTool = buildIndex.returnMassTool();
        InferPTM inferPTM = buildIndex.getInferPTM();

        logger.info("Reading spectra...");
        File spectraFile = new File(spectraPath);
        if ((!spectraFile.exists() || (spectraFile.isDirectory()))) {
            throw new FileNotFoundException("The spectra file not found.");
        }
        String[] temp = spectraPath.split("\\.");
        String ext = temp[temp.length - 1];
        JMzReader spectraParser;
        if (ext.contentEquals("mzXML")) {
            spectraParser = new MzXMLFile(spectraFile);
        } else if (ext.toLowerCase().contentEquals("mgf")) {
            spectraParser = new MgfFile(spectraFile);
        } else {
            throw new Exception(String.format(Locale.US, "Unsupported file format %s. Currently, PIPI only support mzXML and MGF.", ext));
        }

        String sqlPath = "jdbc:sqlite:" + dbName;
        Class.forName("org.sqlite.JDBC").newInstance();

        PreSpectra preSpectra = new PreSpectra(spectraParser, ms1Tolerance, ms1ToleranceUnit, massTool, ext, msLevelSet, sqlPath);

        logger.info("Start searching...");
        int threadNum = Integer.valueOf(parameterMap.get("thread_num"));
        if (threadNum == 0) {
            threadNum = 3 + Runtime.getRuntime().availableProcessors();
        }
        if (debugScanNumArray.length > 0) {
            threadNum = 1;
        }
        ExecutorService threadPool = Executors.newFixedThreadPool(threadNum);
        PrepareSpectrum preSpectrum = new PrepareSpectrum(massTool);
        ArrayList<Future<PIPIWrap.Entry>> taskList = new ArrayList<>(preSpectra.getUsefulSpectraNum() + 10);
        Connection sqlConnection = DriverManager.getConnection(sqlPath);
        Statement sqlStatement = sqlConnection.createStatement();
        ResultSet sqlResultSet = sqlStatement.executeQuery("SELECT scanId, precursorCharge, precursorMass FROM spectraTable");
        ReentrantLock lock = new ReentrantLock();
        Binomial binomial = new Binomial(Integer.valueOf(parameterMap.get("max_peptide_length")) * 2);
        while (sqlResultSet.next()) {
            String scanId = sqlResultSet.getString("scanId");
            int precursorCharge = sqlResultSet.getInt("precursorCharge");
            double precursorMass = sqlResultSet.getDouble("precursorMass");
            taskList.add(threadPool.submit(new PIPIWrap(buildIndex, massTool, ms1Tolerance, leftInverseMs1Tolerance, rightInverseMs1Tolerance, ms1ToleranceUnit, ms2Tolerance, inferPTM.getMinPtmMass(), inferPTM.getMaxPtmMass(), Math.min(precursorCharge > 1 ? precursorCharge - 1 : 1, 3), spectraParser, minClear, maxClear, lock, scanId, precursorCharge, precursorMass, inferPTM, preSpectrum, sqlPath, binomial)));
        }
        sqlResultSet.close();
        sqlStatement.close();

        // check progress every minute, record results,and delete finished tasks.
        PreparedStatement sqlPreparedStatement = sqlConnection.prepareStatement("REPLACE INTO spectraTable (scanNum, scanId, precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, labelling, peptide, theoMass, isDecoy, globalRank, normalizedCorrelationCoefficient, score, deltaLCn, deltaCn, matchedPeakNum, ionFrac, matchedHighestIntensityFrac, explainedAaFrac, otherPtmPatterns, aScore) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
        sqlConnection.setAutoCommit(false);
        int lastProgress = 0;
        int resultCount = 0;
        int totalCount = taskList.size();
        int count = 0;
        while (count < totalCount) {
            // record search results and delete finished ones.
            List<Future<PIPIWrap.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCount - count);
            for (Future<PIPIWrap.Entry> task : taskList) {
                if (task.isDone()) {
                    if (task.get() != null) {
                        PIPIWrap.Entry entry = task.get();
                        sqlPreparedStatement.setInt(1, entry.scanNum);
                        sqlPreparedStatement.setString(2, entry.scanId);
                        sqlPreparedStatement.setInt(3, entry.precursorCharge);
                        sqlPreparedStatement.setDouble(4, entry.precursorMass);
                        sqlPreparedStatement.setString(5, entry.mgfTitle);
                        sqlPreparedStatement.setInt(6, entry.isotopeCorrectionNum);
                        sqlPreparedStatement.setDouble(7, entry.ms1PearsonCorrelationCoefficient);
                        sqlPreparedStatement.setString(8, entry.labelling);
                        sqlPreparedStatement.setString(9, entry.peptide);
                        sqlPreparedStatement.setDouble(10, entry.theoMass);
                        sqlPreparedStatement.setInt(11, entry.isDecoy);
                        sqlPreparedStatement.setInt(12, entry.globalRank);
                        sqlPreparedStatement.setDouble(13, entry.normalizedCorrelationCoefficient);
                        sqlPreparedStatement.setDouble(14, entry.score);
                        sqlPreparedStatement.setDouble(15, entry.deltaLCn);
                        sqlPreparedStatement.setDouble(16, entry.deltaCn);
                        sqlPreparedStatement.setInt(17, entry.matchedPeakNum);
                        sqlPreparedStatement.setDouble(18, entry.ionFrac);
                        sqlPreparedStatement.setDouble(19, entry.matchedHighestIntensityFrac);
                        sqlPreparedStatement.setDouble(20, entry.explainedAaFrac);
                        sqlPreparedStatement.setString(21, entry.otherPtmPatterns);
                        sqlPreparedStatement.setString(22, entry.aScore);
                        sqlPreparedStatement.executeUpdate();
                        ++resultCount;
                    }
                    toBeDeleteTaskList.add(task);
                } else if (task.isCancelled()) {
                    toBeDeleteTaskList.add(task);
                }
            }
            count += toBeDeleteTaskList.size();
            taskList.removeAll(toBeDeleteTaskList);
            taskList.trimToSize();

            sqlConnection.commit();

            int progress = count * 20 / totalCount;
            if (progress != lastProgress) {
                logger.info("Searching {}%...", progress * 5);
                lastProgress = progress;
            }

            if (count == totalCount) {
                break;
            }
            Thread.sleep(6000);
        }

        // shutdown threads.
        threadPool.shutdown();
        if (!threadPool.awaitTermination(60, TimeUnit.SECONDS)) {
            threadPool.shutdownNow();
            if (!threadPool.awaitTermination(60, TimeUnit.SECONDS))
                throw new Exception("Pool did not terminate");
        }

        sqlConnection.commit();
        sqlConnection.setAutoCommit(true);
        sqlConnection.close();
        if (lock.isLocked()) {
            lock.unlock();
        }

        if (resultCount == 0) {
            throw new Exception("There is no useful results.");
        }

        String percolatorInputFileName = spectraPath + "." + labelling + ".input.temp";
        writePercolator(percolatorInputFileName, buildIndex.getPeptide0Map(), sqlPath);
        Map<Integer, PercolatorEntry> percolatorResultMap = null;

        if (parameterMap.get("add_decoy").contentEquals("0")) {
            logger.warn("add_decoy = 0. Don't estimate FDR.");
        } else {
            logger.info("Estimating FDR...");
            String percolatorOutputFileName = spectraPath + "." + labelling + ".output.temp";
            String percolatorProteinOutputFileName = spectraPath + "." + labelling + ".protein.tsv";
            percolatorResultMap = runPercolator(percolatorPath, percolatorInputFileName, percolatorOutputFileName, percolatorProteinOutputFileName, parameterMap.get("db") + ".TD.fasta", parameterMap.get("enzyme_name_1"));
            if (percolatorResultMap.isEmpty()) {
                throw new Exception(String.format(Locale.US, "Percolator failed to estimate FDR. Please check if Percolator is installed and the percolator_path in %s is correct.", parameterPath));
            }
            if (!outputPercolatorInput) {
                (new File(percolatorOutputFileName)).delete();
            }
        }

        if (!outputPercolatorInput) {
            (new File(percolatorInputFileName)).delete();
        }

        logger.info("Saving results...");
        writeFinalResult(percolatorResultMap, spectraPath + "." + labelling + ".pipi.csv", buildIndex.getPeptide0Map(), sqlPath);
        new WritePepXml(spectraPath + "." + labelling + ".pipi.pep.xml", spectraPath, parameterMap, massTool.getMassTable(), percolatorResultMap, buildIndex.getPeptide0Map(), buildIndex.returnFixModMap(), sqlPath);
    }

    private static void help() {
        String helpStr = "PIPI version " + versionStr + "\r\n"
                + "A tool identifying peptides with unlimited PTM.\r\n"
                + "Author: Fengchao Yu\r\n"
                + "Email: fyuab@connect.ust.hk\r\n"
                + "PIPI usage: java -Xmx25g -jar /path/to/PIPI.jar <parameter_file> <data_file>\r\n"
                + "\t<parameter_file>: parameter file. Can be download along with PIPI.\r\n"
                + "\t<data_file>: spectra data file (mzXML)\r\n"
                + "\texample: java -Xmx32g -jar PIPI.jar parameter.def data.mzxml\r\n";
        System.out.print(helpStr);
        System.exit(1);
    }

    private void writePercolator(String resultPath, Map<String, Peptide0> peptide0Map, String sqlPath) throws IOException, SQLException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(resultPath));
        writer.write("id\tlabel\tscannr\tscore\tdelta_c_n\tdelta_L_c_n\tnormalized_cross_corr\tglobal_search_rank\tabs_ppm\tion_frac\tmatched_high_peak_frac\tcharge1\tcharge2\tcharge3\tcharge4\tcharge5\tcharge6\texplained_aa_frac\tpeptide\tprotein\n");
        Connection sqlConnection = DriverManager.getConnection(sqlPath);
        Statement sqlStatement = sqlConnection.createStatement();
        ResultSet sqlResultSet = sqlStatement.executeQuery("SELECT scanNum, precursorCharge, precursorMass, peptide, theoMass, isDecoy, globalRank, normalizedCorrelationCoefficient, score, deltaLCn, deltaCn, matchedPeakNum, ionFrac, matchedHighestIntensityFrac, explainedAaFrac FROM spectraTable");
        while (sqlResultSet.next()) {
            String peptide = sqlResultSet.getString("peptide");
            if (!sqlResultSet.wasNull()) {
                int charge = sqlResultSet.getInt("precursorCharge");
                double theoMass = sqlResultSet.getDouble("theoMass");
                double expMass = sqlResultSet.getDouble("precursorMass");
                double massDiff = getMassDiff(expMass, theoMass, MassTool.C13_DIFF);

                Peptide0 peptide0 = peptide0Map.get(peptide.replaceAll("[^ncA-Z]+", ""));
                TreeSet<String> proteinIdSet = new TreeSet<>();
                for (String protein : peptide0.proteins) {
                    proteinIdSet.add(protein.trim());
                }

                StringBuilder sb = new StringBuilder(20);
                for (int i = 0; i < 6; ++i) {
                    if (i == charge - 1) {
                        sb.append(1);
                    } else {
                        sb.append(0);
                    }
                    sb.append("\t");
                }

                int scanNum = sqlResultSet.getInt("scanNum");
                int isDecoy = sqlResultSet.getInt("isDecoy");
                int globalRank = sqlResultSet.getInt("globalRank");
                double normalizedCorrelationCoefficient = sqlResultSet.getDouble("normalizedCorrelationCoefficient");
                double score = sqlResultSet.getDouble("score");
                double deltaLCn = sqlResultSet.getDouble("deltaLCn");
                double deltaCn = sqlResultSet.getDouble("deltaCn");
                double ionFrac = sqlResultSet.getDouble("ionFrac");
                double matchedHighestIntensityFrac = sqlResultSet.getDouble("matchedHighestIntensityFrac");
                double explainedAaFrac = sqlResultSet.getDouble("explainedAaFrac");

                if (isDecoy == 1) {
                    writer.write(scanNum + "\t-1\t" + scanNum + "\t" + score + "\t" + deltaCn + "\t" + deltaLCn + "\t" + normalizedCorrelationCoefficient + "\t" + globalRank + "\t" + Math.abs(massDiff * 1e6 / theoMass) + "\t" + ionFrac + "\t" + matchedHighestIntensityFrac + "\t" + sb.toString() + explainedAaFrac + "\t" + peptide0.leftFlank + "." + peptide + "." + peptide0.rightFlank + "\t" + String.join("\t", proteinIdSet) + "\n");
                } else {
                    writer.write(scanNum + "\t1\t" + scanNum + "\t" + score + "\t" + deltaCn + "\t" + deltaLCn + "\t" + normalizedCorrelationCoefficient + "\t" + globalRank + "\t" + Math.abs(massDiff * 1e6 / theoMass) + "\t" + ionFrac + "\t" + matchedHighestIntensityFrac + "\t" + sb.toString() + explainedAaFrac + "\t" + peptide0.leftFlank + "." + peptide + "." + peptide0.rightFlank + "\t" + String.join("\t", proteinIdSet) + "\n");
                }
            }
        }
        writer.close();
        sqlResultSet.close();
        sqlStatement.close();
        sqlConnection.close();
    }

    private static Map<Integer, PercolatorEntry> runPercolator(String percolatorPath, String percolatorInputFileName, String percolatorOutputFileName, String percolatorProteinOutputFileName, String tdFastaPath, String enzymeName) throws Exception {
        Map<Integer, PercolatorEntry> percolatorResultMap = new HashMap<>();
        if ((new File(percolatorInputFileName)).exists()) {
            String percolatorEnzyme;
            switch (enzymeName) {
                case "Trypsin": percolatorEnzyme = "trypsin";
                break;
                case "Trypsin/P": percolatorEnzyme = "trypsinp";
                break;
                case "TrypsinR": percolatorEnzyme = "trypsin";
                break;
                case "LysC": percolatorEnzyme = "lys-c";
                break;
                case "ArgC": percolatorEnzyme = "arg-c";
                break;
                case "Chymotrypsin": percolatorEnzyme = "chymotrypsin";
                break;
                case "GluC": percolatorEnzyme = "glu-c";
                break;
                case "LysN": percolatorEnzyme = "lys-n";
                break;
                case "AspN": percolatorEnzyme = "asp-n";
                break;
                default: percolatorEnzyme = "trypsin";
                break;
            }

            String[] commands = new String[]{percolatorPath, "--only-psms", "--verbose", "0", "--no-terminate", "--protein-decoy-pattern", "DECOY_", "--picked-protein", tdFastaPath, "--protein-enzyme", percolatorEnzyme, "--protein-report-fragments", "--protein-report-duplicates", "--results-proteins", percolatorProteinOutputFileName, "--results-psms", percolatorOutputFileName, percolatorInputFileName};

            Process ps = Runtime.getRuntime().exec(commands);
            ps.waitFor();

            if (!(new File(percolatorOutputFileName).exists()) || ps.exitValue() != 0) {
                BufferedReader reader = new BufferedReader(new InputStreamReader(ps.getInputStream()));
                String line;
                while ((line = reader.readLine()) != null) {
                    logger.info("[Percolator info]: {}", line.trim());
                }
                reader.close();
                reader = new BufferedReader(new InputStreamReader(ps.getErrorStream()));
                while ((line = reader.readLine()) != null) {
                    logger.info("[Percolator error]: {}", line.trim());
                }
                reader.close();
                throw new NullPointerException("Percolator didn't exit normally.");
            } else {
                BufferedReader reader = new BufferedReader(new InputStreamReader(ps.getInputStream()));
                String line;
                while ((line = reader.readLine()) != null) {
                    logger.info("[Percolator info]: {}", line.trim());
                }
                reader.close();
                reader = new BufferedReader(new InputStreamReader(ps.getErrorStream()));
                while ((line = reader.readLine()) != null) {
                    logger.info("[Percolator info]: {}", line.trim());
                }
                reader.close();
            }

            logger.info("Percolator finished normally");

            BufferedReader reader = new BufferedReader(new FileReader(percolatorOutputFileName));
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (!line.startsWith("PSMId")) {
                    String[] parts = line.split("\t");
                    percolatorResultMap.put(Integer.valueOf(parts[0]), new PercolatorEntry(Double.valueOf(parts[1]), parts[2], parts[3]));
                }
            }
            reader.close();
        } else {
            logger.error("Cannot find Percolator input file (from {}) for estimating Percolator Q-Value.", percolatorInputFileName);
            return percolatorResultMap;
        }

        return percolatorResultMap;
    }

    private void writeFinalResult(Map<Integer, PercolatorEntry> percolatorResultMap, String outputPath, Map<String, Peptide0> peptide0Map, String sqlPath) throws IOException, SQLException {
        TreeMap<Double, List<String>> tempMap = new TreeMap<>();

        BufferedWriter writer = new BufferedWriter(new FileWriter(outputPath));
        if (percolatorResultMap == null) {
            writer.write("scan_num,peptide,charge,theo_mass,exp_mass,abs_ppm,A_score,protein_ID,score,delta_C_n,other_PTM_patterns,MGF_title,labelling,isotope_correction,MS1_pearson_correlation_coefficient\n");
        } else {
            writer.write("scan_num,peptide,charge,theo_mass,exp_mass,abs_ppm,A_score,protein_ID,score,delta_C_n,percolator_score,posterior_error_prob,q_value,other_PTM_patterns,MGF_title,labelling,isotope_correction,MS1_pearson_correlation_coefficient\n");
        }

        Connection sqlConnection = DriverManager.getConnection(sqlPath);
        Statement sqlStatement = sqlConnection.createStatement();
        ResultSet sqlResultSet = sqlStatement.executeQuery("SELECT scanNum, precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, labelling, peptide, theoMass, isDecoy, score, otherPtmPatterns, aScore, deltaCn FROM spectraTable");
        while (sqlResultSet.next()) {
            int isDecoy = sqlResultSet.getInt("isDecoy");
            if (!sqlResultSet.wasNull()) {
                if (isDecoy == 0) {
                    int scanNum = sqlResultSet.getInt("scanNum");
                    double expMass = sqlResultSet.getDouble("precursorMass");
                    String peptide = sqlResultSet.getString("peptide");
                    double theoMass = sqlResultSet.getDouble("theoMass");
                    double massDiff = getMassDiff(expMass, theoMass, MassTool.C13_DIFF);
                    double ppm = Math.abs(massDiff * 1e6 / theoMass);

                    Peptide0 peptide0 = peptide0Map.get(peptide.replaceAll("[^ncA-Z]+", ""));
                    TreeSet<String> proteinIdSet = new TreeSet<>();
                    for (String protein : peptide0.proteins) {
                        proteinIdSet.add(protein.trim());
                    }

                    String aScore = sqlResultSet.getString("aScore");

                    if (percolatorResultMap == null) {
                        double score = sqlResultSet.getDouble("score");
                        String str = String.format(Locale.US, "%d,%s,%d,%f,%f,%f,%s,%s,%f,%f,%s,\"%s\",%s,%d,%f\n", scanNum, peptide, sqlResultSet.getInt("precursorCharge"), theoMass, expMass, ppm, aScore, String.join(";", proteinIdSet).replaceAll(",", "~"), score, sqlResultSet.getDouble("deltaCn"), sqlResultSet.getString("otherPtmPatterns"), sqlResultSet.getString("mgfTitle"), sqlResultSet.getString("labelling"), sqlResultSet.getInt("isotopeCorrectionNum"), sqlResultSet.getDouble("ms1PearsonCorrelationCoefficient"));
                        if (tempMap.containsKey(score)) {
                            tempMap.get(score).add(str);
                        } else {
                            List<String> tempList = new LinkedList<>();
                            tempList.add(str);
                            tempMap.put(score, tempList);
                        }
                    } else {
                        PercolatorEntry percolatorEntry = percolatorResultMap.get(scanNum);
                        String str = String.format(Locale.US, "%d,%s,%d,%f,%f,%f,%s,%s,%f,%f,%f,%s,%s,%s,\"%s\",%s,%d,%f\n", scanNum, peptide, sqlResultSet.getInt("precursorCharge"), theoMass, expMass, ppm, aScore, String.join(";", proteinIdSet).replaceAll(",", "~"), sqlResultSet.getDouble("score"), sqlResultSet.getDouble("deltaCn"), percolatorEntry.percolatorScore, percolatorEntry.PEP, percolatorEntry.qValue, sqlResultSet.getString("otherPtmPatterns"), sqlResultSet.getString("mgfTitle"), sqlResultSet.getString("labelling"), sqlResultSet.getInt("isotopeCorrectionNum"), sqlResultSet.getDouble("ms1PearsonCorrelationCoefficient"));
                        if (tempMap.containsKey(percolatorResultMap.get(scanNum).percolatorScore)) {
                            tempMap.get(percolatorResultMap.get(scanNum).percolatorScore).add(str);
                        } else {
                            List<String> tempList = new LinkedList<>();
                            tempList.add(str);
                            tempMap.put(percolatorResultMap.get(scanNum).percolatorScore, tempList);
                        }
                    }
                }
            }
        }

        sqlResultSet.close();
        sqlStatement.close();
        sqlConnection.close();

        Double[] tempArray = tempMap.keySet().toArray(new Double[0]);
        for (int i = tempArray.length - 1; i >= 0; --i) {
            List<String> tempList = tempMap.get(tempArray[i]);
            for (String tempStr : tempList) {
                writer.write(tempStr);
            }
        }
        writer.close();
    }

    public static double getMassDiff(double expMass, double theoMass, double C13Diff) {
        double massDiff1 = expMass - theoMass;
        double massDiff2 = expMass - theoMass - C13Diff;
        double massDiff3 = expMass - theoMass - 2 * C13Diff;
        double absMassDiff1 = Math.abs(massDiff1);
        double absMassDiff2 = Math.abs(massDiff2);
        double absMassDiff3 = Math.abs(massDiff3);

        if ((absMassDiff1 <= absMassDiff2) && (absMassDiff1 <= absMassDiff2)) {
            return massDiff1;
        } else if ((absMassDiff2 <= absMassDiff1) && (absMassDiff2 <= absMassDiff3)) {
            return massDiff2;
        } else {
            return massDiff3;
        }
    }
}
