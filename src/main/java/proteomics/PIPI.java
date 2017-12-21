package proteomics;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Output.WritePepXml;
import proteomics.PTM.InferPTM;
import proteomics.Spectrum.PreSpectrum;
import proteomics.Types.*;
import proteomics.Index.BuildIndex;
import proteomics.Parameter.Parameter;
import proteomics.Spectrum.PreSpectra;
import proteomics.TheoSeq.MassTool;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;

import java.io.*;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.sql.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.locks.ReentrantLock;

public class PIPI {

    private static final Logger logger = LoggerFactory.getLogger(PIPI.class);
    public static final String versionStr = "1.4.3";
    public static final boolean useXcorr = true;

    public static final boolean DEV = false;
    public static final int[] debugScanNumArray = new int[]{};

    public static void main(String[] args) {
        // Process inputs
        if (args.length != 2) {
            help();
        }

        // Set parameters
        String parameterPath = args[0].trim();
        String spectraPath = args[1].trim();

        logger.info("Running PIPI version {}.", versionStr);

        try {
            String hostName = InetAddress.getLocalHost().getHostName();
            logger.info("Computer: {}.", hostName);
        } catch (UnknownHostException ex) {
            logger.warn("Cannot get the computer's name.");
        }

        logger.info("Spectra: {}, parameter: {}.", spectraPath, parameterPath);

        if (DEV) {
            logger.info("In DEV mode.");
        }

        new PIPI(parameterPath, spectraPath);
    }

    private PIPI(String parameterPath, String spectraPath) {
        // Get the parameter map
        Parameter parameter = new Parameter(parameterPath);
        Map<String, String> parameterMap = parameter.returnParameterMap();
        float ms2Tolerance = Float.valueOf(parameterMap.get("ms2_tolerance"));
        float ms1Tolerance = Float.valueOf(parameterMap.get("ms1_tolerance"));
        int ms1ToleranceUnit = Integer.valueOf(parameterMap.get("ms1_tolerance_unit"));
        int maxMs2Charge = Integer.valueOf(parameterMap.get("max_ms2_charge"));
        float minClear = Float.valueOf(parameterMap.get("min_clear_mz"));
        float maxClear = Float.valueOf(parameterMap.get("max_clear_mz"));
        String percolatorPath = parameterMap.get("percolator_path");
        boolean outputPercolatorInput = (DEV || (Integer.valueOf(parameterMap.get("output_percolator_input")) == 1));

        String[] tempArray = parameterMap.get("ms_level").split(",");
        Set<Integer> msLevelSet = new HashSet<>(tempArray.length + 1, 1);
        for (String temp : tempArray) {
            msLevelSet.add(Integer.valueOf(temp));
        }

        String labeling = "N14";
        if (parameterMap.get("15N").trim().contentEquals("1")) {
            logger.info("N15 mode is on...");
            labeling = "N15";
        }

        logger.info("Indexing protein database...");
        BuildIndex buildIndexObj = new BuildIndex(parameterMap, labeling);
        MassTool massToolObj = buildIndexObj.returnMassToolObj();

        float minPtmMass = Float.valueOf(parameterMap.get("min_ptm_mass"));
        float maxPtmMass = Float.valueOf(parameterMap.get("max_ptm_mass"));

        logger.info("Reading spectra...");
        JMzReader spectraParser = null;
        String ext = "";
        String sqlPath = "jdbc:sqlite:PIPI.temp.db";
        try {
            File spectraFile = new File(spectraPath);
            if ((!spectraFile.exists() || (spectraFile.isDirectory()))) {
                throw new FileNotFoundException("The spectra file not found.");
            }
            String[] temp = spectraPath.split("\\.");
            ext = temp[temp.length - 1];
            if (ext.contentEquals("mzXML")) {
                spectraParser = new MzXMLFile(spectraFile);
            } else if (ext.toLowerCase().contentEquals("mgf")) {
                spectraParser = new MgfFile(spectraFile);
            } else {
                throw new Exception(String.format(Locale.US, "Unsupported file format %s. Currently, PIPI only support mzXML and MGF.", ext));
            }
            Class.forName("org.sqlite.JDBC").newInstance();
        } catch (Exception ex) {
            ex.printStackTrace();
            logger.error(ex.toString());
            System.exit(1);
        }

        PreSpectra preSpectraObj = new PreSpectra(spectraParser, parameterMap, massToolObj, ext, msLevelSet, sqlPath);

        if (DEV) {
            try (BufferedWriter writer = new BufferedWriter(new FileWriter("spectrum.dev.csv"))) {
                writer.write("scanNum,charge,finalIsotopeCorrectionNum,isotopeCorrectionNum,pearsonCorrelationCoefficient,expMz1,expMz2,expMz3,expInt1,expInt2,expInt3,theoMz1,theoMz2,theoMz3,theoInt1,theoInt2,theoInt3\n");
                Map<Integer, TreeMap<Integer, TreeSet<PreSpectra.DevEntry>>> scanDevEntryMap = preSpectraObj.getScanDevEntryMap();
                for (int scanNum : scanDevEntryMap.keySet()) {
                    TreeMap<Integer, TreeSet<PreSpectra.DevEntry>> chargeDevEntryMap = scanDevEntryMap.get(scanNum);
                    for (int charge : chargeDevEntryMap.keySet()) {
                        for (PreSpectra.DevEntry devEntry : chargeDevEntryMap.get(charge)) {
                            writer.write(String.format(Locale.US, "%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", scanNum, charge, devEntry.isotopeCorrectionNum, devEntry.isotopeCorrectionNum, devEntry.pearsonCorrelationCoefficient, devEntry.expMatrix[0][0], devEntry.expMatrix[1][0], devEntry.expMatrix[2][0], devEntry.expMatrix[0][1], devEntry.expMatrix[1][1], devEntry.expMatrix[2][1], devEntry.theoMatrix[0][0], devEntry.theoMatrix[1][0], devEntry.theoMatrix[2][0], devEntry.theoMatrix[0][1], devEntry.theoMatrix[1][1], devEntry.theoMatrix[2][1]));
                        }
                    }
                }
            } catch (IOException ex) {
                logger.error(ex.toString());
                ex.printStackTrace();
                System.exit(1);
            }
        }

        logger.info("Start searching...");
        int threadNum = Integer.valueOf(parameterMap.get("thread_num"));
        if (threadNum == 0) {
            threadNum = 3 + Runtime.getRuntime().availableProcessors();
        }
        ExecutorService threadPool = Executors.newFixedThreadPool(threadNum);

        InferPTM inferPTM = new InferPTM(massToolObj, maxMs2Charge, buildIndexObj.returnFixModMap(), buildIndexObj.getInference3SegmentObj().getVarModParamSet(), minPtmMass, maxPtmMass, ms2Tolerance);
        PreSpectrum preSpectrumObj = new PreSpectrum(massToolObj);
        List<Future<Boolean>> taskList = new LinkedList<>();

        try {
            Connection sqlConnection = DriverManager.getConnection(sqlPath);
            Statement sqlStatement = sqlConnection.createStatement();
            ResultSet sqlResultSet = sqlStatement.executeQuery("SELECT scanId, precursorCharge, precursorMass FROM spectraTable");
            ReentrantLock lock = new ReentrantLock();
            while (sqlResultSet.next()) {
                String scanId = sqlResultSet.getString("scanId");
                int precursorCharge = sqlResultSet.getInt("precursorCharge");
                float precursorMass = sqlResultSet.getFloat("precursorMass");
                taskList.add(threadPool.submit(new PIPIWrap(buildIndexObj, massToolObj, ms1Tolerance, ms1ToleranceUnit, ms2Tolerance, minPtmMass, maxPtmMass, Math.min(precursorCharge > 1 ? precursorCharge - 1 : 1, maxMs2Charge), spectraParser, minClear, maxClear, lock, scanId, precursorCharge, precursorMass, inferPTM, preSpectrumObj, sqlConnection)));
            }
            sqlResultSet.close();
            sqlStatement.close();

            // check progress every minute, record results,and delete finished tasks.
            int lastProgress = 0;
            int resultCount = 0;
            int totalCount = taskList.size();
            int count = 0;
            while (count < totalCount) {
                // record search results and delete finished ones.
                List<Future<Boolean>> toBeDeleteTaskList = new LinkedList<>();
                for (Future<Boolean> task : taskList) {
                    if (task.isDone()) {
                        if (task.get()) {
                            ++resultCount;
                        }
                        toBeDeleteTaskList.add(task);
                    } else if (task.isCancelled()) {
                        toBeDeleteTaskList.add(task);
                    }
                }
                count += toBeDeleteTaskList.size();
                taskList.removeAll(toBeDeleteTaskList);

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
                    System.err.println("Pool did not terminate");
            }
            threadPool.shutdownNow();
            Thread.currentThread().interrupt();

            sqlConnection.close();
            if (lock.isLocked()) {
                lock.unlock();
            }

            if (resultCount == 0) {
                logger.error("There is no useful results.");
                System.exit(1);
            }
        } catch (Exception ex) {
            ex.printStackTrace();
            logger.error(ex.toString());
            System.exit(1);
        }

        logger.info("Estimating FDR...");
        String percolatorInputFileName = spectraPath + "." + labeling + ".input.temp";
        String percolatorOutputFileName = spectraPath + "." + labeling + ".output.temp";
        writePercolator(percolatorInputFileName, buildIndexObj.getPeptide0Map(), sqlPath);
        Map<Integer, PercolatorEntry> percolatorResultMap = runPercolator(percolatorPath, percolatorInputFileName, percolatorOutputFileName);

        if (percolatorResultMap.isEmpty()) {
            logger.error("Percolator failed to estimate FDR. Please check if Percolator is installed and the percolator_path in {} is correct.", parameterPath);
            (new File("PIPI.temp.db")).delete();
            System.exit(1);
        }

        if (!outputPercolatorInput) {
            (new File(percolatorInputFileName)).delete();
            (new File(percolatorOutputFileName)).delete();
        }

        logger.info("Saving results...");
        writeFinalResult(percolatorResultMap, spectraPath + "." + labeling + ".pipi.csv", buildIndexObj.getPeptide0Map(), sqlPath);
        new WritePepXml(spectraPath + "." + labeling + ".pipi.pep.xml", spectraPath, parameterMap, massToolObj.returnMassTable(), percolatorResultMap, buildIndexObj.getPeptide0Map(), buildIndexObj.returnFixModMap(), sqlPath);

        (new File("PIPI.temp.db")).delete();

        logger.info("Done.");
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

    private void writePercolator(String resultPath, Map<String, Peptide0> peptide0Map, String sqlPath) {
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(resultPath));
            writer.write("id\tlabel\tscannr\tscore\tdelta_c\tdelta_L_c\tnormalized_cross_corr\tglobal_search_rank\tabs_ppm\tion_frac\tmatched_high_peak_frac\tcharge1\tcharge2\tcharge3\tcharge4\tcharge5\tcharge6\texplained_aa_frac\tptm_supporting_peak_frac\tpeptide\tprotein\n");
            Connection sqlConnection = DriverManager.getConnection(sqlPath);
            Statement sqlStatement = sqlConnection.createStatement();
            ResultSet sqlResultSet = sqlStatement.executeQuery("SELECT scanNum, precursorCharge, precursorMass, peptide, theoMass, isDecoy, globalRank, normalizedCorrelationCoefficient, score, deltaLC, deltaC, matchedPeakNum, ionFrac, matchedHighestIntensityFrac, explainedAaFrac, ptmSupportingPeakFrac FROM spectraTable");
            while (sqlResultSet.next()) {
                String peptide = sqlResultSet.getString("peptide");
                if (!sqlResultSet.wasNull()) {
                    int charge = sqlResultSet.getInt("precursorCharge");
                    float theoMass = sqlResultSet.getFloat("theoMass");
                    float expMass = sqlResultSet.getFloat("precursorMass");
                    float massDiff = getMassDiff(expMass, theoMass, MassTool.C13_DIFF);

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
                    double deltaLC = sqlResultSet.getDouble("deltaLC");
                    double deltaC = sqlResultSet.getDouble("deltaC");
                    double ionFrac = sqlResultSet.getDouble("ionFrac");
                    double matchedHighestIntensityFrac = sqlResultSet.getDouble("matchedHighestIntensityFrac");
                    double explainedAaFrac = sqlResultSet.getDouble("explainedAaFrac");
                    double ptmSupportingPeakFrac = sqlResultSet.getDouble("ptmSupportingPeakFrac");

                    if (isDecoy == 1) {
                        writer.write(scanNum + "\t-1\t" + scanNum + "\t" + score + "\t" + deltaC + "\t" + deltaLC + "\t" + normalizedCorrelationCoefficient + "\t" + globalRank + "\t" + Math.abs(massDiff * 1e6f / theoMass) + "\t" + ionFrac + "\t" + matchedHighestIntensityFrac + "\t" + sb.toString() + explainedAaFrac + "\t" + ptmSupportingPeakFrac + "\t" + peptide0.leftFlank + "." + peptide + "." + peptide0.rightFlank + "\t" + String.join(";", proteinIdSet) + "\n");
                    } else {
                        writer.write(scanNum + "\t1\t" + scanNum + "\t" + score + "\t" + deltaC + "\t" + deltaLC + "\t" + normalizedCorrelationCoefficient + "\t" + globalRank + "\t" + Math.abs(massDiff * 1e6f / theoMass) + "\t" + ionFrac + "\t" + matchedHighestIntensityFrac + "\t" + sb.toString() + explainedAaFrac + "\t" + ptmSupportingPeakFrac + "\t" + peptide0.leftFlank + "." + peptide + "." + peptide0.rightFlank + "\t" + String.join(";", proteinIdSet) + "\n");
                    }
                }
            }
            writer.close();
            sqlResultSet.close();
            sqlStatement.close();
            sqlConnection.close();
        } catch (Exception ex) {
            ex.printStackTrace();
            logger.error(ex.toString());
            System.exit(1);
        }
    }

    private static Map<Integer, PercolatorEntry> runPercolator(String percolatorPath, String percolatorInputFileName, String percolatorOutputFileName) {
        Map<Integer, PercolatorEntry> percolatorResultMap = new HashMap<>();
        try {
            if ((new File(percolatorPath)).exists()) {
                Process ps = Runtime.getRuntime().exec(percolatorPath + " --only-psms --verbose 1 --no-terminate --results-psms " + percolatorOutputFileName + " " + percolatorInputFileName);
                ps.waitFor();

                if (!(new File(percolatorOutputFileName).exists())) {
                    logger.error("There is no Percolator output file.");
                    BufferedReader reader = new BufferedReader(new InputStreamReader(ps.getInputStream()));
                    String line;
                    while ((line = reader.readLine()) != null) {
                        logger.info("[Percolator info]: {}", line);
                    }
                    reader.close();
                    reader = new BufferedReader(new InputStreamReader(ps.getErrorStream()));
                    while ((line = reader.readLine()) != null) {
                        logger.error("[Percolator error]: {}", line);
                    }
                    reader.close();
                    return percolatorResultMap;
                } else {
                    BufferedReader reader = new BufferedReader(new InputStreamReader(ps.getInputStream()));
                    String line;
                    while ((line = reader.readLine()) != null) {
                        logger.info("[Percolator info]: {}", line);
                    }
                    reader.close();
                    reader = new BufferedReader(new InputStreamReader(ps.getErrorStream()));
                    while ((line = reader.readLine()) != null) {
                        logger.error("[Percolator error]: {}", line);
                    }
                    reader.close();
                }

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
                logger.error("Cannot find Percolator (from {}) for estimating Percolator Q-Value.", percolatorPath);
                return percolatorResultMap;
            }
        } catch (IOException | InterruptedException ex) {
            ex.printStackTrace();
            logger.error(ex.toString());
            return percolatorResultMap;
        }

        return percolatorResultMap;
    }

    private void writeFinalResult(Map<Integer, PercolatorEntry> percolatorResultMap, String outputPath, Map<String, Peptide0> peptide0Map, String sqlPath) {
        TreeMap<Double, List<String>> tempMap = new TreeMap<>();
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(outputPath));
            writer.write("scan_num,peptide,charge,theo_mass,exp_mass,abs_ppm,ptm_delta_score,ptm_supporting_peak_frac,protein_ID,score,percolator_score,posterior_error_prob,q_value,other_PTM_patterns,MGF_title,labeling,isotope_correction,MS1_pearson_correlation_coefficient\n");
            Connection sqlConnection = DriverManager.getConnection(sqlPath);
            Statement sqlStatement = sqlConnection.createStatement();
            ResultSet sqlResultSet = sqlStatement.executeQuery("SELECT scanNum, precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, labelling, peptide, theoMass, isDecoy, score, ptmSupportingPeakFrac, otherPtmPatterns, ptmDeltaScore FROM spectraTable");
            while (sqlResultSet.next()) {
                int isDecoy = sqlResultSet.getInt("isDecoy");
                if (!sqlResultSet.wasNull()) {
                    if (isDecoy == 0) {
                        int scanNum = sqlResultSet.getInt("scanNum");
                        float expMass = sqlResultSet.getFloat("precursorMass");
                        String peptide = sqlResultSet.getString("peptide");
                        float theoMass = sqlResultSet.getFloat("theoMass");
                        float massDiff = getMassDiff(expMass, theoMass, MassTool.C13_DIFF);
                        float ppm = Math.abs(massDiff * 1e6f / theoMass);

                        Peptide0 peptide0 = peptide0Map.get(peptide.replaceAll("[^ncA-Z]+", ""));
                        TreeSet<String> proteinIdSet = new TreeSet<>();
                        for (String protein : peptide0.proteins) {
                            proteinIdSet.add(protein.trim());
                        }

                        String ptmDeltaScore = sqlResultSet.getString("ptmDeltaScore");

                        PercolatorEntry percolatorEntry = percolatorResultMap.get(scanNum);
                        String str = String.format(Locale.US, "%d,%s,%d,%.4f,%.4f,%.2f,%s,%s,%s,%.4f,%.4f,%s,%s,%s,\"%s\",%s,%d,%f\n", scanNum, peptide, sqlResultSet.getInt("precursorCharge"), theoMass, expMass, ppm, ptmDeltaScore, ptmDeltaScore.contentEquals("-") ? "-" : String.format(Locale.US, "%.4f", sqlResultSet.getDouble("ptmSupportingPeakFrac")), String.join(";", proteinIdSet), sqlResultSet.getDouble("score"), percolatorEntry.percolatorScore, percolatorEntry.PEP, percolatorEntry.qValue, sqlResultSet.getString("otherPtmPatterns"), sqlResultSet.getString("mgfTitle"), sqlResultSet.getString("labelling"), sqlResultSet.getInt("isotopeCorrectionNum"), sqlResultSet.getDouble("ms1PearsonCorrelationCoefficient"));

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

            sqlResultSet.close();
            sqlStatement.close();
            sqlConnection.close();

            Double[] tempArray = tempMap.keySet().toArray(new Double[tempMap.size()]);
            for (int i = tempArray.length - 1; i >= 0; --i) {
                List<String> tempList = tempMap.get(tempArray[i]);
                for (String tempStr : tempList) {
                    writer.write(tempStr);
                }
            }
            writer.close();
        } catch (Exception ex) {
            ex.printStackTrace();
            logger.error(ex.toString());
            System.exit(1);
        }
    }

    public static float getMassDiff(float expMass, float theoMass, float C13Diff) {
        float massDiff1 = expMass - theoMass;
        float massDiff2 = expMass - theoMass - C13Diff;
        float massDiff3 = expMass - theoMass - 2 * C13Diff;
        float absMassDiff1 = Math.abs(massDiff1);
        float absMassDiff2 = Math.abs(massDiff2);
        float absMassDiff3 = Math.abs(massDiff3);

        if ((absMassDiff1 <= absMassDiff2) && (absMassDiff1 <= absMassDiff2)) {
            return massDiff1;
        } else if ((absMassDiff2 <= absMassDiff1) && (absMassDiff2 <= absMassDiff3)) {
            return massDiff2;
        } else {
            return massDiff3;
        }
    }
}
