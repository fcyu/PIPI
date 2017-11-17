package proteomics;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Output.WritePepXml;
import proteomics.Types.*;
import proteomics.Index.BuildIndex;
import proteomics.Parameter.Parameter;
import proteomics.Spectrum.PreSpectra;
import proteomics.TheoSeq.MassTool;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.JMzReaderException;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLParsingException;

import java.io.*;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.util.*;
import java.util.concurrent.*;

public class PIPI {

    private static final Logger logger = LoggerFactory.getLogger(PIPI.class);
    public static final String versionStr = "1.4.2";
    public static final boolean useXcorr = true;

    public static final boolean DEV = false;
    public static final int[] debugScanNumArray = new int[]{};

    private Map<Character, Float> fixModMap = null;

    public static void main(String[] args) {
        // Process inputs
        if (args.length != 2) {
            help();
        }

        // Set parameters
        String parameterPath = args[0].trim();
        String spectraPath = args[1].trim();

        try {
            logger.info("Running PIPI version {}.", versionStr);

            String hostName = InetAddress.getLocalHost().getHostName();
            logger.info("Computer: {}.", hostName);

            logger.info("Spectra: {}, parameter: {}.", spectraPath, parameterPath);

            if (DEV) {
                logger.info("In DEV mode.");
            }

            new PIPI(parameterPath, spectraPath);
        } catch (UnknownHostException ex) {
            logger.warn("Cannot get the computer's name.");
        } catch (Exception ex) {
            ex.printStackTrace();
            logger.error(ex.toString());
            System.exit(1);
        }
    }

    private PIPI(String parameterPath, String spectraPath) {
        // Get the parameter map
        Parameter parameter = new Parameter(parameterPath);
        Map<String, String> parameterMap = parameter.returnParameterMap();
        float ms2Tolerance = Float.valueOf(parameterMap.get("ms2_tolerance"));
        float ms1Tolerance = Float.valueOf(parameterMap.get("ms1_tolerance"));
        int ms1ToleranceUnit = Integer.valueOf(parameterMap.get("ms1_tolerance_unit"));
        int maxMs2Charge = Integer.valueOf(parameterMap.get("max_ms2_charge"));
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
        fixModMap = buildIndexObj.returnFixModMap();

        float minPtmMass = Float.valueOf(parameterMap.get("min_ptm_mass"));
        float maxPtmMass = Float.valueOf(parameterMap.get("max_ptm_mass"));

        logger.info("Reading spectra...");
        JMzReader spectraParser = null;
        String ext = "";
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
                logger.error("Unsupported file format {}.", ext);
                logger.error("Currently, PIPI only support mzXML and MGF.");
            }
        } catch (FileNotFoundException | MzXMLParsingException | JMzReaderException ex) {
            ex.printStackTrace();
            logger.error(ex.toString());
            System.exit(1);
        }

        PreSpectra preSpectraObj = new PreSpectra(spectraParser, parameterMap, massToolObj, ext, msLevelSet);
        Map<Integer, SpectrumEntry> numSpectrumMap = preSpectraObj.returnNumSpectrumMap();
        logger.info("Useful MS/MS spectra number: {}.", numSpectrumMap.size());

        if (DEV) {
            try (BufferedWriter writer = new BufferedWriter(new FileWriter("spectrum.dev.csv"))) {
                writer.write("scanNum,charge,finalIsotopeCorrectionNum,isotopeCorrectionNum,pearsonCorrelationCoefficient,expMz1,expMz2,expMz3,expInt1,expInt2,expInt3,theoMz1,theoMz2,theoMz3,theoInt1,theoInt2,theoInt3\n");
                for (SpectrumEntry spectrumEntry : numSpectrumMap.values()) {
                    for (int charge : spectrumEntry.chargeDevEntryMap.keySet()) {
                        for (SpectrumEntry.DevEntry devEntry : spectrumEntry.chargeDevEntryMap.get(charge)) {
                            writer.write(String.format(Locale.US, "%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", spectrumEntry.scanNum, charge, spectrumEntry.isotopeCorrectionNum, devEntry.isotopeCorrectionNum, devEntry.pearsonCorrelationCoefficient, devEntry.expMatrix[0][0], devEntry.expMatrix[1][0], devEntry.expMatrix[2][0], devEntry.expMatrix[0][1], devEntry.expMatrix[1][1], devEntry.expMatrix[2][1], devEntry.theoMatrix[0][0], devEntry.theoMatrix[1][0], devEntry.theoMatrix[2][0], devEntry.theoMatrix[0][1], devEntry.theoMatrix[1][1], devEntry.theoMatrix[2][1]));
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
            threadNum = 1 + Runtime.getRuntime().availableProcessors();
        }
        ExecutorService threadPool = Executors.newFixedThreadPool(threadNum);

        List<Future<FinalResultEntry>> taskList = new LinkedList<>();
        for (int scanNum : numSpectrumMap.keySet()) {
            SpectrumEntry spectrumEntry = numSpectrumMap.get(scanNum);
            taskList.add(threadPool.submit(new PIPIWrap(buildIndexObj, massToolObj, spectrumEntry, ms1Tolerance, ms1ToleranceUnit, ms2Tolerance, minPtmMass, maxPtmMass, Math.min(spectrumEntry.precursorCharge > 1 ? spectrumEntry.precursorCharge - 1 : 1, maxMs2Charge))));
        }

        // check progress every minute, record results,and delete finished tasks.
        Map<Integer, FinalResultEntry> scanFinalResultMap = new HashMap<>();
        int lastProgress = 0;
        try {
            int totalCount = taskList.size();
            int count = 0;
            while (count < totalCount) {
                // record search results and delete finished ones.
                List<Future<FinalResultEntry>> toBeDeleteTaskList = new LinkedList<>();
                for (Future<FinalResultEntry> task : taskList) {
                    if (task.isDone()) {
                        if (task.get() != null) {
                            FinalResultEntry finalResultEntry = task.get();
                            if (!scanFinalResultMap.containsKey(finalResultEntry.getScanNum()) || (scanFinalResultMap.get(finalResultEntry.getScanNum()).getTopPeptide().getScore() < finalResultEntry.getTopPeptide().getScore()) || (scanFinalResultMap.get(finalResultEntry.getScanNum()).getTopPeptide().getScore() == finalResultEntry.getTopPeptide().getScore() && !finalResultEntry.getTopPeptide().isDecoy()))
                            scanFinalResultMap.put(finalResultEntry.getScanNum(), finalResultEntry);
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
        } catch (InterruptedException | ExecutionException ex) {
            ex.printStackTrace();
            logger.error(ex.toString());
            System.exit(1);
        }

        // shutdown threads.
        threadPool.shutdown();
        try {
            if (!threadPool.awaitTermination(60, TimeUnit.SECONDS)) {
                threadPool.shutdownNow();
                if (!threadPool.awaitTermination(60, TimeUnit.SECONDS))
                    System.err.println("Pool did not terminate");
            }
        } catch (InterruptedException ie) {
            threadPool.shutdownNow();
            Thread.currentThread().interrupt();
            logger.error("Threads were not finished normally.");
            System.exit(1);
        }

        if (scanFinalResultMap.isEmpty()) {
            logger.error("There is no useful results.");
            System.exit(1);
        }

        logger.info("Estimating FDR...");
        String percolatorInputFileName = spectraPath + "." + labeling + ".input.temp";
        String percolatorOutputFileName = spectraPath + "." + labeling + ".output.temp";
        writePercolator(scanFinalResultMap, percolatorInputFileName, buildIndexObj.getPeptide0Map());
        Map<Integer, PercolatorEntry> percolatorResultMap = runPercolator(percolatorPath, percolatorInputFileName, percolatorOutputFileName);

        if (percolatorResultMap.isEmpty()) {
            logger.warn("Percolator failed to estimate FDR. The results won't contain percolator_score, posterior_error_prob, and percolator_q_value.");
        }

        if (!outputPercolatorInput) {
            (new File(percolatorInputFileName)).delete();
            (new File(percolatorOutputFileName)).delete();
        }

        logger.info("Saving results...");
        writeFinalResult(scanFinalResultMap, percolatorResultMap, spectraPath + "." + labeling + ".pipi.csv", buildIndexObj.getPeptide0Map());
        new WritePepXml(scanFinalResultMap, spectraPath + "." + labeling + ".pipi.pep.xml", spectraPath, parameterMap, massToolObj.returnMassTable(), percolatorResultMap, buildIndexObj.getPeptide0Map());

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

    private void writePercolator(Map<Integer, FinalResultEntry> scanFinalResultMap, String resultPath, Map<String, Peptide0> peptide0Map) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(resultPath))) {
            writer.write("id\tlabel\tscannr\tscore\tdelta_c\tdelta_L_c\tnormalized_cross_corr\tglobal_search_rank\tabs_ppm\tion_frac\tmatched_high_peak_frac\tcharge1\tcharge2\tcharge3\tcharge4\tcharge5\tcharge6\texplained_aa_frac\tptm_supporting_peak_frac\tpeptide\tprotein\n");
            for (FinalResultEntry entry : scanFinalResultMap.values()) {
                Peptide peptide = entry.getTopPeptide();
                float theoMass = peptide.getTheoMass();
                float expMass = entry.getCharge() * (entry.getPrecursorMz() - 1.00727646688f);
                float massDiff = getMassDiff(expMass, theoMass, MassTool.C13_DIFF);

                Peptide0 peptide0 = peptide0Map.get(peptide.getPTMFreeSeq());
                StringBuilder sb = new StringBuilder(peptide0.proteins.size()* 10);
                for (String protein : peptide0.proteins) {
                    sb.append(protein);
                    sb.append(";");
                }
                String proteinIdStr = sb.toString();

                sb = new StringBuilder(20);
                int charge = entry.getCharge();
                for (int i = 0; i < 6; ++i) {
                    if (i == charge - 1) {
                        sb.append(1);
                    } else {
                        sb.append(0);
                    }
                    sb.append("\t");
                }

                double deltaLC = (peptide.getScore() - entry.getPeptideSet().last().getScore()) / peptide.getScore();
                double deltaC = 0;
                if (entry.getPeptideSet().size() > 1) {
                    Iterator<Peptide> temp = entry.getPeptideSet().iterator();
                    temp.next();
                    deltaC = (peptide.getScore() - temp.next().getScore()) / peptide.getScore();
                }

                if (peptide.isDecoy()) {
                    writer.write(entry.getScanNum() + "\t-1\t" + entry.getScanNum() + "\t" + peptide.getScore() + "\t" + deltaC + "\t" + deltaLC + "\t" + peptide.getNormalizedCrossCorr() + "\t" + peptide.getGlobalRank() + "\t" + Math.abs(massDiff * 1e6f / theoMass) + "\t" + peptide.getIonFrac() + "\t" + peptide.getMatchedHighestIntensityFrac() + "\t" + sb.toString() + peptide.getExplainedAaFrac() + "\t" + peptide.getPtmSupportingPeakFrac() + "\t" + peptide0.leftFlank + "." + peptide.getPtmContainingSeq(fixModMap) + "." + peptide0.rightFlank + "\t" + proteinIdStr + "\n");
                } else {
                    writer.write(entry.getScanNum() + "\t1\t" + entry.getScanNum() + "\t" + peptide.getScore() + "\t" + deltaC + "\t" + deltaLC + "\t" + peptide.getNormalizedCrossCorr() + "\t" + peptide.getGlobalRank() + "\t" + Math.abs(massDiff * 1e6f / theoMass) + "\t" + peptide.getIonFrac() + "\t" + peptide.getMatchedHighestIntensityFrac() + "\t" + sb.toString() + peptide.getExplainedAaFrac() + "\t" + peptide.getPtmSupportingPeakFrac() + "\t" + peptide0.leftFlank + "." + peptide.getPtmContainingSeq(fixModMap) + "." + peptide0.rightFlank + "\t" + proteinIdStr + "\n");
                }
            }
        } catch (IOException | NullPointerException ex) {
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
                    logger.warn("Error in running Percolator. The results won't contain percolator_score, posterior_error_prob, and percolator_q_value.");
                    BufferedReader reader = new BufferedReader(new InputStreamReader(ps.getInputStream()));
                    String line;
                    while ((line = reader.readLine()) != null) {
                        logger.error("[Percolator info]: {}", line);
                    }
                    reader.close();
                    reader = new BufferedReader(new InputStreamReader(ps.getErrorStream()));
                    while ((line = reader.readLine()) != null) {
                        logger.error("[Percolator info]: {}", line);
                    }
                    reader.close();
                    return percolatorResultMap;
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
                logger.error("Cannot find Percolator (from {}) for estimating Percolator Q-Value. The results won't contain percolator_score, posterior_error_prob, and percolator_q_value.", percolatorPath);
                return percolatorResultMap;
            }
        } catch (IOException | InterruptedException ex) {
            ex.printStackTrace();
            logger.error(ex.toString());
            return percolatorResultMap;
        }

        return percolatorResultMap;
    }

    private void writeFinalResult(Map<Integer, FinalResultEntry> scanFinalResultMap, Map<Integer, PercolatorEntry> percolatorResultMap, String outputPath, Map<String, Peptide0> peptide0Map) {
        TreeMap<Double, List<String>> tempMap = new TreeMap<>();
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputPath))) {
            writer.write("scan_num,peptide,charge,theo_mass,exp_mass,abs_ppm,ptm_delta_score,ptm_supporting_peak_frac,protein_ID,score,percolator_score,posterior_error_prob,q_value,other_PTM_patterns,MGF_title,labeling,isotope_correction,MS1_pearson_correlation_coefficient\n");
            for (FinalResultEntry entry : scanFinalResultMap.values()) {
                Peptide peptide = entry.getTopPeptide();
                if (!peptide.isDecoy()) {
                    int scanNum = entry.getScanNum();
                    float expMass = entry.getCharge() * (entry.getPrecursorMz() - 1.00727646688f);
                    int charge = entry.getCharge();
                    float theoMass = peptide.getTheoMass();
                    float massDiff = getMassDiff(expMass, theoMass, MassTool.C13_DIFF);
                    float ppm = Math.abs(massDiff * 1e6f / theoMass);

                    Peptide0 peptide0 = peptide0Map.get(peptide.getPTMFreeSeq());
                    StringBuilder sb = new StringBuilder(peptide0.proteins.size()* 10);
                    for (String protein : peptide0.proteins) {
                        sb.append(protein);
                        sb.append(";");
                    }

                    StringBuilder otherPtmPatterns = new StringBuilder(300);
                    String ptmDeltaScore;
                    if (entry.getPtmPatterns().containsKey(peptide.getPTMFreeSeq())) {
                        TreeSet<Peptide> tempTreeSet = entry.getPtmPatterns().get(peptide.getPTMFreeSeq());
                        Iterator<Peptide> tempTreeSetIterator = tempTreeSet.iterator();
                        tempTreeSetIterator.next();
                        while (tempTreeSetIterator.hasNext()) {
                            Peptide temp = tempTreeSetIterator.next();
                            otherPtmPatterns.append(String.format(Locale.US, "%s-%.4f;", temp.getPtmContainingSeq(fixModMap), temp.getScore()));
                        }
                        if (tempTreeSet.size() > 1) {
                            Iterator<Peptide> temp = tempTreeSet.iterator();
                            ptmDeltaScore = String.valueOf(temp.next().getScore() - temp.next().getScore());
                        } else {
                            ptmDeltaScore = String.valueOf(tempTreeSet.first().getScore());
                        }
                    } else {
                        otherPtmPatterns.append("-");
                        ptmDeltaScore = "-";
                    }

                    PercolatorEntry percolatorEntry = percolatorResultMap.get(scanNum);
                    String str = String.format(Locale.US, "%d,%s,%d,%.4f,%.4f,%.2f,%s,%s,%s,%.4f,%.4f,%s,%s,%s,\"%s\",%s,%d,%f\n", scanNum, peptide.getPtmContainingSeq(fixModMap), charge, theoMass, expMass, ppm, ptmDeltaScore, ptmDeltaScore.contentEquals("-") ? "-" : String.format("%.4f", peptide.getPtmSupportingPeakFrac()), sb.toString(), peptide.getScore(), percolatorEntry.percolatorScore, percolatorEntry.PEP, percolatorEntry.qValue, otherPtmPatterns.toString(), entry.getMgtTitle(), entry.getLabeling(), entry.getIsotopeCorrectionNum(), entry.getMs1PearsonCorrelationCoefficient());

                    if (tempMap.containsKey(percolatorResultMap.get(scanNum).percolatorScore)) {
                        tempMap.get(percolatorResultMap.get(scanNum).percolatorScore).add(str);
                    } else {
                        List<String> tempList = new LinkedList<>();
                        tempList.add(str);
                        tempMap.put(percolatorResultMap.get(scanNum).percolatorScore, tempList);
                    }
                }
            }

            Double[] tempArray = tempMap.keySet().toArray(new Double[tempMap.size()]);
            for (int i = tempArray.length - 1; i >= 0; --i) {
                List<String> tempList = tempMap.get(tempArray[i]);
                for (String tempStr : tempList) {
                    writer.write(tempStr);
                }
            }
        } catch (IOException | NullPointerException ex) {
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

    private static int finishedFutureNum(List<Future<FinalResultEntry>> futureList) {
        int count = 0;
        for (Future<FinalResultEntry> future : futureList) {
            if (future.isDone()) {
                ++count;
            }
        }
        return count;
    }

    private static boolean hasUnfinishedFuture(List<Future<FinalResultEntry>> futureList) {
        for (Future<FinalResultEntry> future : futureList) {
            if (!future.isDone()) {
                return true;
            }
        }
        return false;
    }
}
