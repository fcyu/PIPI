package proteomics;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.FDR.EstimateFDR;
import proteomics.Segment.InferenceSegment;
import proteomics.Types.*;
import proteomics.Index.BuildIndex;
import proteomics.Parameter.Parameter;
import proteomics.Spectrum.PreSpectra;
import proteomics.TheoSeq.MassTool;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLParsingException;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

public class PIPI {

    private static final Logger logger = LoggerFactory.getLogger(PIPI.class);
    private static final float C13Diff = 1.00335483f;
    private static final String versionStr = "1.2.8";

    public static final boolean DEV = false;

    public static void main(String args[]) {
        // Process inputs
        if (args.length != 2) {
            help();
        }

        // Set parameters
        String parameterPath = args[0].trim();
        String spectraPath = args[1].trim();

        try {
            logger.info("Running PIPI version {}.", versionStr);
            logger.info("Spectra: {}, parameter: {}", spectraPath, parameterPath);

            if (DEV) {
                logger.info("In DEV mode.");
            }

            new PIPI(parameterPath, spectraPath);
        } catch (Exception ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }
    }

    private PIPI(String parameterPath, String spectraPath) throws Exception {
        // Get the parameter map
        Parameter parameter = new Parameter(parameterPath);
        Map<String, String> parameterMap = parameter.returnParameterMap();
        float ms2Tolerance = Float.valueOf(parameterMap.get("ms2_tolerance"));
        float ms1Tolerance = Float.valueOf(parameterMap.get("ms1_tolerance"));
        int ms1ToleranceUnit = Integer.valueOf(parameterMap.get("ms1_tolerance_unit"));
        int maxMs2Charge = Integer.valueOf(parameterMap.get("max_ms2_charge"));
        String percolatorPath = parameterMap.get("percolator_path");
        boolean outputPercolatorInput = (DEV || (Integer.valueOf(parameterMap.get("output_percolator_input")) == 1));
        int minPotentialCharge = Integer.valueOf(parameterMap.get("min_potential_charge"));
        int maxPotentialCharge = Integer.valueOf(parameterMap.get("max_potential_charge"));
        float minPrecursorMass = Float.valueOf(parameterMap.get("min_precursor_mass"));
        float maxPrecursorMass = Float.valueOf(parameterMap.get("max_precursor_mass"));

        logger.info("Indexing protein database...");
        BuildIndex buildIndexObj = new BuildIndex(parameterMap);
        MassTool massToolObj = buildIndexObj.returnMassToolObj();

        InferenceSegment inference3SegmentObj = new InferenceSegment(buildIndexObj, ms2Tolerance, 3);

        // coding peptide beforehand
        Map<String, Peptide0> peptideCodeMap = new HashMap<>();
        Map<String, Float> peptideMassMap = buildIndexObj.returnPepMassMap();
        for (String peptide : peptideMassMap.keySet()) {
            SparseBooleanVector code = inference3SegmentObj.generateSegmentBooleanVector(inference3SegmentObj.cutTheoSegment(peptide));
            peptideCodeMap.put(peptide, new Peptide0(peptideMassMap.get(peptide), code, code.norm2square(), true));
        }
        peptideMassMap = buildIndexObj.returnDecoyPepMassMap();
        for (String peptide : peptideMassMap.keySet()) {
            SparseBooleanVector code = inference3SegmentObj.generateSegmentBooleanVector(inference3SegmentObj.cutTheoSegment(peptide));
            peptideCodeMap.put(peptide, new Peptide0(peptideMassMap.get(peptide), code, code.norm2square(), false));
        }

        logger.info("Reading PTM database...");
        Map<String, TreeSet<Integer>> siteMass1000Map = readPTMDb(parameterMap, buildIndexObj.returnFixModMap());

        int tempMax = 0;
        int tempMin = 99999;
        for (TreeSet<Integer> tempSet : siteMass1000Map.values()) {
            if (tempSet.last() > tempMax) {
                tempMax = tempSet.last();
            }
            if (tempSet.first() < tempMin) {
                tempMin = tempSet.first();
            }
        }

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
            } else if (ext.contentEquals("mgf")) {
                spectraParser = new MgfFile(spectraFile);
            } else {
                logger.error("Unsupported file format {}.", ext);
                logger.error("Currently, PIPI only support mzXML and MGF.");
            }
        } catch (FileNotFoundException | MzXMLParsingException ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }

        PreSpectra preSpectraObj = new PreSpectra(spectraParser, parameterMap, massToolObj, ext);
        Map<Integer, SpectrumEntry> numSpectrumMap = preSpectraObj.returnNumSpectrumMap();
        logger.info("Useful MS/MS spectra number: {}", numSpectrumMap.size());

        logger.info("Start searching...");

        int threadNum = Integer.valueOf(parameterMap.get("thread_num"));
        if (threadNum == 0) {
            threadNum = 1 + Runtime.getRuntime().availableProcessors();
        }
        ExecutorService threadPool = Executors.newFixedThreadPool(threadNum);

        List<Future<FinalResultEntry>> tempResultList = new LinkedList<>();
        for (int scanNum : numSpectrumMap.keySet()) {
            SpectrumEntry spectrumEntry = numSpectrumMap.get(scanNum);
            if (spectrumEntry.precursorCharge > 0) {
                tempResultList.add(threadPool.submit(new PIPIWrap(buildIndexObj, massToolObj, inference3SegmentObj, spectrumEntry, siteMass1000Map, peptideCodeMap, ms1Tolerance, ms1ToleranceUnit, ms2Tolerance, minPtmMass, maxPtmMass, maxMs2Charge)));
            } else {
                for (int potentialCharge = minPotentialCharge; potentialCharge <= maxPotentialCharge; ++potentialCharge) {
                    float potentialPrecursorMass = potentialCharge * (spectrumEntry.precursorMz - 1.00727646688f);
                    if ((potentialPrecursorMass >= minPrecursorMass) && (potentialPrecursorMass <= maxPrecursorMass)) {
                        SpectrumEntry fakeSpectrumEntry = new SpectrumEntry(scanNum, spectrumEntry.precursorMz, potentialPrecursorMass, potentialCharge, new TreeMap<>(spectrumEntry.plMap.subMap(0f, true, potentialPrecursorMass, true)));
                        tempResultList.add(threadPool.submit(new PIPIWrap(buildIndexObj, massToolObj, inference3SegmentObj, fakeSpectrumEntry, siteMass1000Map, peptideCodeMap, ms1Tolerance, ms1ToleranceUnit, ms2Tolerance, minPtmMass, maxPtmMass, maxMs2Charge)));
                    }
                }
            }
        }

        // check progress every minute
        int lastProgress = 0;
        try {
            int count;
            while (((count = finishedFutureNum(tempResultList)) < tempResultList.size()) && hasUnfinishedFuture(tempResultList)) {
                int progress = count * 20 / tempResultList.size();
                if (progress != lastProgress) {
                    logger.info("Searching {}%...", progress * 5);
                    lastProgress = progress;
                }
                Thread.sleep(60000);
            }
        } catch (InterruptedException ex) {
            ex.printStackTrace();
            logger.error(ex.toString());
            System.exit(1);
        }

        logger.info("Searching 100%...");

        // record search results.
        List<FinalResultEntry> tempList = new LinkedList<>();
        try {
            for (Future<FinalResultEntry> tempResult : tempResultList) {
                if (tempResult.isDone() && !tempResult.isCancelled()) {
                    if (tempResult.get() != null) {
                        tempList.add(tempResult.get());
                    }
                } else {
                    logger.error("Threads were not finished normally.");
                    System.exit(1);
                }
            }
        } catch (Exception ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
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

        if (tempList.isEmpty()) {
            logger.error("There is no useful results.");
            System.exit(1);
        }

        // merge different results corresponding to different charges for the same spectrum.
        List<FinalResultEntry> finalScoredPsms = new LinkedList<>();
        for(FinalResultEntry e1 : tempList) {
            boolean keep = true;
            for (FinalResultEntry e2 : tempList) {
                if ((e1.getScanNum() == e2.getScanNum()) && (e1.getNegativeLog10EValue() < e2.getNegativeLog10EValue())) { // Compare e-values because different charge states result in different search space, which result in different e-values even for the same XCorr score.
                    keep = false;
                    break;
                }
            }
            if (keep) {
                finalScoredPsms.add(e1);
            }
        }

        logger.info("Estimating FDR...");
        // estimate T-D FDR
        new EstimateFDR(finalScoredPsms);

        // estimate Percolator FDR
        String percolatorInputFileName = spectraPath + ".input.temp";
        String percolatorOutputFileName = spectraPath + ".output.temp";
        Map<String, Set<String>> peptideProteinMap = buildIndexObj.returnPepProMap();
        Map<String, String> decoyPeptideProteinMap = buildIndexObj.returnDecoyPepProMap();
        writePercolator(finalScoredPsms, peptideProteinMap, decoyPeptideProteinMap, percolatorInputFileName, buildIndexObj.returnFixModMap());
        Map<Integer, PercolatorEntry> percolatorResultMap = runPercolator(percolatorPath, percolatorInputFileName, percolatorOutputFileName);

        if (percolatorResultMap.isEmpty()) {
            logger.warn("Percolator failed to estimate FDR. The results won't contain percolator_score, posterior_error_prob, and percolator_q_value.");
        }

        if (!outputPercolatorInput) {
            (new File(percolatorInputFileName)).delete();
            (new File(percolatorOutputFileName)).delete();
        }

        logger.info("Saving results...");
        writeFinalResult(finalScoredPsms, percolatorResultMap, peptideProteinMap, spectraPath + ".pipi.csv", buildIndexObj.returnFixModMap());

        logger.info("Done.");
    }

    private static void help() {
        String helpStr = "PIPI version " + versionStr + "\r\n"
                + "A tool identifying peptides with unlimited PTM.\r\n"
                + "Author: Fengchao Yu\r\n"
                + "Email: fyuab@connect.ust.hk\r\n"
                + "ECL usage: java -Xmx25g -jar /path/to/PIPI.jar <parameter_file> <data_file>\r\n"
                + "\t<parameter_file>: parameter file. Can be download along with PIPI.\r\n"
                + "\t<data_file>: spectra data file (mzXML)\r\n"
                + "\texample: java -Xmx32g -jar ECL.jar parameter.def data.mzxml\r\n";
        System.out.print(helpStr);
        System.exit(1);
    }

    private static void writePercolator(List<FinalResultEntry> finalScoredResult, Map<String, Set<String>> peptideProteinMap, Map<String, String> decoyPeptideProteinMap, String resultPath, Map<String, Float> fixModMap) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(resultPath))) {
            writer.write("id\tlabel\tscannr\txcorr\tdelta_c\tdelta_L_c\tnegative_log10_e_value\tnormalized_cross_corr\tglobal_search_rank\tabs_ppm\tIonFrac\tmatched_high_peak_frac\tcharge1\tcharge2\tcharge3\tcharge4\tcharge5\tcharge6\tvar_PTM_num\tpeptide\tprotein\n");
            for (FinalResultEntry entry : finalScoredResult) {
                float expMass = entry.getCharge() * (entry.getPrecursorMz() - 1.00727646688f);
                Peptide peptide = entry.getPeptide();
                float theoMass = peptide.getPrecursorMass();
                float massDiff = getMassDiff(expMass, theoMass, C13Diff);
                String proteinIdStr = "";

                if (!entry.isDecoy()) {
                    for (String temp : peptideProteinMap.get(peptide.getPTMFreeSeq())) {
                        proteinIdStr += temp + ";";
                    }
                } else {
                    proteinIdStr = decoyPeptideProteinMap.get(peptide.getPTMFreeSeq());
                }

                StringBuilder sb = new StringBuilder(20);
                int charge = entry.getCharge();
                for (int i = 0; i < 6; ++i) {
                    if (i == charge - 1) {
                        sb.append(1);
                    } else {
                        sb.append(0);
                    }
                    sb.append("\t");
                }

                if (entry.isDecoy()) {
                    writer.write(entry.getScanNum() + "\t-1\t" + entry.getScanNum() + "\t" + entry.getScore() + "\t" + entry.getDeltaC() + "\t" + entry.getDeltaLC() + "\t" + entry.getNegativeLog10EValue() + "\t" + entry.getNormalizedCrossXcorr() + "\t" + entry.getGlobalSearchRank() + "\t" + Math.abs(massDiff * 1e6f / theoMass) + "\t" + entry.getIonFrac() + "\t" + entry.getMatchedHighestIntensityFrac() + "\t" + sb.toString() + peptide.getVarPTMNum() + "\t" + peptide.getLeftFlank() + "." + peptide.getPTMContainedString(fixModMap) + "." + peptide.getRightFlank() + "\t" + proteinIdStr + "\n");
                } else {
                    writer.write(entry.getScanNum() + "\t1\t" + entry.getScanNum() + "\t" + entry.getScore() + "\t" + entry.getDeltaC() + "\t" + entry.getDeltaLC() + "\t" + entry.getNegativeLog10EValue() + "\t" + entry.getNormalizedCrossXcorr() + "\t" + entry.getGlobalSearchRank() + "\t" + Math.abs(massDiff * 1e6f / theoMass) + "\t" + entry.getIonFrac() + "\t" + entry.getMatchedHighestIntensityFrac() + "\t" + sb.toString() + peptide.getVarPTMNum() + "\t" + peptide.getLeftFlank() + "." + peptide.getPTMContainedString(fixModMap) + "." + peptide.getRightFlank() + "\t" + proteinIdStr + "\n");
                }
            }
        } catch (IOException ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
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
                logger.error("Cannot find Percolator for estimating Percolator Q-Value. The results won't contain percolator_score, posterior_error_prob, and percolator_q_value.");
                return percolatorResultMap;
            }
        } catch (IOException | InterruptedException ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            return percolatorResultMap;
        }

        return percolatorResultMap;
    }

    private static void writeFinalResult(List<FinalResultEntry> finalScoredPsms, Map<Integer, PercolatorEntry> percolatorResultMap, Map<String, Set<String>> peptideProteinMap, String outputPath, Map<String, Float> fixModMap) {
        TreeMap<Double, List<String>> tempMap = new TreeMap<>();
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputPath))) {
            writer.write("scan_num,peptide,charge,theo_mass,exp_mass,ppm,protein_ID,xcorr,e_value,naive_q_value,percolator_score,posterior_error_prob,percolator_q_value\n");
            for (FinalResultEntry entry : finalScoredPsms) {
                if (!entry.isDecoy()) {
                    int scanNum = entry.getScanNum();
                    float expMass = entry.getCharge() * (entry.getPrecursorMz() - 1.00727646688f);
                    int charge = entry.getCharge();
                    Peptide peptide = entry.getPeptide();
                    float theoMass = peptide.getPrecursorMass();
                    float massDiff = getMassDiff(expMass, theoMass, C13Diff);
                    float ppm = Math.abs(massDiff * 1e6f / theoMass);


                    String proteinIdStr = "";
                    for (String tempStr : peptideProteinMap.get(peptide.getPTMFreeSeq())) {
                        proteinIdStr += tempStr + ";";
                    }

                    String str;
                    boolean sortedByPercolatorScore = true;
                    if (percolatorResultMap.containsKey(scanNum)) {
                        PercolatorEntry percolatorEntry = percolatorResultMap.get(scanNum);
                        str = String.format("%d,%s,%d,%.4f,%.4f,%.2f,%s,%.4f,%E,%.4f,%.3f,%s,%s\n", scanNum, peptide.getPTMContainedString(fixModMap), charge, theoMass, expMass, ppm, proteinIdStr, entry.getScore(), entry.getEValue(), entry.getQValue(), percolatorEntry.percolatorScore, percolatorEntry.PEP, percolatorEntry.qValue);
                    } else {
                        sortedByPercolatorScore = false;
                        str = String.format("%d,%s,%d,%.4f,%.4f,%.2f,%s,%.4f,%E,%.4f,%s,%s,%s\n", scanNum, peptide.getPTMContainedString(fixModMap), charge, theoMass, expMass, ppm, proteinIdStr, entry.getScore(), entry.getEValue(), entry.getQValue() , "-", "-", "-");
                    }

                    if (sortedByPercolatorScore) {
                        if (tempMap.containsKey(percolatorResultMap.get(scanNum).percolatorScore)) {
                            tempMap.get(percolatorResultMap.get(scanNum).percolatorScore).add(str);
                        } else {
                            List<String> tempList = new LinkedList<>();
                            tempList.add(str);
                            tempMap.put(percolatorResultMap.get(scanNum).percolatorScore, tempList);
                        }
                    } else {
                        if (tempMap.containsKey(entry.getNegativeLog10EValue())) {
                            tempMap.get(entry.getNegativeLog10EValue()).add(str);
                        } else {
                            List<String> tempList = new LinkedList<>();
                            tempList.add(str);
                            tempMap.put(entry.getNegativeLog10EValue(), tempList);
                        }
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
        } catch (IOException ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }
    }

    private static float getMassDiff(float expMass, float theoMass, float C13Diff) {
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

    private Map<String, TreeSet<Integer>> readPTMDb(Map<String, String> parameterMap, Map<String, Float> fixModMap) {
        Map<String, TreeSet<Integer>> siteMass1000Map = new HashMap<>();
        float minPtmMass = Float.valueOf(parameterMap.get("min_ptm_mass"));
        float maxPtmMass = Float.valueOf(parameterMap.get("max_ptm_mass"));

        try (BufferedReader reader = new BufferedReader(new FileReader(parameterMap.get("PTM_db")))) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.startsWith("site") || line.startsWith("#")) {
                    continue;
                }

                String[] parts = line.split("\t");
                String site = parts[0];
                String position = parts[3];
                float mass = Float.valueOf(parts[2]);

                if ((mass > maxPtmMass) || mass < minPtmMass) {
                    continue;
                }

                int mass1000 = (int) Math.floor(mass * 1000);

                String siteString;
                if (site.contentEquals("N-term") && (Math.abs(fixModMap.get("n")) < 1e-6)) {
                    if (position.contentEquals("PROTEIN_N")) {
                        siteString = "PROTEIN_N";
                    } else {
                        siteString = "PEPTIDE_N";
                    }
                } else if (site.contentEquals("C-term")) {
                    if (position.contentEquals("PROTEIN_C")) {
                        siteString = "PROTEIN_C";
                    } else {
                        siteString = "PEPTIDE_C";
                    }
                } else if (Math.abs(fixModMap.get(site)) < 1e-6) { // fix modified amino acid cannot be modified again.
                    if (position.contentEquals("PROTEIN_N")) {
                        siteString = site + "-PROTEIN_N";
                    } else if (position.contentEquals("PROTEIN_C")) {
                        siteString = site + "-PROTEIN_C";
                    } else if (position.contentEquals("PEPTIDE_N")) {
                        siteString = site + "-PEPTIDE_N";
                    } else if (position.contentEquals("PEPTIDE_C")) {
                        siteString = site + "-PEPTIDE_C";
                    } else {
                        siteString = site;
                    }
                } else {
                    continue;
                }

                if (siteMass1000Map.containsKey(siteString)) {
                    siteMass1000Map.get(siteString).add(mass1000);
                } else {
                    TreeSet<Integer> tempSet = new TreeSet<>();
                    tempSet.add(mass1000);
                    siteMass1000Map.put(siteString, tempSet);
                }
            }
        } catch (IOException ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }

        return siteMass1000Map;
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
