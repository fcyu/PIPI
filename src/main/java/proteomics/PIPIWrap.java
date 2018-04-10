package proteomics;

import proteomics.Index.BuildIndex;
import proteomics.PTM.InferPTM;
import proteomics.Search.Binomial;
import proteomics.Search.CalSubscores;
import proteomics.Search.CalScore;
import proteomics.Search.Search;
import proteomics.Segment.InferenceSegment;
import ProteomicsLibrary.PrepareSpectrum;
import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Types.*;
import proteomics.Spectrum.PreSpectra;
import proteomics.Types.*;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;

import java.sql.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantLock;

public class PIPIWrap implements Callable<PIPIWrap.Entry> {

    private final BuildIndex buildIndexObj;
    private final MassTool massToolObj;
    private final double ms1Tolerance;
    private final double leftInverseMs1Tolerance;
    private final double rightInverseMs1Tolerance;
    private final int ms1ToleranceUnit;
    private final double ms2Tolerance;
    private final double minPtmMass;
    private final double maxPtmMass;
    private final int localMaxMs2Charge;
    private final Map<String, Peptide0> peptide0Map;
    private final JMzReader spectraParser;
    private final double minClear;
    private final double maxClear;
    private final ReentrantLock lock;
    private final String scanId;
    private final int precursorCharge;
    private final double precursorMass;
    private final InferPTM inferPTM;
    private final PrepareSpectrum preSpectrum;
    private final String sqlPath;
    private final Binomial binomial;


    public PIPIWrap(BuildIndex buildIndexObj, MassTool massToolObj, double ms1Tolerance, double leftInverseMs1Tolerance, double rightInverseMs1Tolerance, int ms1ToleranceUnit, double ms2Tolerance, double minPtmMass, double maxPtmMass, int localMaxMs2Charge, JMzReader spectraParser, double minClear, double maxClear, ReentrantLock lock, String scanId, int precursorCharge, double precursorMass, InferPTM inferPTM, PrepareSpectrum preSpectrum, String sqlPath, Binomial binomial) {
        this.buildIndexObj = buildIndexObj;
        this.massToolObj = massToolObj;
        this.ms1Tolerance = ms1Tolerance;
        this.leftInverseMs1Tolerance = leftInverseMs1Tolerance;
        this.rightInverseMs1Tolerance = rightInverseMs1Tolerance;
        this.ms1ToleranceUnit = ms1ToleranceUnit;
        this.ms2Tolerance = ms2Tolerance;
        this.minPtmMass = minPtmMass;
        this.maxPtmMass = maxPtmMass;
        this.localMaxMs2Charge = localMaxMs2Charge;
        this.spectraParser = spectraParser;
        this.minClear = minClear;
        this.maxClear = maxClear;
        this.lock = lock;
        this.scanId = scanId;
        this.precursorCharge = precursorCharge;
        this.precursorMass = precursorMass;
        this.inferPTM = inferPTM;
        this.preSpectrum = preSpectrum;
        this.sqlPath = sqlPath;
        this.binomial = binomial;
        peptide0Map = buildIndexObj.getPeptide0Map();
    }

    @Override
    public Entry call() throws Exception {
        Map<Double, Double> rawPLMap;
        try {
            lock.lock();
            // Reading peak list.
            rawPLMap = spectraParser.getSpectrumById(scanId).getPeakList();
        } finally {
            lock.unlock();
        }

        // preprocess peak list
        TreeMap<Double, Double> plMap = preSpectrum.preSpectrumTopNStyle(rawPLMap, precursorMass, precursorCharge, minClear, maxClear, PreSpectra.topN);

        if (plMap.isEmpty()) {
            return null;
        }

        // Coding
        InferenceSegment inference3SegmentObj = buildIndexObj.getInference3SegmentObj();
        List<ThreeExpAA> expAaLists = inference3SegmentObj.inferSegmentLocationFromSpectrum(precursorMass, plMap);
        if (!expAaLists.isEmpty()) {
            SparseVector scanCode = inference3SegmentObj.generateSegmentIntensityVector(expAaLists);

            // Begin search.
            Search searchObj = new Search(buildIndexObj, precursorMass, scanCode, massToolObj, ms1Tolerance, leftInverseMs1Tolerance, rightInverseMs1Tolerance, ms1ToleranceUnit, minPtmMass, maxPtmMass, localMaxMs2Charge);

            // prepare the spectrum
            SparseVector expProcessedPL;
            if (PIPI.useXcorr) {
                expProcessedPL = preSpectrum.prepareXCorr(plMap, false);
            } else {
                expProcessedPL = preSpectrum.digitizePL(plMap);
            }

            double localMS1ToleranceL = -1 * ms1Tolerance;
            double localMS1ToleranceR = ms1Tolerance;
            if (ms1ToleranceUnit == 1) {
                localMS1ToleranceL = (precursorMass * leftInverseMs1Tolerance) - precursorMass;
                localMS1ToleranceR = (precursorMass * rightInverseMs1Tolerance) - precursorMass;
            }

            // infer PTM using the new approach
            TreeSet<Peptide> peptideSet = new TreeSet<>(Collections.reverseOrder());
            Map<String, TreeSet<Peptide>> modSequences = new TreeMap<>();
            for (Peptide peptide : searchObj.getPTMOnlyResult()) {
                Peptide0 peptide0 = peptide0Map.get(peptide.getPTMFreeSeq());
                PeptidePTMPattern peptidePTMPattern = inferPTM.tryPTM(expProcessedPL, plMap, precursorMass, peptide.getPTMFreeSeq(), peptide.isDecoy(), peptide.getNormalizedCrossCorr(), peptide0.leftFlank, peptide0.rightFlank, peptide.getGlobalRank(), localMaxMs2Charge, localMS1ToleranceL, localMS1ToleranceR);
                if (!peptidePTMPattern.getPeptideTreeSet().isEmpty()) {
                    Peptide topPeptide = peptidePTMPattern.getPeptideTreeSet().first();
                    if (peptideSet.size() < 5) {
                        peptideSet.add(topPeptide);
                    } else if (topPeptide.getScore() > peptideSet.last().getScore()) {
                        peptideSet.pollLast();
                        peptideSet.add(topPeptide);
                    }
                    // record scores with different PTM patterns for calculating PTM delta score.
                    modSequences.put(topPeptide.getPTMFreeSeq(), peptidePTMPattern.getPeptideTreeSet());
                }
            }

            // Calculate Score for PTM free peptide
            for (Peptide peptide : searchObj.getPTMFreeResult()) {
                CalScore.calScore(peptide, expProcessedPL, precursorCharge, massToolObj, peptideSet, null);
            }

            if (!peptideSet.isEmpty()) {
                Peptide topPeptide = peptideSet.first();
                TreeSet<Peptide> ptmPatterns = null;
                if (topPeptide.hasVarPTM()) {
                    ptmPatterns = modSequences.get(topPeptide.getPTMFreeSeq());
                }
                new CalSubscores(topPeptide, ms2Tolerance, plMap, precursorCharge, ptmPatterns, binomial);

                Connection sqlConnection = DriverManager.getConnection(sqlPath);
                Statement sqlStatement = sqlConnection.createStatement();
                ResultSet sqlResultSet = sqlStatement.executeQuery(String.format(Locale.US, "SELECT scanNum, scanId, precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, isDecoy, score FROM spectraTable WHERE scanId='%s'", scanId));
                if (sqlResultSet.next()) {
                    int scanNum = sqlResultSet.getInt("scanNum");
                    String scanId = sqlResultSet.getString("scanId");
                    int precursorCharge = sqlResultSet.getInt("precursorCharge");
                    double precursorMass = sqlResultSet.getDouble("precursorMass");
                    String mgfTitle = sqlResultSet.getString("mgfTitle");
                    int isotopeCorrectionNum = sqlResultSet.getInt("isotopeCorrectionNum");
                    double ms1PearsonCorrelationCoefficient = sqlResultSet.getDouble("ms1PearsonCorrelationCoefficient");

                    double deltaLC = 0;
                    if (topPeptide.getScore() > 0) {
                        deltaLC = (topPeptide.getScore() - peptideSet.last().getScore()) / topPeptide.getScore();
                    }
                    double deltaC = 0;
                    if (topPeptide.getScore() > 0) {
                        if (peptideSet.size() > 1) {
                            Iterator<Peptide> temp = peptideSet.iterator();
                            temp.next();
                            deltaC = (topPeptide.getScore() - temp.next().getScore()) / topPeptide.getScore();
                        }
                    }

                    String otherPtmPatterns = "-";
                    if (ptmPatterns != null) {
                        List<String> tempList = new LinkedList<>();
                        Iterator<Peptide> ptmPatternsIterator = ptmPatterns.iterator();
                        ptmPatternsIterator.next();
                        while (ptmPatternsIterator.hasNext()) {
                            Peptide temp = ptmPatternsIterator.next();
                            tempList.add(String.format(Locale.US, "%s-%.4f", temp.getPtmContainingSeq(buildIndexObj.returnFixModMap()), temp.getScore()));
                        }
                        otherPtmPatterns = String.join(";", tempList);
                    }

                    // sqlStatement.executeUpdate(String.format(Locale.US, "REPLACE INTO spectraTable (scanNum, scanId, precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, labelling, peptide, theoMass, isDecoy, globalRank, normalizedCorrelationCoefficient, score, deltaLC, deltaC, matchedPeakNum, ionFrac, matchedHighestIntensityFrac, explainedAaFrac, otherPtmPatterns, aScore) VALUES (%d, '%s', %d, %f, '%s', %d, %f, '%s', '%s', %f, %d, %d, %f, %f, %f, %f, %d, %f, %f, %f, '%s', '%s')", scanNum, scanId, precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, buildIndexObj.getLabelling(), topPeptide.getPtmContainingSeq(buildIndexObj.returnFixModMap()), topPeptide.getTheoMass(), topPeptide.isDecoy() ? 1 : 0, topPeptide.getGlobalRank(), topPeptide.getNormalizedCrossCorr(), topPeptide.getScore(), deltaLC, deltaC, topPeptide.getMatchedPeakNum(), topPeptide.getIonFrac(), topPeptide.getMatchedHighestIntensityFrac(), topPeptide.getExplainedAaFrac(), otherPtmPatterns, topPeptide.getaScore()));
                    Entry entry = new Entry(scanNum, scanId, precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, buildIndexObj.getLabelling(), topPeptide.getPtmContainingSeq(buildIndexObj.returnFixModMap()), topPeptide.getTheoMass(), topPeptide.isDecoy() ? 1 : 0, topPeptide.getGlobalRank(), topPeptide.getNormalizedCrossCorr(), topPeptide.getScore(), deltaLC, deltaC, topPeptide.getMatchedPeakNum(), topPeptide.getIonFrac(), topPeptide.getMatchedHighestIntensityFrac(), topPeptide.getExplainedAaFrac(), otherPtmPatterns, topPeptide.getaScore());

                    sqlResultSet.close();
                    sqlStatement.close();
                    return entry;
                } else {
                    throw new NullPointerException(String.format(Locale.US, "There is no record %s in the spectraTable.", scanId));
                }
            } else {
                return null;
            }
        } else {
            return null;
        }
    }


    public class Entry {

        public final int scanNum;
        public final String scanId;
        public final int precursorCharge;
        public final double precursorMass;
        public final String mgfTitle;
        public final int isotopeCorrectionNum;
        public final double ms1PearsonCorrelationCoefficient;
        public final String labelling;
        public final String peptide;
        public final double theoMass;
        public final int isDecoy;
        public final int globalRank;
        public final double normalizedCorrelationCoefficient;
        public final double score;
        public final double deltaLC;
        public final double deltaC;
        public final int matchedPeakNum;
        public final double ionFrac;
        public final double matchedHighestIntensityFrac;
        public final double explainedAaFrac;
        public final String otherPtmPatterns;
        public final String aScore;

        public Entry(int scanNum, String scanId, int precursorCharge, double precursorMass, String mgfTitle, int isotopeCorrectionNum, double ms1PearsonCorrelationCoefficient, String labelling, String peptide, double theoMass, int isDecoy, int globalRank, double normalizedCorrelationCoefficient, double score, double deltaLC, double deltaC, int matchedPeakNum, double ionFrac, double matchedHighestIntensityFrac, double explainedAaFrac, String otherPtmPatterns, String aScore) {
            this.scanNum = scanNum;
            this.scanId = scanId;
            this.precursorCharge = precursorCharge;
            this.precursorMass = precursorMass;
            this.mgfTitle = mgfTitle;
            this.isotopeCorrectionNum = isotopeCorrectionNum;
            this.ms1PearsonCorrelationCoefficient = ms1PearsonCorrelationCoefficient;
            this.labelling = labelling;
            this.peptide = peptide;
            this.theoMass = theoMass;
            this.isDecoy = isDecoy;
            this.globalRank = globalRank;
            this.normalizedCorrelationCoefficient = normalizedCorrelationCoefficient;
            this.score = score;
            this.deltaLC = deltaLC;
            this.deltaC = deltaC;
            this.matchedPeakNum = matchedPeakNum;
            this.ionFrac = ionFrac;
            this.matchedHighestIntensityFrac = matchedHighestIntensityFrac;
            this.explainedAaFrac = explainedAaFrac;
            this.otherPtmPatterns = otherPtmPatterns;
            this.aScore = aScore;
        }
    }
}
