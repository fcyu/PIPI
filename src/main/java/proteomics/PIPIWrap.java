package proteomics;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Index.BuildIndex;
import proteomics.PTM.InferPTM;
import proteomics.Search.CalSubscores;
import proteomics.Search.CalScore;
import proteomics.Search.Search;
import proteomics.Segment.InferenceSegment;
import proteomics.Spectrum.PreSpectrum;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.JMzReaderException;

import java.sql.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantLock;

public class PIPIWrap implements Callable<Boolean> {

    private static final Logger logger = LoggerFactory.getLogger(PIPIWrap.class);

    private final BuildIndex buildIndexObj;
    private final MassTool massToolObj;
    private final float ms1Tolerance;
    private final int ms1ToleranceUnit;
    private final float ms2Tolerance;
    private final float minPtmMass;
    private final float maxPtmMass;
    private final int maxMs2Charge;
    private final Map<String, Peptide0> peptide0Map;
    private final JMzReader spectraParser;
    private final float minClear;
    private final float maxClear;
    private final ReentrantLock lock;
    private final String scanId;
    private final int precursorCharge;
    private final float precursorMass;
    private final InferPTM inferPTM;
    private final PreSpectrum preSpectrum;
    private final Connection sqlConnection;


    public PIPIWrap(BuildIndex buildIndexObj, MassTool massToolObj, float ms1Tolerance, int ms1ToleranceUnit, float ms2Tolerance, float minPtmMass, float maxPtmMass, int maxMs2Charge, JMzReader spectraParser, float minClear, float maxClear, ReentrantLock lock, String scanId, int precursorCharge, float precursorMass, InferPTM inferPTM, PreSpectrum preSpectrum, Connection sqlConnection) {
        this.buildIndexObj = buildIndexObj;
        this.massToolObj = massToolObj;
        this.ms1Tolerance = ms1Tolerance;
        this.ms1ToleranceUnit = ms1ToleranceUnit;
        this.ms2Tolerance = ms2Tolerance;
        this.minPtmMass = minPtmMass;
        this.maxPtmMass = maxPtmMass;
        this.maxMs2Charge = maxMs2Charge;
        this.spectraParser = spectraParser;
        this.minClear = minClear;
        this.maxClear = maxClear;
        this.lock = lock;
        this.scanId = scanId;
        this.precursorCharge = precursorCharge;
        this.precursorMass = precursorMass;
        this.inferPTM = inferPTM;
        this.preSpectrum = preSpectrum;
        this.sqlConnection = sqlConnection;
        peptide0Map = buildIndexObj.getPeptide0Map();
    }

    @Override
    public Boolean call() throws SQLException {
        Map<Double, Double> rawPLMap = null;
        try {
            lock.lock();
            // Reading peak list.
            rawPLMap = spectraParser.getSpectrumById(scanId).getPeakList();
        } catch (JMzReaderException ex) {
            lock.unlock();
            ex.printStackTrace();
            logger.error(ex.toString());
            System.exit(1);
        } finally {
            lock.unlock();
        }

        // preprocess peak list
        TreeMap<Float, Float> plMap = preSpectrum.preSpectrum(rawPLMap, precursorMass, precursorCharge, ms2Tolerance, minClear, maxClear);

        // Coding
        InferenceSegment inference3SegmentObj = buildIndexObj.getInference3SegmentObj();
        List<ThreeExpAA> expAaLists = inference3SegmentObj.inferSegmentLocationFromSpectrum(precursorMass, plMap);
        if (!expAaLists.isEmpty()) {
            SparseVector scanCode = inference3SegmentObj.generateSegmentIntensityVector(expAaLists);

            // Begin search.
            Search searchObj = new Search(buildIndexObj, precursorMass, scanCode, massToolObj, ms1Tolerance, ms1ToleranceUnit, minPtmMass, maxPtmMass, maxMs2Charge);

            // prepare the spectrum
            SparseVector expProcessedPL;
            if (PIPI.useXcorr) {
                expProcessedPL = preSpectrum.prepareXcorr(plMap, false);
            } else {
                expProcessedPL = preSpectrum.prepareDigitizedPL(plMap, false);
            }

            float localMS1ToleranceL = -1 * ms1Tolerance;
            float localMS1ToleranceR = ms1Tolerance;
            if (ms1ToleranceUnit == 1) {
                localMS1ToleranceL = (precursorMass / (1 + ms1Tolerance * 1e-6f)) - precursorMass;
                localMS1ToleranceR = (precursorMass / (1 - ms1Tolerance * 1e-6f)) - precursorMass;
            }

            // infer PTM using the new approach
            TreeSet<Peptide> peptideSet = new TreeSet<>(Collections.reverseOrder());
            Map<String, TreeSet<Peptide>> modSequences = new TreeMap<>();
            for (Peptide peptide : searchObj.getPTMOnlyResult()) {
                Peptide0 peptide0 = peptide0Map.get(peptide.getPTMFreeSeq());
                PeptidePTMPattern peptidePTMPattern = inferPTM.tryPTM(expProcessedPL, plMap, precursorMass, peptide.getPTMFreeSeq(), peptide.isDecoy(), peptide.getNormalizedCrossCorr(), peptide0.leftFlank, peptide0.rightFlank, peptide.getGlobalRank(), maxMs2Charge, localMS1ToleranceL, localMS1ToleranceR);
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
                new CalSubscores(topPeptide, ms2Tolerance, plMap, precursorCharge);

                Statement sqlStatement = sqlConnection.createStatement();
                ResultSet sqlResultSet = sqlStatement.executeQuery(String.format(Locale.US, "SELECT scanNum, scanId, precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, isDecoy, score FROM spectraTable WHERE scanId='%s'", scanId)); // todo: check
                if (sqlResultSet.next()) {
                    boolean needUpdate = false;
                    int isDecoyOld = sqlResultSet.getInt("isDecoy");
                    if (!sqlResultSet.wasNull()) {
                        double scoreOld = sqlResultSet.getDouble("score");
                        if (topPeptide.getScore() > scoreOld || (topPeptide.getScore() == scoreOld && isDecoyOld == 1 && !topPeptide.isDecoy())) {
                            needUpdate = true;
                        }
                    } else {
                        needUpdate = true;
                    }
                    if (needUpdate) {
                        int scanNum = sqlResultSet.getInt("scanNum");
                        String scanId = sqlResultSet.getString("scanId");
                        int precursorCharge = sqlResultSet.getInt("precursorCharge");
                        float precursorMass = sqlResultSet.getFloat("precursorMass");
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
                        String ptmDeltaScore = "-";
                        if (ptmPatterns != null) {
                            List<String> tempList = new LinkedList<>();
                            Iterator<Peptide> ptmPatternsIterator = ptmPatterns.iterator();
                            ptmPatternsIterator.next();
                            while (ptmPatternsIterator.hasNext()) {
                                Peptide temp = ptmPatternsIterator.next();
                                tempList.add(String.format(Locale.US, "%s-%.4f;", temp.getPtmContainingSeq(buildIndexObj.returnFixModMap()), temp.getScore()));
                            }
                            otherPtmPatterns = String.join(";", tempList);
                            if (ptmPatterns.size() > 1) {
                                Iterator<Peptide> temp = ptmPatterns.iterator();
                                ptmDeltaScore = String.valueOf(temp.next().getScore() - temp.next().getScore());
                            } else {
                                ptmDeltaScore = String.valueOf(ptmPatterns.first().getScore());
                            }
                        }

                        sqlStatement.executeUpdate(String.format(Locale.US, "DELETE FROM spectraTable WHERE scanId=%s", scanId));
                        sqlStatement.executeUpdate(String.format(Locale.US, "INSERT INTO spectraTable (scanNum, scanId, precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, labelling, peptide, theoMass, isDecoy, globalRank, normalizedCorrelationCoefficient, score, deltaLC, deltaC, matchedPeakNum, ionFrac, matchedHighestIntensityFrac, explainedAaFrac, ptmSupportingPeakFrac, otherPtmPatterns, ptmDeltaScore) VALUES (%d, '%s', %d, %f, '%s', %d, %f, '%s', '%s', %f, %d, %d, %f, %f, %f, %f, %d, %f, %f, %f, %f, '%s', '%s')", scanNum, scanId, precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, buildIndexObj.getLabeling(), topPeptide.getPtmContainingSeq(buildIndexObj.returnFixModMap()), topPeptide.getTheoMass(), topPeptide.isDecoy() ? 1 : 0, topPeptide.getGlobalRank(), topPeptide.getNormalizedCrossCorr(), topPeptide.getScore(), deltaLC, deltaC, topPeptide.getMatchedPeakNum(), topPeptide.getIonFrac(), topPeptide.getMatchedHighestIntensityFrac(), topPeptide.getExplainedAaFrac(), topPeptide.getPtmSupportingPeakFrac(), otherPtmPatterns, ptmDeltaScore));
                    }
                } else {
                    throw new NullPointerException(String.format(Locale.US, "There is no record %s in the spectraTable.", scanId));
                }

                sqlResultSet.close();
                sqlStatement.close();
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    }
}
