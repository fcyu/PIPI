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

import java.io.OutputStream;
import java.io.PrintStream;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantLock;

public class PIPIWrap implements Callable<FinalResultEntry> {

    private static final Logger logger = LoggerFactory.getLogger(PIPIWrap.class);

    private final BuildIndex buildIndexObj;
    private final MassTool massToolObj;
    private final SpectrumEntry spectrumEntry;
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


    public PIPIWrap(BuildIndex buildIndexObj, MassTool massToolObj, float ms1Tolerance, int ms1ToleranceUnit, float ms2Tolerance, float minPtmMass, float maxPtmMass, int maxMs2Charge, JMzReader spectraParser, float minClear, float maxClear, ReentrantLock lock, String scanId, int precursorCharge, float precursorMass, InferPTM inferPTM, PreSpectrum preSpectrum, Connection sqlConnection) {
        this.buildIndexObj = buildIndexObj;
        this.massToolObj = massToolObj;
        this.spectrumEntry = spectrumEntry;
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
        peptide0Map = buildIndexObj.getPeptide0Map();
    }

    @Override
    public FinalResultEntry call() {
        Map<Double, Double> rawPLMap = null;
        PrintStream originalStream = System.out;
        try {
            lock.lock();
            PrintStream nullStream = new PrintStream(new OutputStream() {
                @Override
                public void write(int b) { }
            });
            System.setOut(nullStream);

            // Reading peak list.
            rawPLMap = spectraParser.getSpectrumById(scanId).getPeakList();
        } catch (JMzReaderException ex) {
            ex.printStackTrace();
            logger.error(ex.toString());
            System.exit(1);
        } finally {
            System.setOut(originalStream);
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

            FinalResultEntry psm = new FinalResultEntry(spectrumEntry.scanNum, spectrumEntry.precursorCharge, spectrumEntry.precursorMz, spectrumEntry.mgfTitle, buildIndexObj.getLabeling(), spectrumEntry.isotopeCorrectionNum, spectrumEntry.ms1PearsonCorrelationCoefficient);

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

                if (psm.getTopPeptide().hasVarPTM()) {
                    psm.setPtmPatterns(modSequences.get(psm.getTopPeptide().getPTMFreeSeq()));
            if (!peptideSet.isEmpty()) {
                new CalSubscores(topPeptide, ms2Tolerance, plMap, precursorCharge);
                }
                return psm;
            } else {
                return null;
            }
        } else {
            return null;
        }
    }
}
