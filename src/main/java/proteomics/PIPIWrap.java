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

import java.util.*;
import java.util.concurrent.Callable;

public class PIPIWrap implements Callable<FinalResultEntry> {

    private static final Logger logger = LoggerFactory.getLogger(PIPIWrap.class);

    private final BuildIndex buildIndexObj;
    private final MassTool massToolObj;
    private final InferenceSegment inference3SegmentObj;
    private final SpectrumEntry spectrumEntry;
    private final float ms1Tolerance;
    private final int ms1ToleranceUnit;
    private final float ms2Tolerance;
    private final float minPtmMass;
    private final float maxPtmMass;
    private final int maxMs2Charge;
    private final InferPTM inferPTM;


    public PIPIWrap(BuildIndex buildIndexObj, MassTool massToolObj, SpectrumEntry spectrumEntry, float ms1Tolerance, int ms1ToleranceUnit, float ms2Tolerance, float minPtmMass, float maxPtmMass, int maxMs2Charge) {
        this.buildIndexObj = buildIndexObj;
        this.massToolObj = massToolObj;
        this.spectrumEntry = spectrumEntry;
        this.ms1Tolerance = ms1Tolerance;
        this.ms1ToleranceUnit = ms1ToleranceUnit;
        this.ms2Tolerance = ms2Tolerance;
        this.minPtmMass = minPtmMass;
        this.maxPtmMass = maxPtmMass;
        this.maxMs2Charge = maxMs2Charge;
        inference3SegmentObj = buildIndexObj.getInference3SegmentObj();
        inferPTM = new InferPTM(massToolObj, maxMs2Charge, buildIndexObj.returnFixModMap(), inference3SegmentObj.getVarModParamSet(), minPtmMass, maxPtmMass, ms2Tolerance);
    }

    @Override
    public FinalResultEntry call() {
        // Coding
        List<ThreeExpAA> expAaLists = inference3SegmentObj.inferSegmentLocationFromSpectrum(spectrumEntry);
        if (!expAaLists.isEmpty()) {
            SparseVector scanCode = inference3SegmentObj.generateSegmentIntensityVector(expAaLists);

            // Begin search.
            Search searchObj = new Search(buildIndexObj, spectrumEntry, scanCode, massToolObj, ms1Tolerance, ms1ToleranceUnit, minPtmMass, maxPtmMass, maxMs2Charge);

            // prepare the spectrum
            PreSpectrum preSpectrumObj = new PreSpectrum(massToolObj);
            SparseVector expProcessedPL;
            if (PIPI.useXcorr) {
                expProcessedPL = preSpectrumObj.prepareXcorr(spectrumEntry.plMap, false);
            } else {
                expProcessedPL = preSpectrumObj.prepareDigitizedPL(spectrumEntry.plMap, false);
            }

            FinalResultEntry psm = new FinalResultEntry(spectrumEntry.scanNum, spectrumEntry.precursorCharge, spectrumEntry.precursorMz, spectrumEntry.mgfTitle);

            float precursorMass = spectrumEntry.precursorMass;
            float localMS1ToleranceL = -1 * ms1Tolerance;
            float localMS1ToleranceR = ms1Tolerance;
            if (ms1ToleranceUnit == 1) {
                localMS1ToleranceL = (precursorMass / (1 + ms1Tolerance * 1e-6f)) - precursorMass;
                localMS1ToleranceR = (precursorMass / (1 - ms1Tolerance * 1e-6f)) - precursorMass;
            }

            // infer PTM using the new approach
            Map<String, TreeSet<Peptide>> modSequences = new TreeMap<>();
            for (Peptide peptide : searchObj.getPTMOnlyResult()) {
                PeptidePTMPattern peptidePTMPattern = inferPTM.tryPTM(expProcessedPL, spectrumEntry.plMap, precursorMass, peptide.getPTMFreeSeq(), peptide.isDecoy(), peptide.getNormalizedCrossCorr(), peptide.getLeftFlank(), peptide.getRightFlank(), peptide.getGlobalRank(), maxMs2Charge, localMS1ToleranceL, localMS1ToleranceR);
                if (!peptidePTMPattern.getPeptideTreeSet().isEmpty()) {
                    Peptide topPeptide = peptidePTMPattern.getPeptideTreeSet().first();
                    psm.addScore(topPeptide);
                    // record scores with different PTM patterns for calculating PTM delta score.
                    modSequences.put(topPeptide.getPTMFreeSeq(), peptidePTMPattern.getPeptideTreeSet());
                }
            }

            // Calculate Score for PTM free peptide
            for (Peptide peptide : searchObj.getPTMFreeResult()) {
                CalScore.calScore(peptide, expProcessedPL, psm, massToolObj, null);
            }

            if (psm.hasHit()) {
                psm.setPtmPatterns(modSequences);
                for (Peptide peptide : psm.getPeptideSet()) {
                    new CalSubscores(peptide, spectrumEntry, ms2Tolerance);
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
