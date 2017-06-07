package proteomics;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Index.BuildIndex;
import proteomics.PTM.FindPTM;
import proteomics.PTM.GeneratePtmCandidates;
import proteomics.Search.CalSubscores;
import proteomics.Search.CalXcorr;
import proteomics.Search.Search;
import proteomics.Segment.InferenceSegment;
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
    private final Map<String, Peptide0> peptideCodeMap;
    private final float ms1Tolerance;
    private final int ms1ToleranceUnit;
    private final float ms2Tolerance;
    private final float minPtmMass;
    private final float maxPtmMass;
    private final int maxMs2Charge;

    public PIPIWrap(BuildIndex buildIndexObj, MassTool massToolObj, InferenceSegment inference3SegmentObj, SpectrumEntry spectrumEntry, Map<String, Peptide0> peptideCodeMap, float ms1Tolerance, int ms1ToleranceUnit, float ms2Tolerance, float minPtmMass, float maxPtmMass, int maxMs2Charge) {
        this.buildIndexObj = buildIndexObj;
        this.massToolObj = massToolObj;
        this.inference3SegmentObj = inference3SegmentObj;
        this.spectrumEntry = spectrumEntry;
        this.peptideCodeMap = peptideCodeMap;
        this.ms1Tolerance = ms1Tolerance;
        this.ms1ToleranceUnit = ms1ToleranceUnit;
        this.ms2Tolerance = ms2Tolerance;
        this.minPtmMass = minPtmMass;
        this.maxPtmMass = maxPtmMass;
        this.maxMs2Charge = maxMs2Charge;
    }

    @Override
    public FinalResultEntry call() {
        // Coding
        List<ThreeExpAA> expAaLists = inference3SegmentObj.inferSegmentLocationFromSpectrum(spectrumEntry);
        if (!expAaLists.isEmpty()) {
            SparseVector scanCode = inference3SegmentObj.generateSegmentIntensityVector(expAaLists);

            // Begin search.
            Search searchObj = new Search(buildIndexObj, spectrumEntry, scanCode, peptideCodeMap, buildIndexObj.getMassPeptideMap(), massToolObj, ms1Tolerance, ms1ToleranceUnit, minPtmMass, maxPtmMass, maxMs2Charge);

            // Infer PTMs based on DP
            FindPTM findPtmObj = new FindPTM(searchObj.getPTMOnlyResult(), spectrumEntry, expAaLists, inference3SegmentObj.getModifiedAAMassMap(), inference3SegmentObj.getPepNTermPossibleMod(), inference3SegmentObj.getPepCTermPossibleMod(), inference3SegmentObj.getProNTermPossibleMod(), inference3SegmentObj.getProCTermPossibleMod(), minPtmMass, maxPtmMass, ms1Tolerance, ms1ToleranceUnit, ms2Tolerance);

            GeneratePtmCandidates generatePtmCandidates = new GeneratePtmCandidates(inference3SegmentObj, spectrumEntry, massToolObj, buildIndexObj, ms2Tolerance, maxMs2Charge);
            // Generate all candidates based on inferred PTMs and known PTMs.
            Set<Peptide> finalCandidates = generatePtmCandidates.generateAllPtmCandidates(generatePtmCandidates.eliminateMissedCleavageCausedPtms(searchObj.getPTMFreeResult(), findPtmObj.getPeptidesWithPTMs()));

            if (!finalCandidates.isEmpty()) {
                CalXcorr calXcorrObj = new CalXcorr(finalCandidates, spectrumEntry, massToolObj);
                FinalResultEntry scoredPsm = calXcorrObj.getScoredPSM();
                if (scoredPsm != null) {
                    new CalSubscores(scoredPsm, spectrumEntry, ms2Tolerance);
                    return scoredPsm;
                } else {
                    return null;
                }
            } else {
                return null;
            }
        } else {
            return null;
        }
    }
}
