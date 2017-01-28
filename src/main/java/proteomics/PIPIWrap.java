package proteomics;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Index.BuildIndex;
import proteomics.PTM.FindPTM;
import proteomics.Search.CalEValue;
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
    private final Map<String, TreeSet<Integer>> siteMass1000Map;
    private final Map<String, Peptide0> peptideCodeMap;
    private final float ms1Tolerance;
    private final int ms1ToleranceUnit;
    private final float ms2Tolerance;
    private final float minPtmMass;
    private final float maxPtmMass;
    private final int maxMs2Charge;

    public PIPIWrap(BuildIndex buildIndexObj, MassTool massToolObj, InferenceSegment inference3SegmentObj, SpectrumEntry spectrumEntry, Map<String, TreeSet<Integer>> siteMass1000Map, Map<String, Peptide0> peptideCodeMap, float ms1Tolerance, int ms1ToleranceUnit, float ms2Tolerance, float minPtmMass, float maxPtmMass, int maxMs2Charge) {
        this.buildIndexObj = buildIndexObj;
        this.massToolObj = massToolObj;
        this.inference3SegmentObj = inference3SegmentObj;
        this.spectrumEntry = spectrumEntry;
        this.siteMass1000Map = siteMass1000Map;
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

            FindPTM findPtmObj = new FindPTM(searchObj.getPTMOnlyResult(), spectrumEntry, expAaLists, massToolObj, siteMass1000Map, minPtmMass, maxPtmMass, ms1Tolerance, ms1ToleranceUnit, ms2Tolerance);
            List<Peptide> ptmOnlyTemp = findPtmObj.getPeptidesWithPTMs();

            List<Peptide> candidates = new LinkedList<>();
            candidates.addAll(searchObj.getPTMFreeResult());
            candidates.addAll(ptmOnlyTemp);

            if (!candidates.isEmpty()) {
                CalXcorr calXcorrObj = new CalXcorr(candidates, spectrumEntry, massToolObj, buildIndexObj);
                FinalResultEntry scoredPsm = calXcorrObj.getScoredPSM();
                if (scoredPsm != null) {
                    new CalSubscores(scoredPsm, spectrumEntry, ms2Tolerance);
                    new CalEValue(scoredPsm);
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
