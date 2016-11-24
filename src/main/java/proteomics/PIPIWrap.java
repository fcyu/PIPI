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

public class PIPIWrap implements Callable<List<FinalResultEntry>> {

    private static final Logger logger = LoggerFactory.getLogger(PIPIWrap.class);

    private final BuildIndex buildIndexObj;
    private final MassTool massToolObj;
    private final InferenceSegment inference3SegmentObj;
    private final Map<Integer, SpectrumEntry> numSpectrumMap;
    private final NavigableMap<Float, List<Integer>> subMassNumMap;
    private final Map<String, TreeSet<Integer>> siteMass1000Map;
    private final float ms1Tolerance;
    private final int ms1ToleranceUnit;
    private final float ms2Tolerance;
    private final float minPtmMass;
    private final float maxPtmMass;
    private final int maxMs2Charge;
    private final int batchStartIdx;

    public PIPIWrap(BuildIndex buildIndexObj, MassTool massToolObj, InferenceSegment inference3SegmentObj, Map<Integer, SpectrumEntry> numSpectrumMap, NavigableMap<Float, List<Integer>> subMassNumMap, Map<String, TreeSet<Integer>> siteMass1000Map, float ms1Tolerance, int ms1ToleranceUnit, float ms2Tolerance, float minPtmMass, float maxPtmMass, int maxMs2Charge, int batchStartIdx) {
        this.buildIndexObj = buildIndexObj;
        this.massToolObj = massToolObj;
        this.inference3SegmentObj = inference3SegmentObj;
        this.numSpectrumMap = numSpectrumMap;
        this.subMassNumMap = subMassNumMap;
        this.siteMass1000Map = siteMass1000Map;
        this.ms1Tolerance = ms1Tolerance;
        this.ms1ToleranceUnit = ms1ToleranceUnit;
        this.ms2Tolerance = ms2Tolerance;
        this.minPtmMass = minPtmMass;
        this.maxPtmMass = maxPtmMass;
        this.maxMs2Charge = maxMs2Charge;
        this.batchStartIdx = batchStartIdx;
    }

    @Override
    public List<FinalResultEntry> call() {
        try {
            // Inference amino acids.
            Map<Integer, SparseVector> numCodeMap = new HashMap<>();
            Map<Integer, List<ThreeExpAA>> numExp3aaLists = new HashMap<>();
            for (List<Integer> numList : subMassNumMap.values()) {
                for (int scanNum : numList) {
                    SpectrumEntry spectrumEntry = numSpectrumMap.get(scanNum);

                    // Coding
                    List<ThreeExpAA> expAaLists = inference3SegmentObj.inferSegmentLocationFromSpectrum(spectrumEntry);
                    if (!expAaLists.isEmpty()) {
                        numExp3aaLists.put(scanNum, expAaLists);
                        numCodeMap.put(scanNum, inference3SegmentObj.generateSegmentIntensityVector(expAaLists));
                    }
                }
            }

            // Begin search.
            Search searchObj = new Search(buildIndexObj, numCodeMap, inference3SegmentObj, subMassNumMap, ms1Tolerance, ms1ToleranceUnit, minPtmMass, maxPtmMass, maxMs2Charge, batchStartIdx);

            logger.debug("Analyzing PTMs...");
            FindPTM findPtmObj = new FindPTM(searchObj.getPTMOnlyResult(), numSpectrumMap, numExp3aaLists, massToolObj, siteMass1000Map, minPtmMass, maxPtmMass, ms1Tolerance, ms1ToleranceUnit, ms2Tolerance, batchStartIdx);
            Map<Integer, List<Peptide>> ptmOnlyTemp = findPtmObj.getPeptidesWithPTMs();

            logger.debug("Calculating final score...");
            Map<Integer, List<Peptide>> numCandidateMap = new HashMap<>(searchObj.getPTMFreeResult());
            for (int scanNum : ptmOnlyTemp.keySet()) {
                if (numCandidateMap.containsKey(scanNum)) {
                    numCandidateMap.get(scanNum).addAll(ptmOnlyTemp.get(scanNum));
                } else {
                    numCandidateMap.put(scanNum, ptmOnlyTemp.get(scanNum));
                }
            }
            CalXcorr calXcorrObj = new CalXcorr(numCandidateMap, numSpectrumMap, massToolObj, buildIndexObj);
            List<FinalResultEntry> subScoredPsms = calXcorrObj.getScoredPSMs();

            new CalSubscores(subScoredPsms, numSpectrumMap, ms2Tolerance);

            logger.debug("Estimating E-value for each PSM...");
            for (FinalResultEntry psm : subScoredPsms) {
                new CalEValue(psm);
            }

            System.gc();
            return subScoredPsms;
        } catch (Exception ex) {
            logger.error(ex.getMessage());
            ex.printStackTrace();
            System.exit(1);
        }
        return new LinkedList<>();
    }
}
