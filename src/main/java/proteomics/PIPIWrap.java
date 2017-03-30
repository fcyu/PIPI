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

            FindPTM findPtmObj = new FindPTM(searchObj.getPTMOnlyResult(), spectrumEntry, expAaLists, massToolObj, inference3SegmentObj.getModifiedAAMassMap(), inference3SegmentObj.getnTermPossibleMod(), inference3SegmentObj.getcTermPossibleMod(), minPtmMass, maxPtmMass, ms1Tolerance, ms1ToleranceUnit, ms2Tolerance);
            List<Peptide> ptmOnlyTemp = findPtmObj.getPeptidesWithPTMs();

            List<Peptide> tempList = new LinkedList<>();
            tempList.addAll(searchObj.getPTMFreeResult());
            tempList.addAll(ptmOnlyTemp);
            Peptide[] tempArray = tempList.toArray(new Peptide[tempList.size()]);

            // Additional short missed cleavaged amino acid sequence may be canceled out by negative PTM. In this situation, we eliminate the missed cleavaged one.
            List<Peptide> candidates = new LinkedList<>();
            for (int i = 0; i < tempArray.length; ++i) {
                boolean keep = true;
                String tempStr1 = tempArray[i].getNormalizedPeptideString();
                tempStr1 = tempStr1.substring(1, tempStr1.length() - 1);
                for (int j = 0; j < tempArray.length; ++j) {
                    if (i != j) {
                        String tempStr2 = tempArray[j].getNormalizedPeptideString();
                        tempStr2 = tempStr2.substring(1, tempStr2.length() - 1);
                        if (tempStr1.contains(tempStr2) && (tempArray[i].getVarPTMs() != null)) {
                            Map.Entry<Coordinate, Float> tempEntry = tempArray[i].getVarPTMs().firstEntry();
                            if ((tempEntry.getValue() < 0) && (tempEntry.getKey().y - tempEntry.getKey().x > 1)) {
                                keep = false;
                                break;
                            }
                            tempEntry = tempArray[i].getVarPTMs().lastEntry();
                            if ((tempEntry.getValue() < 0) && (tempEntry.getKey().y - tempEntry.getKey().x > 1)) {
                                keep = false;
                                break;
                            }
                        }
                    }
                }
                if (keep) {
                    candidates.add(tempArray[i]);
                }
            }

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
