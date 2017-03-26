package proteomics.Search;


import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Index.BuildIndex;
import proteomics.Spectrum.PreSpectrum;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;

import java.util.*;

public class CalXcorr {

    private static final Logger logger = LoggerFactory.getLogger(CalXcorr.class);
    private static final int minDecoyNum = 3000;
    private static final float evalueTolerance1 = 3;
    private static final float evalueTolerance2 = 20;

    private FinalResultEntry psm;
    private final MassTool massToolObj;

    public CalXcorr(List<Peptide> candidateList, SpectrumEntry spectrum, MassTool massToolObj, BuildIndex buildIndexObj) {
        this.massToolObj = massToolObj;
        PreSpectrum preSpectrumObj = new PreSpectrum(massToolObj);

        // prepare the XCORR vector
        SparseVector expXcorrPl = preSpectrumObj.prepareXcorr(spectrum.plMap);

        // calculate Xcorr
        psm = new FinalResultEntry(spectrum.scanNum, spectrum.precursorCharge, spectrum.precursorMz);
        for (Peptide peptide : candidateList) {
            SparseBooleanVector theoIonVector = massToolObj.buildVector(peptide.getIonMatrix(), spectrum.precursorCharge);
            double xcorr = theoIonVector.dot(expXcorrPl) * 0.25; // scaling the xcorr to original SEQUEST type.
            if (xcorr > 0) {
                psm.addToScoreHistogram(xcorr);
                psm.addScoredPeptide(peptide.getPTMFreeSeq());
                if (psm.noScore() || (xcorr > psm.getScore()) || ((xcorr == psm.getScore()) && psm.isDecoy() && (!peptide.isDecoy()))) {
                    psm.setPeptide(peptide);
                    psm.setGlobalSearchRank(peptide.getGlobalRank());
                    psm.setNormalizedCrossXcorr(peptide.getNormalizedCrossCorr());
                }
                psm.addScore(xcorr);
            }
        }

        if (psm.getPeptide() != null) {
            // generate permutated decoy peptides to fill up the number of total decoy peptides.
            if (psm.getCandidateNum() < minDecoyNum) {
                // generateDecoyScores(psm, expXcorrPl, numCandidateMap.get(scanNum), minDecoyNum - psm.getCandidateNum(), spectrum.precursorCharge);
                generateDecoyScores(psm, expXcorrPl, buildIndexObj.getMassPeptideMap(),  minDecoyNum - psm.getCandidateNum(), spectrum.precursorMass, spectrum.precursorCharge);
            }
        } else {
            psm = null;
        }
    }

    public FinalResultEntry getScoredPSM() {
        return psm;
    }

    private void generateDecoyScores(FinalResultEntry psm, SparseVector expXcorrPl, TreeMap<Float, Set<String>> massPeptideMap, int gapNum, float precursorMass, int precursorCharge) {
        // 3 Da tolerance first
        NavigableMap<Float, Set<String>> subMap = massPeptideMap.subMap(precursorMass - evalueTolerance1, true, precursorMass + evalueTolerance1, true);
        gapNum = accumulate(psm, expXcorrPl, subMap, gapNum, precursorCharge, false);

        // 20 Da tolerance second
        if (gapNum > 0) {
            subMap = new TreeMap<>(massPeptideMap.subMap(precursorMass - evalueTolerance2, true, precursorMass - evalueTolerance1, false));
            subMap.putAll(massPeptideMap.subMap(precursorMass + evalueTolerance1, false, precursorMass + evalueTolerance2, true));
            if (!subMap.isEmpty()) {
                accumulate(psm, expXcorrPl, subMap, gapNum, precursorCharge, true);
            }
        }
    }

    private int accumulate(FinalResultEntry psm, SparseVector expXcorrPl, NavigableMap<Float, Set<String>> subMap, int gapNum, int precursorCharge, boolean notCheckScored) {
        for (Set<String> peptideSet : subMap.values()) {
            for (String peptide : peptideSet) {
                if (notCheckScored || !psm.scored(peptide)) {
                    SparseBooleanVector theoIonVector = massToolObj.buildVector(massToolObj.buildIonArray(peptide, precursorCharge), precursorCharge);
                    double dotProduct = theoIonVector.dot(expXcorrPl) * 0.25; // doesn't to make it larger than 0. for histogram accumulation.
                    if (dotProduct > 0) {
                        psm.addToScoreHistogram(dotProduct);
                        --gapNum;
                        if (gapNum == 0) {
                            return gapNum;
                        }
                    }
                }
            }
        }
        return gapNum;
    }
}
