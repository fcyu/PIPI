package proteomics.Search;


import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;

import java.util.*;

public class CalXcorr {

    private static final Logger logger = LoggerFactory.getLogger(CalXcorr.class);

    public static void calXcorr(Peptide peptide, SparseVector expXcorrPl, FinalResultEntry psm, MassTool massToolObj, Map<String, TreeSet<PeptideScore>> modSequences) {
        double xcorr = massToolObj.buildVector(peptide.getIonMatrix(), psm.getCharge()).fastDot(expXcorrPl) * 0.25; // scaling the xcorr to original SEQUEST type.
        if (xcorr > 0) {
            if (psm.noScore() || (xcorr > psm.getScore()) || ((xcorr == psm.getScore()) && (peptide.getVarPTMNum() < psm.getPeptide().getVarPTMNum()))) {
                psm.setPeptide(peptide);
                psm.setGlobalSearchRank(peptide.getGlobalRank());
                psm.setNormalizedCrossXcorr(peptide.getNormalizedCrossCorr());
            }
            PeptideScore peptideScore = new PeptideScore(xcorr, peptide);
            psm.addScore(peptideScore);

            if (peptide.hasVarPTM()) {
                // record scores with different PTM patterns for calculating PTM delta score.
                if (modSequences.containsKey(peptide.getPTMFreeSeq())) {
                    TreeSet<PeptideScore> temp = modSequences.get(peptide.getPTMFreeSeq());
                    if (temp.size() < 5) {
                        temp.add(peptideScore);
                    } else if (xcorr > temp.last().score) {
                        temp.pollLast();
                        temp.add(peptideScore);
                    }
                } else {
                    TreeSet<PeptideScore> temp = new TreeSet<>(Collections.reverseOrder());
                    temp.add(peptideScore);
                    modSequences.put(peptide.getPTMFreeSeq(), temp);
                }
            }
        }
    }
}
