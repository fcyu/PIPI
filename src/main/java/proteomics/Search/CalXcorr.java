package proteomics.Search;


import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;

import java.util.*;

public class CalXcorr {

    private static final Logger logger = LoggerFactory.getLogger(CalXcorr.class);

    public static void calXcorr(Peptide peptide, SparseVector expXcorrPl, FinalResultEntry psm, MassTool massToolObj, Map<String, LinkedList<PeptideScore>> modSequences) {
        double xcorr = massToolObj.buildVector(peptide.getIonMatrix(), psm.getCharge()).fastDot(expXcorrPl) * 0.25; // scaling the xcorr to original SEQUEST type.
        if (xcorr > 0) {
            if (psm.noScore() || (xcorr > psm.getScore()) || ((xcorr == psm.getScore()) && psm.isDecoy() && (!peptide.isDecoy()))) {
                psm.setPeptide(peptide);
                psm.setGlobalSearchRank(peptide.getGlobalRank());
                psm.setNormalizedCrossXcorr(peptide.getNormalizedCrossCorr());
            }
            psm.addScore(xcorr);

            if (peptide.hasVarPTM()) {
                // record scores with different PTM patterns for calculating PTM delta score.
                if (modSequences.containsKey(peptide.getPTMFreeSeq())) {
                    LinkedList<PeptideScore> temp = modSequences.get(peptide.getPTMFreeSeq());
                    if (temp.size() < 2) {
                        temp.add(new PeptideScore(xcorr, peptide));
                        temp.sort(Collections.reverseOrder());
                    } else if (xcorr > temp.peekLast().score) {
                        temp.pollLast();
                        temp.add(new PeptideScore(xcorr, peptide));
                        temp.sort(Collections.reverseOrder());
                    }
                } else {
                    LinkedList<PeptideScore> temp = new LinkedList<>();
                    temp.add(new PeptideScore(xcorr, peptide));
                    modSequences.put(peptide.getPTMFreeSeq(), temp);
                }
            }
        }
    }
}
