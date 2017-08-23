package proteomics.Search;


import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;

import java.util.*;

public class CalScore {

    private static final Logger logger = LoggerFactory.getLogger(CalScore.class);

    public static void calScore(Peptide peptide, SparseVector expProcessedPL, FinalResultEntry psm, MassTool massToolObj, Map<String, TreeSet<PeptideScore>> modSequences) {
        double score = massToolObj.buildVector(peptide.getIonMatrix(), psm.getCharge()).fastDot(expProcessedPL) * 0.25; // scaling the xcorr to original SEQUEST type.
        if (score > 0) {
            PeptideScore peptideScore = new PeptideScore(score, peptide);
            psm.addScore(peptideScore);

            if (peptide.hasVarPTM()) {
                // record scores with different PTM patterns for calculating PTM delta score.
                if (modSequences.containsKey(peptide.getPTMFreeSeq())) {
                    TreeSet<PeptideScore> temp = modSequences.get(peptide.getPTMFreeSeq());
                    if (temp.size() < 5) {
                        temp.add(peptideScore);
                    } else if (score > temp.last().score) {
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
