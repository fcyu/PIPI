package proteomics.Search;


import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;

import java.util.*;

public class CalScore {

    private static final Logger logger = LoggerFactory.getLogger(CalScore.class);

    public static void calScore(Peptide peptide, SparseVector expProcessedPL, FinalResultEntry psm, MassTool massToolObj, Map<String, TreeSet<Peptide>> modSequences) {
        double score = massToolObj.buildVector(peptide.getIonMatrix(), psm.getCharge()).fastDot(expProcessedPL) * 0.25; // scaling the xcorr to original SEQUEST type.
        if (score > 0) {
            peptide.setScore(score);
            psm.addScore(peptide);

            if (peptide.hasVarPTM()) {
                // record scores with different PTM patterns for calculating PTM delta score.
                if (modSequences.containsKey(peptide.getPTMFreeSeq())) {
                    TreeSet<Peptide> temp = modSequences.get(peptide.getPTMFreeSeq());
                    if (temp.size() < 5) {
                        temp.add(peptide);
                    } else if (score > temp.last().getScore()) {
                        temp.pollLast();
                        temp.add(peptide);
                    }
                } else {
                    TreeSet<Peptide> temp = new TreeSet<>(Collections.reverseOrder());
                    temp.add(peptide);
                    modSequences.put(peptide.getPTMFreeSeq(), temp);
                }
            }
        }
    }
}
