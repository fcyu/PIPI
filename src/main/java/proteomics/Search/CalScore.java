package proteomics.Search;

import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Types.*;
import proteomics.PTM.InferPTM;
import proteomics.Types.*;

import java.util.*;

public class CalScore {

    public static void calScoreForPtmFreePeptide(Peptide peptide, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, int precursorCharge, int localMaxMs2Charge, double ms2Tolerance, MassTool massToolObj, TreeSet<Peptide> peptideSet, Map<String, TreeSet<Peptide>> modSequences) {
        double score = massToolObj.buildVectorAndCalXCorr(peptide.getIonMatrix(), precursorCharge, expProcessedPL); // scaling the xcorr to original SEQUEST type.
        if (score > 0) {
            peptide.setScore(score);
            peptide.setMatchedPeakNum(InferPTM.getMatchedPeakNum(plMap, localMaxMs2Charge, peptide.getIonMatrix(), ms2Tolerance));
            if (peptideSet.size() < 5) {
                peptideSet.add(peptide);
            } else if (peptide.getScore() > peptideSet.last().getScore()) {
                peptideSet.pollLast();
                peptideSet.add(peptide);
            }

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
