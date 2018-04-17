package proteomics.Search;


import ProteomicsLibrary.Binomial;
import ProteomicsLibrary.Score;
import proteomics.Spectrum.PreSpectra;
import proteomics.Types.Peptide;

import java.util.*;

public class CalSubscores {

    public CalSubscores(Peptide peptide, double ms2Tolerance, TreeMap<Double, Double> expPl, int precursorCharge, TreeSet<Peptide> ptmPatterns, Binomial binomial) throws Exception {
        peptide.setIonFrac(Score.calIonFraction(peptide.getIonMatrix(), precursorCharge, expPl, ms2Tolerance));
        peptide.setMatchedHighestIntensityFrac(Score.calMatchedHighestIntensityFraction(peptide.getIonMatrix(), precursorCharge, expPl, ms2Tolerance));
        peptide.setExplainedAaFrac(Score.calExplainedAAFraction(peptide.getIonMatrix(), precursorCharge, expPl, ms2Tolerance));

        // calculate A score
        if (peptide.hasVarPTM()) {
            Peptide[] tempArray = ptmPatterns.toArray(new Peptide[0]);
            peptide.setaScore(String.valueOf(Score.calAScore(expPl, PreSpectra.topN, binomial, peptide.getVarPTMs(), peptide.getIonMatrix(), tempArray.length > 1 ? tempArray[1].getVarPTMs() : null, tempArray.length > 1 ? tempArray[1].getIonMatrix() : null, ms2Tolerance, peptide.length())));
        }
    }
}
