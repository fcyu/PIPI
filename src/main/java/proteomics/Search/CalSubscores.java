package proteomics.Search;


import proteomics.Types.FinalResultEntry;
import proteomics.Types.Peptide;
import proteomics.Types.SpectrumEntry;

import java.util.*;

public class CalSubscores {

    public CalSubscores(FinalResultEntry psm, SpectrumEntry spectrum, float ms2Tolerance) {
        TreeMap<Float, Float> expPl = spectrum.plMap;
        Peptide peptide = psm.getPeptide();
        float[][] ionMatrix = peptide.getIonMatrix();
        int precursorCharge = spectrum.precursorCharge;

        int matchedPeakNum = 0;
        int maxRow = Math.min(ionMatrix.length, 2 * (precursorCharge - 1));
        if (precursorCharge == 1) {
            maxRow = 2;
        }
        int totalIonNum = ionMatrix[0].length * maxRow;
        Float[] intensityArray = expPl.values().toArray(new Float[expPl.size()]);
        Arrays.sort(intensityArray, Collections.reverseOrder());
        float intensityT = 0;
        if (totalIonNum < intensityArray.length) {
            intensityT = intensityArray[totalIonNum];
        }
        int matchedHighestPeakNum = 0;
        for (int i = 0; i < maxRow; ++i) {
            for (int j = 0; j < ionMatrix[0].length; ++j) {
                for (float mz : expPl.keySet()) {
                    if (Math.abs(mz - ionMatrix[i][j]) <= ms2Tolerance) {
                        if (expPl.get(mz) > intensityT) {
                            ++matchedHighestPeakNum;
                        }
                        ++matchedPeakNum;
                        break;
                    }
                }
            }
        }

        psm.setIonFrac((double) matchedPeakNum / (double) totalIonNum);
        psm.setMatchedHighestIntensityFrac((double) matchedHighestPeakNum / (double) totalIonNum);
    }
}
