package proteomics.Search;


import proteomics.Types.FinalResultEntry;
import proteomics.Types.Peptide;
import proteomics.Types.SpectrumEntry;

import java.util.*;

public class CalSubscores {

    public CalSubscores(List<FinalResultEntry> psmList, Map<Integer, SpectrumEntry> numSpectrumMap, float ms2Tolerance) {
        for (FinalResultEntry psm : psmList) {
            SpectrumEntry spectrum = numSpectrumMap.get(psm.getScanNum());
            TreeMap<Float, Float> expPl = spectrum.plMap;
            Peptide peptide = psm.getPeptide();
            float[][] ionMatrix = peptide.getIonMatrix();
            int precursorCharge = spectrum.precursorCharge;

            float matchedPeakIntensity = 0;
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
            float totalIntensity = 0;
            for (int i = 0; i < maxRow; ++i) {
                for (int j = 0; j < ionMatrix[0].length; ++j) {
                    for (float mz : expPl.keySet()) {
                        float intensity = expPl.get(mz);
                        totalIntensity += intensity;
                        if (Math.abs(mz - ionMatrix[i][j]) <= ms2Tolerance) {
                            if (intensity > intensityT) {
                                ++matchedHighestPeakNum;
                            }
                            matchedPeakIntensity += intensity;
                            break;
                        }
                    }
                }
            }

            psm.setIonFrac(matchedPeakIntensity / totalIntensity);
            psm.setMatchedHighestIntensityFrac((double) matchedHighestPeakNum / (double) totalIonNum);
        }
    }
}
