package proteomics.FDR;

import proteomics.Types.FinalResultEntry;

import java.util.List;

public class EstimateFDR {

    private static final float precision = 0.1f;

    public EstimateFDR(List<FinalResultEntry> resultList) {
        double minNegativeLog10Evalue = 9999;
        double maxNegativeLog10Evalue = -9999;
        float[] qvalueArray;

        for (FinalResultEntry entry : resultList) {
            if (entry.getNegativeLog10EValue() > maxNegativeLog10Evalue) {
                maxNegativeLog10Evalue = entry.getNegativeLog10EValue();
            }
            if (entry.getNegativeLog10EValue() < minNegativeLog10Evalue) {
                minNegativeLog10Evalue = entry.getNegativeLog10EValue();
            }
        }

        // accumulate counts
        int arrayLength = (int) ((maxNegativeLog10Evalue - minNegativeLog10Evalue) / precision) + 1;
        long[] decoyCountArray = new long[arrayLength];
        long[] targetCountArray = new long[arrayLength];
        float[] fdrArray = new float[arrayLength];
        qvalueArray = new float[arrayLength];

        for (FinalResultEntry entry : resultList) {
            if (entry.isDecoy()) {
                ++decoyCountArray[(int) ((entry.getNegativeLog10EValue() - minNegativeLog10Evalue) / precision)];
            } else {
                ++targetCountArray[(int) ((entry.getNegativeLog10EValue() - minNegativeLog10Evalue) / precision)];
            }
        }

        // calculate FDR
        int targetCount = 0;
        int decoyCount = 0;
        for (int i = arrayLength - 1; i >= 0; --i) {
            targetCount += targetCountArray[i];
            decoyCount += decoyCountArray[i];
            if (targetCount > 0) {
                fdrArray[i] =  Math.min((float) decoyCount / (float) targetCount, 1);;
            }
        }

        // convert FDR to qValue
        float lastQValue = fdrArray[0];
        qvalueArray[0] = lastQValue;
        for (int i = 1; i < arrayLength; ++i) {
            float qValue = fdrArray[i];
            if (qValue > lastQValue) {
                qvalueArray[i] = lastQValue;
            } else {
                qvalueArray[i] = qValue;
                lastQValue = qValue;
            }
        }

        // record results
        for (FinalResultEntry entry : resultList) {
            entry.setQValue(qvalueArray[(int) ((entry.getNegativeLog10EValue() - minNegativeLog10Evalue) / precision)]);
        }
    }
}
