package proteomics.Spectrum;

import proteomics.TheoSeq.MassTool;
import proteomics.Types.SparseVector;

import java.util.*;

public class PreSpectrum {

    private static final float defaultIntensity = 1; // DO NOT change. Otherwise, change the whole project accordingly.
    private static final float floatZero = 1e-6f;
    private static final int xcorrOffset = 75;
    public static final int topN = 10;
    private static final float removePrecursorPeakTolerance = 1.5f; // this equals the isolation window.

    private final MassTool massToolObj;

    public PreSpectrum(MassTool massToolObj) {
        this.massToolObj = massToolObj;
    }

    public TreeMap<Float, Float> preSpectrum (Map<Double, Double> peaksMap, float precursorMass, int precursorCharge, float ms2Tolerance, float minClear, float maxClear) {
        // remove precursor peak from spectrum
        TreeMap<Float, Float> temp = removeCertainPeaks(peaksMap, precursorMass, precursorCharge, ms2Tolerance, minClear, maxClear);

        temp = new TreeMap<>(temp.subMap(0f, precursorMass));

        return preprocess(temp);
    }

    public SparseVector prepareXcorr(TreeMap<Float, Float> plMap, boolean preprocess) {
        float[] plArray;
        if (preprocess) {
            plArray = digitizeSpec(preprocess(plMap));
        } else {
            plArray = digitizeSpec(plMap);
        }

        SparseVector xcorrPl = new SparseVector();
        int offsetRange = 2 * xcorrOffset + 1;
        double factor = 1 / (double) (offsetRange - 1); // caution: 1/150 rather than 1/151
        double mySum = 0;
        for (int i = 0; i < xcorrOffset; ++i) {
            mySum += plArray[i];
        }

        double[] tempArray = new double[plArray.length];
        for (int i = xcorrOffset; i < plArray.length + xcorrOffset; ++i) {
            if (i < plArray.length) {
                mySum += plArray[i];
            }
            if (i >= offsetRange) {
                mySum -= plArray[i - offsetRange];
            }
            tempArray[i - xcorrOffset] = (mySum - plArray[i - xcorrOffset]) * factor;
        }

        for (int i = 1; i < plArray.length; ++i) {
            double temp = plArray[i] - tempArray[i];
            if (Math.abs(temp) > 1e-6) {
                xcorrPl.put(i, (float) temp);
            }
        }

        return xcorrPl;
    }

    public SparseVector prepareDigitizedPL(TreeMap<Float, Float> plMap, boolean preprocess) {
        float[] plArray;
        if (preprocess) {
            plArray = digitizeSpec(preprocess(plMap));
        } else {
            plArray = digitizeSpec(plMap);
        }

        SparseVector digitizedPL = new SparseVector();

        for (int i = 1; i < plArray.length; ++i) {
            if (Math.abs(plArray[i]) > 1e-6) {
                digitizedPL.put(i, plArray[i]);
            }
        }

        return digitizedPL;
    }

    private TreeMap<Float, Float> removeCertainPeaks(Map<Double, Double> peakMap, float precursorMass, int precursorCharge, float ms2Tolerance, float minClear, float maxClear) {
        TreeMap<Float, Float> mzIntensityMap = new TreeMap<>();
        float precursorMz = precursorMass / precursorCharge + MassTool.PROTON;
        for (double mz : peakMap.keySet()) {
            if (((mz < minClear) || (mz > maxClear)) && (mz > 50)) {
                if ((peakMap.get(mz) > floatZero) && (Math.abs(peakMap.get(mz) - precursorMz) > removePrecursorPeakTolerance)) {
                    mzIntensityMap.put((float) mz, peakMap.get(mz).floatValue());
                }
            }
        }

        return mzIntensityMap;
    }

    private TreeMap<Float, Float> deNoise(TreeMap<Float, Float> plMap) {
        // denoise
        TreeMap<Float, Float> denoisedPlMap = new TreeMap<>();
        float minMz = plMap.firstKey();
        float maxMz = plMap.lastKey();
        float windowSize = (plMap.lastKey() - plMap.firstKey()) / 10 + 1;
        for (int i = 0; i < 10; ++i) {
            float leftMz = Math.min(minMz + i * windowSize, maxMz);
            float rightMz = Math.min(leftMz + windowSize, maxMz);
            NavigableMap<Float, Float> subPlMap;
            if (rightMz < maxMz) {
                subPlMap = plMap.subMap(leftMz, true, rightMz, false);
            } else {
                subPlMap = plMap.subMap(leftMz, true, rightMz, true);
            }

            if (subPlMap.size() > 9) {
                float noiseIntensity = estimateNoiseIntensity(subPlMap);
                for (float mz : subPlMap.keySet()) {
                    if (subPlMap.get(mz) > noiseIntensity) {
                        denoisedPlMap.put(mz, subPlMap.get(mz));
                    }
                }
            } else {
                for (float mz : subPlMap.keySet()) {
                    denoisedPlMap.put(mz, subPlMap.get(mz));
                }
            }
        }

        return denoisedPlMap;
    }

    private float estimateNoiseIntensity(Map<Float, Float> pl) {
        Set<Float> intensitySet = new HashSet<>();
        for (float intensity : pl.values()) {
            intensitySet.add(intensity);
        }
        Float[] uniqueIntensityVector = intensitySet.toArray(new Float[intensitySet.size()]);
        Arrays.sort(uniqueIntensityVector);
        float[] cum = new float[uniqueIntensityVector.length];
        for (int i = 0; i < uniqueIntensityVector.length; ++i) {
            for (float intensity : pl.values()) {
                if (intensity <= uniqueIntensityVector[i]) {
                    ++cum[i];
                }
            }
        }
        float[][] diff = new float[2][uniqueIntensityVector.length - 1];
        for (int i = 0; i < uniqueIntensityVector.length - 1; ++i) {
            diff[0][i] = cum[i + 1] - cum[i];
            diff[1][i] = uniqueIntensityVector[i + 1] - uniqueIntensityVector[i];
        }
        float[] diff2 = new float[uniqueIntensityVector.length - 1];
        for (int i = 0; i < uniqueIntensityVector.length - 1; ++i) {
            diff2[i] = diff[0][i] / (diff[1][i] + floatZero);
        }
        float maxValue = 0;
        int maxIdx = 0;
        for (int i = 0; i < uniqueIntensityVector.length - 1; ++i) {
            if (diff2[i] > maxValue) {
                maxValue = diff2[i];
                maxIdx = i;
            }
        }

        return uniqueIntensityVector[maxIdx];
    }

    private TreeMap<Float, Float> preprocess(TreeMap<Float, Float> plMap) {
        // sqrt the intensity and find the highest intensity.
        TreeMap<Float, Float> sqrtPlMap = new TreeMap<>();
        for (float mz : plMap.keySet()) {
            float sqrtIntensity = (float) Math.sqrt(plMap.get(mz));
            sqrtPlMap.put(mz, sqrtIntensity);
        }

        // normalize the intensity in each 100 Da.
        TreeMap<Float, Float> preprocessedPL = new TreeMap<>();
        float minMz = sqrtPlMap.firstKey();
        float maxMz = sqrtPlMap.lastKey();
        float leftMz = minMz;
        while (leftMz < maxMz) {
            // find the max intensity in each window
            float rightMz = Math.min(leftMz + 100, maxMz);
            NavigableMap<Float, Float> subMap;
            if (rightMz < maxMz) {
                subMap = sqrtPlMap.subMap(leftMz, true, rightMz, false);
            } else {
                subMap = sqrtPlMap.subMap(leftMz, true, rightMz, true);
            }
            if (!subMap.isEmpty()) {
                Float[] intensityArray = subMap.values().toArray(new Float[subMap.size()]);
                Arrays.sort(intensityArray, Comparator.reverseOrder());
                float temp1 = defaultIntensity / intensityArray[0];
                float temp2 = subMap.size() > topN ? intensityArray[topN] : 0;
                for (float mz : subMap.keySet()) {
                    if (subMap.get(mz) > temp2) {
                        preprocessedPL.put(mz, subMap.get(mz) * temp1);
                    }
                }
            }
            leftMz = rightMz;
        }

        return preprocessedPL;
    }

    private float[] digitizeSpec(TreeMap<Float, Float> pl) {
        float[] digitizedPl = new float[massToolObj.mzToBin(pl.lastKey()) + 1];
        for (float mz : pl.keySet()) {
            int idx = massToolObj.mzToBin(mz);
            digitizedPl[idx] = Math.max(pl.get(mz), digitizedPl[idx]);
        }

        return digitizedPl;
    }
}
