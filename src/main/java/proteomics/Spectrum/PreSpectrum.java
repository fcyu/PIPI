package proteomics.Spectrum;

import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Types.*;

import java.util.*;

public class PreSpectrum {

    private static final double defaultIntensity = 1; // DO NOT change. Otherwise, change the whole project accordingly.
    private static final double floatZero = 1e-6;
    private static final int xcorrOffset = 75;
    public static final int topN = 10;
    private static final double removePrecursorPeakTolerance = 1.5; // this equals the isolation window.

    private final MassTool massToolObj;

    public PreSpectrum(MassTool massToolObj) {
        this.massToolObj = massToolObj;
    }

    public TreeMap<Double, Double> preSpectrum (Map<Double, Double> peaksMap, double precursorMass, int precursorCharge, double ms2Tolerance, double minClear, double maxClear) {
        // remove precursor peak from spectrum
        TreeMap<Double, Double> temp = removeCertainPeaks(peaksMap, precursorMass, precursorCharge, ms2Tolerance, minClear, maxClear);

        temp = new TreeMap<>(temp.subMap(0d, precursorMass));

        return preprocess(temp);
    }

    public SparseVector prepareXcorr(TreeMap<Double, Double> plMap, boolean preprocess) {
        double[] plArray;
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
                xcorrPl.put(i, temp);
            }
        }

        return xcorrPl;
    }

    public SparseVector prepareDigitizedPL(TreeMap<Double, Double> plMap, boolean preprocess) {
        double[] plArray;
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

    public static TreeMap<Double, Double> selectTopN(TreeMap<Double, Double> plMap, int localTopN) { // the input plMap must have been normalized.
        // select top N in each 100 Da
        TreeMap<Double, Double> preprocessedPL = new TreeMap<>();
        double minMz = plMap.firstKey();
        double maxMz = plMap.lastKey();
        double leftMz = minMz;
        while (leftMz < maxMz) {
            // find the max intensity in each window
            double rightMz = Math.min(leftMz + 100, maxMz);
            NavigableMap<Double, Double> subMap;
            if (rightMz < maxMz) {
                subMap = plMap.subMap(leftMz, true, rightMz, false);
            } else {
                subMap = plMap.subMap(leftMz, true, rightMz, true);
            }
            if (!subMap.isEmpty()) {
                Double[] intensityArray = subMap.values().toArray(new Double[subMap.size()]);
                Arrays.sort(intensityArray, Comparator.reverseOrder());
                double temp2 = subMap.size() > localTopN ? intensityArray[localTopN] : 0;
                for (double mz : subMap.keySet()) {
                    if (subMap.get(mz) > temp2) {
                        preprocessedPL.put(mz, subMap.get(mz));
                    }
                }
            }
            leftMz = rightMz;
        }

        return preprocessedPL;
    }

    private TreeMap<Double, Double> removeCertainPeaks(Map<Double, Double> peakMap, double precursorMass, int precursorCharge, double ms2Tolerance, double minClear, double maxClear) {
        TreeMap<Double, Double> mzIntensityMap = new TreeMap<>();
        double precursorMz = precursorMass / precursorCharge + MassTool.PROTON;
        for (double mz : peakMap.keySet()) {
            if (((mz < minClear) || (mz > maxClear)) && (mz > 50)) {
                if ((peakMap.get(mz) > floatZero) && (Math.abs(peakMap.get(mz) - precursorMz) > removePrecursorPeakTolerance)) {
                    mzIntensityMap.put(mz, peakMap.get(mz));
                }
            }
        }

        return mzIntensityMap;
    }

    private TreeMap<Double, Double> deNoise(TreeMap<Double, Double> plMap) {
        // denoise
        TreeMap<Double, Double> denoisedPlMap = new TreeMap<>();
        double minMz = plMap.firstKey();
        double maxMz = plMap.lastKey();
        double windowSize = (plMap.lastKey() - plMap.firstKey()) * 0.1 + 1;
        for (int i = 0; i < 10; ++i) {
            double leftMz = Math.min(minMz + i * windowSize, maxMz);
            double rightMz = Math.min(leftMz + windowSize, maxMz);
            NavigableMap<Double, Double> subPlMap;
            if (rightMz < maxMz) {
                subPlMap = plMap.subMap(leftMz, true, rightMz, false);
            } else {
                subPlMap = plMap.subMap(leftMz, true, rightMz, true);
            }

            if (subPlMap.size() > 9) {
                double noiseIntensity = estimateNoiseIntensity(subPlMap);
                for (double mz : subPlMap.keySet()) {
                    if (subPlMap.get(mz) > noiseIntensity) {
                        denoisedPlMap.put(mz, subPlMap.get(mz));
                    }
                }
            } else {
                for (double mz : subPlMap.keySet()) {
                    denoisedPlMap.put(mz, subPlMap.get(mz));
                }
            }
        }

        return denoisedPlMap;
    }

    private double estimateNoiseIntensity(Map<Double, Double> pl) {
        Set<Double> intensitySet = new HashSet<>();
        for (double intensity : pl.values()) {
            intensitySet.add(intensity);
        }
        Double[] uniqueIntensityVector = intensitySet.toArray(new Double[intensitySet.size()]);
        Arrays.sort(uniqueIntensityVector);
        double[] cum = new double[uniqueIntensityVector.length];
        for (int i = 0; i < uniqueIntensityVector.length; ++i) {
            for (double intensity : pl.values()) {
                if (intensity <= uniqueIntensityVector[i]) {
                    ++cum[i];
                }
            }
        }
        double[][] diff = new double[2][uniqueIntensityVector.length - 1];
        for (int i = 0; i < uniqueIntensityVector.length - 1; ++i) {
            diff[0][i] = cum[i + 1] - cum[i];
            diff[1][i] = uniqueIntensityVector[i + 1] - uniqueIntensityVector[i];
        }
        double[] diff2 = new double[uniqueIntensityVector.length - 1];
        for (int i = 0; i < uniqueIntensityVector.length - 1; ++i) {
            diff2[i] = diff[0][i] / (diff[1][i] + floatZero);
        }
        double maxValue = 0;
        int maxIdx = 0;
        for (int i = 0; i < uniqueIntensityVector.length - 1; ++i) {
            if (diff2[i] > maxValue) {
                maxValue = diff2[i];
                maxIdx = i;
            }
        }

        return uniqueIntensityVector[maxIdx];
    }

    private TreeMap<Double, Double> preprocess(TreeMap<Double, Double> plMap) {
        // sqrt the intensity and find the highest intensity.
        TreeMap<Double, Double> sqrtPlMap = new TreeMap<>();
        for (double mz : plMap.keySet()) {
            double sqrtIntensity = Math.sqrt(plMap.get(mz));
            sqrtPlMap.put(mz, sqrtIntensity);
        }

        // normalize the intensity in each 100 Da.
        TreeMap<Double, Double> preprocessedPL = new TreeMap<>();
        double minMz = sqrtPlMap.firstKey();
        double maxMz = sqrtPlMap.lastKey();
        double leftMz = minMz;
        while (leftMz < maxMz) {
            // find the max intensity in each window
            double rightMz = Math.min(leftMz + 100, maxMz);
            NavigableMap<Double, Double> subMap;
            if (rightMz < maxMz) {
                subMap = sqrtPlMap.subMap(leftMz, true, rightMz, false);
            } else {
                subMap = sqrtPlMap.subMap(leftMz, true, rightMz, true);
            }
            if (!subMap.isEmpty()) {
                Double[] intensityArray = subMap.values().toArray(new Double[subMap.size()]);
                Arrays.sort(intensityArray, Comparator.reverseOrder());
                double temp1 = defaultIntensity / intensityArray[0];
                double temp2 = subMap.size() > topN ? intensityArray[topN] : 0;
                for (double mz : subMap.keySet()) {
                    if (subMap.get(mz) > temp2) {
                        preprocessedPL.put(mz, subMap.get(mz) * temp1);
                    }
                }
            }
            leftMz = rightMz;
        }

        return preprocessedPL;
    }

    private double[] digitizeSpec(TreeMap<Double, Double> pl) {
        double[] digitizedPl = new double[massToolObj.mzToBin(pl.lastKey()) + 1];
        for (double mz : pl.keySet()) {
            int idx = massToolObj.mzToBin(mz);
            digitizedPl[idx] = Math.max(pl.get(mz), digitizedPl[idx]);
        }

        return digitizedPl;
    }
}
