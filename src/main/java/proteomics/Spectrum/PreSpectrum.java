package proteomics.Spectrum;

import proteomics.TheoSeq.MassTool;
import proteomics.Types.SparseVector;

import java.util.*;

public class PreSpectrum {

    private static final float defaultIntensity = 1; // DO NOT change. Otherwise, change the whole project accordingly.
    private static final float floatZero = 1e-6f;
    private static final int xcprrOffset = 75;

    private final MassTool massToolObj;
    private final Map<String, Float> massTable;

    public PreSpectrum(MassTool massToolObj) {
        this.massToolObj = massToolObj;
        massTable = massToolObj.returnMassTable();
    }

    public TreeMap<Float, Float> preSpectrum (Map<Double, Double> peaksMap, float precursorMass, int precursorCharge, float ms2Tolerance) {
        // remove precursor peak from spectrum
        TreeMap<Float, Float> temp = removePrecursorPeak(peaksMap, precursorMass, precursorCharge, ms2Tolerance);

        // reduce noise
        TreeMap<Float, Float> deionisedPlMap = deNoise(new TreeMap<>(temp.subMap(0f, precursorMass)));

        // normalize
        return normalizeSpec(deionisedPlMap);
    }

    public SparseVector prepareXcorr(TreeMap<Float, Float> plMap, float precursorMass) {
        float[] plArray = digitizeSpec(plMap, precursorMass);

        SparseVector xcorrPl = new SparseVector();
        float mySum = 0;
        int offsetRange = 2 * xcprrOffset + 1;
        for (int i = 0; i < xcprrOffset; ++i) {
            mySum += plArray[i];
        }

        float factor = 1 / (float) (offsetRange - 1);
        for (int i = xcprrOffset; i < plArray.length + xcprrOffset; ++i) {
            if (i < plArray.length) {
                mySum += plArray[i];
            }
            if (i >= offsetRange) {
                mySum -= plArray[i - offsetRange];
            }
            float temp = (plArray[i - xcprrOffset] - (mySum - plArray[i - xcprrOffset]) * factor);
            if (Math.abs(temp) > floatZero) {
                xcorrPl.put(i - xcprrOffset, temp);
            }
        }

        return xcorrPl;
    }

    private TreeMap<Float, Float> removePrecursorPeak(Map<Double, Double> peakMap, float precursorMass, int precursorCharge, float ms2Tolerance) {
        TreeMap<Float, Float> mzIntensityMap = new TreeMap<>();

        for (double mz : peakMap.keySet()) {
            for (int charge = precursorCharge; charge > 0; --charge) {
                float temp = (precursorMass + charge * massTable.get("PROTON")) / charge;
                if ((peakMap.get(mz) > floatZero) && (Math.abs(peakMap.get(mz) - temp) > ms2Tolerance)) {
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

    private TreeMap<Float, Float> normalizeSpec(TreeMap<Float, Float> plMap) {
        // sqrt the intensity and find the highest intensity.
        TreeMap<Float, Float> sqrtPlMap = new TreeMap<>();
        for (float mz : plMap.keySet()) {
            if (plMap.get(mz) > floatZero) {
                float sqrtIntensity = (float) Math.sqrt(plMap.get(mz));
                sqrtPlMap.put(mz, sqrtIntensity);
            }
        }

        // divide the spectrum into 10 windows and normalize each windows to defaultIntensity
        TreeMap<Float, Float> windowedPlMap = new TreeMap<>();
        float minMz = sqrtPlMap.firstKey();
        float maxMz = sqrtPlMap.lastKey();
        float windowSize = (maxMz - minMz) / 10 + 1;
        for (int i = 0; i < 10; ++i) {
            // find the max intensity in each window
            float leftMz = Math.min(minMz + i * windowSize, maxMz);
            float rightMz = Math.min(leftMz + windowSize, maxMz);
            NavigableMap<Float, Float> subMap;
            if (rightMz < maxMz) {
                subMap = sqrtPlMap.subMap(leftMz, true, rightMz, false);
            } else {
                subMap = sqrtPlMap.subMap(leftMz, true, rightMz, true);
            }
            if (!subMap.isEmpty()) {
                Float[] intensityArray = subMap.values().toArray(new Float[subMap.size()]);
                Arrays.sort(intensityArray);
                float temp1 = defaultIntensity / intensityArray[intensityArray.length - 1];
                float temp2 = (float) 0.05 * intensityArray[intensityArray.length - 1];
                for (float mz : subMap.keySet()) {
                    if (subMap.get(mz) > temp2) {
                        windowedPlMap.put(mz, subMap.get(mz) * temp1);
                    }
                }
            }
        }

        return windowedPlMap;
    }

    private float[] digitizeSpec(TreeMap<Float, Float> pl, float precursorMass) {
        float[] digitizedPl = new float[massToolObj.mzToBin(precursorMass) + 1];
        Set<Float> mzSet = pl.keySet();
        for (float mz : mzSet) {
            int idx = massToolObj.mzToBin(mz);
            digitizedPl[idx] = Math.max(pl.get(mz), digitizedPl[idx]);
        }

        return digitizedPl;
    }
}
