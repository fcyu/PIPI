package proteomics.Segment;


import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Index.BuildIndex;
import proteomics.Types.*;

import java.util.*;

public class InferenceSegment {

    private static final Logger logger = LoggerFactory.getLogger(InferenceSegment.class);
    private static final float precision = 0.0001f;
    private static final int minTagNum = 200;
    private static final float maxTagIntensity = 4; // there are in total 4 peaks whose highest intensity is 1.
    private static final int regionNum = 10;
    private static final int topNumInEachRegion = 20;

    private final Float[] deltaMassArray;
    private final float ms2Tolerance;
    private final boolean lowResolution;
    private TreeMap<Segment, Integer> aaVectorTemplate = new TreeMap<>();
    private Map<Float, Character> massAaMap = new HashMap<>();
    private final int tagLength;
    private final Map<String, Float> massTable;

    public InferenceSegment(BuildIndex buildIndexObj, float ms2Tolerance, int tagLength) throws Exception {
        this.ms2Tolerance = ms2Tolerance;
        this.tagLength = tagLength;
        lowResolution = (ms2Tolerance > 0.1);
        massTable = buildIndexObj.returnMassToolObj().returnMassTable();

        char[] standardAaArray = new char[]{'G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'Q', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W', 'U', 'O'};

        for (char aa : standardAaArray) {
            // # = I/L. $ = Q/K in low resolution.
            if (aa == 'I' || aa == 'L') {
                massAaMap.put(massTable.get(String.valueOf(aa)), '#');
            } else if (((aa == 'Q') || (aa == 'K')) && lowResolution) {
                massAaMap.put(massTable.get("K"), '$');
            } else {
                massAaMap.put(massTable.get(String.valueOf(aa)), aa);
            }
        }
        deltaMassArray = massAaMap.keySet().toArray(new Float[massAaMap.size()]);
        Character[] aaArray = massAaMap.values().toArray(new Character[massAaMap.size()]);

        switch (tagLength) {
            case 1:
                // tag-1
                for (char aa : aaArray) {
                    Segment segment = new Segment(String.valueOf(aa));
                    if (!aaVectorTemplate.containsKey(segment)) {
                        aaVectorTemplate.put(segment, 0);
                    }
                }
                break;
            case 2:
                // tag-2
                for (char aa1 : aaArray) {
                    for (char aa2 : aaArray) {
                        String segmentString = String.valueOf(aa1) + String.valueOf(aa2);
                        Segment segment = new Segment(segmentString);
                        if (!aaVectorTemplate.containsKey(segment)) {
                            aaVectorTemplate.put(segment, 0);
                        }
                    }
                }
                break;
            case 3:
                // tag-3
                for (char aa1 : aaArray) {
                    for (char aa2 : aaArray) {
                        for (char aa3 : aaArray) {
                            String segmentString = String.valueOf(aa1) + String.valueOf(aa2) + String.valueOf(aa3);
                            Segment segment = new Segment(segmentString);
                            if (!aaVectorTemplate.containsKey(segment)) {
                                aaVectorTemplate.put(segment, 0);
                            }
                        }
                    }
                }
                break;
            case 4:
                // tag-4
                for (char aa1 : aaArray) {
                    for (char aa2 : aaArray) {
                        for (char aa3 : aaArray) {
                            for (char aa4 : aaArray) {
                                String segmentString = String.valueOf(aa1) + String.valueOf(aa2) + String.valueOf(aa3) + String.valueOf(aa4);
                                Segment segment = new Segment(segmentString);
                                if (!aaVectorTemplate.containsKey(segment)) {
                                    aaVectorTemplate.put(segment, 0);
                                }
                            }
                        }
                    }
                }
                break;
            default:
                throw new Exception("[ERROR]: The length of tag is larger than 4!");
        }

        int idx = 0;
        for (Segment segment : aaVectorTemplate.keySet()) {
            aaVectorTemplate.put(segment, idx);
            ++idx;
        }
    }

    public List<ThreeExpAA> inferSegmentLocationFromSpectrum(SpectrumEntry spectrumEntry) {
        // Add four virtual peak.
        TreeMap<Float, Float> plMap = spectrumEntry.getPlMapWithVirtualPeaks();
        if (plMap == null) {
            plMap = addVirtualPeaks(spectrumEntry);
            spectrumEntry.setPlMapWithVirtualPeaks(plMap);
        }

        return inferThreeAAFromSpectrum(plMap);
    }

    public Set<Segment> cutTheoSegment(String peptide) {
        String normalizedPeptide = normalizeSequence(lowResolution, peptide);
        Set<Segment> segmentSet = new HashSet<>();
        if (normalizedPeptide.length() == tagLength) {
            segmentSet.add(new Segment(normalizedPeptide));
        } else if (normalizedPeptide.length() > tagLength) {
            for (int i = 0; i <= normalizedPeptide.length() - tagLength; ++i) {
                segmentSet.add(new Segment(normalizedPeptide.substring(i, i + tagLength)));
            }
        }
        return segmentSet;
    }

    public SparseVector generateSegmentIntensityVector(List<ThreeExpAA> inputList) {
        SparseVector finalVector = new SparseVector();
        if (inputList.isEmpty()) {
            return finalVector;
        } else {
            float[] totalIntensityArray = new float[inputList.size()];
            int idx = 0;
            for (ThreeExpAA expAas : inputList) {
                totalIntensityArray[idx] = expAas.getTotalIntensity();
                ++idx;
            }
            Arrays.sort(totalIntensityArray);

            for (ThreeExpAA expAaList : inputList) {
                float totalIntensity = expAaList.getTotalIntensity();
                idx = aaVectorTemplate.get(new Segment(expAaList.getAAString()));
                float value = Math.max(totalIntensity, finalVector.get(idx));
                finalVector.put(idx, value);
            }
            return finalVector;
        }
    }

    public SparseBooleanVector generateSegmentBooleanVector(Set<Segment> cutSegmentSet) {
        SparseBooleanVector finalVector = new SparseBooleanVector();
        if (cutSegmentSet.isEmpty()) {
            return finalVector;
        } else {
            for (Segment segment : cutSegmentSet) {
                finalVector.put(aaVectorTemplate.get(segment));
            }
            return finalVector;
        }
    }

    public static String normalizeSequence(boolean lowResolution, String seq) {
        String normalizedSeq = seq.replaceAll("[IL]", "#");
        if (lowResolution) {
            normalizedSeq = normalizedSeq.replaceAll("[QK]", "$");
        }
        return normalizedSeq;
    }

    public int sparseVectorLength() {
        return aaVectorTemplate.size();
    }

    private List<ThreeExpAA> inferThreeAAFromSpectrum(TreeMap<Float, Float> plMap) { // todo: only support 3-tag
        Float[] mzArray = plMap.keySet().toArray(new Float[plMap.size()]);
        Float[] intensityArray = plMap.values().toArray(new Float[plMap.size()]);
        List<ThreeExpAA> tempList = new LinkedList<>();
        List<ThreeExpAA> outputList = new LinkedList<>();
        for (int i = 0; i < mzArray.length; ++i) {
            float mz1 = mzArray[i];
            float intensity1 = intensityArray[i];
            for (int j = i + 1; j < mzArray.length; ++j) {
                float mz2 = mzArray[j];
                float intensity2 = intensityArray[j];
                char aa1 = inferAA(mz1, mz2);
                if (aa1 != 0) {
                    ExpAA expAa1 = new ExpAA(aa1, mz1, mz2, intensity1, intensity2, -1);
                    List<List<ExpAA>> tempAasList2 = new LinkedList<>();
                    for (int k = j + 1; k < mzArray.length; ++k) {
                        float mz3 = mzArray[k];
                        float intensity3 = intensityArray[k];
                        char aa2 = inferAA(mz2, mz3);
                        if (aa2 != 0) {
                            ExpAA expAa2 = new ExpAA(aa2, mz2, mz3, intensity2, intensity3, -1);
                            List<ExpAA> tempAasList3 = new LinkedList<>();
                            for (int l = k + 1; l < mzArray.length; ++l) {
                                float mz4 = mzArray[l];
                                float intensity4 = intensityArray[l];
                                char aa3 = inferAA(mz3, mz4);
                                if (aa3 != 0) {
                                    ExpAA expAa3 = new ExpAA(aa3, mz3, mz4, intensity3, intensity4, -1);
                                    tempAasList3.add(expAa3);
                                }
                            }
                            for (ExpAA expAas3 : tempAasList3) {
                                List<ExpAA> tempList2 = new LinkedList<>();
                                tempList2.add(expAa2);
                                tempList2.add(expAas3);
                                tempAasList2.add(tempList2);
                            }
                        }
                    }
                    for (List<ExpAA> expAas2 : tempAasList2) {
                        ThreeExpAA threeExpAa = new ThreeExpAA(expAa1, expAas2.get(0), expAas2.get(1));
                        tempList.add(threeExpAa);
                    }
                }
            }
        }


        if (tempList.size() > minTagNum) {
            float minMz = plMap.firstKey();
            float regionWindow = (float) Math.ceil((plMap.lastKey() - minMz) / regionNum);
            for (ThreeExpAA expAa : tempList) {
                expAa.setRegionIdx((int) Math.floor((expAa.getHeadLocation() - minMz) / regionWindow));
            }
            List<List<Float>> regionIntensityList = new ArrayList<>(20);
            for (int i = 0; i < regionNum; ++i) {
                regionIntensityList.add(new ArrayList<>(100));
            }
            for (ThreeExpAA expAa : tempList) {
                regionIntensityList.get(expAa.getRegionIdx()).add(expAa.getTotalIntensity());
            }
            float[] intensityTArray = new float[regionNum];
            for (int i = 0; i < regionNum; ++i) {
                List<Float> intensityList = regionIntensityList.get(i);
                Collections.sort(intensityList, Collections.reverseOrder());
                if (intensityList.size() > topNumInEachRegion) {
                    intensityTArray[i] = intensityList.get(topNumInEachRegion);
                }
            }
            for (ThreeExpAA expAa : tempList) {
                if (expAa.getTotalIntensity() > intensityTArray[expAa.getRegionIdx()]) {
                    outputList.add(expAa);
                }
            }
            return outputList;
        } else {
            return tempList;
        }
    }

    private char inferAA(float mz1, float mz2) {
        float mzDiff = mz2 - mz1;
        for (float mass : deltaMassArray) {
            if (Math.abs(mzDiff - mass) <= ms2Tolerance) {
                return massAaMap.get(mass);
            }
        }
        return 0;
    }

    private TreeMap<Float, Float> addVirtualPeaks(SpectrumEntry spectrumEntry) {
        float totalMass = spectrumEntry.precursorMass + 2 * massTable.get("PROTON");
        TreeMap<Float, Float> plMap = spectrumEntry.plMap;
        TreeMap<Float, Float> finalPlMap = new TreeMap<>();
        for (float mz : plMap.keySet()) {
            finalPlMap.put(mz, plMap.get(mz));
        }
        for (float mz : plMap.keySet()) {
            float anotherMz = totalMass - mz;
            float leftMz = anotherMz - ms2Tolerance;
            float rightMz = anotherMz + ms2Tolerance;
            NavigableMap<Float, Float> temp = null;
            try {
                temp = plMap.subMap(leftMz, true, rightMz, true);
            } catch (IllegalArgumentException ex) {}

            if ((temp == null) || (temp.isEmpty())) {
                finalPlMap.put(anotherMz, plMap.get(mz));
            }
        }

        // Add two virtual peak. Because we have convert all y-ions to b-ions.
        float[] mzList = new float[]{massTable.get("PROTON"), spectrumEntry.precursorMass - massTable.get("H2O") + massTable.get("PROTON")};
        for (float mz : mzList) {
            float leftMz = mz - ms2Tolerance;
            float rightMz = mz + ms2Tolerance;
            NavigableMap<Float, Float> temp = null;
            try {
                temp = plMap.subMap(leftMz, true, rightMz, true);
            } catch (IllegalArgumentException ex) {}

            if ((temp == null) || (temp.isEmpty())) {
                finalPlMap.put(mz, 1f);
            }
        }

        return finalPlMap;
    }

    private double[] calPValue(float[] totalIntensityArray) { // p-values in descending order
        int[] scoreCountArray = new int[(int) (maxTagIntensity / precision) + 1];
        double[] pValueArray = new double[scoreCountArray.length];
        for (float totalIntensity : totalIntensityArray) {
            ++scoreCountArray[intensity2Bin(totalIntensity)];
        }
        int count = scoreCountArray[scoreCountArray.length - 1];
        pValueArray[pValueArray.length - 1] = (double) count / (double) totalIntensityArray.length;
        for (int i = scoreCountArray.length - 2; i >= 0; --i) {
            count += scoreCountArray[i];
            pValueArray[i] = (double) count / (double) totalIntensityArray.length;
        }
        return pValueArray;
    }

    private float getScoreT(double[] pValueArray, double pT) {
        int idx = 0;
        for (int i = 0; i < pValueArray.length; ++i) {
            if (pValueArray[i] <= pT) {
                idx = i;
                break;
            }
        }

        return idx * precision;
    }

    private int intensity2Bin(float totalIntensity) {
        return (int) Math.floor(totalIntensity / precision);
    }
}
