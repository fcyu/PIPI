package proteomics.Segment;


import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Index.BuildIndex;
import proteomics.Types.*;

import java.util.*;

public class InferenceSegment {

    private static final Logger logger = LoggerFactory.getLogger(InferenceSegment.class);
    private static final int minTagNum = 200;
    private static final int regionNum = 10;
    private static final int topNumInEachRegion = 20;

    private final float ms2Tolerance;
    private TreeMap<Segment, Integer> aaVectorTemplate = new TreeMap<>();
    private final Map<String, Float> massTable;
    private Map<Float, String> modifiedAAMap = new HashMap<>();
    private final Float[] deltaMassArray;
    private Map<String, Float> modifiedAAMassMap = new HashMap<>();

    public InferenceSegment(BuildIndex buildIndexObj, float ms2Tolerance, Map<String, String> parameterMap) throws Exception {
        this.ms2Tolerance = ms2Tolerance;
        massTable = buildIndexObj.returnMassToolObj().returnMassTable();

        char[] standardAaArray = new char[]{'G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'Q', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W', 'U', 'O'};

        Map<Float, Character> massAaMap = new HashMap<>();
        for (char aa : standardAaArray) {
            // # = I/L.
            if (aa == 'I' || aa == 'L') {
                massAaMap.put(massTable.get(String.valueOf(aa)), '#');
            } else {
                massAaMap.put(massTable.get(String.valueOf(aa)), aa);
            }
        }

        Character[] aaArray = massAaMap.values().toArray(new Character[massAaMap.size()]);

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

        int idx = 0;
        for (Segment segment : aaVectorTemplate.keySet()) {
            aaVectorTemplate.put(segment, idx);
            ++idx;
        }

        // generate a mass aa map containing modified amino acid
        for (float k : massAaMap.keySet()) {
            modifiedAAMap.put(k, massAaMap.get(k).toString());
        }
        for (String k : parameterMap.keySet()) {
            if (k.startsWith("mod")) {
                String v = parameterMap.get(k);
                if (!v.startsWith("0.0")) {
                    String[] temp = v.split("@");
                    float tempMass = massTable.get(temp[1].substring(0, 1)) + Float.valueOf(temp[0]);
                    // check if the mass has conflict
                    for (float temp2 : modifiedAAMap.keySet()) {
                        if (Math.abs(temp2 - tempMass) <= ms2Tolerance) {
                            logger.error("{} and {} have conflict mass values({} vs {}).", v, modifiedAAMap.get(temp2), tempMass, temp2);
                            System.exit(1);
                        }
                    }
                    if ((temp[1].charAt(0) == 'I') || (temp[1].charAt(0) == 'L')) {
                        modifiedAAMap.put(tempMass, temp[1].replace(temp[1].charAt(0), '#'));
                        modifiedAAMassMap.put(temp[1].replace(temp[1].charAt(0), '#'), Float.valueOf(temp[0]));
                    } else {
                        modifiedAAMap.put(tempMass, temp[1]);
                        modifiedAAMassMap.put(temp[1], Float.valueOf(temp[0]));
                    }
                }
            }
        }
        deltaMassArray = modifiedAAMap.keySet().toArray(new Float[modifiedAAMap.size()]);
    }

    public List<ThreeExpAA> inferSegmentLocationFromSpectrum(SpectrumEntry spectrumEntry) {
        return inferThreeAAFromSpectrum(addVirtualPeaks(spectrumEntry));
    }

    public Set<Segment> cutTheoSegment(String peptide) {
        String normalizedPeptide = normalizeSequence(peptide);
        Set<Segment> segmentSet = new HashSet<>();
        if (normalizedPeptide.length() == 3) {
            segmentSet.add(new Segment(normalizedPeptide));
        } else if (normalizedPeptide.length() > 3) {
            for (int i = 0; i <= normalizedPeptide.length() - 3; ++i) {
                segmentSet.add(new Segment(normalizedPeptide.substring(i, i + 3)));
            }
        }
        return segmentSet;
    }

    public SparseVector generateSegmentIntensityVector(List<ThreeExpAA> inputList) {
        SparseVector finalVector = new SparseVector();
        if (inputList.isEmpty()) {
            return finalVector;
        } else {
            for (ThreeExpAA expAaList : inputList) {
                float totalIntensity = expAaList.getTotalIntensity();
                int idx = aaVectorTemplate.get(new Segment(expAaList.getPtmFreeAAString()));
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

    public static String normalizeSequence(String seq) {
        return seq.replaceAll("[IL]", "#");
    }

    private List<ThreeExpAA> inferThreeAAFromSpectrum(TreeMap<Float, Float> plMap) {
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
                String aa1 = inferAA(mz1, mz2);
                if (aa1 != null) {
                    float mod = 0;
                    if (modifiedAAMassMap.containsKey(aa1)) {
                        mod = modifiedAAMassMap.get(aa1);
                    }
                    ExpAA expAa1 = new ExpAA(aa1, mz1, mz2, intensity1, intensity2, -1, mod);
                    List<List<ExpAA>> tempAasList2 = new LinkedList<>();
                    for (int k = j + 1; k < mzArray.length; ++k) {
                        float mz3 = mzArray[k];
                        float intensity3 = intensityArray[k];
                        String aa2 = inferAA(mz2, mz3);
                        if (aa2 != null) {
                            mod = 0;
                            if (modifiedAAMassMap.containsKey(aa2)) {
                                mod = modifiedAAMassMap.get(aa2);
                            }
                            ExpAA expAa2 = new ExpAA(aa2, mz2, mz3, intensity2, intensity3, -1, mod);
                            List<ExpAA> tempAasList3 = new LinkedList<>();
                            for (int l = k + 1; l < mzArray.length; ++l) {
                                float mz4 = mzArray[l];
                                float intensity4 = intensityArray[l];
                                String aa3 = inferAA(mz3, mz4);
                                if (aa3 != null) {
                                    mod = 0;
                                    if (modifiedAAMassMap.containsKey(aa3)) {
                                        mod = modifiedAAMassMap.get(aa3);
                                    }
                                    ExpAA expAa3 = new ExpAA(aa3, mz3, mz4, intensity3, intensity4, -1, mod);
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

    private String inferAA(float mz1, float mz2) {
        float mzDiff = mz2 - mz1;
        for (float mass : deltaMassArray) {
            if (Math.abs(mzDiff - mass) <= ms2Tolerance) {
                return modifiedAAMap.get(mass);
            }
        }
        return null;
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
        finalPlMap.put(massTable.get("PROTON"), 1f);
        float cTermMz = spectrumEntry.precursorMass - massTable.get("H2O") + massTable.get("PROTON");
        float leftMz = cTermMz - ms2Tolerance;
        float rightMz = cTermMz + ms2Tolerance;
        NavigableMap<Float, Float> temp = null;
        try {
            temp = plMap.subMap(leftMz, true, rightMz, true);
        } catch (IllegalArgumentException ex) {}
        if ((temp == null) || (temp.isEmpty())) {
            finalPlMap.put(cTermMz, 1f);
        }

        return finalPlMap;
    }
}
