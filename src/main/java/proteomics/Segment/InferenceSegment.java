package proteomics.Segment;


import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class InferenceSegment {

    private static final Logger logger = LoggerFactory.getLogger(InferenceSegment.class);
    private static final int minTagNum = 200;
    private static final int regionNum = 10;
    private static final int topNumInEachRegion = 20;
    private static final Pattern pattern = Pattern.compile("([nc][0-9a-i])?([A-Z#$].?)");

    private final float ms2Tolerance;
    private TreeMap<Segment, Integer> aaVectorTemplate = new TreeMap<>();
    private Map<Float, String> modifiedAAMap = new HashMap<>(35, 1);
    private final Float[] deltaMassArray;
    private Map<String, Float> modifiedAAMassMap = new HashMap<>(35, 1);
    private Set<VarModParam> varModParamSet = new HashSet<>();
    private float[] pepNTermPossibleMod = null;
    private float[] pepCTermPossibleMod = null;
    private float[] proNTermPossibleMod = null;
    private float[] proCTermPossibleMod = null;

    public InferenceSegment(Map<Character, Float> massTable, float ms2Tolerance, Map<String, String> parameterMap, Map<Character, Float> fixModMap) throws Exception {
        this.ms2Tolerance = ms2Tolerance;

        char[] standardAaArray = new char[]{'G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'Q', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W', 'U', 'O'};

        Map<Float, Character> massAaMap = new HashMap<>(25, 1);
        for (char aa : standardAaArray) {
            // # = I/L.
            if (aa == 'I' || aa == 'L') {
                massAaMap.put(massTable.get(aa), '#');
            } else {
                massAaMap.put(massTable.get(aa), aa);
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
                    float tempMass = massTable.get(temp[1].charAt(0)) + Float.valueOf(temp[0]);
                    // check if the mass has conflict
                    for (float temp2 : modifiedAAMap.keySet()) {
                        if (Math.abs(temp2 - tempMass) <= ms2Tolerance) {
                            logger.error("{} and {} have conflict mass values({} vs {}).", v, modifiedAAMap.get(temp2), tempMass, temp2);
                            System.exit(1);
                        }
                    }
                    if (Math.abs(fixModMap.get(temp[1].charAt(0))) < 0.1) {
                        // fix modification and var modification cannot be coexist
                        if ((temp[1].charAt(0) == 'I') || (temp[1].charAt(0) == 'L')) {
                            modifiedAAMap.put(tempMass, temp[1].replace(temp[1].charAt(0), '#'));
                            modifiedAAMassMap.put(temp[1].replace(temp[1].charAt(0), '#'), Float.valueOf(temp[0]));
                        } else {
                            modifiedAAMap.put(tempMass, temp[1]);
                            modifiedAAMassMap.put(temp[1], Float.valueOf(temp[0]));
                        }
                        varModParamSet.add(new VarModParam(Float.valueOf(temp[0]), temp[1].charAt(0)));
                    }
                }
            } else if (k.contentEquals("pepNterm")) {
                if (Math.abs(fixModMap.get('n')) < 0.1) {
                    // fix modification and var modification cannot be coexist
                    if (!parameterMap.get(k).startsWith("0.0")) {
                        String[] tempArray = parameterMap.get(k).split(",");
                        pepNTermPossibleMod = new float[tempArray.length];
                        for (int i = 0; i < tempArray.length; ++i) {
                            pepNTermPossibleMod[i] = Float.valueOf(tempArray[i].trim());
                            varModParamSet.add(new VarModParam(Float.valueOf(tempArray[i].trim()), 'n'));
                        }
                    }
                }
            } else if (k.contentEquals("pepCterm")) {
                if (Math.abs(fixModMap.get('c')) < 0.1) {
                    // fix modification and var modification cannot be coexist
                    if (!parameterMap.get(k).startsWith("0.0")) {
                        String[] tempArray = parameterMap.get(k).split(",");
                        pepCTermPossibleMod = new float[tempArray.length];
                        for (int i = 0; i < tempArray.length; ++i) {
                            pepCTermPossibleMod[i] = Float.valueOf(tempArray[i].trim());
                            varModParamSet.add(new VarModParam(Float.valueOf(tempArray[i].trim()), 'c'));
                        }
                    }
                }
            } else if (k.contentEquals("proNterm")) {
                if (Math.abs(fixModMap.get('n')) < 0.1) {
                    if (!parameterMap.get(k).startsWith("0.0")) {
                        String[] tempArray = parameterMap.get(k).split(",");
                        proNTermPossibleMod = new float[tempArray.length];
                        for (int i = 0; i < tempArray.length; ++i) {
                            proNTermPossibleMod[i] = Float.valueOf(tempArray[i].trim());
                            varModParamSet.add(new VarModParam(Float.valueOf(tempArray[i].trim()), 'n')); // to be improve: we don't distinguish protein and peptide N-term
                        }
                    }
                }
            } else if (k.contentEquals("proCterm")) {
                if (Math.abs(fixModMap.get('c')) < 0.1) {
                    if (!parameterMap.get(k).startsWith("0.0")) {
                        String[] tempArray = parameterMap.get(k).split(",");
                        proCTermPossibleMod = new float[tempArray.length];
                        for (int i = 0; i < tempArray.length; ++i) {
                            proCTermPossibleMod[i] = Float.valueOf(tempArray[i].trim());
                            varModParamSet.add(new VarModParam(Float.valueOf(tempArray[i].trim()), 'c')); // to be improve: we don't distinguish protein and peptide C-term
                        }
                    }
                }
            }
        }
        deltaMassArray = modifiedAAMap.keySet().toArray(new Float[modifiedAAMap.size()]);
    }

    public List<ThreeExpAA> inferSegmentLocationFromSpectrum(SpectrumEntry spectrumEntry) {
        return inferThreeAAFromSpectrum(addVirtualPeaks(spectrumEntry), spectrumEntry.precursorMass - MassTool.H2O + MassTool.PROTON);
    }

    public Set<Segment> cutTheoSegment(String peptide) {
        String normalizedPeptide = normalizeSequence(peptide);
        Set<Segment> segmentSet = new HashSet<>(peptide.length() + 1, 1);
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

    public Map<String, Float> getModifiedAAMassMap() {
        return modifiedAAMassMap;
    }

    public Set<VarModParam> getVarModParamSet() {
        return varModParamSet;
    }

    public float[] getPepNTermPossibleMod() {
        return pepNTermPossibleMod;
    }

    public float[] getPepCTermPossibleMod() {
        return pepCTermPossibleMod;
    }

    public float[] getProNTermPossibleMod() {
        return proNTermPossibleMod;
    }

    public float[] getProCTermPossibleMod() {
        return proCTermPossibleMod;
    }

    public static String normalizeSequence(String seq) {
        return seq.replaceAll("[IL]", "#");
    }

    private List<ThreeExpAA> inferThreeAAFromSpectrum(TreeMap<Float, Float> plMap, float cTermMz) {
        Float[] mzArray = plMap.keySet().toArray(new Float[plMap.size()]);
        Float[] intensityArray = plMap.values().toArray(new Float[plMap.size()]);
        Set<ThreeExpAA> tempSet = new HashSet<>();
        List<ThreeExpAA> outputList = new LinkedList<>();
        for (int i = 0; i < mzArray.length; ++i) {
            float mz1 = mzArray[i];
            float intensity1 = intensityArray[i];
            for (int j = i + 1; j < mzArray.length; ++j) {
                float mz2 = mzArray[j];
                float intensity2 = intensityArray[j];
                String aa1 = inferAA(mz1, mz2, Math.abs(mz1 - MassTool.PROTON) <= ms2Tolerance, false, Math.abs(mz1 - MassTool.PROTON) <= ms2Tolerance, false);
                if (aa1 != null) {
                    Matcher matcher = pattern.matcher(aa1);
                    char ptmFreeAA = '\0';
                    float mod = 0;
                    float nTermMod = 0;
                    if (matcher.matches()) {
                        if (modifiedAAMassMap.containsKey(matcher.group(2))) {
                            mod = modifiedAAMassMap.get(matcher.group(2));
                        }
                        ptmFreeAA = matcher.group(2).charAt(0);
                        if (matcher.group(1) != null) {
                            if ((matcher.group(1).charAt(1) - '0' >= 0) && (matcher.group(1).charAt(1) - '0' < 10)) {
                                nTermMod = pepNTermPossibleMod[matcher.group(1).charAt(1) - '0'];
                            } else if ((matcher.group(1).charAt(1) - 'a' >= 0) && (matcher.group(1).charAt(1) - 'a' < 10)) {
                                nTermMod = proNTermPossibleMod[matcher.group(1).charAt(1) - 'a'];
                            } else {
                                logger.error("Something is wrong in inferring tags.");
                                System.exit(1);
                            }
                        }
                    } else {
                        logger.error("Cannot find the PTM free amino acid for {}.", aa1);
                        System.exit(1);
                    }
                    ExpAA expAa1 = new ExpAA(aa1, ptmFreeAA, mz1, mz2, intensity1, intensity2, -1, mod, nTermMod, 0);
                    List<List<ExpAA>> tempAasList2 = new LinkedList<>();
                    for (int k = j + 1; k < mzArray.length; ++k) {
                        float mz3 = mzArray[k];
                        float intensity3 = intensityArray[k];
                        String aa2 = inferAA(mz2, mz3, false, false, false, false);
                        if (aa2 != null) {
                            mod = 0;
                            if (modifiedAAMassMap.containsKey(aa2)) {
                                mod = modifiedAAMassMap.get(aa2);
                            }
                            ExpAA expAa2 = new ExpAA(aa2, aa2.charAt(0), mz2, mz3, intensity2, intensity3, -1, mod, 0, 0);
                            List<ExpAA> tempAasList3 = new LinkedList<>();
                            for (int l = k + 1; l < mzArray.length; ++l) {
                                float mz4 = mzArray[l];
                                float intensity4 = intensityArray[l];
                                String aa3 = inferAA(mz3, mz4, false, Math.abs(mz4 - cTermMz) <= ms2Tolerance, false, Math.abs(mz4 - cTermMz) <= ms2Tolerance);
                                if (aa3 != null) {
                                    matcher = pattern.matcher(aa3);
                                    ptmFreeAA = '\0';
                                    mod = 0;
                                    float cTermMod = 0;
                                    if (matcher.matches()) {
                                        if (modifiedAAMassMap.containsKey(matcher.group(2))) {
                                            mod = modifiedAAMassMap.get(matcher.group(2));
                                        }
                                        ptmFreeAA = matcher.group(2).charAt(0);
                                        if (matcher.group(1) != null) {
                                            if ((matcher.group(1).charAt(1) - '0' >= 0) && (matcher.group(1).charAt(1) - '0' < 10)) {
                                                cTermMod = pepCTermPossibleMod[matcher.group(1).charAt(1) - '0'];
                                            } else if ((matcher.group(1).charAt(1) - 'a' >= 0) && (matcher.group(1).charAt(1) - 'a' < 10)) {
                                                cTermMod = proCTermPossibleMod[matcher.group(1).charAt(1) - 'a'];
                                            } else {
                                                logger.error("Something is wrong in inferring tags.");
                                                System.exit(1);
                                            }
                                        }
                                    } else {
                                        logger.error("Cannot find the PTM free amino acid for {}.", aa3);
                                        System.exit(1);
                                    }
                                    ExpAA expAa3 = new ExpAA(aa3, ptmFreeAA, mz3, mz4, intensity3, intensity4, -1, mod, 0, cTermMod);
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
                        tempSet.add(threeExpAa);
                    }
                }
            }
        }

        // eliminate "overlapped" tags
        ThreeExpAA[] tempArray = tempSet.toArray(new ThreeExpAA[tempSet.size()]);
        List<ThreeExpAA> tempList = new LinkedList<>();
        for (int i = 0; i < tempArray.length; ++i) {
            boolean keep = true;
            for (int j = 0; j < tempArray.length; ++j) {
                if (i != j) {
                    if (tempArray[i].approximateEquals(tempArray[j], 2 * ms2Tolerance)) {
                        if (tempArray[i].getTotalIntensity() < tempArray[j].getTotalIntensity()) {
                            keep = false;
                            break;
                        }
                    }
                }
            }
            if (keep) {
                tempList.add(tempArray[i]);
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

    private String inferAA(float mz1, float mz2, boolean pepNTerm, boolean pepCTerm, boolean proNTerm, boolean proCTerm) {
        float mzDiff = mz2 - mz1;
        for (float mass : deltaMassArray) {
            if (Math.abs(mzDiff - mass) <= 2 * ms2Tolerance) {
                return modifiedAAMap.get(mass);
            }
        }

        if (pepNTerm && (pepNTermPossibleMod != null)) {
            for (float mass : deltaMassArray) {
                for (int i = 0; i < pepNTermPossibleMod.length; ++i) {
                    if (Math.abs(mzDiff - mass - pepNTermPossibleMod[i]) <= 2 * ms2Tolerance) {
                        return "n" + i + modifiedAAMap.get(mass);
                    }
                }
            }
        }

        if (pepCTerm && (pepCTermPossibleMod != null)) {
            for (float mass : deltaMassArray) {
                for (int i = 0; i < pepCTermPossibleMod.length; ++i) {
                    if (Math.abs(mzDiff - mass - pepCTermPossibleMod[i]) <= 2 * ms2Tolerance) {
                        return "c" + i + modifiedAAMap.get(mass);
                    }
                }
            }
        }

        if (proNTerm && (proNTermPossibleMod != null)) {
            for (float mass : deltaMassArray) {
                for (int i = 0; i < proNTermPossibleMod.length; ++i) {
                    if (Math.abs(mzDiff - mass - proNTermPossibleMod[i]) <= 2 * ms2Tolerance) {
                        return "n" + (char) (i + 'a') + modifiedAAMap.get(mass);
                    }
                }
            }
        }

        if (proCTerm && (proCTermPossibleMod != null)) {
            for (float mass : deltaMassArray) {
                for (int i = 0; i < proCTermPossibleMod.length; ++i) {
                    if (Math.abs(mzDiff - mass - proCTermPossibleMod[i]) <= 2 * ms2Tolerance) {
                        return "c" + (char) (i + 'a') + modifiedAAMap.get(mass);
                    }
                }
            }
        }

        return null;
    }

    private TreeMap<Float, Float> addVirtualPeaks(SpectrumEntry spectrumEntry) {
        float totalMass = spectrumEntry.precursorMass + 2 * MassTool.PROTON;
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
        finalPlMap.put(MassTool.PROTON, 1f);
        float cTermMz = spectrumEntry.precursorMass - MassTool.H2O + MassTool.PROTON;
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
