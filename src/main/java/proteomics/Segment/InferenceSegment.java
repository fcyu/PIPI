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

    private final double ms2Tolerance;
    private TreeMap<Segment, Integer> aaVectorTemplate = new TreeMap<>();
    private Map<Double, String> modifiedAAMap = new HashMap<>(35, 1);
    private final Double[] deltaMassArray;
    private Map<String, Double> modifiedAAMassMap = new HashMap<>(35, 1);
    private Set<VarModParam> varModParamSet = new HashSet<>();
    private double[] nTermPossibleMod = null;
    private double[] cTermPossibleMod = null;
    private MassTool massTool;

    public InferenceSegment(MassTool massTool, double ms2Tolerance, Map<String, String> parameterMap, Map<Character, Double> fixModMap) throws Exception {
        this.massTool = massTool;
        this.ms2Tolerance = ms2Tolerance;
        Map<Character, Double> massTable = massTool.returnMassTable();

        char[] standardAaArray = new char[]{'G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'Q', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W', 'U', 'O'};

        Map<Double, Character> massAaMap = new HashMap<>(25, 1);
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
        for (double k : massAaMap.keySet()) {
            modifiedAAMap.put(k, massAaMap.get(k).toString());
        }
        for (String k : parameterMap.keySet()) {
            if (k.startsWith("mod")) {
                String v = parameterMap.get(k);
                if (!v.startsWith("0.0")) {
                    String[] temp = v.split("@");
                    double tempMass = massTable.get(temp[1].charAt(0)) + Double.valueOf(temp[0]);
                    // check if the mass has conflict
                    for (double temp2 : modifiedAAMap.keySet()) {
                        if (Math.abs(temp2 - tempMass) <= ms2Tolerance) {
                            throw new Exception(String.format(Locale.US, "%s and %s have conflict mass values(%f vs %f).", v, modifiedAAMap.get(temp2), tempMass, temp2));
                        }
                    }
                    if (Math.abs(fixModMap.get(temp[1].charAt(0))) < 0.1) {
                        // fix modification and var modification cannot be coexist
                        if ((temp[1].charAt(0) == 'I') || (temp[1].charAt(0) == 'L')) {
                            modifiedAAMap.put(tempMass, temp[1].replace(temp[1].charAt(0), '#'));
                            modifiedAAMassMap.put(temp[1].replace(temp[1].charAt(0), '#'), Double.valueOf(temp[0]));
                        } else {
                            modifiedAAMap.put(tempMass, temp[1]);
                            modifiedAAMassMap.put(temp[1], Double.valueOf(temp[0]));
                        }
                        varModParamSet.add(new VarModParam(Double.valueOf(temp[0]), temp[1].charAt(0), 1, false)); // var mods from the parameter file have the highest priority, those PTM can exist in peptide terminal.
                    }
                }
            } else if (k.contentEquals("Nterm")) {
                if (Math.abs(fixModMap.get('n')) < 0.1) {
                    // fix modification and var modification cannot be coexist
                    if (!parameterMap.get(k).startsWith("0.0")) {
                        String[] tempArray = parameterMap.get(k).split(",");
                        nTermPossibleMod = new double[tempArray.length];
                        for (int i = 0; i < tempArray.length; ++i) {
                            nTermPossibleMod[i] = Double.valueOf(tempArray[i].trim());
                            varModParamSet.add(new VarModParam(Double.valueOf(tempArray[i].trim()), 'n', 1, false)); // var mods from the parameter file have the highest priority
                        }
                    }
                }
            } else if (k.contentEquals("Cterm")) {
                if (Math.abs(fixModMap.get('c')) < 0.1) {
                    // fix modification and var modification cannot be coexist
                    if (!parameterMap.get(k).startsWith("0.0")) {
                        String[] tempArray = parameterMap.get(k).split(",");
                        cTermPossibleMod = new double[tempArray.length];
                        for (int i = 0; i < tempArray.length; ++i) {
                            cTermPossibleMod[i] = Double.valueOf(tempArray[i].trim());
                            varModParamSet.add(new VarModParam(Double.valueOf(tempArray[i].trim()), 'c', 1, false)); // var mods from the parameter file have the highest priority
                        }
                    }
                }
            }
        }
        deltaMassArray = modifiedAAMap.keySet().toArray(new Double[modifiedAAMap.size()]);
    }

    public List<ThreeExpAA> inferSegmentLocationFromSpectrum(double precursorMass, TreeMap<Double, Double> plMap) throws Exception {
        return inferThreeAAFromSpectrum(addVirtualPeaks(precursorMass, plMap), precursorMass - massTool.H2O + MassTool.PROTON);
    }

    public SparseVector generateSegmentIntensityVector(List<ThreeExpAA> inputList) {
        SparseVector finalVector = new SparseVector();
        if (inputList.isEmpty()) {
            return finalVector;
        } else {
            for (ThreeExpAA expAaList : inputList) {
                double totalIntensity = expAaList.getTotalIntensity();
                int idx = aaVectorTemplate.get(new Segment(expAaList.getPtmFreeAAString()));
                double value = Math.max(totalIntensity, finalVector.get(idx));
                finalVector.put(idx, value);
            }
            return finalVector;
        }
    }

    public SparseBooleanVector generateSegmentBooleanVector(String peptide) {
        String normalizedPeptide = normalizeSequence(peptide);
        Set<Integer> tempSet = new HashSet<>(peptide.length() + 1, 1);
        for (int i = 0; i <= normalizedPeptide.length() - 3; ++i) {
            tempSet.add(aaVectorTemplate.get(new Segment(normalizedPeptide.substring(i, i + 3))));
        }
        return new SparseBooleanVector(tempSet);
    }

    public Map<String, Double> getModifiedAAMassMap() {
        return modifiedAAMassMap;
    }

    public Set<VarModParam> getVarModParamSet() {
        return varModParamSet;
    }

    public double[] getNTermPossibleMod() {
        return nTermPossibleMod;
    }

    public double[] getCTermPossibleMod() {
        return cTermPossibleMod;
    }

    public static String normalizeSequence(String seq) {
        return seq.replaceAll("[IL]", "#");
    }

    private List<ThreeExpAA> inferThreeAAFromSpectrum(TreeMap<Double, Double> plMap, double cTermMz) throws Exception {
        Double[] mzArray = plMap.keySet().toArray(new Double[plMap.size()]);
        Double[] intensityArray = plMap.values().toArray(new Double[plMap.size()]);
        Set<ThreeExpAA> tempSet = new HashSet<>();
        List<ThreeExpAA> outputList = new LinkedList<>();
        for (int i = 0; i < mzArray.length; ++i) {
            double mz1 = mzArray[i];
            double intensity1 = intensityArray[i];
            for (int j = i + 1; j < mzArray.length; ++j) {
                double mz2 = mzArray[j];
                double intensity2 = intensityArray[j];
                String aa1 = inferAA(mz1, mz2, Math.abs(mz1 - MassTool.PROTON) <= ms2Tolerance, false);
                if (aa1 != null) {
                    Matcher matcher = pattern.matcher(aa1);
                    char ptmFreeAA = '\0';
                    double mod = 0;
                    double nTermMod = 0;
                    if (matcher.matches()) {
                        if (modifiedAAMassMap.containsKey(matcher.group(2))) {
                            mod = modifiedAAMassMap.get(matcher.group(2));
                        }
                        ptmFreeAA = matcher.group(2).charAt(0);
                        if (matcher.group(1) != null) {
                            if ((matcher.group(1).charAt(1) - '0' >= 0) && (matcher.group(1).charAt(1) - '0' < 10)) {
                                nTermMod = nTermPossibleMod[matcher.group(1).charAt(1) - '0'];
                            } else {
                                throw new Exception("Something is wrong in inferring tags.");
                            }
                        }
                    } else {
                        throw new NullPointerException(String.format(Locale.US, "Cannot find the PTM free amino acid for %s.", aa1));
                    }
                    ExpAA expAa1 = new ExpAA(aa1, ptmFreeAA, mz1, mz2, intensity1, intensity2, -1, mod, nTermMod, 0);
                    List<List<ExpAA>> tempAasList2 = new LinkedList<>();
                    for (int k = j + 1; k < mzArray.length; ++k) {
                        double mz3 = mzArray[k];
                        double intensity3 = intensityArray[k];
                        String aa2 = inferAA(mz2, mz3, false, false);
                        if (aa2 != null) {
                            mod = 0;
                            if (modifiedAAMassMap.containsKey(aa2)) {
                                mod = modifiedAAMassMap.get(aa2);
                            }
                            ExpAA expAa2 = new ExpAA(aa2, aa2.charAt(0), mz2, mz3, intensity2, intensity3, -1, mod, 0, 0);
                            List<ExpAA> tempAasList3 = new LinkedList<>();
                            for (int l = k + 1; l < mzArray.length; ++l) {
                                double mz4 = mzArray[l];
                                double intensity4 = intensityArray[l];
                                String aa3 = inferAA(mz3, mz4, false, Math.abs(mz4 - cTermMz) <= ms2Tolerance);
                                if (aa3 != null) {
                                    matcher = pattern.matcher(aa3);
                                    ptmFreeAA = '\0';
                                    mod = 0;
                                    double cTermMod = 0;
                                    if (matcher.matches()) {
                                        if (modifiedAAMassMap.containsKey(matcher.group(2))) {
                                            mod = modifiedAAMassMap.get(matcher.group(2));
                                        }
                                        ptmFreeAA = matcher.group(2).charAt(0);
                                        if (matcher.group(1) != null) {
                                            if ((matcher.group(1).charAt(1) - '0' >= 0) && (matcher.group(1).charAt(1) - '0' < 10)) {
                                                cTermMod = cTermPossibleMod[matcher.group(1).charAt(1) - '0'];
                                            } else {
                                                throw new Exception("Something is wrong in inferring tags.");
                                            }
                                        }
                                    } else {
                                        throw new NullPointerException(String.format(Locale.US, "Cannot find the PTM free amino acid for %s.", aa3));
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
            double minMz = plMap.firstKey();
            double regionWindow = Math.ceil((plMap.lastKey() - minMz) / regionNum);
            for (ThreeExpAA expAa : tempList) {
                expAa.setRegionIdx((int) Math.floor((expAa.getHeadLocation() - minMz) / regionWindow));
            }
            List<List<Double>> regionIntensityList = new ArrayList<>(20);
            for (int i = 0; i < regionNum; ++i) {
                regionIntensityList.add(new ArrayList<>(100));
            }
            for (ThreeExpAA expAa : tempList) {
                regionIntensityList.get(expAa.getRegionIdx()).add(expAa.getTotalIntensity());
            }
            double[] intensityTArray = new double[regionNum];
            for (int i = 0; i < regionNum; ++i) {
                List<Double> intensityList = regionIntensityList.get(i);
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

    private String inferAA(double mz1, double mz2, boolean nTerm, boolean cTerm) {
        double mzDiff = mz2 - mz1;
        for (double mass : deltaMassArray) {
            if (Math.abs(mzDiff - mass) <= 2 * ms2Tolerance) {
                return modifiedAAMap.get(mass);
            }
        }

        if (nTerm && (nTermPossibleMod != null)) {
            for (double mass : deltaMassArray) {
                for (int i = 0; i < nTermPossibleMod.length; ++i) {
                    if (Math.abs(mzDiff - mass - nTermPossibleMod[i]) <= 2 * ms2Tolerance) {
                        return "n" + i + modifiedAAMap.get(mass);
                    }
                }
            }
        }

        if (cTerm && (cTermPossibleMod != null)) {
            for (double mass : deltaMassArray) {
                for (int i = 0; i < cTermPossibleMod.length; ++i) {
                    if (Math.abs(mzDiff - mass - cTermPossibleMod[i]) <= 2 * ms2Tolerance) {
                        return "c" + i + modifiedAAMap.get(mass);
                    }
                }
            }
        }

        return null;
    }

    private TreeMap<Double, Double> addVirtualPeaks(double precursorMass, TreeMap<Double, Double> plMap) {
        double totalMass = precursorMass + 2 * MassTool.PROTON;
        TreeMap<Double, Double> finalPlMap = new TreeMap<>();
        for (double mz : plMap.keySet()) {
            finalPlMap.put(mz, plMap.get(mz));
        }
        for (double mz : plMap.keySet()) {
            double anotherMz = totalMass - mz;
            double leftMz = anotherMz - ms2Tolerance;
            double rightMz = anotherMz + ms2Tolerance;
            NavigableMap<Double, Double> temp = null;
            try {
                temp = plMap.subMap(leftMz, true, rightMz, true);
            } catch (IllegalArgumentException ex) {}

            if ((temp == null) || (temp.isEmpty())) {
                finalPlMap.put(anotherMz, plMap.get(mz));
            }
        }

        // Add two virtual peak. Because we have convert all y-ions to b-ions.
        finalPlMap.put(MassTool.PROTON, 1d);
        double cTermMz = precursorMass - massTool.H2O + MassTool.PROTON;
        double leftMz = cTermMz - ms2Tolerance;
        double rightMz = cTermMz + ms2Tolerance;
        NavigableMap<Double, Double> temp = null;
        try {
            temp = plMap.subMap(leftMz, true, rightMz, true);
        } catch (IllegalArgumentException ex) {}
        if ((temp == null) || (temp.isEmpty())) {
            finalPlMap.put(cTermMz, 1d);
        }

        return finalPlMap;
    }
}
