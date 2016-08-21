package proteomics.Index;

import java.util.*;
import proteomics.TheoSeq.*;

public class BuildIndex {

    private static final int maxLoopTime = 10;

    private float minPrecursorMass = 0;
    private float maxPrecursorMass = 0;
    private final MassTool massToolObj;
    private Map<String, Float> massTable;
    private Map<String, String> proPeptideMap = new HashMap<>();
    private Map<String, Float> peptideMassMap = new HashMap<>();
    private Map<String, Set<String>> peptideProMap = new HashMap<>();
    private Set<String> forCheckDuplicate = new HashSet<>();
    private Map<String, Float> decoyPeptideMassMap = new HashMap<>();
    private Map<String, String> decoyPeptideProMap = new HashMap<>();
    private Map<String, Float> fixModMap = new HashMap<>();
    private TreeMap<Float, Set<String>> massPeptideMap = new TreeMap<>();

    /////////////////////////////////public methods//////////////////////////////////////////////////////////////////
    public BuildIndex(Map<String, String> parameterMap) {
        // initialize parameters
        minPrecursorMass = Float.valueOf(parameterMap.get("min_precursor_mass"));
        maxPrecursorMass = Float.valueOf(parameterMap.get("max_precursor_mass"));
        String dbPath = parameterMap.get("db");
        int missedCleavage = Integer.valueOf(parameterMap.get("missed_cleavage"));
        float ms2Tolerance = Float.valueOf(parameterMap.get("ms2_tolerance"));
        float oneMinusBinOffset = 1 - Float.valueOf(parameterMap.get("mz_bin_offset"));

        // Read fix modification
        fixModMap.put("G", Float.valueOf(parameterMap.get("G")));
        fixModMap.put("A", Float.valueOf(parameterMap.get("A")));
        fixModMap.put("S", Float.valueOf(parameterMap.get("S")));
        fixModMap.put("P", Float.valueOf(parameterMap.get("P")));
        fixModMap.put("V", Float.valueOf(parameterMap.get("V")));
        fixModMap.put("T", Float.valueOf(parameterMap.get("T")));
        fixModMap.put("C", Float.valueOf(parameterMap.get("C")));
        fixModMap.put("I", Float.valueOf(parameterMap.get("I")));
        fixModMap.put("L", Float.valueOf(parameterMap.get("L")));
        fixModMap.put("N", Float.valueOf(parameterMap.get("N")));
        fixModMap.put("D", Float.valueOf(parameterMap.get("D")));
        fixModMap.put("Q", Float.valueOf(parameterMap.get("Q")));
        fixModMap.put("K", Float.valueOf(parameterMap.get("K")));
        fixModMap.put("E", Float.valueOf(parameterMap.get("E")));
        fixModMap.put("M", Float.valueOf(parameterMap.get("M")));
        fixModMap.put("H", Float.valueOf(parameterMap.get("H")));
        fixModMap.put("F", Float.valueOf(parameterMap.get("F")));
        fixModMap.put("R", Float.valueOf(parameterMap.get("R")));
        fixModMap.put("Y", Float.valueOf(parameterMap.get("Y")));
        fixModMap.put("W", Float.valueOf(parameterMap.get("W")));
        fixModMap.put("U", Float.valueOf(parameterMap.get("U")));
        fixModMap.put("O", Float.valueOf(parameterMap.get("O")));
        fixModMap.put("n", Float.valueOf(parameterMap.get("n")));

        // read protein database
        DbTool dbToolObj = new DbTool(dbPath);
        proPeptideMap = dbToolObj.returnSeqMap();

        // define a new MassTool object
        massToolObj = new MassTool(missedCleavage, fixModMap, parameterMap.get("cleavage_site"), parameterMap.get("protection_site"), Integer.valueOf(parameterMap.get("cleavage_from_c_term")) == 1, ms2Tolerance, oneMinusBinOffset);
        massTable = massToolObj.returnMassTable();

        // build database
        buildPeptideMap();
        buildDecoyPepChainMap();
        buildMassPeptideMap();
    }

    /////////////////////////////////////public methods////////////////////////////////////////////////////////////////////
    public MassTool returnMassToolObj() {
        return massToolObj;
    }

    public Map<String, String> returnProPepMap() {
        return proPeptideMap;
    }

    public Map<String, Float> returnPepMassMap() {
        return peptideMassMap;
    }

    public Map<String, Set<String>> returnPepProMap() {
        return peptideProMap;
    }

    public Map<String, Float> returnDecoyPepMassMap() {
        return decoyPeptideMassMap;
    }

    public Map<String, String> returnDecoyPepProMap() {
        return decoyPeptideProMap;
    }

    public Map<String, Float> returnFixModMap() {
        return fixModMap;
    }

    public TreeMap<Float, Set<String>> getMassPeptideMap() {
        return massPeptideMap;
    }

    private void buildPeptideMap() {
        Set<String> proIdSet = proPeptideMap.keySet();
        for (String proId : proIdSet) {
            String proSeq = proPeptideMap.get(proId);
            Set<String> peptideSet = massToolObj.buildPeptideSet(proSeq);
            for (String peptide : peptideSet) {
                if (peptide.contains("B") || peptide.contains("J") || peptide.contains("X") || peptide.contains("Z")) {
                    continue;
                }

                float massTemp = massToolObj.calResidueMass(peptide) + massTable.get("n") + massTable.get("H2O"); // calMass just calculate the residue mass, so we should add a H2O
                if ((massTemp <= maxPrecursorMass) && (massTemp >= minPrecursorMass)) {
                    peptideMassMap.put(peptide, massTemp);

                    // Add the sequence to the check set for decoy duplicate check
                    String templateSeq = peptide.replace('L', 'I'); // "L" and "I" have the same mass.
                    forCheckDuplicate.add(templateSeq);

                    if (peptideProMap.containsKey(peptide)) {
                        Set<String> proList = peptideProMap.get(peptide);
                        proList.add(proId);
                        peptideProMap.put(peptide, proList);
                    } else {
                        Set<String> proList = new HashSet<>();
                        proList.add(proId);
                        peptideProMap.put(peptide, proList);
                    }
                }
            }
        }
    }

    private void buildDecoyPepChainMap() {
        Set<String> peptideSet = peptideProMap.keySet();
        for (String originalPeptide : peptideSet) {
            String decoyPeptide = shuffleSeq2(originalPeptide);
            if (decoyPeptide.isEmpty()) {
                continue;
            }

            float decoyMassTemp = peptideMassMap.get(originalPeptide);
            decoyPeptideMassMap.put(decoyPeptide, decoyMassTemp);
            String proId = peptideProMap.get(originalPeptide).iterator().next();
            String decoyProId = "DECOY_" + proId;
            decoyPeptideProMap.put(decoyPeptide, decoyProId);
        }
    }

    private void buildMassPeptideMap() {
        for (String peptide : peptideMassMap.keySet()) {
            float mass = peptideMassMap.get(peptide);
            if (massPeptideMap.containsKey(mass)) {
                massPeptideMap.get(mass).add(peptide);
            } else {
                Set<String> temp = new HashSet<>();
                temp.add(peptide);
                massPeptideMap.put(mass, temp);
            }
        }
    }

    private String reverseSeq(String seq) {
        String decoySeq;
        if ((seq.charAt(seq.length() - 1) == 'K') || (seq.charAt(seq.length() - 1) == 'R')) {
            StringBuilder sb = new StringBuilder(seq.substring(0, seq.length() - 1)).reverse();
            decoySeq = sb.toString() + seq.substring(seq.length() - 1);
        } else {
            StringBuilder sb = new StringBuilder(seq).reverse();
            decoySeq = sb.toString();
        }
        if (forCheckDuplicate.contains(decoySeq.replace('L', 'I'))) {
            return "";
        } else {
            return decoySeq;
        }
    }

    private String shuffleSeq(String seq) {
        List<Character> seqList = new ArrayList<>(seq.length() + 2);
        Random randomObj = new Random();
        if ((seq.charAt(seq.length() - 1) == 'K') || (seq.charAt(seq.length() - 1) == 'R')) {
            for (int i = 0; i < seq.length() - 1; ++i) {
                seqList.add(seq.charAt(i));
            }
            StringBuilder sb;
            int loopTime = maxLoopTime;
            do {
                --loopTime;
                Collections.shuffle(seqList, randomObj);
                sb = new StringBuilder(seq.length() + 2);
                for (char aa : seqList) {
                    sb.append(aa);
                }
                sb.append(seq.charAt(seq.length() - 1));
            } while ((forCheckDuplicate.contains(sb.toString().replace('L', 'I'))) && (loopTime > 0));
            return sb.toString();
        } else {
            for (int i = 0; i < seq.length(); ++i) {
                seqList.add(seq.charAt(i));
            }
            StringBuilder sb;
            int loopTime = maxLoopTime;
            do {
                --loopTime;
                Collections.shuffle(seqList, randomObj);
                sb = new StringBuilder(seq.length() + 2);
                for (char aa : seqList) {
                    sb.append(aa);
                }
            } while ((forCheckDuplicate.contains(sb.toString().replace('L', 'I'))) && (loopTime > 0));
            return sb.toString();
        }
    }

    private String shuffleSeq2(String seq) {
        if ((seq.charAt(seq.length() - 1) == 'K') || (seq.charAt(seq.length() - 1) == 'R')) {
            char[] tempArray = seq.substring(0, seq.length() - 1).toCharArray();
            int idx = 0;
            while (idx < tempArray.length - 1) {
                char temp = tempArray[idx];
                tempArray[idx] = tempArray[idx + 1];
                tempArray[idx + 1] = temp;
                idx += 2;
            }
            String decoySeq = String.valueOf(tempArray) + seq.substring(seq.length() - 1, seq.length());
            if (forCheckDuplicate.contains(decoySeq.replace('L', 'I'))) {
                return "";
            } else {
                return decoySeq;
            }
        } else {
            char[] tempArray = seq.toCharArray();
            int idx = 0;
            while (idx < tempArray.length - 1) {
                char temp = tempArray[idx];
                tempArray[idx] = tempArray[idx + 1];
                tempArray[idx + 1] = temp;
                idx += 2;
            }
            String decoySeq = String.valueOf(tempArray);
            if (forCheckDuplicate.contains(decoySeq.replace('L', 'I'))) {
                return "";
            } else {
                return decoySeq;
            }
        }
    }
}
