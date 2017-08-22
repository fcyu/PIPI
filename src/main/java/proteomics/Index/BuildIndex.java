package proteomics.Index;

import java.util.*;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Segment.InferenceSegment;
import proteomics.TheoSeq.*;
import proteomics.Types.Peptide0;
import proteomics.Types.SparseBooleanVector;

public class BuildIndex {

    private static final Logger logger = LoggerFactory.getLogger(BuildIndex.class);

    private final MassTool massToolObj;
    private Map<Character, Float> fixModMap = new HashMap<>(25, 1);
    private float minPeptideMass = 9999;
    private float maxPeptideMass = 0;
    private InferenceSegment inference3SegmentObj;
    private TreeMap<Float, Set<String>> massPeptideMap = new TreeMap<>();
    private Map<String, Peptide0> peptide0Map;

    /////////////////////////////////public methods//////////////////////////////////////////////////////////////////
    public BuildIndex(Map<String, String> parameterMap) {
        // initialize parameters
        int minPeptideLength = Integer.valueOf(parameterMap.get("min_peptide_length"));
        int maxPeptideLength = Integer.valueOf(parameterMap.get("max_peptide_length"));
        String dbPath = parameterMap.get("db");
        int missedCleavage = Integer.valueOf(parameterMap.get("missed_cleavage"));
        float ms2Tolerance = Float.valueOf(parameterMap.get("ms2_tolerance"));
        float oneMinusBinOffset = 1 - Float.valueOf(parameterMap.get("mz_bin_offset"));

        // Read fix modification
        fixModMap.put('G', Float.valueOf(parameterMap.get("G")));
        fixModMap.put('A', Float.valueOf(parameterMap.get("A")));
        fixModMap.put('S', Float.valueOf(parameterMap.get("S")));
        fixModMap.put('P', Float.valueOf(parameterMap.get("P")));
        fixModMap.put('V', Float.valueOf(parameterMap.get("V")));
        fixModMap.put('T', Float.valueOf(parameterMap.get("T")));
        fixModMap.put('C', Float.valueOf(parameterMap.get("C")));
        fixModMap.put('I', Float.valueOf(parameterMap.get("I")));
        fixModMap.put('L', Float.valueOf(parameterMap.get("L")));
        fixModMap.put('N', Float.valueOf(parameterMap.get("N")));
        fixModMap.put('D', Float.valueOf(parameterMap.get("D")));
        fixModMap.put('Q', Float.valueOf(parameterMap.get("Q")));
        fixModMap.put('K', Float.valueOf(parameterMap.get("K")));
        fixModMap.put('E', Float.valueOf(parameterMap.get("E")));
        fixModMap.put('M', Float.valueOf(parameterMap.get("M")));
        fixModMap.put('H', Float.valueOf(parameterMap.get("H")));
        fixModMap.put('F', Float.valueOf(parameterMap.get("F")));
        fixModMap.put('R', Float.valueOf(parameterMap.get("R")));
        fixModMap.put('Y', Float.valueOf(parameterMap.get("Y")));
        fixModMap.put('W', Float.valueOf(parameterMap.get("W")));
        fixModMap.put('U', Float.valueOf(parameterMap.get("U")));
        fixModMap.put('O', Float.valueOf(parameterMap.get("O")));
        fixModMap.put('n', Float.valueOf(parameterMap.get("n")));
        fixModMap.put('c', Float.valueOf(parameterMap.get("c")));

        // read protein database
        DbTool dbToolObj = new DbTool(dbPath);
        Map<String, String> proteinPeptideMap = dbToolObj.returnSeqMap();

        // define a new MassTool object
        massToolObj = new MassTool(missedCleavage, fixModMap, parameterMap.get("cleavage_site"), parameterMap.get("protection_site"), Integer.valueOf(parameterMap.get("cleavage_from_c_term")) == 1, ms2Tolerance, oneMinusBinOffset);

        // build database
        try {
            inference3SegmentObj = new InferenceSegment(massToolObj.returnMassTable(), ms2Tolerance, parameterMap, fixModMap);

            Set<String> forCheckDuplicate = new HashSet<>(500000);
            Map<String, Set<String>> targetPeptideProteinMap = new HashMap<>(500000);
            Map<String, Float> targetPeptideMassMap = new HashMap<>(500000);
            for (String proId : proteinPeptideMap.keySet()) {
                String proSeq = proteinPeptideMap.get(proId);
                Set<String> peptideSet = massToolObj.buildPeptideSet(proSeq);
                for (String peptide : peptideSet) {
                    if (peptide.contains("B") || peptide.contains("J") || peptide.contains("X") || peptide.contains("Z") || peptide.contains("*")) {
                        continue;
                    }

                    if ((peptide.length() <= maxPeptideLength) && (peptide.length() >= minPeptideLength)) {
                        if (!forCheckDuplicate.contains(peptide.replace('L', 'I'))) { // don't record duplicate peptide sequences
                            // Add the sequence to the check set for duplicate check
                            forCheckDuplicate.add(peptide.replace('L', 'I'));

                            float mass = massToolObj.calResidueMass(peptide) + MassTool.H2O;
                            // recode min and max peptide mass
                            if (mass < minPeptideMass) {
                                minPeptideMass = mass;
                            }
                            if (mass > maxPeptideMass) {
                                maxPeptideMass = mass;
                            }

                            targetPeptideMassMap.put(peptide, mass);
                            Set<String> proteins = new HashSet<>(10, 1);
                            proteins.add(proId);
                            targetPeptideProteinMap.put(peptide, proteins);
                        }

                        // considering the case that the sequence has multiple proteins. In the above if clock, such a protein wasn't recorded.
                        if (targetPeptideProteinMap.containsKey(peptide)) {
                            targetPeptideProteinMap.get(peptide).add(proId);
                        }
                    }
                }
            }

            Map<String, Peptide0> tempMap = new HashMap<>();
            for (String targetPeptide : targetPeptideMassMap.keySet()) {
                SparseBooleanVector targetCode = inference3SegmentObj.generateSegmentBooleanVector(inference3SegmentObj.cutTheoSegment(targetPeptide.substring(1, targetPeptide.length() - 1)));

                char leftFlank = '-';
                char rightFlank = '-';
                String peptideString = targetPeptide.substring(1, targetPeptide.length() - 1);
                if (targetPeptideProteinMap.containsKey(targetPeptide)) {
                    String proteinSequence = proteinPeptideMap.get(targetPeptideProteinMap.get(targetPeptide).iterator().next());
                    int startIdx = proteinSequence.indexOf(peptideString);
                    if (startIdx == -1) {
                        logger.warn("Cannot locate {} in protein {}.", targetPeptide, proteinSequence);
                    } else if (startIdx == 0) {
                        int tempIdx = peptideString.length();
                        if (tempIdx < proteinSequence.length()) {
                            rightFlank = proteinSequence.charAt(tempIdx);
                        }
                    } else if (startIdx == proteinSequence.length() - peptideString.length()) {
                        leftFlank = proteinSequence.charAt(startIdx - 1);
                    } else {
                        leftFlank = proteinSequence.charAt(startIdx - 1);
                        rightFlank = proteinSequence.charAt(startIdx + peptideString.length());
                    }
                }

                tempMap.put(targetPeptide, new Peptide0(targetPeptide, targetPeptideMassMap.get(targetPeptide), targetCode, true, new HashSet<>(targetPeptideProteinMap.get(targetPeptide)), leftFlank, rightFlank));

                if (massPeptideMap.containsKey(targetPeptideMassMap.get(targetPeptide))) {
                    massPeptideMap.get(targetPeptideMassMap.get(targetPeptide)).add(targetPeptide);
                } else {
                    Set<String> tempSet = new HashSet<>();
                    tempSet.add(targetPeptide);
                    massPeptideMap.put(targetPeptideMassMap.get(targetPeptide), tempSet);
                }

                // decoy peptides
                String decoyPeptide = shuffleSeq(targetPeptide.substring(1, targetPeptide.length() - 1), forCheckDuplicate);
                if (!decoyPeptide.isEmpty()) {
                    decoyPeptide = "n" + decoyPeptide + "c";
                    forCheckDuplicate.add(decoyPeptide.replace('L', 'I'));
                    SparseBooleanVector decoyCode = inference3SegmentObj.generateSegmentBooleanVector(inference3SegmentObj.cutTheoSegment(decoyPeptide.substring(1, decoyPeptide.length() - 1)));

                    Set<String> decoyProteins = new HashSet<>(targetPeptideProteinMap.get(targetPeptide).size() + 1, 1);
                    for (String proteinId : targetPeptideProteinMap.get(targetPeptide)) {
                        decoyProteins.add("DECOY_" + proteinId);
                    }

                    tempMap.put(decoyPeptide, new Peptide0(decoyPeptide, targetPeptideMassMap.get(targetPeptide), decoyCode, false, decoyProteins, leftFlank, rightFlank));
                    if (massPeptideMap.containsKey(targetPeptideMassMap.get(targetPeptide))) {
                        massPeptideMap.get(targetPeptideMassMap.get(targetPeptide)).add(decoyPeptide);
                    } else {
                        Set<String> tempSet = new HashSet<>();
                        tempSet.add(decoyPeptide);
                        massPeptideMap.put(targetPeptideMassMap.get(targetPeptide), tempSet);
                    }
                }
            }
            peptide0Map = new HashMap<>(tempMap);
        } catch (Exception ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }
    }

    /////////////////////////////////////public methods////////////////////////////////////////////////////////////////////
    public MassTool returnMassToolObj() {
        return massToolObj;
    }

    public float getMinPeptideMass() {
        return minPeptideMass;
    }

    public float getMaxPeptideMass() {
        return maxPeptideMass;
    }

    public Map<Character, Float> returnFixModMap() {
        return fixModMap;
    }

    public InferenceSegment getInference3SegmentObj() {
        return inference3SegmentObj;
    }

    public TreeMap<Float, Set<String>> getMassPeptideMap() {
        return massPeptideMap;
    }

    public Map<String, Peptide0> getPeptide0Map() {
        return peptide0Map;
    }

    private String shuffleSeq2(String seq, Set<String> forCheckDuplicate) {
        Random random = new Random(0);
        char[] tempArray = seq.substring(0, seq.length() - 1).toCharArray();
        String decoySeq;
        int time = 0;
        do {
            // the standard Fisher-Yates shuffle
            for (int i = 0; i < tempArray.length; ++i) {
                int j = random.nextInt(tempArray.length);
                while (j == i) {
                    j = random.nextInt(tempArray.length);
                }
                char temp = tempArray[i];
                tempArray[i] = tempArray[j];
                tempArray[j] = temp;
            }
            decoySeq = String.valueOf(tempArray) + seq.substring(seq.length() - 1, seq.length());
            ++time;
        } while (forCheckDuplicate.contains("n" + decoySeq.replace('L', 'I') + "c") && (time < 10));
        if (forCheckDuplicate.contains("n" + decoySeq.replace('L', 'I') + "c")) {
            return "";
        } else {
            return decoySeq;
        }
    }

    private String shuffleSeq(String seq, Set<String> forCheckDuplicate) {
        char[] tempArray = seq.substring(0, seq.length() - 1).toCharArray();
        int idx = 0;
        while (idx < tempArray.length - 1) {
            char temp = tempArray[idx];
            tempArray[idx] = tempArray[idx + 1];
            tempArray[idx + 1] = temp;
            idx += 2;
        }
        String decoySeq = String.valueOf(tempArray) + seq.substring(seq.length() - 1, seq.length());
        if (forCheckDuplicate.contains("n" + decoySeq.replace('L', 'I') + "c")) {
            return "";
        } else {
            return decoySeq;
        }
    }
}
