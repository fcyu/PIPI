package proteomics.Index;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Segment.InferenceSegment;
import proteomics.TheoSeq.*;
import proteomics.Types.Peptide0;
import proteomics.Types.SparseBooleanVector;

public class BuildIndex {

    private static final Logger logger = LoggerFactory.getLogger(BuildIndex.class);

    private final MassTool massToolObj;
    private Map<Character, Double> fixModMap = new HashMap<>(25, 1);
    private double minPeptideMass = 9999;
    private double maxPeptideMass = 0;
    private InferenceSegment inference3SegmentObj;
    private TreeMap<Double, Set<String>> massPeptideMap = new TreeMap<>();
    private Map<String, Peptide0> peptide0Map;
    private final String labelling;
    private final DbTool dbTool; // this one doesn't contain contaminant proteins.

    public BuildIndex(Map<String, String> parameterMap, String labelling, boolean needCoding, boolean needDecoy) throws Exception {
        // initialize parameters
        int minPeptideLength = Math.max(5, Integer.valueOf(parameterMap.get("min_peptide_length")));
        int maxPeptideLength = Integer.valueOf(parameterMap.get("max_peptide_length"));
        String dbPath = parameterMap.get("db");
        int missedCleavage = Integer.valueOf(parameterMap.get("missed_cleavage"));
        double ms2Tolerance = Double.valueOf(parameterMap.get("ms2_tolerance"));
        double oneMinusBinOffset = 1 - Double.valueOf(parameterMap.get("mz_bin_offset"));
        this.labelling = labelling;

        // Read fix modification
        fixModMap.put('G', Double.valueOf(parameterMap.get("G")));
        fixModMap.put('A', Double.valueOf(parameterMap.get("A")));
        fixModMap.put('S', Double.valueOf(parameterMap.get("S")));
        fixModMap.put('P', Double.valueOf(parameterMap.get("P")));
        fixModMap.put('V', Double.valueOf(parameterMap.get("V")));
        fixModMap.put('T', Double.valueOf(parameterMap.get("T")));
        fixModMap.put('C', Double.valueOf(parameterMap.get("C")));
        fixModMap.put('I', Double.valueOf(parameterMap.get("I")));
        fixModMap.put('L', Double.valueOf(parameterMap.get("L")));
        fixModMap.put('N', Double.valueOf(parameterMap.get("N")));
        fixModMap.put('D', Double.valueOf(parameterMap.get("D")));
        fixModMap.put('Q', Double.valueOf(parameterMap.get("Q")));
        fixModMap.put('K', Double.valueOf(parameterMap.get("K")));
        fixModMap.put('E', Double.valueOf(parameterMap.get("E")));
        fixModMap.put('M', Double.valueOf(parameterMap.get("M")));
        fixModMap.put('H', Double.valueOf(parameterMap.get("H")));
        fixModMap.put('F', Double.valueOf(parameterMap.get("F")));
        fixModMap.put('R', Double.valueOf(parameterMap.get("R")));
        fixModMap.put('Y', Double.valueOf(parameterMap.get("Y")));
        fixModMap.put('W', Double.valueOf(parameterMap.get("W")));
        fixModMap.put('U', Double.valueOf(parameterMap.get("U")));
        fixModMap.put('O', Double.valueOf(parameterMap.get("O")));
        fixModMap.put('n', Double.valueOf(parameterMap.get("n")));
        fixModMap.put('c', Double.valueOf(parameterMap.get("c")));

        // read protein database
        dbTool = new DbTool(dbPath, parameterMap.get("database_type"));
        DbTool contaminantsDb = new DbTool(null, "contaminants");
        Map<String, String> proteinPeptideMap = contaminantsDb.returnSeqMap();
        proteinPeptideMap.putAll(dbTool.returnSeqMap()); // using the target sequence to replace contaminant sequence if there is conflict.

        // define a new MassTool object
        massToolObj = new MassTool(missedCleavage, fixModMap, parameterMap.get("cleavage_site"), parameterMap.get("protection_site"), Integer.valueOf(parameterMap.get("cleavage_from_c_term")) == 1, ms2Tolerance, oneMinusBinOffset, labelling);

        // build database
        inference3SegmentObj = new InferenceSegment(massToolObj, ms2Tolerance, parameterMap, fixModMap);

        Set<String> forCheckDuplicate = new HashSet<>(500000);
        Multimap<String, String> peptideProteinMap = HashMultimap.create();
        Map<String, Double> peptideMassMap = new HashMap<>(500000);
        Map<String, String> targetDecoyProteinSequenceMap = new HashMap<>();
        for (String proId : proteinPeptideMap.keySet()) {
            String proSeq = proteinPeptideMap.get(proId);
            Set<String> peptideSet = massToolObj.buildPeptideSet(proSeq);
            if (proSeq.startsWith("M")) { // Since the digestion doesn't take much time, just digest the whole protein again for easy read.
                peptideSet.addAll(massToolObj.buildPeptideSet(proSeq.substring(1)));
            }

            for (String peptide : peptideSet) {
                if (peptide.contains("B") || peptide.contains("J") || peptide.contains("X") || peptide.contains("Z") || peptide.contains("*")) {
                    continue;
                }

                if ((peptide.length() - 2 <= maxPeptideLength) && (peptide.length() - 2 >= minPeptideLength)) { // caution: there are n and c in the sequence
                    if (!forCheckDuplicate.contains(peptide.replace('L', 'I'))) { // don't record duplicate peptide sequences
                        // Add the sequence to the check set for duplicate check
                        forCheckDuplicate.add(peptide.replace('L', 'I'));

                        double mass = massToolObj.calResidueMass(peptide) + massToolObj.H2O;
                        // recode min and max peptide mass
                        if (mass < minPeptideMass) {
                            minPeptideMass = mass;
                        }
                        if (mass > maxPeptideMass) {
                            maxPeptideMass = mass;
                        }

                        peptideMassMap.put(peptide, mass);
                        peptideProteinMap.put(peptide, proId);
                    }

                    // considering the case that the sequence has multiple proteins. In the above if clock, such a protein wasn't recorded.
                    if (peptideProteinMap.containsKey(peptide)) {
                        peptideProteinMap.put(peptide, proId);
                    }
                }
            }

            targetDecoyProteinSequenceMap.put(proId, proSeq);

            if (needDecoy) {
                // decoy sequence
                String decoyProSeq;
                if (proSeq.startsWith("M")) {
                    decoyProSeq = "M" + shuffleSeq(proSeq.substring(1), massToolObj.getDigestSitePattern());
                } else {
                    decoyProSeq = shuffleSeq(proSeq, massToolObj.getDigestSitePattern());
                }
                peptideSet = massToolObj.buildPeptideSet(decoyProSeq);
                if (decoyProSeq.startsWith("M")) { // Since the digestion doesn't take much time, just digest the whole protein again for easy read.
                    peptideSet.addAll(massToolObj.buildPeptideSet(decoyProSeq.substring(1)));
                }

                for (String peptide : peptideSet) {
                    if (peptide.contains("B") || peptide.contains("J") || peptide.contains("X") || peptide.contains("Z") || peptide.contains("*")) {
                        continue;
                    }

                    if ((peptide.length() - 2 <= maxPeptideLength) && (peptide.length() - 2 >= minPeptideLength)) { // caution: there are n and c in the sequence
                        if (!forCheckDuplicate.contains(peptide.replace('L', 'I'))) { // don't record duplicate peptide sequences
                            // Add the sequence to the check set for duplicate check
                            forCheckDuplicate.add(peptide.replace('L', 'I'));

                            double mass = massToolObj.calResidueMass(peptide) + massToolObj.H2O;
                            // recode min and max peptide mass
                            if (mass < minPeptideMass) {
                                minPeptideMass = mass;
                            }
                            if (mass > maxPeptideMass) {
                                maxPeptideMass = mass;
                            }

                            peptideMassMap.put(peptide, mass);
                            peptideProteinMap.put(peptide, "DECOY_" + proId);
                        }
                    }
                }
                targetDecoyProteinSequenceMap.put("DECOY_" + proId, decoyProSeq);
            }
        }

        if (needDecoy) {
            // writer concatenated fasta
            Map<String, String> proteinAnnotationMap = dbTool.returnAnnotateMap();
            proteinAnnotationMap.putAll(contaminantsDb.returnAnnotateMap());
            BufferedWriter writer = new BufferedWriter(new FileWriter(dbPath + ".TD.fasta"));
            for (String proId : targetDecoyProteinSequenceMap.keySet()) {
                writer.write(String.format(Locale.US, ">%s %s\n", proId, proteinAnnotationMap.getOrDefault(proId, "")));
                writer.write(targetDecoyProteinSequenceMap.get(proId) + "\n");
            }
            writer.close();
        }

        Map<String, Peptide0> tempMap = new HashMap<>();
        for (String peptide : peptideMassMap.keySet()) {
            SparseBooleanVector code = null;
            if (needCoding) {
                code = inference3SegmentObj.generateSegmentBooleanVector(DbTool.getSequenceOnly(peptide));
            }

            Character leftFlank = null;
            Character rightFlank = null;
            String peptideString = DbTool.getSequenceOnly(peptide);
            if (peptideProteinMap.containsKey(peptide)) {
                for (String proteinId : peptideProteinMap.get(peptide)) {
                    String proteinSequence = targetDecoyProteinSequenceMap.get(proteinId);
                    int startIdx = proteinSequence.indexOf(peptideString);
                    while (startIdx >= 0) {
                        if (startIdx == 0 || ((startIdx == 1 && proteinSequence.charAt(0) == 'M'))) { // considering first "M" being cut situation.
                            int tempIdx = startIdx + peptideString.length();
                            if (tempIdx < proteinSequence.length()) {
                                rightFlank = proteinSequence.charAt(tempIdx);
                                if ((parameterMap.get("cleavage_from_c_term").contentEquals("1") && !parameterMap.get("protection_site").contains(rightFlank.toString())) || (parameterMap.get("cleavage_from_c_term").contentEquals("0") && parameterMap.get("cleavage_site").contains(rightFlank.toString()))) {
                                    leftFlank = '-';
                                    break;
                                } else {
                                    rightFlank = null;
                                }
                            } else if (tempIdx == proteinSequence.length()) {
                                leftFlank = '-';
                                rightFlank = '-';
                                break;
                            } else {
                                logger.warn("The peptide {} is longer than its protein {}.", peptideString, proteinSequence);
                            }
                        } else if (startIdx == proteinSequence.length() - peptideString.length()) {
                            leftFlank = proteinSequence.charAt(startIdx - 1);
                            if ((parameterMap.get("cleavage_from_c_term").contentEquals("1") && parameterMap.get("cleavage_site").contains(leftFlank.toString())) || (parameterMap.get("cleavage_from_c_term").contentEquals("0") && !parameterMap.get("protection_site").contains(leftFlank.toString()))) {
                                rightFlank = '-';
                                break;
                            } else {
                                leftFlank = null;
                            }
                        } else {
                            leftFlank = proteinSequence.charAt(startIdx - 1);
                            rightFlank = proteinSequence.charAt(startIdx + peptideString.length());
                            if ((parameterMap.get("cleavage_from_c_term").contentEquals("1") && parameterMap.get("cleavage_site").contains(leftFlank.toString()) && !parameterMap.get("protection_site").contains(rightFlank.toString())) || (parameterMap.get("cleavage_from_c_term").contentEquals("0") && parameterMap.get("cleavage_site").contains(rightFlank.toString()) && !parameterMap.get("protection_site").contains(leftFlank.toString()))) {
                                break;
                            } else {
                                leftFlank = null;
                                rightFlank = null;
                            }
                        }
                        startIdx = proteinSequence.indexOf(peptideString, startIdx + 1);
                    }

                    if (leftFlank != null && rightFlank != null) {
                        break;
                    }
                }

                if (leftFlank != null && rightFlank != null) {
                    tempMap.put(peptide, new Peptide0(code, true, peptideProteinMap.get(peptide).toArray(new String[peptideProteinMap.get(peptide).size()]), leftFlank, rightFlank));

                    if (massPeptideMap.containsKey(peptideMassMap.get(peptide))) {
                        massPeptideMap.get(peptideMassMap.get(peptide)).add(peptide);
                    } else {
                        Set<String> tempSet = new HashSet<>();
                        tempSet.add(peptide);
                        massPeptideMap.put(peptideMassMap.get(peptide), tempSet);
                    }
                }
            }
        }
        peptide0Map = new HashMap<>(tempMap); // Since this map won't be changed any more, using this step to create a HashMap with the capacity exactly equals the actual size.
    }

    public DbTool getDbTool() {
        return dbTool;
    }

    public MassTool returnMassToolObj() {
        return massToolObj;
    }

    public double getMinPeptideMass() {
        return minPeptideMass;
    }

    public double getMaxPeptideMass() {
        return maxPeptideMass;
    }

    public Map<Character, Double> returnFixModMap() {
        return fixModMap;
    }

    public InferenceSegment getInference3SegmentObj() {
        return inference3SegmentObj;
    }

    public TreeMap<Double, Set<String>> getMassPeptideMap() {
        return massPeptideMap;
    }

    public Map<String, Peptide0> getPeptide0Map() {
        return peptide0Map;
    }

    public String getLabelling() {
        return labelling;
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

    private String shuffleSeq(String seq, Pattern digestSitePattern) {
        Set<Integer> cutSiteSet = new HashSet<>();
        Matcher matcher = digestSitePattern.matcher(seq);
        while (matcher.find()) {
            cutSiteSet.add(matcher.start());
        }
        char[] tempArray = seq.toCharArray();
        int idx = 0;
        while (idx < tempArray.length - 1) {
            if (!cutSiteSet.contains(idx) && !cutSiteSet.contains(idx + 1)) {
                char temp = tempArray[idx];
                tempArray[idx] = tempArray[idx + 1];
                tempArray[idx + 1] = temp;
                idx += 2;
            } else {
                ++idx;
            }
        }
        return String.valueOf(tempArray);
    }
}
