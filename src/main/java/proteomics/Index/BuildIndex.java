package proteomics.Index;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.*;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Segment.InferenceSegment;
import ProteomicsLibrary.*;
import ProteomicsLibrary.Types.*;
import proteomics.Types.Peptide0;

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
        Map<String, String> proteinPeptideMap = contaminantsDb.getProteinSequenceMap();
        proteinPeptideMap.putAll(dbTool.getProteinSequenceMap()); // using the target sequence to replace contaminant sequence if there is conflict.

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

            for (String peptide : peptideSet) {
                if (MassTool.containsNonAAAndNC(peptide)) {
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
                    } else if (peptideProteinMap.containsKey(peptide)) {
                        // Considering the case that the sequence has multiple proteins. In the above if block, such a protein ID wasn't recorded. If there are decoy IDs, replace it with the current target ID since the target ID has a higher priority.
                        Set<String> proteinSet = new HashSet<>(peptideProteinMap.get(peptide));
                        peptideProteinMap.get(peptide).clear();
                        for (String protein : proteinSet) {
                            if (!protein.startsWith("DECOY_")) {
                                peptideProteinMap.put(peptide, protein);
                            }
                        }
                        peptideProteinMap.put(peptide, proId);
                    }
                }
            }

            targetDecoyProteinSequenceMap.put(proId, proSeq);

            if (needDecoy) {
                // decoy sequence
                String decoyProSeq = DbTool.shuffleSeq(proSeq, parameterMap.get("cleavage_site"), parameterMap.get("protection_site"), Integer.valueOf(parameterMap.get("cleavage_from_c_term")) == 1);
                peptideSet = massToolObj.buildPeptideSet(decoyProSeq);

                for (String peptide : peptideSet) {
                    if (MassTool.containsNonAAAndNC(peptide)) {
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
            Map<String, String> proteinAnnotationMap = dbTool.getProteinAnnotateMap();
            proteinAnnotationMap.putAll(contaminantsDb.getProteinAnnotateMap());
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

            Character[] leftRightFlank = DbTool.getLeftRightFlank(peptide, peptideProteinMap, targetDecoyProteinSequenceMap, parameterMap.get("cleavage_site"), parameterMap.get("protection_site"), parameterMap.get("cleavage_from_c_term").contentEquals("1"));
            if (leftRightFlank != null) {
                tempMap.put(peptide, new Peptide0(code, isTarget(peptideProteinMap.get(peptide)), peptideProteinMap.get(peptide).toArray(new String[0]), leftRightFlank[0], leftRightFlank[1]));

                if (massPeptideMap.containsKey(peptideMassMap.get(peptide))) {
                    massPeptideMap.get(peptideMassMap.get(peptide)).add(peptide);
                } else {
                    Set<String> tempSet = new HashSet<>();
                    tempSet.add(peptide);
                    massPeptideMap.put(peptideMassMap.get(peptide), tempSet);
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

    private boolean isTarget(Collection<String> proteinIds) {
        for (String protein : proteinIds) {
            if (protein.startsWith("DECOY_")) {
                return false;
            }
        }
        return true;
    }
}
