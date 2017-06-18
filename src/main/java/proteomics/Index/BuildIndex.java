package proteomics.Index;

import java.sql.*;
import java.util.*;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Segment.InferenceSegment;
import proteomics.TheoSeq.*;
import proteomics.Types.SparseBooleanVector;

public class BuildIndex {

    private static final Logger logger = LoggerFactory.getLogger(BuildIndex.class);

    private float minPrecursorMass = 0;
    private float maxPrecursorMass = 0;
    private final MassTool massToolObj;
    private Map<Character, Float> fixModMap = new HashMap<>(25, 1);
    private float minPeptideMass = 9999;
    private float maxPeptideMass = 0;
    private InferenceSegment inference3SegmentObj;

    /////////////////////////////////public methods//////////////////////////////////////////////////////////////////
    public BuildIndex(Map<String, String> parameterMap, String sqlPath) {
        // initialize parameters
        minPrecursorMass = Float.valueOf(parameterMap.get("min_precursor_mass"));
        maxPrecursorMass = Float.valueOf(parameterMap.get("max_precursor_mass"));
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
            // prepare SQL database
            Connection sqlConnection = DriverManager.getConnection("jdbc:sqlite:" + sqlPath);
            Statement sqlStatement = sqlConnection.createStatement();
            sqlStatement.executeUpdate("DROP TABLE IF EXISTS peptideTable");
            sqlStatement.executeUpdate("CREATE TABLE peptideTable (peptideIndex INTEGER PRIMARY KEY ASC , peptideMass REAL NOT NULL, sequence TEXT NOT NULL, peptideCode TEXT NOT NULL, codeNormSquare REAL NOT NULL, isTarget INTEGER NOT NULL, proteins TEXT NOT NULL, leftFlank TEXT NOT NULL, rightFlank TEXT NOT NULL);");
            sqlStatement.executeUpdate("CREATE INDEX peptideMass ON peptideTable (peptideMass)");
            sqlStatement.executeUpdate("CREATE INDEX sequence ON peptideTable (sequence)");

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

                    float mass = massToolObj.calResidueMass(peptide) + MassTool.H2O;
                    if ((mass <= maxPrecursorMass) && (mass >= minPrecursorMass)) { // TODO: fix this bug
                        // Add the sequence to the check set for decoy duplicate check
                        forCheckDuplicate.add(peptide.replace('L', 'I')); // "L" and "I" have the same mass.

                        // recode min and max peptide mass
                        if (mass < minPeptideMass) {
                            minPeptideMass = mass;
                        }
                        if (mass > maxPeptideMass) {
                            maxPeptideMass = mass;
                        }

                        targetPeptideMassMap.put(peptide, mass);
                        if (targetPeptideProteinMap.containsKey(peptide)) {
                            Set<String> proteins = targetPeptideProteinMap.get(peptide);
                            proteins.add(proId);
                            targetPeptideProteinMap.put(peptide, proteins);
                        } else {
                            Set<String> proteins = new HashSet<>(10, 1);
                            proteins.add(proId);
                            targetPeptideProteinMap.put(peptide, proteins);
                        }
                    }
                }
            }

            PreparedStatement sqlPrepareStatement = sqlConnection.prepareStatement("INSERT INTO peptideTable (peptideMass, sequence, peptideCode, codeNormSquare, isTarget, proteins, leftFlank, rightFlank) VALUES (?, ?, ?, ?, ?, ?, ?, ?);");
           sqlConnection.setAutoCommit(false);

            for (String targetPeptide : targetPeptideMassMap.keySet()) {
                sqlPrepareStatement.setFloat(1, targetPeptideMassMap.get(targetPeptide));
                sqlPrepareStatement.setString(2, targetPeptide);

                SparseBooleanVector code = inference3SegmentObj.generateSegmentBooleanVector(inference3SegmentObj.cutTheoSegment(targetPeptide.substring(1, targetPeptide.length() - 1)));
                sqlPrepareStatement.setString(3, code.toString());
                sqlPrepareStatement.setDouble(4, code.norm2square());

                sqlPrepareStatement.setBoolean(5, true);

                StringBuilder sb = new StringBuilder(targetPeptideProteinMap.get(targetPeptide).size() * 10);
                for (String proteinId : targetPeptideProteinMap.get(targetPeptide)) {
                    sb.append(proteinId);
                    sb.append(";");
                }
                sqlPrepareStatement.setString(6, sb.toString());

                String leftFlank = "-";
                String rightFlank = "-";
                String peptideString = targetPeptide.substring(1, targetPeptide.length() - 1);
                if (targetPeptideProteinMap.containsKey(targetPeptide)) {
                    String proteinSequence = proteinPeptideMap.get(targetPeptideProteinMap.get(targetPeptide).iterator().next());
                    int startIdx = proteinSequence.indexOf(peptideString);
                    if (startIdx == -1) {
                        logger.warn("Cannot locate {} in protein {}.", targetPeptide, proteinSequence);
                    } else if (startIdx == 0) {
                        int tempIdx = peptideString.length();
                        if (tempIdx < proteinSequence.length()) {
                            rightFlank = proteinSequence.substring(tempIdx, tempIdx + 1);
                        }
                    } else if (startIdx == proteinSequence.length() - peptideString.length()) {
                        leftFlank = proteinSequence.substring(startIdx - 1, startIdx);
                    } else {
                        leftFlank = proteinSequence.substring(startIdx - 1, startIdx);
                        rightFlank = proteinSequence.substring(startIdx + peptideString.length(), startIdx + peptideString.length() + 1);
                    }
                }
                sqlPrepareStatement.setString(7, leftFlank);
                sqlPrepareStatement.setString(8, rightFlank);
                sqlPrepareStatement.executeUpdate();

                String decoyPeptide = shuffleSeq2(targetPeptide.substring(1, targetPeptide.length() - 1), forCheckDuplicate);
                if (!decoyPeptide.isEmpty()) {
                    decoyPeptide = "n" + decoyPeptide + "c";

                    sqlPrepareStatement.setFloat(1, targetPeptideMassMap.get(targetPeptide));
                    sqlPrepareStatement.setString(2, decoyPeptide);

                    code = inference3SegmentObj.generateSegmentBooleanVector(inference3SegmentObj.cutTheoSegment(decoyPeptide.substring(1, decoyPeptide.length() - 1)));
                    sqlPrepareStatement.setString(3, code.toString());
                    sqlPrepareStatement.setDouble(4, code.norm2square());

                    sqlPrepareStatement.setBoolean(5, false);
                    sb = new StringBuilder(targetPeptideProteinMap.get(targetPeptide).size() * 10);
                    for (String proteinId : targetPeptideProteinMap.get(targetPeptide)) {
                        sb.append("DECOY_");
                        sb.append(proteinId);
                        sb.append(";");
                    }
                    sqlPrepareStatement.setString(6, sb.toString());
                    sqlPrepareStatement.setString(7, leftFlank);
                    sqlPrepareStatement.setString(8, rightFlank);
                    sqlPrepareStatement.executeUpdate();
                }
            }
            sqlConnection.commit();
            sqlConnection.setAutoCommit(true);
            sqlStatement.close();
            sqlPrepareStatement.close();
            sqlConnection.close();
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
}
