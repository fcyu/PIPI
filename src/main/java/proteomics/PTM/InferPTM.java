package proteomics.PTM;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;

import java.io.*;
import java.util.*;

public class InferPTM {

    private static final Logger logger = LoggerFactory.getLogger(InferPTM.class);

    private final MassTool massTool;
    private final int maxMs2Charge;
    private final Map<Character, Float> fixModMap;
    private final float minPtmMass;
    private final float maxPtmMass;
    private final float ms2Tolerance;
    private final float ptmMassTolerance;

    private int highestAasPriority = 0;
    private int matchedPeakNum;
    private double score;

    private Map<Character, Map<VarModParam, Integer>> modMap;

    public InferPTM(MassTool massTool, int maxMs2Charge, Map<Character, Float> fixModMap, Set<VarModParam> varModParamSet, float minPtmMass, float maxPtmMass, float ms2Tolerance) {
        this.massTool = massTool;
        this.maxMs2Charge = maxMs2Charge;
        this.fixModMap = fixModMap;
        this.minPtmMass = minPtmMass;
        this.maxPtmMass = maxPtmMass;
        this.ms2Tolerance = ms2Tolerance;
        ptmMassTolerance = Math.min(2 * ms2Tolerance, 0.02f);

        Map<Character, Float> massTable = massTool.returnMassTable();

        // Building an amino acid substitution matrix.
        Map<Character, Map<VarModParam, Integer>> aasMap = buildAASMap(massTable);
        for (Map<VarModParam, Integer> tempMap : aasMap.values()) {
            for (int priority : tempMap.values()) {
                if (priority > highestAasPriority) {
                    highestAasPriority = priority;
                }
            }
        }

        // Reading Mod table...
        Map<Character, Map<VarModParam, Integer>> ptmMap = readModFile(); // Only those higher than highestAasPriority have a meaning.

        // update ptm table with the high priority mods in the parameter file.
        for (VarModParam varModParam : varModParamSet) {
            if (ptmMap.containsKey(varModParam.aa)) {
                ptmMap.get(varModParam.aa).remove(varModParam); // caution: the put() doesn't update the key.
                ptmMap.get(varModParam.aa).put(varModParam, varModParam.priority);
            } else {
                Map<VarModParam, Integer> tempMap = new HashMap<>();
                tempMap.put(varModParam, varModParam.priority);
                ptmMap.put(varModParam.aa, tempMap);
            }
        }

        // merge PTM map and AA substitution map
        modMap = new HashMap<>(aasMap);
        for (char aa : ptmMap.keySet()) {
            for (VarModParam modEntry : ptmMap.get(aa).keySet()) {
                if (modMap.containsKey(aa)) {
                    if (modMap.get(aa).containsKey(modEntry)) {
                        if (modEntry.priority > highestAasPriority) { // this mod's priority has a meanning and is higher than all the AA substitutions'.
                            modMap.get(aa).remove(modEntry); // caution: the put() doesn't update the key.
                            modMap.get(aa).put(modEntry, modEntry.priority);
                        }
                    } else {
                        modMap.get(aa).put(modEntry, modEntry.priority);
                    }
                } else {
                    Map<VarModParam, Integer> tempMap = new HashMap<>();
                    tempMap.put(modEntry, modEntry.priority);
                    modMap.put(aa, tempMap);
                }
            }
        }
    }

    public PeptidePTMPattern tryPTM(SparseVector expProcessedPL, TreeMap<Float, Float> plMap, float precursorMass, String ptmFreeSequence, boolean isDecoy, double normalizedCrossCorr, char leftFlank, char rightFlank, int globalRank, int localMaxMS2Charge, float localMS1ToleranceL, float localMS1ToleranceR) {
        float ptmFreeMass = massTool.calResidueMass(ptmFreeSequence) + MassTool.H2O;
        float deltaMass = precursorMass - ptmFreeMass;
        float leftMassBound = deltaMass + localMS1ToleranceL;
        float rightMassBound = deltaMass + localMS1ToleranceR;
        Set<Integer> fixModIdxes = getFixModIdxes(ptmFreeSequence, fixModMap);

        // reset values;
        matchedPeakNum = 0;
        score = 0;

        PeptidePTMPattern peptidePTMPattern = new PeptidePTMPattern(ptmFreeSequence);

        // try different PTM combinations
        Set<String> checkedPtmPattern1 = new HashSet<>();
        Set<String> checkedPtmPattern2 = new HashSet<>();
        Set<String> checkedPtmPattern3 = new HashSet<>();
        Set<String> checkedPtmPattern4 = new HashSet<>();
        Set<String> checkedPtmPattern5 = new HashSet<>();

        for (int priorityT : new int[]{highestAasPriority + 1, 0}) { // common PTM first, then all PTMs.
            Map<Integer, Map<VarModParam, Integer>> idxVarModMap = getIdxVarModMap(ptmFreeSequence, fixModIdxes, modMap, priorityT);

            try1PTMs(idxVarModMap, leftMassBound, rightMassBound, ptmFreeSequence, isDecoy, normalizedCrossCorr, leftFlank, rightFlank, globalRank, checkedPtmPattern1, peptidePTMPattern, expProcessedPL, plMap, localMaxMS2Charge);

            if (idxVarModMap.size() > 1) {
                // Try 2 PTMs
                try2PTMs(idxVarModMap, leftMassBound, rightMassBound, ptmFreeSequence, isDecoy, normalizedCrossCorr, leftFlank, rightFlank, globalRank, checkedPtmPattern2, peptidePTMPattern, expProcessedPL, plMap, localMaxMS2Charge);
            }

            if (idxVarModMap.size() > 2) {
                // Try 3 PTMs
                try3PTMs(idxVarModMap, leftMassBound, rightMassBound, ptmFreeSequence, isDecoy, normalizedCrossCorr, leftFlank, rightFlank, globalRank, checkedPtmPattern3, peptidePTMPattern, expProcessedPL, plMap, localMaxMS2Charge);
            }

            if (priorityT > highestAasPriority) {
                if (idxVarModMap.size() > 3) {
                    // Try 4 PTMs
                    try4PTMs(idxVarModMap, leftMassBound, rightMassBound, ptmFreeSequence, isDecoy, normalizedCrossCorr, leftFlank, rightFlank, globalRank, checkedPtmPattern4, peptidePTMPattern, expProcessedPL, plMap, localMaxMS2Charge);
                }

                if (idxVarModMap.size() > 4) {
                    // Try 5 PTMs
                    try5PTMs(idxVarModMap, leftMassBound, rightMassBound, ptmFreeSequence, isDecoy, normalizedCrossCorr, leftFlank, rightFlank, globalRank, checkedPtmPattern5, peptidePTMPattern, expProcessedPL, plMap, localMaxMS2Charge);
                }
            }
        }

        return peptidePTMPattern;
    }

    private Map<Character, Map<VarModParam, Integer>> readModFile() {
        Map<Character, Map<VarModParam, Integer>> siteModMap = new HashMap<>();
        InputStream inputStream = getClass().getClassLoader().getResourceAsStream("modTable.tsv");
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream))) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (!line.isEmpty() && !line.startsWith("site")) {
                    String[] parts = line.split("\t");
                    char site = parts[0].charAt(0);
                    if (parts[0].trim().contentEquals("N-term")) {
                        site = 'n';
                    } else if (parts[0].trim().contentEquals("C-term")) {
                        site = 'c';
                    }
                    float mass = Float.valueOf(parts[2]);
                    if (mass >= minPtmMass && mass <= maxPtmMass) {
                        int priority = Integer.valueOf(parts[3]);
                        VarModParam temp = new VarModParam(mass, site, priority);
                        if (siteModMap.containsKey(site)) {
                            if (siteModMap.get(site).containsKey(temp)) {
                                if (siteModMap.get(site).get(temp) < priority) {
                                    siteModMap.get(site).remove(temp); // caution: the put() doesn't update the key.
                                    siteModMap.get(site).put(temp, priority);
                                }
                            } else {
                                siteModMap.get(site).put(temp, priority);
                            }
                        } else {
                            Map<VarModParam, Integer> tempMap = new HashMap<>();
                            tempMap.put(temp, priority);
                            siteModMap.put(site, tempMap);
                        }
                    }
                }
            }
        } catch (IOException ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }
        return siteModMap;
    }

    private Map<Character, Map<VarModParam, Integer>> buildAASMap(Map<Character, Float> massTable) {
        int[][] pam1Matrix = new int[][]{
                {9867 , 2    , 9    , 10   , 3    , 8    , 17   , 21   , 2    , 6    , 4    , 2    , 6    , 2    , 22   , 35   , 32   , 0    , 2    , 18},
                {1    , 9913 , 1    , 0    , 1    , 10   , 0    , 0    , 10   , 3    , 1    , 19   , 4    , 1    , 4    , 6    , 1    , 8    , 0    , 1},
                {4    , 1    , 9822 , 36   , 0    , 4    , 6    , 6    , 21   , 3    , 1    , 13   , 0    , 1    , 2    , 20   , 9    , 1    , 4    , 1},
                {6    , 0    , 42   , 9859 , 0    , 6    , 53   , 6    , 4    , 1    , 0    , 3    , 0    , 0    , 1    , 5    , 3    , 0    , 0    , 1},
                {1    , 1    , 0    , 0    , 9973 , 0    , 0    , 0    , 1    , 1    , 0    , 0    , 0    , 0    , 1    , 5    , 1    , 0    , 3    , 2},
                {3    , 9    , 4    , 5    , 0    , 9876 , 27   , 1    , 23   , 1    , 3    , 6    , 4    , 0    , 6    , 2    , 2    , 0    , 0    , 1},
                {10   , 0    , 7    , 56   , 0    , 35   , 9865 , 4    , 2    , 3    , 1    , 4    , 1    , 0    , 3    , 4    , 2    , 0    , 1    , 2},
                {21   , 1    , 12   , 11   , 1    , 3    , 7    , 9935 , 1    , 0    , 1    , 2    , 1    , 1    , 3    , 21   , 3    , 0    , 0    , 5},
                {1    , 8    , 18   , 3    , 1    , 20   , 1    , 0    , 9912 , 0    , 1    , 1    , 0    , 2    , 3    , 1    , 1    , 1    , 4    , 1},
                {2    , 2    , 3    , 1    , 2    , 1    , 2    , 0    , 0    , 9872 , 9    , 2    , 12   , 7    , 0    , 1    , 7    , 0    , 1    , 33},
                {3    , 1    , 3    , 0    , 0    , 6    , 1    , 1    , 4    , 22   , 9947 , 2    , 45   , 13   , 3    , 1    , 3    , 4    , 2    , 15},
                {2    , 37   , 25   , 6    , 0    , 12   , 7    , 2    , 2    , 4    , 1    , 9926 , 20   , 0    , 3    , 8    , 11   , 0    , 1    , 1},
                {1    , 1    , 0    , 0    , 0    , 2    , 0    , 0    , 0    , 5    , 8    , 4    , 9874 , 1    , 0    , 1    , 2    , 0    , 0    , 4},
                {1    , 1    , 1    , 0    , 0    , 0    , 0    , 1    , 2    , 8    , 6    , 0    , 4    , 9946 , 0    , 2    , 1    , 3    , 28   , 0},
                {13   , 5    , 2    , 1    , 1    , 8    , 3    , 2    , 5    , 1    , 2    , 2    , 1    , 1    , 9926 , 12   , 4    , 0    , 0    , 2},
                {28   , 11   , 34   , 7    , 11   , 4    , 6    , 16   , 2    , 2    , 1    , 7    , 4    , 3    , 17   , 9840 , 38   , 5    , 2    , 2},
                {22   , 2    , 13   , 4    , 1    , 3    , 2    , 2    , 1    , 11   , 2    , 8    , 6    , 1    , 5    , 32   , 9871 , 0    , 2    , 9},
                {0    , 2    , 0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    , 1    , 0    , 1    , 0    , 9976 , 1    , 0},
                {1    , 0    , 3    , 0    , 3    , 0    , 1    , 0    , 4    , 1    , 1    , 0    , 0    , 21   , 0    , 1    , 1    , 2    , 9945 , 1},
                {13   , 2    , 1    , 1    , 3    , 2    , 2    , 3    , 3    , 57   , 11   , 1    , 17   , 1    , 3    , 2    , 10   , 0    , 2    , 9901},
        };
        char[] aaArray = new char[]{'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};
        Map<Character, Map<VarModParam, Integer>> aasMap = new HashMap<>();
        for (int i = 0; i < aaArray.length; ++i) {
            for (int j = 0; j < aaArray.length; ++j) {
                if (i != j) {
                    if (!(aaArray[i] == 'I' && aaArray[j] == 'L') && !(aaArray[i] == 'L' && aaArray[j] == 'I')) { // "I" and "L" have the same mass. don't consider such a AA substitution.
                        float deltaMass = massTable.get(aaArray[j]) - massTable.get(aaArray[i]);
                        if (deltaMass >= minPtmMass && deltaMass <= maxPtmMass) {
                            VarModParam temp = new VarModParam(deltaMass, aaArray[i], pam1Matrix[i][j]);
                            if (aasMap.containsKey(aaArray[i])) {
                                if (aasMap.get(aaArray[i]).containsKey(temp)) {
                                    if (aasMap.get(aaArray[i]).get(temp) < pam1Matrix[i][j]) {
                                        aasMap.get(aaArray[i]).put(temp, pam1Matrix[i][j]);
                                    }
                                } else {
                                    aasMap.get(aaArray[i]).put(temp, pam1Matrix[i][j]);
                                }
                            } else {
                                Map<VarModParam, Integer> tempMap = new HashMap<>();
                                tempMap.put(temp, pam1Matrix[i][j]);
                                aasMap.put(aaArray[i], tempMap);
                            }
                        }
                    }
                }
            }
        }
        return aasMap;
    }

    private int getMatchedPeakNum(TreeMap<Float, Float> plMap, int localMaxMs2Charge, float[][] ionMatrix) {
        int K1 = 0;
        for (int i = 0; i < localMaxMs2Charge * 2; ++i) {
            for (int j = 0; j < ionMatrix[0].length; ++j) {
                for (float mz : plMap.keySet()) {
                    if (Math.abs(mz - ionMatrix[i][j]) <= ms2Tolerance) {
                        ++K1;
                        break;
                    }
                }
            }
        }
        return K1;
    }

    private Set<Integer> getFixModIdxes(String ptmFreeSequence, Map<Character, Float> fixModMap) {
        Set<Integer> outputSet = new HashSet<>(ptmFreeSequence.length(), 1);
        char[] tempArray = ptmFreeSequence.toCharArray();
        for (int i = 0; i < tempArray.length; ++i) {
            if (Math.abs(fixModMap.get(tempArray[i])) > 0.1) {
                outputSet.add(i);
            }
        }
        return outputSet;
    }

    private Map<Integer, Map<VarModParam, Integer>> getIdxVarModMap(String ptmFreeSequence, Set<Integer> fixModIdxes, Map<Character, Map<VarModParam, Integer>> mods, int priorityT) {
        Map<Integer, Map<VarModParam, Integer>> idxVarModMap = new HashMap<>(ptmFreeSequence.length() + 1, 1);
        for (int i = 0; i < ptmFreeSequence.length(); ++i) {
            if (!fixModIdxes.contains(i)) {
                char aa = ptmFreeSequence.charAt(i);
                if (aa == 'n') {
                    if (mods.containsKey('n')) {
                        Map<VarModParam, Integer> tempMap = new HashMap<>();
                        for (VarModParam modEntry : mods.get('n').keySet()) {
                            if (modEntry.priority >= priorityT) {
                                tempMap.put(modEntry, modEntry.priority);
                            }
                        }
                        if (!tempMap.isEmpty()) {
                            idxVarModMap.put(0, tempMap);
                        }
                    }
                } else if (aa == 'c') {
                    if (mods.containsKey('v')) {
                        Map<VarModParam, Integer> tempMap = new HashMap<>();
                        for (VarModParam modEntry : mods.get('c').keySet()) {
                            if (modEntry.priority >= priorityT) {
                                tempMap.put(modEntry, modEntry.priority);
                            }
                        }
                        if (!tempMap.isEmpty()) {
                            idxVarModMap.put(ptmFreeSequence.length() - 1, tempMap);
                        }
                    }
                } else {
                    if (mods.containsKey(aa)) {
                        Map<VarModParam, Integer> tempMap = new HashMap<>();
                        for (VarModParam modEntry : mods.get(aa).keySet()) {
                            if (modEntry.priority >= priorityT) {
                                tempMap.put(modEntry, modEntry.priority);
                            }
                        }
                        if (!tempMap.isEmpty()) {
                            idxVarModMap.put(i, tempMap);
                        }
                    }
                }
            }
        }
        return idxVarModMap;
    }

    private void try1PTMs(Map<Integer, Map<VarModParam, Integer>> idxVarModMap, float leftMassBound, float rightMassBound, String ptmFreeSequence, boolean isDecoy, double normalizedCrossCorr, char leftFlank, char rightFlank, int globalRank, Set<String> checkedPtmPattern, PeptidePTMPattern peptidePTMPattern, SparseVector expProcessedPL, TreeMap<Float, Float> plMap, int localMaxMS2Charge) { // Sometimes, the precursor mass error may affects the digitized spectrum.
        Integer[] idxArray = idxVarModMap.keySet().toArray(new Integer[idxVarModMap.size()]);
        Arrays.sort(idxArray);
        for (int i = 0; i < idxArray.length - 1; ++i) {
            for (VarModParam modEntry1 : idxVarModMap.get(idxArray[i]).keySet()) {
                if (modEntry1.mass <= rightMassBound && modEntry1.mass >= leftMassBound) {
                    if (!checkedPtmPattern.contains(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000))) {
                        checkedPtmPattern.add(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000));
                        PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap((ptmFreeSequence.length()));
                        positionDeltaMassMap.put(new Coordinate(idxArray[i], idxArray[i] + 1), modEntry1.mass);
                        Peptide peptideObj = new Peptide(ptmFreeSequence, isDecoy, massTool, maxMs2Charge, normalizedCrossCorr, leftFlank, rightFlank, globalRank);
                        peptideObj.setVarPTM(positionDeltaMassMap);
                        peptideObj.setScore(massTool.buildVector(peptideObj.getIonMatrix(), localMaxMS2Charge + 1).fastDot(expProcessedPL) * 0.25);
                        if (peptideObj.getScore() > score) {
                            score = peptideObj.getScore();
                            matchedPeakNum = getMatchedPeakNum(plMap, localMaxMS2Charge, peptideObj.getIonMatrix());
                            peptideObj.setMatchedPeakNum(matchedPeakNum);
                            peptidePTMPattern.update(peptideObj);
                        } else if (peptideObj.getScore() == score) {
                            int localMatchedPeakNum = getMatchedPeakNum(plMap, localMaxMS2Charge, peptideObj.getIonMatrix());
                            if (localMatchedPeakNum > matchedPeakNum) {
                                matchedPeakNum = localMatchedPeakNum;
                                peptideObj.setMatchedPeakNum(matchedPeakNum);
                                peptidePTMPattern.update(peptideObj);
                            }
                        }
                    }
                }
            }
        }
    }

    private void try2PTMs(Map<Integer, Map<VarModParam, Integer>> idxVarModMap, float leftMassBound, float rightMassBound, String ptmFreeSequence, boolean isDecoy, double normalizedCrossCorr, char leftFlank, char rightFlank, int globalRank, Set<String> checkedPtmPattern, PeptidePTMPattern peptidePTMPattern, SparseVector expProcessedPL, TreeMap<Float, Float> plMap, int localMaxMS2Charge) {
        Integer[] idxArray = idxVarModMap.keySet().toArray(new Integer[idxVarModMap.size()]);
        Arrays.sort(idxArray);
        for (int i = 0; i < idxArray.length - 1; ++i) {
            for (VarModParam modEntry1 : idxVarModMap.get(idxArray[i]).keySet()) {
                for (int j = i + 1; j < idxArray.length; ++j) {
                    for (VarModParam modEntry2 : idxVarModMap.get(idxArray[j]).keySet()) {
                        if (modEntry1.mass + modEntry2.mass <= rightMassBound && modEntry1.mass + modEntry2.mass >= leftMassBound) {
                            if (!checkedPtmPattern.contains(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEntry2.mass * 1000))) {
                                checkedPtmPattern.add(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEntry2.mass * 1000));
                                PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap((ptmFreeSequence.length()));
                                positionDeltaMassMap.put(new Coordinate(idxArray[i], idxArray[i] + 1), modEntry1.mass);
                                positionDeltaMassMap.put(new Coordinate(idxArray[j], idxArray[j] + 1), modEntry2.mass);
                                Peptide peptideObj = new Peptide(ptmFreeSequence, isDecoy, massTool, maxMs2Charge, normalizedCrossCorr, leftFlank, rightFlank, globalRank);
                                peptideObj.setVarPTM(positionDeltaMassMap);
                                peptideObj.setScore(massTool.buildVector(peptideObj.getIonMatrix(), localMaxMS2Charge + 1).fastDot(expProcessedPL) * 0.25);
                                if (peptideObj.getScore() > score) {
                                    score = peptideObj.getScore();
                                    matchedPeakNum = getMatchedPeakNum(plMap, localMaxMS2Charge, peptideObj.getIonMatrix());
                                    peptideObj.setMatchedPeakNum(matchedPeakNum);
                                    peptidePTMPattern.update(peptideObj);
                                } else if (peptideObj.getScore() == score) {
                                    int localMatchedPeakNum = getMatchedPeakNum(plMap, localMaxMS2Charge, peptideObj.getIonMatrix());
                                    if (localMatchedPeakNum > matchedPeakNum) {
                                        matchedPeakNum = localMatchedPeakNum;
                                        peptideObj.setMatchedPeakNum(matchedPeakNum);
                                        peptidePTMPattern.update(peptideObj);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    private void try3PTMs(Map<Integer, Map<VarModParam, Integer>> idxVarModMap, float leftMassBound, float rightMassBound, String ptmFreeSequence, boolean isDecoy, double normalizedCrossCorr, char leftFlank, char rightFlank, int globalRank, Set<String> checkedPtmPattern, PeptidePTMPattern peptidePTMPattern, SparseVector expProcessedPL, TreeMap<Float, Float> plMap, int localMaxMS2Charge) {
        Integer[] idxArray = idxVarModMap.keySet().toArray(new Integer[idxVarModMap.size()]);
        Arrays.sort(idxArray);
        for (int i = 0; i < idxArray.length - 2; ++i) {
            for (VarModParam modEntry1 : idxVarModMap.get(idxArray[i]).keySet()) {
                for (int j = i + 1; j < idxArray.length - 1; ++j) {
                    for (VarModParam modEntry2 : idxVarModMap.get(idxArray[j]).keySet()) {
                        if (Math.abs(modEntry1.mass + modEntry2.mass) >= ptmMassTolerance) { // two self cancelled PTM masses are not allowed.
                            for (int k = j + 1; k < idxArray.length; ++k) {
                                for (VarModParam modEntry3 : idxVarModMap.get(idxArray[k]).keySet()) {
                                    if (Math.abs(modEntry1.mass + modEntry3.mass) >= ptmMassTolerance && Math.abs(modEntry2.mass + modEntry3.mass) >= ptmMassTolerance) {
                                        if (modEntry1.mass + modEntry2.mass + modEntry3.mass <= rightMassBound && modEntry1.mass + modEntry2.mass + modEntry3.mass >= leftMassBound) {
                                            if (!checkedPtmPattern.contains(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEntry2.mass * 1000) + "-" + idxArray[k] + "-" + Math.round(modEntry3.mass * 1000))) {
                                                checkedPtmPattern.add(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEntry2.mass * 1000) + "-" + idxArray[k] + "-" + Math.round(modEntry3.mass * 1000));
                                                PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap((ptmFreeSequence.length()));
                                                positionDeltaMassMap.put(new Coordinate(idxArray[i], idxArray[i] + 1), modEntry1.mass);
                                                positionDeltaMassMap.put(new Coordinate(idxArray[j], idxArray[j] + 1), modEntry2.mass);
                                                positionDeltaMassMap.put(new Coordinate(idxArray[k], idxArray[k] + 1), modEntry3.mass);
                                                Peptide peptideObj = new Peptide(ptmFreeSequence, isDecoy, massTool, maxMs2Charge, normalizedCrossCorr, leftFlank, rightFlank, globalRank);
                                                peptideObj.setVarPTM(positionDeltaMassMap);
                                                peptideObj.setScore(massTool.buildVector(peptideObj.getIonMatrix(), localMaxMS2Charge + 1).fastDot(expProcessedPL) * 0.25);
                                                if (peptideObj.getScore() > score) {
                                                    score = peptideObj.getScore();
                                                    matchedPeakNum = getMatchedPeakNum(plMap, localMaxMS2Charge, peptideObj.getIonMatrix());
                                                    peptideObj.setMatchedPeakNum(matchedPeakNum);
                                                    peptidePTMPattern.update(peptideObj);
                                                } else if (peptideObj.getScore() == score) {
                                                    int localMatchedPeakNum = getMatchedPeakNum(plMap, localMaxMS2Charge, peptideObj.getIonMatrix());
                                                    if (localMatchedPeakNum > matchedPeakNum) {
                                                        matchedPeakNum = localMatchedPeakNum;
                                                        peptideObj.setMatchedPeakNum(matchedPeakNum);
                                                        peptidePTMPattern.update(peptideObj);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    private void try4PTMs(Map<Integer, Map<VarModParam, Integer>> idxVarModMap, float leftMassBound, float rightMassBound, String ptmFreeSequence, boolean isDecoy, double normalizedCrossCorr, char leftFlank, char rightFlank, int globalRank, Set<String> checkedPtmPattern, PeptidePTMPattern peptidePTMPattern, SparseVector expProcessedPL, TreeMap<Float, Float> plMap, int localMaxMS2Charge) {
        Integer[] idxArray = idxVarModMap.keySet().toArray(new Integer[idxVarModMap.size()]);
        Arrays.sort(idxArray);
        for (int i = 0; i < idxArray.length - 3; ++i) {
            for (VarModParam modEntry1 : idxVarModMap.get(idxArray[i]).keySet()) {
                for (int j = i + 1; j < idxArray.length - 2; ++j) {
                    for (VarModParam modEntry2 : idxVarModMap.get(idxArray[j]).keySet()) {
                        if (Math.abs(modEntry1.mass + modEntry2.mass) >= ptmMassTolerance) {
                            for (int k = j + 1; k < idxArray.length - 1; ++k) {
                                for (VarModParam modEntry3 : idxVarModMap.get(idxArray[k]).keySet()) {
                                    if (Math.abs(modEntry1.mass + modEntry3.mass) >= ptmMassTolerance && Math.abs(modEntry2.mass + modEntry3.mass) >= ptmMassTolerance) {
                                        for (int l = k + 1; l < idxArray.length; ++l) {
                                            for (VarModParam modEntry4 : idxVarModMap.get(idxArray[l]).keySet()) {
                                                if (Math.abs(modEntry1.mass + modEntry4.mass) >= ptmMassTolerance && Math.abs(modEntry2.mass + modEntry4.mass) >= ptmMassTolerance && Math.abs(modEntry3.mass + modEntry4.mass) >= ptmMassTolerance) {
                                                    if (modEntry1.mass + modEntry2.mass + modEntry3.mass + modEntry4.mass <= rightMassBound && modEntry1.mass + modEntry2.mass + modEntry3.mass + modEntry4.mass >= leftMassBound) {
                                                        if (!checkedPtmPattern.contains(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEntry2.mass * 1000) + "-" + idxArray[k] + "-" + Math.round(modEntry3.mass * 1000) + "-" + idxArray[l] + "-" + Math.round(modEntry4.mass * 1000))) {
                                                            checkedPtmPattern.add(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEntry2.mass * 1000) + "-" + idxArray[k] + "-" + Math.round(modEntry3.mass * 1000) + "-" + idxArray[l] + "-" + Math.round(modEntry4.mass * 1000));
                                                            PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap((ptmFreeSequence.length()));
                                                            positionDeltaMassMap.put(new Coordinate(idxArray[i], idxArray[i] + 1), modEntry1.mass);
                                                            positionDeltaMassMap.put(new Coordinate(idxArray[j], idxArray[j] + 1), modEntry2.mass);
                                                            positionDeltaMassMap.put(new Coordinate(idxArray[k], idxArray[k] + 1), modEntry3.mass);
                                                            positionDeltaMassMap.put(new Coordinate(idxArray[l], idxArray[l] + 1), modEntry4.mass);
                                                            Peptide peptideObj = new Peptide(ptmFreeSequence, isDecoy, massTool, maxMs2Charge, normalizedCrossCorr, leftFlank, rightFlank, globalRank);
                                                            peptideObj.setVarPTM(positionDeltaMassMap);
                                                            peptideObj.setScore(massTool.buildVector(peptideObj.getIonMatrix(), localMaxMS2Charge + 1).fastDot(expProcessedPL) * 0.25);
                                                            if (peptideObj.getScore() > score) {
                                                                score = peptideObj.getScore();
                                                                matchedPeakNum = getMatchedPeakNum(plMap, localMaxMS2Charge, peptideObj.getIonMatrix());
                                                                peptideObj.setMatchedPeakNum(matchedPeakNum);
                                                                peptidePTMPattern.update(peptideObj);
                                                            } else if (peptideObj.getScore() == score) {
                                                                int localMatchedPeakNum = getMatchedPeakNum(plMap, localMaxMS2Charge, peptideObj.getIonMatrix());
                                                                if (localMatchedPeakNum > matchedPeakNum) {
                                                                    matchedPeakNum = localMatchedPeakNum;
                                                                    peptideObj.setMatchedPeakNum(matchedPeakNum);
                                                                    peptidePTMPattern.update(peptideObj);
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    private void try5PTMs(Map<Integer, Map<VarModParam, Integer>> idxVarModMap, float leftMassBound, float rightMassBound, String ptmFreeSequence, boolean isDecoy, double normalizedCrossCorr, char leftFlank, char rightFlank, int globalRank, Set<String> checkedPtmPattern, PeptidePTMPattern peptidePTMPattern, SparseVector expProcessedPL, TreeMap<Float, Float> plMap, int localMaxMS2Charge) {
        Integer[] idxArray = idxVarModMap.keySet().toArray(new Integer[idxVarModMap.size()]);
        Arrays.sort(idxArray);
        for (int i = 0; i < idxArray.length - 4; ++i) {
            for (VarModParam modEntry1 : idxVarModMap.get(idxArray[i]).keySet()) {
                for (int j = i + 1; j < idxArray.length - 3; ++j) {
                    for (VarModParam modEntry2 : idxVarModMap.get(idxArray[j]).keySet()) {
                        if (Math.abs(modEntry1.mass + modEntry2.mass) >= ptmMassTolerance) {
                            for (int k = j + 1; k < idxArray.length - 2; ++k) {
                                for (VarModParam modEntry3 : idxVarModMap.get(idxArray[k]).keySet()) {
                                    if (Math.abs(modEntry1.mass + modEntry3.mass) >= ptmMassTolerance && Math.abs(modEntry2.mass + modEntry3.mass) >= ptmMassTolerance) {
                                        for (int l = k + 1; l < idxArray.length - 1; ++l) {
                                            for (VarModParam modEntry4 : idxVarModMap.get(idxArray[l]).keySet()) {
                                                if (Math.abs(modEntry1.mass + modEntry4.mass) >= ptmMassTolerance && Math.abs(modEntry2.mass + modEntry4.mass) >= ptmMassTolerance && Math.abs(modEntry3.mass + modEntry4.mass) >= ptmMassTolerance) {
                                                    for (int m = l + 1; m < idxArray.length; ++m) {
                                                        for (VarModParam modEntry5 : idxVarModMap.get(idxArray[m]).keySet()) {
                                                            if (Math.abs(modEntry1.mass + modEntry5.mass) >= ptmMassTolerance && Math.abs(modEntry2.mass + modEntry5.mass) >= ptmMassTolerance && Math.abs(modEntry3.mass + modEntry5.mass) >= ptmMassTolerance && Math.abs(modEntry4.mass + modEntry5.mass) >= ptmMassTolerance) {
                                                                if (modEntry1.mass + modEntry2.mass + modEntry3.mass + modEntry4.mass + modEntry5.mass <= rightMassBound && modEntry1.mass + modEntry2.mass + modEntry3.mass + modEntry4.mass + modEntry5.mass + modEntry5.mass >= leftMassBound) {
                                                                    if (!checkedPtmPattern.contains(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEntry2.mass * 1000) + "-" + idxArray[k] + "-" + Math.round(modEntry3.mass * 1000) + "-" + idxArray[l] + "-" + Math.round(modEntry4.mass * 1000) + "-" + idxArray[m] + "-" + Math.round(modEntry5.mass * 1000))) {
                                                                        checkedPtmPattern.add(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEntry2.mass * 1000) + "-" + idxArray[k] + "-" + Math.round(modEntry3.mass * 1000) + "-" + idxArray[l] + "-" + Math.round(modEntry4.mass * 1000) + "-" + idxArray[m] + "-" + Math.round(modEntry5.mass * 1000));
                                                                        PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap((ptmFreeSequence.length()));
                                                                        positionDeltaMassMap.put(new Coordinate(idxArray[i], idxArray[i] + 1), modEntry1.mass);
                                                                        positionDeltaMassMap.put(new Coordinate(idxArray[j], idxArray[j] + 1), modEntry2.mass);
                                                                        positionDeltaMassMap.put(new Coordinate(idxArray[k], idxArray[k] + 1), modEntry3.mass);
                                                                        positionDeltaMassMap.put(new Coordinate(idxArray[l], idxArray[l] + 1), modEntry4.mass);
                                                                        positionDeltaMassMap.put(new Coordinate(idxArray[m], idxArray[m] + 1), modEntry5.mass);
                                                                        Peptide peptideObj = new Peptide(ptmFreeSequence, isDecoy, massTool, maxMs2Charge, normalizedCrossCorr, leftFlank, rightFlank, globalRank);
                                                                        peptideObj.setVarPTM(positionDeltaMassMap);
                                                                        peptideObj.setScore(massTool.buildVector(peptideObj.getIonMatrix(), localMaxMS2Charge + 1).fastDot(expProcessedPL) * 0.25);
                                                                        if (peptideObj.getScore() > score) {
                                                                            score = peptideObj.getScore();
                                                                            matchedPeakNum = getMatchedPeakNum(plMap, localMaxMS2Charge, peptideObj.getIonMatrix());
                                                                            peptideObj.setMatchedPeakNum(matchedPeakNum);
                                                                            peptidePTMPattern.update(peptideObj);
                                                                        } else if (peptideObj.getScore() == score) {
                                                                            int localMatchedPeakNum = getMatchedPeakNum(plMap, localMaxMS2Charge, peptideObj.getIonMatrix());
                                                                            if (localMatchedPeakNum > matchedPeakNum) {
                                                                                matchedPeakNum = localMatchedPeakNum;
                                                                                peptideObj.setMatchedPeakNum(matchedPeakNum);
                                                                                peptidePTMPattern.update(peptideObj);
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
