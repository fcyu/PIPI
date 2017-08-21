package proteomics.PTM;


import org.apache.commons.math3.util.CombinatoricsUtils;
import proteomics.Search.CalScore;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;

import java.util.*;

public class GeneratePtmCandidatesCalculateScore {

    private static final int globalVarModMaxNum = 5; // This value cannot be larger than 5. Otherwise, change generateLocalIdxModMassMap accordingly.
    private static final int maxPermutationNum = 1000;

    private final SpectrumEntry spectrumEntry;
    private final MassTool massToolObj;
    private final int maxMs2Charge;
    private final float tolerance;
    private final Set<VarModParam> varModParamSet;
    private final Map<Character, Float> fixModMap;
    private final SparseVector expProcessedPL;
    private FinalResultEntry psm;
    private Set<Peptide> checkedSequenceSet = new HashSet<>(maxPermutationNum * 2 + 50, 1); // record checked sequence to avoid recording the same sequence twice
    private Map<String, TreeSet<PeptideScore>> modSequences = new HashMap<>(15, 1);

    public GeneratePtmCandidatesCalculateScore(SpectrumEntry spectrumEntry, MassTool massToolObj, Set<VarModParam> varModParamSet, Map<Character, Float> fixModMap, float ms2Tolerance, int maxMs2Charge, SparseVector expProcessedPL) {
        this.spectrumEntry = spectrumEntry;
        this.massToolObj = massToolObj;
        this.maxMs2Charge = maxMs2Charge;
        tolerance = Math.max(ms2Tolerance, 0.1f);
        this.varModParamSet = varModParamSet;
        this.fixModMap = fixModMap;
        this.expProcessedPL = expProcessedPL;
        psm = new FinalResultEntry(spectrumEntry.scanNum, spectrumEntry.precursorCharge, spectrumEntry.precursorMz, spectrumEntry.mgfTitle);
    }

    public Set<Peptide> eliminateMissedCleavageCausedPtms(List<Peptide> ptmFreeCandidates, List<Peptide> ptmOnlyCandidates) {
        // Additional short missed cleavaged amino acid sequence may be canceled out by negative PTM. In this situation, we eliminate the missed cleavaged one.
        List<Peptide> tempList = new LinkedList<>();
        tempList.addAll(ptmFreeCandidates);
        tempList.addAll(ptmOnlyCandidates);
        Peptide[] tempArray = tempList.toArray(new Peptide[tempList.size()]);
        Set<Peptide> candidates = new HashSet<>(11, 1);
        for (int i = 0; i < tempArray.length; ++i) {
            boolean keep = true;
            String tempStr1 = tempArray[i].getNormalizedPeptideString();
            tempStr1 = tempStr1.substring(1, tempStr1.length() - 1);
            for (int j = 0; j < tempArray.length; ++j) {
                if (i != j) {
                    String tempStr2 = tempArray[j].getNormalizedPeptideString();
                    tempStr2 = tempStr2.substring(1, tempStr2.length() - 1);
                    if (tempStr1.contains(tempStr2) && (tempArray[i].getVarPTMs() != null)) {
                        Map.Entry<Coordinate, Float> tempEntry = tempArray[i].getVarPTMs().firstEntry();
                        if ((tempEntry.getValue() < 0) && (tempEntry.getKey().y - tempEntry.getKey().x > 1)) {
                            keep = false;
                            break;
                        }
                        tempEntry = tempArray[i].getVarPTMs().lastEntry();
                        if ((tempEntry.getValue() < 0) && (tempEntry.getKey().y - tempEntry.getKey().x > 1)) {
                            keep = false;
                            break;
                        }
                    }
                }
            }
            if (keep) {
                candidates.add(tempArray[i]);
            }
        }
        return candidates;
    }

    public FinalResultEntry generateAllPtmCandidatesCalculateScore(Set<Peptide> candidates) {
        // Add all PTM sequence given known variable modifications.
        if (!candidates.isEmpty()) {
            for (Peptide candidate : candidates) {
                if (candidate.hasVarPTM()) {
                    Set<Integer> fixModIdxes = getFixModIdxes(candidate.getPTMFreeSeq(), fixModMap);

                    float ptmFreePrecursorMass = massToolObj.calResidueMass(candidate.getPTMFreeSeq()) + MassTool.H2O;
                    // having only one unknown modification
                    float deltaMass = spectrumEntry.precursorMass - ptmFreePrecursorMass;
                    if (Math.abs(deltaMass) >= tolerance) {
                        if (isNewPtmMass(varModParamSet, deltaMass, tolerance)) {
                            for (int i = 1; i < candidate.getPTMFreeSeq().length() - 1; ++i) { // Don't try N-term and C-term because they are the same as the first and the last amino acids in terms of producing spectrum.
                                if (!fixModIdxes.contains(i)) {
                                    PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap(candidate.getPTMFreeSeq().length());
                                    positionDeltaMassMap.put(new Coordinate(i, i + 1), deltaMass);
                                    Peptide peptideObj = new Peptide(candidate.getPTMFreeSeq(), candidate.isDecoy(), massToolObj, maxMs2Charge, candidate.getNormalizedCrossCorr(), candidate.getLeftFlank(), candidate.getRightFlank(), candidate.getGlobalRank());
                                    peptideObj.setVarPTM(positionDeltaMassMap);
                                    peptideObj.setUnknownPtmNum(1);
                                    if (!checkedSequenceSet.contains(peptideObj)) {
                                        CalScore.calScore(peptideObj, expProcessedPL, psm, massToolObj, modSequences);
                                        checkedSequenceSet.add(peptideObj);
                                    }
                                }
                            }
                        }
                    }

                    // having known modifications and maybe an unknown modification.
                    int permutationNum = 0;
                    // summarize all possible modified idxes.
                    String ptmFreeSequence = candidate.getPTMFreeSeq();
                    Map<Integer, List<Float>> idxVarModMassMap = new HashMap<>(ptmFreeSequence.length(), 1);
                    for (int i = 0; i < ptmFreeSequence.length(); ++i) {
                        char aa = ptmFreeSequence.charAt(i);
                        for (VarModParam varModParam : varModParamSet) {
                            if (varModParam.aa == aa) {
                                boolean record = true;
                                if (varModParam.proteinTerminal) {
                                    // check the protein terminal criteria
                                    if (((varModParam.aa == 'n') && (candidate.getLeftFlank() != '-')) || ((varModParam.aa == 'c') && (candidate.getRightFlank() != '-'))) {
                                        record = false;
                                    }
                                }
                                if (record) {
                                    if (idxVarModMassMap.containsKey(i)) {
                                        idxVarModMassMap.get(i).add(varModParam.modMass);
                                    } else {
                                        List<Float> temp = new LinkedList<>();
                                        temp.add(varModParam.modMass);
                                        idxVarModMassMap.put(i, temp);
                                    }
                                }
                            }
                        }
                    }
                    if (!idxVarModMassMap.isEmpty()) {
                        Integer[] allIdxArray = idxVarModMassMap.keySet().toArray(new Integer[idxVarModMassMap.size()]);
                        Arrays.sort(allIdxArray);
                        for (int i = 1; i <= Math.min(idxVarModMassMap.size(), globalVarModMaxNum); ++i) {
                            Iterator<int[]> tempIterator = CombinatoricsUtils.combinationsIterator(allIdxArray.length, i);
                            while (tempIterator.hasNext()) {
                                // generate idx combinations
                                int[] idxIdxArray = tempIterator.next();
                                int[] idxCombination = new int[idxIdxArray.length];
                                for (int j = 0; j < idxIdxArray.length; ++j) {
                                    idxCombination[j] = allIdxArray[idxIdxArray[j]];
                                }
                                Arrays.sort(idxCombination);

                                // generate idx-mass
                                if (idxCombination.length == 5) {
                                    for (int i0 = 0; i0 < idxVarModMassMap.get(idxCombination[0]).size(); ++i0) {
                                        for (int i1 = 0; i1 < idxVarModMassMap.get(idxCombination[1]).size(); ++i1) {
                                            for (int i2 = 0; i2 < idxVarModMassMap.get(idxCombination[2]).size(); ++i2) {
                                                for (int i3 = 0; i3 < idxVarModMassMap.get(idxCombination[3]).size(); ++i3) {
                                                    for (int i4 = 0; i4 < idxVarModMassMap.get(idxCombination[4]).size(); ++i4) {
                                                        Map<Integer, Float> localIdxModMassMap = new HashMap<>(6, 1);
                                                        localIdxModMassMap.put(idxCombination[0], idxVarModMassMap.get(idxCombination[0]).get(i0));
                                                        localIdxModMassMap.put(idxCombination[1], idxVarModMassMap.get(idxCombination[1]).get(i1));
                                                        localIdxModMassMap.put(idxCombination[2], idxVarModMassMap.get(idxCombination[2]).get(i2));
                                                        localIdxModMassMap.put(idxCombination[3], idxVarModMassMap.get(idxCombination[3]).get(i3));
                                                        localIdxModMassMap.put(idxCombination[4], idxVarModMassMap.get(idxCombination[4]).get(i4));
                                                        permutationNum = subFun(candidate, localIdxModMassMap, fixModIdxes, permutationNum);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                } else if (idxCombination.length == 4) {
                                    for (int i0 = 0; i0 < idxVarModMassMap.get(idxCombination[0]).size(); ++i0) {
                                        for (int i1 = 0; i1 < idxVarModMassMap.get(idxCombination[1]).size(); ++i1) {
                                            for (int i2 = 0; i2 < idxVarModMassMap.get(idxCombination[2]).size(); ++i2) {
                                                for (int i3 = 0; i3 < idxVarModMassMap.get(idxCombination[3]).size(); ++i3) {
                                                    Map<Integer, Float> localIdxModMassMap = new HashMap<>(5, 1);
                                                    localIdxModMassMap.put(idxCombination[0], idxVarModMassMap.get(idxCombination[0]).get(i0));
                                                    localIdxModMassMap.put(idxCombination[1], idxVarModMassMap.get(idxCombination[1]).get(i1));
                                                    localIdxModMassMap.put(idxCombination[2], idxVarModMassMap.get(idxCombination[2]).get(i2));
                                                    localIdxModMassMap.put(idxCombination[3], idxVarModMassMap.get(idxCombination[3]).get(i3));
                                                    permutationNum = subFun(candidate, localIdxModMassMap, fixModIdxes, permutationNum);
                                                }
                                            }
                                        }
                                    }
                                } else if (idxCombination.length == 3) {
                                    for (int i0 = 0; i0 < idxVarModMassMap.get(idxCombination[0]).size(); ++i0) {
                                        for (int i1 = 0; i1 < idxVarModMassMap.get(idxCombination[1]).size(); ++i1) {
                                            for (int i2 = 0; i2 < idxVarModMassMap.get(idxCombination[2]).size(); ++i2) {
                                                Map<Integer, Float> localIdxModMassMap = new HashMap<>(4, 1);
                                                localIdxModMassMap.put(idxCombination[0], idxVarModMassMap.get(idxCombination[0]).get(i0));
                                                localIdxModMassMap.put(idxCombination[1], idxVarModMassMap.get(idxCombination[1]).get(i1));
                                                localIdxModMassMap.put(idxCombination[2], idxVarModMassMap.get(idxCombination[2]).get(i2));
                                                permutationNum = subFun(candidate, localIdxModMassMap, fixModIdxes, permutationNum);
                                            }
                                        }
                                    }
                                } else if (idxCombination.length == 2) {
                                    for (int i0 = 0; i0 < idxVarModMassMap.get(idxCombination[0]).size(); ++i0) {
                                        for (int i1 = 0; i1 < idxVarModMassMap.get(idxCombination[1]).size(); ++i1) {
                                            Map<Integer, Float> localIdxModMassMap = new HashMap<>(3, 1);
                                            localIdxModMassMap.put(idxCombination[0], idxVarModMassMap.get(idxCombination[0]).get(i0));
                                            localIdxModMassMap.put(idxCombination[1], idxVarModMassMap.get(idxCombination[1]).get(i1));
                                            permutationNum = subFun(candidate, localIdxModMassMap, fixModIdxes, permutationNum);
                                        }
                                    }
                                } else if (idxCombination.length == 1) {
                                    for (int i0 = 0; i0 < idxVarModMassMap.get(idxCombination[0]).size(); ++i0) {
                                        Map<Integer, Float> localIdxModMassMap = new HashMap<>(2, 1);
                                        localIdxModMassMap.put(idxCombination[0], idxVarModMassMap.get(idxCombination[0]).get(i0));
                                        permutationNum = subFun(candidate, localIdxModMassMap, fixModIdxes, permutationNum);
                                    }
                                }
                            }
                            //  Don't stop in the middle. Check all permutations for a certain i
                            if (permutationNum > maxPermutationNum) {
                                break;
                            }
                        }
                    }
                }
            }

            // permutate modification masses inferred from dynamic programming
            for (Peptide candidate : candidates) {
                if (candidate.hasVarPTM() && (candidate.getVarPTMNum() > 1)) { // One PTM case has already been considered before.
                    Set<Integer> fixModIdxes = getFixModIdxes(candidate.getPTMFreeSeq(), fixModMap);
                    int permutationNum = 0;
                    boolean stop = false;
                    Float[] modMassArray = candidate.getVarPTMs().values().toArray(new Float[candidate.getVarPTMNum()]);
                    List<Float[]> permutedModMassArray = new LinkedList<>();
                    permuteModMassArray(modMassArray, 0, permutedModMassArray);
                    Iterator<int[]> tempIterator = CombinatoricsUtils.combinationsIterator(candidate.getPTMFreeSeq().length(), modMassArray.length);
                    while (tempIterator.hasNext()) {
                        int[] idxArray = tempIterator.next();
                        boolean ok = true;
                        for (int idx : idxArray) {
                            if (fixModIdxes.contains(idx)) {
                                ok = false;
                                break;
                            }
                        }
                        if (ok) {
                            Arrays.sort(idxArray);
                            for (Float[] v : permutedModMassArray) {
                                PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap(candidate.getPTMFreeSeq().length());
                                for (int i = 0; i < idxArray.length; ++i) {
                                    positionDeltaMassMap.put(new Coordinate(idxArray[i], idxArray[i] + 1), v[i]);
                                }
                                Peptide peptideObj = new Peptide(candidate.getPTMFreeSeq(), candidate.isDecoy(), massToolObj, maxMs2Charge, candidate.getNormalizedCrossCorr(), candidate.getLeftFlank(), candidate.getRightFlank(), candidate.getGlobalRank());
                                peptideObj.setVarPTM(positionDeltaMassMap);
                                peptideObj.setUnknownPtmNum(candidate.getUnknownPtmNum());
                                if (!checkedSequenceSet.contains(peptideObj)) {
                                    CalScore.calScore(peptideObj, expProcessedPL, psm, massToolObj, modSequences);
                                    checkedSequenceSet.add(peptideObj);
                                    ++permutationNum;
                                    if (permutationNum > maxPermutationNum) {
                                        stop = true;
                                        break;
                                    }
                                }
                            }
                        }
                        if (stop) {
                            break;
                        }
                    }
                }
            }

            if (psm.getPeptide() != null) {
                // calculate PTM delta score
                double ptmDeltaScore;
                TreeSet<PeptideScore> temp = modSequences.get(psm.getPeptide().getPTMFreeSeq());
                if (temp.size() == 1) {
                    ptmDeltaScore = temp.first().score;
                } else {
                    Iterator<PeptideScore> it = temp.iterator();
                    ptmDeltaScore = it.next().score - it.next().score;
                }
                psm.setPtmDeltasScore(ptmDeltaScore);
                psm.setPtmPatterns(temp);
            }
        }

        return psm;
    }

    private static void permuteModMassArray(Float[] modMassArray, int index, List<Float[]> resultList){
        if(index >= modMassArray.length - 1){
            Float[] temp = Arrays.copyOf(modMassArray, modMassArray.length);
            resultList.add(temp);
            return;
        }

        for(int i = index; i < modMassArray.length; i++){
            float t = modMassArray[index];
            modMassArray[index] = modMassArray[i];
            modMassArray[i] = t;

            permuteModMassArray(modMassArray, index + 1, resultList);

            t = modMassArray[index];
            modMassArray[index] = modMassArray[i];
            modMassArray[i] = t;
        }
    }

    static boolean isNewPtmMass(Set<VarModParam> varModParamSet, float mass, float tolerance) {
        for (VarModParam varModParam : varModParamSet) {
            if (Math.abs(varModParam.modMass - mass) < tolerance) {
                return false;
            }
        }
        return true;
    }

    private boolean isDeltaMassMeaningful(AA[] aaArray, float deltaMass, float tolerance) {
        for (AA aa : aaArray) {
            if (Math.abs(aa.ptmDeltaMass + deltaMass) < tolerance) {
                return false;
            }
        }
        return true;
    }

    private Set<Integer> getFixModIdxes(String ptmFreeSeq, Map<Character, Float> fixModMap) {
        Set<Integer> outputSet = new HashSet<>(ptmFreeSeq.length(), 1);
        char[] tempArray = ptmFreeSeq.toCharArray();
        for (int i = 0; i < tempArray.length; ++i) {
            if (Math.abs(fixModMap.get(tempArray[i])) > 0.1) {
                outputSet.add(i);
            }
        }
        return outputSet;
    }

    private int subFun(Peptide candidate, Map<Integer, Float> localIdxModMassMap, Set<Integer> fixModIdxes, int permutationNum) {
        String ptmFreeSequence = candidate.getPTMFreeSeq();
        StringBuilder sb = new StringBuilder(ptmFreeSequence.length() * 10);
        for (int j = 0; j < ptmFreeSequence.length(); ++j) {
            sb.append(ptmFreeSequence.charAt(j));
            if (localIdxModMassMap.containsKey(j)) {
                sb.append(String.format("(%.1f)", localIdxModMassMap.get(j)));
            }
        }
        String varSeq = sb.toString();

        // Calculate Score
        float seqMass = massToolObj.calResidueMass(varSeq) + MassTool.H2O;
        float deltaMass = spectrumEntry.precursorMass - seqMass;
        if (Math.abs(deltaMass) >= tolerance) {
            if (isNewPtmMass(varModParamSet, deltaMass, tolerance)) {
                // there are one more unknown modification
                AA[] aaArray = MassTool.seqToAAList(varSeq);
                if (isDeltaMassMeaningful(aaArray, deltaMass, tolerance)) {
                    boolean triedNTerm = false;
                    boolean triedLastAA = false;
                    for (int k = 0; k < aaArray.length; ++k) {
                        if ((triedNTerm && (k == 1)) || (triedLastAA && (k == aaArray.length - 1))) {
                            // N-term and the first amino acid are actually the same position. N-term has a higher priority.
                            // C-term and the last amino acid are actually the same position. The last amino acid has a higher priority.
                            continue;
                        }
                        if (!aaArray[k].hasMod() && !fixModIdxes.contains(k)) {
                            PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap(candidate.getPTMFreeSeq().length());
                            for (int l = 0; l < aaArray.length; ++l) {
                                if (aaArray[l].hasMod()) {
                                    positionDeltaMassMap.put(new Coordinate(l, l + 1), aaArray[l].ptmDeltaMass);
                                }
                            }
                            positionDeltaMassMap.put(new Coordinate(k, k + 1), deltaMass);
                            Peptide peptideObj = new Peptide(candidate.getPTMFreeSeq(), candidate.isDecoy(), massToolObj, maxMs2Charge, candidate.getNormalizedCrossCorr(), candidate.getLeftFlank(), candidate.getRightFlank(), candidate.getGlobalRank());
                            peptideObj.setVarPTM(positionDeltaMassMap);
                            peptideObj.setUnknownPtmNum(1);
                            if (!checkedSequenceSet.contains(peptideObj)) {
                                CalScore.calScore(peptideObj, expProcessedPL, psm, massToolObj, modSequences);
                                checkedSequenceSet.add(peptideObj);
                                ++permutationNum;
                                if (k == 0) {
                                    triedNTerm = true;
                                } else if (k == aaArray.length - 2) {
                                    triedLastAA = true;
                                }
                            }
                        }
                    }
                }
            }
        } else {
            PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap(candidate.getPTMFreeSeq().length());
            AA[] aaArray = MassTool.seqToAAList(varSeq);
            for (int k = 0; k < aaArray.length; ++k) {
                if (aaArray[k].hasMod()) {
                    positionDeltaMassMap.put(new Coordinate(k, k + 1), aaArray[k].ptmDeltaMass);
                }
            }
            Peptide peptideObj = new Peptide(candidate.getPTMFreeSeq(), candidate.isDecoy(), massToolObj, maxMs2Charge, candidate.getNormalizedCrossCorr(), candidate.getLeftFlank(), candidate.getRightFlank(), candidate.getGlobalRank());
            peptideObj.setVarPTM(positionDeltaMassMap);
            peptideObj.setUnknownPtmNum(0);
            if (!checkedSequenceSet.contains(peptideObj)) {
                CalScore.calScore(peptideObj, expProcessedPL, psm, massToolObj, modSequences);
                checkedSequenceSet.add(peptideObj);
                ++permutationNum;
            }
        }
        return permutationNum;
    }
}
