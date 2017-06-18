package proteomics.PTM;


import org.apache.commons.math3.util.CombinatoricsUtils;
import proteomics.Search.CalXcorr;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;

import java.util.*;

public class GeneratePtmCandidatesCalculateXcorr {

    private static final int globalVarModMaxNum = 5; // This value cannot be larger than 5. Otherwise, change generateLocalIdxModMassMap accordingly.
    private static final int maxPermutationNum = 1000;

    private final SpectrumEntry spectrumEntry;
    private final MassTool massToolObj;
    private final int maxMs2Charge;
    private final float tolerance;
    private final Set<VarModParam> varModParamSet;
    private final Map<Character, Float> fixModMap;
    private final SparseVector expXcorrPl;
    private FinalResultEntry psm;
    private Set<String> checkedSequenceSet = new HashSet<>(maxPermutationNum * 2 + 50, 1); // record checked sequence to avoid recording the same sequence twice
    private Map<String, LinkedList<PeptideScore>> modSequences = new HashMap<>(15, 1);

    public GeneratePtmCandidatesCalculateXcorr(SpectrumEntry spectrumEntry, MassTool massToolObj, Set<VarModParam> varModParamSet, Map<Character, Float> fixModMap, float ms2Tolerance, int maxMs2Charge, SparseVector expXcorrPl, int scanNum, int precursorCharge, float precursorMz) {
        this.spectrumEntry = spectrumEntry;
        this.massToolObj = massToolObj;
        this.maxMs2Charge = maxMs2Charge;
        tolerance = Math.max(ms2Tolerance, 1);
        this.varModParamSet = varModParamSet;
        this.fixModMap = fixModMap;
        this.expXcorrPl = expXcorrPl;
        psm = new FinalResultEntry(scanNum, precursorCharge, precursorMz);
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

    public FinalResultEntry generateAllPtmCandidatesCalculateXcorr(Set<Peptide> candidates) {
        // Add all PTM sequence given known variable modifications.
        if (!candidates.isEmpty()) {
            for (Peptide candidate : candidates) {
                if (candidate.hasVarPTM()) {
                    Set<Integer> fixModIdxes = getFixModIdxes(candidate.getPTMFreeSeq(), fixModMap);
                    // having only one unknown modification
                    float deltaMass = spectrumEntry.precursorMass - candidate.getPrecursorMass();
                    if (deltaMass >= tolerance) {
                        if (!isKnownPtmMass(varModParamSet, deltaMass, tolerance)) {
                            for (int i = 0; i < candidate.getPTMFreeSeq().length(); ++i) {
                                if (!fixModIdxes.contains(i)) {
                                    PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap();
                                    positionDeltaMassMap.put(new Coordinate(i, i + 1), deltaMass);
                                    Peptide peptideObj = new Peptide(candidate.getPTMFreeSeq(), candidate.isDecoy(), massToolObj, maxMs2Charge, candidate.getNormalizedCrossCorr(), candidate.getLeftFlank(), candidate.getRightFlank(), candidate.getGlobalRank());
                                    peptideObj.setVarPTM(positionDeltaMassMap);
                                    if (!checkedSequenceSet.contains(peptideObj.getVarPtmContainingSeq())) {
                                        CalXcorr.calXcorr(peptideObj, expXcorrPl, psm, massToolObj, modSequences);
                                        checkedSequenceSet.add(peptideObj.getVarPtmContainingSeq());
                                    }
                                }
                            }
                        }
                    }

                    // having known modifications and maybe an unknown modification.
                    int permutationNum = 0;
                    boolean stop = false;
                    String ptmFreeSequence = candidate.getPTMFreeSeq();
                    Map<Integer, List<Float>> idxVarModMassMap = new HashMap<>(ptmFreeSequence.length(), 1);
                    for (int i = 0; i < ptmFreeSequence.length(); ++i) {
                        char aa = ptmFreeSequence.charAt(i);
                        for (VarModParam varModParam : varModParamSet) {
                            if (varModParam.aa == aa) {
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
                    if (!idxVarModMassMap.isEmpty()) {
                        Integer[] allIdxArray = idxVarModMassMap.keySet().toArray(new Integer[idxVarModMassMap.size()]);
                        Arrays.sort(allIdxArray);
                        for (int i = 1; i <= Math.min(idxVarModMassMap.size(), globalVarModMaxNum); ++i) {
                            List<int[]> idxCombinationList = generateIdxCombinations(allIdxArray, i);
                            for (int[] idxCombination : idxCombinationList) {
                                List<Map<Integer, Float>> localIdxModMassMaps = generateLocalIdxModMassMap(idxCombination, idxVarModMassMap);
                                for (Map<Integer, Float> localIdxModMassMap : localIdxModMassMaps) {
                                    StringBuilder sb = new StringBuilder(ptmFreeSequence.length() * 10);
                                    for (int j = 0; j < ptmFreeSequence.length(); ++j) {
                                        sb.append(ptmFreeSequence.charAt(j));
                                        if (localIdxModMassMap.containsKey(j)) {
                                            sb.append(String.format("(%.1f)", localIdxModMassMap.get(j)));
                                        }
                                    }
                                    String varSeq = sb.toString();

                                    // Calculate XCorr
                                    float seqMass = massToolObj.calResidueMass(varSeq) + MassTool.H2O;
                                    deltaMass = spectrumEntry.precursorMass - seqMass;
                                    if (Math.abs(deltaMass) >= tolerance) {
                                        if (!isKnownPtmMass(varModParamSet, deltaMass, tolerance)) {
                                            // there are one more unknown modification
                                            AA[] aaArray = massToolObj.seqToAAList(varSeq);
                                            for (int k = 0; k < aaArray.length; ++k) {
                                                if (!aaArray[k].hasMod() && !fixModIdxes.contains(k)) {
                                                    PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap();
                                                    for (int l = 0; l < aaArray.length; ++l) {
                                                        if (aaArray[l].hasMod()) {
                                                            positionDeltaMassMap.put(new Coordinate(l, l + 1), aaArray[l].ptmDeltaMass);
                                                        }
                                                    }
                                                    positionDeltaMassMap.put(new Coordinate(k, k + 1), deltaMass);
                                                    Peptide peptideObj = new Peptide(candidate.getPTMFreeSeq(), candidate.isDecoy(), massToolObj, maxMs2Charge, candidate.getNormalizedCrossCorr(), candidate.getLeftFlank(), candidate.getRightFlank(), candidate.getGlobalRank());
                                                    peptideObj.setVarPTM(positionDeltaMassMap);
                                                    if (!checkedSequenceSet.contains(peptideObj.getVarPtmContainingSeq())) {
                                                        CalXcorr.calXcorr(peptideObj, expXcorrPl, psm, massToolObj, modSequences);
                                                        checkedSequenceSet.add(peptideObj.getVarPtmContainingSeq());
                                                        ++permutationNum;
                                                        if (permutationNum > maxPermutationNum) {
                                                            stop = true;
                                                            break;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    } else if (!checkedSequenceSet.contains(varSeq)) {
                                        PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap();
                                        AA[] aaArray = massToolObj.seqToAAList(varSeq);
                                        for (int k = 0; k < aaArray.length; ++k) {
                                            if (aaArray[k].hasMod()) {
                                                positionDeltaMassMap.put(new Coordinate(k, k + 1), aaArray[k].ptmDeltaMass);
                                            }
                                        }
                                        Peptide peptideObj = new Peptide(candidate.getPTMFreeSeq(), candidate.isDecoy(), massToolObj, maxMs2Charge, candidate.getNormalizedCrossCorr(), candidate.getLeftFlank(), candidate.getRightFlank(), candidate.getGlobalRank());
                                        peptideObj.setVarPTM(positionDeltaMassMap);
                                        CalXcorr.calXcorr(peptideObj, expXcorrPl, psm, massToolObj, modSequences);
                                        checkedSequenceSet.add(peptideObj.getVarPtmContainingSeq());
                                        ++permutationNum;
                                        if (permutationNum > maxPermutationNum) {
                                            stop = true;
                                        }
                                    }
                                    if (stop) {
                                        break;
                                    }
                                }
                                if (stop) {
                                    break;
                                }
                            }
                            if (stop) {
                                break;
                            }
                        }
                    }
                }
            }

            // permutate modification masses inferred from dynamic programming
            for (Peptide candidate : candidates) {
                if (candidate.hasVarPTM()) {
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
                                PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap();
                                for (int i = 0; i < idxArray.length; ++i) {
                                    positionDeltaMassMap.put(new Coordinate(idxArray[i], idxArray[i] + 1), v[i]);
                                }
                                Peptide peptideObj = new Peptide(candidate.getPTMFreeSeq(), candidate.isDecoy(), massToolObj, maxMs2Charge, candidate.getNormalizedCrossCorr(), candidate.getLeftFlank(), candidate.getRightFlank(), candidate.getGlobalRank());
                                peptideObj.setVarPTM(positionDeltaMassMap);
                                if (!checkedSequenceSet.contains(peptideObj.getVarPtmContainingSeq())) {
                                    CalXcorr.calXcorr(peptideObj, expXcorrPl, psm, massToolObj, modSequences);
                                    checkedSequenceSet.add(peptideObj.getVarPtmContainingSeq());
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
                String secondBestPtmPattern;
                LinkedList<PeptideScore> temp = modSequences.get(psm.getPeptide().getPTMFreeSeq());
                if (temp.size() == 1) {
                    ptmDeltaScore = temp.peekFirst().score;
                    secondBestPtmPattern = "-";
                } else {
                    ptmDeltaScore = temp.peekFirst().score - temp.get(1).score;
                    secondBestPtmPattern = temp.get(1).peptide.getPtmContainingSeq(fixModMap);
                }
                psm.setPtmDeltasScore(ptmDeltaScore);
                psm.setSecondBestPtmPattern(secondBestPtmPattern);
            }
        }

        return psm;
    }

    private List<int[]> generateIdxCombinations(Integer[] allIdxArray, int num) {
        List<int[]> outputList = new LinkedList<>();
        Iterator<int[]> tempIterator = CombinatoricsUtils.combinationsIterator(allIdxArray.length, num);
        while (tempIterator.hasNext()) {
            int[] idxIdxArray = tempIterator.next();
            int[] idxArray = new int[idxIdxArray.length];
            for (int i = 0; i < idxIdxArray.length; ++i) {
                idxArray[i] = allIdxArray[idxIdxArray[i]];
            }
            Arrays.sort(idxArray);
            outputList.add(idxArray);
        }

        return outputList;
    }

    private List<Map<Integer, Float>> generateLocalIdxModMassMap(int[] idxArray, Map<Integer, List<Float>> idxModMassMap) {
        List<Map<Integer, Float>> outputList = new LinkedList<>();
        if (idxArray.length == 5) {
            for (int i0 = 0; i0 < idxModMassMap.get(idxArray[0]).size(); ++i0) {
                for (int i1 = 0; i1 < idxModMassMap.get(idxArray[1]).size(); ++i1) {
                    for (int i2 = 0; i2 < idxModMassMap.get(idxArray[2]).size(); ++i2) {
                        for (int i3 = 0; i3 < idxModMassMap.get(idxArray[3]).size(); ++i3) {
                            for (int i4 = 0; i4 < idxModMassMap.get(idxArray[4]).size(); ++i4) {
                                Map<Integer, Float> localIdxModMassMap = new HashMap<>(6, 1);
                                localIdxModMassMap.put(idxArray[0], idxModMassMap.get(idxArray[0]).get(i0));
                                localIdxModMassMap.put(idxArray[1], idxModMassMap.get(idxArray[1]).get(i1));
                                localIdxModMassMap.put(idxArray[2], idxModMassMap.get(idxArray[2]).get(i2));
                                localIdxModMassMap.put(idxArray[3], idxModMassMap.get(idxArray[3]).get(i3));
                                localIdxModMassMap.put(idxArray[4], idxModMassMap.get(idxArray[4]).get(i4));
                                outputList.add(localIdxModMassMap);
                            }
                        }
                    }
                }
            }
        } else if (idxArray.length == 4) {
            for (int i0 = 0; i0 < idxModMassMap.get(idxArray[0]).size(); ++i0) {
                for (int i1 = 0; i1 < idxModMassMap.get(idxArray[1]).size(); ++i1) {
                    for (int i2 = 0; i2 < idxModMassMap.get(idxArray[2]).size(); ++i2) {
                        for (int i3 = 0; i3 < idxModMassMap.get(idxArray[3]).size(); ++i3) {
                            Map<Integer, Float> localIdxModMassMap = new HashMap<>(5, 1);
                            localIdxModMassMap.put(idxArray[0], idxModMassMap.get(idxArray[0]).get(i0));
                            localIdxModMassMap.put(idxArray[1], idxModMassMap.get(idxArray[1]).get(i1));
                            localIdxModMassMap.put(idxArray[2], idxModMassMap.get(idxArray[2]).get(i2));
                            localIdxModMassMap.put(idxArray[3], idxModMassMap.get(idxArray[3]).get(i3));
                            outputList.add(localIdxModMassMap);
                        }
                    }
                }
            }
        } else if (idxArray.length == 3) {
            for (int i0 = 0; i0 < idxModMassMap.get(idxArray[0]).size(); ++i0) {
                for (int i1 = 0; i1 < idxModMassMap.get(idxArray[1]).size(); ++i1) {
                    for (int i2 = 0; i2 < idxModMassMap.get(idxArray[2]).size(); ++i2) {
                        Map<Integer, Float> localIdxModMassMap = new HashMap<>(4, 1);
                        localIdxModMassMap.put(idxArray[0], idxModMassMap.get(idxArray[0]).get(i0));
                        localIdxModMassMap.put(idxArray[1], idxModMassMap.get(idxArray[1]).get(i1));
                        localIdxModMassMap.put(idxArray[2], idxModMassMap.get(idxArray[2]).get(i2));
                        outputList.add(localIdxModMassMap);
                    }
                }
            }
        } else if (idxArray.length == 2) {
            for (int i0 = 0; i0 < idxModMassMap.get(idxArray[0]).size(); ++i0) {
                for (int i1 = 0; i1 < idxModMassMap.get(idxArray[1]).size(); ++i1) {
                    Map<Integer, Float> localIdxModMassMap = new HashMap<>(3, 1);
                    localIdxModMassMap.put(idxArray[0], idxModMassMap.get(idxArray[0]).get(i0));
                    localIdxModMassMap.put(idxArray[1], idxModMassMap.get(idxArray[1]).get(i1));
                    outputList.add(localIdxModMassMap);
                }
            }
        } else if (idxArray.length == 1) {
            for (int i0 = 0; i0 < idxModMassMap.get(idxArray[0]).size(); ++i0) {
                Map<Integer, Float> localIdxModMassMap = new HashMap<>(2, 1);
                localIdxModMassMap.put(idxArray[0], idxModMassMap.get(idxArray[0]).get(i0));
                outputList.add(localIdxModMassMap);
            }
        }

        return outputList;
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

    private boolean isKnownPtmMass(Set<VarModParam> varModParamSet, float mass, float tolerance) {
        for (VarModParam varModParam : varModParamSet) {
            if (Math.abs(varModParam.modMass - mass) < tolerance) {
                return true;
            }
        }
        return false;
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
}
