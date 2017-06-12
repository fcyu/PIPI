package proteomics.PTM;


import org.apache.commons.math3.util.CombinatoricsUtils;
import proteomics.Index.BuildIndex;
import proteomics.Search.CalXcorr;
import proteomics.Segment.InferenceSegment;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;

import java.util.*;

public class GeneratePtmCandidatesCalculateXcorr {

    private static final int globalVarModMaxNum = 3; // This value cannot be larger than 5. Otherwise, change generateLocalIdxModMassMap accordingly.

    private final InferenceSegment inference3SegmentObj;
    private final SpectrumEntry spectrumEntry;
    private final MassTool massToolObj;
    private final BuildIndex buildIndexObj;
    private final float ms2Tolerance;
    private final int maxMs2Charge;

    public GeneratePtmCandidatesCalculateXcorr(InferenceSegment inference3SegmentObj, SpectrumEntry spectrumEntry, MassTool massToolObj, BuildIndex buildIndexObj, float ms2Tolerance, int maxMs2Charge) {
        this.inference3SegmentObj = inference3SegmentObj;
        this.spectrumEntry = spectrumEntry;
        this.massToolObj = massToolObj;
        this.buildIndexObj = buildIndexObj;
        this.ms2Tolerance = ms2Tolerance;
        this.maxMs2Charge = maxMs2Charge;
    }

    public Set<Peptide> eliminateMissedCleavageCausedPtms(List<Peptide> ptmFreeCandidates, List<Peptide> ptmOnlyCandidates) {
        // Additional short missed cleavaged amino acid sequence may be canceled out by negative PTM. In this situation, we eliminate the missed cleavaged one.
        List<Peptide> tempList = new LinkedList<>();
        tempList.addAll(ptmFreeCandidates);
        tempList.addAll(ptmOnlyCandidates);
        Peptide[] tempArray = tempList.toArray(new Peptide[tempList.size()]);
        Set<Peptide> candidates = new HashSet<>();
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

    public FinalResultEntry generateAllPtmCandidatesCalculateXcorr(Set<Peptide> ptmOnlyCandidates, SparseVector expXcorrPl, int scanNum, int precursorCharge, float precursorMz) {
        FinalResultEntry psm = new FinalResultEntry(scanNum, precursorCharge, precursorMz);
        Map<String, LinkedList<PeptideScore>> modSequences = new HashMap<>();
        Set<String> checkedSequenceSet = new HashSet<>(); // record checked sequence to avoid recording the same sequence twice

        // Add all PTM sequence given known variable modifications.
        Set<VarModParam> varModParamSet = inference3SegmentObj.getVarModParamSet();
        Map<Character, Float> fixModMap = buildIndexObj.returnFixModMap();
        float tolerance = (float) Math.max(ms2Tolerance, 0.1);
        if (!ptmOnlyCandidates.isEmpty()) {
            for (Peptide candidate : ptmOnlyCandidates) {
                // having only on unknown modification
                float deltaMass = spectrumEntry.precursorMass - candidate.getPrecursorMass();
                if (deltaMass > tolerance) {
                    if (!isKnownPtmMass(varModParamSet, deltaMass, tolerance)) {
                        for (int i = 0; i < candidate.getPTMFreeSeq().length(); ++i) {
                            if (Math.abs(fixModMap.get(candidate.getPTMFreeSeq().charAt(i))) <= 0.1) {
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
                Set<String> varSeqSet = generateModSeq(candidate.getPTMFreeSeq(), varModParamSet);
                for (String varSeq : varSeqSet) {
                    float seqMass = massToolObj.calResidueMass(varSeq);
                    deltaMass = spectrumEntry.precursorMass - seqMass;
                    if (Math.abs(deltaMass) > tolerance) {
                        if (!isKnownPtmMass(varModParamSet, deltaMass, tolerance)) {
                            // there are one more unknown modification
                            AA[] aaArray = massToolObj.seqToAAList(varSeq);
                            for (int i = 0; i < aaArray.length; ++i) {
                                if (!aaArray[i].hasMod() && (Math.abs(fixModMap.get(aaArray[i].aa)) <= 0.1)) {
                                    PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap();
                                    for (int j = 0; j < aaArray.length; ++j) {
                                        if (aaArray[j].hasMod()) {
                                            positionDeltaMassMap.put(new Coordinate(j, j + 1), aaArray[j].ptmDeltaMass);
                                        }
                                    }
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
                    } else {
                        PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap();
                        AA[] aaArray = massToolObj.seqToAAList(varSeq);
                        for (int i = 0; i < aaArray.length; ++i) {
                            if (aaArray[i].hasMod()) {
                                positionDeltaMassMap.put(new Coordinate(i, i + 1), aaArray[i].ptmDeltaMass);
                            }
                        }
                        Peptide peptideObj = new Peptide(candidate.getPTMFreeSeq(), candidate.isDecoy(), massToolObj, maxMs2Charge, candidate.getNormalizedCrossCorr(), candidate.getLeftFlank(), candidate.getRightFlank(), candidate.getGlobalRank());
                        peptideObj.setVarPTM(positionDeltaMassMap);
                        if (!checkedSequenceSet.contains(peptideObj.getVarPtmContainingSeq())) {
                            CalXcorr.calXcorr(peptideObj, expXcorrPl, psm, massToolObj, modSequences);
                            checkedSequenceSet.add(peptideObj.getVarPtmContainingSeq());
                        }
                    }
                }
            }

            // permutate modification masses inferred from dynamic programming
            for (Peptide candidate : ptmOnlyCandidates) {
                Float[] modMassArray = candidate.getVarPTMs().values().toArray(new Float[candidate.getVarPTMNum()]);
                List<Float[]> permutedModMassArray = new LinkedList<>();
                permuteModMassArray(modMassArray, 0, permutedModMassArray);
                Iterator<int[]> tempIterator = CombinatoricsUtils.combinationsIterator(candidate.getPTMFreeSeq().length(), modMassArray.length);
                while (tempIterator.hasNext()) {
                    int[] idxArray = tempIterator.next();
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

    private Set<String> generateModSeq(String seq, Set<VarModParam> varModParamSet) {
        Set<String> varSeqSet = new HashSet<>();

        // get all locations' var lists
        Map<Integer, List<Float>> idxVarModMassMap = new HashMap<>();
        for (int i = 0; i < seq.length(); ++i) {
            char aa = seq.charAt(i);
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
            // generate var containing sequences
            Integer[] allIdxArray = idxVarModMassMap.keySet().toArray(new Integer[idxVarModMassMap.size()]);
            Arrays.sort(allIdxArray);
            for (int i = 1; i <= Math.min(idxVarModMassMap.size(), globalVarModMaxNum); ++i) {
                List<int[]> idxCombinationList = generateIdxCombinations(allIdxArray, i);
                for (int[] idxCombination : idxCombinationList) {
                    List<Map<Integer, Float>> localIdxModMassMaps = generateLocalIdxModMassMap(idxCombination, idxVarModMassMap);
                    for (Map<Integer, Float> localIdxModMassMap : localIdxModMassMaps) {
                        StringBuilder sb = new StringBuilder(seq.length() * 10);
                        for (int j = 0; j < seq.length(); ++j) {
                            sb.append(seq.charAt(j));
                            if (localIdxModMassMap.containsKey(j)) {
                                sb.append(String.format("(%.1f)", localIdxModMassMap.get(j)));
                            }
                        }
                        varSeqSet.add(sb.toString());
                    }
                }
            }
        }

        return varSeqSet;
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
                                Map<Integer, Float> localIdxModMassMap = new HashMap<>();
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
                            Map<Integer, Float> localIdxModMassMap = new HashMap<>();
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
                        Map<Integer, Float> localIdxModMassMap = new HashMap<>();
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
                    Map<Integer, Float> localIdxModMassMap = new HashMap<>();
                    localIdxModMassMap.put(idxArray[0], idxModMassMap.get(idxArray[0]).get(i0));
                    localIdxModMassMap.put(idxArray[1], idxModMassMap.get(idxArray[1]).get(i1));
                    outputList.add(localIdxModMassMap);
                }
            }
        } else if (idxArray.length == 1) {
            for (int i0 = 0; i0 < idxModMassMap.get(idxArray[0]).size(); ++i0) {
                Map<Integer, Float> localIdxModMassMap = new HashMap<>();
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
            if (Math.abs(varModParam.modMass - mass) <= tolerance) {
                return true;
            }
        }
        return false;
    }
}
