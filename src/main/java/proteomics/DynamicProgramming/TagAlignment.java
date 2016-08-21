package proteomics.DynamicProgramming;


import proteomics.Types.ExpAA;
import proteomics.Types.Coordinate;
import proteomics.Types.Peptide;
import proteomics.Types.ThreeExpAA;

import java.util.*;

public class TagAlignment {

    private final List<ExpAA> alignedList;
    private static final float ptmPenalty = -2;

    public TagAlignment(ThreeExpAA[] sortedAasArray, Peptide peptide, Map<String, TreeSet<Integer>> siteMass100Map, float NTermMz, float CTermMz, float ms2Tolerance) {
        float[] ptmFreeBIonMz = peptide.getChargeOneBIonArray();
        ms2Tolerance = Math.max(ms2Tolerance, 0.02f); // relax the precision in high resolution situation

        Cell[][] matrix = new Cell[sortedAasArray.length + 2][peptide.length() + 2]; // consider N-term and C-term

        // N-term
        matrix[0][0] = new Cell(1, null);

        for (int i = 0; i < sortedAasArray.length; ++i) {
            ThreeExpAA expAas = sortedAasArray[i];
            ExpAA expAasStart = expAas.get(0);
            int expAasStartIdx = expAasStart.getTheoLocation();
            float expAasStartMz = expAasStart.getHeadLocation();

            // at the head of the tag
            float score = -9999;
            Coordinate from = null;
            // link to N-term
            if ((0 < expAasStartIdx) && (Math.abs(NTermMz - expAasStartMz) > ms2Tolerance)) {
                // there is a gap. but we don't know if there are PTMs in this gap.
                float theoMzGap = ptmFreeBIonMz[expAasStartIdx - 1];
                float scoreTemp;
                if (Math.abs(expAasStartMz - theoMzGap) > ms2Tolerance) {
                    // there are PTMs in the gap
                    String gapSeq = peptide.getPTMFreeSeq().substring(0, expAasStartIdx);
                    float deltaMass = expAasStartMz - theoMzGap;
                    boolean proteinN = false;
                    if (peptide.getLeftFlank().contentEquals("-")) {
                        proteinN = true;
                    }
                    // check PTM type
                    if (checkJumpConstrain(siteMass100Map, gapSeq, true, false, proteinN, false, deltaMass, 2 * ms2Tolerance)) {
                        scoreTemp = matrix[0][0].v + expAasStart.getTotalHalfIntensity() + ptmPenalty;
                        if (scoreTemp > score) {
                            score = scoreTemp;
                            from = new Coordinate(0, 0);
                        }
                    } else {
                        scoreTemp = matrix[0][0].v + expAasStart.getTotalHalfIntensity();
                        if (scoreTemp > score) {
                            score = scoreTemp;
                            from = null;
                        }
                    }
                } else {
                    // there is no PTM in the gap
                    scoreTemp = matrix[0][0].v + expAasStart.getTotalHalfIntensity();
                    if (scoreTemp > score) {
                        score = scoreTemp;
                        from = new Coordinate(0, 0);
                    }
                }
            } else if ((0 == expAasStartIdx) && (Math.abs(NTermMz - expAasStartMz) <= ms2Tolerance)) {
                // there is no pap between.
                float scoreTemp = matrix[0][0].v + expAasStart.getTotalHalfIntensity();
                if (scoreTemp > score) {
                    score = scoreTemp;
                    from = new Coordinate(0, 0);
                }
            }
            // link to another tag
            for (int iFormer = 0; iFormer < i; ++iFormer) {
                ThreeExpAA expAasFormer = sortedAasArray[iFormer];
                ExpAA expAasFormerEnd = expAasFormer.get(expAasFormer.size() - 1);
                int expAasFormerEndIdx = expAasFormerEnd.getTheoLocation();
                float expAasFormerEndMz = expAasFormerEnd.getTailLocation();
                if ((expAasFormerEndIdx < expAasStartIdx - 1) && (expAasFormerEndMz < expAasStartMz)) {
                    // there is a gap. but we don't know if there are PTMs in this gap.
                    float expMzGap = expAasStartMz - expAasFormerEndMz;
                    float theoMzGap = ptmFreeBIonMz[expAasStartIdx - 1] - ptmFreeBIonMz[expAasFormerEndIdx];
                    float scoreTemp;
                    if (Math.abs(expMzGap - theoMzGap) > 2 * ms2Tolerance) {
                        // there are PTMs in the gap
                        String gapSeq = peptide.getPTMFreeSeq().substring(expAasFormerEndIdx + 1, expAasStartIdx);
                        float deltaMass = expMzGap - theoMzGap;
                        // check PTM type
                        if (checkJumpConstrain(siteMass100Map, gapSeq, false, false, false, false, deltaMass, 2 * ms2Tolerance)) {
                            scoreTemp = matrix[iFormer + 1][expAasFormerEndIdx + 1].v + expAasStart.getTotalHalfIntensity() + ptmPenalty;
                            if (scoreTemp > score) {
                                score = scoreTemp;
                                from = new Coordinate(iFormer + 1, expAasFormerEndIdx + 1);
                            }
                        } else {
                            scoreTemp = expAasStart.getTotalHalfIntensity();
                            if (scoreTemp > score) {
                                score = scoreTemp;
                                from = null;
                            }
                        }
                    } else {
                        // there is no PTM in the gap
                        scoreTemp = matrix[iFormer + 1][expAasFormerEndIdx + 1].v + expAasStart.getTotalHalfIntensity();
                        if (scoreTemp > score) {
                            score = scoreTemp;
                            from = new Coordinate(iFormer + 1, expAasFormerEndIdx + 1);
                        }
                    }
                } else if ((expAasFormerEndIdx == expAasStartIdx - 1) && (expAasFormerEndMz == expAasStartMz)) {
                    // it's a AA jump
                    float scoreTemp = matrix[iFormer + 1][expAasFormerEndIdx + 1].v + expAasStart.getTotalHalfIntensity();
                    if (scoreTemp > score) {
                        score = scoreTemp;
                        from = new Coordinate(iFormer + 1, expAasFormerEndIdx + 1);
                    }
                } else if ((expAasFormerEndIdx > expAasStartIdx - 1) && (expAasFormerEndIdx <= expAasStartIdx + expAas.size() - 1) && (expAasFormerEndMz > expAasStartMz)) {
                    // there may be overlaps
                    int overlapIdx = findOverlap(expAasFormer, expAas);
                    if (overlapIdx != -1) {
                        float scoreTemp = matrix[iFormer + 1][expAasFormerEndIdx - overlapIdx + 1].v + expAasStart.getTotalHalfIntensity();
                        if (scoreTemp > score) {
                            score = scoreTemp;
                            from = new Coordinate(iFormer + 1, expAasFormerEndIdx - overlapIdx + 1);
                        }
                    }
                }
            }
            matrix[i + 1][expAasStartIdx + 1] = new Cell(score, from);

            // in the tag.
            for (int j = expAasStartIdx + 2; j < expAasStartIdx + expAas.size() + 1; ++j) {
                matrix[i + 1][j] = new Cell(matrix[i + 1][j - 1].v + expAas.get(j - 1 - expAasStartIdx).getTotalHalfIntensity(), new Coordinate(i + 1, j - 1));
            }

            // after the tag.
            int j = expAasStartIdx + expAas.size() + 1;
            matrix[i + 1][j] = new Cell(matrix[i + 1][j - 1].v, new Coordinate(i + 1, j - 1));
        }

        // C-term
        float score = -9999;
        Coordinate from = null;
        for (int i = 0; i < sortedAasArray.length; ++i) {
            ThreeExpAA expAas = sortedAasArray[i];
            ExpAA expAasEnd = expAas.get(expAas.size() - 1);
            float expAasEndMz = expAasEnd.getTailLocation();
            int expAasEndIdx = expAasEnd.getTheoLocation();
            if ((expAasEndIdx < peptide.length() - 1) && (expAasEndMz <= CTermMz - ms2Tolerance)) {
                // there is a gap. but we don't know if there are PTMs in this gap.
                float expMzGap = CTermMz - expAasEndMz;
                float theoMzGap = ptmFreeBIonMz[ptmFreeBIonMz.length - 1] - ptmFreeBIonMz[expAasEndIdx];
                float scoreTemp;
                if (Math.abs(expMzGap - theoMzGap) > ms2Tolerance) {
                    // there are PTMs in the gap
                    String gapSeq = peptide.getPTMFreeSeq().substring(expAasEndIdx + 1, peptide.length());
                    float deltaMass = expMzGap - theoMzGap;
                    boolean proteinC = false;
                    if (peptide.getRightFlank().contentEquals("-")) {
                        proteinC = true;
                    }
                    // check PTM type
                    if (checkJumpConstrain(siteMass100Map, gapSeq, false, true, false, proteinC, deltaMass, 2 * ms2Tolerance)) {
                        scoreTemp = matrix[i + 1][expAasEndIdx + 1].v + expAasEnd.getTotalHalfIntensity() + ptmPenalty;
                        if (scoreTemp > score) {
                            score = scoreTemp;
                            from = new Coordinate(i + 1, expAasEndIdx + 1);
                        }
                    } else {
                        scoreTemp = expAasEnd.getTotalHalfIntensity();
                        if (scoreTemp > score) {
                            score = scoreTemp;
                            from = null;
                        }
                    }
                } else {
                    // there is no PTM in the gap
                    scoreTemp = matrix[i + 1][expAasEndIdx + 1].v + expAasEnd.getTotalHalfIntensity();
                    if (scoreTemp > score) {
                        score = scoreTemp;
                        from = new Coordinate(i + 1, expAasEndIdx + 1);
                    }
                }
            } else if ((expAasEndIdx == peptide.length() - 1) && (Math.abs(expAasEndMz - CTermMz) <= ms2Tolerance)) {
                // there is no gap in between
                float scoreTemp = matrix[i + 1][expAasEndIdx + 1].v + 1;
                if (scoreTemp > score) {
                    score = scoreTemp;
                    from = new Coordinate(i + 1, expAasEndIdx + 1);
                }
            }
        }
        matrix[matrix.length - 1][matrix[0].length - 1] = new Cell(score, from);

        alignedList = backTrackTagAlignment(matrix, sortedAasArray, matrix.length - 1, matrix[0].length - 1);
    }

    public List<ExpAA> getAlignedList() {
        return alignedList;
    }

    private List<ExpAA> backTrackTagAlignment(Cell[][] matrix, ThreeExpAA[] a, int i, int j) {
        if ((i == 0) || (j == 0)) {
            return new ArrayList<>();
        } else if ((i > 0) && (j > 0) && (matrix[i][j].direction == null)) {
            return null;
        } else if (i > a.length) {
            return backTrackTagAlignment(matrix, a, matrix[i][j].direction.x, matrix[i][j].direction.y);
        } else {
            List<ExpAA> temp = backTrackTagAlignment(matrix, a, matrix[i][j].direction.x, matrix[i][j].direction.y);

            if (temp == null) {
                return null;
            }

            ThreeExpAA temp1 = a[i - 1];
            if ((j - 1 - temp1.get(0).getTheoLocation() < temp1.size()) && (j - 1 >= temp1.get(0).getTheoLocation())) {
                temp.add(temp1.get(j - 1 - temp1.get(0).getTheoLocation()));
            }
            return temp;
        }
    }

    private int findOverlap(ThreeExpAA list1, ThreeExpAA list2) {
        int idxGap = list1.get(list1.size() - 1).getTheoLocation() - list2.get(0).getTheoLocation() + 1;
        for (int i = 0; i < idxGap; ++i) {
            if (!list1.get(list1.size() - idxGap + i).equals(list2.get(i))) {
                return -1;
            }
        }
        return idxGap;
    }

    private boolean checkJumpConstrain(Map<String, TreeSet<Integer>> siteMass1000Map, String seq, boolean peptideN, boolean peptideC, boolean proteinN, boolean proteinC, float deltaMass, float tolerance) { // todo: check
        Map<Integer, int[]> localSiteMass1000Map = new HashMap<>();
        for (int i = 0; i < seq.length(); ++i) {
            String aa = seq.substring(i, i + 1);
            TreeSet<Integer> tempSet = new TreeSet<>();
            if (siteMass1000Map.containsKey(aa)) {
                tempSet.addAll(siteMass1000Map.get(aa));
            } else if ((i == 0) && peptideN && siteMass1000Map.containsKey(aa + "-PEPTIDE_N")) {
                tempSet.addAll(siteMass1000Map.get(aa + "-PEPTIDE_N"));
            } else if ((i == 0) && proteinN && siteMass1000Map.containsKey(aa + "-PROTEIN_N")) {
                tempSet.addAll(siteMass1000Map.get(aa + "-PROTEIN_N"));
            } else if ((i == seq.length() - 1) && peptideC && siteMass1000Map.containsKey(aa + "-PEPTIDE_C")) {
                tempSet.addAll(siteMass1000Map.get(aa + "-PEPTIDE_C"));
            } else if ((i == seq.length() - 1) && proteinC && siteMass1000Map.containsKey(aa + "-PROTEIN_C")) {
                tempSet.addAll(siteMass1000Map.get(aa + "-PROTEIN_C"));
            } else if ((i == 0) && peptideN && siteMass1000Map.containsKey("PEPTIDE_N")) {
                tempSet.addAll(siteMass1000Map.get("PEPTIDE_N"));
            } else if ((i == 0) && proteinN && siteMass1000Map.containsKey("PROTEIN_N")) {
                tempSet.addAll(siteMass1000Map.get("PROTEIN_N"));
            } else if ((i == seq.length() - 1) && peptideC && siteMass1000Map.containsKey("PEPTIDE_C")) {
                tempSet.addAll(siteMass1000Map.get("PEPTIDE_C"));
            } else if ((i == seq.length() - 1) && proteinC && siteMass1000Map.containsKey("PEPTIDE_C")) {
                tempSet.addAll(siteMass1000Map.get("PEPTIDE_C"));
            }

            if (!tempSet.isEmpty()) {
                int[] tempArray = new int[tempSet.size()];
                int tempIdx = 0;
                for (int mass1000 : tempSet) {
                    tempArray[tempIdx] = mass1000;
                    ++tempIdx;
                }
                localSiteMass1000Map.put(i, tempArray);
            }
        }

        int leftBoard = (int) ((deltaMass - tolerance) * 1000);
        int rightBoard = (int) ((deltaMass + tolerance) * 1000);

        int[] idxArray = new int[localSiteMass1000Map.size()];
        int tempIdx = 0;
        for (int idx : localSiteMass1000Map.keySet()) {
            idxArray[tempIdx] = idx;
            ++tempIdx;
        }
        Arrays.sort(idxArray);

        boolean ok = false;

        // one PTM
        for (int idx : idxArray) {
            int[] mass1000Array = localSiteMass1000Map.get(idx);
            for (int mass1000 : mass1000Array) {
                if ((mass1000 >= leftBoard) && (mass1000 <= rightBoard)) {
                    ok = true;
                    break;
                }
            }
            if (ok) {
                break;
            }
        }

        // two PTM
        if ((!ok) && (idxArray.length > 1)) {
            for (int i = 0; i < idxArray.length - 1; ++i) {
                int[] mass1000Array1 = localSiteMass1000Map.get(idxArray[i]);
                for (int j = i + 1; j < idxArray.length; ++j) {
                    int[] mass1000Array2 = localSiteMass1000Map.get(idxArray[j]);
                    for (int mass1000_1 : mass1000Array1) {
                        for (int mass1000_2 : mass1000Array2) {
                            int mass1000T = mass1000_1 + mass1000_2;
                            if ((mass1000T >= leftBoard) && (mass1000T <= rightBoard)) {
                                ok = true;
                                break;
                            }
                        }
                        if (ok) {
                            break;
                        }
                    }
                    if (ok) {
                        break;
                    }
                }
                if (ok) {
                    break;
                }
            }
        }

        return ok;
    }
}
