package proteomics.DynamicProgramming;


import proteomics.Types.ExpAA;
import proteomics.Types.Coordinate;
import proteomics.Types.Peptide;
import proteomics.Types.ThreeExpAA;

import java.util.*;

public class TagAlignment {

    private final List<ExpAA> alignedList;
    private static final float ptmPenalty = -2;
    private static final float negativePtmPenalty = 2 * ptmPenalty;

    public TagAlignment(ThreeExpAA[] sortedAasArray, Peptide peptide, float NTermMz, float CTermMz, float ms2Tolerance) {
        float[] ptmFreeBIonMz = peptide.getChargeOneBIonArray();
        ms2Tolerance = Math.max(ms2Tolerance, 0.02f); // relax the precision in high resolution situation

        Cell[][] matrix = new Cell[sortedAasArray.length + 2][peptide.length()];

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
            if ((1 < expAasStartIdx) && (Math.abs(NTermMz - expAasStartMz) > ms2Tolerance)) {
                // there is a gap. but we don't know if there are PTMs in this gap.
                float theoMzGap = ptmFreeBIonMz[expAasStartIdx - 2];
                float scoreTemp;
                if (Math.abs(expAasStartMz - theoMzGap) > ms2Tolerance) {
                    // there are PTMs in the gap
                    scoreTemp = matrix[0][0].v + expAasStart.getTotalHalfIntensity() + (expAasStartMz - theoMzGap > 0 ? ptmPenalty : negativePtmPenalty);
                    if (scoreTemp > score) {
                        score = scoreTemp;
                        from = new Coordinate(0, 0);
                    }
                } else {
                    // there is no PTM in the gap
                    scoreTemp = matrix[0][0].v + expAasStart.getTotalHalfIntensity();
                    if (scoreTemp > score) {
                        score = scoreTemp;
                        from = new Coordinate(0, 0);
                    }
                }
            } else if ((1 == expAasStartIdx) && (Math.abs(NTermMz - expAasStartMz) <= ms2Tolerance)) {
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
                    float theoMzGap = ptmFreeBIonMz[expAasStartIdx - 2] - ptmFreeBIonMz[expAasFormerEndIdx - 1];
                    float scoreTemp;
                    if (Math.abs(expMzGap - theoMzGap) > 2 * ms2Tolerance) {
                        // there are PTMs in the gap
                        scoreTemp = matrix[iFormer + 1][expAasFormerEndIdx].v + expAasStart.getTotalHalfIntensity() + (expMzGap - theoMzGap > 0 ? ptmPenalty : negativePtmPenalty);
                        if (scoreTemp > score) {
                            score = scoreTemp;
                            from = new Coordinate(iFormer + 1, expAasFormerEndIdx);
                        }
                    } else {
                        // there is no PTM in the gap
                        scoreTemp = matrix[iFormer + 1][expAasFormerEndIdx].v + expAasStart.getTotalHalfIntensity();
                        if (scoreTemp > score) {
                            score = scoreTemp;
                            from = new Coordinate(iFormer + 1, expAasFormerEndIdx);
                        }
                    }
                } else if ((expAasFormerEndIdx == expAasStartIdx - 1) && (expAasFormerEndMz == expAasStartMz)) {
                    // it's a AA jump
                    float scoreTemp = matrix[iFormer + 1][expAasFormerEndIdx].v + expAasStart.getTotalHalfIntensity();
                    if (scoreTemp > score) {
                        score = scoreTemp;
                        from = new Coordinate(iFormer + 1, expAasFormerEndIdx);
                    }
                } else if ((expAasFormerEndIdx > expAasStartIdx - 1) && (expAasFormerEndIdx <= expAasStartIdx + expAas.size() - 1) && (expAasFormerEndMz > expAasStartMz)) {
                    // there may be overlaps
                    int overlapIdx = findOverlap(expAasFormer, expAas);
                    if (overlapIdx != -1) {
                        float scoreTemp = matrix[iFormer + 1][expAasFormerEndIdx - overlapIdx].v + expAasStart.getTotalHalfIntensity();
                        if (scoreTemp > score) {
                            score = scoreTemp;
                            from = new Coordinate(iFormer + 1, expAasFormerEndIdx - overlapIdx);
                        }
                }
                }
            }
            matrix[i + 1][expAasStartIdx] = new Cell(score, from);

            // in the tag.
            for (int j = expAasStartIdx + 1; j < expAasStartIdx + expAas.size(); ++j) {
                matrix[i + 1][j] = new Cell(matrix[i + 1][j - 1].v + expAas.get(j - expAasStartIdx).getTotalHalfIntensity(), new Coordinate(i + 1, j - 1));
            }

            // after the tag.
            int j = expAasStartIdx + expAas.size();
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
            if ((expAasEndIdx < peptide.length() - 2) && (expAasEndMz <= CTermMz - ms2Tolerance)) {
                // there is a gap. but we don't know if there are PTMs in this gap.
                float expMzGap = CTermMz - expAasEndMz;
                float theoMzGap = ptmFreeBIonMz[ptmFreeBIonMz.length - 1] - ptmFreeBIonMz[expAasEndIdx - 1];
                float scoreTemp;
                if (Math.abs(expMzGap - theoMzGap) > ms2Tolerance) {
                    // there are PTMs in the gap
                    scoreTemp = matrix[i + 1][expAasEndIdx].v + (expMzGap - theoMzGap > 0 ? ptmPenalty : negativePtmPenalty) + 1;
                    if (scoreTemp > score) {
                        score = scoreTemp;
                        from = new Coordinate(i + 1, expAasEndIdx);
                    }
                } else {
                    // there is no PTM in the gap
                    scoreTemp = matrix[i + 1][expAasEndIdx].v + 1;
                    if (scoreTemp > score) {
                        score = scoreTemp;
                        from = new Coordinate(i + 1, expAasEndIdx);
                    }
                }
            } else if ((expAasEndIdx == peptide.length() - 2) && (Math.abs(expAasEndMz - CTermMz) <= ms2Tolerance)) {
                // there is no gap in between
                float scoreTemp = matrix[i + 1][expAasEndIdx].v + 1;
                if (scoreTemp > score) {
                    score = scoreTemp;
                    from = new Coordinate(i + 1, expAasEndIdx);
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
            if ((j < temp1.get(0).getTheoLocation() + temp1.size()) && (j >= temp1.get(0).getTheoLocation())) {
                temp.add(temp1.get(j - temp1.get(0).getTheoLocation()));
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
}
