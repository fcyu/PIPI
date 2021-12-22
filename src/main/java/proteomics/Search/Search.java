/*
 * Copyright 2016-2019 The Hong Kong University of Science and Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package proteomics.Search;

import proteomics.Index.BuildIndex;
import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Types.*;
import proteomics.Types.*;

import java.util.*;

public class Search {

    private static final int rankNum = 5;

    private List<Peptide> ptmOnlyResult = new LinkedList<>();
    private List<Peptide> ptmFreeResult = new LinkedList<>();


    public Search(BuildIndex buildIndex, double precursorMass, SparseVector scanCode, MassTool massTool, double ms1Tolerance, double leftInverseMs1Tolerance, double rightInverseMs1Tolerance, int ms1ToleranceUnit, double minPtmMass, double maxPtmMass, int localMaxMs2Charge) {
        PriorityQueue<ResultEntry> ptmFreeQueue = new PriorityQueue<>(rankNum * 2);
        PriorityQueue<ResultEntry> ptmOnlyQueue = new PriorityQueue<>(rankNum * 2);
        double scanNormSquare = scanCode.norm2square();
        double leftTol = ms1Tolerance;
        double rightTol = ms1Tolerance;
        if (ms1ToleranceUnit == 1) {
            leftTol = precursorMass - (precursorMass * leftInverseMs1Tolerance);
            rightTol = (precursorMass * rightInverseMs1Tolerance) - precursorMass;
        }
        double leftMass = Math.max(precursorMass + minPtmMass - leftTol, buildIndex.getMinPeptideMass());
        double rightMass = Math.min(precursorMass + maxPtmMass + rightTol, buildIndex.getMaxPeptideMass());

        if (leftMass >= rightMass) {
            return;
        }

        Map<String, Peptide0> peptide0Map = buildIndex.getPeptide0Map();
        TreeMap<Double, Set<String>> massPeptideMap = buildIndex.getMassPeptideMap();

        NavigableMap<Double, Set<String>> subMassPeptideMap = massPeptideMap.subMap(leftMass, true, rightMass, true);

        if (!subMassPeptideMap.isEmpty()) {
            for (double mass : subMassPeptideMap.keySet()) {
                for (String sequence : massPeptideMap.get(mass)) {
                    Peptide0 peptide0 = peptide0Map.get(sequence);
                    double score = 0;
                    double temp1 = Math.sqrt(peptide0.code.norm2square() * scanNormSquare);
                    if (temp1 > 1e-6) {
                        score = peptide0.code.dot(scanCode) / temp1;
                    }
                    double deltaMass = mass - precursorMass; // caution: the order matters under ms1ToleranceUnit == 1 situation

                    if (peptide0.isTarget) {
                        if ((deltaMass <= rightTol) && (deltaMass >= -1 * leftTol)) {
                            // PTM-free
                            if (ptmFreeQueue.size() < rankNum) {
                                ptmFreeQueue.add(new ResultEntry(score, sequence, false));
                            } else {
                                if (score > ptmFreeQueue.peek().score) {
                                    ptmFreeQueue.poll();
                                    ptmFreeQueue.add(new ResultEntry(score, sequence, false));
                                }
                            }
                        }

                        if ((deltaMass > rightTol) || (deltaMass < -1 * leftTol)) {
                            // PTM-only
                            if (ptmOnlyQueue.size() < rankNum) {
                                ptmOnlyQueue.add(new ResultEntry(score, sequence, false));
                            } else {
                                if (score > ptmOnlyQueue.peek().score) {
                                    ptmOnlyQueue.poll();
                                    ptmOnlyQueue.add(new ResultEntry(score, sequence, false));
                                }
                            }
                        }
                    } else {
                        if ((deltaMass <= rightTol) && (deltaMass >= -1 * leftTol)) {
                            // PTM-free
                            if (ptmFreeQueue.size() < rankNum) {
                                ptmFreeQueue.add(new ResultEntry(score, sequence, true));
                            } else {
                                if (score > ptmFreeQueue.peek().score) {
                                    ptmFreeQueue.poll();
                                    ptmFreeQueue.add(new ResultEntry(score, sequence, true));
                                }
                            }
                        }

                        if ((deltaMass > rightTol) || (deltaMass < -1 * leftTol)) {
                            // PTM-only
                            if (ptmOnlyQueue.size() < rankNum) {
                                ptmOnlyQueue.add(new ResultEntry(score, sequence, true));
                            } else {
                                if (score > ptmOnlyQueue.peek().score) {
                                    ptmOnlyQueue.poll();
                                    ptmOnlyQueue.add(new ResultEntry(score, sequence, true));
                                }
                            }
                        }
                    }
                }
            }
        }

        if (!(ptmFreeQueue.isEmpty() && ptmOnlyQueue.isEmpty())) {
            ptmFreeResult = convertResult(ptmFreeQueue, massTool, localMaxMs2Charge);
            ptmOnlyResult = convertResult(ptmOnlyQueue, massTool, localMaxMs2Charge);
        }
    }

    private List<Peptide> convertResult(PriorityQueue<ResultEntry> inputQueue, MassTool massTool, int localMaxMs2Charge) {
        List<Peptide> peptideList = new LinkedList<>();
        int globalRank = inputQueue.size();
        while (!inputQueue.isEmpty()) {
            ResultEntry temp = inputQueue.poll();
            peptideList.add(new Peptide(temp.peptide, temp.isDecoy(), massTool, localMaxMs2Charge, temp.score, globalRank));
            --globalRank;
        }

        return peptideList;
    }

    public List<Peptide> getPTMOnlyResult() {
        return ptmOnlyResult;
    }

    public List<Peptide> getPTMFreeResult() {
        return ptmFreeResult;
    }
}
