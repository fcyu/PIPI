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


    public Search(BuildIndex buildIndexObj, double precursorMass, SparseVector scanCode, MassTool massToolObj, double ms1Tolerance, double leftInverseMs1Tolerance, double rightInverseMs1Tolerance, int ms1ToleranceUnit, double minPtmMass, double maxPtmMass, int localMaxMs2Charge) {
        PriorityQueue<ResultEntry> ptmFreeQueue = new PriorityQueue<>(rankNum * 2);
        PriorityQueue<ResultEntry> ptmOnlyQueue = new PriorityQueue<>(rankNum * 2);
        double scanNormSquare = scanCode.norm2square();
        double leftTol = ms1Tolerance;
        double rightTol = ms1Tolerance;
        if (ms1ToleranceUnit == 1) {
            leftTol = precursorMass - (precursorMass * leftInverseMs1Tolerance);
            rightTol = (precursorMass * rightInverseMs1Tolerance) - precursorMass;
        }
        double leftMass = Math.max(precursorMass + minPtmMass, buildIndexObj.getMinPeptideMass());
        double rightMass = Math.min(precursorMass + maxPtmMass, buildIndexObj.getMaxPeptideMass());

        if (leftMass >= rightMass) {
            return;
        }

        Map<String, Peptide0> peptide0Map = buildIndexObj.getPeptide0Map();
        TreeMap<Double, Set<String>> massPeptideMap = buildIndexObj.getMassPeptideMap();

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
            ptmFreeResult = convertResult(ptmFreeQueue, massToolObj, localMaxMs2Charge);
            ptmOnlyResult = convertResult(ptmOnlyQueue, massToolObj, localMaxMs2Charge);
        }
    }

    private List<Peptide> convertResult(PriorityQueue<ResultEntry> inputQueue, MassTool massToolObj, int localMaxMs2Charge) {
        List<Peptide> peptideList = new LinkedList<>();
        int globalRank = inputQueue.size();
        while (!inputQueue.isEmpty()) {
            ResultEntry temp = inputQueue.poll();
            peptideList.add(new Peptide(temp.peptide, temp.isDecoy(), massToolObj, localMaxMs2Charge, temp.score, globalRank));
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
