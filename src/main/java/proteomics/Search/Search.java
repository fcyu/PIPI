package proteomics.Search;


import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Index.BuildIndex;
import proteomics.PIPI;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.*;
import java.util.*;

public class Search {

    private static final Logger logger = LoggerFactory.getLogger(Search.class);
    private static final int rankNum = 10;

    private List<Peptide> ptmOnlyResult = new LinkedList<>();
    private List<Peptide> ptmFreeResult = new LinkedList<>();


    public Search(BuildIndex buildIndexObj, SpectrumEntry spectrumEntry, SparseVector scanCode, String sqlPath, MassTool massToolObj, float ms1Tolerance, int ms1ToleranceUnit, float minPtmMass, float maxPtmMass, int maxMs2Charge) {
        PriorityQueue<ResultEntry> ptmFreeQueue = new PriorityQueue<>(rankNum * 2);
        PriorityQueue<ResultEntry> ptmOnlyQueue = new PriorityQueue<>(rankNum * 2);
        try {
            double scanNormSquare = scanCode.norm2square();
            float leftTol = ms1Tolerance;
            float rightTol = ms1Tolerance;
            if (ms1ToleranceUnit == 1) {
                leftTol = spectrumEntry.precursorMass - (spectrumEntry.precursorMass / (1 + ms1Tolerance * 1e-6f));
                rightTol = (spectrumEntry.precursorMass / (1 - ms1Tolerance * 1e-6f)) - spectrumEntry.precursorMass;
            }
            float leftMass = Math.max(spectrumEntry.precursorMass + minPtmMass, buildIndexObj.getMinPeptideMass());
            float rightMass = Math.min(spectrumEntry.precursorMass + maxPtmMass, buildIndexObj.getMaxPeptideMass());

            if (leftMass >= rightMass) {
                return;
            }

            Connection sqlConnection = DriverManager.getConnection(sqlPath);
            Statement sqlStatement = sqlConnection.createStatement();
            ResultSet sqlResultSet = sqlStatement.executeQuery(String.format("SELECT peptideMass, sequence, peptideCode, codeNormSquare, isTarget, leftFlank, rightFlank FROM peptideTable WHERE peptideMass BETWEEN %f AND %f", leftMass, rightMass));
            while (sqlResultSet.next()) {
                float peptideMass = sqlResultSet.getFloat(1);
                String peptide = sqlResultSet.getString(2);
                SparseBooleanVector theoVector = SparseBooleanVector.toSparseBooleanVector(sqlResultSet.getString(3));
                double theoNormSquare = sqlResultSet.getDouble(4);
                boolean isTarget = sqlResultSet.getBoolean(5);
                char leftFlank = sqlResultSet.getString(6).charAt(0);
                char rightFlank = sqlResultSet.getString(7).charAt(0);

                double score = 0;
                double temp1 = Math.sqrt(theoNormSquare * scanNormSquare);
                if (temp1 > 1e-6) {
                    score = theoVector.dot(scanCode) / temp1;
                }
                float deltaMass = peptideMass - spectrumEntry.precursorMass; // caution: the order matters under ms1ToleranceUnit == 1 situation

                if (isTarget) {
                    if ((deltaMass <= rightTol) && (deltaMass >= -1 * leftTol)) {
                        // PTM-free
                        if (ptmFreeQueue.size() < rankNum) {
                            ptmFreeQueue.add(new ResultEntry(score, peptide, leftFlank, rightFlank, false, false));
                        } else {
                            if (score > ptmFreeQueue.peek().score) {
                                ptmFreeQueue.poll();
                                ptmFreeQueue.add(new ResultEntry(score, peptide, leftFlank, rightFlank, false, false));
                            }
                        }
                    }

                    if ((deltaMass > rightTol) || (deltaMass < -1 * leftTol)) {
                        // PTM-only
                        if (ptmOnlyQueue.size() < rankNum) {
                            ptmOnlyQueue.add(new ResultEntry(score, peptide, leftFlank, rightFlank, false, true));
                        } else {
                            if (score > ptmOnlyQueue.peek().score) {
                                ptmOnlyQueue.poll();
                                ptmOnlyQueue.add(new ResultEntry(score, peptide, leftFlank, rightFlank, false, true));
                            }
                        }
                    }
                } else {
                    if ((deltaMass <= rightTol) && (deltaMass >= -1 * leftTol)) {
                        // PTM-free
                        if (ptmFreeQueue.size() < rankNum) {
                            ptmFreeQueue.add(new ResultEntry(score, peptide, leftFlank, rightFlank, true, false));
                        } else {
                            if (score > ptmFreeQueue.peek().score) {
                                ptmFreeQueue.poll();
                                ptmFreeQueue.add(new ResultEntry(score, peptide, leftFlank, rightFlank, true, false));
                            }
                        }
                    }

                    if ((deltaMass > rightTol) || (deltaMass < -1 * leftTol)) {
                        // PTM-only
                        if (ptmOnlyQueue.size() < rankNum) {
                            ptmOnlyQueue.add(new ResultEntry(score, peptide, leftFlank, rightFlank, true, true));
                        } else {
                            if (score > ptmOnlyQueue.peek().score) {
                                ptmOnlyQueue.poll();
                                ptmOnlyQueue.add(new ResultEntry(score, peptide, leftFlank, rightFlank, true, true));
                            }
                        }
                    }
                }
            }
            sqlStatement.close();
            sqlConnection.close();
            // (new File(sqlPath)).delete();
        } catch (SQLException ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }

        if (!(ptmFreeQueue.isEmpty() && ptmOnlyQueue.isEmpty())) {
            ptmFreeResult = convertResult(ptmFreeQueue, massToolObj, maxMs2Charge);
            ptmOnlyResult = convertResult(ptmOnlyQueue, massToolObj, maxMs2Charge);

            if (PIPI.DEV) {
                try {
                    BufferedWriter writer = new BufferedWriter(new FileWriter("PTMFreeContainingGlobalScoreBound" + "." + spectrumEntry.scanNum + "." + spectrumEntry.precursorCharge + ".csv"));
                    writer.write("PTM_free_low,PTM_free_high,PTM_containing_low,PTM_containing_high\n");
                    if (!ptmFreeQueue.isEmpty()) {
                        double highScore = 0;
                        double lowScore = 1;
                        for (ResultEntry temp : ptmFreeQueue) {
                            if (temp.score > highScore) {
                                highScore = temp.score;
                            }
                            if (temp.score < lowScore) {
                                lowScore = temp.score;
                            }
                        }
                        writer.write(String.format("%f,%f,", lowScore, highScore));
                    } else {
                        writer.write("-,-,");
                    }
                    if (!ptmOnlyQueue.isEmpty()) {
                        double highScore = 0;
                        double lowScore = 1;
                        for (ResultEntry temp : ptmOnlyQueue) {
                            if (temp.score > highScore) {
                                highScore = temp.score;
                            }
                            if (temp.score < lowScore) {
                                lowScore = temp.score;
                            }
                        }
                        writer.write(String.format("%f,%f\n", lowScore, highScore));
                    } else {
                        writer.write("-,-\n");
                    }
                    writer.close();

                    writer = new BufferedWriter(new FileWriter("ptm_only_candidates" + "." + spectrumEntry.scanNum + "." + spectrumEntry.precursorCharge + ".csv"));
                    writer.write("peptide,globalRank,is_decoy\n");
                    for (Peptide peptide : ptmOnlyResult) {
                        if (peptide.isDecoy()) {
                            writer.write(peptide.getPTMFreeSeq() + "," + peptide.getGlobalRank() + ",1\n");
                        } else {
                            writer.write(peptide.getPTMFreeSeq() + "," + peptide.getGlobalRank() + ",0\n");
                        }
                    }
                    writer.close();
                } catch (IOException ex) {
                    ex.printStackTrace();
                    System.exit(1);
                }
            }
        }
    }

    private List<Peptide> convertResult(PriorityQueue<ResultEntry> inputQueue, MassTool massToolObj, int maxMs2Charge) {
        List<Peptide> peptideList = new LinkedList<>();
        int globalRank = inputQueue.size();
        while (!inputQueue.isEmpty()) {
            ResultEntry temp = inputQueue.poll();
            peptideList.add(new Peptide(temp.peptide, temp.isDecoy(), massToolObj, maxMs2Charge, temp.score, temp.leftFlank, temp.rightFlank, globalRank));
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
