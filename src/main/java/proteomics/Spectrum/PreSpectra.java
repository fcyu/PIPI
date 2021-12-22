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

package proteomics.Spectrum;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.PIPI;
import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.IsotopeDistribution;
import static ProteomicsLibrary.Utilities.*;
import uk.ac.ebi.pride.tools.jmzreader.*;
import uk.ac.ebi.pride.tools.jmzreader.model.*;
import uk.ac.ebi.pride.tools.mgf_parser.model.Ms2Query;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.Statement;
import java.util.*;

public class PreSpectra {

    private static final Logger logger = LoggerFactory.getLogger(PreSpectra.class);
    public static final int topN = 6;

    private final IsotopeDistribution isotopeDistribution;

    private int usefulSpectraNum = 0;

    public PreSpectra(JMzReader spectraParser, double ms1Tolerance, int ms1ToleranceUnit, MassTool massTool, String ext, Set<Integer> msLevelSet, String sqlPath) throws Exception {
        isotopeDistribution = new IsotopeDistribution(massTool.getElementTable(), 0, massTool.getLabelling());

        // prepare SQL database
        Connection sqlConnection = DriverManager.getConnection(sqlPath);
        Statement sqlStatement = sqlConnection.createStatement();
        sqlStatement.executeUpdate("PRAGMA journal_mode=WAL");
        sqlStatement.executeUpdate("DROP TABLE IF EXISTS spectraTable");
        sqlStatement.executeUpdate("CREATE TABLE spectraTable (scanNum INTEGER NOT NULL, scanId TEXT PRIMARY KEY, precursorCharge INTEGER NOT NULL, precursorMass REAL NOT NULL, mgfTitle TEXT NOT NULL, isotopeCorrectionNum INTEGER NOT NULL, ms1PearsonCorrelationCoefficient REAL NOT NULL, labelling TEXT, peptide TEXT, theoMass REAL, isDecoy INTEGER, globalRank INTEGER, normalizedCorrelationCoefficient REAL, score REAL, deltaLCn REAL, deltaCn REAL, matchedPeakNum INTEGER, ionFrac REAL, matchedHighestIntensityFrac REAL, explainedAaFrac REAL, otherPtmPatterns TEXT, aScore TEXT)");
        sqlStatement.close();

        PreparedStatement sqlPrepareStatement = sqlConnection.prepareStatement("INSERT INTO spectraTable (scanNum, scanId, precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient) VALUES (?, ?, ?, ?, ?, ?, ?)");
        sqlConnection.setAutoCommit(false);

        Iterator<Spectrum> spectrumIterator = spectraParser.getSpectrumIterator();
        String parentId = null;
        while (spectrumIterator.hasNext()) {
            try {
                Spectrum spectrum = spectrumIterator.next();

                if (ext.toLowerCase().contentEquals("mzxml")) {
                    if (!msLevelSet.contains(spectrum.getMsLevel())) {
                        parentId = spectrum.getId();
                        continue;
                    }
                }

                if (spectrum.getPeakList().size() < 5) {
                    continue;
                }

                int scanNum;
                double precursorMz = spectrum.getPrecursorMZ();
                int precursorCharge = -1;
                double precursorMass;
                int isotopeCorrectionNum = 0;
                double pearsonCorrelationCoefficient = -1;
                String mgfTitle = "";
                if (ext.toLowerCase().contentEquals("mgf")) {
                    mgfTitle = ((Ms2Query) spectrum).getTitle();
                    scanNum = getScanNum(mgfTitle);

                    if (PIPI.debugScanNumArray.length > 0) {
                        if (Arrays.binarySearch(PIPI.debugScanNumArray, scanNum) < 0) {
                            continue;
                        }
                    }

                    if (spectrum.getPrecursorCharge() == null) {
                        logger.warn("Scan {} does not contain charge information.", scanNum);
                        continue;
                    } else {
                        precursorCharge = spectrum.getPrecursorCharge();
                        precursorMass = precursorMz * precursorCharge - precursorCharge * MassTool.PROTON;
                    }
                } else {
                    scanNum = Integer.valueOf(spectrum.getId());

                    if (PIPI.debugScanNumArray.length > 0) {
                        if (Arrays.binarySearch(PIPI.debugScanNumArray, scanNum) < 0) {
                            continue;
                        }
                    }

                    TreeMap<Double, Double> parentPeakList = new TreeMap<>(spectraParser.getSpectrumById(parentId).getPeakList());
                    if (spectrum.getPrecursorCharge() == null) {
                        // We have to infer the precursor charge.
                        for (int charge = 2; charge <= 4; ++charge) {
                            IsotopeDistribution.Entry entry = isotopeDistribution.getIsotopeCorrectionNum(precursorMz, ms1Tolerance, ms1ToleranceUnit, charge, parentPeakList);
                            if (entry.pearsonCorrelationCoefficient > pearsonCorrelationCoefficient) {
                                pearsonCorrelationCoefficient = entry.pearsonCorrelationCoefficient;
                                isotopeCorrectionNum = entry.isotopeCorrectionNum;
                                precursorCharge = charge;
                            }
                        }
                        if (precursorCharge > 0) {
                            precursorMass = (precursorMz - MassTool.PROTON) * precursorCharge + isotopeCorrectionNum * MassTool.C13_DIFF;
                        } else {
                            logger.warn("Cannot infer the precursor charge for scan {}.", scanNum);
                            continue;
                        }
                    } else {
                        // We do not try to correct the precursor charge if there is one.
                        precursorCharge = spectrum.getPrecursorCharge();
                        IsotopeDistribution.Entry entry = isotopeDistribution.getIsotopeCorrectionNum(precursorMz, ms1Tolerance, ms1ToleranceUnit, precursorCharge, parentPeakList);
                        if (entry.pearsonCorrelationCoefficient >= 0.7) { // If the Pearson correlation coefficient is smaller than 0.7, there is not enough evidence to change the original precursor mz.
                            isotopeCorrectionNum = entry.isotopeCorrectionNum;
                            pearsonCorrelationCoefficient = entry.pearsonCorrelationCoefficient;
                        }
                        precursorMass = (precursorMz - MassTool.PROTON) * precursorCharge + isotopeCorrectionNum * MassTool.C13_DIFF;
                    }
                }

                sqlPrepareStatement.setInt(1, scanNum);
                sqlPrepareStatement.setString(2, spectrum.getId());
                sqlPrepareStatement.setInt(3, precursorCharge);
                sqlPrepareStatement.setDouble(4, precursorMass);
                sqlPrepareStatement.setString(5, mgfTitle);
                sqlPrepareStatement.setInt(6, isotopeCorrectionNum);
                sqlPrepareStatement.setDouble(7, pearsonCorrelationCoefficient);
                sqlPrepareStatement.executeUpdate();
                ++usefulSpectraNum;
            } catch (RuntimeException ex) {
                logger.error(ex.toString());
            }
        }
        sqlConnection.commit();
        sqlConnection.setAutoCommit(true);
        sqlPrepareStatement.close();
        sqlConnection.close();
        logger.info("Useful MS/MS spectra number: {}.", usefulSpectraNum);
    }

    public int getUsefulSpectraNum() {
        return usefulSpectraNum;
    }
}
