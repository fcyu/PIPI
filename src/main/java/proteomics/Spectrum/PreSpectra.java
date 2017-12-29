package proteomics.Spectrum;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.PIPI;
import proteomics.TheoSeq.MassTool;
import uk.ac.ebi.pride.tools.jmzreader.*;
import uk.ac.ebi.pride.tools.jmzreader.model.*;
import uk.ac.ebi.pride.tools.mgf_parser.model.Ms2Query;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.Statement;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class PreSpectra {

    private static final Logger logger = LoggerFactory.getLogger(PreSpectra.class);
    private static final Pattern scanNumPattern1 = Pattern.compile("Scan:([0-9]+) ", Pattern.CASE_INSENSITIVE);
    private static final Pattern scanNumPattern2 = Pattern.compile("scan=([0-9]+)", Pattern.CASE_INSENSITIVE);
    private static final Pattern scanNumPattern3 = Pattern.compile("^[^.]+\\.([0-9]+)\\.[0-9]+\\.[0-9]");
    private static final int[] isotopeCorrectionArray = new int[]{-2, -1, 0}; // do not change it

    private final float ms1Tolerance;
    private final int ms1ToleranceUnit;
    private final IsotopeDistribution isotopeDistribution;

    private Map<Integer, TreeMap<Integer, TreeSet<DevEntry>>> scanDevEntryMap = new HashMap<>();

    public PreSpectra(JMzReader spectraParser, Map<String, String> parameterMap, MassTool massToolObj, String ext, Set<Integer> msLevelSet, String sqlPath) throws Exception {
        int minMs1Charge = Integer.valueOf(parameterMap.get("min_ms1_charge"));
        int maxMs1Charge = Integer.valueOf(parameterMap.get("max_ms1_charge"));
        int minPeakNum = Integer.valueOf(parameterMap.get("min_peak_num"));
        ms1Tolerance = Float.valueOf(parameterMap.get("ms1_tolerance"));
        ms1ToleranceUnit = Integer.valueOf(parameterMap.get("ms1_tolerance_unit"));
        isotopeDistribution = new IsotopeDistribution(massToolObj.elementTable, 0, massToolObj.getLabeling());

        // prepare SQL database
        Connection sqlConnection = DriverManager.getConnection(sqlPath);
        Statement sqlStatement = sqlConnection.createStatement();
        sqlStatement.executeUpdate("DROP TABLE IF EXISTS spectraTable");
        sqlStatement.executeUpdate("CREATE TABLE spectraTable (scanNum INTEGER NOT NULL, scanId TEXT PRIMARY KEY, precursorCharge INTEGER NOT NULL, precursorMass REAL NOT NULL, mgfTitle TEXT NOT NULL, isotopeCorrectionNum INTEGER NOT NULL, ms1PearsonCorrelationCoefficient REAL NOT NULL, labelling TEXT, peptide TEXT, theoMass REAL, isDecoy INTEGER, globalRank INTEGER, normalizedCorrelationCoefficient REAL, score REAL, deltaLC REAL, deltaC REAL, matchedPeakNum INTEGER, ionFrac REAL, matchedHighestIntensityFrac REAL, explainedAaFrac REAL, ptmSupportingPeakFrac REAL, otherPtmPatterns TEXT, ptmDeltaScore TEXT)");
        sqlStatement.close();

        PreparedStatement sqlPrepareStatement = sqlConnection.prepareStatement("INSERT INTO spectraTable (scanNum, scanId, precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient) VALUES (?, ?, ?, ?, ?, ?, ?)");
        sqlConnection.setAutoCommit(false);
        int usefulSpectraNum = 0;

        Iterator<Spectrum> spectrumIterator = spectraParser.getSpectrumIterator();
        String parentId = null;
        while (spectrumIterator.hasNext()) {
            Spectrum spectrum = spectrumIterator.next();

            if (!msLevelSet.contains(spectrum.getMsLevel())) {
                parentId = spectrum.getId();
                continue;
            }

            if (spectrum.getPeakList().size() < minPeakNum) {
                logger.debug("Scan {} doesn't contain enough peak number ({}). Skip.", spectrum.getId(), minPeakNum);
                continue;
            }

            int scanNum = -1;
            float precursorMz = spectrum.getPrecursorMZ().floatValue();
            int precursorCharge = -1;
            float precursorMass = -1;
            int isotopeCorrectionNum = 0;
            double pearsonCorrelationCoefficient = -1;
            String mgfTitle = "";
            TreeMap<Integer, TreeSet<DevEntry>> chargeDevEntryMap = new TreeMap<>();
            if (ext.toLowerCase().contentEquals("mgf")) {
                mgfTitle = ((Ms2Query) spectrum).getTitle();
                Matcher matcher1 = scanNumPattern1.matcher(mgfTitle);
                Matcher matcher2 = scanNumPattern2.matcher(mgfTitle);
                Matcher matcher3 = scanNumPattern3.matcher(mgfTitle);
                if (matcher1.find()) {
                    scanNum = Integer.valueOf(matcher1.group(1));
                } else if (matcher2.find()) {
                    scanNum = Integer.valueOf(matcher2.group(1));
                } else if (matcher3.find()) {
                    scanNum = Integer.valueOf(matcher3.group(1));
                } else {
                    throw new Exception("Cannot get scan number from the MGF title " + mgfTitle + ". Please report your MGF title to fyuab@connect.ust.hk.");
                }

                if (PIPI.debugScanNumArray.length > 0) {
                    if (Arrays.binarySearch(PIPI.debugScanNumArray, scanNum) < 0) {
                        continue;
                    }
                }

                if (spectrum.getPrecursorCharge() == null) {
                    throw new Exception("MGF file does not contain charge information.");
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
                    for (int charge = minMs1Charge; charge <= maxMs1Charge; ++charge) {
                        Entry entry = getIsotopeCorrectionNum(precursorMz, charge, parentPeakList, chargeDevEntryMap);
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
                    Entry entry = getIsotopeCorrectionNum(precursorMz, precursorCharge, parentPeakList, chargeDevEntryMap);
                    if (entry.pearsonCorrelationCoefficient >= 0.7) { // If the Pearson correlation coefficient is smaller than 0.7, there is not enough evidence to change the original precursor mz.
                        isotopeCorrectionNum = entry.isotopeCorrectionNum;
                        pearsonCorrelationCoefficient = entry.pearsonCorrelationCoefficient;
                    }
                    precursorMass = (precursorMz - MassTool.PROTON) * precursorCharge + isotopeCorrectionNum * MassTool.C13_DIFF;
                }
            }

            if ((precursorCharge < minMs1Charge) || (precursorCharge > maxMs1Charge) || (precursorMass < 400)) {
                continue;
            }

            sqlPrepareStatement.setInt(1, scanNum);
            sqlPrepareStatement.setString(2, spectrum.getId());
            sqlPrepareStatement.setInt(3, precursorCharge);
            sqlPrepareStatement.setFloat(4, precursorMass);
            sqlPrepareStatement.setString(5, mgfTitle);
            sqlPrepareStatement.setInt(6, isotopeCorrectionNum);
            sqlPrepareStatement.setDouble(7, pearsonCorrelationCoefficient);
            sqlPrepareStatement.executeUpdate();
            ++usefulSpectraNum;

            if (PIPI.DEV) {
                scanDevEntryMap.put(scanNum, chargeDevEntryMap);
            }
        }
        sqlConnection.commit();
        sqlConnection.setAutoCommit(true);
        sqlPrepareStatement.close();
        sqlConnection.close();
        logger.info("Useful MS/MS spectra number: {}.", usefulSpectraNum);
    }

    public Map<Integer, TreeMap<Integer, TreeSet<DevEntry>>> getScanDevEntryMap() {
        return scanDevEntryMap;
    }

    private Entry getIsotopeCorrectionNum(double precursorMz, int charge, TreeMap<Double, Double> parentPeakList, TreeMap<Integer, TreeSet<DevEntry>> chargeDevEntryMap) {
        Entry entry = new Entry(0, 0);
        double leftTol = ms1Tolerance * 2;
        double rightTol = ms1Tolerance * 2;
        if (ms1ToleranceUnit == 1) {
            leftTol = (precursorMz - precursorMz / (1 + ms1Tolerance * 1e-6)) * 2;
            rightTol = (precursorMz / (1 - ms1Tolerance * 1e-6) - precursorMz) * 2;
        }
        for (int isotopeCorrectionNum : isotopeCorrectionArray) {
            double[][] expMatrix = new double[3][2];
            for (int i = 0; i < 3; ++i) {
                expMatrix[i][0] = precursorMz + (isotopeCorrectionNum + i) * MassTool.C13_DIFF / charge;
                NavigableMap<Double, Double> subMap = parentPeakList.subMap(expMatrix[i][0] - leftTol, true, expMatrix[i][0] + rightTol, true);
                for (double intensity : subMap.values()) {
                    expMatrix[i][1] = Math.max(expMatrix[i][1], intensity);
                }
            }

            if (Math.abs(expMatrix[0][1]) > 1) { // In bottom-up proteomics, the precursor mass won't be so large that we cannot observe the monoisotopic peak.
                Map<String, Integer> elementMap = isotopeDistribution.getElementMapFromMonoMass((expMatrix[0][0] - MassTool.PROTON) * charge);
                List<IsotopeDistribution.Peak> theoIsotopeDistribution = isotopeDistribution.calculate(elementMap);
                double pearsonCorrelationCoefficient = scaleAndCalPearsonCorrelationCoefficient(expMatrix, theoIsotopeDistribution, charge, isotopeCorrectionNum);
                if (pearsonCorrelationCoefficient > entry.pearsonCorrelationCoefficient) {
                    entry.pearsonCorrelationCoefficient = pearsonCorrelationCoefficient;
                    entry.isotopeCorrectionNum = isotopeCorrectionNum;
                }
                if (PIPI.DEV) {
                    double[][] theoMatrix = new double[expMatrix.length][2];
                    for (int i = 0; i < expMatrix.length; ++i) {
                        IsotopeDistribution.Peak peak = theoIsotopeDistribution.get(i);
                        theoMatrix[i][0] = peak.mass / charge + MassTool.PROTON;
                        theoMatrix[i][1] = peak.realArea;
                    }
                    if (chargeDevEntryMap.containsKey(charge)) {
                        chargeDevEntryMap.get(charge).add(new DevEntry(isotopeCorrectionNum, pearsonCorrelationCoefficient, expMatrix, theoMatrix));
                    } else {
                        TreeSet<DevEntry> tempSet = new TreeSet<>();
                        tempSet.add(new DevEntry(isotopeCorrectionNum, pearsonCorrelationCoefficient, expMatrix, theoMatrix));
                        chargeDevEntryMap.put(charge, tempSet);
                    }
                }
            }
        }
        return entry;
    }

    private double scaleAndCalPearsonCorrelationCoefficient(double[][] expMatrix, List<IsotopeDistribution.Peak> theoIsotopeDistribution, int precursorCharge, int isotopeCorrection) {
        // get theo peaks.
        int peakNum = Math.min(expMatrix.length, theoIsotopeDistribution.size());
        double[][] theoMatrix = new double[peakNum][2];
        for (int i = 0; i < peakNum; ++i) {
            IsotopeDistribution.Peak peak = theoIsotopeDistribution.get(i);
            theoMatrix[i][0] = peak.mass / precursorCharge + MassTool.PROTON;
            theoMatrix[i][1] = peak.realArea;
        }

        // scale theo peaks
        double scale = expMatrix[-1 * isotopeCorrection][1] / theoMatrix[-1 * isotopeCorrection][1];
        if (Math.abs(scale) > 1e-6) {
            for (int i = 0; i < peakNum; ++i) {
                theoMatrix[i][1] *= scale;
            }
            // calculate Pearson correlation coefficient.
            double theoIntensityMean = 0;
            double expIntensityMean = 0;
            for (int i = 0; i < peakNum; ++i) {
                theoIntensityMean += theoMatrix[i][1];
                expIntensityMean += expMatrix[i][1];
            }
            theoIntensityMean /= peakNum;
            expIntensityMean /= peakNum;
            double temp1 = 0;
            double temp2 = 0;
            double temp3 = 0;
            for (int i = 0; i < peakNum; ++i) {
                temp1 += (theoMatrix[i][1] - theoIntensityMean) * (expMatrix[i][1] - expIntensityMean);
                temp2 += Math.pow(theoMatrix[i][1] - theoIntensityMean, 2);
                temp3 += Math.pow(expMatrix[i][1] - expIntensityMean, 2);
            }
            return (temp1 == 0 || temp2 == 0) ? 0 : temp1 / (Math.sqrt(temp2 * temp3));
        } else {
            return 0;
        }
    }

    private class Entry {

        double pearsonCorrelationCoefficient;
        int isotopeCorrectionNum;

        Entry(double pearsonCorrelationCoefficient, int isotopeCorrectionNum) {
            this.pearsonCorrelationCoefficient = pearsonCorrelationCoefficient;
            this.isotopeCorrectionNum = isotopeCorrectionNum;
        }
    }

    public static class DevEntry implements Comparable<DevEntry> {

        public final int isotopeCorrectionNum;
        public final double pearsonCorrelationCoefficient;
        public final double[][] expMatrix;
        public final double[][] theoMatrix;
        private final String toString;
        private final int hashCode;

        public DevEntry(int isotopeCorrectionNum, double pearsonCorrelationCoefficient, double[][] expMatrix, double[][] theoMatrix) {
            this.isotopeCorrectionNum = isotopeCorrectionNum;
            this.pearsonCorrelationCoefficient = pearsonCorrelationCoefficient;
            this.expMatrix = expMatrix;
            this.theoMatrix = theoMatrix;
            toString = isotopeCorrectionNum + "-" + pearsonCorrelationCoefficient;
            hashCode = toString.hashCode();
        }

        public String toString() {
            return toString;
        }

        public int hashCode() {
            return hashCode;
        }

        public boolean equals(Object other) {
            if (other instanceof DevEntry) {
                DevEntry temp = (DevEntry) other;
                return isotopeCorrectionNum == temp.isotopeCorrectionNum && pearsonCorrelationCoefficient == temp.pearsonCorrelationCoefficient;
            } else {
                return false;
            }
        }

        public int compareTo(DevEntry other) {
            if (isotopeCorrectionNum > other.isotopeCorrectionNum) {
                return 1;
            } else if (isotopeCorrectionNum < other.isotopeCorrectionNum) {
                return -1;
            } else {
                return 0;
            }
        }
    }
}
