package proteomics.Spectrum;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.SpectrumEntry;
import uk.ac.ebi.pride.tools.jmzreader.*;
import uk.ac.ebi.pride.tools.jmzreader.model.*;
import uk.ac.ebi.pride.tools.mgf_parser.model.Ms2Query;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class PreSpectra {

    private static final Logger logger = LoggerFactory.getLogger(PreSpectra.class);
    private static final Pattern scanNumPattern1 = Pattern.compile("Scan:([0-9]+) ");
    private static final Pattern scanNumPattern2 = Pattern.compile("^[^.]+\\.([0-9]+)\\.[0-9]+\\.[0-9]");

    private Map<Integer, SpectrumEntry> numSpectrumMap = new HashMap<>();

    public PreSpectra(JMzReader spectraParser, Map<String, String> parameterMap, MassTool massToolObj, String ext, Set<Integer> msLevelSet) {
        int minMs1Charge = Integer.valueOf(parameterMap.get("min_ms1_charge"));
        int maxMs1Charge = Integer.valueOf(parameterMap.get("max_ms1_charge"));
        float minPrecursorMass =  Float.valueOf(parameterMap.get("min_precursor_mass"));
        float maxPrecursorMass = Float.valueOf(parameterMap.get("max_precursor_mass"));
        int minPeakNum = Integer.valueOf(parameterMap.get("min_peak_num"));
        float ms2Tolerance = Float.valueOf(parameterMap.get("ms2_tolerance"));
        float minClear = Float.valueOf(parameterMap.get("min_clear_mz"));
        float maxClear = Float.valueOf(parameterMap.get("max_clear_mz"));

        PreSpectrum preSpectrumObj = new PreSpectrum(massToolObj);
        PrintStream originalStream = System.out;
        PrintStream nullStream = new PrintStream(new OutputStream() {
            @Override
            public void write(int b) throws IOException {}
        });
        System.setOut(nullStream);

        Iterator<Spectrum> spectrumIterator = spectraParser.getSpectrumIterator();
        while (spectrumIterator.hasNext()) {
            Spectrum spectrum = spectrumIterator.next();

            if (!msLevelSet.contains(spectrum.getMsLevel())) {
                continue;
            }

            Map<Double, Double> rawMzIntensityMap = spectrum.getPeakList();
            if (rawMzIntensityMap.size() < minPeakNum) {
                logger.debug("Scan {} doesn't contain enough peak number ({}). Skip.", spectrum.getId(), minPeakNum);
                continue;
            }

            int scanNum = -1;
            String mgfTitle = "";
            try {
                if (ext.toLowerCase().contentEquals("mgf")) {
                    mgfTitle = ((Ms2Query) spectrum).getTitle();
                    Matcher matcher1 = scanNumPattern1.matcher(mgfTitle);
                    Matcher matcher2 = scanNumPattern2.matcher(mgfTitle);
                    if (matcher1.find()) {
                        scanNum = Integer.valueOf(matcher1.group(1));
                    } else if (matcher2.find()) {
                        scanNum = Integer.valueOf(matcher2.group(1));
                    } else {
                        logger.error("Cannot get scan number from the MGF title {}. PIPI only support the MGF files converted from ProteoWizard or ReAdw4Mascot4.", mgfTitle);
                        System.exit(1);
                    }
                } else {
                    scanNum = Integer.valueOf(spectrum.getId());
                }
            } catch (Exception ex) {
                ex.printStackTrace();
                logger.error(ex.getMessage());
                System.exit(1);
            }

            float precursorMz = spectrum.getPrecursorMZ().floatValue();

            SpectrumEntry spectrumEntry;
            int precursorCharge;
            float precursorMass;
            if (spectrum.getPrecursorCharge() == null) {
                logger.debug("Scan {} doesn't have charge information.", spectrum.getId());
                precursorCharge = 0;
                precursorMass = 0;
                TreeMap<Float, Float> plMap = preSpectrumObj.preSpectrum(rawMzIntensityMap, precursorMass, precursorCharge, ms2Tolerance, minClear, maxClear);
                if (plMap.size() <= minPeakNum) {
                    continue;
                }
                TreeMap<Float, Float> unprocessedPlMap = new TreeMap<>();
                for (double mz : rawMzIntensityMap.keySet()) {
                    if (mz >= 50) {
                        unprocessedPlMap.put((float) mz, rawMzIntensityMap.get(mz).floatValue());
                    }
                }
                spectrumEntry = new SpectrumEntry(scanNum, precursorMz, precursorMass, precursorCharge, plMap, unprocessedPlMap, mgfTitle);
            } else {
                precursorCharge = spectrum.getPrecursorCharge();
                if ((precursorCharge < minMs1Charge) || (precursorCharge > maxMs1Charge)) {
                    continue;
                }
                precursorMass = precursorMz * precursorCharge - precursorCharge * MassTool.PROTON;
                if ((precursorMass > maxPrecursorMass) || (precursorMass < minPrecursorMass)) {
                    continue;
                }
                TreeMap<Float, Float> plMap = preSpectrumObj.preSpectrum(rawMzIntensityMap, precursorMass, precursorCharge, ms2Tolerance, minClear, maxClear);
                if (plMap.size() <= minPeakNum) {
                    continue;
                }
                TreeMap<Float, Float> unprocessedPlMap = new TreeMap<>();
                for (double mz : rawMzIntensityMap.keySet()) {
                    if ((mz >= 50) && (mz <= precursorMass)) {
                        unprocessedPlMap.put((float) mz, rawMzIntensityMap.get(mz).floatValue());
                    }
                }
                spectrumEntry = new SpectrumEntry(scanNum, precursorMz, precursorMass, precursorCharge, plMap, unprocessedPlMap, mgfTitle);
            }

            numSpectrumMap.put(scanNum, spectrumEntry);
        }

        System.setOut(originalStream);
    }

    public Map<Integer, SpectrumEntry> returnNumSpectrumMap() {
        return numSpectrumMap;
    }
}
