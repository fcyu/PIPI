package proteomics.Spectrum;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.ChargeMassTuple;
import proteomics.Types.SpectrumEntry;
import uk.ac.ebi.pride.tools.jmzreader.*;
import uk.ac.ebi.pride.tools.jmzreader.model.*;
import uk.ac.ebi.pride.tools.mgf_parser.model.Ms2Query;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.*;

public class PreSpectra {

    private static final Logger logger = LoggerFactory.getLogger(PreSpectra.class);

    private Map<Integer, ChargeMassTuple> numChargeMassMap = new HashMap<>();
    private Map<Integer, SpectrumEntry> numSpectrumMap = new HashMap<>();
    private TreeMap<Float, List<Integer>> massNumMap = new TreeMap<>();

    public PreSpectra(JMzReader spectraParser, Map<String, String> parameterMap, MassTool massToolObj, String ext, FilterSpectra.MassScan[] scanNumArray, int startIdx, int endIdx) {
        int minPeakNum = Integer.valueOf(parameterMap.get("min_peak_num"));
        float ms2Tolerance = Float.valueOf(parameterMap.get("ms2_tolerance"));
        Map<String, Float> massTable = massToolObj.returnMassTable();

        PreSpectrum preSpectrumObj = new PreSpectrum(massToolObj);
        PrintStream originalStream = System.out;
        PrintStream nullStream = new PrintStream(new OutputStream() {
            @Override
            public void write(int b) throws IOException {}
        });
        System.setOut(nullStream);

        int idx = startIdx;
        try {
            while ((idx < endIdx) && (idx < scanNumArray.length)) {
                Spectrum spectrum = spectraParser.getSpectrumById(scanNumArray[idx].scanId);
                int precursorCharge = spectrum.getPrecursorCharge();
                float precursorMz = spectrum.getPrecursorMZ().floatValue();
                float precursorMass = precursorMz * precursorCharge - precursorCharge * massTable.get("PROTON");
                Map<Double, Double> rawMzIntensityMap = spectrum.getPeakList();
                TreeMap<Float, Float> plMap = preSpectrumObj.preSpectrum(rawMzIntensityMap, precursorMass, precursorCharge, ms2Tolerance);
                if (plMap.size() <= minPeakNum) {
                    ++idx;
                    continue;
                }

                int scanNum = Integer.valueOf(spectrum.getId());
                try {
                    if (ext.contentEquals("mgf")) {
                        String title = ((Ms2Query) spectrum).getTitle();
                        String[] temp = title.split("\\.");
                        scanNum = Integer.valueOf(temp[temp.length - 2]);
                    }
                } catch (Exception ex) {}

                SpectrumEntry spectrumEntry = new SpectrumEntry(scanNum, precursorMz, precursorMass, precursorCharge, plMap, null);

                if (massNumMap.containsKey(precursorMass)) {
                    List<Integer> spectrumList = massNumMap.get(precursorMass);
                    spectrumList.add(scanNum);
                    massNumMap.put(precursorMass, spectrumList);
                } else {
                    List<Integer> scanNumList = new LinkedList<>();
                    scanNumList.add(scanNum);
                    massNumMap.put(precursorMass, scanNumList);
                }

                numSpectrumMap.put(scanNum, spectrumEntry);
                numChargeMassMap.put(scanNum, new ChargeMassTuple(precursorCharge, precursorMass));

                ++idx;
            }
        } catch (JMzReaderException ex) {
            logger.error(ex.getMessage());
            ex.printStackTrace();
            System.exit(1);
        }

        System.setOut(originalStream);
    }

    public Map<Integer, SpectrumEntry> returnNumSpectrumMap() {
        return numSpectrumMap;
    }

    public TreeMap<Float, List<Integer>> returnMassNumMap() {
        return massNumMap;
    }

    public Map<Integer, ChargeMassTuple> getNumChargeMassMap() {
        return numChargeMassMap;
    }
}
