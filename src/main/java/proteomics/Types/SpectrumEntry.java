package proteomics.Types;

import java.util.*;

public class SpectrumEntry {
    public final int scanNum;
    public final float precursorMz;
    public final float precursorMass;
    public final int precursorCharge;
    public final TreeMap<Float, Float> plMap;
    private TreeMap<Float, Float> plMapWithVirtualPeaks = null;
    public final SparseVector plMapXcorr;
    private final String toString;

    public SpectrumEntry(int scanNum, float precursorMz, float precursorMass, int precursorCharge, TreeMap<Float, Float> plMap, SparseVector plMapXcorr) {
        this.scanNum = scanNum;
        this.precursorMz = precursorMz;
        this.precursorMass = precursorMass;
        this.precursorCharge = precursorCharge;
        this.plMap = plMap;
        this.plMapXcorr = plMapXcorr;
        toString = this.scanNum + " (charge = " + this.precursorCharge + ", mass = " + this.precursorMass + ", peak_num = " + this.plMap.size() + ")";
    }

    public String toString() {
        return toString;
    }

    public void setPlMapWithVirtualPeaks(TreeMap<Float, Float> input) {
        plMapWithVirtualPeaks = input;
    }

    public TreeMap<Float, Float> getPlMapWithVirtualPeaks() {
        return plMapWithVirtualPeaks;
    }
}