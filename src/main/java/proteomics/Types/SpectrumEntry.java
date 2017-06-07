package proteomics.Types;

import java.util.*;

public class SpectrumEntry {
    public final int scanNum;
    public final float precursorMz;
    public final float precursorMass;
    public final int precursorCharge;
    public final TreeMap<Float, Float> plMap;
    public final TreeMap<Float, Float> unprocessedPlMap;
    private final String toString;

    public SpectrumEntry(int scanNum, float precursorMz, float precursorMass, int precursorCharge, TreeMap<Float, Float> plMap, TreeMap<Float, Float> unprocessedPlMap) {
        this.scanNum = scanNum;
        this.precursorMz = precursorMz;
        this.precursorMass = precursorMass;
        this.precursorCharge = precursorCharge;
        this.plMap = plMap;
        this.unprocessedPlMap = unprocessedPlMap;
        toString = this.scanNum + " (charge = " + this.precursorCharge + ", mass = " + this.precursorMass + ", peak_num = " + this.plMap.size() + ")";
    }

    public String toString() {
        return toString;
    }
}