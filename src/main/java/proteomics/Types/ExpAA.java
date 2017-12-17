package proteomics.Types;

public class ExpAA implements Comparable<ExpAA> {
    private final String aa;
    private final char ptmFreeAA;
    private final float headLocation;
    private final float tailLocation;
    private final float headIntensity;
    private final float tailIntensity;
    private final float totalHalfIntensity;
    private final float mod;
    private final float nTermMod;
    private final float cTermMod;
    private int theoLocation; // starts from 0, include N/C-terminal
    private int hashCode;

    public ExpAA(String aa, char ptmFreeAA, float headLocation, float tailLocation, float headIntensity, float tailIntensity, int theoLocation, float mod, float nTermMod, float cTermMod) {
        this.aa = aa;
        this.ptmFreeAA = ptmFreeAA;
        this.headLocation = headLocation;
        this.tailLocation = tailLocation;
        this.headIntensity = headIntensity;
        this.tailIntensity = tailIntensity;
        this.totalHalfIntensity = (headIntensity + tailIntensity) / 2;
        this.mod = mod;
        this.nTermMod = nTermMod;
        this.cTermMod = cTermMod;
        this.theoLocation = -1;
        String toString = headLocation + "." + aa + "." + theoLocation + "." + tailLocation;
        hashCode = toString.hashCode();
    }

    public String getAA() {
        return aa;
    }

    public char getPtmFreeAA() {
        return ptmFreeAA;
    }

    public float getMod() {
        return mod;
    }

    public float getnTermMod() {
        return nTermMod;
    }

    public float getcTermMod() {
        return cTermMod;
    }

    void setTheoLocation(int theo) {
        theoLocation = theo;
        // update toString and hashCode
        String toString = headLocation + "." + aa + "." + theoLocation + "." + tailLocation;
        hashCode = toString.hashCode();
    }

    public int getTheoLocation() {
        return theoLocation;
    }

    public int compareTo(ExpAA other) {
        if (headLocation > other.getHeadLocation()) {
            return 1;
        } else if (headLocation < other.getHeadLocation()) {
            return -1;
        } else {
            return 0;
        }
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        if (other instanceof ExpAA) {
            ExpAA temp = (ExpAA) other;
            return this.hashCode() == temp.hashCode();
        } else {
            return false;
        }
    }

    public boolean approximateEquals(ExpAA other, float tolerance) {
        return (this.aa.contentEquals(other.aa) && (this.theoLocation == other.theoLocation) && (Math.abs(this.headLocation - other.headLocation) <= tolerance));
    }

    public ExpAA clone() {
        return new ExpAA(aa, ptmFreeAA, headLocation, tailLocation, headIntensity, tailIntensity, theoLocation, mod, nTermMod, cTermMod);
    }

    public float getHeadLocation() {
        return headLocation;
    }

    public float getTailLocation() {
        return tailLocation;
    }

    public float getHeadIntensity() {
        return headIntensity;
    }

    public float getTailIntensity() {
        return tailIntensity;
    }

    public float getTotalHalfIntensity() {
        return totalHalfIntensity;
    }
}