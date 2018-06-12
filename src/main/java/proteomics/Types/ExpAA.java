package proteomics.Types;

public class ExpAA implements Comparable<ExpAA> {
    private final String aa;
    private final char ptmFreeAA;
    private final double headLocation;
    private final double tailLocation;
    private final double headIntensity;
    private final double tailIntensity;
    private final double mod;
    private final double nTermMod;
    private final double cTermMod;
    private int theoLocation; // starts from 0, include N/C-terminal
    private int hashCode;

    public ExpAA(String aa, char ptmFreeAA, double headLocation, double tailLocation, double headIntensity, double tailIntensity, int theoLocation, double mod, double nTermMod, double cTermMod) {
        this.aa = aa;
        this.ptmFreeAA = ptmFreeAA;
        this.headLocation = headLocation;
        this.tailLocation = tailLocation;
        this.headIntensity = headIntensity;
        this.tailIntensity = tailIntensity;
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

    public double getMod() {
        return mod;
    }

    public double getnTermMod() {
        return nTermMod;
    }

    public double getcTermMod() {
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
        if (headLocation > other.headLocation) {
            return 1;
        } else if (headLocation < other.headLocation) {
            return -1;
        } else {
            return Double.compare(tailIntensity, other.tailIntensity);
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

    public boolean approximateEquals(ExpAA other, double tolerance) {
        return (this.aa.contentEquals(other.aa) && (this.theoLocation == other.theoLocation) && (Math.abs(this.headLocation - other.headLocation) <= tolerance));
    }

    public ExpAA clone() throws CloneNotSupportedException {
        super.clone();
        return new ExpAA(aa, ptmFreeAA, headLocation, tailLocation, headIntensity, tailIntensity, theoLocation, mod, nTermMod, cTermMod);
    }

    double getHeadLocation() {
        return headLocation;
    }

    double getTailLocation() {
        return tailLocation;
    }

    double getHeadIntensity() {
        return headIntensity;
    }

    double getTailIntensity() {
        return tailIntensity;
    }
}