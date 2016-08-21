package proteomics.Types;

public class ExpAA implements Comparable<ExpAA> {

    private final char aa;
    private final float headLocation;
    private final float tailLocation;
    private final float headIntensity;
    private final float tailIntensity;
    private final float totalHalfIntensity;
    private int theoLocation; // starts from 0
    private String toString;
    private int hashCode;

    public ExpAA(char aa, float headLocation, float tailLocation, float headIntensity, float tailIntensity, int theoLocation) {
        this.aa = aa;
        this.headLocation = headLocation;
        this.tailLocation = tailLocation;
        this.headIntensity = headIntensity;
        this.tailIntensity = tailIntensity;
        this.totalHalfIntensity = (headIntensity + tailIntensity) / 2;
        this.theoLocation = -1;
        toString = headLocation + "." + aa + "." + theoLocation + "." + tailLocation;
        hashCode = toString.hashCode();
    }

    public char getAA() {
        return aa;
    }

    protected void setTheoLocation(int theo) {
        theoLocation = theo;
        // update toString and hashCode
        toString = headLocation + "." + aa + "." + theoLocation + "." + tailLocation;
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

    public String toString() {
        return toString;
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

    public ExpAA clone() {
        return new ExpAA(aa, headLocation, tailLocation, headIntensity, tailIntensity, theoLocation);
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