package proteomics.Types;

public class AA {

    public final char aa;
    public final double ptmDeltaMass;
    private final int hashCode;

    public AA(char aa, double ptmDeltaMass) {
        this.aa = aa;
        this.ptmDeltaMass = ptmDeltaMass;
        String toString = aa + "-" + ptmDeltaMass;
        hashCode = toString.hashCode();
    }

    public boolean hasMod() {
        return Math.abs(ptmDeltaMass) > 0.1;
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        if (other instanceof AA) {
            return ((AA) other).hashCode == hashCode;
        } else {
            return false;
        }
    }

    public AA clone() {
        return new AA(aa, ptmDeltaMass);
    }
}
