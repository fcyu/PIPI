package proteomics.Types;

public class AA {

    public final char aa;
    public final float ptmDeltaMass;
    private final String toString;
    private final int hashCode;

    public AA(char aa, float ptmDeltaMass) {
        this.aa = aa;
        this.ptmDeltaMass = ptmDeltaMass;
        toString = aa + "-" + ptmDeltaMass;
        hashCode = toString.hashCode();
    }

    public String toString() {
        return toString;
    }

    public int hashCode() {
        return hashCode;
    }

    public AA clone() {
        return new AA(aa, ptmDeltaMass);
    }
}
