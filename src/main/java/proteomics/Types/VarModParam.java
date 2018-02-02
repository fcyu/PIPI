package proteomics.Types;


public class VarModParam {

    public final double mass;
    public final char aa;
    public final int priority; // 1 = high; 0 = low.
    public final boolean onlyProteinTerminalIfnc;

    private final int hashCode;

    public VarModParam(double mass, char aa, int priority, boolean onlyProteinTerminalIfnc) {
        this.mass = mass;
        this.aa = aa;
        this.priority = priority;
        this.onlyProteinTerminalIfnc = onlyProteinTerminalIfnc;

        String toString = Math.round(mass * 1000) + "@" + aa;
        hashCode = toString.hashCode();
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        if (other instanceof VarModParam) {
            VarModParam temp = (VarModParam) other;
            return temp.hashCode == hashCode;
        } else {
            return false;
        }
    }
}
