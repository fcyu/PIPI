package proteomics.Types;


public class VarModParam {

    public final float modMass;
    public final char aa;
    public final boolean proteinTerminal;
    public final String toString;

    public VarModParam(float modMass, char aa, boolean proteinTerminal) {
        this.modMass = modMass;
        this.aa = aa;
        this.proteinTerminal = proteinTerminal;
        toString = modMass + "@" + aa + proteinTerminal;
    }

    public String toString() {
        return toString;
    }

    public int hashCode() {
        return toString.hashCode();
    }

    public boolean equals(Object other) {
        if (other instanceof VarModParam) {
            VarModParam temp = (VarModParam) other;
            return (Math.abs(temp.modMass - modMass) <= 0.01) && (temp.aa == aa) && (temp.proteinTerminal == proteinTerminal);
        } else {
            return false;
        }
    }
}
