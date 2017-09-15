package proteomics.Types;


public class VarModParam {

    public final float mass;
    public final char aa;
    public final int priority;

    private final String toString;
    private final int hashCode;

    public VarModParam(float mass, char aa, int priority) {
        this.mass = mass;
        this.aa = aa;
        this.priority = priority;

        toString = Math.round(mass * 100) + "@" + aa;
        hashCode = toString.hashCode();
    }

    public String toString() {
        return toString;
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
