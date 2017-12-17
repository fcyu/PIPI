package proteomics.Types;

public class Coordinate implements Comparable<Coordinate> {

    public final int x;
    public final int y;
    private final int hashCode;

    public Coordinate(int x, int y) {
        this.x = x;
        this.y = y;
        String toString = "(" + x + "-" + y + ")";
        hashCode = toString.hashCode();
    }

    public int compareTo(Coordinate other) {
        if (x > other.x) {
            return 1;
        } else if (x < other.x) {
            return -1;
        } else {
            if (y > other.y) {
                return 1;
            } else if (y < other.y) {
                return -1;
            } else {
                return 0;
            }
        }
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        if (other instanceof Coordinate) {
            Coordinate temp = (Coordinate) other;
            return temp.hashCode() == this.hashCode();
        } else {
            return false;
        }
    }
}
