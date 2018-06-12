package proteomics.Types;


public class Segment implements Comparable<Segment> {

    private final String segmentString;
    private final int hashCode;

    public Segment(String segmentString) {
        String temp = reverseString(segmentString);
        // Compare the string and its reversed version. Kept the smaller one.
        int compareResult = segmentString.compareTo(temp);
        if (compareResult > 0) {
            this.segmentString = temp;
        } else {
            this.segmentString = segmentString;
        }

        hashCode = segmentString.hashCode();
    }

    public int length() {
        return segmentString.length();
    }

    public  boolean equals(Object other) {
        return (other instanceof Segment) && this.segmentString.contentEquals(((Segment) other).segmentString);
    }

    public String toString() {
        return segmentString;
    }

    public int hashCode() {
        return hashCode;
    }

    public int compareTo(Segment other) {
        return this.segmentString.compareTo(other.segmentString);
    }

    private String reverseString(String seq) {
        return new StringBuilder(seq).reverse().toString();
    }
}
