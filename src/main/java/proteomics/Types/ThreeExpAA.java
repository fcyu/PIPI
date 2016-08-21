package proteomics.Types;


public class ThreeExpAA implements Comparable<ThreeExpAA> {

    private final ExpAA[] threeExpAa;
    private String toString;
    private int hashCode;
    private final String aaString;
    private final float totalIntensity;
    private int regionIdx;

    public ThreeExpAA(ExpAA aa1, ExpAA aa2, ExpAA aa3) {
        threeExpAa = new ExpAA[]{aa1, aa2, aa3};
        toString = threeExpAa[0].toString() + "-" + threeExpAa[1].toString() + "-" + threeExpAa[2].toString();
        hashCode = toString.hashCode();

        StringBuilder sb = new StringBuilder(5);
        for (ExpAA aa : threeExpAa) {
            sb.append(aa.getAA());
        }
        aaString = sb.toString();

        float intensity = threeExpAa[0].getHeadIntensity();
        for (ExpAA aa : threeExpAa) {
            intensity += aa.getTailIntensity();
        }
        totalIntensity = intensity;
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        return (other instanceof ThreeExpAA) && (this.hashCode() == other.hashCode());
    }

    public void setTheoLocation(int i, int theoLoc) {
        threeExpAa[i].setTheoLocation(theoLoc);
        // update toString and hashCode
        toString = threeExpAa[0].toString() + "-" + threeExpAa[1].toString() + "-" + threeExpAa[2].toString();
        hashCode = toString.hashCode();
    }

    public String toString() {
        return toString;
    }

    public int compareTo(ThreeExpAA other) {
        if (this.threeExpAa[0].getHeadLocation() > other.getExpAAs()[0].getHeadLocation()) {
            return 1;
        } else if (this.threeExpAa[0].getHeadLocation() < other.getExpAAs()[0].getHeadLocation()) {
            return -1;
        } else {
            return 0;
        }
    }

    public ExpAA[] getExpAAs() {
        return threeExpAa;
    }

    public String getAAString() {
        return aaString;
    }

    public float getTotalIntensity() {
        return totalIntensity;
    }

    public float getHeadLocation() {
        return threeExpAa[0].getHeadLocation();
    }

    public float getTailLocation() {
        return threeExpAa[threeExpAa.length - 1].getTailLocation();
    }

    public ThreeExpAA clone() {
        return new ThreeExpAA(threeExpAa[0].clone(), threeExpAa[1].clone(), threeExpAa[2].clone());
    }

    public int size() {
        return threeExpAa.length;
    }

    public ExpAA get(int i) {
        return threeExpAa[i];
    }

    public void setRegionIdx(int regionIdx) {
        this.regionIdx = regionIdx;
    }

    public int getRegionIdx() {
        return regionIdx;
    }
}
