package proteomics.Types;


import java.util.TreeMap;

public class PositionDeltaMassMap extends TreeMap<Coordinate, Float> {

    public PositionDeltaMassMap() {
        super();
    }

    public String toString() {
        StringBuilder sb = new StringBuilder(1000);
        for (Coordinate co : this.keySet()) {
            sb.append(this.get(co).toString());
            sb.append("@");
            sb.append(co.toString());
            sb.append(";");
        }
        return sb.toString();
    }

    public int hashCode() {
        return this.toString().hashCode();
    }

    public boolean equals(Object other) {
        if (other instanceof PositionDeltaMassMap) {
            PositionDeltaMassMap temp = (PositionDeltaMassMap) other;
            return temp.hashCode() == this.hashCode();
        } else {
            return false;
        }
    }

    public PositionDeltaMassMap clone() {
        PositionDeltaMassMap other = new PositionDeltaMassMap();
        other.clear();
        for (Coordinate co : this.keySet()) {
            other.put(co, this.get(co));
        }
        return other;
    }
}