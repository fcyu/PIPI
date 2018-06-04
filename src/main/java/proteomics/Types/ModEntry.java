package proteomics.Types;

public class ModEntry {

    public final String id;
    public final String name;
    public final double mass;

    public ModEntry(String id, String name, double mass) {
        this.id = id;
        this.name = name;
        this.mass = mass;
    }
}
