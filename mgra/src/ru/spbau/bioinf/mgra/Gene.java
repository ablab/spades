package ru.spbau.bioinf.mgra;

import org.jdom.Element;

public class Gene {

    private int id;
    private Direction direction;

    public Gene(int id, Direction direction) {
        this.id = id;
        this.direction = direction;
    }

    public int getId() {
        return id;
    }

    public Element toXml() {
        Element gene = new Element("gene");
        XmlUtil.addElement(gene, "id", id);
        XmlUtil.addElement(gene, "direction", direction.toString());
        return gene;
    }
}
