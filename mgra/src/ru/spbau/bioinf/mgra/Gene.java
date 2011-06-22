package ru.spbau.bioinf.mgra;

import org.jdom.Element;

public class Gene {

    private int id;
    private Direction direction;
    private End end = null;

    public Gene(int id, Direction direction) {
        this.id = id;
        this.direction = direction;
    }

    public int getId() {
        return id;
    }

    public void setEnd(End end) {
        this.end = end;
    }

    public void clearEnd() {
        this.end = null;
    }

    public Element toXml() {
        Element gene = new Element("gene");
        XmlUtil.addElement(gene, "id", id);
        XmlUtil.addElement(gene, "direction", direction.toString());
        if (end !=null) {
            gene.addContent(end.toXml());
        }

        return gene;
    }
}
