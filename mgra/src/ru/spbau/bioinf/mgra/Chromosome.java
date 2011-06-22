package ru.spbau.bioinf.mgra;

import org.jdom.Element;

import java.util.ArrayList;
import java.util.List;

public class Chromosome {

    private int id;

    private List<Gene> genes = new ArrayList<Gene>();

    public Chromosome(int id, String s) {
        this.id = id;
        String[] data = s.split(" ");
        for (String v : data) {
           if (!v.startsWith("$")) {
               Gene gene = new Gene(Integer.parseInt(v.substring(1)), Direction.getDirection(v.charAt(0)));
               genes.add(gene);
           }
        }
    }

    public boolean contains(End end) {
        for (Gene gene : genes) {
            if (gene.getId() == end.getId()) {
                return true;
            }
        }
        return false;
    }

    public void mark(End end) {
        for (Gene gene : genes) {
            if (gene.getId() == end.getId()) {
                gene.setEnd(end);
            }
        }
    }

    public void clearEnds() {
        for (Gene gene : genes) {
            gene.clearEnd();
        }
    }

    public Element toXml() {
        Element chr = new Element("chromosome");
        XmlUtil.addElement(chr, "id", id);
        for (Gene gene : genes) {
            chr.addContent(gene.toXml());
        }
        return chr;
    }
}
