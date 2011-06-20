package ru.spbau.bioinf.mgra;

import org.jdom.Element;

import java.util.ArrayList;
import java.util.List;

public class Chromosome {

    List<Gene> genes = new ArrayList<Gene>();

    public Chromosome(String s) {
        String[] data = s.split(" ");
        for (String v : data) {
           if (!v.startsWith("$")) {
               Gene gene = new Gene(Integer.parseInt(v.substring(1)), Direction.getDirection(v.charAt(0)));
               genes.add(gene);
           }
        }
    }

    public Element toXml() {
        Element chr = new Element("chromosome");
        for (Gene gene : genes) {
            chr.addContent(gene.toXml());
        }
        return chr;
    }
}
