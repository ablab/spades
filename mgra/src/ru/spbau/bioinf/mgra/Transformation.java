package ru.spbau.bioinf.mgra;

import org.jdom.Element;

import java.util.ArrayList;
import java.util.List;

public class Transformation {
    private List<Chromosome> chromosomes = new ArrayList<Chromosome>();
    private Element before = new Element("before");
    private Element after = new Element("after");

    End[] ends = new End[4];

    public Transformation(String text) {
        String[] data = text.split("[ \t]");
        for (int i = 0; i < ends.length; i++) {
            ends[i] = new End(data[i]);
        }
    }

    public void update(Genome genome) {
        List<Chromosome> all = genome.getChromosomes();
        for (Chromosome chromosome : all) {
            for (End end : ends) {
                if (chromosome.contains(end)) {
                    this.chromosomes.add(chromosome);
                    before.addContent(chromosome.toXml());
                    break;
                }
            }
        }
    }

    public Element toXml() {
        Element tr = new Element("transformation");
        tr.addContent(before);
        for (End end : ends) {
            tr.addContent(end.toXml());
        }
        tr.addContent(after);
        return tr;
    }
}
