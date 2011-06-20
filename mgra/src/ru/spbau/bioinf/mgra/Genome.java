package ru.spbau.bioinf.mgra;

import org.jdom.Element;

import java.util.ArrayList;
import java.util.List;

public class Genome {
    List<Chromosome> chromosomes = new ArrayList<Chromosome>();

    public void addChromosome(Chromosome chr) {
        chromosomes.add(chr);
    }

    public Element toXml() {
        Element genome = new Element("genome");
        for (Chromosome chromosome : chromosomes) {
            genome.addContent(chromosome.toXml());
        }
        return genome;
    }
}
