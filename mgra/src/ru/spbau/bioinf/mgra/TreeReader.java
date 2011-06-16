package ru.spbau.bioinf.mgra;

import org.jdom.Document;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;

public class TreeReader {

    public static void main(String[] args) throws Exception {
        String fileName = "data/mam6.cfg";
        if (args.length > 0) {
            fileName = args[0];
        }
        new TreeReader(new File(fileName));
    }

    public TreeReader(File cfg) throws IOException{
        BufferedReader input = getBufferedInputReader(cfg);
        String s;
        while (!(s = input.readLine()).startsWith("[Trees]")) {}
        Tree tree = new Tree(null, input.readLine());
        Document doc = new Document();
        doc.setRootElement(tree.toXml(cfg.getParentFile()));
        XmlUtil.saveXml(doc, "tree.xml");
    }

    public static BufferedReader getBufferedInputReader(File file) throws FileNotFoundException {
        return new BufferedReader(new InputStreamReader(new FileInputStream(file)));
    }
}
