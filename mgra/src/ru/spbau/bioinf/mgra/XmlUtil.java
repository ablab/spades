package ru.spbau.bioinf.mgra;

import org.apache.log4j.Logger;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.output.Format;
import org.jdom.output.XMLOutputter;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Set;

public class XmlUtil {

    private static final Logger log = Logger.getLogger(XmlUtil.class);

    public static XMLOutputter outputter = new XMLOutputter(Format.getPrettyFormat());

    public static Element addElement(Element element, String tag) {
        element.addContent(new Element(tag));
        return element;
    }

    public static Element addElement(Element element, String tag, String value) {
        element.addContent(new Element(tag).addContent(value));
        return element;
    }

    public static Element addElement(Element element, String tag, StringBuilder value) {
        element.addContent(new Element(tag).addContent(value.toString()));
        return element;
    }

    public static Element addElement(Element element, String tag, char value) {
        element.addContent(new Element(tag).addContent(Character.toString(value)));
        return element;
    }

    public static Element addElement(Element element, String tag, int value) {
        return  addElement(element, tag, Integer.toString(value));
    }

    public static Element addElement(Element element, String tag, boolean value) {
        return  addElement(element, tag, Boolean.toString(value));
    }

    public static Element addElement(Element element, String tag, float value) {
        return  addElement(element, tag, Float.toString(value));
    }

    public static Element addElement(Element element, String tag, double value) {
        return  addElement(element, tag, Double.toString(value));
    }

    public static void addTags(Element prsm, Set<String> tags) {
        for (String tag : tags) {
            Element tagElement = new Element("tag");
            addElement(tagElement, "value", tag);
            prsm.addContent(tagElement);
        }
    }

    public static void saveXml(Document doc, String fileName) throws IOException {
        FileWriter writer = null;
        try {
            writer = new FileWriter(fileName);
            outputter.output(doc, writer);
        } catch (IOException e) {
            log.error("Error saving xml", e);
            throw new RuntimeException(e);
        } finally {
            if (writer != null)
                writer.close();
        }
    }
}
