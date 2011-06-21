package ru.spbau.bioinf.mgra;

import org.jdom.Element;

public class End {

    private int id;
    private EndType type;

    public End(String s) {
        if ("oo".equals(s)) {
            id = -1;
            type = EndType.OO;
        } else {
            int lastChar = s.length() - 1;
            id = Integer.parseInt(s.substring(0, lastChar));
            type = EndType.getType(s.charAt(lastChar));
        }
    }

    public int getId() {
        return id;
    }

    public Element toXml() {
        Element end = new Element("end");
        XmlUtil.addElement(end, "id", id > -1 ? Integer.toString(id) : "");
        XmlUtil.addElement(end, "type", type.toString());
        return end;
    }
}
