package ru.spbau.bioinf.mgra;

import org.jdom.Element;

import java.util.ArrayList;
import java.util.List;

public class Tree {

    private String root = "";

    private Tree parent;

    private List<Tree> children = new ArrayList<Tree>();

    private int depth = 0;
    private int width = 0;

    @Override
    public String toString() {
        return root;
    }

    public Tree(Tree parent, String s) {
        this.parent = parent;
        for (int i = 0;  i < s.length(); i++) {
             char ch = s.charAt(i);
             if (Character.isLetter(ch)) {
                 root += ch;
             }
        }
        if (s.length() < 2)
            return;
        if (s.startsWith("(")) {
            s = s.substring(1, s.length() - 1);
        }
        int stat = 0;
        int last = 0;
        for (int i = 0;  i < s.length(); i++) {
            char ch = s.charAt(i);
            if (ch == '(') {
                stat++;
            }
            if (ch == ')') {
                stat--;
            }
            if (stat == 0 && ch != ',') {
                children.add(new Tree(this, s.substring(last, i+1)));
                last = i+2;
            }
        }
    }

    public Element toXml() {
        Element tree = new Element("tree");
        int maxDepth = getDepth();
        Element row = new Element("row");
        addCell(row, this, 1);
        tree.addContent(row);


        for (int level = 1; level < maxDepth; level++) {
            row = new Element("row");
            for (Tree child : children) {
                child.addCells(row, level, maxDepth - level);
            }
            tree.addContent(row);

        }
        return tree;
    }

    private void addCell(Element row, Tree t, int height) {
        Element cell = new Element("cell");
        XmlUtil.addElement(cell, "width", t.getWidth());
        XmlUtil.addElement(cell, "height", children.size() > 0 ? 1 : height);
        XmlUtil.addElement(cell, "text", t.root);
        row.addContent(cell);
    }

    public void addCells(Element row, int level, int depth) {
        if (level == 1) {
            addCell(row, this, depth);
        } else {
            for (Tree child : children) {
                child.addCells(row, level - 1, depth);
            }
        }
    }
    public int getDepth() {
        if (depth == 0) {
            depth = 1;
            for (Tree child : children) {
                int cd = child.getDepth();
                if (cd >= depth)
                    depth = cd + 1;
            }
        }
        return depth;
    }

    public int getWidth() {
        if (width == 0) {
            if (children.size() == 0) {
                width = 1;
            } else {
                for (Tree child : children) {
                    width += child.getWidth();
                }
            }
        }
        return width;
    }
}
