package ru.spbau.bioinf.mgra;

public enum EndType {
    HEAD("h"), TAIL("t"), OO("oo");

    private String text;

    EndType(String text) {
        this.text = text;
    }

    public static EndType getType(char ch) {
        if (ch == 'h') return HEAD;
        if (ch == 't') return TAIL;
        return OO;
    }

    @Override
    public String toString() {
        return text;
    }
}
