package ru.spbau.bioinf.mgra;


public enum Direction {
    MINUS("minus"), PLUS("plus");

    private String text;

    Direction(String text) {
        this.text = text;
    }

    public static Direction getDirection(char ch) {
        return ch == '-' ? MINUS : PLUS;
    }

    @Override
    public String toString() {
        return text;
    }
}
