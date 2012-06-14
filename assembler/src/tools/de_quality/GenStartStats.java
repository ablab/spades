import java.io.*;
import java.util.*;
import java.math.*;
import static java.lang.Math.*;
import utils.*;

public class GenStartStats implements Runnable {

    private class Pair<T extends Comparable> /*implements Comparable<Pair<T>>*/ {
        private T first;
        private T second;

        Pair(T x, T y) {
            first = x;
            second = y;
        }
        
        public T first() {
            return first;   
        }
        
        public T second() {
            return second;   
        }

        @Override
        public int hashCode() {
            return first.hashCode() + second.hashCode();
        }

        public boolean equals(Object o) {
            Pair<T> another = ((Pair<T>) o);
            return (first.compareTo(another.first()) == 0) && (second.compareTo(another.second()) == 0);   
        }
        
        public String toString() {
            if (first instanceof Double) {
                return String.format("%.2f %.2f", (Double) first, (Double) second);   
            } else 
                return first.toString() + " " + second.toString();
        }
    }
    
    private class Three<T extends Comparable> {
        private T first;
        private T second;
        private T third;

        Three(T x, T y, T z) {
            first = x;
            second = y;
            third = z;
        }
        
        public T first() {
            return first;   
        }
        
        public T second() {
            return second;   
        }

        public T third() {
            return third;   
        }

        public boolean equals(Object o) {
            Three<T> another = ((Three<T>) o);
            return (first.compareTo(another.first()) == 0);// && (second.compareTo(another.second()) == 0) && (third.compareTo(another.third()) == 0);   
        }

        @Override
        public int hashCode() {
            return first.hashCode();// + second.hashCode() + third.hashCode();
        }

        public String toString() {
            if (first instanceof Double) {
                return String.format("%.2f %.2f %.2f", (Double) first, (Double) second, (Double) third);   
            } else 
                return first.toString() + " " + second.toString() + " " + third.toString();
        }
    }

    public void run() {
        try {
            Locale.setDefault(Locale.US);
            MyScanner etalon;
            MyScanner clustered;
            etalon = new MyScanner("etalon.prd");
            clustered = new MyScanner("clustered.prd");
            
            PrintWriter fpr = new PrintWriter("test_fp.prd");
            PrintWriter tpr = new PrintWriter("test_tp.prd");

            HashMap<Pair<Integer>, HashSet<Three<Double>>> edges = new HashMap<Pair<Integer>, HashSet<Three<Double>>>(); 

            while (etalon.hasMoreTokens()) {
                Pair<Integer> edgePair = new Pair<Integer>(etalon.nextInt(), etalon.nextInt());
                Three<Double> estThree = new Three<Double>(etalon.nextDouble(), etalon.nextDouble(), etalon.nextDouble());
                etalon.nextToken();
                //if (round(estThree.second()) == 0) continue;
                if (round(estThree.first()) <= 0) continue;

                if (!edges.containsKey(edgePair)) {
                    HashSet<Three<Double>> set = new HashSet<Three<Double>>();
                    set.add(estThree);
                    edges.put(edgePair, set);
                } else {
                    //System.out.println("Warn! Contains already : " + edgePair + " " + estThree);
                    HashSet<Three<Double>> set = edges.get(edgePair);
                    set.add(estThree);

                    //edges.put(edgePair, set);
                }
            }

            while (clustered.hasMoreTokens()) {
                Pair<Integer> edgePair = new Pair<Integer>(clustered.nextInt(), clustered.nextInt());
                Three<Double> estThree = new Three<Double>(clustered.nextDouble(), clustered.nextDouble(), clustered.nextDouble());
                double est_dist = estThree.first();
                double est_weight = estThree.second();
                double est_var = estThree.third();
                clustered.nextToken();  
                if (round(est_dist) <= 0) continue;
                if (edges.containsKey(edgePair)) {
                    boolean okay = false;
                    for (Three<Double> three : edges.get(edgePair)) {
                        double dist = three.first();
                        double weight = three.second();
                        double var = three.third();
                        if (est_dist + est_var + var >= dist && est_dist - est_var - var <= dist) {
                            tpr.println(edgePair + " " + estThree + " ."); 
                            okay = true;
                            break;
                        }
                    }
                    if (!okay) 
                        fpr.println(edgePair + " " + estThree + " .");
                } else {
                    fpr.println(edgePair + " " + estThree + " .");   
                }


            }


            etalon.close();
            clustered.close();
            fpr.close();
            tpr.close();

        } 
        catch (Exception e) {
            e.printStackTrace();   
        }
    }


    public static void main(String[] args) { 
        new Thread(new GenStartStats()).start();
    }
}
