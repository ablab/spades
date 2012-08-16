import utils.*;
import java.io.*;
import java.util.*;
import java.math.*;
import javax.swing.*;
import java.awt.*;
import java.awt.Graphics;
import static java.lang.Math.*;


public class PlotFPR implements Runnable {

	private static String filename = "";

    private final static String output_dir = "data/output/";

	private static double Threshold;

    private static boolean output = true;

    private class Pair implements Comparable<Pair> {
        int a;
        int b;

        public Pair(int x, int y) {
            a = x;
            b = y;
        }

        public int compareTo(Pair p) {
            if (a == p.a) 
                return (b - p.b);
            else 
                return (a - p.a);
        }

        public String toString() {
            return a + " " + b;
        }

    }

	public static void main(String[] args) {
        if (args.length > 0 && args[0].equals("-s")) 
            output = false;
		new Thread(new PlotFPR()).start();
	}

	private void debug(Object obj) {
		System.out.println(obj);
	}

    private static final String dir_name = "data/input/";


	public void run() {
		try {
			MyScanner in_tp, in_fpr, in_fnr, in_et, in_cl;
			Locale.setDefault(Locale.US);
			in_tp = new MyScanner(dir_name + "tp.prd");
			in_fpr = new MyScanner(dir_name + "fp.prd");
			in_fnr = new MyScanner(dir_name + "fn.prd");
			in_cl = new MyScanner(dir_name + "clustered.prd");
			in_et = new MyScanner(dir_name + "etalon.prd");
			PrintWriter out = new PrintWriter(System.out);
            int size_fn = 0;
            double maxfnr = -1;
			while (in_fnr.hasMoreTokens()) {
				int a = in_fnr.nextInt();
				int b = in_fnr.nextInt();

				double x = in_fnr.nextDouble();
				double y = in_fnr.nextDouble();
				double z = in_fnr.nextDouble();
                maxfnr = Math.max(maxfnr, x);
                in_fnr.nextToken();
                if (x > 1e-9) size_fn++;
            }
            in_fnr.close();
            int size_et = 0;
			//in_et.nextToken();
			while (in_et.hasMoreTokens()) {
				int a = in_et.nextInt();
				int b = in_et.nextInt();

				double x = in_et.nextDouble();
				double y = in_et.nextDouble();
				double z = in_et.nextDouble();
				in_et.nextToken();
                if (x > 1e-9) size_et++;
            }
            in_et.close();

            int size_cl = 0;
			//in_cl.nextToken();
			while (in_cl.hasMoreTokens()) {
				int a = in_cl.nextInt();
				int b = in_cl.nextInt();

				double x = in_cl.nextDouble();
				double y = in_cl.nextDouble();
				double z = in_cl.nextDouble();
				in_cl.nextToken();
                if (x > 1e-9) size_cl++;
            }
            in_cl.close();
//          getting fpr info
            
            double maxfpr = -1;
			int fpr = 0;
			int tp = 0;
            int size_fp = 0;
            StringBuffer buf = new StringBuffer("");
            double weight = -1;
            double total = 0;
            double cur = 0;
            int ind = 0;
            TreeMap<Double, Double> fpr_total = new TreeMap<Double, Double>();
            TreeMap<Double, Double> fnr_total = new TreeMap<Double, Double>();
            TreeSet<Double> all_positives = new TreeSet<Double>();
            all_positives.add(0.);
            fpr_total.put(Double.MAX_VALUE, 0.);
            fnr_total.put(Double.MAX_VALUE, 0.);
            
			while (in_fpr.hasMoreTokens()) {
				int a = in_fpr.nextInt();
				int b = in_fpr.nextInt();
				double x = in_fpr.nextDouble();
				double y = in_fpr.nextDouble();
				double z = in_fpr.nextDouble();
				in_fpr.nextToken();
                maxfpr = Math.max(maxfpr, x);
                if (x == 0) 
                    continue;
                if (weight > y + 1e-9) {
                    total += cur;
                    cur = 0;
                    fpr_total.put(weight, total);
                    all_positives.add(weight);
                    weight = y;
                } else if (weight < -1e-9) {
                    cur = 0;
                    weight = y;
                }
                cur++;
                size_fp++;
            }
            total += cur;
            cur = 0;
            fpr_total.put(weight, total);
            all_positives.add(weight);
            weight = 0;
            fpr_total.put(0.0, total);
            
            weight = -1;
            cur = 0;
            total = 0;
            ind = 0;
            int size_tp = 0;
			while (in_tp.hasMoreTokens()) {
				int a = in_tp.nextInt();
				int b = in_tp.nextInt();
				double x = in_tp.nextDouble();
				double y = in_tp.nextDouble();
				double z = in_tp.nextDouble();
				in_tp.nextToken();
                if (x == 0) 
                    continue;
                if (weight > y + 1e-9) {
                    total += cur;
                    cur = 0;
                    fnr_total.put(weight, total);
                    all_positives.add(weight);
                    weight = y;
                } else if (weight < 0) {
                    cur = 0;
                    weight = y;
                }
                cur++;
                size_tp++;
            }
            total += cur;
            cur = 0;
            fnr_total.put(weight, total);
            all_positives.add(weight);
            weight = 0;
            fnr_total.put(0.0, total);

            if (!output) {
                System.out.println("False positive rate now is " + size_fp * 1. / size_cl);  
                System.out.println("False negative rate now is " + size_fn * 1. / size_et);   
            }
            StringBuffer buffer_fpr = new StringBuffer("");
            StringBuffer buffer_fnr = new StringBuffer("");

            double lastKey = Math.max(maxfpr, maxfnr);
            double intersectionPoint = -1.;

            for (double thr : all_positives) {
                //System.out.println("Threshold " + thr);
                double size_fpr_ = fpr_total.ceilingEntry(thr).getValue();
                double size_fnr = size_fn + size_tp - fnr_total.ceilingEntry(thr).getValue();
                //System.out.println("Sizes are " + size_fpr_ + " " + size_fnr);
                double fpr_for_threshold = 0.;
                if (size_fpr_ < 1e-9) 
                    fpr_for_threshold = 0.;
                else 
                    fpr_for_threshold = (size_fpr_ * 100.) / (size_fpr_ + fnr_total.ceilingEntry(thr).getValue());
                
                double fnr_for_threshold = (size_fnr * 100.) / size_et;

                buffer_fpr.append(thr + " " + fpr_for_threshold + "\n");
                buffer_fnr.append(thr + " " + fnr_for_threshold + "\n");


                if (fpr_for_threshold > fnr_for_threshold) 
                    intersectionPoint = fpr_for_threshold;
            }

            System.out.println("Intersection of graphics is at fpr = fnr = " + intersectionPoint);

            PrintWriter out_fnr = new PrintWriter(output_dir + "plot/plot_fnr.out");
            out_fnr.println(buffer_fnr);
            out_fnr.close();
            PrintWriter out_fpr = new PrintWriter(output_dir + "plot/plot_fpr.out");
            out_fpr.println(buffer_fpr);
            out_fpr.close();
            
            PrintWriter out_config = new PrintWriter(output_dir + "plot/fpnr_plot.conf");
                        String text =
                        "#!/usr/bin/gnuplot -persist\n" + 
                        "set term x11 0\n" +
                        "plot \"plot_fpr.out\" with linespoints, \"plot_fnr.out\" with linespoints\n" + 
                        "pause -1 \"press any key to continue\"\n";
            out_config.print(text);
            out_config.close();
            out.close();
            in_fpr.close();
            in_tp.close();

		} catch(Exception e) {
			e.printStackTrace();
		}
	}


}
