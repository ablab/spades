import java.io.*;
import java.util.*;
import java.math.*;
import javax.swing.*;
import java.awt.*;
import java.awt.Graphics;
import static java.lang.Math.*;


public class PlotFPR implements Runnable{
	
	private static String filename = "";

	private static double Threshold;

    private int N = 100000;

    private class Pair implements Comparable<Pair>{
        int a;
        int b;

        public Pair(int x, int y){
            a = x;
            b = y;
        }

        public int compareTo(Pair p){
              if (a == p.a) return (b - p.b);
              return (a - p.a);
        }

        public String toString(){
              
            return a + " " + b;
        }

    }

	public static void main(String[] args){
		new Thread(new PlotFPR()).start();
	}

	private void debug(Object obj){
		System.out.println(obj);
	}

	public void run(){
		try{
			MyScanner in_tp, in_fpr, in_fnr, in_et;
			String folder0 = "data_zero/";
			Locale.setDefault(Locale.US);
			in_tp = new MyScanner(folder0 + "tp.prd");
			in_fpr = new MyScanner(folder0 + "fpr.prd");
			in_fnr = new MyScanner(folder0 + "fn.prd");
			in_et = new MyScanner("distance_estimation_cl.prd");
			PrintWriter out = new PrintWriter(System.out);
            //debug(Threshold);
            int size_fnr = 0;
			in_fnr.nextToken();
			while (in_fnr.hasMoreTokens()){
				int a = in_fnr.nextInt();
				int b = in_fnr.nextInt();

				double x = in_fnr.nextDouble();
				double y = in_fnr.nextDouble();
				double z = in_fnr.nextDouble();
				in_fnr.nextToken();
                if (!(abs(x) < 1e-9 && abs(y) < 1e-9)) size_fnr++;
            }
            in_fnr.close();
            int size_et = 0;
			in_et.nextToken();
			while (in_et.hasMoreTokens()){
				int a = in_et.nextInt();
				int b = in_et.nextInt();

				double x = in_et.nextDouble();
				double y = in_et.nextDouble();
				double z = in_et.nextDouble();
				in_et.nextToken();
                if (!(abs(x) < 1e-9 && abs(y) < 1e-9)) size_et++;
            }
            in_et.close();
//          getting fpr info
            
			int fpr = 0;
			int tp = 0;
            //int size_fpr = 0;
            StringBuffer buf = new StringBuffer("");
            double weight = -1;
            double total = 0;
            double cur = 0;
            int ind = 0;
            TreeMap<Double, Double> fpr_total = new TreeMap<Double, Double>();
            TreeMap<Double, Double> fnr_total = new TreeMap<Double, Double>();
			while (in_fpr.hasMoreTokens()){
				int a = in_fpr.nextInt();
				int b = in_fpr.nextInt();

				double x = in_fpr.nextDouble();
				double y = in_fpr.nextDouble();
				double z = in_fpr.nextDouble();
				in_fpr.nextToken();
                if (x == 0) continue;
                if (weight > y + 1e-9){
                    total += cur;
                    cur = 0;
                    fpr_total.put(weight, total);
                    weight = y;
                }else if (weight<0){
                    cur = 0;
                    weight = y;
                }
                cur++;
                //size_fpr++;
            }
            total += cur;
            fpr_total.put(0.0, total);
            N = ind;
            
            weight = -1;
            cur = 0;
            total = 0;
            ind = 0;
            int size_tp = 0;
			while (in_tp.hasMoreTokens()){
				int a = in_tp.nextInt();
				int b = in_tp.nextInt();
				double x = in_tp.nextDouble();
				double y = in_tp.nextDouble();
				double z = in_tp.nextDouble();
				in_tp.nextToken();
                if (x == 0) continue;
                if (weight > y + 1e-9){
                    total += cur;
                    cur = 0;
                    fnr_total.put(weight, total);
                    weight = y;
                }else if (weight<0){
                    cur = 0;
                    weight = y;
                }
                cur++;
                size_tp++;
            }
            total += cur;
            fnr_total.put(0.0, total);
            Set<Double> fpr_set = fpr_total.keySet();
            for (double w : fpr_set){
                double size_fpr = fpr_total.get(w);
            
                double size_tpv = 0;
                if (w > fnr_total.lastKey() + 1e-9) break;
                else size_tpv =  fnr_total.ceilingEntry(w).getValue();
                buf.append(w + " " + (size_fpr)*100.0 / (size_tpv + size_fpr) + "\n");
            }
            PrintWriter out_fpr = new PrintWriter("plot_fpr.out");
            out_fpr.println(buf);
            out_fpr.close();
            buf = new StringBuffer("");
            Set<Double> fnr_set = fnr_total.keySet();
            for (double w : fnr_set){
                double size_fpr = 0;
                if ( w > fpr_total.lastKey() + 1e-9) break;
                else size_fpr = fpr_total.ceilingEntry(w).getValue();
                double size_tpv =  fnr_total.get(w);
                double size_fn = size_fnr + size_tp - size_tpv;
                //debug(size_tp - size_tpv);
                buf.append(w + " " + (size_fn)*100.0 / (size_et) + "\n");
            }
            PrintWriter out_tp = new PrintWriter("plot_fnr.out");
            out_tp.println(buf);
            out_tp.close();
            //out.println("With threshold = " + Threshold + ":");
            //out.println("False positive rate is going to be " + (fpr*100.0/(fpr + tp)));
            //if (fpr <= size_fpr) out.println("FPR will be decreased by " + ( -fpr + size_fpr)*100.0 / size_fpr + " percent");
            //if (tp <= size_tp) out.println("Perfect match will be decreased by " + (-tp + size_tp)*100.0 / size_tp + " percent");
            
            PrintWriter out1 = new PrintWriter("plot_fpnr.conf");
                        String text =
                        "#!/usr/bin/gnuplot -persist\n" + 
                        "set term x11 0\n" +
                        "plot \"plot_fpr.out\" with linespoints, \"plot_fnr.out\" with linespoints\n" + 
                        "pause -1 \"press any key to continue\"\n";
            out1.print(text);
            out1.close();
            out.close();
            in_fpr.close();
            in_tp.close();

		}catch(Exception e){
			e.printStackTrace();
		}
	}

}
