import java.io.*;
import java.util.*;
import java.math.*;
import javax.swing.*;
import java.awt.*;
import java.awt.Graphics;


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
		if (args != null){
            Threshold = Double.parseDouble(args[0]);
		}
		new Thread(new PlotFPR()).start();
	}

	private void debug(Object obj){
		System.out.println(obj);
	}

	public void run(){
		try{
			MyScanner in_tp, in_fpr;
			
			Locale.setDefault(Locale.US);
			in_tp = new MyScanner("tp.prd");
			in_fpr = new MyScanner("fpr.prd");
			PrintWriter out = new PrintWriter(System.out);
            debug(Threshold);

//          getting fpr info
            
			int fpr = 0;
			int tp = 0;
            int size_fpr = 0;
            StringBuffer buf = new StringBuffer("");
            double weight = 100000;
            double total = 0;
            double cur = 0;
            int ind = 0;
            double[] fpr_weight = new double[N];
            double[] fpr_total = new double[N];
            double[] tp_weight = new double[N];
            double[] tp_total = new double[N];
			while (in_fpr.hasMoreTokens()){
				int a = in_fpr.nextInt();
				int b = in_fpr.nextInt();

				double x = in_fpr.nextDouble();
				double y = in_fpr.nextDouble();
				double z = in_fpr.nextDouble();
				in_fpr.nextToken();
                if (weight > y + 1e-9){
                    total += cur;
                    cur = 0;
                    fpr_weight[ind] = weight;
                    fpr_total[ind++] = total;
                    weight = y;
                }
                cur++;
                size_fpr++;
            }
            N = ind;
            for (int i = 0; i<N; i++){
                double w = fpr_total[i];
                buf.insert(0, fpr_weight[i] + " " + (-w + size_fpr)*100.0 / size_fpr + "\n");
            }
            PrintWriter out_fpr = new PrintWriter("plot_fpr.out");
            out_fpr.println(buf);
            out_fpr.close();
            buf = new StringBuffer("");
            weight = 100000;
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
                if (weight > y + 1e-9){
                    total += cur;
                    cur = 0;
                    tp_weight[ind] = weight;
                    tp_total[ind++] = total;
                    weight = y;
                }
                cur++;
                size_tp++;
            }
            int NN = ind;
            for (int i = 0; i<NN; i++){
                double w = tp_total[i];
                buf.insert(0, tp_weight[i] + " " + (-w + size_tp)*100.0 / size_tp + "\n");
            }
            PrintWriter out_tp = new PrintWriter("plot_tp.out");
            out_tp.println(buf);
            out_tp.close();
            //out.println("With threshold = " + Threshold + ":");
            //out.println("False positive rate is going to be " + (fpr*100.0/(fpr + tp)));
            //if (fpr <= size_fpr) out.println("FPR will be decreased by " + ( -fpr + size_fpr)*100.0 / size_fpr + " percent");
            //if (tp <= size_tp) out.println("Perfect match will be decreased by " + (-tp + size_tp)*100.0 / size_tp + " percent");
            
            PrintWriter out1 = new PrintWriter("plot_fpr.conf");
                        String text =
                        "#!/usr/bin/gnuplot -persist\n" + 
                        "set term x11 0\n" +
                        "plot \"plot_fpr.prd\" with linespoints, \"plot_tp.prd\" with linespoints,\n" + 
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
