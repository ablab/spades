import java.io.*;
import java.util.*;
import java.math.*;
import javax.swing.*;
import java.awt.*;
import java.awt.Graphics;


public class PlotFPR implements Runnable{
	
	private static String filename = "";

	private static double Threshold;

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
            int weight = 100000;
            int total = 0;
            int cur = 0;
            TreeMap<Integer, Integer> map = new TreeMap<Integer, Integer>();
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
                    map.put(weight, total);
                }
                cur++;
                size_fpr++;
            }
            for (Map.Entry<Integer, Integer> v : map){
                int w = v.second();
                buf.append(v.first + " " + (-w + size_fpr)*100.0 / size_fpr + "\n");
            }
            weight = 100000;
            cur = 0;
            total = 0;
            int size_tp = 0;
            map = new TreeMap<Integer, Integer>();
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
                    map.put(weight, total);
                }
                cur++;
                size_tp++;
            }
            for (Map.Entry<Integer, Integer> v : map){
                int w = v.second();
                buf.append(v.first + " " + (-w + size_tp)*100.0 / size_tp + "\n");
            }
            out.println("With threshold = " + Threshold + ":");
            out.println("False positive rate is going to be " + (fpr*100.0/(fpr + tp)));
            if (fpr <= size_fpr) out.println("FPR will be decreased by " + ( -fpr + size_fpr)*100.0 / size_fpr + " percent");
            if (tp <= size_tp) out.println("Perfect match will be decreased by " + (-tp + size_tp)*100.0 / size_tp + " percent");
            
            PrintWriter out1 = new PrintWriter(folder1 + "plot.conf");
                        String text =
                        "#!/usr/bin/gnuplot -persist\n" + 
                        "set term x11 0\n" +
                        "plot \"unclustered.prd\" with linespoints, \"clustered.prd\" with impulses," + 
                        "\"fpr.prd\" with points lt 1 lc 4 pt 7 ps 2," + " \"fnr.prd\" with points lt 1 lc 3 pt 7 ps 2 \n" + 
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
