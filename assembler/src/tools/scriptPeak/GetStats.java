import java.io.*;
import java.util.*;
import java.math.*;
import javax.swing.*;
import java.awt.*;
import java.awt.Graphics;


public class GetStats implements Runnable{
	
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
		new Thread(new GetStats()).start();
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
			while (in_fpr.hasMoreTokens()){
				int a = in_fpr.nextInt();
				int b = in_fpr.nextInt();

				double x = in_fpr.nextDouble();
				double y = in_fpr.nextDouble();
				double z = in_fpr.nextDouble();
				in_fpr.nextToken();
                if (x == 0) continue;
			    if (y > Threshold - 1e-9) fpr++;
                size_fpr++;
            }
            int size_tp = 0;
			while (in_tp.hasMoreTokens()){
				int a = in_tp.nextInt();
				int b = in_tp.nextInt();
				double x = in_tp.nextDouble();
				double y = in_tp.nextDouble();
				double z = in_tp.nextDouble();
				in_tp.nextToken();
                if (x == 0) continue;
			    if (y > Threshold - 1e-9) tp++;
                size_tp++;
            }
            //debug(tp + " " + size_tp);
            out.println("With threshold = " + Threshold + ":");
            out.println("False positive rate is " + (size_fpr*100.0/(size_fpr + size_tp)));
            out.println("True positive rate is " + (size_tp*100.0/(size_fpr + size_tp)));
            out.println("False positive rate is going to be " + (fpr*100.0/(fpr + tp)));
            out.println("True positive rate is going to be" + (tp*100.0/(fpr + tp)));
            if (fpr <= size_fpr) out.println("FPR will be decreased by " + ( -fpr + size_fpr)*100.0 / size_fpr + " percent");
            if (tp <= size_tp) out.println("Perfect match will be decreased by " + (-tp + size_tp)*100.0 / size_tp + " percent");
            out.close();
            in_fpr.close();
            in_tp.close();

		}catch(Exception e){
			e.printStackTrace();
		}
	}

}
