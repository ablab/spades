import utils.*;
import java.io.*;
import java.util.*;
import java.math.*;
import javax.swing.*;
import java.awt.*;
import java.awt.Graphics;


public class GetErrors implements Runnable{
	
	private static String filename = "";

	private final static String input_dir = "data/input/";

	private final static String output_dir = "data/output/hists";

	private static int upperBound;

	private static int lowerBound;

	private static boolean symmetric = false;

    private static int Nmax = 100000000;
	
	private static int mask = 0; 

	//	2^0 = symmetric
    //	2^1 = non-symmetric
	//  2^2 = morethan lowerBound
	//	2^3 = lessthan upperBound
	//

    private class Pair implements Comparable<Pair> {
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

	public static void main(String[] args) {
		if (args != null){
			for (int i = 0; i<args.length; i++) if (args[i].equals("-s")) mask |= 1<<0;
			for (int i = 0; i<args.length; i++) if (args[i].equals("-ns")) mask |= 1<<1;
			for (int i = 0; i<args.length; i++) if (args[i].equals("-m")) {
				lowerBound = Integer.parseInt(args[i+1]);
				mask |= 1<<2;
			}
			for (int i = 0; i<args.length; i++) if (args[i].equals("-l")) {
				upperBound = Integer.parseInt(args[i+1]);
				mask |= 1<<3;
			}
			for (int i = 0; i<args.length; i++) if (args[i].equals("-f")) filename = args[i+1];
		}
		new Thread(new GetErrors()).start();
	}

	private void debug(Object obj) {
		System.out.println(obj);
	}

	private boolean isAcceptable(int len1, int len2) {
		int b = 0;
        if (len1 == len2) b |= 1;
        if (len1 != len2) b |= 2;
        if (len1 > lowerBound && len2 > lowerBound) b |= 4;
        if (len1 < upperBound && len2 < upperBound) b |= 8;
		return ((b & mask) == mask);
	}



	
	public void run() {
		try {
			MyScanner in, in2, incl, fnrin, fprin, inpaths;
			
			Locale.setDefault(Locale.US);
			if (filename.length() > 0) {
				int a = 0;
				while (a < filename.length() && filename.charAt(a)!='.') a++;
				filename = filename.substring(0, a);
				in = new MyScanner(input_dir + filename + ".prd");
				in2 = new MyScanner(input_dir + filename + ".grp");
                incl = new MyScanner(input_dir + filename + "_cl.prd");
                inpaths = new MyScanner(input_dir + "paths.prd");
                fnrin = new MyScanner(input_dir + "fn.prd");
                fprin = new MyScanner(input_dir + "fp.prd");
			}
			else throw new IOException("no input data");
			debug(filename);
			//output_dir = ((mask & 1) == 0) ? ("testfiles_" + filename + "/!sym/") : 
			//("testfiles_" + filename + "/sym/");
			File file = new File(output_dir);
			file.mkdirs();
			PrintWriter out = null;
			in2.nextInt();
			int[] edges = new int[Nmax];
			int ind = 0;
			while (in2.hasMoreTokens()) {
				String s = in2.nextToken();
				if (s.equals("Edge"))
                    ind = in2.nextInt();					
				else if (s.equals("=")) 
                    edges[ind-1] = in2.nextInt();
			}

            TreeSet<Pair> set = new TreeSet<Pair>();
            StringBuffer buf = new StringBuffer("");
			int len1 = 0;
			int len2 = 0;
			int e1 = 0;
			int e2 = 0;
			boolean filtering = false;

//          getting fpr info
            
            debug("Fpr info");
            buf = new StringBuffer("");
			len1 = 0;
		    len2 = 0;
			e1 = 0;
			e2 = 0;
			filtering = false;
			while (fprin.hasMoreTokens()) {
				int a = fprin.nextInt();
				int b = fprin.nextInt();
				double x = fprin.nextDouble();
				double y = fprin.nextDouble();
				double z = fprin.nextDouble();
				fprin.nextToken();
				if (!(a == e1 && b == e2)) {
					if (out != null && filtering) {
						out.println(buf);
						out.close();
						buf = new StringBuffer("");
					}	
					filtering = (isAcceptable(edges[a-1], edges[b-1]));
					if (filtering) {
                        set.add(new Pair(e1, e2));
						debug("Current edge is processing : " + a + " " + b + " " + edges[a-1] + " " + edges[b-1]);
                        String folder1 = output_dir + "/" + a + "_" + b + "_" + edges[a-1] + "_" + edges[b-1] + "/";
            			file = new File(folder1);
			            file.mkdirs();
						out = new PrintWriter(folder1 + "fpr.prd");
					}
				}
				if (filtering) if (x > 0) {
                    buf.append(x + " " + 0 + "\n");
                }
				e1 = a;
				e2 = b;
			}
		    if (out != null) {
                if (filtering) 
                    out.println(buf);
                out.close();
            }
            fprin.close();


//          getting fnr info
            debug("Fnr info");
            buf = new StringBuffer("");
			len1 = 0;
			len2 = 0;
			e1 = 0;
			e2 = 0;
            out = null;
			filtering = false;
			while (fnrin.hasMoreTokens()) {
				int a = fnrin.nextInt();
				int b = fnrin.nextInt();
				double x = fnrin.nextDouble();
				double y = fnrin.nextDouble();
				double z = fnrin.nextDouble();
				fnrin.nextToken();
				if (!(a == e1 && b == e2)) {
					if (out != null && filtering) {
						out.println(buf);
						out.close();
						buf = new StringBuffer("");
					}	
					filtering = (isAcceptable(edges[a-1], edges[b-1]));
					if (filtering) {

						//debug("Current edge is processing : " + a + " " + b + " " + edges[a-1] + " " + edges[b-1]);
                        set.add(new Pair(e1, e2));
                        String folder1 = output_dir + "/" + a + "_" + b + "_" + edges[a-1] + "_" + edges[b-1] + "/";
                        file = new File(folder1);
                        file.mkdirs();
						out = new PrintWriter(folder1 + "fnr.prd");
					}
				}
				if (filtering) if (x > 0) {
                    buf.append(x + " " + 0 + "\n");
                }
				e1 = a;
				e2 = b;
			}
		    if (out != null) {
                if (filtering) 
                    out.println(buf);
                out.close();
            }
            fnrin.close();
                        
            e1 = e2 = 0;
//          getting paired info
			debug("Paired Info");
			in.nextInt();
			while (in.hasMoreTokens()) {
				int a = in.nextInt(); //edge1
				int b = in.nextInt(); //edge2
				double x = in.nextDouble(); //distance
				double y = in.nextDouble(); //weight
				double z = in.nextDouble(); //variance
				in.nextToken();
				if (!(a == e1 && b == e2)) {
					if (out != null && filtering) {
						out.println(buf);
						out.close();
						buf = new StringBuffer("");
					}	
					filtering = (isAcceptable(edges[a - 1], edges[b - 1]));
					if (filtering) {

						if (!set.contains(new Pair(a, b))) 
                            continue;
                        String folder1 = output_dir + "/" + a + "_" + b + "_" + edges[a - 1] + "_" + edges[b - 1] + "/";
                        out = new PrintWriter(folder1 + "unclustered.prd");
                        //generating config
						debug("Current edge is processing : " + a + " " + b + " " + edges[a - 1] + " " + edges[b - 1]);
                        debug("Generating config files");
                        PrintWriter out1 = new PrintWriter(folder1 + "plot.conf");
                        String text =
                        "#!/usr/bin/gnuplot -persist\n" + 
                        "set term x11 0\n" +
                        "plot \"unclustered.prd\" with linespoints, \"clustered.prd\" with impulses," + 
                        "\"fpr.prd\" with points lt 1 lc 4 pt 7 ps 1," + " \"fnr.prd\" with points lt 1 lc 3 pt 7 ps 1, " + " \"paths.prd\" with impulses lt 1 lw 2 lc 5\n" +
                        "pause -1 \"press any key to continue\"\n";
                        out1.print(text);
                        out1.close();

                        out1 = new PrintWriter(folder1 + "png_plot.conf");
                        text =
                        "#!/usr/bin/gnuplot -persist\n" + 
                        "set term png enhanced size 1920, 1080 14\n" +
                        "plot \"unclustered.prd\" with linespoints, \"clustered.prd\" with impulses," + 
                        "\"fpr.prd\" with points lt 1 lc 4 pt 7 ps 1," + " \"fnr.prd\" with points lt 1 lc 3 pt 7 ps 1, " + " \"paths.prd\" with impulses lt 1 lw 2 lc 5\n" +
                        "set output \"" + a + "_" + b + "_" + edges[a - 1] + "_" + edges[b - 1] + ".png\"\n" +
                        "replot\n";
                        out1.print(text);
                        out1.close();
					}
				}
				if (filtering) if (x > 0. && y != 0.) {
                    buf.append(x + " " + y + "\n");
                }
				e1 = a;
				e2 = b;
			}
		    if (out != null) {
                if (filtering) 
                    out.println(buf);
                out.close();
            }
			in.close();

//          getting clustered info
            
            debug("Clustered info");
            buf = new StringBuffer("");
			incl.nextInt();
			len1 = 0;
		
			len2 = 0;
			e1 = 0;
			e2 = 0;
			filtering = false;
			while (incl.hasMoreTokens()) {
				int a = incl.nextInt();
				int b = incl.nextInt();
				double x = incl.nextDouble();
				double y = incl.nextDouble();
				double z = incl.nextDouble();
				incl.nextToken();
				if (!(a == e1 && b == e2)) {
					if (out != null && filtering) {
						out.println(buf);
						out.close();
						buf = new StringBuffer("");
					}	
					filtering = (isAcceptable(edges[a-1], edges[b-1]));
					if (filtering) {

						debug("Current edge is processing : " + a + " " + b + " " + edges[a-1] + " " + edges[b-1]);
						if (!set.contains(new Pair(a, b))) continue;
                        String folder1 = output_dir + "/" + a + "_" + b + "_" + edges[a-1] + "_" + edges[b-1] + "/";
						out = new PrintWriter(folder1 + "clustered.prd");
					}
				}
				if (filtering) if (x > 0.) {
                    buf.append(x + " " + y + "\n");
                }
				e1 = a;
				e2 = b;
			}
		    if (out != null) {
                if (filtering) 
                    out.println(buf);
                out.close();
            }
            incl.close();

//          getting paths info
            
            debug("Paths info");
            buf = new StringBuffer("");
			inpaths.nextInt();
			len1 = 0;
			len2 = 0;
			e1 = 0;
			e2 = 0;
			filtering = false;
			while (inpaths.hasMoreTokens()) {
				int a = inpaths.nextInt();
				int b = inpaths.nextInt();
				double x = inpaths.nextDouble();
				double y = inpaths.nextDouble();
				double z = inpaths.nextDouble();
				inpaths.nextToken();
				if (!(a == e1 && b == e2)) {
					if (out != null && filtering) {
						out.println(buf);
						out.close();
						buf = new StringBuffer("");
					}	
					filtering = (isAcceptable(edges[a-1], edges[b-1]));
					if (filtering) {

						debug("Current edge is processing : " + a + " " + b + " " + edges[a-1] + " " + edges[b-1]);
						if (!set.contains(new Pair(a, b))) continue;
                        String folder1 = output_dir + "/" + a + "_" + b + "_" + edges[a-1] + "_" + edges[b-1] + "/";
						out = new PrintWriter(folder1 + "paths.prd");
					}
				}
				if (filtering) {
                        buf.append(x + " " + y + "\n");
                }
				e1 = a;
				e2 = b;
			}
		    if (out != null) {
                if (filtering) 
                    out.println(buf);
                out.close();
            }
            inpaths.close();
             
			
			                                              
		} catch(Exception e) {
			e.printStackTrace();
		}
	}

}
