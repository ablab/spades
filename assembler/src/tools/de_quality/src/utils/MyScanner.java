package utils;
import java.io.*;
import java.util.*;

public class MyScanner{
	private BufferedReader in;
	private StringTokenizer st = null;

	public String nextToken() {
		if (hasMoreTokens()){
			return st.nextToken();
		}
		return null;
	}

	public MyScanner(String file){
	        try{
		in = new BufferedReader(new FileReader(new File(file)));
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	public MyScanner(InputStream inp){
		try{                                   
                	in = new BufferedReader(new InputStreamReader(inp));
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public boolean hasMoreTokens() {
		String s = null;
		try{
			while ((st==null || !st.hasMoreTokens()) && ((s=in.readLine()) != null)) 
                st = new StringTokenizer(s);
			if ((st==null || !st.hasMoreTokens()) && s==null) 
                return false;
	    } catch(IOException e) {
	        e.printStackTrace();
	    }
		return true;
	}

	public int nextInt() {
		return Integer.parseInt(nextToken());
	}

	public long nextLong() {
		return Long.parseLong(nextToken());
	}

	public double nextDouble() {
		return Double.parseDouble(nextToken());
	}


	public String nextString() {
		return nextToken();
	}

    public void close() {
		try{
			in.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
}

