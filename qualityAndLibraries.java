import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class QHC2 {
	private static int cp = -1;
	private static int cpdiff = -1;
	private static int n = 30;
	private static double th = 0.4;
	final private static String afeesfn = "X:\\Data\\Collaborations\\Afees\\SSw\\ValidateGC\\AfeesGCTest.txt";
	private static int var = -1;
	private static int varindex = -1;
	private static int todo = -1;
	private static String store_ss = "";

	public static void main(String args[]) {
		try {
			
				TestRig2();
			
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public static void TestRig2() throws FileNotFoundException {
		PrintStream prt = new PrintStream(new File("./sampleResults/04.txt"));
		PrintStream console = System.out;
		System.setOut(prt);
		System.out.println("Name of File" + " ");
		System.out.println(" k " + "," + " repid " + "," + "fitbasic " + "," + "timebasic " + "," + "cpb " + "," + "fd " + "," + "td " + "," + "cpd " + ","
				+ "Et " + "," + "bestExf " + "," + "ssex " + ","+"ssb "+","+"ssd ");
		for ( int count =2 ; count <=32; count++) {
			String fileId = (count < 10 ? "0" : "") + count;
			String fileName = String.format("./DataWK2/04_%s.csv", fileId);
			double wk[][] = readCSVFile(fileName);
			th = ComputeThreshHold(wk);
			long before = System.nanoTime();
			double bestExf = Exhaustive(wk, 30);
			//double bestExf = 0.0;
			long Et = System.nanoTime() - before;
			String ssex = store_ss;
			for (int i = 0; i < 100; ++i) {
				before = System.nanoTime();
				double fb = DoBasicHC(wk);
				int cpb = cp;
				long tb = System.nanoTime() - before;
				String ssb = store_ss;
				before = System.nanoTime();
				double fd = DoBasicHCDiff(wk);
				int cpd = cpdiff;
				long td = System.nanoTime() - before;
				String ssd = store_ss;
				System.setOut(prt);
				System.out.println(count + "," + i + "," + fb + "," + tb + "," + cpb + "," + fd + "," + td + "," + cpd + ","
						+ Et + "," + bestExf+ "," + ssex + ","+ssb+","+ssd);
		}
	}
		prt.close();
	}

	private static double Exhaustive(double[][] wk, int n) {
		double bestf = Double.MIN_VALUE;
		String bestss = "";
		for (int i = 0; i < (int) Math.pow(2, n); ++i) {
			String bin = Integer.toBinaryString(i);
			if (bin.length() < n) {
				bin = "0".repeat(n - bin.length()) + bin;
			}
			// System.out.println(i + " " + bin);
			ArrayList<Integer> ss = new ArrayList<>();
			for (int j = 0; j < bin.length(); ++j)
				if (bin.charAt(j) == '1')
					ss.add(j);
			if (ss.size() > 1) {
				double f = Fitness(ss, wk);
				// System.out.println(i + " " + bin + " " + f);
				if (f > bestf) {
					
					bestf = f;
					bestss = bin;
				}
			}
			store_ss = bestss;
		}
		// System.out.println("Best: " + bestss + " " + bestf);
		return (bestf);
	}

	private static void DoFitnessSim(double[][] wk) {
		int rep = 10000;

		for (int i = 0; i < rep; ++i) {
			ArrayList<Integer> s = RandomStart();
			double f = Fitness(s, wk);
			System.out.println(i + " " + f);
		}
	}

	private static double ComputeThreshHold(double[][] wk) {
		double av = 0.0;
		double c = 0.0;
		int k = wk.length;
		for (int i = 0; i < k - 1; ++i) {
			for (int j = i + 1; j < k; ++j) {
				av += wk[i][j];
				c += 1.0;
			}
		}
		return (av / c);
	}

	@SuppressWarnings("unchecked")
	private static double DoBasicHC(double[][] wk) {
		cp = 0;
		ArrayList<Integer> s = RandomStart();
		double f = Fitness(s, wk);
		int iter = 10000000;

		for (int i = 0; i < iter; ++i) {
			ArrayList<Integer> sdash = SmallChange(s);
			double fdash = Fitness(sdash, wk);
			if (fdash >= f) {
				s = (ArrayList<Integer>) sdash.clone();
				f = fdash;
				cp = i;
			}
		}
		store_ss = ALtoString(s,wk.length);
		return (f);
	}

	private static double DoBasicHCDiff(double[][] wk) {
		cpdiff = 0;
		ArrayList<Integer> s = RandomStart();
		double f = Fitness(s, wk);
		int iter = 10000000;

		for (int i = 0; i < iter; ++i) {
			double fdash = SmallChangeDiff(s, f, wk);
			if (fdash >= f) {
				DoSmallChange(s);
				f = fdash;
				cpdiff = i;
			}
		}
		store_ss = ALtoString(s,wk.length);
		return (f);
	}

	private static String ALtoString(ArrayList<Integer> s,int n) 
	{
		String ss = "";
		for(int i=0;i<n;++i)
		{
			if (s.contains(i) == true)
			{
				ss = ss + "1";
			}
			else
			{
				ss = ss + "0";
			}
		}
		return(ss);
	}

	private static void DoSmallChange(ArrayList<Integer> s) {
		if (todo == 1) {
			s.remove(varindex);
		} else {
			s.add(var);
		}
	}

	private static double SmallChangeDiff(ArrayList<Integer> s, double fold, double wk[][]) {
		double fnew = fold;
		int k = s.size();
		todo = UI(1, 2);
		if (s.size() == n)
			todo = 1;
		if (s.size() == 2)
			todo = 2;
		if (todo == 1) {
			varindex = UI(0, k - 1);
			var = s.get(varindex);
			for (int i = 0; i < k; ++i) {
				if (i != varindex) {
					int x = s.get(i);
					fnew -= (wk[x][var] - th);
				}
			}
		} else {
				do {
				var = UI(0, n - 1);
			} while (s.contains(var));
			varindex = -1;
			for (int i = 0; i < k; ++i) {
				int x = s.get(i);
				fnew += (wk[x][var] - th);
			}
		}
		return (fnew);
	}

	private static ArrayList<Integer> SmallChange(ArrayList<Integer> s) {
		@SuppressWarnings("unchecked")
		ArrayList<Integer> sdash = (ArrayList<Integer>) (s.clone());
		int k = sdash.size();
		todo = UI(1, 2);
		if (sdash.size() == n)
			todo = 1;
		if (sdash.size() == 2)
			todo = 2;
		if (todo == 1) {
			// Remove an item
			int i = UI(0, k - 1);
			sdash.remove(i);
		} else {
			int x = -1;
			// Add an item
			do {
				x = UI(0, n - 1);
			} while (sdash.contains(x));
			sdash.add(x);
		}
		return (sdash);
	}

	private static double Fitness(ArrayList<Integer> s, double wk[][]) {
		double f = 0.0;
		int k = s.size();
		for (int i = 0; i < k - 1; ++i) {
			for (int j = i + 1; j < k; ++j) {
				int si = s.get(i);
				int sj = s.get(j);
				f += (wk[si][sj] - th);
			}
		}
		return (f);
	}

	private static double[][] AfeesWK() {
		double wk[][] = {
				{ 0, 0.495297453, 0.500513802, 0.329461302, 0.610565284, 0.158452849, 0.228240307, 0.409039334,
						0.810661144, 0.512720926 },
				{ 0.495297453, 0, 0.185909895, 0.24446673, 0.726895921, 0.203204702, 0.596463194, 0.23341616,
						0.74563442, 0.291514941 },
				{ 0.500513802, 0.185909895, 0, 0.520274984, 0.701755535, 0.876365997, 0.759703806, 0.50563971,
						0.669195343, 0.127061511 },
				{ 0.329461302, 0.24446673, 0.520274984, 0, 0.253250901, 0.585403556, 0.592234721, 0.339351459,
						0.483439122, 0.299736668 },
				{ 0.610565284, 0.726895921, 0.701755535, 0.253250901, 0, 0.548054823, 0.774874784, 0.446209111,
						0.484820151, 0.630068214 },
				{ 0.158452849, 0.203204702, 0.876365997, 0.585403556, 0.548054823, 0, 0.689319186, 0.370919698,
						0.784223826, 0.422548543 },
				{ 0.228240307, 0.596463194, 0.759703806, 0.592234721, 0.774874784, 0.689319186, 0, 0.610019975,
						0.49962605, 0.3672312 },
				{ 0.409039334, 0.23341616, 0.50563971, 0.339351459, 0.446209111, 0.370919698, 0.610019975, 0,
						0.956316672, 0.237204032 },
				{ 0.810661144, 0.74563442, 0.669195343, 0.483439122, 0.484820151, 0.784223826, 0.49962605, 0.956316672,
						0, 0.489938696 },
				{ 0.512720926, 0.291514941, 0.127061511, 0.299736668, 0.630068214, 0.422548543, 0.3672312, 0.237204032,
						0.489938696, 0 } };

		return (wk);

	}

	private static ArrayList<Integer> RandomStart() {
		// ArrayList<Integer> m = new ArrayList<Integer>();
		ArrayList<Integer> s = new ArrayList<>();
		for (int i = 0; i < n; ++i) {
			if (UI(0, 1) == 1)
				s.add(i);
		}
		while (s.size() < 2) {
			int x = UI(0, n - 1);
			if (!s.contains(x))
				s.add(x);
		}
		return (s);
	}

	private static int UI(int a, int b) {
		double d = Math.abs(a - b) + 1.0;
		int min = Math.min(a, b);
		int max = Math.max(a, b);
		double x = Math.floor(Math.random() * d) + min;
		if (x < min)
			x = min;
		if (x > max)
			x = max;
		return ((int) x);
	}

	private static void PrintArray(double[][] wk) {
		for (double[] element : wk) {
			for (int j = 0; j < element.length; ++j) {
				System.out.print(element[j] + " ");
			}
			System.out.println();
		}
	}

	private static double[][] CreateWK(int n) {
		double wk[][] = new double[n][n];
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				wk[i][j] = Math.abs(Math.random());
				wk[j][i] = wk[i][j];
			}
			wk[i][i] = 1.0;
		}
		return (wk);
	}


	// Find the index of the second smallest value in an array
	private static int SecondSmallest(double x[]) {
		int i1 = -1, i2 = -1;
		double m1 = 0.0, m2 = 0.0;

		if (x[0] < x[1]) {
			m1 = x[0];
			m2 = x[1];
			i1 = 0;
			i2 = 1;
		} else {
			m1 = x[1];
			m2 = x[0];
			i1 = 1;
			i2 = 0;
		}
		for (int i = 2; i < x.length; ++i) {
			double y = x[i];

			// y < m1 < m2
			if (y < m1) {
				m2 = m1;
				m1 = y;
				i2 = i1;
				i1 = i;
			}
			// m1 < y < m2
			if (m1 < y && y < m2) {
				m2 = y;
				i2 = i;
			}
		}
		return (i2);
	}

	@SuppressWarnings("deprecation")
	public static double[][] readCSVFile(String fileName) {
		try {
			String thisLine;
			FileInputStream fis = new FileInputStream(fileName);
			DataInputStream myInput = new DataInputStream(fis);

			List<double[]> lines = new ArrayList<>();
			while ((thisLine = myInput.readLine()) != null) {

				String[] tuples = thisLine.split(",");
				double[] dataset = new double[tuples.length];

				for (int i = 0; i < tuples.length; i++) {
					dataset[i] = Double.parseDouble(tuples[i]);
				}
				lines.add(dataset);
			}

			// convert our list to a String array.
			double[][] array = new double[lines.size()][0];
			lines.toArray(array);
			return array;
		} catch (Exception exp) {
			exp.printStackTrace();
		}
		return null;
	}

	static public double[][] ReadArrayFile(String filename, String sep) {
		double res[][] = null;
		try {
			BufferedReader input = null;
			input = new BufferedReader(new FileReader(filename));
			String line = null;
			int ncol = 0;
			int nrow = 0;

			while ((line = input.readLine()) != null) {
				++nrow;
				String[] columns = line.split(sep);
				ncol = Math.max(ncol, columns.length);
			}
			res = new double[nrow][ncol];
			input.close();
			input = new BufferedReader(new FileReader(filename));
			int i = 0, j = 0;
			while ((line = input.readLine()) != null) {

				String[] columns = line.split(sep);
				for (j = 0; j < columns.length; ++j) {
					res[i][j] = Double.parseDouble(columns[j]);
				}
				++i;
			}
			input.close();
		} catch (Exception E) {
			// System.out.println("SSS");
			System.out.println(E.getMessage());
		}
		return (res);
	}
}