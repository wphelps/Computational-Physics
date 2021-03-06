import java.awt.Color;
import java.util.ArrayList;
import java.util.Random;

import javax.swing.JFrame;
import javax.swing.JPanel;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYErrorRenderer;
import org.jfree.data.xy.IntervalXYDataset;
import org.jfree.data.xy.XYIntervalSeries;
import org.jfree.data.xy.XYIntervalSeriesCollection;

public class project4MultiThreading {
	int num = 1000000;
	int Dimensions = 25;
	double err[] = new double[Dimensions];
	double vol[] = new double[Dimensions];
	double dim[] = new double[Dimensions];
	
	public static void main(String args[]) throws InterruptedException {
		new project4MultiThreading();
	}
	ArrayList<ThreadWithResult> threads = new ArrayList<ThreadWithResult>();
	project4MultiThreading() throws InterruptedException{
		for (int D = 1; D < Dimensions; D++) {
			ResultSetter setter = new ResultSetter() {
				public void setResult(int dimension, double volume, double error) {
					err[dimension-1] = error;
					vol[dimension-1] = volume;
					dim[dimension-1] = (double) dimension;
					//System.out.println("Setter:"+volume+" "+ error + " " + dimension);
				}
			};
			ThreadWithResult thread = new ThreadWithResult();
			threads.add(thread);
			thread.setResultSetter(setter);
			thread.setDimension(D);
			thread.setNum(num);
			thread.start();
		}
         for(ThreadWithResult thread: threads){
             thread.join();  
         }
		createPlot(vol, dim, err, num);
	}

	private void createPlot(double[] vol, double[] dim, double[] err, int num) {
		JFrame frame = new JFrame();
		JPanel chartPanel = new ChartPanel(createChart(createDataset(vol, dim, err), num));
		frame.add(chartPanel);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
	}

	private IntervalXYDataset createDataset(double[] Volume, double[] Dimensions, double[] error) {
		XYIntervalSeriesCollection dataset = new XYIntervalSeriesCollection();
		ArrayList<XYIntervalSeries> series = new ArrayList<XYIntervalSeries>();
		series.add(new XYIntervalSeries("Set 1"));
		for (int i = 0; i < Volume.length; i++) {
			series.get(0).add(Dimensions[i], Dimensions[i], Dimensions[i], Volume[i], Volume[i] - error[i],
					(Volume[i] + error[i]));
		}
		dataset.addSeries(series.get(0));
		return dataset;
	}

	private JFreeChart createChart(IntervalXYDataset dataset, int n) {

		NumberAxis xAxis = new NumberAxis("Dimension");
		NumberAxis yAxis = new NumberAxis("Volume");
		XYErrorRenderer renderer = new XYErrorRenderer();
		renderer.setDrawXError(false);
		renderer.setBaseLinesVisible(false);
		renderer.setBaseShapesVisible(true);
		XYPlot plot = new XYPlot(dataset, xAxis, yAxis, renderer);

		plot.setBackgroundPaint(Color.lightGray);
		plot.setDomainGridlinePaint(Color.white);
		plot.setRangeGridlinePaint(Color.white);

		JFreeChart chart = new JFreeChart("Monte Carlo Spherical Volume Estimator (N=" + n + ")", plot);
		chart.setBackgroundPaint(Color.white);
		return chart;
	}
	
	interface ResultSetter {
		public void setResult(int dimension, double volume, double error);
	}
	
	class ThreadWithResult extends Thread {
		private ResultSetter setter;
		private int dimension;
		private int num;
		public void setResultSetter(ResultSetter setter) {
			this.setter = setter;
		}
		
		public void setDimension(int dimension) {
			this.dimension = dimension;
		}
		
		public void setNum(int num) {
			this.num = num;
		}
		
		public void run() {
			double error = 0;
			double volume = 0;
			int N, N0;
			double err1,err2;
			double x[] = new double[dimension];
			N = N0 = 0;
			err1 = err2 = error = 0.0;
			Random rand = new Random();
			for (long j = 0; j < num; j++) {
				double sum = 0.0;
				for (int k = 0; k < x.length; k++) {
					x[k] = (rand.nextDouble() % 2.0) - 1.0;
					sum += x[k] * x[k];
				}
				sum = Math.sqrt(sum);
				N0 = (Math.pow(sum, 2) < 1) ? N0 + 1 : N0;
				N++;

				err1 = (err1 * ((double) N - 1.0) + Math.pow((Math.pow(2.0, dimension) * (double) N0 / (double) N), 2)) / N;
				err2 = (err2 * ((double) N - 1.0) + (Math.pow(2.0, dimension) * (double) N0 / (double) N)) / N;
			}
			error = err1 - Math.pow(err2, 2);
			volume = Math.pow(2.0, dimension) * (double) N0 / (double) N;
			System.out.println("V" + dimension + ":[" + Math.pow(2.0, dimension) * (double) N0 / (double) N + "] Sigma^2:[" + error
					+ "] Sigma:[" + Math.pow(error, .5) + "]");

			setter.setResult(dimension, volume, error);
		}
	}
}

