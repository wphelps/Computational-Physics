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

public class testChart {
	public static void main(String args[]) {
		int nPoints =10;
		double[] xValues =  new double[nPoints];
		double[] yValues =  new double[nPoints];
		double[] xError =  new double[nPoints];
		Random rand = new Random();

		for(int i=0; i<nPoints; i++){
			xValues[i] = Math.abs((rand.nextDouble()%4.0));
			yValues[i] = (double)i;
			xError[i] = Math.abs((rand.nextDouble()%.1));
		}
		createPlot(xValues, yValues,xError,nPoints);
		
	}
	
	
	private static void createPlot(double[] vol, double[] dim, double[] err, int num) {
		JFrame frame = new JFrame();
		JPanel chartPanel = new ChartPanel(createChart(createDataset(vol, dim, err), num));
		frame.add(chartPanel);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
	}

	private static IntervalXYDataset createDataset(double[] Volume, double[] Dimensions, double[] error) {
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

	private static JFreeChart createChart(IntervalXYDataset dataset, int n) {

		NumberAxis xAxis = new NumberAxis("Rand x");
		NumberAxis yAxis = new NumberAxis("Rand y");
		XYErrorRenderer renderer = new XYErrorRenderer();
		renderer.setDrawXError(false);
		renderer.setBaseLinesVisible(false);
		renderer.setBaseShapesVisible(true);
		XYPlot plot = new XYPlot(dataset, xAxis, yAxis, renderer);

		plot.setBackgroundPaint(Color.lightGray);
		plot.setDomainGridlinePaint(Color.white);
		plot.setRangeGridlinePaint(Color.white);

		JFreeChart chart = new JFreeChart("Random plot (N=" + n + ")", plot);
		chart.setBackgroundPaint(Color.white);
		return chart;
	}
}
