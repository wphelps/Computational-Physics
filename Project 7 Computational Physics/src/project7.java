import java.awt.Color;

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

public class project7 {
	static double sInDay = 8.64E+04;
	static double m;
	static double m1 = 2.0E+30;
	static double m2 = 1.0E+4;
	static double au = 1.4959787061E+11; // in meters
	static vector r = new vector(0.3170787, -1.887907, 0.2396305);
	static vector v = new vector(0.008180196, 0.01188029, 0.0006081217);
	static double G = 6.67E-11;
	static double l = Math.sqrt(r.cross(v).magnitudeSquared()*Math.pow(au, 4)/Math.pow(sInDay, 2));
	static double radius = r.magnitude()*au;
	static double phi = 3.0*Math.PI/2.0+ r.phi();
	static double theta = r.phi();
	static double pr = v.dot(r)/radius * au/sInDay;
	static double dt = sInDay;
	static double x = radius/au * Math.cos(phi);
	static double y = radius/au * Math.sin(phi);
	
	public static void main(String[] args) {
		m = (m1 * m2) / (m1 + m2);
		XYIntervalSeries position = new XYIntervalSeries("Comet Position");
		
		XYIntervalSeries trajectory = new XYIntervalSeries("Trajectory");

		for (int i = 0; i < 365*10; i++) {
			//System.out.println(i+" " + radius + " " +pr+ " "+phi);
			radius += pr * dt;
			pr     += Uprime(radius) * dt;
			
			phi += l / Math.pow(radius, 2) * dt;
			phi = phi % (2.0 * Math.PI);
			x = radius/au * Math.cos(phi);
			y = radius/au * Math.sin(phi);
			
			position.add(i, i, i,radius/au, radius/au, radius/au);
			trajectory.add(x,x,x,y,y,y);
		}
		
		
		XYIntervalSeriesCollection dataset = new XYIntervalSeriesCollection();
		dataset.addSeries(position);
		createWindow(createXYChart(dataset,  "Time (Days)","Position (AU)", "Comet Radius"), "Project #7");
		
		XYIntervalSeriesCollection dataset2 = new XYIntervalSeriesCollection();
		dataset2.addSeries(trajectory);
		createWindow(createXYChart(dataset2,  "X (AU)","Y (AU)", "Comet Position"), "Project #7");

	}
	public static double Uprime(double r) {
		return Math.pow(l, 2) / ( Math.pow(r, 3)) - G * m1 / Math.pow(r, 2);
	}

	public static double U(double r) {
		return Math.pow(l, 2) / ( 2.0*m*Math.pow(r, 2)) - G * m1*m2 / Math.pow(r, 1);
	}

	private static void createWindow(JFreeChart chart, String frameTitle) {
		JFrame frame = new JFrame(frameTitle);
		JPanel chartPanel = new ChartPanel(chart);
		frame.add(chartPanel);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
	}

	private static JFreeChart createXYChart(IntervalXYDataset dataset, String xLabel, String yLabel, String title) {
		NumberAxis xAxis = new NumberAxis(xLabel);
		NumberAxis yAxis = new NumberAxis(yLabel);
		XYErrorRenderer renderer = new XYErrorRenderer();
		renderer.setDrawXError(false);
		renderer.setBaseLinesVisible(false);
		renderer.setBaseShapesVisible(true);
		XYPlot plot = new XYPlot(dataset, xAxis, yAxis, renderer);

		plot.setBackgroundPaint(Color.lightGray);
		plot.setDomainGridlinePaint(Color.white);
		plot.setRangeGridlinePaint(Color.white);

		JFreeChart chart = new JFreeChart(title, plot);
		chart.setBackgroundPaint(Color.white);
		return chart;
	}

}

// Vector object
class vector {
	double x, y, z;
	vector() {
		x = 0.0;
		y = 0.0;
		z = 0.0;
	}
	vector(double x, double y, double z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}
	public double magnitudeSquared() {
		return Math.pow(x, 2) + Math.pow(y, 2) + Math.pow(z, 2);
	}

	public double magnitude() {
		return Math.sqrt(magnitudeSquared());
	}

	public double difference(vector v) {
		vector diff = new vector(x - v.x, y - v.y, z - v.z);
		return diff.magnitude();
	}

	public double[] toDoubleArray() {
		double[] ret = new double[3];
		ret[0] = x;
		ret[1] = y;
		ret[2] = z;
		return ret;
	}

	public vector cross(vector A) {
		return new vector(y * A.z - z * A.y, z * A.x - x * A.z, x * A.y - y * A.x);
	}

	public double r() {
		return magnitude();
	}

	public double cosTheta() {
		double ret = z / r();
		return (ret);
	}

	public double theta() {
		double ret = Math.acos(this.cosTheta());
		return (ret);
	}

	public double phi() {
		double ret = Math.atan2(y, x);
		return (ret);
	}
	public double dot(vector v) {
		return x * v.x + y * v.y + z * v.z;
	}
}