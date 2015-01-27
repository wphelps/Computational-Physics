import java.awt.Color;
import java.util.ArrayList;

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

public class project6 {
	static double h,l,g,m,q0,p0,a,w,t;
	static ArrayList<XYIntervalSeries> series;
	public static void main(String[] args) {
		l = g = m = 1.0; 
		
		//Initial conditions
		int n = 20000;
		h  = 0.001;
		w = .5;
		a = 0.5;
		t  = 0.0;
		q0 = 0.2;
		p0 = 0.0;
		
		XYIntervalSeriesCollection dataset = new XYIntervalSeriesCollection();
		series = new ArrayList<XYIntervalSeries>();
		series.add(new XYIntervalSeries("Phase Space"));
		series.add(new XYIntervalSeries("Position"));

		series = recursion(q0, p0, n);
		dataset.addSeries(series.get(0));
		createWindow(createXYChart(dataset,"Momentum","Position","Pendulum Phase Space (N="+n+" h="+h+" a="+a+" w="+w+")"),"Project #6");	
		
		XYIntervalSeriesCollection dataset2 = new XYIntervalSeriesCollection();
		dataset2.addSeries(series.get(1));
		createWindow(createXYChart(dataset2,"Time","Position","Pendulum Position (N="+n+" h="+h+" a="+a+" w="+w+")"),"Project #6");	
	}
	
	public static ArrayList<XYIntervalSeries> recursion(double p, double q, int i){
		t += h;
		if(t/h < i){
			double pph = 2.0*h*F(q,t) + p;
			double qph = 2.0*h*p/m + q;
			series.get(0).add(p,p,p,q,q,q);
			series.get(1).add(t,t,t,q,q,q);
			recursion(pph, qph, i);
		}
		return series;
	}
	
	public static double F(double q, double t){
		return a*Math.cos(w*t) - q;
	}
	
	private static void createWindow(JFreeChart chart,String frameTitle){
		JFrame frame = new JFrame(frameTitle);
		JPanel chartPanel = new ChartPanel(chart);
		frame.add(chartPanel);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
	}

	private static JFreeChart createXYChart(IntervalXYDataset dataset,String xLabel, String yLabel, String title) {
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
