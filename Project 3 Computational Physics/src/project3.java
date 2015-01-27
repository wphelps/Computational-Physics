import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

import javax.swing.ButtonGroup;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYErrorRenderer;
import org.jfree.data.xy.IntervalXYDataset;
import org.jfree.data.xy.XYIntervalSeries;
import org.jfree.data.xy.XYIntervalSeriesCollection;

public class project3 {
	static ArrayList<element> elements;
	static ArrayList<JRadioButton> buttons;
	static JRadioButton int_button;
	static JRadioButton c_button;
	static JRadioButton rms_button;
	static JPanel chartPanel = new JPanel();
	static int mode;
	static int elementNumber;
	static JPanel mainPanel;
	static JFrame frame;
	public static void main(String[] args) throws FileNotFoundException {
		//Settings
		int n = 16;
		String dataFile = "gauss_laguerre.dat";
		
		//Initialize ArrayList of elements
		 elements = new ArrayList<element>();
		elements.add(new element("Ca",  40, 20));
		elements.add(new element("Fe",  56, 26));
		elements.add(new element("Cu",  63, 29));
		elements.add(new element("Cu",  65, 29));
		elements.add(new element("Au", 197, 79));
		elements.add(new element("Pb", 208, 82));
		elements.add(new element("U",  238, 92));
		
		//Read in file to 3 dimensional array
		Scanner scanner = new Scanner(new File(dataFile));
		double[][][] gauss_laguerre =  new double[n][][];
		for (int j = 0; scanner.hasNext(); j++) {
			if (scanner.hasNextInt()) {
				int nLines = scanner.nextInt();
				gauss_laguerre[j] = new double[nLines][2];
				for (int i = 0; i < nLines; i++) {
					gauss_laguerre[j][i][0] = scanner.nextDouble();
					gauss_laguerre[j][i][1] = scanner.nextDouble(); 
					//double x = gauss_laguerre[j][i][0];
					//double w = gauss_laguerre[j][i][1];
					//System.out.println("N:["+(j+1)+ "] w:["+gauss_laguerre[j][i][0] + "] x:[" + gauss_laguerre[j][i][1]+"]");
				}
			}
		}
		
		//Calculate C for each element
		double[][] C = new double[elements.size()][n];
		double[][] Integral = new double[elements.size()][n];
		for(element e: elements){
			for(int i = 0; i<gauss_laguerre.length; i++){
				double I = 0.0;
				for(int j = 0; j<gauss_laguerre[i].length; j++){
					double x = gauss_laguerre[i][j][0];
					double w = gauss_laguerre[i][j][1];
					//System.out.println("N:"+i+ " x "+x + " w " + w);
					I += w * Math.pow(x, 2)/(Math.pow(Math.E, -x) + Math.pow(Math.E, -2.0 * Math.pow(e.A, (1.0 / 3.0))));
				}
				//System.out.println(e.name + " I:"+I);
				Integral[elements.indexOf(e)][i]  = I;
				C[elements.indexOf(e)][i] = e.Q / (4.0 * Math.PI * Math.pow(.55, 3) * I);
			}
		}
		
		//Calculate charge radii
		double[][] rms = new double[elements.size()][n];
		for(element e: elements){
			for(int i = 0; i<gauss_laguerre.length; i++){
				double I = 0.0;
				for(int j = 0; j<gauss_laguerre[i].length; j++){
					double x = gauss_laguerre[i][j][0];
					double w = gauss_laguerre[i][j][1];
					I += w* Math.pow(x, 4)/(Math.pow(Math.E, -x) + Math.pow(Math.E, -2.0 * Math.pow(e.A, (1.0 / 3.0))));
				}
				//rms[elements.indexOf(e)][i] = Math.sqrt(4.0*Math.PI*Math.pow(.55,5)*C[elements.indexOf(e)][i]*I);
				rms[elements.indexOf(e)][i] = Math.sqrt((Math.pow(.55,2)*I)/Integral[elements.indexOf(e)][i]);
			}
		}
		/*
		//Print out results
		for(element e: elements){
			for(int i = 0; i<gauss_laguerre.length; i++){
				System.out.println("Element:["+e.name+"] C:["+C[elements.indexOf(e)][i]+"] sqrt(<r^2>):["+rms[elements.indexOf(e)][i]+"]");
			}
		}
		*/
		
		initiateGUI(Integral,C,rms);
	}
	
	private static void initiateGUI(final double[][] Integral, final double[][] C, final double[][] rms){
		mode = 0;
		elementNumber = 0;
		
		frame = new JFrame();
		mainPanel = new JPanel(new BorderLayout());
		JPanel elementPanel =  new JPanel(new GridLayout(elements.size(), 1, 200, 20));
		buttons = new ArrayList<JRadioButton>();
		ButtonGroup elementGroup = new ButtonGroup();
		for(element e: elements){
			buttons.add(new JRadioButton(e.name));
		}
		buttons.get(0).setSelected(true);
		for(JRadioButton b : buttons){
			elementGroup.add(b);
			elementPanel.add(b);
			b.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent event) {
					update(Integral,C, rms);
				}
			});
		}
		
		JPanel controlPanel = new JPanel();
		ButtonGroup controlGroup = new ButtonGroup();
		int_button = new JRadioButton("Integral");
		int_button.setSelected(true);
		c_button = new JRadioButton("C");
		rms_button = new JRadioButton("RMS");
		controlGroup.add(int_button);
		controlGroup.add(c_button);
		controlGroup.add(rms_button);
		controlPanel.add(int_button);
		controlPanel.add(c_button);
		controlPanel.add(rms_button);
		update(Integral,C, rms);
		int_button.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent event) {
				update(Integral,C, rms);
			}
		});
		c_button.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent event) {
				update(Integral,C, rms);
			}
		});
		rms_button.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent event) {
				update(Integral,C, rms);
			}
		});
		
		mainPanel.add(controlPanel, BorderLayout.NORTH);
		mainPanel.add(elementPanel, BorderLayout.WEST);
		mainPanel.add(chartPanel,BorderLayout.EAST);
		frame.add(mainPanel);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
	}
	
	

	protected static void update(final double[][] Integral, final double[][] C, final double[][] rms) {
		elementNumber = 0;
		for(JRadioButton b : buttons){
			if(b.isSelected()){
				elementNumber = buttons.indexOf(b);
			}
		}
		mode = 0;
		double[] n1 = new double[16];
		double[] err = new double[16];
		for(int i = 1; i<=16; i++){
			n1[i-1] = (double) i;
			err[i-1] = 0.0;
		}
		
		mainPanel.remove(chartPanel);
		if(int_button.isSelected()){
			chartPanel = new ChartPanel((createChart(createDataset(Integral[elementNumber],n1,err),elementNumber,"Integral","Integral")));
			mode = 0;
		}else if(c_button.isSelected()){
			chartPanel = new ChartPanel((createChart(createDataset(C[elementNumber],n1,err),elementNumber,"C","Integration Constant")));
			mode = 1;
		}else if(rms_button.isSelected()){
			chartPanel = new ChartPanel((createChart(createDataset(rms[elementNumber],n1,err),elementNumber,"RMS (fm)","RMS")));
			mode = 2;
		}
		mainPanel.add(chartPanel,BorderLayout.EAST);
		frame.pack();
		//frame.setVisible(true);
		//System.out.println("elementNumber:["+elementNumber+"] mode:["+mode+"]");
	}
	
	
	
	private static IntervalXYDataset createDataset(double[] Volume,double[] Dimensions, double[] error) {
		XYIntervalSeriesCollection dataset = new XYIntervalSeriesCollection();
		ArrayList<XYIntervalSeries> series = new ArrayList<XYIntervalSeries>();
		series.add(new XYIntervalSeries("Set 1"));
		for (int i = 0; i < Volume.length; i++) {
			series.get(0).add(Dimensions[i], Dimensions[i], Dimensions[i], Volume[i],
					Volume[i] - error[i], (Volume[i] + error[i]));
		}
		dataset.addSeries(series.get(0));
		return dataset;
	}

	private static JFreeChart createChart(IntervalXYDataset dataset, int n,String yLabel, String title) {
		
		NumberAxis xAxis = new NumberAxis("N");
		NumberAxis yAxis = new NumberAxis(yLabel);
		XYErrorRenderer renderer = new XYErrorRenderer();
		renderer.setDrawXError(false);
		renderer.setBaseLinesVisible(false);
		renderer.setBaseShapesVisible(true);
		XYPlot plot = new XYPlot(dataset, xAxis, yAxis, renderer);

		plot.setBackgroundPaint(Color.lightGray);
		plot.setDomainGridlinePaint(Color.white);
		plot.setRangeGridlinePaint(Color.white);

		JFreeChart chart = new JFreeChart(title+" Plot for element ["+elements.get(n).name+"]", plot);
		chart.setBackgroundPaint(Color.white);
		return chart;
	}
}

class element {
	String name;
	int A;
	int Q;
	element(String name, int A, int Q) {
		this.name = name;
		this.A = A;
		this.Q = Q;
	}
	
}