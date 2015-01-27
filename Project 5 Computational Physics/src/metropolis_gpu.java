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

public class metropolis_gpu {
	final static int M = 100; // N steps
	final static int L = 1;
	final static int N = 9*50;

	public static void main(String args[]) {
		ensemble_gpu ensemble = new ensemble_gpu(N);
		double[] kineticEnergy = new double[M];
		double[] potentialEnergy = new double[M];
		double[] kineticEnergySquared = new double[M];
		double[] potentialEnergySquared = new double[M];
		double[] interactionPotential = new double[M];
		double[] total = new double[M];
		/*
		 * double totalEnergy = 0.0; double totalEnergySquared = 0.0; double[]
		 * specificHeat = new double[M];
		 */

		for (int i = 0; i < M; i++) {
			ensemble.oneCycle();
			if (i % L == 0) {
				ensemble.k++;
				System.out.printf("%6d%7.4f%4d%14.6e complete:[%5.2f]\n", ensemble.accepted, ensemble.acceptance(),
						ensemble.k, (ensemble.kineticEnergy() + ensemble.potentialEnergy()) / (double) ensemble.k,
						((double) i / (double) M) * 100.0);
			}
			kineticEnergy[i] = ensemble.kineticEnergy();
			kineticEnergySquared[i] = Math.pow(kineticEnergy[i], 2);
			potentialEnergy[i] = ensemble.potentialEnergy();
			potentialEnergySquared[i] = Math.pow(potentialEnergy[i], 2);
			interactionPotential[i] = ensemble.interactioPotential();
			if (interactionPotential[i] > 1000000) {
				interactionPotential[i] = 1000000;
			}
			total[i] = potentialEnergy[i] + kineticEnergy[i] + interactionPotential[i];
			/*
			 * if (i > 500) { totalEnergySquared = (totalEnergySquared *
			 * ((double) i - 501.0) + Math.pow(kineticEnergySquared[i] +
			 * potentialEnergySquared[i], 2)) / i - 501.0; totalEnergy =
			 * (totalEnergy * ((double) i - 501.0) + (kineticEnergy[i] +
			 * potentialEnergy[i])) / i - 501.0; specificHeat[i] =
			 * (totalEnergySquared - Math.pow(totalEnergy, 2));
			 * System.out.println("Specific Heat:[" + specificHeat[i] + "]"); }
			 */
		}

		double[] zero = new double[M];
		for (int i = 0; i < M; i++) {
			zero[i] = 0.0;
		}
		double[] steps = new double[M];
		for (int i = 0; i < M; i++) {
			steps[i] = (double) i;
		}
		createPlot(kineticEnergy, potentialEnergy, total, interactionPotential, steps, zero, N);
	}

	private static void createPlot(double[] vol, double[] potential, double[] total, double[] interaction,
			double[] dim, double[] err, int num) {
		JFrame frame = new JFrame();
		JPanel chartPanel = new ChartPanel(
				createChart(createDataset(vol, potential, total, interaction, dim, err), num));
		frame.add(chartPanel);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
	}

	private static IntervalXYDataset createDataset(double[] kinetic, double[] potential, double[] total,
			double[] interaction, double[] Dimensions, double[] error) {
		XYIntervalSeriesCollection dataset = new XYIntervalSeriesCollection();
		ArrayList<XYIntervalSeries> series = new ArrayList<XYIntervalSeries>();
		series.add(new XYIntervalSeries("Kinetic"));
		for (int i = 0; i < kinetic.length; i++) {
			series.get(0).add(Dimensions[i], Dimensions[i], Dimensions[i], kinetic[i], kinetic[i], kinetic[i]);
		}
		series.add(new XYIntervalSeries("Gravitational Potential"));
		for (int i = 0; i < potential.length; i++) {
			series.get(1).add(Dimensions[i], Dimensions[i], Dimensions[i], potential[i], potential[i], potential[i]);
		}
		series.add(new XYIntervalSeries("Interaction Potential"));
		for (int i = 0; i < potential.length; i++) {
			series.get(2).add(Dimensions[i], Dimensions[i], Dimensions[i], interaction[i], interaction[i],
					interaction[i]);
		}
		series.add(new XYIntervalSeries("Total"));
		for (int i = 0; i < potential.length; i++) {
			series.get(3).add(Dimensions[i], Dimensions[i], Dimensions[i], total[i], total[i], total[i]);
		}
		dataset.addSeries(series.get(0));
		dataset.addSeries(series.get(1));
		dataset.addSeries(series.get(2));
		dataset.addSeries(series.get(3));
		return dataset;
	}

	private static JFreeChart createChart(IntervalXYDataset dataset, int n) {
		NumberAxis xAxis = new NumberAxis("Number of iterations");
		NumberAxis yAxis = new NumberAxis("Energy");
		XYErrorRenderer renderer = new XYErrorRenderer();
		renderer.setDrawXError(false);
		renderer.setBaseLinesVisible(false);
		renderer.setBaseShapesVisible(true);
		XYPlot plot = new XYPlot(dataset, xAxis, yAxis, renderer);

		plot.setBackgroundPaint(Color.lightGray);
		plot.setDomainGridlinePaint(Color.white);
		plot.setRangeGridlinePaint(Color.white);

		JFreeChart chart = new JFreeChart("Metropolis (N particles=" + n + ")", plot);
		chart.setBackgroundPaint(Color.white);
		return chart;
	}
}

class ensemble_gpu {
	ArrayList<particle> particles = new ArrayList<particle>();
	long accepted = 0;
	long rejected = 0;
	long k = 0;
	int nParticles = 0;

	public ensemble_gpu(int N) {
		double z = .00;
		for (int i = 0; i < N; i += 9) {
			particles.add(new particle(.1,.1,z));
			particles.add(new particle(.1,.5,z));
			particles.add(new particle(.1,.9,z));
			particles.add(new particle(.5,.1,z));
			particles.add(new particle(.5,.5,z));
			particles.add(new particle(.5,.9,z));
			particles.add(new particle(.9,.1,z));
			particles.add(new particle(.9,.5,z));
			particles.add(new particle(.9,.9,z));
			z+=0.45;
		}
		nParticles = particles.size();
	}

	public double kineticEnergy() {
		double kineticEnergy = 0.0;
		for (particle p : particles) {
			kineticEnergy += p.kineticEnergy();
		}
		return kineticEnergy;
	}
	public boolean tooCloseTrial(int j) {
		boolean tooClose = false;
		double distance = 0.2;
		for (int i = 0; i < particles.size(); i++) {
			if ((i != j)) {
				if (Math.abs(particles.get(i).position.difference(particles.get(j).future_position)) < distance) {
					return true;
				}
			}

		}
		return tooClose;
	}

	public boolean tooClose(int j) {
		boolean tooClose = false;
		double distance = 0.2;
		for (int i = 0; i < particles.size(); i++) {
			if ((i != j)) {
				if (Math.abs(particles.get(i).position.difference(particles.get(j).position)) < distance) {
					return true;
				}
			}

		}
		return tooClose;
	}

	public double interactioPotential() {
		OpenCL_Potential opencl = new OpenCL_Potential();
		int n = particles.size();
		int N = (n * (n - 1)) / 2;
		float[] x1 = new float[N];
		float[] y1 = new float[N];
		float[] z1 = new float[N];

		float[] x2 = new float[N];
		float[] y2 = new float[N];
		float[] z2 = new float[N];

		int k = 0;
		for (int i = 0; i < particles.size(); i++) {
			for (int j = i + 1; j < particles.size(); j++) {
				x1[k] = (float) particles.get(i).position.x;
				y1[k] = (float) particles.get(i).position.y;
				z1[k] = (float) particles.get(i).position.z;

				x2[k] = (float) particles.get(j).position.x;
				y2[k] = (float) particles.get(j).position.y;
				z2[k] = (float) particles.get(j).position.z;
				k++;
			}
		}
		opencl.setN(N);
		opencl.x_input = x1;
		opencl.y_input = y1;
		opencl.z_input = z1;

		opencl.x2_input = x2;
		opencl.y2_input = y2;
		opencl.z2_input = z2;

		opencl.Run();

		return opencl.result;
	}
	public double interactioPotentialTrial(int p) {
		OpenCL_Potential opencl = new OpenCL_Potential();
		int n = particles.size();
		//int N = (n * (n - 1)) / 2;
		float[] x1 = new float[n-1];
		float[] y1 = new float[n-1];
		float[] z1 = new float[n-1];

		float[] x2 = new float[n-1];
		float[] y2 = new float[n-1];
		float[] z2 = new float[n-1];

		int k = 0;
		float xval = (float) particles.get(p).future_position.x;
		float yval = (float) particles.get(p).future_position.y;
		float zval = (float) particles.get(p).future_position.z;

		for (int i = 0; i < particles.size() - 1; i++) {
			if (i != p) {
				x1[k] = xval;
				y1[k] = yval;
				z1[k] = zval;
				x2[k] = (float) particles.get(i).position.x;
				y2[k] = (float) particles.get(i).position.y;
				z2[k] = (float) particles.get(i).position.z;
				k++;
			}
		}
		opencl.setN(n-1);
		opencl.x_input = x1;
		opencl.y_input = y1;
		opencl.z_input = z1;

		opencl.x2_input = x2;
		opencl.y2_input = y2;
		opencl.z2_input = z2;

		opencl.Run();

		return opencl.result;
	}

	public double kineticEnergyTrial(int j) {
		double kineticEnergy = 0.0;
		for (int i = 0; i < particles.size(); i++) {
			if (i == j) {
				kineticEnergy += particles.get(i).kineticEnergyTrial();
			} else {
				kineticEnergy += particles.get(i).kineticEnergy();
			}
		}
		return kineticEnergy;
	}

	public double potentialEnergy() {
		double potentialEnergy = 0.0;
		for (particle p : particles) {
			potentialEnergy += p.position.z;
		}
		return potentialEnergy;
	}
	public double potentialEnergyTrial(int j) {
		double potentialEnergy = 0.0;
		for (int i = 0; i < particles.size(); i++) {
			if (i == j) {
				potentialEnergy += particles.get(i).future_position.z;
			} else {
				potentialEnergy += particles.get(i).position.z;
			}
		}
		return potentialEnergy;
	}

	public void oneCycle() {
		for (particle p : particles) {
			/*
			 * while (true) { p.newConfiguration();
			 * //System.out.println(p.position.x + " " + p.position.y + " " +
			 * p.position.z); //System.out.println(p.future_position.x + " " +
			 * p.future_position.y + " " + p.future_position.z);
			 * 
			 * if (p.accepted()) { accepted++; p.storePrevious(); break; } else
			 * { rejected++; } }
			 */

			p.newConfiguration();
			if (p.accepted(this, this.particles.indexOf(p))) {
				accepted++;
				p.storePrevious();
			} else {
				rejected++;
			}

			// System.out.println("Accepted:["+accepted+"] Rejected:["+rejected+"]");
		}
	}

	public double acceptance() {
		return ((double) accepted / (double) (accepted + rejected));
	}
}
