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

public class metropolis {
	final static int M = 10000; // N steps
	final static int L = 1;
	final static int N = 9*2;

	public static void main(String args[]) {
		ensemble ensemble = new ensemble(N);
		double[] kineticEnergy = new double[M];
		double[] potentialEnergy = new double[M];
		double[] kineticEnergySquared = new double[M];
		double[] potentialEnergySquared = new double[M];
		double[] interactionPotential = new double[M];
		double[] total = new double[M];
		/*double totalEnergy = 0.0;
		double totalEnergySquared = 0.0;
		double[] specificHeat = new double[M];*/
		
		for (int i = 0; i < M; i++) {
			ensemble.oneCycle();
			if (i % L == 0) {
				ensemble.k++;
				 System.out.printf("%6d%7.4f%4d%14.6e complete:[%5.2f]\n", ensemble.accepted,ensemble.acceptance(), ensemble.k,(ensemble.kineticEnergy() + ensemble.potentialEnergy()) /(double) ensemble.k,((double)i/(double)M)*100.0);
			}
			kineticEnergy[i] = ensemble.kineticEnergy();
			kineticEnergySquared[i] = Math.pow(kineticEnergy[i], 2);
			potentialEnergy[i] = ensemble.potentialEnergy();
			potentialEnergySquared[i] = Math.pow(potentialEnergy[i], 2);
			interactionPotential[i] = ensemble.interactioPotential();
			if (interactionPotential[i] > 100000) {
				interactionPotential[i] = 100000;
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

class ensemble {
	ArrayList<particle> particles = new ArrayList<particle>();
	long accepted = 0;
	long rejected = 0;
	long k = 0;
	int nParticles = 0;
	
	public ensemble(int N) {
		double z=.00;
		for (int i = 0; i < N; i+=9) {
			particles.add(new particle(.1,.1,z));
			particles.add(new particle(.1,.5,z));
			particles.add(new particle(.1,.9,z));
			particles.add(new particle(.5,.1,z));
			particles.add(new particle(.5,.5,z));
			particles.add(new particle(.5,.9,z));
			particles.add(new particle(.9,.1,z));
			particles.add(new particle(.9,.5,z));
			particles.add(new particle(.9,.9,z));
			z+=0.4;
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
		double distance = 0.05;
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
		double potential = 0.0;
		for (int i = 0; i < particles.size(); i++) {
			double radius = 0.0;
			for (int j = i + 1; j < particles.size(); j++) {
				radius = Math.abs(particles.get(i).position.difference(particles.get(j).position));
				potential += .001 / Math.pow(radius, 12) - .001 / Math.pow(radius, 6);
			}

		}
		return potential;
	}

	public double interactioPotentialTrial(int k) {
		double potential = 0.0;
		for (int i = 0; i < particles.size(); i++) {
			double radius = 0.0;
			for (int j = i + 1; j < particles.size(); j++) {
				if (i == k) {
					radius = Math.abs(particles.get(i).future_position.difference(particles.get(j).position));
				} else if (j == k) {
					radius = Math.abs(particles.get(i).position.difference(particles.get(j).future_position));
				} else {
					radius = Math.abs(particles.get(i).position.difference(particles.get(j).position));
				}
				potential += .001 / Math.pow(radius, 12) - .001 / Math.pow(radius, 6);

			}

		}
		return potential;
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

class particle {
	vector momentum;
	vector position;
	vector future_momentum;
	vector future_position;	
	double eps = .1;
	particle() {
		Random rand = new Random();
		momentum = new vector(0, 0, 0);
		position = new vector((rand.nextDouble() % 1.0), (rand.nextDouble() % 1.0), (rand.nextDouble() % 10.0) + 1);
		future_momentum = new vector();
		future_position = new vector();
	}
	particle(double z) {
		momentum = new vector(0, 0, 0);
		position = new vector(.5,.5, z);
		future_momentum = new vector();
		future_position = new vector();
	}

	public particle(double x, double y, double z) {
		momentum = new vector(0, 0, 0);
		position = new vector(x,y, z);
		future_momentum = new vector();
		future_position = new vector();	
		}
	public double kineticEnergyTrial() {
		return future_momentum.magnitudeSquared() / 2.0;
	}

	public void storePrevious() {
		momentum.x = future_momentum.x;
		momentum.y = future_momentum.y;
		momentum.z = future_momentum.z;

		position.x = future_position.x;
		position.y = future_position.y;
		position.z = future_position.z;
	}

	public void newConfiguration() {
		Random rand = new Random();
		future_momentum.x = momentum.x + ((rand.nextDouble() % 1.0) - .5) * eps;
		future_momentum.y = momentum.y + ((rand.nextDouble() % 1.0) - .5) * eps;
		future_momentum.z = momentum.z + ((rand.nextDouble() % 1.0) - .5) * eps;
		// Position
		future_position.x = position.x + ((rand.nextDouble() % 1.0) - .5) * eps;
		future_position.y = position.y + ((rand.nextDouble() % 1.0) - .5) * eps;
		future_position.z = position.z + ((rand.nextDouble() % 1.0) - .5) * eps;
	}

	public boolean isNotInBox(vector position) {
		if (position.x < 0.0 || position.x > 1.0) {
			return true;
		} else if (position.y < 0.0 || position.y > 1.0) {
			return true;
		} else if (position.z < 0.0) {
			return true;
		}
		return false;
	}

	public double kineticEnergy() {
		return momentum.magnitudeSquared() / 2.0;
	}

	public boolean accepted(ensemble ensemble, int i) {
		Random rand = new Random();
		// System.out.println("p_future():["+p_future()+"] p():["+p()+"]");
		double p_future = p_future(ensemble, i);
		double p = p(ensemble);
		boolean re = ((p_future >= p) || (p_future > (rand.nextDouble() % 1.0) * p));
		// System.out.println("p_future:["+p_future+"] p:["
		// +p+"] return:["+re+"]");
		return re;
	}

	public double p(ensemble ensemble) {
		double probability = 0.0;
		double beta = 4.0;;

		// Bounds Check
		if (isNotInBox(position))
			return probability;

		/*double hamiltonian = (ensemble.kineticEnergy() + ensemble.interactioPotential() + ensemble.potentialEnergy())
				/ ensemble.nParticles;*/
		double hamiltonian = kineticEnergy() + (ensemble.interactioPotential())/ensemble.nParticles + position.z;
		//System.out.println("p Hamiltonian:["+hamiltonian+"]");
		// Boltzmann
		probability = Math.pow(Math.E, (-beta * hamiltonian));
		return probability;
	}
	public double p_future(ensemble_gpu ensemble, int i) {
		double s = 0.0;
		double beta = 4.0;
		// Bounds Check
		if (isNotInBox(future_position))
			return s;
		if (ensemble.tooCloseTrial(i))
			return s;
		
		/*double hamiltonian = (ensemble.kineticEnergyTrial(i) + ensemble.interactioPotentialTrial(i) + ensemble
				.potentialEnergyTrial(i)) / ensemble.nParticles;*/
		double hamiltonian = kineticEnergyTrial() + ensemble.interactioPotentialTrial(i)/ensemble.nParticles + future_position.z;
		 //System.out.println("p_trial Hamiltonian:["+hamiltonian+"] ensemble.kineticEnergyTrial(i):["+ensemble.kineticEnergyTrial(i)+"] ensemble.interactioPotentialTrial(i):["+ensemble.interactioPotentialTrial(i)+"] ensemble.potentialEnergyTrial(i)"+ensemble.potentialEnergyTrial(i)+"] ");

		// Boltzmann
		s = Math.pow(Math.E, (-beta * hamiltonian));
		return s;
	}
	
	public double p(ensemble_gpu ensemble) {
		double probability = 0.0;
		double beta = 4.0;;

		// Bounds Check
		if (isNotInBox(position))
			return probability;

		/*double hamiltonian = (ensemble.kineticEnergy() + ensemble.interactioPotential() + ensemble.potentialEnergy())
				/ ensemble.nParticles;*/
		double hamiltonian = kineticEnergy() + (ensemble.interactioPotential())/ensemble.nParticles + position.z;
		//System.out.println("p Hamiltonian:["+hamiltonian+"]");
		// Boltzmann
		probability = Math.pow(Math.E, (-beta * hamiltonian));
		return probability;
	}
	public double p_future(ensemble ensemble, int i) {
		double s = 0.0;
		double beta = 4.0;
		// Bounds Check
		if (isNotInBox(future_position))
			return s;
		if (ensemble.tooCloseTrial(i))
			return s;
		
		/*double hamiltonian = (ensemble.kineticEnergyTrial(i) + ensemble.interactioPotentialTrial(i) + ensemble
				.potentialEnergyTrial(i)) / ensemble.nParticles;*/
		double hamiltonian = kineticEnergyTrial() + ensemble.interactioPotentialTrial(i)/ensemble.nParticles + future_position.z;
		 //System.out.println("p_trial Hamiltonian:["+hamiltonian+"] ensemble.kineticEnergyTrial(i):["+ensemble.kineticEnergyTrial(i)+"] ensemble.interactioPotentialTrial(i):["+ensemble.interactioPotentialTrial(i)+"] ensemble.potentialEnergyTrial(i)"+ensemble.potentialEnergyTrial(i)+"] ");

		// Boltzmann
		s = Math.pow(Math.E, (-beta * hamiltonian));
		return s;
	}
	public boolean accepted(ensemble_gpu ensemble_gpu, int i) {
		Random rand = new Random();
		// System.out.println("p_future():["+p_future()+"] p():["+p()+"]");
		double p_future = p_future(ensemble_gpu, i);
		double p = p(ensemble_gpu);
		boolean re = ((p_future >= p) || (p_future > (rand.nextDouble() % 1.0) * p));
		// System.out.println("p_future:["+p_future+"] p:["
		// +p+"] return:["+re+"]");
		return re;
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
}
