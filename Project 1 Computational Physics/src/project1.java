public class project1 {
	//Golden Mean
	public static void main(String[] args) {
		int m, n;
		double  g, phi, phi0, phi1, phina, phinb, phinc;
		g = (Math.sqrt(5.0) - 1.0) / 2.0;

		phi0 = 1.0;
		phi1 = g;
		phinc = g;

		m = 100;

		for (n = 2; n <= m; n++) {
			phi = phi0 - phi1;

			phina = Math.pow(g, n);
			phinb = Math.pow(g, (float) n);
			phinc = g * phinc;

			System.out.printf("%4d%15.7e%15.7e%15.7e%15.7e\n", n, phi, phina, phinb, phinc);

			phi0 = phi1;
			phi1 = phi;
		}

	}
}
