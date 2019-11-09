/**
 * Computer Modelling, Nbody Simulation project: 
 *
 * This uses velocity Verlet time integration algorithm to model 
 * the orbits of N bodies interacting with one another
 * via the gravitational force F = -(m1m2/|r|^2)r
 *
 * @author D. Xiong
 * @author T. Bui
 * @version "16/03/2016"
 */

// IO package is rquired for file writing
import java.io.* ;
import java.util.Scanner ; 
import java.lang.Math ;

public class NbodyVerlet  {
    
    /* **********************************************************************
     * Methods to compute gravitational force vector and energy of the system
     * **********************************************************************
     */

    /** Method to calculate gravitational force acting on particle p1 by particle p2
     * @param GRAV_CONSTANT the gravitational constant in the chosen units
     * @param p1 the first particle
     * @param p2 the second particle
     * @return force a Vector3D representing the gravitational force vector from p2 to p1
     */
    static Vector3D gravitationalForce(double GRAV_CONSTANT, Particle3D p1, Particle3D p2){
	Vector3D separation = Particle3D.getSeparation(p1, p2) ;
	double norm = separation.norm();
	double massProduct = p1.getMass() * p2.getMass() ;
	double k = -GRAV_CONSTANT * massProduct / Math.pow(separation.norm(), 3) ;
	Vector3D force = Vector3D.multVector(separation, k) ;
	return force ;
    }

    /** Method to calculate potential energy between two particles
     * @param GRAV_CONSTANT the gravitational constant in the chosen units
     * @param p1 the first particle
     * @param p2 the second particle
     * @return potentialEnergy a double representing the gravitational potential energy between two particles
     */
    static double potentialEnergy(double GRAV_CONSTANT, Particle3D p1, Particle3D p2){
	double distance = Particle3D.getDistance(p1, p2) ;
	double massProduct = p1.getMass() * p2.getMass() ;
	double PE = - GRAV_CONSTANT * massProduct / distance ;
	return PE;
    }

    /** Method to calculate the total energy of the whole system
     * @param GRAV_CONSTANT the gravitational constant in the chosen units
     * @param particles an array holding all particles of the system
     * @return totalEnergy a double representing the sum of kinetic and potential energy of the system
     */
    static double totalEnergy(double GRAV_CONSTANT, Particle3D[] particles) {
	double totalKinetic = 0;
	double totalPotential = 0;
	for (int i = 0; i < particles.length; i++) {
	    totalKinetic += particles[i].kineticEnergy();
	}
	for (int j = 0; j < particles.length-1; j++) {
	    for (int k = j + 1; k < particles.length; k++) {
		totalPotential += potentialEnergy(GRAV_CONSTANT, particles[j], particles[k]);
	    }
	}
	double totalEnergy = totalKinetic + totalPotential;
	return totalEnergy;
    }

    /** Method to calculate the total gravitational force acting on one particle p1 due to all pairwise interactions
     * between p1 and other particles
     * @param GRAV_CONSTANT the gravitational constant in the chosen units
     * @param label an integer suggesting the target particle
     * @param particles an array holding all particles in the system
     * @return force a Vector3D representing the total gravitational force vector on p1
     */
    static Vector3D totalForce (double GRAV_CONSTANT, int label, Particle3D[] particles) {
	Vector3D totalForce = new Vector3D(0,0,0);
	for (int i = 0; i < particles.length; i++) {
	    if (label == i) continue;
	    totalForce = Vector3D.addVector(totalForce, gravitationalForce(GRAV_CONSTANT, particles[label], particles[i]));
	}

	return totalForce;
    }

    /** Method to update the total forces for the whole array of particles
     *
     * @param GRAV_CONSTANT the gravitational constant in the chosen units
     * @param forces an array holding the current forces on each particle
     * @param particles an array holding all particles in the system
     */
    static void totalForces (double GRAV_CONSTANT, Vector3D[] forces, Particle3D[] particles) {
	for (int i = 0; i < forces.length; i++) {
	    forces[i] = totalForce(GRAV_CONSTANT, i, particles);
	}
    }

    /** Method to calculate the Earth-Moon mass centre
     *
     * @param Earth the particle representing the Earth
     * @param Moon the particle representing the Moon
     * @return centreEarthMoon the position vector of Earth-Moon mass centre
     */
    //static Vector3D centreEarthMoon (Particle3D Earth, Particle3D Moon) {
    //	Vector3D temp1 = Vector3D.multVector(Earth.getPosition(), Earth.getMass());
    //	Vector3D temp2 = Vector3D.multVector(Moon.getPosition(), Moon.getMass());
    //	Vector3D temp3 = Vector3D.subVector(temp1, temp2);
    //	Vector3D centreEarthMoon = Vector3D.divideVector(temp3, (Earth.getMass() - Moon.getMass()));
    //	return centreEarthMoon;
    //}

    /* *****************************************************
     * Main method
     * As we are doing file IO we need to throw an exception
     * *****************************************************
     */
    public static void main (String[] argv) throws IOException {

        // Open the trajectory output file
        String trajectoryOut = argv[0];
        PrintWriter out1 = new PrintWriter(new FileWriter(trajectoryOut));

	// Open the energy ouput file
	String energyOut = argv[1];
	PrintWriter out2 = new PrintWriter(new FileWriter(energyOut));

        /*
         * Read in the parameters of the simulations and 
	 * number of particles, N, and initial states of N particles from a file 
	 * using a Scanner.
         */
	String settings = argv[2];
	BufferedReader settingsFile = new BufferedReader(new FileReader(settings));
	Scanner scan1 = new Scanner(settingsFile);

	// Read in number of steps, size of time step and  gravitational constant
	int step = scan1.nextInt();
	double dt = scan1.nextDouble();
	final double GRAV_CONSTANT = scan1.nextDouble();

	// Open up the file of the entries' information
	String entries = argv[3] ;
	BufferedReader entriesFile = new BufferedReader(new FileReader(entries));

	// Attach a Scanner to the file to parse the strings/ numbers
	Scanner scan2 = new Scanner(entriesFile);

	// Read in the number of particles
	int N = scan2.nextInt();

	//Create a new array to hold N particles
	Particle3D particles[] = new Particle3D[N];

	//Read in the N particles from the input file 2
	for (int i = 0; i < N; i++){
	    particles[i] = new Particle3D();
	    particles[i].scanner(scan2) ;
	}

	double t = 0;
	double totalEnergy = totalEnergy(GRAV_CONSTANT, particles);
	Vector3D forces[] = new Vector3D[N];
	for (int w = 0; w < forces.length; w++) {
	    forces[w] = new Vector3D(0,0,0);
        }
	totalForces(GRAV_CONSTANT, forces, particles);

	// Print out the initial position and energy
	out1.printf("%d\n" + "Point = " + "%d\n", N, 1);
	for (int j = 0; j < particles.length; j++) {
	    out1.println(particles[j]);
	}
	out2.printf("%f %10.5f\n", t, totalEnergy);


	/* ****************************************
	 * Algorithm to evolve the system with time
	 * and count number of orbits completed
	 ******************************************
	 */

	// Create an array to hold the number of orbits completed by N-bodies
	double orbits[] = new double[N];
	double periods[] = new double[N];

	// Create an array to hold the postions of N-bodies in the previous time step
	Vector3D prePositions[] = new Vector3D[N];
	double perihelions[] = new double[N];
	double aphelions[] = new double[N];
	for (int p = 0; p < N; p++){
	    prePositions[p] = new Vector3D(particles[p].getPosition());
	    perihelions[p] = Vector3D.subVector(prePositions[p], prePositions[0]).norm();
	    aphelions[p] = 0;
	}
	
	// Create an array to accumulate the angles that N-bodies have moved through
	double anglesInTime[] = new double[N];
	double twoPi[] = new double[N];
	Vector3D forcesNew[] = new Vector3D[N];
	for (int f1 = 0; f1 < N; f1++) {
	    forcesNew[f1] = new Vector3D(0,0,0);
	}

	// Time integration
	int endCount = 0;
	int periodCount[] = new int[N];
	int stepCount = 0;

	while (endCount != N-1 || stepCount < step) {
	    Particle3D.leapPositionArray(dt, particles, forces);
	    totalForces(GRAV_CONSTANT, forcesNew, particles);
	    Particle3D.leapVelocityArray(dt, particles, forces, forcesNew);
	    for (int i = 0; i < forcesNew.length; i++){
		forces[i] = forcesNew[i];
            }
	   
	    for (int l = 1; l < particles.length; l++){
		double tempAngle = 0;
		if (l == 10) {		   
		    Vector3D tempMoonEarth = Vector3D.subVector(particles[3].getPosition(), particles[l].getPosition());
		    Vector3D tempMoonEarthPre = Vector3D.subVector(prePositions[3], prePositions[l]);
		    if (tempMoonEarth.norm() > aphelions[l]) aphelions[l] = tempMoonEarth.norm();
		    if (tempMoonEarth.norm() < perihelions[l]) perihelions[l] = tempMoonEarth.norm();
		    tempAngle = Vector3D.getAngle(tempMoonEarth, tempMoonEarthPre);
		} else {
		    Vector3D tempSunBody = Vector3D.subVector(particles[0].getPosition(), particles[l].getPosition());
		    Vector3D tempSunBodyPre = Vector3D.subVector(particles[0].getPosition(), prePositions[l]);
		    if (tempSunBody.norm() > aphelions[l]) aphelions[l] = tempSunBody.norm();
		    if (tempSunBody.norm() < perihelions[l]) perihelions[l] = tempSunBody.norm();
		    tempAngle = Vector3D.getAngle(tempSunBody, tempSunBodyPre);
		   
		}
	
		twoPi[l] += tempAngle;
	
		if (stepCount < step) { 
		    anglesInTime[l] += tempAngle;
		}
		if ((twoPi[l] / (periodCount[l] + 1)) < (2 * Math.PI) && (2 * Math.PI) - (twoPi[l] / (periodCount[l] + 1)) < tempAngle) {
		    if (periodCount[l] == 0) { 
			endCount++;
		    }
		    periods[l] = t / (periodCount[l] + 1);
		    periodCount[l]++;
		}
		
	    }
	    for (int o = 1; o < particles.length; o++) {
		prePositions[o] = particles[o].getPosition();
	    }
	    if (stepCount < step) {
		out1.printf("%d\n" + "Point = " + "%d\n", N, stepCount+2);
		for (int j = 0; j < particles.length; j++) {
		    out1.println(particles[j]);
		}
		stepCount++;
		totalEnergy = totalEnergy(GRAV_CONSTANT, particles);
		out2.printf("%f %10.5f\n", t, totalEnergy);  
	    }
	    t += dt;
	}	
	out1.close();
	out2.close();

	// Create 2 arrays to hold the semimajor distance and the periods calculated using Kepler's third law
	double semimajors[] = new double[N];
	double periodsKepler[] = new double[N];
	double c;

	//Calculate the number of orbits turned by dividing the accumulated angles by 2*pi
	for (int z = 0; z < N; z++){
	    orbits[z] = anglesInTime[z]/(2*Math.PI);
	    //Calculate semimajor distance by averaging perihelion and aphelion
	    semimajors[z] = (perihelions[z] + aphelions[z])/2;
	    if (z == 10) {
		c = (2*Math.PI) / Math.sqrt(GRAV_CONSTANT*particles[3].getMass());
	    }
	    else {
		c = (2*Math.PI) / Math.sqrt(GRAV_CONSTANT*particles[0].getMass());
	    }
	    periodsKepler[z] = c * Math.pow(semimajors[z], 1.5);
	    System.out.println(particles[z].getLabel() + "\t" + orbits[z] + "\t" + periods[z] + "\t" + perihelions[z] + "\t" + aphelions[z] + "\t" + periodsKepler[z]);
	}
	
    }
}



