/**
 * Computer Modelling, Exercise 3: Particle3D class.
 *
 * @author D. Xiong
 * @author T. Bui
 * @version "16/03/2016"
 *
 */
import java.util.Scanner;

public class Particle3D {

    private String label;
    private double mass;
    private Vector3D position;
    private Vector3D velocity;

    /* ***********************
     * Setters and Getters
     * ***********************/    
    /** Get the position of a particle.
     *
     * @return a Vector3D representing the position.
     */
    public Vector3D getPosition() { return position; }

    /** Get the velocity of a particle.
     *
     * @return a double representing the velocity.
     */
    public Vector3D getVelocity() { return velocity; }

    /** Get the mass of a particle.
     *
     * @return a double representing the mass.
     */
    public double getMass()     { return mass; }

    /** Get the label of a particle.
     *
     * @return a String representing the label.
     */
    public String getLabel()     { return label; }

    /** Set the position of a particle.
     *
     * @param r a Vector3D representing the position.
     */    public void setPosition(Vector3D r) {
	this.position = new Vector3D(r);
    }
    
    /** Set the velocity of a particle.
     *
     * @param v a Vector3D representing the velocity.
     */
    public void setVelocity(Vector3D v) {
	this.velocity = new Vector3D(v);
    }
    
    /** Set the mass of a particle.
     *
     * @param m a double representing the mass.
     */
    public void setMass(double m)     { this.mass = m; }

    /** Set the label of a particle.
     *
     * @param l a String representing the label.
     */
 
    public void setLabel(String l)     { this.label = l; }

    /* ******************************************
     * Constructors
     ********************************************/
    
    /** Default constructor. Sets all properties to "not a number" 
     * to indicate that they are uninitialised.
     */
    public Particle3D() {
	this.setLabel("");
        this.setMass(Double.NaN);
        this.setPosition(new Vector3D());
        this.setVelocity(new Vector3D());
    }

    /** Explicit constructor. Constructs a new Particle3D with
     * explicitly given position, velocity, and mass.
     *
     * @param l a String that defines the label.
     * @param m a double that defines the mass.
     * @param r a Vector3D  that defines the position.
     * @param v a Vector3D that defines the velocity.
     */
    public Particle3D(String l, double m, Vector3D r, Vector3D v) {
	this.setLabel(l);
        this.setMass(m);
        this.setPosition(r);
        this.setVelocity(v);
    }

    /* ******************************************
     * toString Method
     ********************************************/
    
    /** Returns a String representation of Particle3D.
     * Used to print a Particle3D instance using the "%s"
     * format identifier.
     * Used to print a Particle3D instance in the form
     * "[label], position: [position]" 
     */
    
    public String toString() {
        return label + " " + position.toString();
    }

    /** Scan the input file given in the command line
     * to get the properties of a Particle3D.
     *
     * @return a Particle3D with properties as
     * given in the input file.
     */
    public void scanner(Scanner scan) {
	this.label = scan.next();
	this.mass = scan.nextDouble();

	this.position.setX(scan.nextDouble());
       	this.position.setY(scan.nextDouble());
       	this.position.setZ(scan.nextDouble());
       
	this.velocity.setX(scan.nextDouble());
        this.velocity.setY(scan.nextDouble());
        this.velocity.setZ(scan.nextDouble());
    }
    

    /* ******************************************
     * Instance Methods
     ********************************************/
    
    /** Returns the kinetic energy of a Particle3D,
     * calculated as 1/2*m*|v|^2.
     *
     * @return a double that is the kinetic energy.
     */
    public double kineticEnergy() { 
	double KE = 0.5*mass*velocity.normSquared();
	return KE;
    }

    /** Time integration support: evolve the velocity
     * according to dv = f/m * dt.
     *
     * @param dt a double that is the timestep.
     * @param force a Vector3D that is the current force on the Particle3D.
     * @param forceNew a Vector3D that is the updated force on the Particle3D after this timestep.
     */
    public void leapVelocity(double dt, Vector3D force, Vector3D forceNew) {
	Vector3D sumForce = new Vector3D(Vector3D.addVector(force, forceNew));
	Vector3D dv = new Vector3D(Vector3D.multVector(sumForce, 0.5*dt/mass));
	velocity = Vector3D.addVector(velocity, dv);
    }
    
    /** Time integration support: evolve the position
     * according to dr = v * dt.
     *
     * @param dt a double that is the timestep.
     */
    public void leapPosition(double dt) {
	Vector3D dr = new Vector3D(Vector3D.multVector(velocity, dt));
        position = Vector3D.addVector(position, dr);
    }

    /** Time integration support: evolve the position
     * according to dr = v * dt + 0.5 * (f/m) * dt**2.
     *
     * @param dt a double that is the timestep.
     * @param force a Vector3D that is the current force.
     */
    public void leapPosition(double dt, Vector3D force) {
       Vector3D dr1 = new Vector3D(Vector3D.multVector(velocity, dt));
       Vector3D dr2 = new Vector3D(Vector3D.multVector(force, (dt*dt)/(2*mass)));
       Vector3D dr = new Vector3D(Vector3D.addVector(dr1, dr2));
       position = Vector3D.addVector(position, dr);
    } 
    
  /* ******************************************
     * Static Methods
     ********************************************/


    /** Manipulate the separation vector between two particles.
     *
     * @param p1 a Particle3D instance representing particle 1.
     * @param p2 a Particle3D instance representing particle 2.
     *
     * @return separation a Vector3D that is the separation vector
     * between two particles, pointing from p2 to p1
     */
    public static Vector3D getSeparation(Particle3D p1, Particle3D p2) {
	Vector3D r1 = new Vector3D(p1.getPosition());
	Vector3D r2 = new Vector3D(p2.getPosition());
	Vector3D separation = new Vector3D(Vector3D.subVector(r1, r2));
	return separation;
    }

    /** Calculate the distance between two particles.
     *
     * @param p1 a Particle3D instance representing particle 1.
     * @param p2 a Particle3D instance representing particle 2.
     * @return distance a double that is the distance
     * between two particles
     */
    public static double getDistance(Particle3D p1, Particle3D p2){
	double distance = getSeparation(p1, p2).norm();
	return distance;
    }

   /** Time integration support: evolve the velocities of an array of Particle3Ds 
     * according to dv = f/m * dt.
     *
     * @param dt a double that is the timestep.
     * @param Nbody a Particle3D array that holds all particles of the system.
     * @param forces a Vector3D array that holds the current forces on each Particle3D of the system.
     * @param forcesNew a Vector3D array that holds the updated forces on each Particle3D after this timestep.
     */
    public static void leapVelocityArray (double dt, Particle3D[] Nbody, Vector3D[] forces, Vector3D[] forcesNew) {
	for (int i = 0; i < Nbody.length; i++) {
	    Nbody[i].leapVelocity(dt, forces[i], forcesNew[i]);
	}
    }

   /** Time integration support: evolve the positions of an array of Particle3Ds
     * according to dr = v * dt + 0.5 * (f/m) * dt**2.
     *
     * @param dt a double that is the timestep.
     * @param Nbody a Particle3D array that holds all particles of the system.
     * @param forces a Vector3D array that holds the current forces on each Particle3D of the system.
     */    
    public static void leapPositionArray (double dt, Particle3D[] Nbody, Vector3D[] forces) {
	for (int i = 0; i < Nbody.length; i++) {
	    Nbody[i].leapPosition(dt, forces[i]);
	}
    }
    
};
