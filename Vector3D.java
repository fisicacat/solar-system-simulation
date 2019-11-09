/** 
 * A class for 3D vectors, complete with constructors,
 * setters and getters, instance methods to calculate scalar multiplication and division, the magnitude squared
 * and the magnitude, and static methods to implement addition, subtraction, dot product and cross product between
 * two 3D vectors. 
 *
 * @author Doreen Xiong 
 * @author Tien Bui
 * @version "16/03/2016"
 *
 */
import java.lang.Math;

public class Vector3D {

    private double X;
    private double Y;
    private double Z;

    /* ************************
     *  Constructors
     * ************************/

    /** Default constructor. Constructs a new Vector3D, with uninitialised
     *  X, Y and Z components.
     */
    public Vector3D () {
        // Set to "not-a-number" to indicate uninitialised
	this.setVector(Double.NaN, Double.NaN, Double.NaN);
    }

    /** Copy constructor. Constructs a new Vector3D by copying the three components
     *  of another Vector3D instance.
     *
     * @param original the Vector3D to be copied
     */
    public Vector3D (Vector3D original) {
	this.setVector(original.getX(), original.getY(), original.getZ());
    }

    /** Explicit constructor. Constructs a new Vector3D from explicitly
     *  given X, Y and Z components.
     *
     * @param xx a double giving the x-component of the new Vector3D
     * @param yy a double giving the y-component of the new Vector3D
     * @param zz a double giving the z-component of the new Vector3D
     */
    public Vector3D(double xx, double yy, double zz) {
	this.setVector(xx, yy, zz);
    }

    /* **********************
     * Setters and getters
     * **********************/

    /** Convenient set method to set all three components.
     *
     * @param xx a double to set the x-component
     * @param yy a double to set the y-component
     * @param zz a double to set the z-component
     */
    public void setVector(double xx, double yy, double zz) {
	this.setX(xx);
	this.setY(yy);
	this.setZ(zz);
    }


    /** Sets the x-component of Vector3D only.
     *
     * @param xx a double to set the x-component
     */
    public void setX(double xx) { 
	this.X = xx; 
	}

    /** Sets the y-component of Vector3D only.
     *
     * @param yy a double to set the y-component
     */
    public void setY(double yy) { 
	this.Y = yy; 
	}

    /** Sets the z-component of Vector3D only.
     *
     * @param zz a double to set the z-component
     */
    public void setZ(double zz) { 
	this.Z = zz; 
	}
 
    /** Gets the x-component of a Vector3D.
     *
     * @return a double instance representing the Vector3D's x-component.
     */
    public double getX() { 
	return this.X; 
	}

    /** Gets the y-component of a Vector3D.
     *
     * @return a double instance representing the Vector3D's y-component.
     */
    public double getY() { 
	return this.Y; 
	}
    
     /** Gets the z-component of a Vector3D.
     *
     * @return a double instance representing the Vector3D's z-component.
     */
    public double getZ() { 
	return this.Z; 
	}

    /* **********************
     * toString method
     * **********************/

    /** Returns a String representation of the Vector3D. Methods 
     * called 'toString' are automagically called when an object
     * is printed.<br>
     *
     * @return a string representation of the Vector3D instance
     */
    public String toString() { 
        double xx = this.getX();
        double yy = this.getY();
	double zz = this.getZ();
        return xx + " " + yy + " " + zz;
    } 

    /* ******************
     * Instance methods
     * ******************/

    /** Calculates the square modulus |A|^2 (norm squared) of Vector3D A.
     * 
     * @return a double representing the Vector3D's square modulus.
     */ 
    public double normSquared() { 
	return this.getX() * this.getX() + 
	       this.getY() * this.getY() + 
	       this.getZ() * this.getZ();  
	}

    /** Calculates the modulus (norm) |A| of the Vector3D A.
     *
     * @return a double representing the Vector3D's modulus.
     */
    public double norm() { 
	return Math.sqrt(this.normSquared());  
	}

    /* *********************
     * Static methods
     * *********************/

    /** Multiplies a Vector3D A = (X,Y,Z) with a double b as A*b = (X*b,Y*b,Z*b).
     *
     * @param A a Vector3D
     * @param b the scalar factor
     * @return the Vector3D A scaled by the real number b.
     */
    public static Vector3D multVector(Vector3D A, double b) {
	return new Vector3D(A.getX()*b, A.getY()*b, A.getZ()*b);
    }
    
    /** Divides Vector3D A = (X, Y, Z) by double b as A/b = (X/b, Y/b, Z/b).
     *
     * @param A a Vector3D
     * @param b the scalar divisor
     * @return the Vector3D divided by a scalar
     */
    public static Vector3D divideVector(Vector3D A, double b) {
	return new Vector3D(A.getX()/b, A.getY()/b, A.getZ()/b);
    }

    /** Adds two Vector3Ds.
     *
     * @param A the first Vector3D to be added
     * @param B the second Vector3D to be added
     * @return the Vector3D sum of A and B, A+B.
     */
    public static Vector3D addVector(Vector3D A, Vector3D B) { 
	return new Vector3D(A.getX() + B.getX(), 
			    A.getY() + B.getY(), 
			    A.getZ() + B.getZ());
    }
    
    /** Subtracts two Vector3Ds.
     *
     * @param A the subtrahend
     * @param B the subtractor
     * @return the Vector3D difference between A and B, A-B.
     */
    public static Vector3D subVector(Vector3D A, Vector3D B) { 
	return new Vector3D(A.getX() - B.getX(), 
			    A.getY() - B.getY(), 
			    A.getZ() - B.getZ());
    }
    
    /** Do dot product of two Vector3Ds A = (X1, Y1, Z1) and B = (X2, Y2, Z2)
     *
     * @param A the first Vector3D in the dot product
     * @param B the second Vector3D in the dot product
     * @return the dot product of A and B.
     */
    public static double doDotProduct(Vector3D A, Vector3D B) {
	double dotProduct = A.getX()*B.getX() + A.getY()*B.getY() + A.getZ()*B.getZ();
	return dotProduct;
    }
    
    
    /** Do cross product of 2 Vector3Ds A = (X1, Y1, Z1) and B = (X2, Y2, Z2)
     *
     * @param A the first Vector3D in the cross product
     * @param B the second Vector3D in the cross product
     * @return angle the cross product of A and B.
     */
    public static Vector3D doCrossProduct(Vector3D A, Vector3D B) {
	return new Vector3D(A.getY()*B.getZ() - A.getZ()*B.getY(),
			    A.getZ()*B.getX() - A.getX()*B.getZ(),
			    A.getX()*B.getY() - A.getY()*B.getX());
    }
 
  /** Compare 2 Vector3Ds A = (X1, Y1, Z1) and B = (X2, Y2, Z2)
     *
     * @param A the first Vector3D
     * @param B the second Vector3D
     * @return a Boolean type to indicate if A and B are equal.
     */
    public static Boolean compareVectors(Vector3D A, Vector3D B) {
	double deltaX = Math.abs(A.getX() - B.getX());
	double deltaY = Math.abs(A.getY() - B.getY());
	double deltaZ = Math.abs(A.getZ() - B.getZ());
	if (deltaX < 0.000000001 && deltaY < 0.000000001 && deltaZ < 0.000000001 ) {
	    return true;
	}
	else {
	    return false;
	}
	
    }

  /** Find the angle between two Vector3Ds
     *
     * @param A the first Vector3D
     * @param B the second Vector3D
     * @return the acute angle between A and B.
     */
    public static double getAngle(Vector3D A, Vector3D B){
	double dotProduct = doDotProduct(A, B);
	double angle = Math.acos(dotProduct/(A.norm()*B.norm()));
	return angle;
    }
    
}

