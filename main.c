#include <stdio.h>
#include <stdlib.h>		//for the random function, which returns a long int in the range [0, RAND_MAX]
#include <math.h>		//for trig functions
#include <string.h>		//for memset

//random number in range [-1, 1] and (0, 1], respectively
#define rand11() (2.0*((double)random())/((double)(RAND_MAX))-1.0)
#define rand01() (((double)random())/((double)(RAND_MAX)))

typedef struct {
    double x, y, z;
} vec;

typedef struct {
    vec pos, vel, acc;
    //mass assumed = 1
    //speed of light c also assumed = 1, for simplicity's sake
} obj;

//Generate a vector in a random direction, of random length ≤ L
vec randVec(double L)
{
    double x = L, y = L, z = L;
    while (sqrt(x * x + y * y + z * z) > L) {
	x = rand11() * L;
	y = rand11() * L;
	z = rand11() * L;
    }
    vec vec = {x, y, z};
    return vec;
}

//Generate a unit vector in x-y plane tangent to point (x,y).
//analagous to a flat rotation curve,
//as the magnitude of the vector doesn't depend on the distance (x, y) is from (0,0)
vec tangentUnitVec(double x, double y)
{
    vec v = { .z = 0 };
    double theta = atan((y<0?-y:y)/(x<0?-x:x));
    if (x > 0) {
	if (y > 0) { //1st quad.
	    x = -cos(theta);
	    y = sin(theta);
	} else { //4th quad.
	    x = cos(theta);
	    y = sin(theta);
	}
    } else {
	if (y > 0) { //2nd quad.
	    x = -cos(theta);
	    y = -sin(theta);
	} else { //3rd quad.
	    x = cos(theta);
	    y = -sin(theta);
	}
    }
    v.x = x;
    v.y = y;
    return v;
}

//Generate a random distribution of numObjects objects in sphere of radius rad.
//They will have 1 unit of rotational velocity ∥ to z axis.
obj* generate(int numObjects, double rad)
{
    obj *objects = calloc(numObjects, sizeof(obj));

    for (int i = 0; i < numObjects; i++) {
	objects[i].pos = randVec(rad);
	objects[i].vel = tangentUnitVec(objects[i].pos.x,
					objects[i].pos.y);
    }
    return objects;
}

//Compute the norm (magnitude) of a vector.
double vecNorm(vec *v) {
	return sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
}

void printPhysicsStats(obj *objects, int numObjects, FILE *outFile) {
	vec CoM = {0}; //center of mass position vector
	double totalVel = 0, avgVel = 0, totalAcc = 0, avgAcc = 0;
	for (int i = 0; i < numObjects; i++) {
		totalVel += vecNorm(&objects[i].vel);
		totalAcc += vecNorm(&objects[i].acc);
	}
	avgVel = totalVel / (double)numObjects;
	avgVel = totalVel / (double)numObjects;
	fprintf(outFile, "# Total Velocity: % e  Average Velocity: % e   "
			 "Total Acceleration: % e  Average Acceleration: % e\n",
			 totalVel, avgVel, totalAcc, avgAcc);
}

void printOut(obj *objects, int numObjects, FILE *outFile) {
	vec pos, vel, acc;
	//print physics stats (total KE, RMS vel./acc., etc.) for time-snapshot
	printPhysicsStats(objects, numObjects, outFile);
	//print header
	fprintf(outFile, "%3s  %13s %13s %13s  %13s %13s %13s  %13s %13s %13s\n",
			 "#ID",
			 "x", "y", "z",
			 "vx", "vy", "vz",
			 "ax", "ay", "az");
	//print object data
	for (int i = 0; i < numObjects; i++) {
		pos = objects[i].pos;
		vel = objects[i].vel;
		acc = objects[i].acc;
		fprintf(outFile, "%3d  % e % e % e  % e % e % e  % e % e % e\n",
				i, //id of object
				pos.x, pos.y, pos.z,
				vel.x, vel.y, vel.z,
				acc.x, acc.y, acc.z);
	}
}

int main() {
	int n = 8; //number of objects
	double radius = 128; //radius of "Plummer" sphere

    	srandom(23483920); //Seed the RNG.

	obj *objects = generate(n, radius);

	double t_0, t_f, dt;


	printOut(objects, n, stdout);

	return 0;
}
