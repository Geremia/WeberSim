//JMJ
//AMDG

#include <stdio.h>
#include <stdlib.h>		//for the random function, which returns a long int in the range [0, RAND_MAX]
#include <math.h>		//for trig functions
#include <string.h>		//for memset
#include <time.h>		//for seeding the RNG with the current time

//random number in range [-1, 1] and (0, 1], respectively
#define rand11() (2.0*((double)random())/((double)(RAND_MAX))-1.0)
#define rand01() (((double)random())/((double)(RAND_MAX)))

typedef struct {
    double x, y, z;
} vec;

typedef struct {
    vec pos, vel, acc;
    //mass assumed = 1
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

//Generate a vector of magnitude magnitude in x-y plane tangent to point (x,y).
//analagous to a flat rotation curve,
//as the magnitude of the vector doesn't depend on the distance (x, y) is from (0,0)
vec tangentVec(double x, double y, double magnitude)
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
    v.x = x*magnitude;
    v.y = y*magnitude;
    return v;
}

//Generate a random distribution of numObjects objects in sphere of radius rad.
//They will have velMag worth of rotational velocity ∥ to z axis.
obj* generate(int numObjects, double rad, double velMag)
{
    obj *objects = calloc(numObjects, sizeof(obj));
    for (int i = 0; i < numObjects; i++) {
	objects[i].pos = randVec(rad);
	objects[i].vel = tangentVec(objects[i].pos.x,
		       		    objects[i].pos.y,
				    velMag);
    }
    return objects;
}

//Compute the norm (magnitude) of a vector.
double vecNorm(vec *v) {
	return sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
}

//Compute the norm².
double vecNormSqrd(vec *v) {
	return v->x * v->x + v->y * v->y + v->z * v->z;
}

void printPhysicsStats(obj *objects, int numObjects, FILE *outFile) {
	int i;
	vec CoMpos = {0}, CoMvel = {0}, CoMacc = {0}; //center of mass pos., vel., & acc. vectors
	double totalVel = 0, totalAcc = 0,
	       avgVel = 0, avgVelStd = 0, rmsVel = 0, rmsVelStd = 0,
	       avgAcc = 0, avgAccStd = 0, rmsAcc = 0, rmsAccStd = 0;
	//compute CoMs
	for (i = 0; i < numObjects; i++) {
		//position
		CoMpos.x += objects[i].pos.x;
		CoMpos.y += objects[i].pos.y;
		CoMpos.z += objects[i].pos.z;
		CoMpos.x /= (double)numObjects;
		CoMpos.y /= (double)numObjects;
		CoMpos.z /= (double)numObjects;
		//velocity
		CoMvel.x += objects[i].vel.x;
		CoMvel.y += objects[i].vel.y;
		CoMvel.z += objects[i].vel.z;
		CoMvel.x /= (double)numObjects;
		CoMvel.y /= (double)numObjects;
		CoMvel.z /= (double)numObjects;
		//acceleration
		CoMacc.x += objects[i].acc.x;
		CoMacc.y += objects[i].acc.y;
		CoMacc.z += objects[i].acc.z;
		CoMacc.x /= (double)numObjects;
		CoMacc.y /= (double)numObjects;
		CoMacc.z /= (double)numObjects;
	}
	//compute total vel & acc:
	for (i = 0; i < numObjects; i++) {
		totalVel += vecNorm(&objects[i].vel);
		totalAcc += vecNorm(&objects[i].acc);
	}
	//compute average vel & acc:
	avgVel = totalVel / (double)numObjects;
	avgAcc = totalAcc / (double)numObjects;
	//compute averages' standard dev.s:
	for (i = 0; i < numObjects; i++) {
		rmsVel += vecNormSqrd(&objects[i].vel);
		rmsVel = sqrt(rmsVel);
		rmsAcc += vecNormSqrd(&objects[i].acc);
		rmsAcc = sqrt(rmsAcc);
	}
	//compute RMSs' standard dev.s:
	for (i = 0; i < numObjects; i++) {
		rmsVelStd += pow(vecNormSqrd(&objects[i].vel) - rmsVel, 2.0);
		rmsVelStd = sqrt(rmsVelStd);
		rmsAccStd += pow(vecNormSqrd(&objects[i].acc) - rmsAcc, 2.0);
		rmsAccStd = sqrt(rmsAccStd);
	}
	fprintf(outFile, "# Center of Mass:\n#\t pos.: (% e, % e, % e)\n#\t vel.: (% e, % e, % e)\n#\t acc.: (% e, % e, % e)\n"
			 "# Speed:\n#\ttotal: %e\n#\t mean: %e ± %e\n#\t  RMS: %e ± %e\n"
			 "# Acceleration:\n#\ttotal: %e\n#\t mean: %e ± %e\n#\t  RMS: %e ± %e\n",
			 CoMpos.x, CoMpos.y, CoMpos.z,
			 CoMvel.x, CoMvel.y, CoMvel.z,
			 CoMacc.x, CoMacc.y, CoMacc.z,
			 totalVel, avgVel, avgVelStd, rmsVel, rmsVelStd,
			 totalAcc, avgAcc, avgAccStd, rmsAcc, rmsAccStd);
}

void printOut(obj *objects, int numObjects, FILE *outFile, double t) {
	vec pos, vel, acc;
	//print timestep
	fprintf(outFile, "# t = %f\n", t);
	//print physics stats (total KE, RMS vel./acc., etc.) for time-snapshot
	printPhysicsStats(objects, numObjects, outFile);
	//print header
	fprintf(outFile, "%6s  %13s %13s %13s  %13s %13s %13s  %13s %13s %13s\n",
			 "#   ID",
			 "x", "y", "z",
			 "vel_x", "vel_y", "vel_z",
			 "acc_x", "acc_y", "acc_z");
	//print object data
	for (int i = 0; i < numObjects; i++) {
		pos = objects[i].pos;
		vel = objects[i].vel;
		acc = objects[i].acc;
		fprintf(outFile, "%6d  % e % e % e  % e % e % e  % e % e % e\n",
				i, //id of object
				pos.x, pos.y, pos.z,
				vel.x, vel.y, vel.z,
				acc.x, acc.y, acc.z);
	}
}

void writeDataFile(obj *objects, int numObjects, int stepNumber, double t) {
	FILE *f;
	char str[10];
	sprintf(str, "%06d.txt", stepNumber);
	if ((f = fopen(str, "w"))) {
		printOut(objects, numObjects, f, t);
		fclose(f);
	} else
		fprintf(stderr, "Could not open file %s for writing\n"
				"Perhaps there is no space left on the device.", str);
}

//find the Euclidean/Pythagorean distance r between two points represented by pos. vectors
double distance(vec *a, vec *b) {
	return sqrt((a->x - b->x)*(a->x - b->x)+
		    (a->y - b->y)*(a->y - b->y)+
		    (a->z - b->z)*(a->z - b->z));
}

void integrate(obj *objects, int numObjects, double dt) {
	int i, j;
	double r = 0, r3 = 0, dr_dt = 0, dr_dt_sqrd, d2r_dt2 = 0;
	obj *a, *b;
	for (i = 0; i < numObjects; i++) {
		//compute force on object i due to all objects ≠ i
		for (j = 0; j < numObjects; j++) { // 2 for loops → O(n²)
			a = objects + i;
			if (j != i) { //exclude self-interactions
				b = objects + j;
				r = distance(&a->pos, &b->pos);
					r3 = r*r*r;
				//units: c = e = e' = 1
				dr_dt = ((a->pos.x - b->pos.x)*(a->vel.x - b->vel.x)+
					 (a->pos.y - b->pos.y)*(a->vel.y - b->vel.y)+
					 (a->pos.z - b->pos.z)*(a->vel.z - b->vel.z))/r;
					dr_dt_sqrd = dr_dt * dr_dt;
				d2r_dt2 = ((a->vel.x - b->vel.x)*(a->vel.x - b->vel.x) + (a->pos.x - b->pos.x)*(a->acc.x - b->acc.x) +
					   (a->vel.y - b->vel.y)*(a->vel.y - b->vel.y) + (a->pos.y - b->pos.y)*(a->acc.y - b->acc.y) +
					   (a->vel.z - b->vel.z)*(a->vel.z - b->vel.z) + (a->pos.z - b->pos.z)*(a->acc.z - b->acc.z)
					   - dr_dt_sqrd)/r;
				//update acceleration
				//use Weber's force law here:
				a->acc.x = ((1 - 0.5*dr_dt_sqrd + r*d2r_dt2) * (b->pos.x - a->pos.x))/r3;
				a->acc.y = ((1 - 0.5*dr_dt_sqrd + r*d2r_dt2) * (b->pos.y - a->pos.y))/r3;
				a->acc.z = ((1 - 0.5*dr_dt_sqrd + r*d2r_dt2) * (b->pos.z - a->pos.z))/r3;
				//update position using current velocity and next time step's acceleration
				a->pos.x += 0.5*a->acc.x*dt*dt + a->vel.x*dt;
				a->pos.y += 0.5*a->acc.y*dt*dt + a->vel.y*dt;
				a->pos.z += 0.5*a->acc.z*dt*dt + a->vel.z*dt;
				//update velocity using next time step's acceleration
				a->vel.x += a->acc.x * dt;
				a->vel.y += a->acc.y * dt;
				a->vel.z += a->acc.z * dt;
			}
		}
	}
}

void simulate(double t_0, double t_f, double dt, obj *objects, int numObjects) {
	double t; //current simulation time
	int stepNumber = 0;
	//output initial condition:
	writeDataFile(objects, numObjects, stepNumber++, 0.);
	for (t = t_0+dt; t < t_f; t += dt) {
		//integrate 
		integrate(objects, numObjects, dt);	
		//write data file
		writeDataFile(objects, numObjects, stepNumber++, t);
	}
}

int main(int argc, char *argv[]) {
	int n = 0; //number of objects
	double radius = 0, //radius of "Plummer" sphere
	       velMag = 0; //magnitude of initial velocity for each object
	double t_0 = 0, t_f = 0, dt = 0;
    	srandom(time(NULL)); //Seed the RNG with the time.
	for (int i = 1; i < argc; i++) {
		if (strcmp(&argv[i][1], "h")==0)
			goto usage;
		if (argc != 13) {
			printf("Wrong number of arguments\n");
			goto usage;
		}
        	if (argv[i][0] == '-') {
        		switch (argv[i][1]) {
        		case 'n':
        			n = atoi(argv[++i]);
				if (n <= 0) {
					printf("(n = %d) ≯ 0\n", n);
					return 1;
				}
        			break;
        		case 'r':
        			radius = atof(argv[++i]);
				if (radius <= 0) {
					printf("(r = %f) ≯ 0\n", radius);
					return 1;
				}
        			break;
        		case 'v':
        			velMag = atof(argv[++i]);
				if (velMag <= 0) {
					printf("(v = %f) ≯ 0\n", velMag);
					return 1;
				}
        			break;
        		case 'i':
        			t_0 = atof(argv[++i]);
				if (t_0 < 0) {
					printf("(t_0 = %f) must be > 0.\n", t_0);
					return 1;
				}
        			break;
        		case 'f':
        			t_f = atof(argv[++i]);
				if (t_f <= 0) {
					printf("(f = %f) ≯ 0\n", t_f);
					return 1;
				}
				if (t_f <= t_0) {
					printf("(f = %f) ≮ (i = %f)\n", t_0, t_f);
					return 1;
				}
        			break;
        		case 'd':
        			dt = atof(argv[++i]);
				if (dt <= 0) {
					printf("(d = %f) ≯ 0\n", dt);
					return 1;
				}
        			break;
			case 'h':
usage:
				printf("Usage:\n"
				       "\t-n: number of objects\n"
				       "\t-r: radius of \"Plummer\" sphere\n"
				       "\t-v: initial velocity for each object\n"
				       "\t\tThe initial velocity is in the x-y plane.\n"
				       "\t-i: initial time\n"
				       "\t-f: final time\n"
				       "\t-d: time step\n"
				       "Units are defined in terms of c = 1 and e = e' = 1.\n"
				       "\t-h: this help message\n");
				return 0;
			}
		}
	}
	obj *objects = generate(n, radius, velMag);
	simulate(t_0, t_f, dt, objects, n);
	return 0;
}
