//JMJ
//AMDG

#include <stdio.h>
#include <unistd.h>		//for access()
#include <stdlib.h>		//for the random function, which returns a long int in the range [0, RAND_MAX]
#include <math.h>		//for trig functions
#include <string.h>		//for memset
#include <time.h>		//for seeding the RNG with the current time

#include <omp.h>		//for multithreading

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
inline vec randVec(double L)
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

//Compute the norm (magnitude) of a vector.
inline double vecNorm(vec *v) {
	return sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
}

//Generate a vector, of magnitude ∝ r from z-axis, in x-y plane tangent to point (x,y).
//analagous to a flat rotation curve,
//as the magnitude of the vector doesn't depend on the distance (x, y) is from (0,0)
inline vec tangentVec(double x, double y, double magnitude, double radius)
{
    vec v = { .z = 0 };
    double theta = atan((y<0?-y:y)/(x<0?-x:x)),
	   r = sqrt(x*x+y*y);
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
    v.x = x*magnitude*r/radius;
    v.y = y*magnitude*r/radius;
    return v;
}

//Generate a random distribution of numObjects objects in sphere of radius rad.
//They will have velMag worth of rotational velocity ∥ to z axis.
obj* generate(int numObjects, double rad, double velMag)
{
    obj *objects = calloc(numObjects, sizeof(obj));
#pragma omp parallel for
    for (int i = 0; i < numObjects; i++) {
	objects[i].pos = randVec(rad);
	objects[i].vel = tangentVec(objects[i].pos.x,
		       		    objects[i].pos.y,
				    velMag, rad);
    }
    return objects;
}

//Compute the norm².
inline double vecNormSqrd(vec *v) {
	return v->x * v->x + v->y * v->y + v->z * v->z;
}

inline void vecPlusEqual(vec *result, vec *addend) {
	result->x += addend->x;
	result->y += addend->y;
	result->z += addend->z;
}

inline void vecDivEqual(vec *result, double divend) {
	result->x /= divend;
	result->y /= divend;
	result->z /= divend;
}

inline void crossProdPlusEqual(vec *result, vec *a, vec *b) {
	result->x += a->y * b->z - a->z * b->y;
	result->y += a->z * b->x - a->x * b->z;
	result->z += a->x * b->z - a->z * b->x;
}

inline vec crossProd(vec *a, vec *b) {
	vec result;
	result.x = a->y * b->z - a->z * b->y;
	result.y = a->z * b->x - a->x * b->z;
	result.z = a->x * b->z - a->z * b->x;
	return result;
}

void printPhysicsStats(obj *objects, int numObjects, FILE *outFile, double t /*= current time step*/) {
	int i;
	vec CoMpos = {0}, CoMvel = {0}, CoMacc = {0}, //center of mass pos., vel., & acc. vectors
	    angularMo = {0}, torque = {0}; //angular momentum (r×v) and torque (r×a)
	double avgD = 0, avgDstd = 0, rmsD = 0, rmsDstd = 0, //distance-from-origin stats
	       totalVel = 0, totalAcc = 0,
	       avgVel = 0, avgVelStd = 0, rmsVel = 0, rmsVelStd = 0,
	       avgAcc = 0, avgAccStd = 0, rmsAcc = 0, rmsAccStd = 0;
#pragma omp parallel
	{ //compute CoMs
		for (i = 0; i < numObjects; i++) {
			//position
			vecPlusEqual(&CoMpos, &objects[i].pos);
			//velocity
			vecPlusEqual(&CoMvel, &objects[i].vel);
			//acceleration
			vecPlusEqual(&CoMacc, &objects[i].acc);
		}
		vecDivEqual(&CoMpos, numObjects);
		vecDivEqual(&CoMvel, numObjects);
		vecDivEqual(&CoMacc, numObjects);
	}
#pragma omp parallel
	{ //compute total pos (distance-from-origin), vel, & acc:
		for (i = 0; i < numObjects; i++) {
			avgD += vecNorm(&objects[i].pos);
			totalVel += vecNorm(&objects[i].vel);
			totalAcc += vecNorm(&objects[i].acc);
		}
		//compute average vel & acc:
		avgD /= (double)numObjects;
		avgVel = totalVel / (double)numObjects;
		avgAcc = totalAcc / (double)numObjects;
	}
#pragma omp parallel
	{ //compute RMSs:
		for (i = 0; i < numObjects; i++) {
			rmsD += vecNormSqrd(&objects[i].pos);
			rmsVel += vecNormSqrd(&objects[i].vel);
			rmsAcc += vecNormSqrd(&objects[i].acc);
		}
		rmsD = sqrt(rmsD/(double)numObjects);
		rmsVel = sqrt(rmsVel/(double)numObjects);
		rmsAcc = sqrt(rmsAcc/(double)numObjects);
		//compute standard dev.s:
		for (i = 0; i < numObjects; i++) {
			//for means
			avgDstd += pow(vecNormSqrd(&objects[i].pos) - rmsD, 2.0);
			avgVelStd += pow(vecNormSqrd(&objects[i].vel) - rmsVel, 2.0);
			avgAccStd += pow(vecNormSqrd(&objects[i].acc) - rmsAcc, 2.0);
			//for RMSs
			rmsDstd += pow(vecNormSqrd(&objects[i].pos) - rmsD, 2.0);
			rmsVelStd += pow(vecNormSqrd(&objects[i].vel) - rmsVel, 2.0);
			rmsAccStd += pow(vecNormSqrd(&objects[i].acc) - rmsAcc, 2.0);
		}
		avgDstd = sqrt(rmsDstd/(double)numObjects);
		avgVelStd = sqrt(rmsVelStd/(double)numObjects);
		avgAccStd = sqrt(rmsAccStd/(double)numObjects);
		rmsDstd = sqrt(rmsDstd/(double)numObjects);
		rmsVelStd = sqrt(rmsVelStd/(double)numObjects);
		rmsAccStd = sqrt(rmsAccStd/(double)numObjects);
	}
#pragma omp parallel
	{ //compute angular mo. and torque
		for (i = 0; i < numObjects; i++) {
			crossProdPlusEqual(&angularMo, &objects[i].pos, &objects[i].vel); // angular mo. = r × v
			crossProdPlusEqual(&torque, &objects[i].pos, &objects[i].acc); // torque = r × a
			vecDivEqual(&angularMo, numObjects);
			vecDivEqual(&torque, numObjects);
		}
	}
	fprintf(outFile, "% e "
			 "% e % e % e "
			 "% e % e % e "
			 "% e % e % e "
			 "% e % e % e "
			 "% e % e % e "
			 "% e % e % e % e "
			 "% e % e % e % e % e "
			 "% e % e % e % e % e\n",
			 t,
			 CoMpos.x, CoMpos.y, CoMpos.z,
			 CoMvel.x, CoMvel.y, CoMvel.z,
			 CoMacc.x, CoMacc.y, CoMacc.z,
			 angularMo.x, angularMo.y, angularMo.z,
			 torque.x, torque.y, torque.z,
			 avgD, avgDstd, rmsD, rmsDstd,
			 totalVel, avgVel, avgVelStd, rmsVel, rmsVelStd,
			 totalAcc, avgAcc, avgAccStd, rmsAcc, rmsAccStd);
}

void printOut(obj *objects, int numObjects, FILE *outFile, double t) {
	vec pos, vel, acc, angmo, torque;
	//print timestep
	fprintf(outFile, "# t = %f\n", t);
	//print physics stats (total KE, RMS vel./acc., etc.) for time-snapshot
	//print header
	fprintf(outFile, "%6s "
			 "%13s %13s %13s %13s "
			 "%13s %13s %13s %13s "
			 "%13s %13s %13s %13s "
			 "%13s %13s\n"
			 "#%5d "
			 "%13d %13d %13d %13d "
			 "%13d %13d %13d %13d "
			 "%13d %13d %13d %13d "
			 "%13d %13d\n",
			 "#   ID",
			 "x", "y", "z", "dist",
			 "vel_x", "vel_y", "vel_z", "|vel|",
			 "acc_x", "acc_y", "acc_z", "|acc|",
			 "ang_mo", "torque",
			 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);
	//print object data
	for (int i = 0; i < numObjects; i++) {
		pos = objects[i].pos;
		vel = objects[i].vel;
		acc = objects[i].acc;
		angmo = crossProd(&pos, &vel), // =r×v
		torque = crossProd(&pos, &acc); // =r×a
		fprintf(outFile,
				"%6d "
			 	"% e % e % e % e "
			 	"% e % e % e % e "
			 	"% e % e % e % e "
			 	"% e % e\n",
				i, //id of object
				pos.x, pos.y, pos.z, vecNorm(&pos),
				vel.x, vel.y, vel.z, vecNorm(&vel),
				acc.x, acc.y, acc.z, vecNorm(&acc),
				vecNorm(&angmo), // =|r×v|
				vecNorm(&torque)); // =|r×a|
	}
}

void writeDataFile(obj *objects, int numObjects, int stepNumber, double t) {
	FILE *f, *s; //frame (f) and stats (s) files
	char str[10];
	//Output "frame" file
	sprintf(str, "%06d.txt", stepNumber);
	if ((f = fopen(str, "w")))
		printOut(objects, numObjects, f, t);
	else
		goto error;
	//Append to stats file
	sprintf(str, "stats.txt");
	if(access(str, F_OK)==0) { //stats.txt exists.
		if (!(s = fopen(str, "a"))) //a = append
			goto error;
	} else { //stats.txt doesn't exist; prepend header first
		if ((s = fopen(str, "a"))) //a = append
			//header:
			fprintf(s, //column headers:
				   "%13s "
				   "%13s %13s %13s "
				   "%13s %13s %13s "
				   "%13s %13s %13s "
				   "%13s %13s %13s "
				   "%13s %13s %13s "
				   "%13s %13s %13s %13s "
				   "%13s %13s %13s %13s %13s "
				   "%13s %13s %13s %13s %13s\n"
				   //column numbers:
				   "#%12d "
				   "%13d %13d %13d "
				   "%13d %13d %13d "
				   "%13d %13d %13d "
				   "%13d %13d %13d "
				   "%13d %13d %13d "
				   "%13d %13d %13d %13d "
				   "%13d %13d %13d %13d %13d "
				   "%13d %13d %13d %13d %13d\n",
				   "#           t",
				   "CoM x", "CoM y", "CoM z",
				   "CoM vel x", "CoM vel y", "CoM vel z",
				   "CoM acc x", "CoM acc y", "CoM acc z",
				   "ang. mo. x", "ang. mo. y", "ang. mo. z",
				   "torque x", "torque y", "torque z",
				   "mean dist", "+/-",  "RMS dist", "+/-",
				   "total vel", "mean vel", "+/-", "RMS vel", "+/-",
				   "total acc", "mean  acc", "+/-", "RMS acc", "+/-",
				   1,
				   2,3,4,
				   5,6,7,
				   8,9,10,
				   11,12,13,
				   14,15,16,
				   17,18,19,20,
				   21,22,23,24,25,
				   26,27,28,29,30);
		else
			goto error;
	}
	printPhysicsStats(objects, numObjects, s, t);
	fclose(f);
	fclose(s);
	return;
error:
	fprintf(stderr, "Could not open file %s for writing\n"
		"Perhaps there is no space left on the device.", str);
}

//find the Euclidean/Pythagorean distance r between two points represented by pos. vectors
inline double distance(vec *a, vec *b) {
	return sqrt((a->x - b->x)*(a->x - b->x)+
		    (a->y - b->y)*(a->y - b->y)+
		    (a->z - b->z)*(a->z - b->z));
}

void integrateWeber(obj *objects, int numObjects, double dt) {
	int i, j;
	double r = 0, r3 = 0, dr_dt = 0, dr_dt_sqrd, d2r_dt2 = 0;
	obj *a, *b;
	obj *objectsOld = malloc(sizeof(obj)*numObjects);
	objectsOld = memcpy(objectsOld, objects, sizeof(obj)*numObjects); //backup old one
	//compute force on object i due to all objects ≠ i
#pragma omp parallel for
	for (i = 0; i < numObjects; i++) {
		a = objects + i; //"a" is the object we're updating, so take it from "objects"
		//rezero acceleration because acceleration is only for current time step 
		a->acc.x = 0;
		a->acc.y = 0;
		a->acc.z = 0;
		for (j = 0; j < numObjects; j++) { // 2 for loops → O(n²)
			if (j != i) { //exclude self-interactions
				b = objectsOld + j; //"b" is the "source" object, which we're not updating, so take it from "objectsOld"
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
				//use Weber's force law here:                   V---- Negative sign here should make the force attractive
				a->acc.x += ((1 - 0.5*dr_dt_sqrd + r*d2r_dt2) * (b->pos.x - a->pos.x))/r3;
				a->acc.y += ((1 - 0.5*dr_dt_sqrd + r*d2r_dt2) * (b->pos.y - a->pos.y))/r3;
				a->acc.z += ((1 - 0.5*dr_dt_sqrd + r*d2r_dt2) * (b->pos.z - a->pos.z))/r3;
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
	free(objectsOld);
}

void integrateNewton(obj *objects, int numObjects, double dt) {
	int i, j;
	double r = 0, r3 = 0;
	obj *a, *b;
	obj *objectsOld = malloc(sizeof(obj)*numObjects);
	objectsOld = memcpy(objectsOld, objects, sizeof(obj)*numObjects); //backup old one
	//compute force on object i due to all objects ≠ i
#pragma omp parallel for
	for (i = 0; i < numObjects; i++) {
		a = objects + i; //"a" is the object we're updating, so take it from "objects"
		//rezero acceleration because acceleration is only for current time step 
		a->acc.x = 0;
		a->acc.y = 0;
		a->acc.z = 0;
		for (j = 0; j < numObjects; j++) { // 2 for loops → O(n²)
			if (j != i) { //exclude self-interactions
				b = objectsOld + j; //"b" is the "source" object, which we're not updating, so take it from "objectsOld"
				r = distance(&a->pos, &b->pos);
					r3 = r*r*r;
				//units: c = e = e' = 1
				//update acceleration
				//use Newton's force law here:
				a->acc.x += (b->pos.x - a->pos.x)/r3;
				a->acc.y += (b->pos.y - a->pos.y)/r3;
				a->acc.z += (b->pos.z - a->pos.z)/r3;
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
	free(objectsOld);
}

void simulate(double t_0, double t_f, double dt, obj *objects, int numObjects, void (*integrate)(obj*, int, double)) {
	double t = t_0; //current time
	int stepNumber = 0, totalSteps = (int)((t_f - t_0)/dt);
	//output initial condition:
	writeDataFile(objects, numObjects, stepNumber++, t);
	for (stepNumber = 1; stepNumber <= totalSteps; stepNumber++) {
		//integrate
		integrate(objects, numObjects, dt);
		//write data file
		writeDataFile(objects, numObjects, stepNumber, t += dt);
	}
}

int main(int argc, char *argv[]) {
	int n = 0; //number of objects
	double radius = 0, //radius of "Plummer" sphere
	       velMag = 0; //magnitude of initial velocity for each object
	double t_0 = 0, t_f = 0, dt = 0;
	_Bool weber = 1; // ≠1 → Newton's force law 
	omp_set_num_threads(omp_get_num_procs()); //# threads == # CPUs
    	srandom(999); //Seed the RNG with time
	for (int i = 1; i < argc; i++) {
		if (strcmp(&argv[i][1], "h")==0)
			goto usage;
		if (!(argc >= 13 && argc <= 16)) {
			fprintf(stderr, "Wrong number of arguments\n");
			goto usage;
		}
        	if (argv[i][0] == '-') {
        		switch (argv[i][1]) {
        		case 'n':
        			n = atoi(argv[++i]);
				if (n <= 0) {
					fprintf(stderr, "(n = %d) ≯ 0\n", n);
					return 1;
				}
        			break;
        		case 'r':
        			radius = atof(argv[++i]);
				if (radius <= 0) {
					fprintf(stderr, "(r = %f) ≯ 0\n", radius);
					return 1;
				}
        			break;
        		case 'v':
        			velMag = atof(argv[++i]);
				if (velMag < 0) {
					fprintf(stderr, "(v = %f) must be ≥ 0.\n", velMag);
					return 1;
				}
        			break;
        		case 'i':
        			t_0 = atof(argv[++i]);
				if (t_0 < 0) {
					fprintf(stderr, "(t_0 = %f) must be > 0.\n", t_0);
					return 1;
				}
        			break;
        		case 'f':
        			t_f = atof(argv[++i]);
				if (t_f <= 0) {
					fprintf(stderr, "(f = %f) ≯ 0\n", t_f);
					return 1;
				}
				if (t_f <= t_0) {
					fprintf(stderr, "(f = %f) ≮ (i = %f)\n", t_0, t_f);
					return 1;
				}
        			break;
        		case 'd':
        			dt = atof(argv[++i]);
				if (dt <= 0) {
					fprintf(stderr, "(d = %f) ≯ 0\n", dt);
					return 1;
				}
        			break;
			case 'w':
				weber = atoi(argv[++i]);
				if (!(weber == 0 || weber == 1)) {
					fprintf(stderr, "(w = %d) ≠ (0 or 1)\n", weber);
					return 1;
				}
				break;
			case 'h':
usage:
				printf("Usage:\n"
				       "\t-n: number of objects\n"
				       "\t-r: radius of \"Plummer\" sphere\n"
				       "\t-v: initial velocity for object on equator\n"
				       "\t\tThe initial velocity is in the x-y plane and ⟂ r.\n"
				       "\t\tIt is ∝ (dist. from z-axis)×(value given).\n"
				       "\t-i: initial time\n"
				       "\t-f: final time\n"
				       "\t-d: time step\n"
				       "Units are defined in terms of c = 1 and e = e' = 1.\n"
				       "Optional:\n\t-w: 0 = Newton's force (default: 1 = Weber's)\n"
				       "\t-h: this help message\n");
				return 0;
			}
		}
	}
	printf("Generating random distribution of"
	       "\n\t%e objects each with"
	       "\n\tvelocity %e ∝ (dist. from z-axis)"
	       "\n\tin sphere of radius %e…\n",
	       (double)n, velMag, radius);
	obj *objects = generate(n, radius, velMag);
	void (*integrationFunction)(obj*, int, double) = (weber==1)?integrateWeber:integrateNewton;
	printf("Integrating %s's force law using %d OpenMP procs…"
	       "\n\tt = %e → %e\n\tΔt = %e\n",
	       (weber==1)?"Weber":"Newton", omp_get_num_procs(),
	       t_0, t_f, dt);
	simulate(t_0, t_f, dt, objects, n, integrationFunction);
	printf("%e files written.\n", (double)((int)((t_f - t_0)/dt)));
	free(objects);
	return 0;
}
