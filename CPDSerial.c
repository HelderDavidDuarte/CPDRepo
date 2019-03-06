#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.01

struct particle
{
	double x;
	double y;
	double vx;
	double vy;
	double m;	
};

void init_particles(long seed, long ncside, long n_part, long particle_t)
{
    long i;

    srandom(seed);

    struct particle par[n_part];

    for(i = 0; i < n_part; i++)
    {
        par[i].x = RND0_1;
        par[i].y = RND0_1;
        par[i].vx = RND0_1 / ncside / 10.0;
        par[i].vy = RND0_1 / ncside / 10.0;

        par[i].m = RND0_1 * ncside / (G * 1e6 * n_part);
    }
}

double* forca (double mass, double accelx, double accely){
	double *force = (double*) malloc(sizeof(double)*2);
	*force = mass*accelx;
	*(force+1) = mass*accely;
	return force;
}

double* velocidade (double prevspeedx, double prevspeedy, double accelx, double accely, int tstep){
	double *speed = (double*) malloc(sizeof(double)*2);
	*speed = prevspeedx + accelx*tstep;
	*(speed+1) = prevspeedy + accely*tstep;
	return speed;
}

double* movement (double xpos, double ypos, double accelx, double accely, double speedx, double speedy, int tstep){
	double *newpos = (double*) malloc(sizeof(double)*2);
	*newpos = xpos + speedx*tstep + accelx*tstep*tstep;
	*(newpos+1) = ypos + speedy*tstep + accely*tstep*tstep;
	return newpos;
}
/*
double centerofmass (){

	return centerofmass;
}*/

void main(){
	int seed = scanf("%d", &seed);
	int ncside = scanf("%d", &ncside);
	int n_part = scanf("%d", &n_part);
	int tstep = scanf("%d", &tstep);
	init_particles(seed, ncside, n_part, tstep);
	printf("Finished!\n");

}