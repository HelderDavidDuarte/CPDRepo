#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.01

typedef struct particle
{
	double x;
	double y;
	double vx;
	double vy;
	double m;
	long xi;
	long yj;

}PARTICLE;

PARTICLE par[sizeof(long)];

void init_particles(long seed, long ncside, long n_part, long particle_t)
{
    srandom(seed);
    for(long i=0; i < n_part; i++)
    {
        par[i].x = RND0_1;
        par[i].y = RND0_1;
        par[i].vx = RND0_1 / ncside / 10.0;
        par[i].vy = RND0_1 / ncside / 10.0;
        par[i].m = RND0_1 * ncside / (G * 1e6 * n_part);
    }
}

void particle_pos(){
	for(long i=0; i < n_part; i++){
		par[i].xi=floor(par[i].x*ncside);
		par[i].yj=floor(par[i].y*ncside);
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

double* centerofmass (double totalmass){
	double* centerofmass = (double*) malloc(sizeof(double)*2);
		double* M;
		for(long i=0; i<ncside; i++){
			for(long j=0; j<ncside; j++){
				for(long k=0; k<n_part; k++){
					if(par[k].xi==i && par[k].yj==j) M[k]+=par[k].m;
				}
			}
		} //ALTERAR AQUI
	for(long i=0;i<n_part:i++){
		*centerofmass=(par[i].m*par[i].x)/totalmass;
		*(centerofmass+1)=(par[i].m*par[i].y)/totalmass;
	}
	return centerofmass;
}

void main(){

	int seed = scanf("%d", &seed);
	int ncside = scanf("%d", &ncside);
	int n_part = scanf("%d", &n_part);
	int tstep = scanf("%d", &tstep);
	init_particles(seed, ncside, n_part, tstep);
	printf("Finished!\n");

}