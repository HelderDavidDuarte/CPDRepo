#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <math.h>

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

typedef struct matrix
{
	double mass;
	double cmx;
	double cmy;
	long xi;
	long yj;

}MATRIX;

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

void particle_pos(long ncside, long n_part){
	for(long i=0; i < n_part; i++){
		par[i].xi=floor(par[i].x*ncside);
		par[i].yj=floor(par[i].y*ncside);
	}
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

void centerofmass (long ncside, long n_part){
	double M[ncside][ncside];
	double x,y;
	MATRIX m;
	for(long i=0; i<ncside; i++){
		for(long j=0; j<ncside; j++){
			for(long k=0; k<n_part; k++){
				if(par[k].xi==i && par[k].yj==j) M[i][j]+=par[k].m;
			}
			x=0, y=0;
			for(long k=0;k<n_part;k++){
				if(par[k].xi==i && par[k].yj==j){
					x+=(par[k].m*par[k].x)/M[i][j]; //centro de massa em x, para uma dada celula
					y+=(par[k].m*par[k].y)/M[i][j]; //centro de massa em y, para uma dada celula
				}
			}
			for(long k=0;k<n_part;k++){
				if(par[k].xi==i && par[k].yj==j){
					m.xi=par[k].xi;
					m.yj=par[k].yj;
					m.cmx=x; //centro de massa em x, para uma dada celula
					m.cmy=y; //centro de massa em y, para uma dada celula
					m.mass=M[i][j];
					printf("%f\n", m.cmx);
					printf("%f\n", m.cmy);
				}
			}
		}
	}
}

double* forca (long k){
	double *force = (double*) malloc(sizeof(double)*2);
	MATRIX n;
	for(long i=0; i<ncside; i++){
		for(long j=0; j<ncside; j++){
			*force = G*n.mass*particle[k].m/((n.cmx-particle[k].x)*(n.cmx-particle[k].x));
			*(force+1) = G*n.mass*particle[k].m/((n.cmy-particle[k].y)*(n.cmy-particle[k].y));
		}
	}	
	return force;
}

double* accel (double mass, double accelx, double accely){
	double *force = (double*) malloc(sizeof(double)*2);
	*force = mass*accelx;
	*(force+1) = mass*accely;
	return force;
}

void main(int argc, char** argv){
	long seed = atoi(argv[1]);
	long ncside = atoi(argv[2]);
	long n_part = atoi(argv[3]);
	long tstep = atoi(argv[4]);

	init_particles(seed, ncside, n_part, tstep);
	particle_pos(ncside, n_part);
	centerofmass(ncside, n_part);
	printf("Finished!\n");

}