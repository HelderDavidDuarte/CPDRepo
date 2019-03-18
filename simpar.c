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

}MATRIX;

PARTICLE par[sizeof(long)];
MATRIX mtr[sizeof(long)][sizeof(long)];

void init_particles(long seed, long ncside, long n_part, long particle_t)
{
	long long i;
    srandom(seed);
    for(i=0; i < n_part; i++)
    {
        par[i].x = RND0_1;
        par[i].y = RND0_1;
        par[i].vx = RND0_1 / ncside / 10.0;
        par[i].vy = RND0_1 / ncside / 10.0;
        par[i].m = RND0_1 * ncside / (G * 1e6 * n_part);
    }
}

void centerofmass (long ncside, long n_part){
	for(long i=0; i<ncside; i++){
		for(long j=0; j<ncside; j++){
			for(long k=0; k<n_part; k++){
				par[k].xi=floor(par[k].x*ncside);
				par[k].yj=floor(par[k].y*ncside);
				if(par[k].xi==i && par[k].yj==j) mtr[i][j].mass+=par[k].m;
			}
			for(long k=0;k<n_part;k++){
				if(par[k].xi==i && par[k].yj==j){
					mtr[i][j].cmx+=(par[k].m*par[k].x)/mtr[i][j].mass; //centro de massa em x, para uma dada celula
					mtr[i][j].cmy+=(par[k].m*par[k].y)/mtr[i][j].mass; //centro de massa em y, para uma dada celula
				}
			}
		}
	}
}

void globalcenterofmass (long n_part){
	double masssum, xcm, ycm;
	for(long k=0; k<n_part; k++) masssum += par[k].m;
	for(long k=0; k<n_part; k++){
		xcm+=(par[k].m*par[k].x)/masssum;
		ycm+=(par[k].m*par[k].y)/masssum;
	}
	printf("%f %f\n", xcm, ycm);
}

double accel (long i, long j, long k, int c){//utilizar na func avgforce
	double accel;
	if(!c) accel = G*mtr[i][j].mass/pow((mtr[i][j].cmx-par[k].x),2);
	else accel = G*mtr[i][j].mass/pow((mtr[i][j].cmy-par[k].y),2);
	return accel;
}

double avgaccel(long i, long j, long k, int c){//utilizar na func accel
	double avgaccel;
	avgaccel = (accel(i,j,k,c)+accel(i+1,j,k,c)+accel(i-1,j,k,c)+accel(i,j+1,k,c)+accel(i,j-1,k,c)+accel(i+1,j+1,k,c)+accel(i-1,j-1,k,c),accel(i+1,j-1,k,c)+accel(i-1,j+1,k,c))/8;
}

void velocidade (long k, long tstep){//utilizar na func movement
	par[k].vx+=avgaccel(par[k].xi, par[k].yj, k, 0)*tstep;
	par[k].vy+=avgaccel(par[k].xi, par[k].yj, k, 1)*tstep;
}

void movement (long k, long tstep){
	par[k].x+= par[k].vx*tstep + avgaccel(par[k].xi, par[k].yj,k,0)*tstep;
	par[k].y+= par[k].vy*tstep + avgaccel(par[k].xi, par[k].yj,k,1)*tstep;
}

void updater(long ncside, long n_part){
	long tstep=1;
	for(long k=0;k<n_part;k++){
		velocidade(k, tstep);
		movement(k, tstep);
	}
}

void loop(long ncside, long n_part, long particle_t){
	centerofmass(ncside, n_part);
	for(long k=0; k<particle_t; k++){
		updater(ncside,n_part);
		centerofmass(ncside, n_part);
	}
	printf("%f %f\n", par[0].x, par[0].y);
	globalcenterofmass(n_part);
}

void main(int argc, char** argv){
	const long seed = atoi(argv[1]);
	const long ncside = atoi(argv[2]);
	const long n_part = atoi(argv[3]);
	const long particle_t = atoi(argv[4]);

	init_particles(seed, ncside, n_part, particle_t);
	loop(ncside, n_part, particle_t);
	printf("Finished!\n");

}

