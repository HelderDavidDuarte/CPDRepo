#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <time.h> 

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

PARTICLE *par;
MATRIX **mtr;
double masssum=0;

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
        masssum+=par[i].m;
    }
}

void check(long k){
	if(par[k].x>=1) par[k].x-=1;
	if(par[k].x<0)  par[k].x+=1;
	if(par[k].y>=1) par[k].y-=1;
	if(par[k].y<0)  par[k].y+=1;
}

void centerofmass (long ncside, long n_part){
	for(long k=0; k<n_part; k++){
		for(long i=0; i<ncside; i++){
			for(long j=0; j<ncside; j++){
				check(k);
				par[k].xi=floor(par[k].x*ncside);
				par[k].yj=floor(par[k].y*ncside);
				if(par[k].xi==i && par[k].yj==j) mtr[i][j].mass+=par[k].m;
			}
		}
	}
	for(long k=0;k<n_part;k++){
		for(long i=0; i<ncside; i++){
			for(long j=0; j<ncside; j++){
				if(par[k].xi==i && par[k].yj==j){
					mtr[i][j].cmx+=(par[k].m*par[k].x)/mtr[i][j].mass; //centro de massa em x, para uma dada celula
					if(mtr[i][j].cmx>=1) mtr[i][j].cmx-=1;
					if(mtr[i][j].cmx<=1) mtr[i][j].cmx+=1;
					mtr[i][j].cmy+=(par[k].m*par[k].y)/mtr[i][j].mass; //centro de massa em y, para uma dada celula
					if(mtr[i][j].cmy>=1) mtr[i][j].cmy-=1;
					if(mtr[i][j].cmy<=1) mtr[i][j].cmy+=1;
				}
			}
		}
	}
}

double accel (long i, long j, long k, int c){//utilizar na func avgforce
	double accel, rx=mtr[i][j].cmx-par[k].x, ry=mtr[i][j].cmy-par[k].y;
	check(k);
	if(!c) accel = G*mtr[i][j].mass/pow((rx),2);
	else accel = G*mtr[i][j].mass/pow((ry),2);
	if(rx<0.01||ry<0.01) accel=0;
	return accel;
}

double avgaccel(long k, long ncside, int c){//utilizar na func accel
	double avgaccel=0;
	long i=par[k].xi;
	long j=par[k].yj;
	long p=i+1,q=i-1,r=j+1,s=j-1;
	if(p>=ncside) p=0;
	if(q<0) q=ncside-1;
	if(r>=ncside) r=0;
	if(s<0) s=ncside-1;
	avgaccel = (accel(i,j,k,c)+accel(p,j,k,c)+accel(q,j,k,c)+accel(i,r,k,c)+accel(i,s,k,c)+accel(p,r,k,c)+accel(q,s,k,c),accel(p,s,k,c)+accel(q,r,k,c))/9;
	return avgaccel;
}

void velocidade (long k, long tstep, long ncside){//utilizar na func movement
	par[k].vx+=avgaccel(k, ncside, 0)*tstep;
	par[k].vy+=avgaccel(k, ncside, 1)*tstep;
}

void movement (long k, long tstep, long ncside){
	par[k].x+= par[k].vx*tstep + (avgaccel(k,ncside,0)*tstep*tstep)/2;
	par[k].y+= par[k].vy*tstep + (avgaccel(k,ncside,1)*tstep*tstep)/2;
}

void run(long ncside, long n_part, long particle_t){
	long tstep=1;
	centerofmass(ncside, n_part);
	for(long l=0; l<particle_t; l++){
		for(long k=0;k<n_part;k++){
			velocidade(k, tstep, ncside);
			movement(k, tstep, ncside);
		}
		centerofmass(ncside, n_part);
	}
	printf("%.2f %.2f\n", par[0].x, par[0].y);
}

void globalcenterofmass (long n_part){
	double xcm=0, ycm=0;
	for(long k=0; k<n_part; k++){
		xcm+=(par[k].m*par[k].x)/masssum;
		ycm+=(par[k].m*par[k].y)/masssum;
	}
	printf("%.2f %.2f\n", xcm, ycm);
}

void main(int argc, char** argv){
	const long seed = atoi(argv[1]);
	const long ncside = atoi(argv[2]);
	const long n_part = atoi(argv[3]);
	const long particle_t = atoi(argv[4]);
	clock_t start, end;
    double cpu_time_used;
    start = clock();

	par = calloc(n_part,sizeof(PARTICLE));
	mtr = (MATRIX**)calloc(ncside,sizeof(MATRIX*));
	for (int l=0; l<ncside; l++){
		mtr[l]=(MATRIX*)calloc(ncside,sizeof(MATRIX));
	}
	
	init_particles(seed, ncside, n_part, particle_t);
	run(ncside, n_part, particle_t);
	globalcenterofmass(n_part);

	end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("%f\n", cpu_time_used);
}

