#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <omp.h> 

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

}particle_t;

typedef struct matrix
{
	double mass;
	double cmx;
	double cmy;
	long ix;
	long jy;

}MATRIX;

particle_t *par;
MATRIX *mtr;
double masssum=0, xcm=0, ycm=0;
int w=0;

void init_particles(long seed, long ncside, long long n_part, particle_t *par)
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

double accel (long i, long j, long k, int c){
	double accel, rx=mtr[i+j].cmx-par[k].x, ry=mtr[i+j].cmy-par[k].y;
	if(!c) accel = G*mtr[i+j].mass/(rx*rx);
	else accel = G*mtr[i+j].mass/(ry*ry);
	if(rx<0.01||ry<0.01) accel=0;
	return accel;
}

double avgaccel(long i, long j, long p, long q, long r, long s, long k, int c){
	return (accel(i,j,k,c)+accel(p,j,k,c)+accel(q,j,k,c)+accel(i,r,k,c)+accel(i,s,k,c)+accel(p,r,k,c)+accel(q,s,k,c),accel(p,s,k,c)+accel(q,r,k,c))/9;
}

void wrapcalc(long ncside, long n_part){
	long tstep=1,i,j,p,q,r,s;
	double ax, ay;
	for(long k=0; k<n_part; k++){
		i=par[k].xi,j=par[k].yj;
		p=i+1,q=i-1,r=j+1,s=j-1;
		if(p>=ncside) p=0;
		if(q<0) q=ncside-1;
		if(r>=ncside) r=0;
		if(s<0) s=ncside-1;
		ax=avgaccel(i,j,p,q,r,s,k,0);
		ay=avgaccel(i,j,p,q,r,s,k,1);
		par[k].vx+= ax*tstep;
		par[k].x+= par[k].vx*tstep + (ax*tstep*tstep)/2;
		par[k].vy+= ay*tstep;
		par[k].y+= par[k].vy*tstep + (ay*tstep*tstep)/2;
	}
}

void centerofmass (long ncside, long n_part){
	for(long k=0; k<n_part; k++){
		if(par[k].x>=1) par[k].x-=1;
		else if(par[k].x<0) par[k].x+=1;
		if(par[k].y>=1) par[k].y-=1;
		else if(par[k].y<0) par[k].y+=1;
		par[k].xi=floor(par[k].x*ncside);
		par[k].yj=floor(par[k].y*ncside);
		for(long i=0; i<ncside*ncside; i++){
			if(par[k].xi==mtr[i].ix && par[k].yj==mtr[i].jy){
				mtr[i].mass+=par[k].m;
				mtr[i].cmx+=(par[k].m*par[k].x)/mtr[i].mass; //centro de massa em x, para uma dada celula
				if(mtr[i].cmx>=1) mtr[i].cmx-=1;
				else if(mtr[i].cmx<=1) mtr[i].cmx+=1;
				mtr[i].cmy+=(par[k].m*par[k].y)/mtr[i].mass; //centro de massa em y, para uma dada celula
				if(mtr[i].cmy>=1) mtr[i].cmy-=1;
				else if(mtr[i].cmy<=1) mtr[i].cmy+=1;
			}
		}
		if(w){
			xcm+=(par[k].m*par[k].x)/masssum;
			ycm+=(par[k].m*par[k].y)/masssum;
		}
	}
}

void run(long ncside, long n_part, long particle_iter){
	centerofmass(ncside, n_part);
	for(long l=0; l<particle_iter; l++){
		if(l==particle_iter-1) w=1;
		wrapcalc(ncside,n_part);
		centerofmass(ncside, n_part);
	}
}

void main(int argc, char** argv){
	const long seed = atoi(argv[1]);
	const long ncside = atoi(argv[2]);
	const long long n_part = atoi(argv[3]);
	const long particle_iter = atoi(argv[4]);
	clock_t start, end;
    double cpu_time_used;
    start = clock();

	par = (particle_t*)calloc(n_part,sizeof(particle_t));
	mtr = (MATRIX*)calloc(ncside*ncside,sizeof(MATRIX));
	
	init_particles(seed, ncside, n_part, par);
	run(ncside, n_part, particle_iter);

	printf("%.2f %.2f\n", par[0].x, par[0].y);
	printf("%.2f %.2f\n", xcm, ycm);

	end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("%f\n", cpu_time_used);
}

