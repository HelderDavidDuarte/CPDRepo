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
double masssum=0;

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

void init_matrix(long ncside){
	for(long i=0;i<ncside;i++){
		mtr[i].ix=i;
		for(long j=0;j<ncside;j++){
			mtr[j].jy=j;
		}
	}
}

double accelx (long t, long k){
	double rx=mtr[t].cmx-par[k].x;
	if(rx<0.01) return 0;
	return G*mtr[t].mass/(rx*rx);
}

double accely (long t, long k){
	double ry=mtr[t].cmy-par[k].y;
	if(ry<0.01) return 0;
	return G*mtr[t].mass/(ry*ry);
}

void centerofmassinit (long ncside, long n_part){
	for(long k=0; k<n_part; k++){
		for(long i=0; i<ncside*ncside && floor(par[k].x*ncside)==mtr[i].ix && floor(par[k].y*ncside)==mtr[i].jy; i++){
			mtr[i].mass+=par[k].m;
			mtr[i].cmx+=(par[k].m*par[k].x)/mtr[i].mass; //centro de massa em x, para uma dada celula
			mtr[i].cmy+=(par[k].m*par[k].y)/mtr[i].mass; //centro de massa em y, para uma dada celula
		}
	}
}

void wrapcalc(long ncside, long n_part){
	long tstep=1,i,j,p,q,r,s;
	double compvx, compvy;
	for(long k=0; k<n_part; k++){
		i=floor(par[k].x*ncside),j=floor(par[k].y*ncside);
		p=i+1,q=i-1,r=j+1,s=j-1;
		if(p>=ncside) p=0;
		if(q<0) q=ncside-1;
		if(r>=ncside) r=0;
		if(s<0) s=ncside-1;
		compvx=((accelx(i+j,k)+accelx(p+j,k)+accelx(q+j,k)+accelx(i+r,k)+accelx(i+s,k)+accelx(p+r,k)+accelx(q+s,k),accelx(p+s,k)+accelx(q+r,k))/9)*tstep;
		compvy=((accely(i+j,k)+accely(p+j,k)+accely(q+j,k)+accely(i+r,k)+accely(i+s,k)+accely(p+r,k)+accely(q+s,k),accely(p+s,k)+accely(q+r,k))/9)*tstep;
		par[k].vx+= compvx;
		par[k].x+= par[k].vx*tstep + (compvx*tstep)/2;
		par[k].vy+= compvy;
		par[k].y+= par[k].vy*tstep + (compvy*tstep)/2;
		if(par[k].x>=1) par[k].x-=1;
		else if(par[k].x<0) par[k].x+=1;
		if(par[k].y>=1) par[k].y-=1;
		else if(par[k].y<0) par[k].y+=1;
		for(long i=0; i<ncside*ncside && floor(par[k].x*ncside)==mtr[i].ix && floor(par[k].y*ncside)==mtr[i].jy; i++){
			mtr[i].mass+=par[k].m;
			mtr[i].cmx+=(par[k].m*par[k].x)/mtr[i].mass; //centro de massa em x, para uma dada celula
			mtr[i].cmy+=(par[k].m*par[k].y)/mtr[i].mass; //centro de massa em y, para uma dada celula
		}
	}
}

void centerofmassfinal (long ncside, long n_part){
	double xcm=0, ycm=0;
	for(long k=0; k<n_part; k++){
		xcm+=(par[k].m*par[k].x)/masssum;
		ycm+=(par[k].m*par[k].y)/masssum;
	}
	printf("%.2f %.2f\n", par[0].x, par[0].y);
	printf("%.2f %.2f\n", xcm, ycm);
}

void run(long ncside, long n_part, long particle_iter){
	centerofmassinit(ncside, n_part);
	for(long l=0; l<particle_iter; l++){
		wrapcalc(ncside,n_part);
	}
	centerofmassfinal(ncside, n_part);
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
	init_matrix(ncside);
	run(ncside, n_part, particle_iter);

	end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("%f\n", cpu_time_used);
}

