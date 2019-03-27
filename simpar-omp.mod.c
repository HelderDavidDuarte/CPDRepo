#ifdef _POMP
#  undef _POMP
#endif
#define _POMP 200110

#include "simpar-omp.c.opari.inc"
#line 1 "simpar-omp.c"
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
	double cmx;
	double cmy;

}particle_t;

typedef struct matrix
{
	double mass;
	long ix;
	long jy;

}MATRIX;

particle_t *par;
MATRIX *mtr;
double masssum=0;

void init_particles(long seed, long ncside, long long n_part, particle_t *par){
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

void init_matrix(long ncside){//funcao que inicializa o vetor de estruturas, associando valores i,j para a matriz
POMP_Parallel_fork(&omp_rd_1);
#line 53 "simpar-omp.c"
	#pragma omp parallel    
{ POMP_Parallel_begin(&omp_rd_1);
POMP_For_enter(&omp_rd_1);
#line 53 "simpar-omp.c"
 #pragma omp          for nowait
	for(long i=0;i<ncside;i++){
		mtr[i].ix=i;
		for(long j=0;j<ncside;j++) mtr[j].jy=j;
	}
POMP_Barrier_enter(&omp_rd_1);
#pragma omp barrier
POMP_Barrier_exit(&omp_rd_1);
POMP_For_exit(&omp_rd_1);
POMP_Parallel_end(&omp_rd_1); }
POMP_Parallel_join(&omp_rd_1);
#line 58 "simpar-omp.c"
}

double accelx (long t, long long k){//calculo da aceleracao de uma particula a um dado centro de massa, em x
	double rx;
	if((rx=par[k].cmx-par[k].x)<0.0005) return 0;
	return G*mtr[t].mass/(rx*rx);
}

double accely (long t, long long k){//calculo da aceleracao de uma particula a um dado centro de massa, em y
	double ry;
	if((ry=par[k].cmy-par[k].y)<0.0005) return 0;
	return G*mtr[t].mass/(ry*ry);
}

void centerofmassinit (long ncside, long long n_part){//calcula a primeira iteracao dos centros de massa, necessaria aos calculos seguintes
POMP_Parallel_fork(&omp_rd_2);
#line 73 "simpar-omp.c"
	#pragma omp parallel    
{ POMP_Parallel_begin(&omp_rd_2);
POMP_For_enter(&omp_rd_2);
#line 73 "simpar-omp.c"
 #pragma omp          for nowait
	for(long n=0; n<ncside*ncside; n++){
			for(long long k=0; k<n_part && floor(par[k].x*ncside)==mtr[n].ix && floor(par[k].y*ncside)==mtr[n].jy; k++)
				mtr[n].mass+=par[k].m;
		}
POMP_Barrier_enter(&omp_rd_2);
#pragma omp barrier
POMP_Barrier_exit(&omp_rd_2);
POMP_For_exit(&omp_rd_2);
POMP_Parallel_end(&omp_rd_2); }
POMP_Parallel_join(&omp_rd_2);
#line 78 "simpar-omp.c"
	for(long long k=0; k<n_part; k++){
		masssum+=par[k].m;
		for(long n=0; floor(par[k].x*ncside)==mtr[n].ix && floor(par[k].y*ncside)==mtr[n].jy && n<ncside*ncside; n++){
			par[k].cmx+=(par[k].m*par[k].x)/mtr[n].mass; //centro de massa em x, para uma dada celula
			par[k].cmy+=(par[k].m*par[k].y)/mtr[n].mass; //centro de massa em y, para uma dada celula
			break;
		}
	}
POMP_Parallel_fork(&omp_rd_3);
#line 86 "simpar-omp.c"
	#pragma omp parallel    
{ POMP_Parallel_begin(&omp_rd_3);
POMP_For_enter(&omp_rd_3);
#line 86 "simpar-omp.c"
 #pragma omp          for nowait
	for(long n=0; n<ncside*ncside; n++) mtr[n].mass=0;
POMP_Barrier_enter(&omp_rd_3);
#pragma omp barrier
POMP_Barrier_exit(&omp_rd_3);
POMP_For_exit(&omp_rd_3);
POMP_Parallel_end(&omp_rd_3); }
POMP_Parallel_join(&omp_rd_3);
#line 88 "simpar-omp.c"
}

void wrapcalc(long ncside, long long n_part, long particle_iter){
	long i,j,p,q,r,s; //timestep = 1
	double compvx, compvy, xcm=0, ycm=0;
	for(long l=0; l<particle_iter; l++){
POMP_Parallel_fork(&omp_rd_4);
#line 94 "simpar-omp.c"
		#pragma omp parallel    
{ POMP_Parallel_begin(&omp_rd_4);
POMP_For_enter(&omp_rd_4);
#line 94 "simpar-omp.c"
  #pragma omp          for nowait
		for(long n=0; n<ncside*ncside; n++){
			for(long long k=0; k<n_part && floor(par[k].x*ncside)==mtr[n].ix && floor(par[k].y*ncside)==mtr[n].jy; k++)
				mtr[n].mass+=par[k].m;
		}
POMP_Barrier_enter(&omp_rd_4);
#pragma omp barrier
POMP_Barrier_exit(&omp_rd_4);
POMP_For_exit(&omp_rd_4);
POMP_Parallel_end(&omp_rd_4); }
POMP_Parallel_join(&omp_rd_4);
#line 99 "simpar-omp.c"
		for(long long k=0; k<n_part; k++){
			i=par[k].x*ncside,j=par[k].y*ncside;
			p=i+1,q=i-1,r=j+1,s=j-1;
			if(p>=ncside) p=0;
			else if(q<0) q=ncside-1;
			if(r>=ncside) r=0;
			else if(s<0) s=ncside-1;
			//update de velocidade e posicao em x
			compvx=(accelx(i+j,k)+accelx(p+j,k)+accelx(q+j,k)+accelx(i+r,k)+accelx(i+s,k)+accelx(p+r,k)+accelx(q+s,k),accelx(p+s,k)+accelx(q+r,k))/9;
			par[k].vx+= compvx;
			par[k].x+= par[k].vx + compvx*0.5;
			if(par[k].x>=1) par[k].x-=1;
			else if(par[k].x<0) par[k].x+=1;
			//update de velocidade e posicao em y
			compvy=(accely(i+j,k)+accely(p+j,k)+accely(q+j,k)+accely(i+r,k)+accely(i+s,k)+accely(p+r,k)+accely(q+s,k),accely(p+s,k)+accely(q+r,k))/9;
			par[k].vy+= compvy;
			par[k].y+= par[k].vy + compvy*0.5;
			if(par[k].y>=1) par[k].y-=1;
			else if(par[k].y<0) par[k].y+=1;
			
			if(l==particle_iter-1){//na ultima iteracao, calcula o centro de massa de todas as particulas
				xcm+=(par[k].m*par[k].x)/masssum;
				ycm+=(par[k].m*par[k].y)/masssum;
			}
			else{
				for(long n=0; floor(par[k].x*ncside)==mtr[n].ix && floor(par[k].y*ncside)==mtr[n].jy && n<ncside*ncside; n++){//update do centro de massa
					par[k].cmx+=(par[k].m*par[k].x)/mtr[n].mass; //centro de massa em x, para uma dada celula
					par[k].cmy+=(par[k].m*par[k].y)/mtr[n].mass; //centro de massa em y, para uma dada celula
					break;
				}
			}
		}
POMP_Parallel_fork(&omp_rd_5);
#line 131 "simpar-omp.c"
		#pragma omp parallel    
{ POMP_Parallel_begin(&omp_rd_5);
POMP_For_enter(&omp_rd_5);
#line 131 "simpar-omp.c"
  #pragma omp          for nowait
		for(long n=0; n<ncside*ncside; n++) mtr[n].mass=0;
POMP_Barrier_enter(&omp_rd_5);
#pragma omp barrier
POMP_Barrier_exit(&omp_rd_5);
POMP_For_exit(&omp_rd_5);
POMP_Parallel_end(&omp_rd_5); }
POMP_Parallel_join(&omp_rd_5);
#line 133 "simpar-omp.c"
	}
	printf("%.2f %.2f\n", par[0].x, par[0].y);
	printf("%.2f %.2f\n", xcm, ycm);
}

void usage(){
	printf("Usage: simpar <random generator seed> <grid size> <particle no> <time-step no>\n");
	exit(0);
}

void main(int argc, char** argv){

	if(argc!=5 || argv[1]<=0 || argv[2]<=0 || argv[3]<=0 || argv[4]<=0) usage();

	char *ptr1, *ptr2, *ptr3, *ptr4;
	const long seed = strtol(argv[1], &ptr1, 10);
	const long ncside = strtol(argv[2], &ptr2, 10);
	const long long n_part = strtol(argv[3], &ptr3, 10);
	const long particle_iter = strtol(argv[4], &ptr4, 10);
	if (*ptr1!=0 || *ptr2!=0 || *ptr3!=0 || *ptr4!=0) usage();

	clock_t start, end;
    double cpu_time_used;
    start = clock();

	if ((par = (particle_t*)calloc(n_part,sizeof(particle_t)))==NULL) exit (0);

	if ((mtr = (MATRIX*)calloc(ncside*ncside,sizeof(MATRIX)))==NULL) exit (0);
	
	init_particles(seed, ncside, n_part, par);
	init_matrix(ncside);
	centerofmassinit(ncside, n_part);
	wrapcalc(ncside, n_part, particle_iter);

	end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("%f\n", cpu_time_used);
}

