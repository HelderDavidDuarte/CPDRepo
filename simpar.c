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
#define EPSLON 0.0005

typedef struct particle
{
	double x;
	double y;
	double vx;
	double vy;
	double m;
	long ix;
	long jy;

}particle_t;

typedef struct matrix
{
	double mass;
	double cmx;
	double cmy;

}MATRIX;

particle_t *par;
MATRIX **mtr;
double masssum=0;
double compvx=0, compvy=0, wwx, wwy;

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

double accelx (long i, long j, long long k){//calculo da aceleracao de uma particula a um dado centro de massa, em x
	double tempx=mtr[i][j].cmx, rx;
	if(wwx) tempx=1-tempx;
	rx=tempx-par[k].x;
	if(rx<0 && (-rx)>EPSLON) compvx-=G*mtr[i][j].mass/(rx*rx*9);
	else if(rx>EPSLON) compvx+=G*mtr[i][j].mass/(rx*rx*9);
}

double accely (long i, long j, long long k){//calculo da aceleracao de uma particula a um dado centro de massa, em y
	double tempy=mtr[i][j].cmy, ry;
	if(wwy) tempy=1-tempy;
	ry=tempy-par[k].y;
	if(ry<0 && (-ry)>EPSLON) compvy-=G*mtr[i][j].mass/(ry*ry*9);
	else if(ry>EPSLON) compvy+=G*mtr[i][j].mass/(ry*ry*9);
}

void centerofmassinit (long ncside, long long n_part){//calcula a primeira iteracao dos centros de massa, necessaria aos calculos seguintes
	for(long long k=0;k<n_part;k++){
		par[k].ix=par[k].x*ncside;
        par[k].jy=par[k].y*ncside;
		for(long i=0; i<ncside; i++){
			for(long j=0; j<ncside; j++){
				mtr[i][j].mass+=par[k].m;
			}
		}
	}
	for(long long k=0; k<n_part; k++){
		masssum+=par[k].m;
		for(long i=0; i<ncside; i++){
			for(long j=0; j<ncside; j++){
				mtr[i][j].cmx+=(par[k].m*par[k].x)/mtr[i][j].mass; //centro de massa em x, para uma dada celula
				mtr[i][j].cmy+=(par[k].m*par[k].y)/mtr[i][j].mass; //centro de massa em y, para uma dada celula
			}
		}
	}
}

void wrapcalc(long ncside, long long n_part, long particle_iter){
	long i,j,p,q,r,s,t,u; //timestep = 1
	double xcm=0, ycm=0;
	for(long l=0; l<particle_iter; l++){
		for(long i=0; i<ncside; i++){
			for(long j=0; j<ncside; j++) mtr[i][j].mass=0;
		}
		for(long long k=0;k<n_part;k++){
			for(long i=0; i<ncside; i++){
				for(long j=0; j<ncside; j++){
					mtr[i][j].mass+=par[k].m;
				}
			}
		}
		for(long long k=0; k<n_part; k++){
			compvx=0, compvy=0;
			t=par[k].ix=par[k].x*ncside;
        	u=par[k].jy=par[k].y*ncside;
			p=t+1,q=t-1,r=u+1,s=u-1;
			if(p>=ncside) {p=0; wwx=1;}
			if(q<0) {q=ncside-1; wwx=1;}
			if(r>=ncside) {r=0; wwy=1;}
			if(s<0) {s=ncside-1; wwy=1;}
			accelx(t,u,k);accelx(p,u,k);accelx(q,u,k);accelx(t,r,k);accelx(t,s,k);accelx(p,r,k);accelx(q,s,k);accelx(p,s,k);accelx(q,r,k);
			accely(t,u,k);accely(p,u,k);accely(q,u,k);accely(t,r,k);accely(t,s,k);accely(p,r,k);accely(q,s,k);accely(p,s,k);accely(q,r,k);
			
			//update de velocidade e posicao em x
			par[k].vx+= compvx;
			par[k].x+= par[k].vx + compvx*0.5;
			while(par[k].x>=1) par[k].x-=1;
			while(par[k].x<0) par[k].x+=1;
			//update de velocidade e posicao em y
			par[k].vy+= compvy;
			par[k].y+= par[k].vy + compvy*0.5;
			while(par[k].y>=1) par[k].y-=1;
			while(par[k].y<0) par[k].y+=1;
			par[k].ix=par[k].x*ncside;
        	par[k].jy=par[k].y*ncside;
			if(l!=particle_iter-1){
				for(long i=0; i<ncside; i++){
					for(long j=0; j<ncside; j++){mtr[i][j].cmx=0;mtr[i][j].cmy=0;}
				}
				for(long i=0; i<ncside; i++){
					for(long j=0; j<ncside; j++){
						mtr[i][j].cmx+=(par[k].m*par[k].x)/mtr[i][j].mass; //centro de massa em x, para uma dada celula
						mtr[i][j].cmy+=(par[k].m*par[k].y)/mtr[i][j].mass; //centro de massa em y, para uma dada celula
					}
				}
			}
			else{//na ultima iteracao, calcula o centro de massa de todas as particulas
				xcm+=(par[k].m*par[k].x)/masssum;
				ycm+=(par[k].m*par[k].y)/masssum;
			}
		}
	}
	printf("%.2f %.2f\n", par[0].x, par[0].y);
	printf("%.2f %.2f\n", xcm, ycm);
}

void usage(){
	printf("Usage: simpar <random generator seed> <grid size> <particle no> <time-step no>\n");
	exit(0);
}

void main(int argc, char** argv){

	if(argc!=5) usage();

	char *ptr1, *ptr2, *ptr3, *ptr4;
	const long seed = strtol(argv[1], &ptr1, 10);
	const long ncside = strtol(argv[2], &ptr2, 10);
	const long long n_part = strtol(argv[3], &ptr3, 10);
	const long particle_iter = strtol(argv[4], &ptr4, 10);
	if (*ptr1!=0 || *ptr2!=0 || *ptr3!=0 || *ptr4!=0 || seed <=0 || ncside<=0 || n_part<=0 || particle_iter<=0) usage();

	clock_t start, end;
    double cpu_time_used;
    start = clock();

	if ((par = (particle_t*)calloc(n_part,sizeof(particle_t)))==NULL) exit (0);

	if ((mtr = (MATRIX**)calloc(ncside,sizeof(MATRIX*)))==NULL) exit (0);

	for (long l=0; l<ncside; l++){
		if ((mtr[l]=(MATRIX*)calloc(ncside,sizeof(MATRIX)))==NULL) exit (0);
	}
	
	init_particles(seed, ncside, n_part, par);
	centerofmassinit(ncside, n_part);
	wrapcalc(ncside, n_part, particle_iter);

	end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("%f\n", cpu_time_used);
}

