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
	double rx=mtr[i][j].cmx-par[k].x;
	if(rx<EPSLON) return 0;
	return G*mtr[i][j].mass/(rx*rx*9);
}
/////////////////////////////////////////////////// as funcoes de calculo das aceleracoes nao aparentam ser fontes de erro
double accely (long i, long j, long long k){//calculo da aceleracao de uma particula a um dado centro de massa, em y
	double ry=mtr[i][j].cmy-par[k].y;
	if(ry<EPSLON) return 0;
	return G*mtr[i][j].mass/(ry*ry*9);
}

void centerofmassinit (long ncside, long long n_part){//calcula a primeira iteracao dos centros de massa, necessaria aos calculos seguintes
	for(long long k=0;k<n_part;k++){
		par[k].ix=par[k].x*ncside;
        par[k].jy=par[k].y*ncside;
		for(long i=0; i<ncside; i++){
			for(long j=0; j<ncside && i==par[k].ix && j==par[k].jy; j++){
				//if(i==par[k].ix && j==par[k].jy)
				mtr[i][j].mass+=par[k].m; ///////////////////////PODE SER FONTE DE ERRO!
			}
		}
	}
	for(long long k=0; k<n_part; k++){
		masssum+=par[k].m;
		for(long i=0; i<ncside; i++){
			for(long j=0; j<ncside && mtr[i][j].mass!=0; j++){  ///////////////////////PODE SER FONTE DE ERRO!
				mtr[i][j].cmx+=(par[k].m*par[k].x)/mtr[i][j].mass; //centro de massa em x, para uma dada celula
				mtr[i][j].cmy+=(par[k].m*par[k].y)/mtr[i][j].mass; //centro de massa em y, para uma dada celula
			}
		}
	}
}

void wrapcalc(long ncside, long long n_part, long particle_iter){
	long i,j,p,q,r,s,m[3],n[3]; //timestep = 1
	double compvx, compvy, xcm=0, ycm=0;
	for(long l=0; l<particle_iter; l++){
		for(long i=0; i<ncside; i++){
			for(long j=0; j<ncside; j++) mtr[i][j].mass=0;
		}
		for(long long k=0;k<n_part;k++){
			for(long i=0; i<ncside; i++){
				for(long j=0; j<ncside && i==par[k].ix && j==par[k].jy; j++){
					//if(i==par[k].ix && j==par[k].jy)
					mtr[i][j].mass+=par[k].m; ///////////////////////PODE SER FONTE DE ERRO!
				}
			}
		}
		for(long long k=0; k<n_part; k++){
			compvx=0, compvy=0;
			m[1]=par[k].ix=par[k].x*ncside; //m e n sao vetores que constroem a grelha auxiliar de 3*3
        	n[1]=par[k].jy=par[k].y*ncside;
			m[2]=m[1]+1,m[0]=m[1]-1,n[2]=n[1]+1,n[0]=n[1]-1;
			if(m[2]>=ncside) m[2]=0;
			if(m[0]<0) m[0]=ncside-1;
			if(n[2]>=ncside) n[2]=0;
			if(n[0]<0) n[0]=ncside-1;
			for(long v=0; v<3;v++){
				for (long w=0; w<3;w++){
					compvx+=accelx(m[v],n[w],k);
					compvy+=accely(m[v],n[w],k);
				}
			}
			//printf("%g\n", compvx);
			//update de velocidade e posicao em x
			par[k].vx+= compvx;
			par[k].x+= par[k].vx + compvx*0.5;
			if(par[k].x>=1) par[k].x-=1;
			else if(par[k].x<0) par[k].x+=1;
			//update de velocidade e posicao em y
			par[k].vy+= compvy;
			par[k].y+= par[k].vy + compvy*0.5;
			if(par[k].y>=1) par[k].y-=1;
			else if(par[k].y<0) par[k].y+=1;
			par[k].ix=par[k].x*ncside;
        	par[k].jy=par[k].y*ncside;
			if(l!=particle_iter-1){
				for(long i=0; i<ncside; i++){
					for(long j=0; j<ncside; j++){mtr[i][j].cmx=0;mtr[i][j].cmy=0;}
				}
				for(long i=0; i<ncside; i++){
					for(long j=0; j<ncside && mtr[i][j].mass!=0; j++){  ///////////////////////PODE SER FONTE DE ERRO!
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

