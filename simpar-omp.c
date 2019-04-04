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



//#pragma omp threadprivate(n_part)
//#pragma omp threadprivate(par)
//#pragma omp threadprivate(mtr)
//#pragma omp threadprivate(masssum)

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

/*void accelx (long i, long j, long long k, int flag){//calculo da aceleracao de uma particula a um dado centro de massa, em x
	double rx=mtr[i][j].cmx;
	if(flag) rx=(-rx);
	if(rx<0 && (-rx)>EPSLON) compvx-=G*mtr[i][j].mass/(rx*rx*9);
	else if(rx>EPSLON) compvx+=G*mtr[i][j].mass/(rx*rx*9);
}

void accely (long i, long j, long long k, int flag){//calculo da aceleracao de uma particula a um dado centro de massa, em y
	double ry=mtr[i][j].cmy;
	if(flag) ry=(-ry);
	if(ry<0 && (-ry)>EPSLON) compvy-=G*mtr[i][j].mass/(ry*ry*9);
	else if(ry>EPSLON) compvy+=G*mtr[i][j].mass/(ry*ry*9);
}*/

void centerofmassinit (long ncside, long long n_part){//calcula a primeira iteracao dos centros de massa, necessaria aos calculos seguintes
	//#pragma omp parallel for reduction(+:masssum)
	for(long long k=0;k<n_part;k++){
		par[k].ix=par[k].x*ncside;
		par[k].jy=par[k].y*ncside;
		mtr[par[k].ix][par[k].jy].mass+=par[k].m;
		masssum+=par[k].m;
	}
	#pragma omp parallel for
	for(long long k=0; k<n_part; k++){
		mtr[par[k].ix][par[k].jy].cmx+=(par[k].m*par[k].x)/mtr[par[k].ix][par[k].jy].mass; //centro de massa em x, para uma dada celula
		mtr[par[k].ix][par[k].jy].cmy+=(par[k].m*par[k].y)/mtr[par[k].ix][par[k].jy].mass; //centro de massa em y, para uma dada celula
	}
}

void wrapcalc(long ncside, long long n_part, long particle_iter){
	int wwx, wwy;
	int t, u;
	long long k;
	double xcm=0, ycm=0;
	long m[]={0,0,0}, n[]={0,0,0};
	long nc_side;
	double local_vx=0, local_vy=0, local_x=0, local_y=0;
	double rx, ry;
	double compvx=0, compvy=0;

	for(long l=0; l<particle_iter; l++){
		

				for(long i=0; i<ncside; i++){
				for(long j=0; j<ncside; j++) mtr[i][j].mass=0;
			}

			#pragma omp parallel for schedule(static, 2)
			for(long long k=0;k<n_part;k++){
			par[k].ix=par[k].x*ncside;
			par[k].jy=par[k].y*ncside;
			mtr[par[k].ix][par[k].jy].mass+=par[k].m;
			}

			//#pragma omp parallel
			//{
			//#pragma omp single
			//#pragma omp parallel for private(p,q,r,s,t,u,compvx,compvy, wwx, wwy, ncside) schedule(dynamic,2)
			//#pragma omp for private(t,u,compvx,compvy, wwx, wwy, m, n) schedule(dynamic, 2) nowait
				/*for(k=0; k<n_part; k++){
					/*local_vx=par[k].vx;
					local_vy=par[k].vy;
					local_x=par[k].x;
					local_y=par[k].y;
					//printf("%lld ", n_part);
					compvx=0, compvy=0;
					wwx=0, wwy=0;
					m[1]=par[k].x*ncside;
		        	n[1]=par[k].y*ncside;
					m[2]=m[1]+1,m[0]=m[1]-1,n[2]=n[1]+1,n[0]=n[1]-1;
					if(m[2]>=ncside) {m[2]=0; wwx=1;}
					else if(m[0]<0) {m[0]=ncside-1; wwx=1;}
					if(n[2]>=ncside) {n[2]=0; wwy=1;}
					else if(n[0]<0) {n[0]=ncside-1; wwy=1;}
				}*/
				#pragma omp parallel
			{
				#pragma omp for private(t,u,compvx, compvy, wwy, wwx, m, n, rx) schedule(static, 2) nowait
				for(k=0; k<n_part; k++){
					compvx=0, compvy=0;
					wwx=0, wwy=0;
					m[1]=par[k].x*ncside;
		        	n[1]=par[k].y*ncside;
					m[2]=m[1]+1,m[0]=m[1]-1,n[2]=n[1]+1,n[0]=n[1]-1;
					if(m[2]>=ncside) {m[2]=0; wwx=1;}
					else if(m[0]<0) {m[0]=ncside-1; wwx=1;}
					if(n[2]>=ncside) {n[2]=0; wwy=1;}
					else if(n[0]<0) {n[0]=ncside-1; wwy=1;}
					for(t=0; t<3; t++){
						for(u=0; u<3; u++){
							rx=mtr[m[t]][n[u]].cmx;
							//printf("%f ", rx);
							if(wwx) rx=(-rx);
							if(rx<0 && (-rx)>EPSLON) compvx-=G*mtr[m[t]][n[u]].mass/(rx*rx*9);
							else if(rx>EPSLON) compvx+=G*mtr[m[t]][n[u]].mass/(rx*rx*9);
						}
					}
					par[k].vx+= compvx;
						//#pragma omp atomic
						par[k].x+= par[k].vx + compvx*0.5;
						//for (local_x; local_x>=1;local_x--){}
						//for (local_x;local_x<0;local_x++){}
						while(par[k].x>=1) par[k].x-=1;
						while(par[k].x<0) par[k].x+=1;
				}
				#pragma omp for private(t,u, compvy, compvx, wwx, wwy, m, n, ry) schedule(static, 2) nowait
				for(k=0; k<n_part; k++){
					compvx=0, compvy=0;
					wwx=0, wwy=0;
					m[1]=par[k].x*ncside;
		        	n[1]=par[k].y*ncside;
					m[2]=m[1]+1,m[0]=m[1]-1,n[2]=n[1]+1,n[0]=n[1]-1;
					if(m[2]>=ncside) {m[2]=0; wwx=1;}
					else if(m[0]<0) {m[0]=ncside-1; wwx=1;}
					if(n[2]>=ncside) {n[2]=0; wwy=1;}
					else if(n[0]<0) {n[0]=ncside-1; wwy=1;}
					for(t=0; t<3; t++){
						for(u=0; u<3; u++){
							ry=mtr[m[t]][n[u]].cmy;
							if(wwy) ry=(-ry);
							if(ry<0 && (-ry)>EPSLON) compvy-=G*mtr[m[t]][n[u]].mass/(ry*ry*9);
							else if(ry>EPSLON) compvy+=G*mtr[m[t]][n[u]].mass/(ry*ry*9);
						}
					}
						//accelx(m[1],n[1],k,wwx);accelx(m[2],n[1],k,wwx);accelx(m[0],n[1],k,wwx);accelx(m[1],n[2],k,wwx);accelx(m[1],n[0],k,wwx);accelx(m[2],n[2],k,wwx);accelx(m[0],n[0],k,wwx);accelx(m[2],n[0],k,wwx);accelx(q,r,k,wwx);
						//update de velocidade e posicao em x
						//#pragma omp atomic
						
						//accely(t,u,k,wwy);accely(p,u,k,wwy);accely(q,u,k,wwy);accely(t,r,k,wwy);accely(t,s,k,wwy);accely(p,r,k,wwy);accely(q,s,k,wwy);accely(p,s,k,wwy);accely(q,r,k,wwy);
						//par[k].vx=local_vx;
						//par[k].x=local_x;
					//update de velocidade e posicao em y
					//#pragma omp atomic
					par[k].vy+= compvy;
					//#pragma omp atomic
					par[k].y+= par[k].vy + compvy*0.5;
					while(par[k].y>=1) par[k].y-=1;
					while(par[k].y<0) par[k].y+=1;
					//par[k].vy=local_vy;
					//par[k].y=local_y;

				}
			

				/*for(k=0;k<n_part;k++){
					par[k].vx=local_vx[k];
					par[k].vy=local_vy[k];
				}*/
				//#pragma omp for
				for(long i=0; i<ncside; i++){
					for(long j=0; j<ncside; j++){mtr[i][j].cmx=0;mtr[i][j].cmy=0;}
				}
		//}
		
		#pragma omp for private(k) reduction(+:xcm) reduction(+:ycm) schedule(dynamic, 4) nowait
		for(k=0; k<n_part; k++){
			par[k].ix=par[k].x*ncside;
			par[k].jy=par[k].y*ncside;
			mtr[par[k].ix][par[k].jy].cmx+=(par[k].m*par[k].x)/mtr[par[k].ix][par[k].jy].mass; //centro de massa em x, para uma dada celula
			mtr[par[k].ix][par[k].jy].cmy+=(par[k].m*par[k].y)/mtr[par[k].ix][par[k].jy].mass; //centro de massa em y, para uma dada celula
			if(l==particle_iter-1){//na ultima iteracao, calcula o centro de massa de todas as particulas
				xcm+=(par[k].m*par[k].x)/masssum;
				ycm+=(par[k].m*par[k].y)/masssum;
			}
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
	long ncside = strtol(argv[2], &ptr2, 10);
	long long n_part = strtol(argv[3], &ptr3, 10);
	const long particle_iter = strtol(argv[4], &ptr4, 10);
	if (*ptr1!=0 || *ptr2!=0 || *ptr3!=0 || *ptr4!=0 || seed <=0 || ncside<=0 || n_part<=0 || particle_iter<=0) usage();

	if ((par = (particle_t*)calloc(n_part,sizeof(particle_t)))==NULL) exit (0);

	if ((mtr = (MATRIX**)calloc(ncside,sizeof(MATRIX*)))==NULL) exit (0);

	long l;
	//#pragma omp parallel for
	for (l=0; l<ncside; l++){
		if ((mtr[l]=(MATRIX*)calloc(ncside,sizeof(MATRIX)))==NULL) exit (0);
	}
	
	init_particles(seed, ncside, n_part, par);
	centerofmassinit(ncside, n_part);
	wrapcalc(ncside, n_part, particle_iter);

}

