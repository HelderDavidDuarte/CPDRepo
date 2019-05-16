#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <stddef.h>

#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.0005
#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)

#define TAG_MTR 1
#define TAG_PAR 2
#define TAG_MAS 3


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

}MATRIX;

particle_t *par;
MATRIX *mtr;

int rank, p;
MPI_Datatype MPI_particle_t;
	MPI_Datatype MPI_matrix_cell;
	MPI_Datatype matrix_1D;
	MPI_Datatype matrix_2D;


/******************************************************************
void init_particles 
long seed:			seed for random generator (input)
long long n_part:	number of particles (input)
particle_t *par:	pointer to particle vector

PURPOSE : 	Initializes particle position, velocity and mass from 
			random value
RETURN :  	void
********************************************************************/
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

/******************************************************************
double centerofmassinit
long ncside:		    size of grid(input)
long long n_part:	    number of particles (input)
double masssum:		    sum of mass of all particles
MATRIX *mat:            pointer to array which represents grid
particle_t * particle:  pointer to array of particles
int par_size:           number of particles handled in this process

PURPOSE : 	Calculates first iteration of grid cell center of mass, sum 
			of mass of all particles and which cell particles are in
RETURN :  	masssum
********************************************************************/
double centerofmassinit (long ncside, long long n_part, double masssum, MATRIX *mat, particle_t *part, int par_size){
	int ix, jy, idx;
	for(long long k=0;k<par_size;k++){ 
        ix=part[k].x*ncside;
		jy=part[k].y*ncside;
		mat[ix*ncside+jy].mass+=part[k].m; 
		masssum+=part[k].m;
	}

	for(long long k=0; k<par_size; k++){
        ix=part[k].x*ncside;
		jy=part[k].y*ncside;
        idx=ix*ncside+jy;
		mat[idx].cmx+=(part[k].m*part[k].x)/mat[idx].mass; //centre of mass for x, for each grid cell
		mat[idx].cmy+=(part[k].m*part[k].y)/mat[idx].mass; //centre of mass for y, for each grid cell
	}

	for(int k=0;k<ncside*ncside;k++){
			MPI_Allreduce(&(mat[k].mass), &(mat[k].mass), 1, MPI_DOUBLE,
				MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&(mat[k].cmx), &(mat[k].cmx), 1, MPI_DOUBLE,
				MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&(mat[k].cmy), &(mat[k].cmy), 1, MPI_DOUBLE,
				MPI_SUM, MPI_COMM_WORLD);
	}

	MPI_Allreduce(&(masssum), &(masssum), 1, MPI_DOUBLE,
				MPI_SUM, MPI_COMM_WORLD);
	
	return masssum;
}

/******************************************************************
void wrapcalc
long ncside:		    size of grid(input)
long long n_part:	    number of particles (input)
long particle_iter:	    number of iterations to run (input)
double masssum:		    sum of mass of all particles
particle_t * particle:  pointer to array of particles
MATRIX *mat:            pointer to array which represents grid

PURPOSE : 	Bulk of calculations for particle position and velocity,
			as well as grid cell center of mass
RETURN :  	void
********************************************************************/
void wrapcalc(long ncside, long long n_part, long particle_iter, double masssum, particle_t * particle, MATRIX *mat){
	int wwx, wwy;
	int t, u, par_size;
	long long k;
	double xcm=0, ycm=0;
	long m[]={0,0,0}, n[]={0,0,0};
	double local_vx=0, local_vy=0, local_x=0, local_y=0;
	double rx, ry, compvx, compvy;
	int ix, jy, idx;

	if(n_part%p==0){
		par_size=n_part/p;
	}
	else{
		if(rank<n_part%p)
			par_size=(int)(n_part/p)+1;
		else
			par_size=(int)n_part/p;
	}

	masssum=centerofmassinit (ncside, n_part,masssum, mat, particle, par_size);

	for(long l=0; l<particle_iter; l++){

		for(k=0; k<par_size; k++){
			compvx=0, compvy=0;
			wwx=0, wwy=0;
			m[1]=particle[k].x*ncside;
	       	n[1]=particle[k].y*ncside;
			m[2]=m[1]+1,m[0]=m[1]-1,n[2]=n[1]+1,n[0]=n[1]-1;
			if(m[2]>=ncside) {m[2]=0; wwx=1;}
			else if(m[0]<0) {m[0]=ncside-1; wwx=1;}
			if(n[2]>=ncside) {n[2]=0; wwy=1;}
			else if(n[0]<0) {n[0]=ncside-1; wwy=1;}
			for(t=0; t<3; t++){
				for(u=0; u<3; u++){
					rx=mat[m[t]*ncside+n[u]].cmx;
					if(wwx) rx=(-rx);
					if(rx<0 && (-rx)>EPSLON) compvx-=G*mat[m[t]*ncside+n[u]].mass/(rx*rx*9);
					else if(rx>EPSLON) compvx+=G*mat[m[t]*ncside+n[u]].mass/(rx*rx*9);
					ry=mat[m[t]*ncside+n[u]].cmy;
					if(wwy) ry=(-ry);
					if(ry<0 && (-ry)>EPSLON) compvy-=G*mat[m[t]*ncside+n[u]].mass/(ry*ry*9);
					else if(ry>EPSLON) compvy+=G*mat[m[t]*ncside+n[u]].mass/(ry*ry*9);
				}
			}
			particle[k].vx+= compvx;
			particle[k].x+= particle[k].vx + compvx*0.5;
			while(particle[k].x>=1) particle[k].x-=1;
			while(particle[k].x<0) particle[k].x+=1;
			particle[k].vy+= compvy;
			particle[k].y+= particle[k].vy + compvy*0.5;
			while(particle[k].y>=1) particle[k].y-=1;
			while(particle[k].y<0) particle[k].y+=1;	
		}
		
		for(k=0; k<ncside*ncside; k++){
			mat[k].cmx=0;mat[k].cmy=0; mat[k].mass=0;
		}
		
		for(k=0; k<par_size; k++){
            ix=particle[k].x*ncside;
            jy=particle[k].y*ncside;
            idx=ix*ncside+jy;
			mat[idx].mass+=particle[k].m;
			mat[idx].cmx+=(particle[k].m*particle[k].x)/mat[idx].mass; //centre of mass for x, for a given cell
			mat[idx].cmy+=(particle[k].m*particle[k].y)/mat[idx].mass; //centre of mass for y, for a given cell
			if(l==particle_iter-1){//on the last iteration, calculates centre of mass of all particles
				xcm+=(particle[k].m*particle[k].x)/masssum;
				ycm+=(particle[k].m*particle[k].y)/masssum;	
			}
		}
		for(k=0;k<ncside*ncside;k++){
			MPI_Allreduce(&(mat[k].cmx), &(mat[k].cmx), 1, MPI_DOUBLE,
				MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&(mat[k].cmy), &(mat[k].cmy), 1, MPI_DOUBLE,
				MPI_SUM, MPI_COMM_WORLD);
		}
	}
	
	MPI_Allreduce(&xcm, &xcm, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&ycm, &ycm, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
	if(rank==0){
		printf("%.2f %.2f\n", particle[0].x, particle[0].y);
		printf("%.2f %.2f\n", xcm, ycm);
	}
}

/******************************************************************
void usage

PURPOSE : 	To be used in case of input error, prints proper usage
			and exits
RETURN :  	void
********************************************************************/
void usage(){
	printf("Usage: simpar <random generator seed> <grid size> <particle no> <time-step no>\n");
	exit(0);
}

/******************************************************************
particle_t * alloc_particles
long long n_part:	number of particles (input)

PURPOSE : 	Allocate array of particles according to the number of
            tasks (particles) to be handled in each process
RETURN :  	part (pointer to array of particles)
********************************************************************/
particle_t * alloc_particles(long long n_part){

	particle_t * part;
	if(n_part%p==0){
		part=(particle_t*)calloc((n_part/p),sizeof(particle_t));
	}
	else{
		if(rank<n_part%p)
			part=(particle_t*)calloc((int)(n_part/p)+1,sizeof(particle_t));
		else
			part=(particle_t*)calloc((int)(n_part/p),sizeof(particle_t));
	}
	return part;
}

/******************************************************************
void SendParticles
long long n_part:	number of particles (input)

PURPOSE : 	Distributes (sends) particles among all running processes
RETURN :  	void
********************************************************************/
void SendParticles(long long n_part){

	particle_t * prt;
	for (int i=0;i<n_part;i++){
		prt=&par[i];
		MPI_Send(prt, 5, MPI_DOUBLE, i%p, TAG_PAR, MPI_COMM_WORLD);
	}
}

/******************************************************************
particle_t * ReceiveParticle
long long n_part:	number of particles (input)
long ncside:		size of grid(input)

PURPOSE : 	Allocates and receives particles from process 0, which
            initializes them
RETURN :  	part(pointer to particle array)
********************************************************************/
particle_t * ReceiveParticle(long long n_part, long ncside){
	particle_t *part=alloc_particles(n_part);
	int flag = 1, i=0;
	while(1){
		MPI_Iprobe(MPI_ANY_SOURCE, TAG_PAR, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
		if (!flag)
			break;
		MPI_Recv(&part[i], 5, MPI_DOUBLE, MPI_ANY_SOURCE, TAG_PAR, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		i++;
	}
	return part;
}

/******************************************************************
void main
int arcg
char** argv

PURPOSE : 	Gets user inputs, allocates memory for matrix and particles,
			and runs calculation functions
RETURN :  	void
********************************************************************/
void main(int argc, char** argv){

	double secs;
	MATRIX *mat;
	particle_t *part;
    
	MPI_Init (&argc, &argv);
 	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
 	MPI_Comm_size (MPI_COMM_WORLD, &p);

	secs = - MPI_Wtime();
	if(argc!=5) usage();
		double masssum=0;
		long l;
		char *ptr1, *ptr2, *ptr3, *ptr4;
		const long seed = strtol(argv[1], &ptr1, 10);
		long ncside = strtol(argv[2], &ptr2, 10);
		long long n_part = strtol(argv[3], &ptr3, 10);
		const long particle_iter = strtol(argv[4], &ptr4, 10);
		if (*ptr1!=0 || *ptr2!=0 || *ptr3!=0 || *ptr4!=0 || seed <=0 || ncside<=0 || n_part<=0 || particle_iter<=0) usage();
	
	if(rank==0){
		
		if ((par = (particle_t*)calloc(n_part,sizeof(particle_t)))==NULL) exit (0);

		init_particles(seed, ncside, n_part, par);
		SendParticles(n_part);

	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	mat = (MATRIX*)calloc(ncside*ncside,sizeof(MATRIX));

	part=ReceiveParticle(n_part, ncside);

	MPI_Barrier(MPI_COMM_WORLD);
	wrapcalc(ncside, n_part, particle_iter, masssum, part, mat);

	MPI_Barrier(MPI_COMM_WORLD);
    secs += MPI_Wtime();
	if(!rank) printf("Time:%12.6f sec,\n", secs);
    
	MPI_Finalize();

	exit(0);

}

