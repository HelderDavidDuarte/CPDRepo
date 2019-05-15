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

#define ECHO_ON 1

#define TAG_MTR 1
#define TAG_PAR 2


typedef struct particle
{
	double x;
	double y;
	double vx;
	double vy;
	double m;
	double ix;
	double jy;

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


void Create_MPI_particle_t(){

	int count = 7;
	int array_of_blocklengths[] = { 1, 1, 1, 1, 1, 1, 1 };
	MPI_Aint array_of_displacements[] = { offsetof(particle_t,x), offsetof(particle_t,y), offsetof(particle_t,vx), 
										offsetof(particle_t,vy), offsetof(particle_t,m), offsetof(particle_t,ix), offsetof(particle_t,jy)};
	MPI_Datatype array_of_types[] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_LONG, MPI_LONG};
	MPI_Datatype MPI_particle_t;

	MPI_Type_create_struct( count, array_of_blocklengths, array_of_displacements, array_of_types, &MPI_particle_t );
	MPI_Type_commit( &MPI_particle_t );

}


void Create_MPI_matrix(long ncside){

	int count = 3;
	int array_of_blocklengths[] = { 1, 1, 1};
	MPI_Aint array_of_displacements[] = { offsetof(MATRIX,mass), offsetof(MATRIX, cmx), offsetof(MATRIX, cmy)};
	MPI_Datatype array_of_types[] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
	MPI_Datatype MPI_matrix_cell;

	MPI_Type_create_struct( count, array_of_blocklengths, array_of_displacements, array_of_types, &MPI_matrix_cell );
	MPI_Type_commit( &MPI_matrix_cell);

	MPI_Type_contiguous(ncside, MPI_matrix_cell, &matrix_1D);
    MPI_Type_commit(&matrix_1D);

    MPI_Type_contiguous(ncside, matrix_1D, &matrix_2D);
    MPI_Type_commit(&matrix_2D);

}

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
long ncside:		size of grid(input)
long long n_part:	number of particles (input)
double masssum:		sum of mass of all particles

PURPOSE : 	Calculates first iteration of grid cell center of mass, sum 
			of mass of all particles and which cell particles are in
RETURN :  	masssum
********************************************************************/
double centerofmassinit (long ncside, long long n_part, double masssum){
	int ix=0, jy=0;
	for(long long k=0;k<n_part;k++){
		par[k].ix=par[k].x*ncside;
		par[k].jy=par[k].y*ncside;
		ix=(int)par[k].ix;
		jy=(int)par[k].jy;
		mtr[ix*ncside+jy].mass+=par[k].m;
		masssum+=par[k].m;
	}

	for(long long k=0; k<n_part; k++){
		mtr[ix*ncside+jy].cmx+=(par[k].m*par[k].x)/mtr[ix*ncside+jy].mass; //centre of mass for x, for each grid cell
		mtr[ix*ncside+jy].cmy+=(par[k].m*par[k].y)/mtr[ix*ncside+jy].mass; //centre of mass for y, for each grid cell
	}
	return masssum;
}

/******************************************************************
void wrapcalc
long ncside:		size of grid(input)
long long n_part:	number of particles (input)
long particle_iter:	number of iterations to run (input)
double masssum:		sum of mass of all particles

PURPOSE : 	Bulk of calculations for particle position and velocity,
			as well as grid cell center of mass
RETURN :  	void
********************************************************************/
/*void wrapcalc(long ncside, long long n_part, long particle_iter, double masssum){
	int wwx, wwy;
	int t, u;
	long long k;
	double xcm=0, ycm=0;
	long m[]={0,0,0}, n[]={0,0,0};
	long j;
	double local_vx=0, local_vy=0, local_x=0, local_y=0;
	double rx, ry;
	double compvx=0, compvy=0;

	for(long l=0; l<particle_iter; l++){
		
		for(long i=0; i<ncside; i++){	//set mass in each cell to 0, to be calculated for this iteration
			for(long j=0; j<ncside; j++) mtr[i][j].mass=0;
			}

		for(long long k=0;k<n_part;k++){
			par[k].ix=par[k].x*ncside;
			par[k].jy=par[k].y*ncside;
			mtr[par[k].ix][par[k].jy].mass+=par[k].m;
		}


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
						if(wwx) rx=(-rx);
						if(rx<0 && (-rx)>EPSLON) compvx-=G*mtr[m[t]][n[u]].mass/(rx*rx*9);
						else if(rx>EPSLON) compvx+=G*mtr[m[t]][n[u]].mass/(rx*rx*9);
					}
				}
				par[k].vx+= compvx;
				par[k].x+= par[k].vx + compvx*0.5;
				while(par[k].x>=1) par[k].x-=1;
				while(par[k].x<0) par[k].x+=1;
			}


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
				par[k].vy+= compvy;
				par[k].y+= par[k].vy + compvy*0.5;
				while(par[k].y>=1) par[k].y-=1;
				while(par[k].y<0) par[k].y+=1;	
			}
		
		
		for(long i=0; i<ncside; i++){
			for(j=0; j<ncside; j++){mtr[i][j].cmx=0;mtr[i][j].cmy=0;}
		}
		

		for(k=0; k<n_part; k++){
			par[k].ix=par[k].x*ncside;
			par[k].jy=par[k].y*ncside;
			mtr[par[k].ix][par[k].jy].cmx+=(par[k].m*par[k].x)/mtr[par[k].ix][par[k].jy].mass; //centre of mass for x, for a given cell
			mtr[par[k].ix][par[k].jy].cmy+=(par[k].m*par[k].y)/mtr[par[k].ix][par[k].jy].mass; //centre of mass for y, for a given cell
			if(l==particle_iter-1){//on the last iteration, calculates centre of mass of all particles
				xcm+=(par[k].m*par[k].x)/masssum;
				ycm+=(par[k].m*par[k].y)/masssum;
			}
		}
	}
	printf("%.2f %.2f\n", par[0].x, par[0].y);
	printf("%.2f %.2f\n", xcm, ycm);
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

void Echo(char * msg) {
	if (ECHO_ON)
	{
		printf("Node %d: %s\n", rank, msg);
	}
}

void SendMatrix(long ncside, long long n_part){

	MATRIX *m1=&mtr[0];
	
	for (int i=0; i<p;i++){
		Echo("Sending matrix...");
		MPI_Send(m1, 3*ncside*ncside, MPI_DOUBLE, i, TAG_MTR,  MPI_COMM_WORLD);

		int row, columns;
		for (int row=0; row<ncside*ncside; row++)
		{
        		 printf("%f     ", mtr[row].mass);

 		}
	}

}

void ReceiveMatrix(long ncside){

		MATRIX *mat;
		mat = (MATRIX*)calloc(ncside,sizeof(MATRIX));

		Echo("Receiving matrix...");
		MPI_Recv(&(mat[0]), ncside*ncside*3, MPI_DOUBLE, 0, TAG_MTR, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    	Echo("Matrix Received!");
    	for (int row=0; row<ncside*ncside; row++)
		{
        		 printf("%f     ", mat[row].mass);

 		}
 		printf("\n");
    	
	
}

particle_t * alloc_particles(long long n_part){

	particle_t * part;
	if(n_part%p==0){
		part=(particle_t*)calloc((n_part/p),sizeof(particle_t));
	}
	else{
		if(rank<n_part%p)
			part=(particle_t*)calloc((n_part/p)+1,sizeof(particle_t));
		else
			part=(particle_t*)calloc((n_part/p),sizeof(particle_t));
	}

	return part;

}

void SendParticles(long long n_part){

	particle_t * prt;
	for (int i=0;i<n_part;i++){
		prt=&par[i];
		printf("HEY THERE %f %f\n", par[i].x, par[i].y);
		Echo("Sending particle...");
		MPI_Send(prt, 7, MPI_DOUBLE, i%p, TAG_PAR, MPI_COMM_WORLD);
	}
}

void ReceiveParticle(long long n_part){
	particle_t *part=alloc_particles(n_part);
	int flag = 1, i=0;
	while(1){
		MPI_Iprobe(MPI_ANY_SOURCE, TAG_PAR, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
		if (!flag)
			break;
		Echo("Receiving particle...");
		MPI_Recv(&part[i], 7, MPI_DOUBLE, MPI_ANY_SOURCE, TAG_PAR, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		printf("HEY THERE RANK %d: %f %f\n", rank, part[i].x, part[i].y);
		i++;
    	Echo("Particle Received!");

	}		
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

	

	MPI_Init (&argc, &argv);
 	//MPI_Barrier(MPI_COMM_WORLD);
 	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
 	MPI_Comm_size (MPI_COMM_WORLD, &p);

 	printf("P IS %d\n", p);
	
	if(argc!=5) usage();
		double masssum=0;
		long l;
		char *ptr1, *ptr2, *ptr3, *ptr4;
		const long seed = strtol(argv[1], &ptr1, 10);
		long ncside = strtol(argv[2], &ptr2, 10);
		long long n_part = strtol(argv[3], &ptr3, 10);
		const long particle_iter = strtol(argv[4], &ptr4, 10);
		if (*ptr1!=0 || *ptr2!=0 || *ptr3!=0 || *ptr4!=0 || seed <=0 || ncside<=0 || n_part<=0 || particle_iter<=0) usage();
	
		Create_MPI_particle_t();
		Create_MPI_matrix(ncside);
	//SERIAL
	if(rank==0){
		
		if ((par = (particle_t*)calloc(n_part,sizeof(particle_t)))==NULL) exit (0);

		if ((mtr = (MATRIX*)calloc(ncside*ncside,sizeof(MATRIX)))==NULL) exit (0);

		/*for (l=0; l<ncside; l++){
			if ((mtr[l]=(MATRIX*)calloc(ncside,sizeof(MATRIX)))==NULL) exit (0);
		}*/
	
		/*for(int i=0; i<ncside;i++){
			for(int j=0; j<ncside;j++){
				mtr[i][j].mass=0.0;
				mtr[i][j].cmx=0.0;
				mtr[i][j].cmy=0.0;
			}
		}*/
		init_particles(seed, ncside, n_part, par);
		masssum=centerofmassinit(ncside, n_part, masssum);

		SendMatrix(ncside, n_part);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);

	ReceiveMatrix(ncside);

	MPI_Barrier(MPI_COMM_WORLD);

	if(rank==0)
		SendParticles(n_part);

	MPI_Barrier(MPI_COMM_WORLD);



	ReceiveParticle(n_part);
	

	//RECEIVE INFO
	//wrapcalc(ncside, n_part, particle_iter, masssum);*/

	Echo("Finished!");

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();

	exit(0);

}

