/** \file quantfdtd.c 
    This program calculates eigenenergies and eigenfunctions of the nanoscale
	QD on the surface.
	
	Author: T.Ivas
	Date  : 08.02.13
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include "quantfdtd.h"

/**
  Potential of the free particle
  */		
void pot_free(double **V, int N, int M){
	for(int i = 0; i < N; i++){
		for(int j = 0; j < M; j++){
			//printf("i: %i, j: %i\n",i,j);
			V[i][j] = 0.0;
		}
	}
}

void pot_step(double **V, int N, int M){
	for(int i = 0; i < N; i++){
		for(int j = 0; j < M; j++){
			V[i][j] = 0.0;
		}
	}
}

void pot_slanted(double **V, int N, int M){
		
}
/**
 Potential of the harmonic oscillator
 */
void pot_harmonic(double **V, int N, int M, int NC, int MC, double k0, double delX2){
	for( int i = 0; i < N; i++) {
		for ( int j = 0; j < M; j++) {
			V[i][j] = k0 * ((NC-i)*(NC-i)+(MC-j)*(MC-j))*delX2;
		}
	}
}
/**
 Potential for the quantum corral as presented by Cromie et. al 
 Science, 282, 181, (1993)
 */
void pot_corral(double **V, int N, int M, int NC, int MC, double radius, double V0, double DX){
	for( int i = 0; i < N; i++) {
		for ( int j = 0; j < M; j++) {
			if (((NC-i)*(NC-i)+(MC-j)*(MC-j))*DX > radius*radius )
				V[i][j] = V0;
		}
	}
}

/* function calculates if the given point is inside the polygon */
int pnpoly(int nvert, float *vertx, float *verty, float testx, float testy) {
	int i,j, c=0;
	for (i = 0, j = nvert-1; i < nvert; j=i++) {
		printf ("%i,%i\n",i,j);
		if( ((verty[i]>testy) != (verty[j]>testy)) &&
			(testx < (vertx[j]-vertx[i])* (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]))
				c = !c; 
	}
	return c;
}
/**
Reads input from the file descriptor.

\param input_parameters_s represents input parameters to read.
\param input if this set to \c NULL, read from the stdin.
*/
int readinput(input_data_s *inppar, FILE *input){
	// Read input
	FILE *inpd = input ? input : stdin;
	
	printf("# Number of cells in x direction (int)?:"); 
	fscanf(inpd,"%i", &(inppar->NN));
	printf("# Number of cells in y direction (int)?:");
	fscanf(inpd,"%i", &(inppar->MM));
	printf("# Size of the cell in x direction (double)?:");
	fscanf(inpd,"%lf", &(inppar->del_x));
	printf("# Conversion factor for size of the cell in nm (double)?:");
	fscanf(inpd,"%lf", &(inppar->DX));
	printf("# Time step duration (double)?:"); 
	fscanf(inpd,"%lf", &(inppar->dt));
	printf("# Effective mass of an electron (double)?:");
	fscanf(inpd,"%lf", &(inppar->meff));
	printf("# Potential (int)?:");
	fscanf(inpd,"%i", &(inppar->potential));
	printf("# Number of timesteps (int):?");
	fscanf(inpd,"%i", &(inppar->n_step));
	return 0;
}

/**
This function represents the test function which is used to
find eigenenergies and eigenfunctions of the QD.

\param N  number of discrete cells in X direction.
\param M  number of discrete cells in Y direction.
\param win2D Hanning window used to smooth the FFT results. 
\param double ** prl real part of the wavefunction.

\return Initialize and returns real wavefunction.
*/

int test_function( int N, int M, int NC, int MC, double **win2D, double **prl ) {
	double dist,dist1,dist2;
	double sigma = 6.0;
	for(int i = 0; i < N; i++){
		for(int j = 0; j < M; j++){
			win2D[i][j] = 0.25 * (1. - cos(2*M_PI*i/N)) * (1. - cos(2*M_PI*j/M));
			
			// Single point
			dist = sqrt((NC-i)*(NC-i) + (MC-j)*(MC-j));
			prl[i][j]  = win2D[i][j] * exp(-(dist/sigma)*(dist/sigma));
			
			// Double point
			dist1 = sqrt( (MC-j)*(MC-j) + (NC+10-i)*(NC-10-i) );
			dist2 = sqrt( (MC-j)*(MC-j) + (NC-10-i)*(NC-10-i) );								
		}
	}
	return 0;
}

/**
Calculate the total sum of wavefunctions which are used to normalize the probability.

\param N  number of discrete cells in X direction.
\param M  number of discrete cells in Y direction.
\param double ** prl real part of the wavefunction.
\param double ** pim imaginary part of the wavefunction.

\return sum of the wavefunction over the whole space.  
*/
double check_normalization( int N, int M, double **prl, double **pim ) {
	double ptot = 0.0;
	for(int i = 0; i < N; i++){
		for(int j = 0; j < M; j++) {
			ptot += prl[i][j]*prl[i][j]+pim[i][j]*pim[i][j]; 
		}
	}
	return ptot;
}

/**
Do normalization of the real and imaginary part of the wavefunction.

\param N  number of discrete cells in X direction.
\param M  number of discrete cells in Y direction.
\param double **prl real part of the wavefunction.
\param double **pim imaginary part of the wavefunction.

\return Normalized real and imaginary parts of the wavefunction.
*/
double normalization( int N, int M, double **prl, double **pim ) {
	double ptot = check_normalization( N, M, prl, pim );
	double pnorm = sqrt(ptot);
	for(int i = 0; i < N; i++){
		for(int j = 0; j < M; j++) {
			prl[i][j] /= pnorm;
			pim[i][j] /= pnorm; 
		}
	}
	return pnorm;
}
/**
  This is the core of the program it calculates the time step for Schrodinger
  equation using forward stepping.
  
  \param N  number of discrete cells in X direction.
  \param M  number of discrete cells in Y direction.
  \param nsource x point where to collect time-domain data.
  \param msource y point where to collect time-domain data.
  \param input_data_s structure containing all parameters data for calculations.
  \param **prl  real part of the wavefunction.
  \param **pim  imaginary part of the wavefunction.
  \param **V    potential for the particle 
  \param n_step  number of timestep to calculate.
 */
int calculate_FDTD( int N, int M, int nsource, int msource, input_data_s *pars, fft_parameters_s *fft_pars, 
  	double **prl, double **pim, double **V, int n_step ){
	
	int t, i, j;
	double dt, hbar, ra;
	dt = pars->dt;
	hbar = pars->hbar;
	ra = (0.5*hbar/pars->melec)*(dt/(pow(pars->del_x,2))); // ra must be < .1
		
	// ------ This is the core FDTD program ---------------
	for (t=1; t < n_step; t++) {
    	//first real part
		for (i = 1; i < N-1; i++) {
        	for(j = 1; j < M-1; j++) {
            	prl[i][j] = prl[i][j] - ra * (-4.*pim[i][j] + pim[i][j-1] + pim[i][j+1]
            		+ pim[i-1][j] + pim[i+1][j]) + (dt/hbar)*V[i][j]*pim[i][j];
			}
		}
	
        //now imaginary part
    	for (i = 1; i < N-1; i++) {
			for(j = 1; j < M-1; j++) {
            	pim[i][j] = pim[i][j] + ra * (-4*prl[i][j] + prl[i][j-1] + prl[i][j+1]
            		+ prl[i-1][j] + prl[i+1][j]) - (dt/hbar)*V[i][j]*prl[i][j];
			}
		}
		
    	fft_pars->Ptime[t] = prl[msource][nsource]-I*pim[msource][nsource];
	}
	return t;
}

/**
	Creates Hanning window to smooth the FFT data.
	
	\param fft_parameters parameters which are used for FFT.
	\return corrected array Pwin which represents corrected time-domain data
*/
int create_Hanning(fft_parameters_s *pars, int n_step) {
	
	for(int n = 0; n < n_step; n++) {
		pars->win[n] = 0.5*(1.0-cos(2*M_PI*n/n_step));
		pars->Pwin[n] = pars->win[n] * pars->Ptime[n];
	}
	return 0;
}
/**
  FFT routine take discrete fft from time-domain data which are recorded
  on the particular point in space-domain.
  \param fft_parameters contains all parameters for FFT.
  \return FFT transform in fft_parameters.PF array.
 */
int take_fft(fft_parameters_s *pars, int n_step) {
	//fftw_complex *in, *out;
	fftw_plan plan;
	complex double out[Ntime];
	
	//in = (fftw_complex*) fft_malloc(sizeof(fftw_complex) * pars->Ntime);
	//out = (fftw_complex*) fft_malloc(sizeof(fftw_complex) * pars->Ntime);
	
	plan = fftw_plan_dft_1d(Ntime, pars->Pwin, pars->PF, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);

	fftw_destroy_plan(plan);
	
	plan = fftw_plan_dft_1d(Ntime, pars->PF, &out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
	FILE *ftest = fopen("test_fft", "w");
	for(int i = 0; i < n_step; i++)
		fprintf(ftest,"%g \n", out[i]/Ntime);
	
	fftw_destroy_plan(plan);
	//fftw_free(in);
	//fftw_free(out); 
	return 0;
}
int initialize_calculation(input_data_s *pars, fft_parameters_s *fft_pars){
	pars->melec = pars->m0 * pars->meff; // mass of electron
	for (int i = 1; i < pars->NN; i++){
			pars->XX[i] = pars->XX[i-1] + pars->del_x * pars->DX; 
	}
	fft_pars->del_F = 1./(Ntime * pars->dt);
	fft_pars->del_E = pars->J2eV * 2 * M_PI * pars->hbar * fft_pars->del_F;
	fft_pars->delT = pars->dt * 1e12;
	for (int i = 1; i < Ntime; i++) {
		fft_pars->FF[i] += fft_pars->FF[i-1]+fft_pars->del_E;
	}
	for (int i = 1; i < Ntime; i++) {
		fft_pars->PS[i] = fft_pars->PS[i-1] + fft_pars->delT; 
	}
		
	return 0;
}
int plot_data(FILE *output, int N, int M, double *x, double *y, double **prl, double **pim, double **V){
	FILE *fd = output ? output : stdout;
	for (int i=0; i < N; i++){
		for (int j = 0; j < M; j++){
			fprintf(fd, "%g %g %g\n", x[i], y[j], prl[i][j]);
		}
		fprintf(fd,"\n");
	}
	
	return 0;
}

int plot_fft(FILE *output, fft_parameters_s *fft_pars, int n_step) {
	FILE *fd = output ? output : stdout;
	double nnorm = 1./sqrt(Ntime);
	for (int i = 1; i < n_step; i++)
		fprintf(fd, "%g %g\n", fft_pars->PS[i],creal(fft_pars->Pwin[i]));
	fprintf(fd,"\n\n");
	for (int i = 1; i < n_step; i++) {
		fft_pars->PF[i] = nnorm * abs(fft_pars->PF[i]);
		fprintf(fd, "%g %g\n", 1e3*fft_pars->FF[i], creal(fft_pars->PF[i])); 
	}
	return 0;
}

double **matrix(int sizeX, int sizeY){
	double **m;
	printf("Allocating matrix of size: (%i,%i)\n", sizeX, sizeY);
	m = (double **)malloc(sizeof(double*)*sizeY);
	m[0] = (double*)malloc(sizeof(double)*sizeY*sizeX);
	if (m == NULL) 
		fprintf(stderr, "Problem with allocation of memory for matrix\n");
	if (m[0] == NULL)
		fprintf(stderr, "Can't allocate enough memory for matrix\n");
	for (int i = 1; i < sizeY; i++) m[i] = m[i-1]+sizeX; 
	
	return m;
}

void free_matrix(double **m, int sizeX, int sizeY){
	for(int i = 0; i < sizeY; i++)
		free(m[i]);
	free(m);
}	

int allocate_memory(input_data_s *pars, double ***prl, double ***pim, double ***V, double ***win2D){
	pars->XX = (double *)malloc(sizeof(double)*pars->NN);
	// allocating 2D arrays
	*prl = matrix(pars->NN,pars->MM);
	*pim = matrix(pars->NN,pars->MM);
	*V   = matrix(pars->NN,pars->MM);
	*win2D = matrix(pars->NN,pars->MM);
	return 0;
}

int deallocate_memory(input_data_s *pars, double **prl, double **pim, double **V, double **win2D){
	free(pars->XX);
	free_matrix(prl,pars->NN,pars->MM);
	free_matrix(pim,pars->NN,pars->MM);
	free_matrix(V,pars->NN,pars->MM);
	free_matrix(win2D,pars->NN,pars->MM);
	
	return 0;
}

int main(int argc, char **argv) {
	input_data_s inppars = {
		.m0   = 9.1e-31,   /* mass of electron */
		.hbar = 1.054e-34, /* Planck's constant */
		.eV2J = 1.6e-19,   /* Energy conversion factors */
		.J2eV = 1/1.6e-19};

	fft_parameters_s fft_pars; /* FFT data and parameters */
	double k0=0, V0=0, radius=0;
	double **prl = NULL, **pim = NULL, **V = NULL, **win2D = NULL;  /* real and imag part of the wavefunction, potential*/	
	FILE *file = NULL;
	
	// simple processing of command line
	// if (argc == 2) {
	// 		file = argv[1];
	// 		// Open the file to read parameters
	// 		FILE *fileinput = fopen(file,'r');
	// 		Stopif(!fileinput,return -1, "File not found!\n");
	// 	} 
		
	readinput(&inppars, file);
	
	// allocate memory
	allocate_memory(&inppars,&prl,&pim,&V,&win2D);
	initialize_calculation(&inppars, &fft_pars);
	printf("Mass of electron: %g\n",inppars.melec);
	int N = inppars.NN;
	int M = inppars.MM;
	int n_step = inppars.n_step;
	printf("n_steps: %i", inppars.n_step);
	//printf("I got here\n");
	// Specify the potentials
	switch(inppars.potential) {
		case 1:
		pot_free(V, N, M);
		break;
		case 2:
		pot_step(V, N, M);
		break;
		case 3:
		pot_slanted(V, N, M);
		break;
		case 4:
		pot_harmonic(V, N, M, N/2, M/2, k0, pow(inppars.del_x,2));
		break;
		case 5:
		pot_corral(V, N, M, N/2, M/2, radius, V0, inppars.DX);
		break;
		default:
		printf("Unrecognized potential type: %i\n",inppars.potential);
		break;
	}
	//printf("I got here");
	// Test function
	test_function( N, M, N/2, M/2, win2D, prl );
	// Normalize and check 
	normalization( N, M, prl, pim );	
	// Core of FDTD 
	calculate_FDTD( N, M, N/2, M/2-10, &inppars, &fft_pars, prl, pim, V, n_step );
	// Check normalization
	double normalize = check_normalization(N, M, prl, pim);
	printf("Normalization constant after %i steps: %g\n", n_step, normalize);
	// Plot the time domain data and FFT.
	plot_data(NULL, N, M, inppars.XX, inppars.XX, prl, pim, V);
	// Create the Hanning window for the time-domain data
	create_Hanning(&fft_pars, n_step);
	// Take the FFT of the windowed time-domain data 
	take_fft(&fft_pars, n_step);
	// Plot the time-domain data together with the FFT.
	FILE *fdtd = fopen("timedomain","w");
	plot_fft(fdtd, &fft_pars, n_step);
	//deallocate_memory(&inppars,prl,pim,V,win2D);
	return 0;
}