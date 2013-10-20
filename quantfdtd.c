/** \file quantfdtd.c 
    This program calculates eigenenergies and eigenfunctions of the nanoscale
	QD on the surface.
	
	Author: ∃teslos 
	Date  : 08.02.13
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <assert.h>
#include "quantfdtd.h"
int **imatrix(int, int);
void free_imatrix(int **,int,int);
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
	radius /=DX;
	for( int i = 0; i < N; i++) {
		for ( int j = 0; j < M; j++) {
			if ( ((NC-i)*(NC-i)+(MC-j)*(MC-j)) > (radius*radius) )
				V[i][j] = V0*1.6e-19;
		}
	}
}

/* function calculates if the given point is inside the polygon */
int pnpoly(int nvert, float *vertx, float *verty, float testx, float testy) {
	int i,j, c=0;
	for (i = 0, j = nvert-1; i < nvert; j=i++) {
		if( ((verty[i]>testy) != (verty[j]>testy)) &&
			(testx < (vertx[j]-vertx[i])* (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]))
				c = !c; 
	}
	return c;
}

/**
 Potential for the quantum dot as presented by Lobo-Checa et al. 
 */
void pot_hexagon(double **V, int N, int M, int NC, int MC, double radius, double V0, double *XX){
	//float vertx[] ={4.615,13.845,18.46,13.845,4.615,0.1};
	//float verty[] ={0.1,0.1,8.0,16.0,16.0,8.0};
	float r = radius;
	float vertx[6], verty[6];
	for (int i = 0; i < 6; i++) {
		vertx[i] = r*cos(2*M_PI * i / 6.0) + XX[NC];
		verty[i] = r*sin(2*M_PI * i / 6.0) + XX[MC];
	} 

	float testx, testy;
	int nvert = 6;
	for (int i = 0; i < M; i++){
		 testy = XX[i];
		 for (int j = 0; j < N; j++){
			testx = XX[j];
			if (!pnpoly(nvert,vertx,verty,testx,testy)){
				 V[i][j] = V0*1.6e-19;
			 }
		 }
	 }	
}

/**
 Potential which simulates honeycomb network, contains 7 hexagon units arange
 in honeycomb structure.
 
 
 */
void pot_honeycomb(double **V, int N, int M, int NC, int MC, double radius, double t, double V0, double *XX) {
	int **flag;
	flag = (int**)malloc(sizeof(int*)*M);
	flag[0] = (int*)malloc(sizeof(int)*N*M);
	if (flag == NULL)
		fprintf(stderr, "Problem with allocation of memory for flag array\n");
	if (flag[0] == NULL)
		fprintf(stderr, "Can't allocate enough memory for flag array\n");
	for (int i=1; i < M; i++) flag[i] = flag[i-1]+N;

	// create honeycomb structure
	float comb_distance = 2*radius*cos(M_PI / 6.0) + t;
	// calculate six distances from the central hexagon
	float comb_vertx[7], comb_verty[7];
	// starting point in the center of the domain
	comb_vertx[0] = XX[NC];
	comb_verty[0] = XX[MC];
	for (int i = 1; i < 7; i++) {
		comb_vertx[i] = comb_distance*cos(2 * M_PI * i / 6.0)+comb_vertx[0];
		comb_verty[i] = comb_distance*sin(2 * M_PI * i / 6.0)+comb_verty[0];
	}
	// first flag all cells which have non-zero potential
	// we use the same procedure as in the hexagon case
	// but add additional comb cells around.
	for (int k = 0; k < 7; k++) {
		float r = radius;
		float vertx[6], verty[6];
		for (int i = 0; i < 6; i++) {
			vertx[i] = r*cos(2*M_PI * i / 6.0+M_PI/6.0) + comb_vertx[k];
			verty[i] = r*sin(2*M_PI * i / 6.0+M_PI/6.0) + comb_verty[k];
		} 

		float testx, testy;
		int nvert = 6;
		for (int i = 0; i < M; i++){
			testy = XX[i];
			for (int j = 0; j < N; j++){
				testx = XX[j];
				if (pnpoly(nvert,vertx,verty,testx,testy)){
					flag[i][j] = 1;
				}
			}
		}	
	}
	
	// go finally over the complete space and set potential
	// for flagged cells.
	for (int i = 0; i < M; i++){
		for (int j = 0; j < N; j++){
			if (!flag[i][j])
				V[i][j] = V0*1.6e-19;
		}
	} 
	// free memory
	free(flag[0]);
	free(flag);
}
/**
 Reads potential from the 8-bit based image in PGM format. 
 \param double **V - potential to set 
 \param int N - number of discrete cells in X-direction
 \param int M - number of discrete cells in Y-direction
 \param char *filename - contains file name of the PGM image to read.
 \return The potential V set by the pixels of the image.
 */
void pot_image(double **V, int N, int M, char *filename){
    int numcols, numrows;
    int maxgrey;
    
    char *buffer = NULL;
    size_t len = 0;
    int **array;
    // open file for reading
    FILE *fp = fopen(filename,"rb");
    if (fp == NULL)
        fprintf(stderr, " Problem with opening file: %s\n", filename);
    // read the magic number
    getline(&buffer, &len, fp);
    if (strncmp(buffer,"P5",2) != 0) 
        fprintf(stderr, " The image %s is not in PGM format\n", filename);
    // read comment line
    getline(&buffer, &len, fp);
    printf("#Comment:%s \n", buffer);
    
    // read the size of image
    fscanf(fp,"%i %i", &numcols, &numrows);
    printf("#Num cols: %i, num rows: %i\n", numcols, numrows);
    
    fscanf(fp,"%i", &maxgrey);
    printf("#Max grey level: %i\n", maxgrey);

    // allocate memory for image
    array = imatrix(numrows, numcols);
    // read the data in 
    for (int row = 0; row < numrows; ++row){
        for (int col = 0; col < numcols; ++col){
            array[row][col] = (int) fgetc(fp);
        }
    }

    if ((N % numrows != 0) || (M % numcols != 0)) {
        fprintf(stderr, "Size of the domain (%i,%i) is not multiple image size (%i,%i)\n"
            "Please change the cells number to be multiple of image size:\n"
            " example 70x70 image size; 140x140 size of the domain\n", N,M, numrows,numcols);
        exit(0);
    }
    // finally set the potential array with the pixel values from the
    // image. Maximum grey component corresponds to potential value V0.
    for (int i = 0; i < N; ++i){
        for (int j = 0; j < M; ++j){
            V[i][j] = array[i][j]/((float)maxgrey)*1.6e-19*2.16;  // this is workround
        }
    }
    free(buffer);
    free_imatrix(array,numcols,numrows);
} 

/** 
Reads the line until it finds keyword in the file.

\param input is file input
\param keyword to search in file 
\param type of parameter
\param parameter to read
\param comment char that starts comment line 

*/
int readline(FILE *input, char *keyword, char type, void* parameter, char comment){
	char line[512];
	int val;
	double fval;
	char *result;
	
	while (fgets(line, 512, input) != NULL) {
		if (strchr(line,comment) != NULL) continue; // it is comment, next line
		if ((strstr(line, keyword)) != NULL) {
			printf("%s:", keyword);
			result = strtok(line, "=");
			if (result != NULL) {
				result = strtok(NULL, "=");
				if (type == 'I') {
					val = atoi(result);
					int *p = (int *) parameter;
					*p = val;
					printf("%i\n", val);
				} else if (type == 'D') {
					fval = atof(result);
					double *pv = (double *)parameter;
					*pv = fval;
					printf("%e\n",fval);
				} else
					printf("Error: unknown type\n");
			}
			// you found the keyword rewind the file to begin
			rewind(input);
			break;
		}
	}
	return 0;
}

/**
Reads input from the file descriptor.

\param input_parameters_s represents input parameters to read.
\param input if this set to \c NULL, read from the stdin.
*/
int readinput(input_data_s *inppar, FILE *input){
	
	// Read input
	FILE *inpd = input ? input : stdin;
	
	readline(inpd, "CELLSX", 'I', &(inppar->NN), '#');
	readline(inpd, "CELLSY", 'I', &(inppar->MM), '#');
	readline(inpd, "SIZEX", 'D', &(inppar->del_x), '#');
	readline(inpd, "CONV", 'D', &(inppar->DX), '#');
	readline(inpd, "TIMESTEP", 'D', &(inppar->dt), '#');
	readline(inpd, "EMASS", 'D', &(inppar->meff), '#');
	readline(inpd, "POTENTIAL", 'I', &(inppar->potential), '#');
	readline(inpd, "NSTEPS", 'I', &(inppar->n_step), '#');	
	readline(inpd, "EIGEN", 'D', &(inppar->Ein), '#');
	readline(inpd, "TESTX", 'I', &(inppar->testx), '#');
	readline(inpd, "TESTY", 'I', &(inppar->testy), '#');
	
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
	double sigma = 2.0;
	for(int i = 0; i < N; i++){
		for(int j = 0; j < M; j++){
			win2D[i][j] = 0.25 * (1. - cos(2*M_PI*i/N)) * (1. - cos(2*M_PI*j/M));
			
			// Single point
			dist = sqrt((NC-i)*(NC-i) + (MC-j)*(MC-j));
			prl[i][j]  = win2D[i][j] * exp(-(dist/sigma)*(dist/sigma));
			
			// Double point
			dist1 = sqrt( (MC-j)*(MC-j) + (NC-10-i)*(NC-10-i) );
			dist2 = sqrt( (MC-j)*(MC-j) + (NC+10-i)*(NC+10-i) );
			//prl[i][j] = /*win2D[i][j]*/( exp(-(dist1/sigma)*(dist1/sigma))-exp(-(dist2/sigma)*(dist2/sigma)));								
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
int calculate_FDTD( int N, int M, int nsource, int msource, input_data_s *pars, fft_parameters_s *fft_pars, eigenfunction_s *eigen,
  	double **prl, double **pim, double **V, int n_step ){
	
	int t, i, j;
	double dt, hbar, ra;
	dt = pars->dt;
	hbar = pars->hbar;
	ra = (0.5*hbar/pars->melec)*(dt/(pow(pars->del_x,2))); // ra must be < .1
	printf("# Value of Ra constant: %g\n", ra);
	FILE *phi_fd = fopen("phi_test.txt","w");	
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
		
		if (!isnan(pars->Ein)){
		    for (i = 1; i < N-1; i++) {
		        for (j = 1; j < M-1; j++) {
		            eigen->psi[i][j]  = prl[i][j] + I*pim[i][j];
		            eigen->phi[i][j] += fft_pars->win[t]* (cos(eigen->arg*t)-I*sin(eigen->arg*t))*(eigen->psi[i][j]);
					//printf("win: %g\n", fft_pars->win[t]);
				}
				
			}
			// print test
		    // for (i = 0; i < N; i++) {
// 				for (j = 0; j < M; j++) {
// 					fprintf(phi_fd,"%g ",creal(eigen->phi[i][j]));
// 				}
// 				fprintf(phi_fd,"\n");
// 			}
// 			fclose(phi_fd);
// 			exit(0);
		}
		
		
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
	
	FILE *fw = fopen("test_win","w");
	for (int n = 0; n < n_step; n++)
		fprintf(fw, "%i %g \n", n, creal(pars->win[n]));
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
	
	//in = (fftw_complex*) fft_malloc(sizeof(fftw_complex) * pars->Ntime);
	//out = (fftw_complex*) fft_malloc(sizeof(fftw_complex) * pars->Ntime);
	
	plan = fftw_plan_dft_1d(Ntime, pars->Pwin, pars->PF, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);

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
/**
 This function writes results to a file used by gnuplot to plot the data.
 
 \param output file descriptor if set \c NULL then write to the stdout.
 \param N  number of cells in x direc.
 \param M  number of cells in y direc.
 \param *x coordinates to write out in x direc.
 \param *y coordinates to write out in y direc.
 \param **prl  real part of the wavefunction or eigenfunction.
 \param **pim  imag part of the wavefunction or eigenfunction.
 \param **V    potential used for calculation.
 */
int plot_data(FILE *output, int N, int M, double *x, double *y, double **prl, double **pim, double **V){
	FILE *fd = output ? output : stdout;
	if (prl != NULL) {
		for (int i=0; i < N; i++){
			for (int j = 0; j < M; j++){
				fprintf(fd, "%g %g %g\n", x[i], y[j], prl[i][j]);
			}
			fprintf(fd,"\n");
		}
	}
	fprintf(fd,"\n\n");
	if (pim != NULL) {
		for (int i=0; i < N; i++){
			for (int j = 0; j < M; j++){
				fprintf(fd, "%g %g %g\n", x[i], y[j], pim[i][j]);
			}
			fprintf(fd,"\n");
		}
	}
	if (V != NULL)
	for (int i=0; i < N; i++){
		for (int j = 0; j < M; j++){
			fprintf(fd, "%g %g %g\n", x[i], y[j], V[i][j]);
		}
		fprintf(fd,"\n");
	}
	return 0;
}
/**
 Save FFT data to be ploted by gnuplot.
 \param FILE *output - file handle of the output file
 \param fft_parameters_s *fft_pars - results of fft analysis
 \param int n_step - number of the stime steps to save (display)
 */
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

/**
 Allocate memory for the matrix 2d array.
 \param int sizeX size of the matrix in x 
 \param int sizeY size of the matrix in y
 \return returns allocate memory
 */
double **matrix(int sizeX, int sizeY){
	double **m;
	printf("#Allocating matrix of size: (%i,%i)\n", sizeX, sizeY);
	m = (double **)malloc(sizeof(double*)*sizeY);
	m[0] = (double*)malloc(sizeof(double)*sizeY*sizeX);
	if (m == NULL) 
		fprintf(stderr, "Problem with allocation of memory for matrix\n");
	if (m[0] == NULL)
		fprintf(stderr, "Can't allocate enough memory for matrix\n");
	for (int i = 1; i < sizeY; i++) m[i] = m[i-1]+sizeX; 
	
	return m;
}
/**
 Allocate memory for the matrix 2d array of ints.
 \param int sizeX size of the matrix in x 
 \param int sizeY size of the matrix in y
 \return returns allocate memory
 */
int **imatrix(int sizeX, int sizeY){
	int **m;
	printf("#Allocating matrix of size: (%i,%i)\n", sizeX, sizeY);
	m = (int **)malloc(sizeof(int*)*sizeY);
	m[0] = (int*)malloc(sizeof(int)*sizeY*sizeX);
	if (m == NULL) 
		fprintf(stderr, "Problem with allocation of memory for matrix\n");
	if (m[0] == NULL)
		fprintf(stderr, "Can't allocate enough memory for matrix\n");
	for (int i = 1; i < sizeY; i++) m[i] = m[i-1]+sizeX; 
	
	return m;
}

/**
 Allocate memory for the matrix 2d array of complex numbers.
 \param int sizeX size of the matrix in x 
 \param int sizeY size of the matrix in y
 \return returns allocated memory
*/
complex double **cmatrix(int sizeX, int sizeY){
	complex double **m;
	printf("#Allocating matrix of size: (%i,%i)\n", sizeX, sizeY);
	m = (complex double **)malloc(sizeof(complex double*)*sizeY);
	m[0] = (complex double*)malloc(sizeof(complex double)*sizeY*sizeX);
	if (m == NULL) 
		fprintf(stderr, "Problem with allocation of memory for matrix\n");
	if (m[0] == NULL)
		fprintf(stderr, "Can't allocate enough memory for matrix\n");
	for (int i = 1; i < sizeY; i++) m[i] = m[i-1]+sizeX; 
	
	return m;
}
/**
 Free previous memory allocated by \c matrix.
 \param double **m is matrix to be freed.
 \param int sizeX size of the matrix in X direction 
 \param int sizeY size of the matrix in Y direction
 */
void free_matrix(double **m, int sizeX, int sizeY){
	free(m[0]);
	free(m);
}	

/**
 Free previous memory allocated by \c imatrix.
 \param double **m is matrix to be freed.
 \param int sizeX size of the matrix in X direction 
 \param int sizeY size of the matrix in Y direction
 */
void free_imatrix(int **m, int sizeX, int sizeY){
	free(m[0]);
	free(m);
}	
/**
 Free previous memory allocated by \c cmatrix.
 \param double **m is matrix to be freed.
 \param int sizeX size of the matrix in X direction 
 \param int sizeY size of the matrix in Y direction
 */
void free_cmatrix(complex double **m, int sizeX, int sizeY){
	free(m[0]);
	free(m);
}
/**
 Do allocation of the memory for the simulation.
 \param input_data_s *pars  input para
 \param fft_parameters_s *fft_pars
 \param eigenfunction_s *eigen - contains data and parameters for the eigenfunction calculations 
 \param double ***prl - 2D array representing real part of the wave function
 \param double ***pim - 2D array representing imaginary part of the wave function
 \param double ***V - 2D array representing arbitary potential for Schrodinger equation
 \param double ***win2D - 2D array representing Hanning window 
 */
long int allocate_memory(input_data_s *pars, fft_parameters_s *fft_pars, eigenfunction_s *eigen, double ***prl, double ***pim, double ***V, double ***win2D){
	pars->XX = (double *)malloc(sizeof(double)*pars->NN);
	// allocating 2D arrays
	//printf("NN: %i, MM: %i\n", pars->NN, pars->MM);
	
	*prl = matrix(pars->NN,pars->MM);
	*pim = matrix(pars->NN,pars->MM);
	*V   = matrix(pars->NN,pars->MM);
	*win2D = matrix(pars->NN,pars->MM);
	if (!isnan(pars->Ein)) {
		eigen->psi = cmatrix(pars->NN,pars->MM);
		eigen->phi = cmatrix(pars->NN,pars->MM);
		eigen->phi_m = matrix(pars->NN,pars->MM);
		eigen->angle = matrix(pars->NN,pars->MM);
		eigen->phi0_rl = matrix(pars->NN,pars->MM);
	}
	// allocate memory for FFT parameters
	fft_pars->Ptime = (complex double*) malloc(sizeof(complex double)*Ntime);
	fft_pars->Pwin  = (complex double*) malloc(sizeof(complex double)*Ntime);
	fft_pars->win   = (double *) malloc(sizeof(double)*Ntime);
	fft_pars->PF    = (complex double *) malloc(sizeof(complex double)*Ntime);
	fft_pars->FF    = (double *)malloc(sizeof(double)*Ntime);
	fft_pars->PS    = (double *)malloc(sizeof(double)*Ntime);
	
	// for (int i = 0; i < pars->NN; i++) {
// 		for (int j = 0; j < pars->MM; j++) {
// 			eigen->psi[i][j] = 0.0;
// 			eigen->phi[i][j] = 0.0;
// 			eigen->phi_m[i][j] = 0.0;
// 			eigen->angle[i][j] = 0.0;
// 			eigen->phi0_rl[i][j] = 0.0;
// 		}
// 	}
	long int memory_allocated = sizeof(double)*(pars->NN + 7*pars->NN*pars->MM +3*Ntime) +
		sizeof(complex double)*(2*pars->NN*pars->MM+3*Ntime);
	return memory_allocated;
}
/**
 Deallocate memory previously allocated by allocate function.
 */
int deallocate_memory(input_data_s *pars, fft_parameters_s *fft_pars, eigenfunction_s *eigen, double **prl, double **pim, double **V, double **win2D){
	free(pars->XX);
	free_matrix(prl,pars->NN,pars->MM);
	free_matrix(pim,pars->NN,pars->MM);
	free_matrix(V,pars->NN,pars->MM);
	free_matrix(win2D,pars->NN,pars->MM);
	if (!isnan(pars->Ein)) {
		free_cmatrix(eigen->psi,pars->NN,pars->MM);
		free_cmatrix(eigen->phi,pars->NN,pars->MM);
		free_matrix(eigen->phi_m,pars->NN,pars->MM);
		free_matrix(eigen->angle,pars->NN,pars->MM);
		free_matrix(eigen->phi0_rl,pars->NN,pars->MM);
	}
	
	free(fft_pars->Ptime);
	free(fft_pars->Pwin);
	free(fft_pars->win);
	free(fft_pars->PF);
	free(fft_pars->FF);
	free(fft_pars->PS);
	
	return 0;
}
/**
 Intialize the parameters like frequency, period -- for calculating eigenfunction
 of the Schrodinger equation.
 \param input_data_s *pars represents input data (eigenenergy, physical constants).
 \param eigenfunction_s *eigen eigenfunction parameters.
 */
int eigenfunction(input_data_s *pars, eigenfunction_s *eigen) {
	if ( !isnan(pars->Ein) ) {
		eigen->freq = pars->Ein/(pars->J2eV*2*M_PI*pars->hbar);
		eigen->omega = 2*M_PI*eigen->freq;
		eigen->arg = eigen->omega * pars->dt;
		eigen->T_period = 1/(eigen->freq * pars->dt);
		printf("#freq: %g, omega: %g, arg: %g, T_period = %g\n", eigen->freq, eigen->omega,
			eigen->arg, eigen->T_period);
	}
	return 0;
}
/**
 Finds eigenfunction of the problem.
 
 \param int NC  x position of the source  
 \param int MC  y position of the source 
 \param input_data_s contains input data parameters (geometry) 
 \param eigenfunction_s *eigen eigenfunction parameters
 */
int find_eigenfunction(int NC, int MC, input_data_s *pars, eigenfunction_s *eigen){
	complex double ptot_phi = 0 + I*0;
	// shorthand writing
	complex double **phi = eigen->phi;
	double **phi_m = eigen->phi_m;
	double **angle = eigen->angle;
	double **phi0_rl = eigen->phi0_rl;
	
	for (int i = 0; i < pars->NN; i++) {
		for (int j = 0; j < pars->MM; j++) {
			ptot_phi += phi[i][j]*conj(phi[i][j]);   //complex conjugate
		}
	}
	//printf("Ptot_phi:%g + I%g\n", creal(ptot_phi),cimag(ptot_phi));
	for (int i = 0; i < pars->NN; i++) {
		for (int j = 0; j < pars->MM; j++) {
			phi_m[i][j] = phi[i][j] / ptot_phi;
			//printf("phi_m: %g \n", phi_m[i][j]);
			assert(!isnan(phi_m[i][j]));
			angle[i][j] = atan2(cimag(phi[i][j]),creal(phi[i][j]));
			//printf("i: %i, j: %i\n",i,j);
			//printf("phi: %g + I%g, angle: %g\n", creal(phi[i][j]),cimag(phi[i][j]), angle[i][j]);
			assert(!isnan(angle[i][j]));
		}
	}
	
	FILE *fdtest = fopen("test_phi1_m", "w");
	for (int i = 0; i < pars->NN; i++){
		for (int j = 0; j < pars->MM; j++){
			fprintf(fdtest, "%i %i %g\n", i,j, creal(phi_m[i][j]));
		}
		fprintf(fdtest,"\n");
	}
	fprintf(fdtest, "\n\n");
	for (int i = 0; i < pars->NN; i++){
		for (int j = 0; j < pars->MM; j++){
			fprintf(fdtest, "%i %i %g\n", i,j, cimag(phi_m[i][j]));
		}
		fprintf(fdtest,"\n");
	}
	fclose(fdtest);
	
	double ang0 = angle[MC][NC];      // The angle at the source point
	//double ang0 = 0.0;
	printf("#ang0: %g\n", ang0);
	double ptot0 = 0.;
	for (int i = 0; i < pars->NN; i++) {
		for (int j = 0; j < pars->MM; j++) {	
			angle[i][j] = angle[i][j] - ang0;
			assert(!isnan(angle[i][j]));
			//printf("abs(phi): %g, cos(angle): %g\n", fabs(phi_m[i][j]),cos(angle[i][j]));
			phi0_rl[i][j] = fabs(phi_m[i][j])*cos(angle[i][j]);
			assert(!isnan(phi0_rl[i][j]));
			//printf("phi0_rl: %g\n", phi0_rl[i][j]);
			ptot0 += phi0_rl[i][j]*phi0_rl[i][j];
		}
	}
	//printf("Ptot: %g\n", ptot0);
	double nptot0 = sqrt(ptot0);

	for (int i = 0; i < pars->NN; i++) {
		for (int j = 0; j < pars->MM; j++) {	
			phi0_rl[i][j] = phi0_rl[i][j] / nptot0;
		}
	}
	return 0;
}
/**
 Calculates critical timestep 
 \return critical time-step in sec.
 */
double critical_dt(input_data_s *pars, double Vmax){
	return pars->hbar/(pow(pars->hbar,2)/pars->melec *(2.0/pars->del_x)+Vmax);
}

int main(int argc, char **argv) {
	input_data_s inppars = {
		.m0   = 9.1e-31,   /* mass of electron */
		.hbar = 1.054e-34, /* Planck's constant */
		.eV2J = 1.6e-19,   /* Energy conversion factors */
		.J2eV = 1/1.6e-19};

	fft_parameters_s fft_pars; /* FFT data and parameters */
	eigenfunction_s  eigen;    /* eigenfunction parameters and data */
	int MC, NC;                /* source point for the finding eigenfunction */
    double k0=0, V0=0.8, radius=13.8/2.0, t= 0.7;
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
	long int mem_allocated = allocate_memory(&inppars,&fft_pars,&eigen,&prl,&pim,&V,&win2D);
	printf("#Current Memory Allocated: %li\n", mem_allocated);
	initialize_calculation(&inppars, &fft_pars);
	eigenfunction(&inppars,&eigen);
	printf("#Mass of electron: %g\n",inppars.melec);
	printf("#Critical timestep: %g\n", critical_dt(&inppars,V0*inppars.eV2J));
	printf("#Timestep: %g\n",inppars.dt);
	printf("#Timesteps: %i\n",inppars.n_step);
	printf("#Eigenenergy: %g\n",inppars.Ein);
	int N = inppars.NN;
	int M = inppars.MM;
	int n_step = inppars.n_step;
	printf("#n_steps: %i", inppars.n_step);
	// test point
    NC = inppars.testx;
    MC = inppars.testy;
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
		pot_corral(V, N, M, N/2, M/2, radius, V0, inppars.del_x * inppars.DX);
		break;
		case 6:
		pot_hexagon(V, N, M, N/2, M/2, radius, V0, inppars.XX);
		break;
		case 7:
		pot_honeycomb(V, N, M, N/2, M/2, radius, t, V0, inppars.XX);
		break;
        case 8:
        pot_image(V, N, M, "TV151.pgm");
        break;
		default:
		printf("Unrecognized potential type: %i\n",inppars.potential);
		break;
	}
	// Test function
	test_function( N, M, NC, MC, win2D, prl );
	
	// Normalize and check 
	normalization( N, M, prl, pim );	
	
	if (!isnan(inppars.Ein)) 
		create_Hanning(&fft_pars, n_step);
	
	// Core of FDTD 
	calculate_FDTD( N, M, N/2, M/2, &inppars, &fft_pars, &eigen, prl, pim, V, n_step );
	
	// Check normalization
	double normalize = check_normalization(N, M, prl, pim);
	printf("#Normalization constant after %i steps: %g\n", n_step, normalize);
	
	if (!isnan(inppars.Ein)) {
		FILE *find_phi = fopen("Eigenfunction","w");
		find_eigenfunction(NC, MC, &inppars, &eigen);
		plot_data(find_phi, N, M, inppars.XX, inppars.XX, prl, pim, NULL);
		fclose(find_phi);
		plot_data(NULL, N, M, inppars.XX, inppars.XX, eigen.phi0_rl, NULL, NULL);	
	} else {
		// Plot the time domain data and FFT.
		FILE *poten = fopen("Potential","w");
		plot_data(NULL, N, M, inppars.XX, inppars.XX, prl, pim, NULL);
		plot_data(poten, N, M, inppars.XX, inppars.XX, NULL, NULL, V);
		fclose(poten);
		// Create the Hanning window for the time-domain data
		create_Hanning(&fft_pars, n_step);
		// Take the FFT of the windowed time-domain data 
		take_fft(&fft_pars, n_step);
		// Plot the time-domain data together with the FFT.
		FILE *fdtd = fopen("timedomain","w");
		plot_fft(fdtd, &fft_pars, n_step);
	}
	deallocate_memory(&inppars,&fft_pars,&eigen,prl,pim,V,win2D);
	return 0;
}