//#include <glib.h>
#include "quantfdtd.h"
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

input_data_s data = {
		.m0   = 9.1e-31,   /* mass of electron */
		.hbar = 1.054e-34, /* Planck's constant */
		.eV2J = 1.6e-19,   /* Energy conversion factors */
		.J2eV = 1/1.6e-19};

fft_parameters_s fft_data;

int readinput(input_data_s *inppar, FILE *input){
	// Read input
	FILE *inpd = input ? input : stdin;
	
	printf("# Number of cells in x direction (int)?:"); 
	fscanf(inpd,"%i", &(inppar->NN));
	printf("# Number of cells in y direction (int)?:");
	fscanf(inpd,"%i", &(inppar->MM));
	printf("# Size of the cell in x direction (double)?:");
	fscanf(inpd,"%lf", &(inppar->del_x));
	printf("# Time step duration (double)?:"); 
	fscanf(inpd,"%lf", &(inppar->dt));
	printf("# Effective mass of an electron (double)?:");
	fscanf(inpd,"%lf", &(inppar->meff));
	
	printf("Number of cells in x: %i\n", inppar->NN);
	printf("Number of cells in y: %i\n", inppar->MM);
	printf("Size of cell in x: %g\n", inppar->del_x);
	printf("Time step duration: %g\n", inppar->dt);
	printf("Effective mass of an electron: %g\n", inppar->meff);
	return 0;
}

double **matrix(int sizeX, int sizeY){
	double **m;
	m = (double **)malloc(sizeof(double*)*sizeY);
	m[0] = (double*)malloc(sizeof(double)*sizeY*sizeX);
	
	for (int i = 1; i < sizeY; i++) m[i] = m[i-1]+sizeX; 
	
	return m;
}
	
int allocate_memory(input_data_s *pars, double ***prl, double ***pim, double ***V,double ***win2D){
	pars->XX = (double *)malloc(sizeof(double)*pars->NN);
	// allocating 2D arrays
	*prl = matrix(pars->NN,pars->MM);
	*pim = matrix(pars->NN,pars->MM);
	*V   = matrix(pars->NN,pars->MM);
	*win2D = matrix(pars->NN,pars->MM);
	return 0;
}

int check_allocate(input_data_s *pars, double **prl, double **pim) {
	printf("Data size: %li\n",sizeof(prl));
	assert(pars->NN == 10);
	assert(pars->MM == 10);
	// write in memory
	for (int j = 0; j < pars->MM; j++ ) {
		for (int i = 0; i  < pars->NN; i++ ) {
			printf("i: %i, j: %i\n",i,j);
			prl[j][i] = 0xBABA;
		}
	}
	assert(prl[9][9] == 0xBABA);
	assert(prl[0][0] == 0xBABA);
	return 0;
}

void pot_free(double **V, int N, int M){
	for(int i = 0; i < N; i++){
		for(int j = 0; j < M; j++){
			V[i][j] = 0.0;
		}
	}
}

int check_potential(double **V){
	int potential = 0;
	switch(potential){
		case 0:
		assert(V[0][0]==0.0);
		assert(V[9][9]==0.0);
		break;
		default:
		break;
	}
	return 0;
}

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

int check_function(int NC, int MC, double **win2D, double **prl){
	assert(win2D[NC][MC] == prl[NC][MC]);
	return 0;
}

double check_normalization( int N, int M, double **prl, double **pim ) {
	double ptot = 0.0;
	for(int i = 0; i < N; i++){
		for(int j = 0; j < M; j++) {
			ptot += prl[i][j]*prl[i][j]+pim[i][j]*pim[i][j]; 
		}
	}
	return ptot;
}

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

int calculate_FDTD( int N, int M, int nsource, int msource, input_data_s *pars, fft_parameters_s *fft_pars, 
  	double **prl, double **pim, double **V, int n_step ){
	
	int t, i, j;
	double dt, hbar, ra;
	dt = pars->dt;
	hbar = pars->hbar;
	ra = (0.5*hbar/pars->melec)*(dt/(pow(pars->del_x,2))); // ra must be < .1
	printf("dt: %g, hbar: %g, melec: %g, ra:%g\n", dt, hbar, pars->melec, ra);	
	// ------ This is the core FDTD program ---------------
	for (t=0; t < n_step; t++) {
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

int check_FDTD(double **prl, double **pim, fft_parameters_s *fft_data){
	int t = 0;
	printf("Real:%g, Imag: %g\n", creal(fft_data->Ptime[t]),cimag(fft_data->Ptime[t]));
	return 0;
}

int create_Hanning(fft_parameters_s *pars, int n_step) {
	
	for(int n = 1; n < n_step; n++) {
		pars->win[n] = 0.5*(1.0-cos(2*M_PI*n/n_step));
		pars->Pwin[n] = pars->win[n] * pars->Ptime[n];
	}
	return 0;
}

int check_Hanning(fft_parameters_s *fft_data){
	assert(fft_data->Pwin[1] == 0.0);
	return 0;
}

int take_fft(fft_parameters_s *pars) {
	//fftw_complex *in, *out;
	fftw_plan plan;
	
	//in = (fftw_complex*) fft_malloc(sizeof(fftw_complex) * pars->Ntime);
	//out = (fftw_complex*) fft_malloc(sizeof(fftw_complex) * pars->Ntime);
	
	plan = fftw_plan_dft_1d(Ntime, pars->Pwin, pars->PF, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
	
	fftw_destroy_plan(plan);
	return 0;
	//fftw_free(in);
	//fftw_free(out); 
}

int check_fft(fft_parameters_s *pars){
	for (int i = 0; i < 10; i++)
		printf("PF:%g + I%g\n", creal(pars->PF[i]), cimag(pars->PF[i]));
}

int main() {
	FILE *fd = NULL;
	double **prl=NULL,**pim=NULL,**V=NULL,**win2d=NULL;
	readinput(&data, fd);
	// mass of electron
	data.melec = data.m0 * data.meff;
	
	allocate_memory(&data,&prl,&pim,&V,&win2d);
	check_allocate(&data,prl,pim);
	pot_free(V,data.NN, data.MM);
	check_potential(V);
	test_function( data.NN, data.MM, data.NN/2, data.MM/2, win2d, prl );
	check_function( data.NN/2, data.MM/2, win2d, prl );
	printf("Normalization constant: %g\n", check_normalization( data.NN, data.MM, prl, pim ));
	normalization( data.NN, data.MM, prl, pim );
	printf("Normalization constant: %g\n", check_normalization( data.NN, data.MM, prl, pim ));
    calculate_FDTD(data.NN, data.MM, data.NN/2, data.MM/2, &data, &fft_data, prl, pim, V, 1);
	check_FDTD(prl, pim, &fft_data);
	printf("Normalization constant after one step: %g\n", check_normalization( data.NN, data.MM, prl, pim ));
	create_Hanning(&fft_data,1);
	take_fft(&fft_data);
	check_fft(&fft_data);	
	return 0;
}