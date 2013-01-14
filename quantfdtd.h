/* quantfdtd.h - contains all constants and definitions of data structures
   for the quantfdtd program.
*/
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

char error_mode;

FILE *error_log;

#define Stopif(assertion, error_action, ...)                      \
	   if (assertion){                                            \
		   fprintf(error_log ? error_log : stderr, __VA_ARGS__);  \
			   fprintf(error_log? error_log: stderr, "\n");       \
				   if (error_mode=='s') abort();                  \
	               else {error_action;}                           \			   
	   }

typedef struct {
   int NN;                            // Number of points in the space.
   int MM;
   const double hbar;                 // Planck's constant
   const double m0;                   // Mass of an free electron
   double meff;         	          // Effective mass of an electron for Cu(111)
   double melec;                      // Mass of an electron
   const double eV2J;                 // Energy conversion factors
   const double J2eV;

   double del_x;       // The cells size
   double dt;          // Time steps
   double ra;          // Stability criteria
   double DX;          // Delta X in real units          
   double *XX;         // Mesh of points in X direc.
   int potential;      // potential to calculate.
   int n_step;         // Number of time steps.
   }input_data_s;
	   
#define  Ntime  (2<<15)        // Size of the FFT buffer  

typedef struct {
   complex double Ptime[Ntime];          // wavefunction saved at particular point    
   complex double Pwin[Ntime];           // wavefunction corrected by Hanning
   complex double win[Ntime];            // Size of Hanning window
   complex double PF[Ntime];             // FFT transform of the wavefunction
   double FF[Ntime]; 
   double PS[Ntime];
   double del_F;               // delta in Fourier space
   double del_E;               // delta in Energy
   double delT;                // Delta time scaled
   }fft_parameters_s;