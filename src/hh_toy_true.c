/*
 *   homogeneous half space
 *   last update 02.11.02, T. Bohlen
 */

#include "fd.h"

void model(float  ***  rho, float ***  pi, float ***  u,
           float ***  taus, float ***  taup, float   *eta) {

	/*-----------------------------------------------------------------------*/
	/* extern variables */

	extern float DT, *FL, TAU;
	extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], L;
	extern FILE *FP;
	/* local variables */
	float muv, piv, ws;
	float Vp, Vs, Rho;
	float *pts=NULL, ts=0.0, tp=0.0, sumu, sumpi;
	int i, j, k, l, ii, jj, kk;
	int kasten=1;

	/* parameters for layer 1 */
	const float vp1=6200.0, vs1=3600.0, rho1=2800.0;
	/*const float vp2=6200.0, vs2=3600.0, rho2=2800.0;*/
	/*const float vp2=7000.0, vs2=3900.0, rho2=2800.0;*/

	/*-----------------------------------------------------------------------*/
	sumu=0.0;
	sumpi=0.0;

	fprintf(FP," start creation of toy-model");

	if (L) {
		/* vector for maxwellbodies */
		pts=vector(1,L);

		for (l=1; l<=L; l++) {
			pts[l]=1.0/(2.0*PI*FL[l]);
			eta[l]=DT/pts[l];
		}

		ts=TAU;
		tp=TAU;
		ws=2.0*PI*FL[1];

		for (l=1; l<=L; l++) {
			sumu=sumu+((ws*ws*pts[l]*pts[l]*ts)/(1.0+ws*ws*pts[l]*pts[l]));
			sumpi=sumpi+((ws*ws*pts[l]*pts[l]*tp)/(1.0+ws*ws*pts[l]*pts[l]));
		}
	}

	/* loop over global grid */
	for (j=1; j<=NYG; j++) {
		for (i=1; i<=NXG; i++) {
			for (k=1; k<=NZG; k++) {


				Vp=vp1;
				Vs=vs1;
				Rho=rho1;

				muv=Vs*Vs*Rho/(1.0+sumu);
				piv=Vp*Vp*Rho/(1.0+sumpi);

				/* only the PE which belongs to the current global gridpoint
				  is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) &&
				        (POS[2]==((j-1)/NY)) &&
				        (POS[3]==((k-1)/NZ))) {
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;
					kk=k-POS[3]*NZ;

					if (L) {
						taus[jj][ii][kk]=ts;
						taup[jj][ii][kk]=tp;
					}

					u[jj][ii][kk]=muv;
					rho[jj][ii][kk]=Rho;
					pi[jj][ii][kk]=piv;
				}
			}
		}
	}

	if (kasten) {
		for (j=50; j<=65; j++) {	/*vertical*/
			for (i=50; i<=80; i++) {
				for (k=50; k<=90; k++) {

					Vp=0.0;
					Vs=0.0;
					Rho=0.0;
					Vp=vp1+500.0;
					Vs=vs1+300;
					Rho=rho1;

					muv=Vs*Vs*Rho/(1.0+sumu);
					piv=Vp*Vp*Rho/(1.0+sumpi);

					/* only the PE which belongs to the current global gridpoint
					  is saving model parameters in his local arrays */
					if ((POS[1]==((i-1)/NX)) &&
					        (POS[2]==((j-1)/NY)) &&
					        (POS[3]==((k-1)/NZ))) {
						ii=i-POS[1]*NX;
						jj=j-POS[2]*NY;
						kk=k-POS[3]*NZ;

						if (L) {
							taus[jj][ii][kk]=ts;
							taup[jj][ii][kk]=tp;
						}

						u[jj][ii][kk]=muv;
						rho[jj][ii][kk]=Rho;
						pi[jj][ii][kk]=piv;
					}
				}

			}
		}

		for (j=50; j<=95; j++) {
			for (i=81; i<=110; i++) {
				for (k=50; k<=90; k++) {

					Vp=0.0;
					Vs=0.0;
					Rho=0.0;
					Vp=vp1-250;
					Vs=vs1-100;
					Rho=rho1;

					muv=Vs*Vs*Rho/(1.0+sumu);
					piv=Vp*Vp*Rho/(1.0+sumpi);

					/* only the PE which belongs to the current global gridpoint
					  is saving model parameters in his local arrays */
					if ((POS[1]==((i-1)/NX)) &&
					        (POS[2]==((j-1)/NY)) &&
					        (POS[3]==((k-1)/NZ))) {
						ii=i-POS[1]*NX;
						jj=j-POS[2]*NY;
						kk=k-POS[3]*NZ;

						if (L) {
							taus[jj][ii][kk]=ts;
							taup[jj][ii][kk]=tp;
						}

						u[jj][ii][kk]=muv;
						rho[jj][ii][kk]=Rho;
						pi[jj][ii][kk]=piv;
					}
				}

			}
		}

		for (j=66; j<=95; j++) {
			for (i=50; i<=80; i++) {
				for (k=50; k<=90; k++) {

					Vp=0.0;
					Vs=0.0;
					Rho=0.0;
					Vp=vp1+300;
					Vs=vs1+200;
					Rho=rho1;

					muv=Vs*Vs*Rho/(1.0+sumu);
					piv=Vp*Vp*Rho/(1.0+sumpi);

					/* only the PE which belongs to the current global gridpoint
					  is saving model parameters in his local arrays */
					if ((POS[1]==((i-1)/NX)) &&
					        (POS[2]==((j-1)/NY)) &&
					        (POS[3]==((k-1)/NZ))) {
						ii=i-POS[1]*NX;
						jj=j-POS[2]*NY;
						kk=k-POS[3]*NZ;

						if (L) {
							taus[jj][ii][kk]=ts;
							taup[jj][ii][kk]=tp;
						}

						u[jj][ii][kk]=muv;
						rho[jj][ii][kk]=Rho;
						pi[jj][ii][kk]=piv;
					}
				}

			}
		}

		for (j=50; j<=95; j++) {
			for (i=50; i<=110; i++) {
				for (k=91; k<=110; k++) {

					Vp=0.0;
					Vs=0.0;
					Rho=0.0;
					Vp=vp1-400;
					Vs=vs1-150;
					Rho=rho1;

					muv=Vs*Vs*Rho/(1.0+sumu);
					piv=Vp*Vp*Rho/(1.0+sumpi);

					/* only the PE which belongs to the current global gridpoint
					  is saving model parameters in his local arrays */
					if ((POS[1]==((i-1)/NX)) &&
					        (POS[2]==((j-1)/NY)) &&
					        (POS[3]==((k-1)/NZ))) {
						ii=i-POS[1]*NX;
						jj=j-POS[2]*NY;
						kk=k-POS[3]*NZ;

						if (L) {
							taus[jj][ii][kk]=ts;
							taup[jj][ii][kk]=tp;
						}

						u[jj][ii][kk]=muv;
						rho[jj][ii][kk]=Rho;
						pi[jj][ii][kk]=piv;
					}
				}

			}
		}

	}


	MPI_Barrier(MPI_COMM_WORLD);

	if (L)	{
		free_vector(pts,1,L);
	}
}



