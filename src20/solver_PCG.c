/*
 * solver_PCG
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <omp.h>

#include "solver_PCG.h"
#include "pcg_ext.h"

extern int
solve_PCG (int N, int NL, int NU, int *indexL, int *itemL, int *indexU, int *itemU,
		double *D, double *B, double *X, double *AL, double *AU,
		int PEsmpTOT, int *SMPindexG,
  	        double EPS, int *ITR, int *IER, int N2)
{
	double **W;
	double VAL, BNRM2, WVAL, SW, RHO, BETA, RHO1, C1, DNRM2, ALPHA, ERR;
	double Stime, Etime;
	int i, j, ic, ip, L, ip1, N3;
	int R = 0;
	int Z = 1;
	int Q = 1;
	int P = 2;
	int DD = 3;

/*********
 * INIT. *
 *********/
        N3= N;
	W = (double **)malloc(sizeof(double *)*4);
	if(W == NULL) {
		fprintf(stderr, "Error: %s\n", strerror(errno));
		return -1;
	}
	for(i=0; i<4; i++) {
  	    W[i] = (double *)malloc(sizeof(double)*N3);
		if(W[i] == NULL) {
			fprintf(stderr, "Error: %s\n", strerror(errno));
			return -1;
		}
	}

if(NFLAG == 0){
	for(i=0; i<N; i++) {
		X[i] = 0.0;
		W[1][i] = 0.0;
		W[2][i] = 0.0;
		W[3][i] = 0.0;
	}
}else{
#pragma omp parallel for private (ip, i)
	for(ip=0; ip<PEsmpTOT; ip++) {
		for(i=SMPindexG[ip]; i<SMPindexG[ip+1]; i++) {
			X[i] = 0.0;
			W[1][i] = 0.0;
			W[2][i] = 0.0;
			W[3][i] = 0.0;
		}
	}
}
/**************************
 * {r0} = {b} - {A}{xini} *
 **************************/
#pragma omp parallel for private (ip, i, VAL, j)
	for(ip=0; ip<PEsmpTOT; ip++) {
		for(i=SMPindexG[ip]; i<SMPindexG[ip+1]; i++) {
			VAL = D[i] * X[i];
			for(j=indexL[i]; j<indexL[i+1]; j++) {
				VAL += AL[j] * X[itemL[j]-1];
			}
			for(j=indexU[i]; j<indexU[i+1]; j++) {
				VAL += AU[j] * X[itemU[j]-1];
			}
			W[R][i] = B[i] - VAL;
		}
	}

	BNRM2 = 0.0;
#pragma omp parallel for private (ip, i) reduction (+:BNRM2)
	for(ip=0; ip<PEsmpTOT; ip++) {
		for(i=SMPindexG[ip]; i<SMPindexG[ip+1]; i++) {
			BNRM2 += B[i]*B[i];
		}
	}

#pragma omp parallel for private(ip, i)
		for(ip=0; ip<PEsmpTOT; ip++) {
			for(i=SMPindexG[ip]; i<SMPindexG[ip+1]; i++) {
	  W[DD][i]= 1.e0/D[i];
	}
	}

/************************************************************** ITERATION */
	*ITR = N;

	Stime = omp_get_wtime();
	for(L=0; L<(*ITR); L++) {

/*******************
 * {z} = [Minv]{r} *
 *******************/
 #pragma omp parallel for private(ip, i)
		for(ip=0; ip<PEsmpTOT; ip++) {
			for(i=SMPindexG[ip]; i<SMPindexG[ip+1]; i++) {
				W[Z][i] = W[R][i]*W[DD][i];
			}
		}
/****************
 * RHO = {r}{z} *
 ****************/
		RHO = 0.0;
#pragma omp parallel for private (ip, i) reduction (+:RHO)
		for(ip=0; ip<PEsmpTOT; ip++) {
			for(i=SMPindexG[ip]; i<SMPindexG[ip+1]; i++) {
				RHO += W[R][i] * W[Z][i];
			}
		}

/********************************
 * {p}  = {z} if      ITER=0    *
 * BETA = RHO / RHO1  otherwise *
 ********************************/
		if(L == 0) {
#pragma omp parallel for private (ip, i)
			for(ip=0; ip<PEsmpTOT; ip++) {
				for(i=SMPindexG[ip]; i<SMPindexG[ip+1]; i++) {
					W[P][i] = W[Z][i];
				}
			}
		} else {
			BETA = RHO / RHO1;
#pragma omp parallel for private (ip, i)
			for(ip=0; ip<PEsmpTOT; ip++) {
				for(i=SMPindexG[ip]; i<SMPindexG[ip+1]; i++) {
					W[P][i] = W[Z][i] + BETA * W[P][i];
				}
			}
		}

/****************
 * {q} = [A]{p} *
 ****************/
#pragma omp parallel for private (ip, i, VAL, j)
		for(ip=0; ip<PEsmpTOT; ip++) {
			for(i=SMPindexG[ip]; i<SMPindexG[ip+1]; i++) {
				VAL = D[i] * W[P][i];
				for(j=indexL[i]; j<indexL[i+1]; j++) {
					VAL += AL[j] * W[P][itemL[j]-1];
				}
				for(j=indexU[i]; j<indexU[i+1]; j++) {
					VAL += AU[j] * W[P][itemU[j]-1];
				}
				W[Q][i] = VAL;
			}
		}

/************************
 * ALPHA = RHO / {p}{q} *
 ************************/
		C1 = 0.0;
#pragma omp parallel for private (ip, i) reduction (+:C1)
		for(ip=0; ip<PEsmpTOT; ip++) {
			for(i=SMPindexG[ip]; i<SMPindexG[ip+1]; i++) {
				C1 += W[P][i] * W[Q][i];
			}
		}

		ALPHA = RHO / C1;

/***************************
 * {x} = {x} + ALPHA * {p} *
 * {r} = {r} - ALPHA * {q} *
 ***************************/
#pragma omp parallel for private (ip, i)
		for(ip=0; ip<PEsmpTOT; ip++) {
			for(i=SMPindexG[ip]; i<SMPindexG[ip+1]; i++) {
				X[i]    += ALPHA * W[P][i];
			}
		}

#pragma omp parallel for private (ip, i)
		for(ip=0; ip<PEsmpTOT; ip++) {
			for(i=SMPindexG[ip]; i<SMPindexG[ip+1]; i++) {
				W[R][i] -= ALPHA * W[Q][i];
			}
		}

		DNRM2 = 0.0;
#pragma omp parallel for private (ip, i) reduction (+:DNRM2)
		for(ip=0; ip<PEsmpTOT; ip++) {
			for(i=SMPindexG[ip]; i<SMPindexG[ip+1]; i++) {
			  DNRM2 += W[R][i]*W[R][i];
			}
		}

		ERR = sqrt(DNRM2/BNRM2);
                if( (L+1)%100 ==1) {
                        fprintf(stdout, "%5d%16.6e\n", L+1, ERR);
                }

		if(ERR < EPS) {
			*IER = 0;
			goto N900;
		} else {
			RHO1 = RHO;
		}
	}
	*IER = 1;

N900:
	Etime = omp_get_wtime();
	fprintf(stdout, "%5d%16.6e\n", L+1, ERR);
	fprintf(stdout, "%16.6e sec. (solver)\n", Etime - Stime);
	*ITR = L;

	free(W);

	return 0;
}
