#ifndef __H_PCG
#define __H_PCG

	static int N2 = 128;
	int NUmax, NLmax, NCOLORtot, NCOLORk, NU, NL;
	int METHOD, ORDER_METHOD, NFLAG;

	double EPSICCG;

	double *D, *PHI, *BFORCE;
	//double **AL, **AU;
	double *AL, *AU;

	int *INL, *INU, *COLORindex;
	int *indexL, *indexU;
	int *indexLorg, *indexUorg;
	int *SMPindex, *SMPindexG;
	int *OLDtoNEW, *NEWtoOLD;

	int **IAL, **IAU;
	int *itemL, *itemU;
	int *itemLorg, *itemUorg;

	int NPL, NPU;

#endif /* __H_PCG */
