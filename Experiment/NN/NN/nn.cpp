﻿/*******************************************************************************
*    nn.c   1.0                                       � JOHN BULLINARIA  2004  *
*******************************************************************************/

/*      To compile use "cc nn.c -O -lm -o nn" and then run using "./nn"       */
/*      For explanations see:  http://www.cs.bham.ac.uk/~jxb/NN/nn.html       */

#include <stdio.h>
#include <stdlib.h> 
#include <time.h>
#include <math.h>
#include <fcntl.h>
#include <windows.h>

#define NUMPAT 9
#define NUMIN  2
#define NUMHID 5
#define NUMOUT 1

#define rando() ((double)rand()/(RAND_MAX+1))

int    NumPattern = NUMPAT, NumInput = NUMIN, NumHidden = NUMHID, NumOutput = NUMOUT;
double WeightIH[NUMIN + 1][NUMHID + 1], WeightHO[NUMHID + 1][NUMOUT + 1];

void Compute(double *testInput)
{
	double testSumH[NUMHID + 1];
	double testSumO[NUMOUT + 1];
	double testHidden[NUMHID + 1];
	double testOutput[NUMOUT + 1];
	//double testInput[NUMIN + 1];
	int i, j, k;
	char str[256];

	//testInput = input;

	for (j = 1; j <= NumHidden; j++) {    /* compute hidden unit activations */
		testSumH[j] = WeightIH[0][j];
		for (i = 1; i <= NumInput; i++) {
			testSumH[j] += testInput[i] * WeightIH[i][j];
		}
		testHidden[j] = 1.0 / (1.0 + exp(-testSumH[j]));
	}
	for (k = 1; k <= NumOutput; k++) {    /* compute output unit activations and errors */
		testSumO[k] = WeightHO[0][k];
		for (j = 1; j <= NumHidden; j++) {
			testSumO[k] += testHidden[j] * WeightHO[j][k];
		}
		//testOutput[k] = 1.0 / (1.0 + exp(-testSumO[k]));   /* Sigmoidal Outputs */
		testOutput[k] = testSumO[k]; /* Linear Outputs */
		sprintf_s(str, "%f\t", testOutput[k]);
		OutputDebugString(str);   /* print network outputs */
	}
}

int main() {
	int    i, j, k, p, np, op, ranpat[NUMPAT + 1], epoch;
	char str[256];
	double Input[NUMPAT + 1][NUMIN + 1] = { 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 2, 1, 0, 1, 2, 0, 2, 2, 0, 4, 2, 0, 2, 4 };
	double Target[NUMPAT + 1][NUMOUT + 1] = { 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 3, 0, 3, 0, 0, 0, 6, 0, 6 };
	double SumH[NUMPAT + 1][NUMHID + 1], Hidden[NUMPAT + 1][NUMHID + 1];
	double SumO[NUMPAT + 1][NUMOUT + 1], Output[NUMPAT + 1][NUMOUT + 1];
	double DeltaO[NUMOUT + 1], SumDOW[NUMHID + 1], DeltaH[NUMHID + 1];
	double DeltaWeightIH[NUMIN + 1][NUMHID + 1], DeltaWeightHO[NUMHID + 1][NUMOUT + 1];
	double Error, eta = 0.05, alpha = 0.6, smallwt = 0.2;

	for (j = 1; j <= NumHidden; j++) {    /* initialize WeightIH and DeltaWeightIH */
		for (i = 0; i <= NumInput; i++) {
			DeltaWeightIH[i][j] = 0.0;
			WeightIH[i][j] = 2.0 * (rando() - 0.5) * smallwt;
		}
	}
	for (k = 1; k <= NumOutput; k++) {    /* initialize WeightHO and DeltaWeightHO */
		for (j = 0; j <= NumHidden; j++) {
			DeltaWeightHO[j][k] = 0.0;
			WeightHO[j][k] = 2.0 * (rando() - 0.5) * smallwt;
		}
	}

	for (epoch = 0; epoch < 1000000; epoch++) {    /* iterate weight updates */
		for (p = 1; p <= NumPattern; p++) {    /* randomize order of individuals */
			ranpat[p] = p;
		}
		for (p = 1; p <= NumPattern; p++) {
			np = p + rando() * (NumPattern + 1 - p);
			op = ranpat[p]; ranpat[p] = ranpat[np]; ranpat[np] = op;
		}
		Error = 0.0;
		for (np = 1; np <= NumPattern; np++) {    /* repeat for all the training patterns */
			p = ranpat[np];
			for (j = 1; j <= NumHidden; j++) {    /* compute hidden unit activations */
				SumH[p][j] = WeightIH[0][j];
				for (i = 1; i <= NumInput; i++) {
					SumH[p][j] += Input[p][i] * WeightIH[i][j];
				}
				Hidden[p][j] = 1.0 / (1.0 + exp(-SumH[p][j]));
			}
			for (k = 1; k <= NumOutput; k++) {    /* compute output unit activations and errors */
				SumO[p][k] = WeightHO[0][k];
				for (j = 1; j <= NumHidden; j++) {
					SumO[p][k] += Hidden[p][j] * WeightHO[j][k];
				}
				/* Output[p][k] = 1.0 / (1.0 + exp(-SumO[p][k])); */  /* Sigmoidal Outputs */
				Output[p][k] = SumO[p][k];   /*   Linear Outputs */
				Error += 0.5 * (Target[p][k] - Output[p][k]) * (Target[p][k] - Output[p][k]);   /* SSE */
				//Error -= ( Target[p][k] * log( Output[p][k] ) + ( 1.0 - Target[p][k] ) * log( 1.0 - Output[p][k] ) ) ; /*    Cross-Entropy Error */
				//DeltaO[k] = (Target[p][k] - Output[p][k]) * Output[p][k] * (1.0 - Output[p][k]);   /* Sigmoidal Outputs, SSE */
				/*              DeltaO[k] = Target[p][k] - Output[p][k];     Sigmoidal Outputs, Cross-Entropy Error */
				DeltaO[k] = Target[p][k] - Output[p][k];    /*  Linear Outputs, SSE */
			}
			for (j = 1; j <= NumHidden; j++) {    /* 'back-propagate' errors to hidden layer */
				SumDOW[j] = 0.0;
				for (k = 1; k <= NumOutput; k++) {
					SumDOW[j] += WeightHO[j][k] * DeltaO[k];
				}
				DeltaH[j] = SumDOW[j] * Hidden[p][j] * (1.0 - Hidden[p][j]);
			}
			for (j = 1; j <= NumHidden; j++) {     /* update weights WeightIH */
				DeltaWeightIH[0][j] = eta * DeltaH[j] + alpha * DeltaWeightIH[0][j];
				WeightIH[0][j] += DeltaWeightIH[0][j];
				for (i = 1; i <= NumInput; i++) {
					DeltaWeightIH[i][j] = eta * Input[p][i] * DeltaH[j] + alpha * DeltaWeightIH[i][j];
					WeightIH[i][j] += DeltaWeightIH[i][j];
				}
			}
			for (k = 1; k <= NumOutput; k++) {    /* update weights WeightHO */
				DeltaWeightHO[0][k] = eta * DeltaO[k] + alpha * DeltaWeightHO[0][k];
				WeightHO[0][k] += DeltaWeightHO[0][k];
				for (j = 1; j <= NumHidden; j++) {
					DeltaWeightHO[j][k] = eta * Hidden[p][j] * DeltaO[k] + alpha * DeltaWeightHO[j][k];
					WeightHO[j][k] += DeltaWeightHO[j][k];
				}
			}
		}
		if (epoch % 100 == 0) fprintf(stdout, "\nEpoch %-5d :   Error = %f", epoch, Error);
		if (Error < 0.00002) 
		{ 
			sprintf_s(str, "\nEpoch %-5d :   Error = %f", epoch, Error);
			OutputDebugString(str);   /* print network outputs */
			break;  /* stop learning when 'near enough' */
		}
	}
	
	sprintf_s(str, "\n\nNETWORK DATA - EPOCH %d\n\nPat\t", epoch);
	OutputDebugString(str);   /* print network outputs */
	for (i = 1; i <= NumInput; i++) {
		//printf("Input%-4d\t", i);
		sprintf_s(str, "Input%-4d\t", i);
		OutputDebugString(str);   /* print network outputs */
	}
	for (k = 1; k <= NumOutput; k++) {
		//printf("Target%-4d\tOutput%-4d\t", k, k);
		sprintf_s(str, "Target%-4d\tOutput%-4d\t", k, k);
		OutputDebugString(str);   /* print network outputs */
	}
	for (p = 1; p <= NumPattern; p++) {
		//fprintf(stdout, "\n%d\t", p);
		sprintf_s(str, "\n%d\t", p);
		OutputDebugString(str);   /* print network outputs */
		for (i = 1; i <= NumInput; i++) {
			//printf("%f\t", Input[p][i]);
			sprintf_s(str, "%f\t", Input[p][i]);
			OutputDebugString(str);   /* print network outputs */
		}
		for (k = 1; k <= NumOutput; k++) {
			//printf("%f\t%f\t", Target[p][k], Output[p][k]);
			sprintf_s(str, "%f\t%f\t", Target[p][k], Output[p][k]);
			OutputDebugString(str);   /* print network outputs */
		}
	}
	//printf("\n\nGoodbye!\n\n");
	sprintf_s(str, "\n\Network Trained!\n\n");
	OutputDebugString(str);   /* print network outputs */

	for (i = 1; i < NUMHID + 1; i++)
	{
		sprintf_s(str, "Bias %d: %f\n", i, WeightIH[0][i]);
		OutputDebugString(str);
	}
	for (i = 1; i < NUMIN + 1; i++)
	{
		for (j = 1; j < NUMHID + 1; j++)
		{
			sprintf_s(str, "Input %d: Hidden %d: Weight: %f\n", i, j, WeightIH[i][j]);
			OutputDebugString(str);
		}
	}

	for (i = 1; i < NUMOUT + 1; i++)
	{
		sprintf_s(str, "Bias %d: %f\n", i, WeightHO[0][i]);
		OutputDebugString(str);
	}
	for (i = 1; i < NUMOUT + 1; i++)
	{
		for (j = 1; j < NUMHID + 1; j++)
		{
			sprintf_s(str, "Output %d: Hidden %d: Weight: %f\n", i, j, WeightHO[i][j]);
			OutputDebugString(str);
		}
	}

	/* test */

	double testInput[NUMIN + 1];
	testInput[1] = 0;
	testInput[2] = 1;
	Compute(testInput);

	testInput[1] = 1;
	testInput[2] = 1;
	Compute(testInput);

	testInput[1] = 1;
	testInput[2] = 2;
	Compute(testInput);

	testInput[1] = 2;
	testInput[2] = 1;
	Compute(testInput);

	testInput[1] = 2;
	testInput[2] = 2;
	Compute(testInput);

	testInput[1] = 2;
	testInput[2] = 4;
	Compute(testInput);

	return 1;
}



/*******************************************************************************/