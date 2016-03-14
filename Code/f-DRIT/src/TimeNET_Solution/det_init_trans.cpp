/*
 * det_init_trans.cpp
 *
 *  Created on: Mar 14, 2016
 *      Author: sushma_n
 */

/*
 * det_init_trans.c
 *
 *  Created on: Mar 14, 2016
 *      Author: sushma_n
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <iostream>
#include <cstddef>

#include "num_utils.h"
#include "stat_utils.h"
#include "sim_slave_types.h"
#include "det_init_trans.h"
#include "stat_utils.c"
#include "num_utils.c"

/* global variables used by the modul */
static int TRANS_n_max = 0;
static int TRANS_n0_max = 0;
static int paw_crosses = 0;

static SampSeq	*TempDetSeq	= NULL;

void mylog(const char* format, ...)
{
	FILE *fff = fopen("C:\\log.txt", "a");
	va_list argptr;
	va_start(argptr, format);
	vfprintf(fff, format, argptr);
	fprintf(fff, "\n");
	va_end(argptr);
	fclose(fff);
}


/*************************************************************************

  function:	DetInitTrans_data_out

  parameters:	MeasureInfo -> pointer to a reward measure statistic

  purpose:	performs the output of results of the detection of the initial
		transient

  author:	Christian Kelling
		Technical University of Berlin

  date:		August 1993, Janauar 1994

 *************************************************************************/
void DetInitTrans_data_out(StatInfo *MeasureInfo)
{
	printf("1");
	printf("\nDetection Initial Transient Phase:\n");
	printf("==================================\n\n");
	printf("TransSamples: %d\t",MeasureInfo->InitLength);
	printf("TransTime: %f\t",MeasureInfo->TransTime);
	printf("Tryings: %d\n",MeasureInfo->InitTryings);
	printf("StationaryTest: %f\t",MeasureInfo->TestStationary);
	printf("TestLength: %d\t",MeasureInfo->TestLength );
	printf("t_student: %f\n",MeasureInfo->InitStudent);
} /* DetInitTrans_data_out */

/*************************************************************************

  function:	DetInitTrans_Init

  parameters:	maxsamples -> the maximum allowed length of the simulation
			run measured in the number of recorded observations.
			The calling function uses the number of messages
			to transmit.

  purpose:	init local data of this modul and global variables used
                by all statistics

  author:	Christian Kelling
		Technical University of Berlin

  date:		Nov 1992, January 1994

 *************************************************************************/
void DetInitTrans_Init (long maxsamples)
{
	TRANS_n_max = maxsamples;
	TRANS_n0_max = TRANS_n_max / 5;
	paw_crosses = CROSSES_DEF;
	printf("2");
} /* DetInitTrans_Init */

/*************************************************************************

  function:	Init_One_DetInitTrans

  purpose:	init local data used by any single statistic
                of this modul to perform determining the
		start point of steady-state simulation.

 return:   	MeasureInfo -> a pointer to a measure info structure

 author:	Christian Kelling
		Technical University of Berlin

  date:		May 1991, Nov 1992, January 1994

 *************************************************************************/
char *er_malloc (int bytes);

StatInfo *Init_One_DetInitTrans ()
{
	StatInfo *MeasureInfo;
	int i;
	MeasureInfo = (StatInfo *) er_malloc (sizeof (StatInfo));

	MeasureInfo->StopSimulation = false; /* Renate, 3.7.97 */
	MeasureInfo->InitTryings = 0;
	MeasureInfo->InitCrossLength = 0;
	MeasureInfo->InitLength = 0;
	MeasureInfo->SeqLength = 0;
	MeasureInfo->SampleCount = 0;
	MeasureInfo->SendCount= 0;
	MeasureInfo->sum = 0.0;
	MeasureInfo->LastUpdateTime = 0.0;
	MeasureInfo->TransTime = 0.0;
	MeasureInfo->StatStatus = TRANS_CROSSING;
	MeasureInfo->ReferenceSequence_first = NULL;
	MeasureInfo->ReferenceSequence_last = NULL;
	MeasureInfo->CVList_first = NULL;
	MeasureInfo->SampSequence =
			(sample_type *) er_malloc(SEND_BLOCK_DEF*sizeof(sample_type));
	MeasureInfo->CVSequence =
			(double *) er_malloc(MAX_CVS*SEND_BLOCK_DEF*sizeof(double));
	//	    StatSema ++;

	for(i = 0; i < SEND_BLOCK_DEF; i++)
	{
		MeasureInfo->SampSequence[i].value = 0.0;
		MeasureInfo->SampSequence[i].s_time = 0.0;
	}
	for(i = 0; i < MAX_CVS*SEND_BLOCK_DEF; i++)
	{
		MeasureInfo->CVSequence[i] = 0.0;
	}


	return (MeasureInfo);
}



/*************************************************************************

  function:	DetectInitialTransient

  parameters:	MeasureInfo   -> pointer to a reward measure statistic
		sample-> sampled item of a random variable

  returns:	true  -> the point of steady-state is passed
		false -> the initial transient period isn't finished yet.

  purpose:	determines the length of the initial transient period.
		The first approximation of the number of observations to be
		deleted	is evaluated by the heuristic rule R5
		([Pawlikowsky]: p.145).
		After that first approximation the function uses the
		stationary test	proposed by Schruben et al. [1983].
		This function works with a input data stream asynchronously.
		It collects the single sample item and performs tests as soon
		as sufficient sampled data are collected.

  author:	Christian Kelling
		Technical University of Berlin

  date:		August 93, Jan 1994

 *************************************************************************/

bool DetectInitialTransient (StatInfo *MeasureInfo, double sample, bool trans_detection, double current_time, int	conf_level)
{
	printf("4");
	static int varseqlen = VARSEQLEN_DEF;
	static int safecoeff = SAFE_DEF_TRANS;
	static double excoeff = EXCOEFF_DEF;
	static double alphat1, alphat;

	bool oldflag = false, newflag = false;
	int crosscount = 0, i = 0, free = 0, tryings = 0;
	double variance = 0.0, teststat = 0.0, sum = 0.0, temp_sum = 0.0;
	double *testsamplesarray = NULL, nt_mean = 0.0;

	alphat = (1.0-(double)conf_level/100.0); //replaced def for ALPHAT_DEF

	if(!(*MeasureInfo == null))
	{
		mylog("Measure Info");
	}

	if (!trans_detection)
	{
		MeasureInfo->InitLength = 0;//0.0;
		MeasureInfo->TransTime = current_time;
		MeasureInfo->StatStatus = STEADY_STATE;
		MeasureInfo->SampleCount = 1;
		MeasureInfo->sum = 0.0;

		mylog("No transient");

		return (true);
	}
	switch (MeasureInfo->StatStatus)
	{
	case TRANS_CROSSING:
		/* STEP 1 */
		TempDetSeq = (SampSeq*)er_malloc(sizeof(SampSeq));
		TempDetSeq->Value = sample;
		if (MeasureInfo->ReferenceSequence_first == NULL)
		{
			TempDetSeq->Next = NULL;
			MeasureInfo->ReferenceSequence_first = TempDetSeq;
		}
		else
			TempDetSeq->Next = MeasureInfo->ReferenceSequence_first;

		MeasureInfo->ReferenceSequence_first = TempDetSeq;
				MeasureInfo->sum += sample;
				MeasureInfo->SampleCount++;
				if (MeasureInfo->SampleCount >= paw_crosses)
				{
					if (MeasureInfo->SampleCount > TRANS_n0_max)
					{
						printf("Initial transient period too long\n");
						printf("or measure with (approx.) deterministic behaviour\n");
						printf("Samples tested: %d\n",TRANS_n0_max );
						exit(FAIL);
					}
					MeasureInfo->MeanValue = MeasureInfo->sum / MeasureInfo->SampleCount;
								crosscount = 0;
								tryings = 0;
								TempDetSeq = MeasureInfo->ReferenceSequence_first;
								oldflag = (MeasureInfo->MeanValue > TempDetSeq->Value);
								while (TempDetSeq->Next != NULL)
								{
									newflag = (MeasureInfo->MeanValue > TempDetSeq->Value);
									if (newflag != oldflag)
										crosscount++;
									oldflag = newflag;
									TempDetSeq = TempDetSeq->Next;
									tryings ++;
								}
								if (crosscount >= paw_crosses)
											{
												MeasureInfo->InitCrossLength = MeasureInfo->SampleCount;
												if ((safecoeff * varseqlen) >=
														(excoeff * MeasureInfo->SampleCount))
													MeasureInfo->TestLength = safecoeff * varseqlen;
												else
													MeasureInfo->TestLength = (int)(excoeff * MeasureInfo->SampleCount);
												MeasureInfo->DeltaLength = MeasureInfo->TestLength;
												if (MeasureInfo->InitCrossLength + MeasureInfo->TestLength
														> TRANS_n0_max)
												{
													printf("Initial transient period too long\n");
													printf("TRANS n0 max: %d\n",TRANS_n0_max );
													printf("CrossL + TestL: %d\n",
															MeasureInfo->InitCrossLength+
															MeasureInfo->TestLength);
													exit(FAIL);
												}
												MeasureInfo->StatStatus = TRANS_TESTING;
											}
											if ((tryings >= EQU_MAX) && (crosscount < 3))
											{
												printf("NOTE: A statistic seems to have no variance,\n");
												printf("%d successive equal samples (value: %f) recorded.\n",
														tryings, sample);
												printf
												("Check the result definitions for deterministic behaviour.\n");
												MeasureInfo->InitLength = 0;
												MeasureInfo->TransTime = current_time;
												MeasureInfo->StatStatus = STEADY_STATE;
												MeasureInfo->SampleCount = 1;
												MeasureInfo->sum = 0.0;
												return(true);
											}
										}
										break;

	case TRANS_TESTING:
			/* STEP 2 */
			MeasureInfo->SeqLength++;
			TempDetSeq = (SampSeq*)er_malloc(sizeof(SampSeq));
			TempDetSeq->Value = sample;
			TempDetSeq->Next = MeasureInfo->ReferenceSequence_first;
			MeasureInfo->ReferenceSequence_first = TempDetSeq;
			if (MeasureInfo->SeqLength >= MeasureInfo->DeltaLength)
			{
				MeasureInfo->InitTryings++;
				/* STEP 3 */
				testsamplesarray = dvector(1,varseqlen);
				for (i=varseqlen; i > 0; i--)
				{
					testsamplesarray[i] = TempDetSeq->Value;
					TempDetSeq = TempDetSeq->Next;
				}
				variance = SpectralVarAnalysis(testsamplesarray, varseqlen);
				/* STEP 4 */
				TempDetSeq = MeasureInfo->ReferenceSequence_first;
				nt_mean = means(MeasureInfo->TestLength, TempDetSeq);
				sum = 0.0;for (i=MeasureInfo->TestLength; i>0; i--)
				{
					temp_sum = nt_mean - means(i,TempDetSeq);
					temp_sum *=
						(double)i * (1.0 - (double)i / MeasureInfo->TestLength);
					sum += temp_sum;
					TempDetSeq = TempDetSeq->Next;
				}
				teststat = sum / sqrt(variance);
				teststat *= (sqrt(45.0) / sqrt((double)varseqlen));
				teststat /= pow((double)MeasureInfo->TestLength,1.5);
				/* use two sided test because sign of initial bias is unknown */
				teststat = fabs(teststat);
				alphat1 = alphat / 2.0;
				free = FREEDOM_DEF_TRANS;
				/* STEP 5 */
				free_dvector(testsamplesarray,1,varseqlen);
				if (teststat <= student_distr(free,1.0-alphat1))
				{
					/* initial transient period is determined */
					/* initial transient period = n0 + nt !? */
					MeasureInfo->TestStationary = teststat;
					MeasureInfo->InitStudent = student_distr(free,1.0-alphat1);
					MeasureInfo->InitLength =
						MeasureInfo->InitCrossLength + MeasureInfo->TestLength;
					MeasureInfo->TransTime = current_time;
					MeasureInfo->StatStatus = STEADY_STATE;
					MeasureInfo->SampleCount = 1;
					MeasureInfo->sum = 0.0;
					/*	      DetInitTrans_data_out (MeasureInfo);*/
					return(true);
				}
				else
				{
					/* T > t(k,1-alpha_t1) -> no steady-state ! */
					/* n0 und delta_n : see Pawlikowsky ! */
					MeasureInfo->InitLength +=
						(int)(excoeff * MeasureInfo->SampleCount);
					MeasureInfo->SeqLength = 0;
					if (MeasureInfo->InitLength + MeasureInfo->TestLength
							> TRANS_n0_max)
					{
						printf("Initial transient period too long\n");
						printf("TRANS n0 max: %d\n",TRANS_n0_max );
						printf("InitL + TestL: %d\n",
								MeasureInfo->InitLength+
								MeasureInfo->TestLength);
						exit(FAIL);
					}
				}
			}
			break;

		default:
			printf("DetectInitialTransient: unknown Mode\n");
			exit (FAIL);
			break;
		}
		return(false);
	} /* DetectInitialTransient */



