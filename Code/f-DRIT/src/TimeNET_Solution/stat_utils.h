/*
 * stat_utils.h
 *
 *  Created on: Mar 14, 2016
 *      Author: sushma_n
 */

#ifndef TIMENET_SOLUTION_STAT_UTILS_H_
#define TIMENET_SOLUTION_STAT_UTILS_H_


/*NEUE_VERSION (BM)*/
/*************************************************************************

  filename:	stat_utils.h

  purpose:	Constants and types will be defined for stat_utils.c

  author:	Christian Kelling
		Technical University of Berlin
		Computer Science Department
		Institute for Technical Computer Science

  date:		May 1993

  Datei im Ramen des PDV-Projekts im WS 95/96 ge�ndert.

  �nderung:       STAT_BLOCK_DEF und RECEIVE_BLOCK_DEF definiert.

  	          STAT_BLOCK_DEF defineirt die kleinste Blockgr��e im
		  Ma�speicher. Es ist sinnvoll den Wert 1 beizu behalten, da
		  dies die gr��te Performance-Verbesserung mit sich bringt.

		  RECEIVE_BLOCK_DEF definiert die maximale Empfangspuffergr��e.
		  Es ist zu beachten, da� die RECEIVE-Blockgr��e nicht kleiner
		  als die SEND-Blockgr��e sein darf.
		  Dieses MACRO ist nur der Verst�ndlichkeit wegen eingef�hrt
		  worden


  Autoren:	Mathias Neumann und Benjamin Hohnh�user

  Datum:   	April 1996

*************************************************************************/

#define  SUCCEED    0
#define  FAIL       1

#define EQU_MAX		100     /* max. number of successive equal samples */
				/* before assuming: no variance */

#define MIN_FIRINGS	50 	/* each transition should have been fired
				   MIN_FIRINGS times to ensure
				   valid statistics*/

#define LOCAL_BATCH_SIZE_DEF	20	/*-20- number of samples to batch at
					 *slave before record in sequence */

#define SEND_BLOCK_DEF	10	/*maximale Sendeblockgr��e (BM)*/

#define STAT_BLOCK_DEF	1	/* Blockgr��e im Ma�speicher (BM)
				 * Es ist sinnvoll den Wert 1 beizu-
				 * behalten, da dies die gr��te
				 * Performance-Verbesserung
				 * mit sich bringt*/

#define RECEIVE_BLOCK_DEF SEND_BLOCK_DEF   /* max. Empfangspuffergr��e (BM)
					    * Es ist zu beachten, da� die
					    * RECEIVE-Blockgr��e gleich der
					    * SEND-Blockgr��e sein mu�.
					    * Dieses MACRO ist nur der
					    * Verst�ndlichkeit wegen eingef�hrt
					    * worden */

#define MAX_SAMPLES_DEF 0x7FFFFFFF  /* 500000: a useful value,
				  0x7FFFFFFF: max possible value */

/* node of linked list to collect sampled data */
typedef struct t_SampleSequence
	{
	double Value;
	struct t_SampleSequence *Next;
	} SampSeq;


double means ();
 double normal_distr();
 double student_distr();
double periodogram ();
double SpectralVarAnalysis();

#define MAX_THRESHOLDS 30

/*#define DEBUG_RST */      /* If this activated you get some debugging Information */
                       /* Simulation with RESTART, see 'receive_samples()' below */






#endif /* TIMENET_SOLUTION_STAT_UTILS_H_ */
