/*
 * det_init_trans.h
 *
 *  Created on: Mar 14, 2016
 *      Author: sushma_n
 */

#ifndef TIMENET_SOLUTION_DET_INIT_TRANS_H_
#define TIMENET_SOLUTION_DET_INIT_TRANS_H_

/*************************************************************************

  filename:	det_init_trans.h

  purpose:	Constants and types will be defined for init_trans.c

  author:       Christian Kelling
  		Technical University of Berlin
		Computer Science Department
		Institute for Technical Computer Science

  date:		May 1991, November 1992, Aug 93, Jan 1994

*************************************************************************/

/* Modes for function DetectInitialTransient */
#define STEADY_STATE 	0
#define TRANS_CROSSING	1	/* Test for k times crossing of mean value */
#define TRANS_TESTING	2	/* Test for stationary data */

/* Constant values for initializing parameters */
#define VARSEQLEN_DEF	100
#define CROSSES_DEF	25
#define SAFE_DEF_TRANS	2
#define EXCOEFF_DEF	0.5
#define ALPHAT_DEF	(1.0-(double)conf_level/100.0)
#define FREEDOM_DEF_TRANS	7

/* node of linked list to access a sequence of CVs */
typedef struct t_ControlSequenceList
  {
    SampSeq  *Seq_first;
    SampSeq  *Seq_last;
    SampSeq  *SeqPointer;
    SampSeq  *CurrentPosition;
    struct t_ControlSequenceList *Next;
  } ContrSeqList;


/* Environment used by the detection initial transient procedure for each opened STATISTICS */
typedef struct t_StatisticsEnvironment
	{
	bool		StopSimulation;
	int		StatStatus;
	SampSeq		*ReferenceSequence_first;
	SampSeq		*ReferenceSequence_last;
        double  	MeanValue;
        int             UpdatedYet;
        double          sum;
	double          InitStudent;
	double          TestStationary;
	int             InitCrossLength;
	int             DeltaLength;
	int             TestLength;
	int             SeqLength;
	int             InitTryings;
	int            	InitLength;
        long     	SampleCount;
        int     	SendCount;
        int     	ControlVariates;
	double		LastUpdateTime;
	double          TransTime;
	ContrSeqList    *CVList_first;
	sample_type	*SampSequence;
	str_sample_type *SampPhasis;
	double		*CVSequence;
      } StatInfo;





#endif /* TIMENET_SOLUTION_DET_INIT_TRANS_H_ */
