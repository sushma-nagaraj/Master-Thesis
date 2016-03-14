/*
 * sim_slave_types.h
 *
 *  Created on: Mar 14, 2016
 *      Author: sushma_n
 */

#ifndef TIMENET_SOLUTION_SIM_SLAVE_TYPES_H_
#define TIMENET_SOLUTION_SIM_SLAVE_TYPES_H_

/******************************************************************************

  filename :	sim_slave_types.h

  purpose  :	type definitions for DSPN simulation slave

  author   :    Christian Kelling
  		Technical University Berlin
		Computer Science Department

  date     :	June 1993


******************************************************************************/
#ifndef __SIM_SLAVE_TYPES
#define __SIM_SLAVE_TYPES

#define false 0
#define true 1

#define  SUCCEED    0
#define  FAIL       1

#define MASTER_STOP 0
#define MASTER_RDY  1
#define MEASURE_STOP 2

#define MAX_TRIES_ROOT_FINDING  50
#define ACCURACY              1.0e-10

#define MAX_TOKEN	100000  /* max. number of tokens in a place */
#define MAX_CVS		10	/* max. number of control variates */
#define NO_TOK		0
#define ONE_TOK		1
#define FILE_NAME	1000    /* max. length of a filename */
#define	NEXT_LINE_BUF	500

/* constants for measure_flag in places */
#define NO_MEASURE	0      /* no measure to estimate for this place */

/* constant for function check_enabled*/
#define NO_TRANS	-1	/* no transition is enabled -> dead */

/* constants for transition types */
#define IMM             0
#define EXP		1
#define DET		2
#define GEN		3

/* constants for firing policies */
#define RS             0
#define RA		1
#define RE		2

#define INF_SERV	0
#define SING_SERV	1

/* constants for arc types */
#define INP_ARC		0
#define OUTP_ARC	1
#define INH_ARC		2
#define INP_OR_INH	3

/* constants for other connectivities */
#define ENA_CONN	4  /* no arc connection but influence via enabling
			      function*/
#define MDA_CONN	5  /* no arc connection but influence via m.-d.
			      arc multiplicity */
/* constants for selection of measure output in _ERC.c */
#define SEL_MEAN      	0
#define SEL_UP		1
#define SEL_LO		2

/* constants for File-Interface (generell Transition) */
#define SEQUENCE 0
#define SEQU_REV 1
#define RAND_UNI 2


//#ifdef RST
#define RESTART_BATCH_SIZE_DEF 20           /* 25.02.1997 */
/*#define RESTART_BATCH_SIZE_DEF 100 */     /* alter Wert */
#define MAX_RESETTING_OF_THRESHOLDS 20
//#endif

/* type definitions */

typedef struct sample_type	/* type for a sample containing current time */
{
  double value;		/* the value itself */
  double s_time;	/* the sampling time of the value or the batching
   			   time if batches are used */
} sample_type;

typedef double str_sample_type;	/* type for a sample  */

typedef struct par_name		/* type for marking or delay parameter */
{
  char	name[30];
  float value;
} par_name;

typedef struct conn_list	/* type for list of connections */
{				/* of a transition*/
  int	multiplicity;
  int	element;
  int	id;
} conn_list;

typedef struct cv_list_elem
{
  int meas_number;
  int cv_number;
  int recorded;
} cv_list_elem;

typedef struct enabled_trans	/* type for list of enabled transitions, */
{				/* both immediate and timed */
  int	number;
  double sample;
  double elapsed_time;
  struct enabled_trans *next;
} enabled_trans;

typedef struct place		/* type for all places */
{
  char	name[30];
  int	token;
  int	no_measures;  /* denotes the number of measures this place */
  		       /* is involved in */
  int	*meas_list;
  		       /* a list of identifiers used as arguments when
		       calling the select_measure-functions from <net>_ERC.c*/
} place;
/** shorter version of Place needed for in sim trans */
typedef struct hiplace
{

	int tokens;
} hiplace;
/** for every hiPlace arrays there is the history info */
typedef struct history
{
	int passed; //0 not passed 1 passed
	double time;
} history;

typedef struct transition	/* type for all transitions */
{
  char	name[30];
  int	type;
  int  serv_type;
  double  delay;
  double  sample;
  int	prio;
  int 	fire_policy;	/* firing policy */
  int	ena_id;		/* identifikation of enabling function */
  int	mark_dep_id; 	/* identification
			   of marking dependent delay or weight */
  int	ok;             /* flag checks wether the m-d delay or weight could
			   determined by parsing .DEF-file */
  int noEntries;        /* for File-Interface */
  int lastEntry;        /* for File-Interface */
  double *buffer;       /* size dependended from noEntries => look at random.c */
  long	no_firings;
  int no_outputs;
  int no_inputs;
  int no_inh;
  int no_aff;		/* number of affected transitions */
  int no_cv_meas;	/* number of measures the transitions is involved in
			   as control variate*/
  int just_enabled;
  conn_list  *out_to_place;
  conn_list  *in_from_place;
  conn_list  *inh_from_place;
  int no_in_list;
#ifdef WIN32
  char affected;
#else
  bool affected;
#endif
  int *affected_trans;
  void (*add_ttf)();
  void (*change_ttf)();
  void (*del_ttf)();
  cv_list_elem **cv_list;
  double samp_rnd;
  double (*get_sample)();
} transition;

// Datatype for Breakpoint-Management of Tokengame
typedef struct breakpoint {
	char *name;
	char relation;
	int value;
} breakpoint;

/******************************* sim_list ****************************/

#define E_IMM 0
#define E_TIM 1
#define E_RAG 2


#define NO_SAVESTATES 30

/* ttf: eventuell feuernde Transition */

typedef struct ttf
{ struct ttf *prev;
  struct ttf *next;
  struct ttf *t_prev;
  struct ttf *t_next;
  int    tra_nr;
  int    ttf_type;
  int    prio;
  double samp_rnd;
  double sample;
  double rem_delay;
  double time_in_list;
} ttf;


/* exported functions */
void init_ttf();

void add_imm_ttf();
void add_tim_ttf();
void add_rag_ttf();

void del_ttf();
void del_fir_tra();
void off_rag_ttf();

void start_changing();
double get_samp_rnd();
void change_tim_ttf();
void replace_RS_ttf();
void replace_imm_ttf();

int trans_to_fire();

double get_CVsample();

void save_list();
void load_list();
#endif





#endif /* TIMENET_SOLUTION_SIM_SLAVE_TYPES_H_ */
