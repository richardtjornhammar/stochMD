#include <stdio.h>

#define XX     0
#define YY     1
#define ZZ     2
#define DIM    3
#define STRLEN 1000

typedef enum a_type {i,f,d,c} type;
typedef double vector[DIM];
typedef vector *t_vector;

typedef struct abstr{
  union {
       int i;
       float f;
       double d;
       char c;
       vector v;
    } dat;
  type t;
} abstract;

typedef abstract *t_abstract;

// RUN PARAMETERS
typedef struct param{
  int p;
  double NB[STRLEN][10];
} parameters;
typedef parameters *t_parameters;

// INTERACTION PARAMETERS NOT REALLY USED
typedef struct top{
  double **nb_mat;
  double **bon_mat;
  double **ang_mat;
} topology;
typedef topology *t_topology;

// CLUSTER & CLUSTERS
typedef struct clust{
  int    nQM;
  int    nCL;
  int    *qmndx;
  int    *clndx;
} cluster;
typedef cluster *t_cluster;

typedef struct clus{
  int       nQMc;
  int       nQMa;
  vector    *rQM;
  t_cluster *qmcluster;
  int       *indx_mat;
  double    *dist_mat;
  int       *bool_mat;
} clusters;
typedef clusters *t_clusters;

typedef struct molekyl{
  char   name[STRLEN];
  int    uid,pid;
  double *k_bond, *k_angle,*k_improp,*k_dihed;
  int    nr_t, nr_b, nr_a, nr_i, nr_d, bZ;
  int    *z_seq;   // atomic sequence (representation)
  int    *t_seq;   // type sequence
  int    *b_seq;   // bond sequence
  int    *a_seq;   // angl sequence
  int    *i_seq;   // impr sequence
  int    *d_seq;   // dihe sequence
} molecule;
typedef molecule *t_molecule;

// FRAME VALUES
typedef struct fram{
  double T,Tref,taut;
  double p,pref,taup,vpres,kpres,pcorr;
  double rc2;
  int    cmd,bGenv,bDebug,bSystem,bIR;
  double dt;
  int    NSTEPS;
  double E[10],virial,C6av;
  double mu[3];
  double deBroglie, deBroglie_vol, pc2;
  double box[DIM][DIM];
  int    nr_parts, nr_alloc,*nr_mols,nr_umols;
  t_molecule mol[STRLEN];
  int    mol_seq[STRLEN],n_start[STRLEN];
  vector *r;
  vector *ir;
  vector *v;
  vector *f;
  vector vcom,rcom;
  double W,fipimi,fifimi;
  int    *restyp;
  int    *partyp,*partyp_key,nr_unique;        //  Z 
  double *C6,*C12,*Q,*M;
  int    *resnr;
} frame;
typedef frame *t_frame;

typedef struct filer{
  FILE *fil[10];
  int  bOpen[10];
  char ftype[10];
  int nrfiles;
} files;
typedef files *t_files;

typedef struct quantum_interaction{
  int nrQMc;
} qmmm;

struct arguments{
  char *args[1];                        /* ARG1 and ARG2 */
  int  verbose,keep;                    /* The -v flag */
  char *outfile,*outfile_X,*outfile_E;
  char *infile,*infile2,*infile3;       /* Argument for -o and -i*/
  char *solitp;
  char *insert;                         /* Argument for -j */
  int  frames,start;
  int  cmd;
};

typedef struct arguments arguments;
typedef struct arguments *t_arguments;

