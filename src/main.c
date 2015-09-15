/*
*	THIS MOLECULAR DYNAMICS PACKAGE SUPPORT GRAND CANONICAL SIMULATIONS
*	OF SYSTEM WITH BONDED AND NON-BONDED INTERACTIONS. THE VELOCITY RE-
*	SCALING FOR THE TEMPERATURE COUPLING IS STOCHASTIC.
*       CODE WAS WRITTEN BY: "Richard Tjoernhammar"
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <argp.h>
#include <math.h>
#include <time.h>
#include "define.h"
#include "structs.h"
#include "parse.h"
#include "clean.h"
#include "iofunctions.h"
#include "force.h"
#include "thermo.h"

int main (int argc, char **argv)
{
  //GENERAL
  int          i,j,k,nr;
  arguments    args;
  t_arguments  pargs;
  files        inputfiler;
  t_files      filepkg;
  frame        frm;
  t_frame      pframe;
  parameters   params;
  t_parameters pparams;
  topology     top;
  t_topology   ptop;
  clusters     clusts;
  t_clusters   pclusters;

  //TIMING 
  time_t tid[11];
  double sec[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double stp[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  char   what[10][100]={"null","TEMP","CHECK","QM","NB","BON","CHECK","THERMO","IO","null"};
  int    itid;

  // ASSIGN
  pclusters = &clusts;
  pargs     = &args;
  filepkg   = &inputfiler;
  pframe    = &frm;
  pparams   = &params;
  ptop      = &top;

  tid[0]    = time(NULL);

  // PARSE INPUT
  parse(argc,argv,pargs,filepkg);

  // GET RUN DATA
  read_params(pframe,pargs,pparams,filepkg);
  init_particles(pframe,pparams,filepkg);
  pframe->mu[0]=1.0;  pframe->mu[1]=1.0;

  if(pframe->bGenv)
    gen_velocities(pframe);

  if(pframe->cmd==cmdMD){
    fprintf(filepkg->fil[0],"GOING TO DO MOLECULAR DYNAMICS\n");
  }
  if(pframe->cmd==cmdMC){
    fprintf(filepkg->fil[0],"GOING TO DO GRAND CANONICAL SIMULATION\n");
  }

  if(pframe->bDebug)
    fprintf(filepkg->fil[0],"DEBUG IS TURNED ON\n");

  // COMMAND LOOP
  for(i=0;i<pframe->NSTEPS;i++){
    tid[1] = time(NULL);
    
    if( pframe->cmd == cmdMC && !fmod(i,1)){
      compute_temperature(pframe,filepkg);
      if(!fmod(i,10))
        grand_canonical(pframe,filepkg);
    }

    tid[2] = time(NULL);
    if(pframe->bDebug)
      check_frame(pframe,i+1);
    tid[3] = time(NULL);
    // CLEAN 

    // BEGIN FORCE
    tid[4] = time(NULL);
    force_calc_nb(pframe,filepkg);    // NONBONDED
    tid[5] = time(NULL);
    force_calc_bon(pframe,filepkg);   // BONDED
    tid[6] = time(NULL);
    // END FORCE

    if(pframe->bDebug)
      check_frame(pframe,(i+1)*(-1));
    tid[7] = time(NULL);
    // CLEAN

    // BEGIN THERMODYNAMICS
    compute_temperature(pframe,filepkg);
    compute_pressure(pframe,filepkg);
    propagate_thermostat(pframe);
    propagate_barostat(pframe);
    update_conf(pframe);
    // END THERMO

    tid[8] = time(NULL);
    if(pframe->bDebug)
      check_frame(pframe,(i+1)*(-1));

    write_data(pframe,filepkg,i*(pframe->dt));

    tid[9] = time(NULL);
    for(itid=1;itid<9;itid++){
      sec[itid]+=difftime(tid[itid+1],tid[itid]);
      stp[itid]+=1.0;
    }
  }
  // FINISH COMMAND

  // OUTPUT RESULTS
  print_coord(pframe,filepkg);
  
  tid[10] = time(NULL);

  // CLEAN UP
  dealloc_data(pframe,pclusters,0);

  fprintf(filepkg->fil[0],"EXECUTION END\n");
  
  close_streams(filepkg);

  fprintf(stdout,"\n==========================\nTIMING::\n");
  for(itid=1;itid<9;itid++){
    fprintf(stdout,"%10s TOOK %10lf \n",what[itid],sec[itid]);///stp[itid]
  }
  fprintf(stdout,"%10s TOOK %10lf \n==========================\n\n","TOTAL",difftime(tid[10],tid[0]));
  fflush(stdout);

  return 0;
}
