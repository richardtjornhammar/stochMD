#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <argp.h>
#include <math.h>
#include "define.h"
#include "structs.h"
#include "parse.h"

const char *argp_program_version =
"stochMD 1.0";

const char *argp_program_bug_address =
"<richard.tjornhammar@gmail.com>";

void fatal(char errstr[])
{
  fprintf(stderr,"ERROR::  %s\n",errstr);
  exit(0);
}

int file_ext(char filenm[], char filext[])
{
  char* ptr;
  
  ptr=strpbrk(filenm,".");
  ptr;
  if(ptr==NULL)
    return(-1);
  else{
    ptr++;
    return(strcmp(ptr,filext));
  }
}

/*
   OPTIONS.  Field 1 in ARGP.
   Order of fields: {NAME, KEY, ARG, FLAGS, DOC}.
*/
static struct argp_option options[] =
{
  {"verbose", 'v', 0, 0, "Produce verbose output"},
  {"keep", 'k', 0, 0, "Store the positions of the first molecule in the first frame for grand canonical simulations."},
  {"coord",  'c', "COORDFILE", 0,
   "Coordinate information in coordfile"},
  {"topol",  'p', "TOPFILE", 0,
   "System information from topology"},
  {"input",  'i', "INPUTFILE", 0,
   "System information from input"},
  {"output",  'o', "OUTFILE", 0,
   "Output to OUTFILE instead of to standard output"},
  {"outx",  'x', "OUTX", 0,
   "Writes output coordinates to outx"},
  {"oute",  'e', "OUTE", 0,
   "Writes output energy to oute"},
  {0}
};
/*
   PARSER. Field 2 in ARGP.
   Order of parameters: KEY, ARG, STATE.
*/
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  t_arguments pargs = state->input;

  switch (key)
    {
    case 'v':
      pargs->verbose = 1;
      break;
    case 'k':
      pargs->keep = 1;
      break;
    case 'o':
      pargs->outfile = arg;
      break;
    case 'i':
      pargs->infile3 = arg;
      break;
    case 'c':
      pargs->infile  = arg;
      break;
    case 'p':
      pargs->infile2 = arg;
      break;
    case ARGP_KEY_ARG:
      if (state->arg_num >= 1){
	  argp_usage(state);
      }
      pargs->args[state->arg_num] = arg;
      break;
    case ARGP_KEY_END:
      if (state->arg_num < 1)
	{
	  argp_usage (state);
	}
      break;
    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

/*
   ARGS_DOC. Field 3 in ARGP.
   A description of the non-option command-line arguments
     that we accept.
*/
static char args_doc[] = "COMMAND";//"ARG1 ARG2";

/*
  DOC.  Field 4 in ARGP.
  Program documentation.
*/
static char doc[] =
"\n\nstochMD (C) 2015 Richard Tjornhammar is a small molecular dynamics program written by Richard Tjoernhammar.\n\n\tThis program is free software: you can redistribute it and/or modify\n\tit under the terms of the GNU General Public License as published by\n\tthe Free Software Foundation, either version 3 of the License, or\n\t(at your option) any later version.\n\n\tThis program is distributed in the hope that it will be useful,\n\tbut WITHOUT ANY WARRANTY; without even the implied warranty of\n\tMERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\tGNU General Public License for more details.\n\tYou should have received a copy of the GNU General Public License\n\talong with this program.  If not, see <http://www.gnu.org/licenses/> ";
/*
   The ARGP structure itself.
*/
static struct argp argp = {options, parse_opt, args_doc, doc};

/*
   The main function.
   Notice how now the only function call needed to process
   all command-line options and arguments nicely
   is argp_parse.
*/
int parse (int argc, char **argv, t_arguments pargs, t_files fpkg)
{
  //GENERAL
  int i,j,k,nr;
  arguments args;

  //IO
  int typ;
  char *exten[]={"GRO","PDB","XYZ","TOP","DAT","INP"};
  char *commands[]={"MD","CE","MC"};
  int  NRCMDS=3,icmd;
  char *cmd;

  // DEFAULTS 
  pargs->outfile = NULL;
  pargs->outfile_E = NULL;
  pargs->outfile_X = NULL;
  pargs->infile  = NULL;
  pargs->infile2 = NULL;
  pargs->infile3 = NULL;
  pargs->solitp  = NULL;
  pargs->verbose = 0;
  pargs->keep    = 0;
  pargs->frames  = 1;

  // PARSER 
  argp_parse (&argp, argc, argv, 0, 0, pargs);

  cmd=pargs->args[0];
  fpkg->nrfiles=10;
  for(i=0;i<fpkg->nrfiles;i++)
    fpkg->bOpen[i]=0;

  i=0;icmd=-1;
  while(!(cmd[i]=='\0')){
    cmd[i]=toupper(cmd[i]);
    i++;
  }
  for(i=1;i<=NRCMDS;i++){
    if(!(strcmp(cmd,commands[i-1]))){
      icmd=i;
    }
  }
  if(icmd==-1){
    fatal("COULD NOT FIND AN APPROPRIATE COMMAND");
  }else{
    pargs->cmd=icmd;
  }
  if (pargs->verbose)
    fprintf(stdout,"FOUND COMMAND = %d\n",pargs->cmd);

  if(  pargs->keep == 1  )
    if( pargs->cmd != 3)
      fatal("THE KEEP FLAG MUST BE COMBINED WITH A GRAND CANONICAL SIMULATION");

  /* OUTPUT */
  if (pargs->outfile){
    fpkg->fil[0] = fopen (pargs->outfile, "w");
    if(fpkg->fil[0]==NULL)
      fatal("COULD NOT OPEN OUTPUT FILE");
    fpkg->bOpen[0]=1;
    fpkg->ftype[0]='w';
  }
  else
    fpkg->fil[0]=stdout; 

  // If in verbose mode 
  if(pargs->verbose)
    fprintf(stdout,"%s\n","VERBOSE");
  
  if(pargs->infile){
    if(!file_ext(pargs->infile,"gro"))
      typ=etGRO;
    if(!file_ext(pargs->infile,"pdb"))
      typ=etPDB;
    if(!file_ext(pargs->infile,"xyz"))
      typ=etPDB;
    if(file_ext(pargs->infile,"gro") & file_ext(pargs->infile,"pdb") & file_ext(pargs->infile,"xyz"))
      fatal("INPUT FILE IS NOT A RECONGIZED FILETYPE");
    fpkg->fil[1]=fopen(pargs->infile,"r");
    if(fpkg->fil[1]==NULL)
      fatal("FAILED TO OPEN INPUT FILE");
    fpkg->bOpen[1]=1;
    fpkg->ftype[1]='r';
  }
  else
    fatal("NO COORDINATE FILE GIVEN");
  
  if(pargs->infile2){
    if(!file_ext(pargs->infile2,"top")){
      fpkg->fil[2]=fopen(pargs->infile2,"r");
      if(fpkg->fil[2] == NULL)
	fatal("FAILED TO OPEN INPUT TOP FILE");
      fpkg->bOpen[2] = 1 ;
      fpkg->ftype[2] ='r';
    }
    else
      fatal("NOT A RECOGNIZED TOP FILE");
  }else{
    fatal("NO SPECIFIC TOPOLOGY INFO GIVEN\n");
  }

  if(pargs->infile3){
    if(!file_ext(pargs->infile3,"inp")){
      fpkg->fil[3] = fopen(pargs->infile3,"r");
      if( fpkg->fil[3] == NULL )
	fatal("FAILED TO OPEN INPUT INPUT FILE");
      fpkg->bOpen[3] = 1 ;
      fpkg->ftype[3] ='r';
    }else{
      fatal("NOT A RECOGNIZED INP FILE");
    }
  }else{
    fatal("NO SPECIFIC INPUT INFO GIVEN\n");
  }

  if(pargs->cmd==2){
    fpkg->fil[4]=fopen("charge.dat","w");
    fpkg->bOpen[4]=1;
    fpkg->ftype[4]='w';
  } 

  if(pargs->outfile_E){
    fpkg->fil[5]=fopen(pargs->outfile_E,"w");
    fpkg->bOpen[5]=1;
    fpkg->ftype[5]='w';
  } else {
    fpkg->fil[5]=fopen("energy.dat","w");
    fpkg->bOpen[5]=1;
    fpkg->ftype[5]='w';
  }
  if(pargs->outfile_X){
    fpkg->fil[6]=fopen(pargs->outfile_X,"w");
    fpkg->bOpen[6]=1;
    fpkg->ftype[6]='w';
  } else {
    fpkg->fil[6]=fopen("outcoord.xyz","w");
    fpkg->bOpen[6]=1;
    fpkg->ftype[6]='w';
  }

  if (pargs->verbose)
    fprintf(stdout,"DONE WITH FILES\n");

  return 0;
}
