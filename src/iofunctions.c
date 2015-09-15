#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "structs.h"
#include "define.h"
#include "units.h"
#include "tools.h"

int get_input(char *cptr,int nrid, int nrtot, int *tmp_nr, double *dvec, int *ivec)
{
  char *array[50];
  int   loop;

  array[0]=strtok(cptr," ");

  if(array[0]==NULL)
    fatal("NO TEST TO SEARCH");
  for(loop=1;loop<50;loop++) { 
    array[loop]=strtok(NULL," ");
    if(array[loop]==NULL)
      break;
  }

  for(loop=0;loop<50;loop++) {
    if(array[loop]==NULL|| *tmp_nr == nrtot)
      break;
    if(loop<nrid){
      ivec[loop]=atoi(array[loop]);
    }else{
      dvec[loop-nrid]=atof(array[loop]);
    }
    *tmp_nr+=1;
  }

  return(0);
}

int assign_mol(molecule *mol, int nrid, int nrtot, double *dvec, int *ivec,int typ)
{
  int dnr=0,i;

  if(typ==0){
    mol->b_seq = realloc(mol->b_seq,sizeof(int)*2*mol->nr_b);
    for(i=0;i<nrid;i++){
      mol->b_seq[2*(mol->nr_b - 1) + i ]=ivec[i];
    }
    mol->k_bond = realloc(mol->k_bond,sizeof(double)*2*mol->nr_b);
    for(i=0;i<(nrtot-nrid);i++){
      mol->k_bond[2*(mol->nr_b - 1) + i ]=dvec[i];
    }
  }
  if(typ==1){
    mol->a_seq = realloc(mol->a_seq,sizeof(int)*3*mol->nr_a);
    for(i=0;i<nrid;i++){
      mol->a_seq[3*(mol->nr_a - 1) + i ]=ivec[i];
    }
    mol->k_angle = realloc(mol->k_angle,sizeof(double)*2*mol->nr_a);
    for(i=0;i<(nrtot-nrid);i++){
      mol->k_angle[2*(mol->nr_a - 1) + i ]=dvec[i];
    }
  }
  if(typ==2){
    mol->i_seq = realloc(mol->i_seq,sizeof(int)*4*mol->nr_i);
    for(i=0;i<nrid;i++){
      mol->i_seq[4*(mol->nr_i - 1) + i ]=ivec[i];
    }
    mol->k_improp = realloc(mol->k_improp,sizeof(double)*4*mol->nr_i);
    for(i=0;i<(nrtot-nrid);i++){
      mol->k_improp[4*(mol->nr_i - 1) + i ]=dvec[i];
    }
  }
  if(typ==3){
    mol->d_seq = realloc(mol->d_seq,sizeof(int)*4*mol->nr_d);
    for(i=0;i<nrid;i++){
      mol->d_seq[nrid*(mol->nr_d - 1) + i ]=ivec[i];
    }
    mol->k_dihed = realloc(mol->k_dihed,sizeof(double)*4*mol->nr_d);
    for(i=0;i<(nrtot-nrid);i++){
      mol->k_dihed[4*(mol->nr_d - 1) + i ]=dvec[i];
    }
  }

  return(0);
}

int init_particles(t_frame pframe,t_parameters params, t_files filepkg)
{
  int  bHasVelocities = 0;
  int  i,j,k,l,f = 1;
  int  typ,nr,resnr,restyp,bSame=0,bSE=0;
  char line[STRLEN],partyp[STRLEN];
  double count;

  // SETS THE INITIAL COORDINATES AND FINDS CORRECT C6 C12
  // FIRST GET HEADER
  fgets(line,STRLEN,filepkg->fil[1]);
  sscanf(line,"%d %lf ",&nr,&(pframe->box[XX][XX]));
  pframe->box[XX][XX]*=0.1;
  pframe->box[YY][YY]=pframe->box[XX][XX];
  pframe->box[ZZ][ZZ]=pframe->box[XX][XX];

  if(pframe->nr_parts!=nr && pframe->bSystem == 1){
    fprintf(stderr,"NR_top = %d NR_xyz = %d \n",pframe->nr_parts,nr);
    fatal("TOPOLOGY AND COORDINATES DO NOT MATCH ( NR )");
  }
  pframe->nr_parts = nr;

  i=(pframe->nr_parts)*((pframe->nr_parts>=0)?1.0:-1.0);
  j=pframe->box[XX][XX]*((pframe->box[XX][XX]>=0)?1.0:-1.0);
  fprintf(filepkg->fil[0],"INIT INFO :: GOT %d PARTICLES\n",(pframe->nr_parts));

  if( !(i<10000000 && j<10000000) && pframe->nr_parts>0 && pframe->box[XX][XX] > 0 && bHasVelocities >= 0){
    fatal("WRONG WITH COORDINATE INPUT");
  }

  // ALLOCATE MEMORY
  pframe->r        = malloc(sizeof(vector)*(pframe->nr_parts));
  pframe->v        = malloc(sizeof(vector)*(pframe->nr_parts)*2); // leap frog needs half step velocities
  pframe->f        = malloc(sizeof(vector)*(pframe->nr_parts));
  pframe->restyp   = malloc(sizeof(int)*(pframe->nr_parts));
  pframe->resnr    = malloc(sizeof(int)*(pframe->nr_parts));
  pframe->partyp   = malloc(sizeof(int)*(pframe->nr_parts));
  pframe->nr_alloc = pframe->nr_parts;

  fgets(line,STRLEN,filepkg->fil[f]); //title

  count=0.0; pframe->C6av=0.0;
  for(i=0;i<pframe->nr_parts;i++){
    if(!feof(filepkg->fil[f])){
      if(fgets(line,STRLEN,filepkg->fil[f])==NULL){
        fatal("GOT EMPTY LINE FROM COORDINATE FILE");
      }
      sscanf(line," %s %lf %lf %lf ",partyp, &(pframe->r[i][XX]),
                                   &(pframe->r[i][YY]) , &(pframe->r[i][ZZ]));
      pframe->partyp[i]=atoi(partyp);
      if(pframe->partyp[i]<=0||pframe->partyp[i]>100){ // means that partyp is not a valid number so treat it as a string
        for(k=0;k<10;k++){
            if(partyp[k]=='\0')
              break;
        }
        for(j=1;j<100;j++){
          if(!strncmp(partyp,atom.name[j],k)){
            pframe->partyp[i]=j;
            break;
          }
        }
      }
      pframe->r[i][XX]*=0.1;
      pframe->r[i][YY]*=0.1;
      pframe->r[i][ZZ]*=0.1;
      k = pframe->partyp_key[pframe->partyp[i]];
      if( pframe->M[pframe->nr_unique*k+k] <= 0.0 ){
        fprintf(stderr,"GOT NR = %d, ID = %d, MASS = % lf\n",i,k,pframe->M[pframe->nr_unique*k+k]);
        fatal("NEGATIVE OR ZERO MASS ENCOUNTERED");
      }
      if( pframe->partyp[i] <= 0 || pframe->partyp[i] > 99 ){
        fprintf(stderr,"GOT Z = %d\n",pframe->partyp[i]);
        fatal("NON REASONABLE ATOM ENCOUNTERED");
      }
      if( pframe->C6[pframe->nr_unique*k+k] < 0.0 || pframe->C12[pframe->nr_unique*k+k] < 0.0 ){
        fatal("NEGATIVE OR ZERO C6 or C12 ENCOUNTERED");
      }
      if( pframe->C6[pframe->nr_unique*k+k] > 0.0 && pframe->C12[pframe->nr_unique*k+k] > 0.0 ){
        pframe->C6av += pframe->C6[pframe->nr_unique*k+k];
        count += 1.0;
      }
    }else{
      break;
    }
  }
  if( i < (pframe->nr_parts)-1 ){
    fatal("FILE ENDED BEFORE ALL PARTICLES WERE READ");
  }
  pframe->C6av/=count;
  if(pframe->bDebug)
    fprintf(filepkg->fil[0],"INIT INFO :: < C6 > = %f %f\n",pframe->C6av,count);

  // TYPE CHECKING BLOCK
  if(pframe->bSystem){
    for(i=0;i<pframe->nr_umols;i++){
      for(j=0;j<pframe->nr_mols[i];j++){
        if(pframe->mol[i]->bZ==0){
          for(k=0;k<pframe->mol[pframe->mol_seq[i]]->nr_t;k++){
            pframe->resnr [j*pframe->mol[i]->nr_t+k+pframe->n_start[i]] = j;
            pframe->restyp[j*pframe->mol[i]->nr_t+k+pframe->n_start[i]] = i;
            if(pframe->partyp[j*pframe->mol[pframe->mol_seq[i]]->nr_t+k+pframe->n_start[i]] != pframe->mol[pframe->mol_seq[i]]->t_seq[k]){
              printf("TYPES %d %d ( %d :: %d ) %d %s %d %d\n",pframe->partyp[j*pframe->mol[pframe->mol_seq[i]]->nr_t+k+pframe->n_start[i]],pframe->mol[pframe->mol_seq[i]]->t_seq[k],j*pframe->mol[pframe->mol_seq[i]]->nr_t+k+pframe->n_start[i],i,j,pframe->mol[i]->name,pframe->mol_seq[i],pframe->n_start[i]);
              fatal("TYPES IN TOPOLOGY AND COORDINATE FILES DO NOT MATCH");
            }
          }

          for(k=0;k<pframe->mol[pframe->mol_seq[i]]->nr_b;k++){
            if( pframe->partyp[j*pframe->mol[pframe->mol_seq[i]]->nr_t+pframe->mol[pframe->mol_seq[i]]->b_seq[2*k]+pframe->n_start[i]] != pframe->mol[pframe->mol_seq[i]]->t_seq[pframe->mol[pframe->mol_seq[i]]->b_seq[2*k]] )
              fatal("TYPES IN TOPOLOGY AND COORDINATE FILES DO NOT MATCH BOND P1");
            if( pframe->partyp[j*pframe->mol[pframe->mol_seq[i]]->nr_t+pframe->mol[pframe->mol_seq[i]]->b_seq[2*k+1]+pframe->n_start[i]] != pframe->mol[pframe->mol_seq[i]]->t_seq[pframe->mol[pframe->mol_seq[i]]->b_seq[2*k+1]] )
              fatal("TYPES IN TOPOLOGY AND COORDINATE FILES DO NOT MATCH BOND P2");
          }

          for(k=0;k<pframe->mol[pframe->mol_seq[i]]->nr_a;k++){
            if( pframe->partyp[j*pframe->mol[pframe->mol_seq[i]]->nr_t+pframe->mol[pframe->mol_seq[i]]->a_seq[2*k]+pframe->n_start[i]] != pframe->mol[pframe->mol_seq[i]]->t_seq[pframe->mol[pframe->mol_seq[i]]->a_seq[2*k]] )
              fatal("TYPES IN TOPOLOGY AND COORDINATE FILES DO NOT MATCH ANGLE P1");
            if( pframe->partyp[j*pframe->mol[pframe->mol_seq[i]]->nr_t+pframe->mol[pframe->mol_seq[i]]->a_seq[2*k+1]+pframe->n_start[i]] != pframe->mol[pframe->mol_seq[i]]->t_seq[pframe->mol[pframe->mol_seq[i]]->a_seq[2*k+1]] )
              fatal("TYPES IN TOPOLOGY AND COORDINATE FILES DO NOT MATCH ANGLE P2");
            if( pframe->partyp[j*pframe->mol[pframe->mol_seq[i]]->nr_t+pframe->mol[pframe->mol_seq[i]]->a_seq[2*k+2]+pframe->n_start[i]] != pframe->mol[pframe->mol_seq[i]]->t_seq[pframe->mol[pframe->mol_seq[i]]->a_seq[2*k+2]] )
              fatal("TYPES IN TOPOLOGY AND COORDINATE FILES DO NOT MATCH ANGLE P3");
          }
        }
        if(pframe->mol[i]->bZ==1){
          for(k=0;k<pframe->mol[pframe->mol_seq[i]]->nr_t;k++){
            pframe->resnr [j*pframe->mol[i]->nr_t+k+pframe->n_start[i]] = j;
            pframe->restyp[j*pframe->mol[i]->nr_t+k+pframe->n_start[i]] = i;
            if(pframe->partyp[j*pframe->mol[pframe->mol_seq[i]]->nr_t+k+pframe->n_start[i]] != pframe->mol[pframe->mol_seq[i]]->z_seq[k]){
              printf("TYPES %d %d ( %d :: %d ) %d %s %d %d| %d,%d,%d\n",pframe->partyp[j*pframe->mol[pframe->mol_seq[i]]->nr_t+k+pframe->n_start[i]],pframe->mol[pframe->mol_seq[i]]->z_seq[k],j*pframe->mol[pframe->mol_seq[i]]->nr_t+k+pframe->n_start[i],i,j,pframe->mol[i]->name,pframe->mol_seq[i],pframe->n_start[i],i,j,k);
              fatal("TYPES IN TOPOLOGY AND COORDINATE FILES DO NOT MATCH");
            }
          }

          for(k=0;k<pframe->mol[pframe->mol_seq[i]]->nr_b;k++){
            if( pframe->partyp[j*pframe->mol[pframe->mol_seq[i]]->nr_t+pframe->mol[pframe->mol_seq[i]]->b_seq[2*k]+pframe->n_start[i]] != pframe->mol[pframe->mol_seq[i]]->z_seq[pframe->mol[pframe->mol_seq[i]]->b_seq[2*k]] )
              fatal("TYPES IN TOPOLOGY AND COORDINATE FILES DO NOT MATCH BOND P1");
            if( pframe->partyp[j*pframe->mol[pframe->mol_seq[i]]->nr_t+pframe->mol[pframe->mol_seq[i]]->b_seq[2*k+1]+pframe->n_start[i]] != pframe->mol[pframe->mol_seq[i]]->z_seq[pframe->mol[pframe->mol_seq[i]]->b_seq[2*k+1]] )
              fatal("TYPES IN TOPOLOGY AND COORDINATE FILES DO NOT MATCH BOND P2");
          }

          for(k=0;k<pframe->mol[pframe->mol_seq[i]]->nr_a;k++){
            if( pframe->partyp[j*pframe->mol[pframe->mol_seq[i]]->nr_t+pframe->mol[pframe->mol_seq[i]]->a_seq[2*k]+pframe->n_start[i]] != pframe->mol[pframe->mol_seq[i]]->z_seq[pframe->mol[pframe->mol_seq[i]]->a_seq[2*k]] )
              fatal("TYPES IN TOPOLOGY AND COORDINATE FILES DO NOT MATCH ANGLE P1");
            if( pframe->partyp[j*pframe->mol[pframe->mol_seq[i]]->nr_t+pframe->mol[pframe->mol_seq[i]]->a_seq[2*k+1]+pframe->n_start[i]] != pframe->mol[pframe->mol_seq[i]]->z_seq[pframe->mol[pframe->mol_seq[i]]->a_seq[2*k+1]] )
              fatal("TYPES IN TOPOLOGY AND COORDINATE FILES DO NOT MATCH ANGLE P2");
            if( pframe->partyp[j*pframe->mol[pframe->mol_seq[i]]->nr_t+pframe->mol[pframe->mol_seq[i]]->a_seq[2*k+2]+pframe->n_start[i]] != pframe->mol[pframe->mol_seq[i]]->z_seq[pframe->mol[pframe->mol_seq[i]]->a_seq[2*k+2]] )
              fatal("TYPES IN TOPOLOGY AND COORDINATE FILES DO NOT MATCH ANGLE P3");
          }
          for( k = 0 ; k < pframe->mol[pframe->mol_seq[i]]->nr_t ; k++ ){
            pframe->partyp[j*pframe->mol[pframe->mol_seq[i]]->nr_t+k+pframe->n_start[i]] = pframe->mol[pframe->mol_seq[i]]->t_seq[k]; //exchange z with type sequence
          }
        }        
      }
    }
  }

  if(pframe->bIR){
    i=pframe->nr_umols-1;
    j=pframe->nr_mols[i]-1;
    pframe->ir       = malloc( sizeof( vector )*( pframe->mol[pframe->mol_seq[i]]->nr_t ) );
    for(k=0;k<pframe->mol[pframe->mol_seq[i]]->nr_t;k++){
      l=j*pframe->mol[pframe->mol_seq[i]]->nr_t+k+pframe->n_start[i];
      pframe->ir[k][XX]=pframe->r[l][XX];
      pframe->ir[k][YY]=pframe->r[l][YY];
      pframe->ir[k][ZZ]=pframe->r[l][ZZ];
      if(pframe->bDebug)
        fprintf(stderr,"[%lf %lf %lf]\n",pframe->ir[k][XX],pframe->ir[k][YY],pframe->ir[k][ZZ]);
    }
  }

  return 0;
}

int print_coord(t_frame pframe, t_files filepkg)
{
  int i,j,k,I;

  fprintf(filepkg->fil[6],"% d % f \n OUTPUT\n",pframe->nr_parts,pframe->box[XX][XX]*10.0);
  for(i=0;i<pframe->nr_umols;i++){
    for(j=0;j<pframe->nr_mols[i];j++){
      for( k = 0 ; k < pframe->mol[pframe->mol_seq[i]]->nr_t ; k++ ){
        I = j*pframe->mol[pframe->mol_seq[i]]->nr_t+k+pframe->n_start[i];
        if( pframe->mol[i]->bZ == 0 ){
          fprintf(filepkg->fil[6],"%d %lf %lf %lf\n",(pframe->partyp[I]), (pframe->r[I][XX]*10.0),
                                               (pframe->r[I][YY]*10.0) ,(pframe->r[I][ZZ]*10.0));
        }
        if( pframe->mol[i]->bZ == 1 ){
          fprintf(filepkg->fil[6],"%d %lf %lf %lf\n",pframe->mol[pframe->mol_seq[i]]->z_seq[k], (pframe->r[I][XX]*10.0),
                                               (pframe->r[I][YY]*10.0) ,(pframe->r[I][ZZ]*10.0));
        }
      }
    }
  }
  return 0;
}

int print_charge(t_frame pframe, t_files filepkg)
{
  int i;

  fprintf(filepkg->fil[3],"\n" );
  for(i=0;i<pframe->nr_parts;i++){
    fprintf(filepkg->fil[3],"% 10.6lf ",(pframe->Q[i]) );
  }
  return 0;
}

int read_topology(t_topology ptop, t_files filepkg)
{
  return 0;
}

int read_params(t_frame pframe, t_arguments pargs, t_parameters params, t_files filepkg)
{
  int      bHasVelocities = 0, nrid, nrtot;
  int      i, j, k, l, f  = 1,sys_nr=0;
  int      typ, nr, resnr, restyp, bSame = 0, bSE = 0;
  int      TYPES[STRLEN], type_flag, interaction_flag = 0;
  int      SEQ[STRLEN], nr_mols = 0;
  int      tmp_flag, tmp_nr, prop_nr,ivec[STRLEN];
  molecule   mole;
  t_molecule mol;
  double   C[2],Q,M,s6,s12,nrc6=0,tmp_real,dvec[STRLEN];
  char     line[STRLEN],key[STRLEN];
  char     value[STRLEN],delimiter[STRLEN];
  char     *cptr;

  pframe->cmd = pargs->cmd;
  mol=&mole;
  
  // FIRST SET INPUT DEFAULTS
  pframe->NSTEPS        =  0;
  pframe->dt            =  0;
  pframe->pref          = -1;
  pframe->Tref          = -1;
  pframe->bDebug        =  0;
  pframe->rc2           = 1.5*1.5;
  pframe->bGenv         =  0;
  pframe->mu[3]         = -1;
  pframe->C6av          =  0.0;
  pframe->deBroglie     =  1.0; 
  pframe->deBroglie_vol =  1.0; 
  pframe->bIR           = pargs->keep;
  pframe->bSystem       = 0;
  while(!feof(filepkg->fil[3])){
    fgets(line,STRLEN,filepkg->fil[3]); 
    if(line[0]=='#'||line[0]=='&'||line[0]=='@'|| line[0]=='['){
      continue;
    }
    sscanf(line,"%s %s %s",key,delimiter,value);
    if(!strcmp(key,"NSTEPS")){
      pframe->NSTEPS=atoi(value);
    }
    if(!strcmp(key,"DT")){
      pframe->dt=atof(value);
    }
    if(!strcmp(key,"TREF")){
      pframe->Tref=atof(value);
    }
    if(!strcmp(key,"TAUT")){
      pframe->taut=atof(value);
    }
    if(!strcmp(key,"TAUP")){
      pframe->taup=atof(value);
    }
    if(!strcmp(key,"PREF")){
      pframe->pref=atof(value);
    }
    if(!strcmp(key,"DEBUG")){
      pframe->bDebug=1;
    }
    if(!strcmp(key,"CUTOFF")){
      pframe->rc2=atof(value)*atof(value);   
    }
    if(!strcmp(key,"GENV")){
      pframe->bGenv=1;
    }
    if(!strcmp(key,"MU")){ // has its value given in [ K ]
      pframe->mu[2]=atof(value);
    }
  }
  
  if( pframe->mu[2]==-1 && pframe->cmd == cmdMC )
    fatal("FAILED IN SETTING UP GRAND CANONICAL SIMULATION");

  if(pframe->NSTEPS==0||pframe->dt==0)
    fatal("FAILED IN SETTING UP THE SYSTEM 1");
  if(pframe->Tref==-1&&pframe->bGenv==1)
    fatal("FAILED IN SETTING UP TEMPERATURE GENERATION ");
  if(pframe->Tref!=-1&&(pframe->taut<0||pframe->taut>1000.0))
    fatal("FAILED IN SETTING UP TEMPERATURE BATH FOR THE SYSTEM");
  if(pframe->pref!=-1&&(pframe->taup<0||pframe->taup>1000.0))
    fatal("FAILED IN SETTING UP PRESSURE BATH FOR THE SYSTEM 3");

  // NOW LOAD TOPOLOGY
  mol->t_seq    = malloc(sizeof(int));
  mol->z_seq    = malloc(sizeof(int));
  mol->a_seq    = malloc(sizeof(int));
  mol->b_seq    = malloc(sizeof(int));
  mol->i_seq    = malloc(sizeof(int));
  mol->d_seq    = malloc(sizeof(int));
  mol->k_bond   = malloc(sizeof(double));
  mol->k_angle  = malloc(sizeof(double));
  mol->k_improp = malloc(sizeof(double));
  mol->k_dihed  = malloc(sizeof(double));
  pframe->nr_parts = 0;

  mol->nr_t = 0;
  mol->nr_b = 0;
  mol->nr_a = 0;
  mol->nr_i = 0;
  mol->nr_d = 0;
  mol->bZ   = 0;

  TYPES[0]=0;
  while(!feof(filepkg->fil[2])){
    sprintf(line," ");
    fgets(line,STRLEN,filepkg->fil[2]); 
    if( line[0]=='#' || line[0]=='[' || line[0]=='@' || line[0]=='\n' ){
      continue;
    }
    if(!strncmp(line,"&nonbonded",6)){
      interaction_flag=0;
      if(pframe->bDebug){
        fprintf(filepkg->fil[0],"INIT INFO :: NONBOND\n");
      }    
      continue;
    }
    if(!strncmp(line,"&molecule",6)){
      interaction_flag=1;
      prop_nr=0;
      nr_mols++;
      if(nr_mols>(STRLEN-2)){
        fatal("CURRENTLY ONLY SUPPORT A LIMITED NUMBER OF MOLECULES IN THE SYSTEM");
      }
      if(pframe->bDebug){
        fprintf(filepkg->fil[0],"INIT INFO :: MOLECULE\n");
      }    
      continue;
    }
    if(!strncmp(line,"&system",6)){
      interaction_flag=2;
      sys_nr=0;
      pframe->nr_mols      = malloc(sizeof(int)*nr_mols);
      if(pframe->bDebug){
        fprintf(filepkg->fil[0],"INIT INFO :: SYSTEM\n");
      }    
      pframe->bSystem=1;
      continue;
    }
    if(line[0]=='&' && interaction_flag != -1) {
      if(interaction_flag==1 && prop_nr < 2) {
        fatal("BADLY FORMATED MOLECULE BLOCK");
      }
      if( interaction_flag == 1 ){ // END OF MOLECULE BLOCK
        // ALLOC MOLECULES
        pframe->mol[nr_mols-1]          = malloc(sizeof(molecule)*1);

        sprintf(pframe->mol[nr_mols-1]->name,"%s",mol->name);
        pframe->mol[nr_mols-1]->t_seq   = malloc(sizeof(int)*mol->nr_t);
        pframe->mol[nr_mols-1]->z_seq   = malloc(sizeof(int)*mol->nr_t);
        pframe->mol[nr_mols-1]->nr_t    = mol->nr_t;
        for(i=0;i<mol->nr_t;i++){
          pframe->mol[nr_mols-1]->t_seq[i]=mol->t_seq[i];
          if(mol->bZ){
            pframe->mol[nr_mols-1]->z_seq[i]=mol->z_seq[i];
          }
        }
        pframe->mol[nr_mols-1]->bZ      = mol->bZ;
        pframe->mol[nr_mols-1]->b_seq   = malloc(sizeof(int)*2*mol->nr_b);
        pframe->mol[nr_mols-1]->k_bond  = malloc(sizeof(double)*2*mol->nr_b);
        pframe->mol[nr_mols-1]->nr_b    = mol->nr_b;
        for(j=0;j<mol->nr_b;j++){
          for(i=0;i<2;i++){
            pframe->mol[nr_mols-1]->b_seq[2*j + i ]  = mol->b_seq[2*j + i ];
          }
          for(i=0;i<2;i++){
            pframe->mol[nr_mols-1]->k_bond[2*j + i ] = mol->k_bond[2*j + i ];
          }
        }
        pframe->mol[nr_mols-1]->a_seq   = malloc(sizeof(int)*3*mol->nr_a);
        pframe->mol[nr_mols-1]->k_angle = malloc(sizeof(double)*2*mol->nr_a);
        pframe->mol[nr_mols-1]->nr_a    = mol->nr_a;
        for(j=0;j<mol->nr_a;j++){
          for(i=0;i<3;i++){
            pframe->mol[nr_mols-1]->a_seq[3*j + i ] = mol->a_seq[3*j + i];
          }
          for(i=0;i<2;i++){
             pframe->mol[nr_mols-1]->k_angle[2*j + i ] = mol->k_angle[2*j + i ];
          }
        }  
        pframe->mol[nr_mols-1]->i_seq    = malloc(sizeof(int)*4*mol->nr_i);
        pframe->mol[nr_mols-1]->k_improp = malloc(sizeof(double)*4*mol->nr_i);
        pframe->mol[nr_mols-1]->nr_i     = mol->nr_i;
        for(j=0;j<mol->nr_i;j++){
          for(i=0;i<4;i++){
            pframe->mol[nr_mols-1]->i_seq[4*j + i ] = mol->i_seq[4*j + i];
          }
          for(i=0;i<2;i++){
             pframe->mol[nr_mols-1]->k_improp[4*j + i ] = mol->k_improp[4*j + i ];
          }
        }  
        pframe->mol[nr_mols-1]->d_seq    = malloc(sizeof(int)*4*mol->nr_d);
        pframe->mol[nr_mols-1]->k_dihed  = malloc(sizeof(double)*4*mol->nr_d);
        pframe->mol[nr_mols-1]->nr_d     = mol->nr_d;
        for(j=0;j<mol->nr_d;j++){
          for(i=0;i<4;i++){
            pframe->mol[nr_mols-1]->d_seq[4*j + i ] = mol->d_seq[4*j + i];
          }
          for(i=0;i<3;i++){
             pframe->mol[nr_mols-1]->k_dihed[4*j + i ] = mol->k_dihed[4*j + i ];
          }
        }
        if(pframe->bDebug){
          print_molecule(mol,filepkg);
          print_molecule(pframe->mol[nr_mols-1],filepkg);
        }
        mol->nr_t=0;
        mol->nr_b=0;
        mol->nr_a=0;
        mol->nr_i=0;
        mol->nr_d=0;
      }
      interaction_flag=-1;
    }
    if(interaction_flag==-1){
      continue;
    }
    switch(interaction_flag){
      case 0:
        sscanf(line," %d %d %lf %lf %lf %lf ",&nr,&typ,&C[0],&C[1],&Q,&M);
        if(nr<=0 || nr>=1000){
          continue;
        }
        if(typ==1){
          type_flag=0;
          for(i=0;i<TYPES[0];i++){
            if(TYPES[i+1] == nr){
              type_flag=1;
            }
          }
          if( type_flag == 0 ) { // UNIQUE
            TYPES[0]+=1;
            TYPES[TYPES[0]]=nr;
            s6 = C[0]*C[0]*C[0]*C[0]*C[0]*C[0];
            params->NB[nr][0] = 4.0*C[1]*s6*RGAS/KILO; // kJ mol^-1
            params->NB[nr][1] = 4.0*C[1]*s6*s6*RGAS/KILO;
            if( params->NB[nr][0] > 0.0 ){
              pframe->C6av     += params->NB[nr][0];
              nrc6             += 1.0;
            }
            params->NB[nr][2] = Q;
            params->NB[nr][3] = M;
            if(pframe->bDebug)
              fprintf(filepkg->fil[0],"INIT INFO :: %d %d | %f %f | %20.17f %20.17f | %20.17f %20.17f\n",nr,TYPES[0],C[0],C[1],params->NB[nr][0],params->NB[nr][1],params->NB[nr][2],params->NB[nr][3]);
          }
        }else{
          fatal("ONLY LENNARD JONES PARAMETERS GIVEN BY SIGMA IN [ NM ] AND EPSILON IN [ K ] ");
        }
        break;
      case 1:
        if(line[0]=='S'){
          prop_nr++;
          if(pframe->bDebug){
            fprintf(filepkg->fil[0],"INIT INFO :: FOUND SEQUENCE > ");
          }
          tmp_nr=0; tmp_flag=0; cptr=&line[1];
          for(i=2;i<STRLEN;i++){
            if((line[i] != ' '||line[i]!='\n'||line[i]!='\0') && tmp_flag==0){
              tmp_flag = 1;
            }
            if((line[i]==' '||line[i]=='\n'||line[i]=='\0') && tmp_flag==1){
              tmp_flag = atoi(cptr);
              if(tmp_flag<=0)
                continue;
              if(pframe->bDebug){
                fprintf(filepkg->fil[0]," %d ",tmp_flag);
              }
              ivec[tmp_nr] = tmp_flag;
              cptr=&line[i]; tmp_nr++;
              tmp_flag = 0;
            }
          }
          if(pframe->bDebug){
            fprintf(filepkg->fil[0]," NR = %d\n",tmp_nr);
          }
          mol->nr_t=tmp_nr;
          mol->t_seq=realloc(mol->t_seq,sizeof(int)*tmp_nr);
          for(i=0;i<mol->nr_t;i++){
            mol->t_seq[i]=ivec[i];
          }
        }
        if(line[0]=='Z'){
          if(prop_nr==0){
            fatal("YOU MUST FIRST SPECIFY THE PARTICLE LABELS");
          }
          if(pframe->bDebug){
            fprintf(filepkg->fil[0],"INIT INFO :: FOUND ATOMIC INDECES > ");
          }
          tmp_nr=0; tmp_flag=0; cptr=&line[1];
          for(i=2;i<STRLEN;i++){
            if((line[i] != ' '||line[i]!='\n'||line[i]!='\0') && tmp_flag==0){
              tmp_flag = 1;
            }
            if((line[i]==' '||line[i]=='\n'||line[i]=='\0') && tmp_flag==1){
              tmp_flag = atoi(cptr);
              if(tmp_flag<=0)
                continue;
              if(pframe->bDebug){
                fprintf(filepkg->fil[0]," %d ",tmp_flag);
              }
              ivec[tmp_nr] = tmp_flag;
              cptr=&line[i]; tmp_nr++;
              tmp_flag = 0;
            }
          }
          if( mol->nr_t != tmp_nr ){
            fatal("MUST HAVE THE SAME NUMBER OF ATOMIC INDECES AS PARTICLES");
          }
          mol->bZ = 1;
          mol->z_seq=realloc(mol->z_seq,sizeof(int)*tmp_nr);
          for(i=0;i<mol->nr_t;i++){
            mol->z_seq[i]=ivec[i];
          }
        }
        if(line[0]=='M'){
          prop_nr++;
          cptr=&line[2];
          for(i=2;i<STRLEN;i++){
            line[i]=line[i]=='\n'?'\0':line[i];
          } 
          if(pframe->bDebug){
            fprintf(filepkg->fil[0],"INIT INFO :: FOUND NAME %s\n",cptr);
          }
          sprintf(mol->name,"%s",cptr);
        }
        if(line[0]=='B'){
          if(prop_nr<2){
            fprintf(filepkg->fil[0],"\n");
            fatal("BAD MOLECULE FORMAT");
          }
          prop_nr++;
          tmp_nr=0; cptr=&line[1];
          nrid=2; nrtot=4;
          get_input(cptr,nrid,nrtot,&tmp_nr,dvec,ivec);
          mol->nr_b  += 1;
          assign_mol(mol, nrid, nrtot, dvec, ivec,0);
          if(tmp_nr<nrtot){
            fprintf(filepkg->fil[0],"\n");
            fatal("ERROR IN BOND FORMAT");
          }
        }
        if(line[0]=='A'){
          if(prop_nr<2){
            fatal("BAD MOLECULE FORMAT");
          }
          prop_nr++;
          tmp_nr=0; cptr=&line[1];
          nrid=3; nrtot=5;
          get_input(cptr,nrid,nrtot,&tmp_nr,dvec,ivec);
          mol->nr_a  += 1;
          assign_mol(mol, nrid, nrtot, dvec, ivec,1);
          if(tmp_nr<5){
            fprintf(filepkg->fil[0],"\n");
            fatal("ERROR IN ANGLE FORMAT");
          }
        }
        if(line[0]=='I'){ // IMPROPER
          if(prop_nr<2){
            fatal("BAD MOLECULE FORMAT");
          }
          prop_nr++;
          tmp_nr=0; tmp_flag=0; cptr=&line[1];
          nrid=4; nrtot=6;
          get_input(cptr,nrid,nrtot,&tmp_nr,dvec,ivec);
          mol->nr_i  += 1;
          assign_mol(mol, nrid, nrtot, dvec, ivec,2);
          if(tmp_nr<6){
            fprintf(filepkg->fil[0],"\n");
            fatal("ERROR IN IMPROPER FORMAT");
          }
        }
        if(line[0]=='D'){ // DIHEDRAL
          if(prop_nr<2){
            fatal("BAD MOLECULE FORMAT");
          }
          prop_nr++;
          tmp_nr=0; cptr=&line[1];
          nrid=4; nrtot=7;
          get_input(cptr,nrid,nrtot,&tmp_nr,dvec,ivec);
          mol->nr_d  += 1;
          assign_mol(mol, nrid, nrtot, dvec, ivec,3);
          if(tmp_nr<7){
            fprintf(filepkg->fil[0],"\n");
            fatal("ERROR IN DIHEDRAL FORMAT");
          }
        }
        break;
      case 2:
          sscanf(line,"%s %d",delimiter,&tmp_flag);
          if(line[0]==';'){
            continue;
          }
          fprintf(filepkg->fil[0],"INIT INFO :: %s %d\n",delimiter,tmp_flag);
          pframe->n_start[0]=0;
          for(i=0;i<nr_mols;i++){
            if(!strcmp(delimiter,pframe->mol[i]->name)) // find matching molecule
              j=i;
          }
          if( j == -1 ){
            fatal("MOLECULE NOT SPECIFIED");
          }
          pframe->nr_mols[sys_nr]      = tmp_flag;
          if(sys_nr>STRLEN-2)
            fatal("TOO MANY MOLECULES SPECIFIED");
          pframe->mol_seq[sys_nr] = j;
          pframe->mol[j]->uid=sys_nr;
          pframe->mol[j]->pid=j;
          sys_nr++;
          pframe->nr_umols=sys_nr;
          pframe->n_start[sys_nr] = tmp_flag*pframe->mol[j]->nr_t;
          pframe->nr_parts       += tmp_flag*pframe->mol[j]->nr_t;
          fprintf(filepkg->fil[0],"INIT NUMR :: %d %d %d\n",pframe->nr_parts,tmp_flag,mol->nr_t);
          fprintf(filepkg->fil[0],"INIT MSEQ :: %d %d (->%d) [ %d ]\n",sys_nr-1,pframe->mol_seq[sys_nr-1], pframe->n_start[sys_nr],sys_nr);
        break;
    }
  }  

  // WE RECALCULATE THIS AFTER THE COORDINATES HAVE BEEN READ
  pframe->C6av /= nrc6;
  if(pframe->bDebug){
    fprintf(filepkg->fil[0],"INIT INFO :: < C6_top > = %f  NR UNIQUE TYPES = %d\n",pframe->C6av,TYPES[0]);
  }

  // BUILDING INTERACTION RULES
  k=0;
  for(i=0;i<TYPES[0];i++){
    k=TYPES[i+1]>k?TYPES[i+1]:k;
  }
  pframe->partyp_key = malloc(sizeof(double)*(k+1));
  for(i=0;i<TYPES[0];i++){ 
    k=TYPES[i+1];
    pframe->partyp_key[k]=i;  // MAP TYPE TO MATRIX INDEX
    if(pframe->bDebug){
      fprintf(filepkg->fil[0],"INIT INFO :: %d -> %d ( %d )\n",k,pframe->partyp_key[k],TYPES[0]);
    }
  }
  pframe->C6       = malloc(sizeof(double)*(TYPES[0]*TYPES[0]));   //THESE ARE FIXED EVEN IN A GC SIMULATION
  pframe->C12      = malloc(sizeof(double)*(TYPES[0]*TYPES[0]));
  pframe->Q        = malloc(sizeof(double)*(TYPES[0]*TYPES[0]));
  pframe->M        = malloc(sizeof(double)*(TYPES[0]*TYPES[0]));
  pframe->nr_unique= TYPES[0];
// should implement a list override here
  for(i=0;i<TYPES[0];i++){ 
    if(pframe->bDebug)
      fprintf(filepkg->fil[0],"INIT INFO %d %d:: ",i,TYPES[i+1]);
    for(j=0;j<TYPES[0];j++){
      pframe->C6 [ i*TYPES[0]+j ] = sqrt( params->NB[ TYPES[i+1] ][0]*params->NB[ TYPES[j+1] ][0] );
      pframe->C12[ i*TYPES[0]+j ] = sqrt( params->NB[ TYPES[i+1] ][1]*params->NB[ TYPES[j+1] ][1] );
      pframe->Q  [ i*TYPES[0]+j ] =       params->NB[ TYPES[i+1] ][2]*params->NB[ TYPES[j+1] ][2]  ;
      pframe->M  [ i*TYPES[0]+j ] = sqrt( params->NB[ TYPES[i+1] ][3]*params->NB[ TYPES[j+1] ][3] );
      if(pframe->bDebug){
        fprintf(filepkg->fil[0]," ( % f,% f, %f, %f ) ",pframe->C6 [ i*TYPES[0]+j ],pframe->C12 [ i*TYPES[0]+j ],pframe->Q  [ i*TYPES[0]+j ],pframe->M  [ i*TYPES[0]+j ]);
      }
    }
    if(pframe->bDebug)
      fprintf(filepkg->fil[0],"\n");
  }
  free(mol->t_seq);
  free(mol->z_seq);
  free(mol->a_seq);
  free(mol->b_seq);
  free(mol->i_seq);
  free(mol->d_seq);
  free(mol->k_bond);
  free(mol->k_angle);
  free(mol->k_improp);
  free(mol->k_dihed);

  if(!pframe->bSystem)
    fatal("MUST SPECIFY SYSTEM SECTION");

  return 0;
}
