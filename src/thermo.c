#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "define.h"
#include "structs.h"
#include "units.h"
#include "tools.h"
#include "stochastic_temp.h"

int propagate_thermostat(t_frame pframe)
{
  int    i,j,k,l,nx,ny,nz;
  double dt,hdt;

  if(pframe->Tref!=-1){
    pframe->mu[0] = sqrt(1+(pframe->dt/pframe->taut)*(pframe->Tref/pframe->T-1.0) );
    if(pframe->pref!=-1){
      pframe->mu[1] *= pframe->mu[0];
    }
  }
  else{
    pframe->mu[0] = 1.0;
  }

  return(0);
}

int propagate_barostat(t_frame pframe)
{
  double V;

  V = pframe->box[XX][XX]*pframe->box[YY][YY]*pframe->box[ZZ][ZZ];
  pframe->mu[1] = 1.00000000000;

  if(pframe->pref!=-1){
    pframe->mu[1] += 3.0*(V*(pframe->p-pframe->pref)+2.0/BOLTZ/pframe->Tref)
                     + pframe->fipimi*pframe->dt*0.5 + pframe->fifimi/3.0*SQUARE(pframe->dt*0.5);
    pframe->mu[1] *= pframe->dt*0.5/pframe->W;
    pframe->mu[1]  = exp(-pframe->mu[1]);
  }

  return(0);
}

int propagate_velocity(t_frame pframe)
{
  int i,k;
  double half=0.5, dt;

  dt  = (pframe->dt)*half;

  for(i=0;i<pframe->nr_parts;i++){
    k = pframe->partyp_key[pframe->partyp[i]];
    pframe->v[i][XX]=(pframe->v[i][XX]+(pframe->f[i][XX])/(pframe->M[pframe->nr_unique*k+k])*dt);
    pframe->v[i][YY]=(pframe->v[i][YY]+(pframe->f[i][YY])/(pframe->M[pframe->nr_unique*k+k])*dt);
    pframe->v[i][ZZ]=(pframe->v[i][ZZ]+(pframe->f[i][ZZ])/(pframe->M[pframe->nr_unique*k+k])*dt);
  }

  return(0);
}

int update_conf(t_frame pframe)
{
  int i,k,bRCOM=0;
  double dt,rc,nx,ny,nz;
  double hboxx,hboxy,hboxz;
  vector vcom;
  char   err_string[STRLEN],tmp_str[STRLEN];

  dt=pframe->dt;
  pframe->vcom[XX]=0.0; pframe->vcom[YY]=0.0; pframe->vcom[ZZ]=0.0;
  for(i=0;i<pframe->nr_parts;i++){

    if(pframe->bDebug){
      sprintf(tmp_str,"\n R0 = ( % 10.7lf , % 10.7lf , % 10.7lf )",pframe->r[i][XX],pframe->r[i][YY],pframe->r[i][ZZ]);
    }

    nx = (floor(pframe->r[i][XX]/pframe->box[XX][XX]));
    pframe->r[i][XX] -= pframe->box[XX][XX]*nx;
    ny = (floor(pframe->r[i][YY]/pframe->box[YY][YY]));
    pframe->r[i][YY] -= pframe->box[YY][YY]*ny;
    nz = (floor(pframe->r[i][ZZ]/pframe->box[ZZ][ZZ]));
    pframe->r[i][ZZ] -= pframe->box[ZZ][ZZ]*nz;
    k = pframe->partyp_key[pframe->partyp[i]];

    if(pframe->bDebug && (nx*nx>1||ny*ny>1 ||nz*nz>1)){
      sprintf(err_string,"ONE PARTICLE 1 ( NR % d ) MOVED MORE THAN ONE BOXLENGTH \n N =  ( % lf , % lf , % lf )",i,nx,ny,nz);
      strcat(err_string,tmp_str);
      sprintf(tmp_str,"\n R  = ( % 10.7lf , % 10.7lf , % 10.7lf )",pframe->r[i][XX],pframe->r[i][YY],pframe->r[i][ZZ]);
      strcat(err_string,tmp_str);
      sprintf(tmp_str,"\n M  = ( % 10.7lf , % 10.7lf  )",pframe->M[pframe->nr_unique*k+k],pframe->Q[pframe->nr_unique*k+k]);
      strcat(err_string,tmp_str);
      sprintf(tmp_str,"\n B  = ( % 10.7lf , % 10.7lf , % 10.7lf )",pframe->box[XX][XX],pframe->box[YY][YY],pframe->box[ZZ][ZZ]);
      strcat(err_string,tmp_str);
      sprintf(tmp_str,"\nMU  = ( % 10.7lf , % 10.7lf ) W = % 10.7lf",pframe->mu[0],pframe->mu[1],pframe->W);
      strcat(err_string,tmp_str);
      fatal(err_string);
    }

    pframe->v[i][XX] = pframe->mu[0]*(pframe->v[i][XX]+(pframe->f[i][XX])/(pframe->M[pframe->nr_unique*k+k])*dt);
    pframe->v[i][YY] = pframe->mu[0]*(pframe->v[i][YY]+(pframe->f[i][YY])/(pframe->M[pframe->nr_unique*k+k])*dt);
    pframe->v[i][ZZ] = pframe->mu[0]*(pframe->v[i][ZZ]+(pframe->f[i][ZZ])/(pframe->M[pframe->nr_unique*k+k])*dt);

    pframe->vcom[XX] += pframe->v[i][XX]/pframe->nr_parts;
    pframe->vcom[YY] += pframe->v[i][YY]/pframe->nr_parts;
    pframe->vcom[ZZ] += pframe->v[i][ZZ]/pframe->nr_parts;

    if(pframe->mu[1]>0){
      pframe->box[XX][XX] *= pframe->mu[1];
      pframe->box[YY][YY] *= pframe->mu[1];
      pframe->box[YY][YY] *= pframe->mu[1];
    }

    pframe->r[i][XX] = pframe->mu[1] * ( pframe->r[i][XX] + pframe->v[i][XX]*dt );
    pframe->r[i][YY] = pframe->mu[1] * ( pframe->r[i][YY] + pframe->v[i][YY]*dt );
    pframe->r[i][ZZ] = pframe->mu[1] * ( pframe->r[i][ZZ] + pframe->v[i][ZZ]*dt );

    nx = (floor(pframe->r[i][XX]/pframe->box[XX][XX]));
    pframe->r[i][XX] -= pframe->box[XX][XX]*nx;
    ny = (floor(pframe->r[i][YY]/pframe->box[YY][YY]));
    pframe->r[i][YY] -= pframe->box[YY][YY]*ny;
    nz = (floor(pframe->r[i][ZZ]/pframe->box[ZZ][ZZ]));
    pframe->r[i][ZZ] -= pframe->box[ZZ][ZZ]*nz;
    k = pframe->partyp_key[pframe->partyp[i]];

    if(pframe->bDebug && (nx*nx>1||ny*ny>1 ||nz*nz>1)){
      sprintf(err_string,"ONE PARTICLE 2 ( NR % d ) MOVED MORE THAN ONE BOXLENGTH \n N =  ( % lf , % lf , % lf )",i,nx,ny,nz);
      strcat(err_string,tmp_str);
      sprintf(tmp_str,"\n R  = ( % 10.7lf , % 10.7lf , % 10.7lf )",pframe->r[i][XX],pframe->r[i][YY],pframe->r[i][ZZ]);
      strcat(err_string,tmp_str);
      sprintf(tmp_str,"\n M  = ( % 10.7lf , % 10.7lf  )",pframe->M[pframe->nr_unique*k+k],pframe->Q[pframe->nr_unique*k+k]);
      strcat(err_string,tmp_str);
      sprintf(tmp_str,"\n B  = ( % 10.7lf , % 10.7lf , % 10.7lf )",pframe->box[XX][XX],pframe->box[YY][YY],pframe->box[ZZ][ZZ]);
      strcat(err_string,tmp_str);
      sprintf(tmp_str,"\nMU  = ( % 10.7lf , % 10.7lf ) W = % 10.7lf",pframe->mu[0],pframe->mu[1],pframe->W);
      strcat(err_string,tmp_str);
      fatal(err_string);
    }

    pframe->f[i][XX] = 0.0;
    pframe->f[i][YY] = 0.0;
    pframe->f[i][ZZ] = 0.0;
  }

  return 0;
}

int compute_pressure(t_frame pframe, t_files filepkg)
{
  double v,kpres,vpres,pcorr,PV=0.0,pres_unit,Ndf;
  double n,cutoff=1.0,V;
  int    i,k;

  // FROM KJ MOL^-1 TO BAR 
  pres_unit=1.0/AVOGADRO*1E25; // 1E27*1E3*1E-5;
  V = pframe->box[XX][XX]*pframe->box[YY][YY]*pframe->box[ZZ][ZZ];

  cutoff= sqrt(pframe->rc2);
  n     = pframe->nr_parts;
  v     = pframe->box[XX][XX]*pframe->box[YY][YY]*pframe->box[ZZ][ZZ];

  pframe->fipimi = 0.0;
  pframe->fifimi = 0.0;

  for(i=0;i<pframe->nr_parts;i++){
    k               = pframe->partyp_key[pframe->partyp[i]];
    PV             += iprod(pframe->r[i],pframe->f[i]);
    pframe->fipimi  = iprod(pframe->f[i],pframe->v[i]);
    pframe->fifimi  = iprod(pframe->f[i],pframe->f[i])/pframe->M[pframe->nr_unique*k+k];
  }

  kpres        =  2.0 * pframe->E[0] / (3.0 * v);
  vpres        =  2.0 * PV / (3.0 * v);
  pcorr        = -4.0 * M_PI * n * n / (3.0 * v*v * cutoff*cutoff*cutoff)*pframe->C6av; 

  Ndf          =  3.0 * ( pframe->nr_parts-1 );
  pframe->W    = BOLTZMANN*SQUARE(pframe->pref/pframe->deBroglie)/AMU*1E-6*pframe->Tref*Ndf; // PISTON MASS 
  pframe->p    = (kpres + vpres + pcorr)*pres_unit;
  
  pframe->kpres= kpres*pres_unit;
  pframe->vpres= vpres*pres_unit;
  pframe->pcorr= pcorr*pres_unit;

  return 0;
}

int compute_temperature(t_frame pframe, t_files filepkg) {
  int i,k;
  double dt,rc,Ekin=0.0, Ndf=0.0,v2,unit,Ekin_ref,m0=0.0;
  
  unit=AMU*1E6; // nm / ps = 1E3 m/s squared AMU*1E6 IN [ J ] SO AMU*1E3 IN [ kJ ]
  dt=pframe->dt;
  pframe->pc2=0;

  for(i=0;i<pframe->nr_parts;i++){
    v2           = iprod(pframe->v[i],pframe->v[i]);
    k            = pframe->partyp_key[pframe->partyp[i]];
    Ekin        += v2*0.5*pframe->M[pframe->nr_unique*k+k];
    pframe->pc2 += v2*pframe->M[pframe->nr_unique*k+k]*AMU*(SI_c)*1E6*pframe->M[pframe->nr_unique*k+k]*AMU*(SI_c);
    m0          += pframe->M[pframe->nr_unique*k+k]/(pframe->nr_parts);
  }
  pframe->E[0]          = Ekin*unit*AVOGADRO*1E-3;
  Ndf                   = pframe->nr_parts>1?(3.0*(pframe->nr_parts-1)):(3.0);
  Ekin_ref              = pframe->Tref*Ndf*BOLTZ*0.5;

  pframe->deBroglie     = SI_hc/sqrt(pframe->Tref*Ndf*BOLTZMANN*AMU*m0*SI_c*SI_c)*1E9; // in nm
  pframe->deBroglie_vol = pframe->deBroglie*pframe->deBroglie*pframe->deBroglie;
  Ekin                  = rekin(Ekin, Ekin_ref, Ndf, pframe->taut/pframe->dt );
  pframe->T             = 2.0*Ekin*unit/Ndf/BOLTZMANN;

  if(0){
    fprintf(stderr,"\n%f %f %f %f | %f %f \n",pframe->E[0],Ndf,Ekin_ref,pframe->deBroglie,Ekin,pframe->T);
  }

  return 0;
}

