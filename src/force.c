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

int force_calc_bon(t_frame pframe, t_files fpkg)
{
  int     i,j,k,l,N,ci,cj,nqi,id1,id2,id3,id4;
  vector  *x,*v,dx,dy;
  t_vector distx,disty;
  double  dr,r,r2,r6,r12,ir,ir2,ir6,ir12,E=0.0,rc2;
  double  ar21,ar22,air21,air22,air11,air12,ang;
  double  C6ij,C12ij,FORCE,f,df,kc,ddr,fij;
  double  q2,dfq,fq,fe,phi,cosval,cosv2,dV;
  int     ist,isl,ds,m;

  distx=&dx; disty=&dy;
  for(i=0;i<pframe->nr_umols;i++){
    for(j=0;j<pframe->nr_mols[i];j++){
      // BONDED
      for(k=0;k<pframe->mol[pframe->mol_seq[i]]->nr_b;k++){
        id1  = j*pframe->mol[pframe->mol_seq[i]]->nr_t + pframe->mol[pframe->mol_seq[i]]->b_seq[2*k  ]+pframe->n_start[i];
        id2  = j*pframe->mol[pframe->mol_seq[i]]->nr_t + pframe->mol[pframe->mol_seq[i]]->b_seq[2*k+1]+pframe->n_start[i];
        r2   = pbc_dx(pframe,id1,id2,dx);
        dr   = r2*InvSqrt(r2);
        // BEGIN HARMONIC
        kc   = pframe->mol[pframe->mol_seq[i]]->k_bond[2*k];
        ddr  = dr-pframe->mol[pframe->mol_seq[i]]->k_bond[2*k+1];
        df   = kc*ddr;
        dV   = 0.5*df*ddr;
        // END HARMONIC
        E   += dV;
        df  *= InvSqrt(r2);
        for (m=XX; m<=ZZ; m++) {
          fij=-df*dx[m];
          pframe->f[id1][m]+=fij;
          pframe->f[id2][m]-=fij;
        }
      }

      // ANGLES
      for(k=0;k<pframe->mol[pframe->mol_seq[i]]->nr_a;k++){
        id1   = j*pframe->mol[pframe->mol_seq[i]]->nr_t + pframe->mol[pframe->mol_seq[i]]->a_seq[2*k  ]+pframe->n_start[i];
        id2   = j*pframe->mol[pframe->mol_seq[i]]->nr_t + pframe->mol[pframe->mol_seq[i]]->a_seq[2*k+1]+pframe->n_start[i];
        id3   = j*pframe->mol[pframe->mol_seq[i]]->nr_t + pframe->mol[pframe->mol_seq[i]]->a_seq[2*k+2]+pframe->n_start[i];
        ar21  = pbc_dx(pframe,id1,id2,dx);
        ar22  = pbc_dx(pframe,id3,id2,dy);
        phi   = calc_angle(dx,dy,&cosval);
        cosv2 = cosval*cosval;

        // BEGIN HARMONIC
        kc    = pframe->mol[pframe->mol_seq[i]]->k_angle[2*k];
        ddr   = phi-pframe->mol[pframe->mol_seq[i]]->k_angle[2*k+1]*DEG2RAD;
        df    = kc*ddr;
        dV    = 0.5*df*ddr;
        // END HARMONIC
        E    += dV;

        if(cosv2<1){
          double st,sth;
          double cik,cii,ckk;
          double nrkj2,nrij2;
          vector f_i,f_j,f_k;

          st    = -df*InvSqrt(1.0 - cosv2);
          sth   = st*cosval;
          nrkj2 = iprod(dy,dy);
          nrij2 = iprod(dx,dx);
          cik   = st*InvSqrt(nrkj2*nrij2);
          cii   = sth/nrij2;	
          ckk   = sth/nrkj2;	
      
          for (m=XX; (m<=ZZ); m++) {	
	    f_i[m]=-(cik*dy[m]-cii*dx[m]);
	    f_k[m]=-(cik*dx[m]-ckk*dy[m]);
	    f_j[m]=-f_i[m]-f_k[m];
	    pframe->f[id1][m]+=f_i[m];
	    pframe->f[id2][m]+=f_j[m];
	    pframe->f[id3][m]+=f_k[m];
          }
        }
      }
      // PROPER DIHEDRALS
      E += pdihs(pframe,i,j);

      // IMPROPER DIHEDRALS
      E += idihs(pframe,i,j);
    }
  }
  pframe->E[2]  = E;
  pframe->E[0] += E;

  return 0;
}

int force_calc_nb(t_frame pframe, t_files fpkg)
{
  int     i,j,k,N,ci,cj,nqi;
  vector  *x,*v,distx;
  double   dr,r,r2,r6,r12,ir,ir2,ir6,ir12,E=0.0,rc2;
  double   C6ij,C12ij,FORCE,f,df;
  double   q2,dfq,fq,fe;
  int      ui,uj,mi,mj,ip,jp;

  N   = pframe->nr_parts;
  pframe->virial = 0.0;
  rc2 = pframe->rc2;
  nqi = pframe->nr_unique;

  for(i=0;i<N;i++){
    for(j=i+1;j<N;j++){
      mi   = pframe->resnr [i];
      ui   = pframe->restyp[i];
      mj   = pframe->resnr [j];
      uj   = pframe->restyp[j];
      if( (mi==mj && ui==uj) )
        continue;
      r2   = pbc_dx(pframe,i,j,distx);
      if( r2>rc2  ) //if outside cutoff or same residue
        continue;
      ir2  = 1.0/r2;
      ci   = pframe->partyp_key[pframe->partyp[i]];
      cj   = pframe->partyp_key[pframe->partyp[j]];
      if(pframe->C6 [ci*nqi+cj] > 0 && pframe->C12[ci*nqi+cj]>0){
        ir6   = ir2*ir2*ir2;
        ir12  = ir6*ir6;
        C6ij = pframe->C6 [ci*nqi+cj]*ir6;
        C12ij= pframe->C12[ci*nqi+cj]*ir12;
        E   += C12ij-C6ij;
        df   = (12.0*C12ij - 6.0*C6ij);
        f    = df*ir2;
        for(k=XX;k<=ZZ;k++){
          FORCE             = f * distx[k] ;
          pframe->f[i][k]  += FORCE;
          pframe->f[j][k]  -= FORCE;
        }
      }
      if( pframe->Q [ci*nqi+cj] != 0.0 ){
        ir  = InvSqrt(r2);
        fe  = ELUNIT*pframe->Q[ci*nqi+cj]*ir; 
        E  +=  fe;
        f   = fe*ir2;
        for(k=XX;k<=ZZ;k++){
          FORCE             = f * distx[k] ;
          pframe->f[i][k]  += FORCE;
          pframe->f[j][k]  -= FORCE;
        }
      }
    }
  }

  if(!isfinite(E))
    fatal("NAN ENERGY");
  pframe->E[1]=E;

  return 0;
}

