#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "define.h"
#include "structs.h"
#include "units.h"
#include "tools.h"

double rnd1(){
	return ((float)rand())/((float)RAND_MAX);
}

double gauss(){
  return((rnd1()+rnd1()+rnd1()+rnd1()+rnd1()+rnd1()+rnd1()+rnd1()+rnd1()+rnd1()+rnd1()+rnd1()-6));
}

double InvSqrt (double y) 
// Inverted Square root is not more accurate than float
{
    float x;
    x = (float)y;
    {
      float xhalf = 0.5f*x;
      int i = *(int*)&x;

      i = 0x5f3759df - (i >> 1);
      x = *(float*)&i;
      x = x*(1.5f - xhalf*x*x);
      y = (double)x;
    }
    return y;
}

int quaternion(vector x, vector v, vector origo, double fi, vector nx)
// x is the old vector , v is the axis of rotation, origo is that, fi is the angle in radians, nx is the new vector
{
  double norm,q[4];
  vector xo;

  norm   = InvSqrt(v[XX]*v[XX]+v[YY]*v[YY]+v[ZZ]*v[ZZ]);
  xo[XX] = x[XX];
  xo[YY] = x[YY];
  xo[ZZ] = x[ZZ];

  q[0] = cos(fi*0.5);
  q[1] = v[0] * norm * sin(fi*0.5) ;
  q[2] = v[1] * norm * sin(fi*0.5) ;
  q[3] = v[2] * norm * sin(fi*0.5) ;

  nx[XX] = (q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3])*xo[XX] + (2*q[1]*q[2] - 2*q[0]*q[3])*xo[YY] + (2*q[1]*q[3] + 2*q[0]*q[2])*xo[ZZ]+origo[XX];
  nx[YY] = (2*q[1]*q[2] + 2*q[0]*q[3])*xo[XX] + (q[0]*q[0]-q[1]*q[1] + q[2]*q[2]-q[3]*q[3])*xo[YY] + (2*q[2]*q[3]-2*q[0]*q[1])*xo[ZZ]+origo[YY];
  nx[ZZ] = (2*q[1]*q[3] - 2*q[0]*q[2])*xo[XX] + (2*q[2]*q[3] + 2*q[0]*q[1])*xo[YY] + (q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3])*xo[ZZ]+origo[ZZ];

  return 0;
}

double pbc_dx(t_frame pframe, int i, int j, vector dist)
{
  dist[XX]=pframe->r[i][XX]-pframe->r[j][XX];
  dist[YY]=pframe->r[i][YY]-pframe->r[j][YY];
  dist[ZZ]=pframe->r[i][ZZ]-pframe->r[j][ZZ];

  dist[XX] -= pframe->box[XX][XX]*((int)nearbyint((dist[XX])/pframe->box[XX][XX]));
  dist[YY] -= pframe->box[YY][YY]*((int)nearbyint((dist[YY])/pframe->box[YY][YY]));
  dist[ZZ] -= pframe->box[ZZ][ZZ]*((int)nearbyint((dist[ZZ])/pframe->box[ZZ][ZZ]));

  return dist[XX]*dist[XX]+dist[YY]*dist[YY]+dist[ZZ]*dist[ZZ];
}


int print_molecule(t_molecule mol, t_files filepkg){
  int i;

  fprintf(filepkg->fil[0],"MOLECULE  :: %s\n",mol->name);
  fprintf(filepkg->fil[0],"MOLECULE  :: SEQ(%d) = {",mol->nr_t);
  for(i=0;i<mol->nr_t;i++)
    fprintf(filepkg->fil[0]," %d ",mol->t_seq[i]);
  fprintf(filepkg->fil[0],"}\n");
  fprintf(filepkg->fil[0],"MOLECULE  :: BON(%d) \n",mol->nr_b);
  for(i=0;i<mol->nr_b;i++)
    fprintf(filepkg->fil[0],"MOLECULE  :: BON(%d) = { %d %d %lf %lf }\n",i,mol->b_seq[2*i],mol->b_seq[2*i+1],mol->k_bond[2*i],mol->k_bond[2*i+1]);
  fprintf(filepkg->fil[0],"MOLECULE  :: ANG(%d) \n",mol->nr_a);
  for(i=0;i<mol->nr_a;i++)
    fprintf(filepkg->fil[0],"MOLECULE  :: ANG(%d) = { %d %d %d %lf %lf }\n",i,mol->a_seq[3*i],mol->a_seq[3*i+1],mol->a_seq[3*i+2],mol->k_angle[2*i],mol->k_angle[2*i+1]);
  fprintf(filepkg->fil[0],"MOLECULE  :: IMP(%d) \n",mol->nr_i);
  for(i=0;i<mol->nr_i;i++)
    fprintf(filepkg->fil[0],"MOLECULE  :: IMP(%d) = { %d %d %d %d %lf %lf %lf %lf  }\n",i,mol->i_seq[4*i],mol->i_seq[4*i+1],mol->i_seq[4*i+2],mol->i_seq[4*i+3],mol->k_improp[4*i],mol->k_improp[4*i+1],mol->k_improp[4*i+2],mol->k_improp[4*i+3]);
  fprintf(filepkg->fil[0],"MOLECULE  :: DIH(%d) \n",mol->nr_d);
  for(i=0;i<mol->nr_d;i++)
    fprintf(filepkg->fil[0],"MOLECULE  :: DIH(%d) = { %d %d %d %d %lf %lf %lf  }\n",i,mol->d_seq[4*i],mol->d_seq[4*i+1],mol->d_seq[4*i+2],mol->d_seq[4*i+3],mol->k_dihed[4*i],mol->k_dihed[4*i+1],mol->k_dihed[4*i+2]);

  return 0;
}

double calc_angle(vector x, vector y, double *cosval){
  double phi=0.0,arg,ry,rx,irx,iry;

  irx   = InvSqrt(x[XX]*x[XX]+x[YY]*x[YY]+x[ZZ]*x[ZZ]);
  iry   = InvSqrt(y[XX]*y[XX]+y[YY]*y[YY]+y[ZZ]*y[ZZ]);

  *cosval  = (x[XX]*y[XX]+x[YY]*y[YY]+x[ZZ]*y[ZZ])*irx*iry;
  phi      = acos( *cosval );

  return(phi);
}

double spec_angle(vector a, vector b)
{
    vector w;
    double wlen,s;
    
    wlen  = cprod(a,b,w);
    wlen *= InvSqrt(wlen);    
    s     = iprod(a,b);
    
    return atan2(wlen,s);
}

double dih_angle(int i, int j, int k, int l, t_frame pframe,
               vector r_ij,vector r_kj,vector r_kl,vector m,vector n,
               double *sign) //, int *t1, int *t2, int *t3)
{
  double ipr,phi;

  pbc_dx(pframe,i,j,r_ij);
  pbc_dx(pframe,k,j,r_kj);
  pbc_dx(pframe,k,l,r_kl);

  cprod(r_ij,r_kj,m); 			
  cprod(r_kj,r_kl,n);			
  phi=spec_angle(m,n); 			
  ipr=iprod(r_ij,n); 			
  (*sign)=(ipr<0.0)?-1.0:1.0;
  phi=(*sign)*phi; 			
					
  return phi;
}

int vector_inc(vector a,vector b)
{  
  a[XX]+=b[XX];
  a[YY]+=b[YY];
  a[ZZ]+=b[ZZ]; 
  return(0);
}

int vector_dec(vector a,vector b)
{
  a[XX]-=b[XX];
  a[YY]-=b[YY];
  a[ZZ]-=b[ZZ];
  return(0);
}

int vector_add(vector a,vector b,vector c)
{
  c[XX]=a[XX]+b[XX];
  c[YY]=a[YY]+b[YY];
  c[ZZ]=a[ZZ]+b[ZZ];
  return(0);  
}

int vector_sub(vector a,vector b,vector c)
{
  c[XX]=a[XX]-b[XX];
  c[YY]=a[YY]-b[YY];
  c[ZZ]=a[ZZ]-b[ZZ];
  return(0);  
}

int svmul(double a,vector v1,vector v2)
{
  v2[XX]=a*v1[XX];
  v2[YY]=a*v1[YY];
  v2[ZZ]=a*v1[ZZ];
  return(0);
}

int do_dih_fup(int i,int j,int k,int l,double ddphi,
		vector r_ij,vector r_kj,vector r_kl,
		vector m, vector n, t_frame pframe )
{
  vector f_i,f_j,f_k,f_l;
  vector uvec,vvec,svec,dx_jl;
  double iprm,iprn,nrkj,nrkj2;
  double a,p,q,toler;
  vector *f,*x;

  f=pframe->f;
  x=pframe->r;  

  iprm  = iprod(m,m);		
  iprn  = iprod(n,n);		
  nrkj2 = iprod(r_kj,r_kj);	
  toler = nrkj2*FLOAT_EPS;
  if ((iprm > toler) && (iprn > toler)) {
    nrkj  = nrkj2*InvSqrt(nrkj2);	
    a     = -ddphi*nrkj/iprm;	
    svmul(a,m,f_i);		
    a     = ddphi*nrkj/iprn;	
    svmul(a,n,f_l);		
    p     = iprod(r_ij,r_kj);	
    p    /= nrkj2;		
    q     = iprod(r_kl,r_kj);	
    q    /= nrkj2;		
    svmul(p,f_i,uvec);		
    svmul(q,f_l,vvec);		
    vector_sub(uvec,vvec,svec);	
    vector_sub(f_i,svec,f_j);	
    vector_add(f_l,svec,f_k);	
    vector_inc(f[i],f_i);   	
    vector_dec(f[j],f_j);   	
    vector_dec(f[k],f_k);   	
    vector_inc(f[l],f_l);   	
  }

  return(0);
}

double dopdihs(double cpA,double cpB ,double phiA,double phiB, int mult,
	       double phi,double lambda,double *V,double *F)
{
  double v,dvdl,mdphi,v1,sdphi,ddphi;
  double L1   = 1.0 - lambda;
  double ph0  = (L1*phiA + lambda*phiB)*DEG2RAD;
  double dph0 = (phiB - phiA)*DEG2RAD;
  double cp   = L1*cpA + lambda*cpB;
  
  mdphi = mult*phi - ph0;
  sdphi = sin(mdphi);
  ddphi = -cp*mult*sdphi;
  v1    = 1.0 + cos(mdphi);
  v     = cp*v1;
  
  dvdl  = (cpB - cpA)*v1 + cp*dph0*sdphi;
  
  *V = v;
  *F = ddphi;
  
  return dvdl;
}

double pdihs(t_frame pframe,int pi,int pj)
{
  int  i,type,ai,aj,ak,al;
  int  t1,t2,t3,k;
  vector r_ij,r_kj,r_kl,m,n;
  double phi,sign,ddphi,vpd,vtot,lambda=0.0;

  vtot = 0.0;

  for(k=0;k<pframe->mol[pframe->mol_seq[pi]]->nr_d;k++){

    ai   = pj*pframe->mol[pframe->mol_seq[pi]]->nr_t + pframe->mol[pframe->mol_seq[pi]]->d_seq[2*k  ]+pframe->n_start[pi];
    aj   = pj*pframe->mol[pframe->mol_seq[pi]]->nr_t + pframe->mol[pframe->mol_seq[pi]]->d_seq[2*k+1]+pframe->n_start[pi];
    ak   = pj*pframe->mol[pframe->mol_seq[pi]]->nr_t + pframe->mol[pframe->mol_seq[pi]]->d_seq[2*k+2]+pframe->n_start[pi];
    al   = pj*pframe->mol[pframe->mol_seq[pi]]->nr_t + pframe->mol[pframe->mol_seq[pi]]->d_seq[2*k+2]+pframe->n_start[pi];

    phi=dih_angle(ai,aj,ak,al,pframe,r_ij,r_kj,r_kl,m,n,&sign);

    dopdihs(pframe->mol[pframe->mol_seq[pi]]->k_dihed[2*k  ], 0.0, pframe->mol[pframe->mol_seq[pi]]->k_dihed[2*k+1], 0.0,
            pframe->mol[pframe->mol_seq[pi]]->k_dihed[2*k+2], phi, 0.0, &vpd, &ddphi);

    vtot += vpd;
    do_dih_fup(ai,aj,ak,al,ddphi,r_ij,r_kj,r_kl,m,n,pframe);
  }

  return vtot;
}

double idihs(t_frame pframe, int pi, int pj)
{
  int    i,type,ai,aj,ak,al;
  int    t1,t2,t3,k;
  double phi,phi0,dphi0,ddphi,sign,vtot;
  vector r_ij,r_kj,r_kl,m,n;
  double L1,kk,dp,dp2,kA,kB,pA,pB,dvdl,lambda=0.0;

  L1 = 1.0-lambda;
  dvdl = 0;

  vtot = 0.0;

  for(k=0;k<pframe->mol[pframe->mol_seq[pi]]->nr_d;k++){

    ai   = pj*pframe->mol[pframe->mol_seq[pi]]->nr_t + pframe->mol[pframe->mol_seq[pi]]->i_seq[2*k  ]+pframe->n_start[pi];
    aj   = pj*pframe->mol[pframe->mol_seq[pi]]->nr_t + pframe->mol[pframe->mol_seq[pi]]->i_seq[2*k+1]+pframe->n_start[pi];
    ak   = pj*pframe->mol[pframe->mol_seq[pi]]->nr_t + pframe->mol[pframe->mol_seq[pi]]->i_seq[2*k+2]+pframe->n_start[pi];
    al   = pj*pframe->mol[pframe->mol_seq[pi]]->nr_t + pframe->mol[pframe->mol_seq[pi]]->i_seq[2*k+2]+pframe->n_start[pi];

    phi  = dih_angle(ai,aj,ak,al,pframe,r_ij,r_kj,r_kl,m,n, &sign);
    
    kA = pframe->mol[pframe->mol_seq[pi]]->k_improp[2*k];
    kB = 0.0;
    pA = pframe->mol[pframe->mol_seq[pi]]->k_improp[2*k+1];
    pB = 0.0;

    kk    = L1*kA + lambda*kB;
    phi0  = (L1*pA + lambda*pB)*DEG2RAD;
    dphi0 = (pB - pA)*DEG2RAD;
    dp = phi-phi0;  
    if (dp >= M_PI)
      dp -= 2*M_PI;
    else if(dp < -M_PI)
      dp += 2*M_PI;
    
    dp2 = dp*dp;
    vtot += 0.5*kk*dp2;
    ddphi = -kk*dp;
    dvdl += 0.5*(kB - kA)*dp2 - kk*dphi0*dp;
    do_dih_fup(ai,aj,ak,al,(double)(-ddphi),r_ij,r_kj,r_kl,m,n,pframe);			
  }

  return vtot;
}

double iprod(vector x, vector y){
  return(x[XX]*y[XX]+x[YY]*y[YY]+x[ZZ]*y[ZZ]);
}

double cprod(vector x, vector y, vector out)
{
  out[XX]=x[YY]*y[ZZ]-x[ZZ]*y[YY];
  out[YY]=x[ZZ]*y[XX]-x[XX]*y[ZZ];
  out[ZZ]=x[XX]*y[YY]-x[YY]*y[XX];

  return(out[XX]*out[XX]+out[YY]*out[YY]+out[ZZ]*out[ZZ]);
}

double calc_energy(t_frame pframe, int I){ 
  // CALCULATE ENERGY OF PARTICLE OR RESIDUE I
  double v2,ekin=0,epot=0,rc2,r,r2,r6,r12,C6ij,C12ij,ir,ir2,ir6,ir12,fe;
  double e_tmp,unit,ndf,temp=1.0,scale;
  vector distx,dx,dy;
  int i,j,ci,cj,nqi;
  int ri,nrid,iq,k,kk,ui,uj,mj,ip,jp,pid;
  int id1,id2,id3,pjd,N,mi;
  double ar21,ar22,phi,cosv,cosv2,cosval,dr,dV;

  nqi  = pframe->nr_unique;
  rc2  = pframe->rc2;
  i    = I;
  mi   = pframe->resnr [i];
  ui   = pframe->restyp[i];
  ri   = pframe->resnr[I];
  ui   = pframe->restyp[I];
  pid  = pframe->mol_seq[ui];
  N    = pframe->nr_parts;

  for(j=0;j<N;j++){
    mj   = pframe->resnr [j];
    uj   = pframe->restyp[j];
    if( ( mi == mj && ui == uj ) )
      continue;
    if( i==j )
      fatal("KILL!");
    r2   = pbc_dx(pframe,i,j,distx);
    if( r2>rc2 )
      continue;
    ir2  = 1.0/r2;
    ci   = pframe->partyp_key[pframe->partyp[i]];
    cj   = pframe->partyp_key[pframe->partyp[j]];
    if(pframe->C6 [ci*nqi+cj] > 0 && pframe->C12[ci*nqi+cj]>0){
      ir6     = ir2*ir2*ir2;
      ir12    = ir6*ir6;
      C6ij    = pframe->C6 [ci*nqi+cj]*ir6;
      C12ij   = pframe->C12[ci*nqi+cj]*ir12;
      epot   += (C12ij-C6ij);
    }
    if( pframe->Q [ci*nqi+cj] != 0.0 ){
      ir      = InvSqrt(r2);
      fe      = ELUNIT*pframe->Q[ci*nqi+cj]*ir; 
      epot   +=  fe;
    }
  }

  return (epot);
}

int grand_canonical(t_frame pframe, t_files filepkg){
// PERFORM GC UPDATE
  int rc_val,i,nrcyc,k,l,m;
  int ri,rj,rk,rl,uid,pid,dN;
  double arg,etot,zz,vol,fi;
  vector translate,axis,origo,x,x0;

  zz     = exp(pframe->mu[2]/pframe->T)/(pframe->deBroglie_vol);
  nrcyc  = 1;
  vol    = (pframe->box[XX][XX]*pframe->box[YY][YY]*pframe->box[ZZ][ZZ]);
  uid    = pframe->nr_umols-1;
  pid    = pframe->mol_seq[uid];

  if(1){
    int rid,nrid,iq;
    for(k=0;k<nrcyc;k++){
      rc_val = rnd1()>0.5; 
      dN = pframe->nr_parts-pframe->n_start[uid];
      if( rc_val ){   // ANHILATION ( WE ONLY DO THIS WITH THE MOLECULE SPECIES THAT WAS ADDED LAST TO THE SYSTEM )
        fprintf(stderr,"-");
        if( dN <= pframe->mol[pid]->nr_t )
          return rc_val;
        i    = ceil((pframe->nr_parts-pframe->n_start[uid])*rnd1())-1;
        i   += pframe->n_start[uid];
        rid  = pframe->resnr[i];
        if( i < pframe->n_start[uid] || i >= pframe->nr_parts ){
          fatal("ALGORITHM FAIL");
        }
        etot = calc_energy(pframe,i)*1E3; // kJ mol^-1 -> J mol^-1
        arg  = exp(1.0*etot/pframe->T/RGAS)/(zz*vol)*dN;
        if( rnd1() < arg ){
          nrid = pframe->resnr[pframe->nr_parts-1];
          for( iq = 0; iq < pframe->mol[uid]->nr_t; iq++ ){
            pframe->r[pframe->n_start[uid]+rid*pframe->mol[pid]->nr_t+iq][XX]=pframe->r[pframe->n_start[uid]+nrid*pframe->mol[pid]->nr_t+iq][XX];
            pframe->r[pframe->n_start[uid]+rid*pframe->mol[pid]->nr_t+iq][YY]=pframe->r[pframe->n_start[uid]+nrid*pframe->mol[pid]->nr_t+iq][YY];
            pframe->r[pframe->n_start[uid]+rid*pframe->mol[pid]->nr_t+iq][ZZ]=pframe->r[pframe->n_start[uid]+nrid*pframe->mol[pid]->nr_t+iq][ZZ];
          }
          pframe->nr_parts-=pframe->mol[pid]->nr_t;
          pframe->nr_mols[uid]-=1;
          fprintf(stderr,"S");
          //update resnr
          for(rj=0;rj<pframe->nr_mols[uid];rj++)
            for(rk=0;rk<pframe->mol[pid]->nr_t;rk++)
              pframe->resnr[pframe->n_start[uid]+rj*pframe->mol[pid]->nr_t+rk]  = rj;
        }else{
          fprintf(stderr,"F");
        }
      }else{      // CREATION
        vector ivec;
        fprintf(stderr,"+");
        if( pframe->nr_parts + 100 > pframe->nr_alloc){
          renew_frame(pframe,200);
        }

        pframe->resnr[pframe->nr_parts]  = 1 + pframe->resnr[pframe->nr_parts-1];
        pframe->restyp[pframe->nr_parts] = pframe->resnr[pframe->nr_parts-1];

        rid  = pframe->resnr[pframe->nr_parts];
        i    = ceil((pframe->nr_parts-pframe->n_start[uid])*rnd1())-1;
        i   += pframe->n_start[uid];
        nrid  = pframe->resnr[i];

        if( i < pframe->n_start[uid] || i >= pframe->nr_parts ){
          fatal("ALGORITHM FAIL");
        }

        // BELOW IS FOR TRANSLATION
        translate[XX] = rnd1()*pframe->box[XX][XX];
        translate[YY] = rnd1()*pframe->box[YY][YY];
        translate[ZZ] = rnd1()*pframe->box[ZZ][ZZ];

        // BELOW IS FOR ROTATION
        axis[XX]=rnd1();
        axis[YY]=rnd1();
        axis[ZZ]=rnd1();
        fi = rnd1()*M_PI;

        l         = pframe->n_start[uid]+ rid*pframe->mol[pid]->nr_t;
        m         = pframe->n_start[uid]+nrid*pframe->mol[pid]->nr_t;

        origo[XX] = translate[XX];
        origo[YY] = translate[YY];
        origo[ZZ] = translate[ZZ];

        if(pframe->bIR){
          x0[XX] = pframe->ir[0][XX];
          x0[YY] = pframe->ir[0][YY];
          x0[ZZ] = pframe->ir[0][ZZ];
        }else{
          x0[XX] = pframe->r[m][XX];
          x0[YY] = pframe->r[m][YY];
          x0[ZZ] = pframe->r[m][ZZ];
        }

        for( iq = 0; iq < pframe->mol[pid]->nr_t; iq++ ){
          l=pframe->n_start[uid]+ rid*pframe->mol[pid]->nr_t+iq;
          m=pframe->n_start[uid]+nrid*pframe->mol[pid]->nr_t+iq;

          pframe->resnr [l] = pframe->resnr [pframe->nr_parts];
          pframe->restyp[l] = pframe->restyp[pframe->nr_parts];

          if(pframe->bIR){
            x[XX]             =  (pframe->ir[iq][XX] - x0[XX]) ;
            x[YY]             =  (pframe->ir[iq][YY] - x0[YY]) ;
            x[ZZ]             =  (pframe->ir[iq][ZZ] - x0[ZZ]) ;
          }else{
            x[XX]             =  (pframe->r[m][XX] - x0[XX]) ;
            x[YY]             =  (pframe->r[m][YY] - x0[YY]) ;
            x[ZZ]             =  (pframe->r[m][ZZ] - x0[ZZ]) ;
          }
// DO THE ROTATION
          quaternion(x,axis,origo,fi,pframe->r[l]);
// INIT VELOCITY
          pframe->v[l][XX]  = CUBE(gauss()); //pframe->v[m][XX];
          pframe->v[l][YY]  = CUBE(gauss()); //pframe->v[m][YY];
          pframe->v[l][ZZ]  = CUBE(gauss()); //pframe->v[m][ZZ];

          pframe->partyp[l]     = pframe->mol[pid]->t_seq[iq];
        }
        etot = calc_energy(pframe,pframe->nr_parts)*1E3;
        arg  = exp(-1.0*etot/pframe->T/RGAS)*(zz*vol)/(dN+1);
        if( rnd1() < arg ) {
          pframe->nr_parts+=pframe->mol[uid]->nr_t;
          pframe->nr_mols[uid]+=1;
          fprintf(stderr,"S");
        }else{
          fprintf(stderr,"F");
        }
      }
    }
  }

  return rc_val;
}

int renew_frame(t_frame pframe, int nr){
  // REALLOCATE MEMORY
  int i;

  pframe->r        = realloc(pframe->r     , sizeof(vector)*(pframe->nr_alloc+nr)  );
  pframe->v        = realloc(pframe->v     , sizeof(vector)*(pframe->nr_alloc+nr)*2);
  pframe->f        = realloc(pframe->f     , sizeof(vector)*(pframe->nr_alloc+nr)  );
  pframe->restyp   = realloc(pframe->restyp, sizeof(  int )*(pframe->nr_alloc+nr)  );
  pframe->resnr    = realloc(pframe->resnr , sizeof(  int )*(pframe->nr_alloc+nr)  );
  pframe->partyp   = realloc(pframe->partyp, sizeof(  int )*(pframe->nr_alloc+nr)  );
  pframe->nr_alloc = pframe->nr_parts;

  if(pframe->r     == NULL){
    fatal("REALLOC FAILED");
  }
  if(pframe->v     == NULL){
    fatal("REALLOC FAILED");
  }
  if(pframe->f     == NULL){
    fatal("REALLOC FAILED");
  }
  if(pframe->restyp== NULL){
    fatal("REALLOC FAILED");
  }
  if(pframe->resnr == NULL){
    fatal("REALLOC FAILED");
  }
  pframe->nr_alloc+=nr;

  for(i=0;i<nr;i++){
    pframe->partyp[i + pframe->nr_parts] = pframe->partyp[ pframe->nr_parts - 1];
    pframe->restyp[i + pframe->nr_parts] = pframe->restyp[ pframe->nr_parts - 1];
    pframe->resnr [i + pframe->nr_parts] = pframe->resnr [ pframe->nr_parts - 1];
  }

  return 0;
}

int check_frame(t_frame pframe, int step){
  int      j,i,ret_code=0,flag=0;
  int      ui,uj,mi,mj,ip,jp;

  if(!isfinite(pframe->box[XX][XX])){
    ret_code=-1;
    fprintf(stderr,"WARNING # S[%10d]::BOX IS NOT FINITE\n",step);
  }
  if(!isfinite(pframe->box[YY][YY])){
    ret_code=-1;
    fprintf(stderr,"WARNING # S[%10d]::BOX IS NOT FINITE\n",step);
  }
  if(!isfinite(pframe->box[ZZ][ZZ])){
    ret_code=-1;
    fprintf(stderr,"WARNING # S[%10d]::BOX IS NOT FINITE\n",step);
  }
  if(!isfinite(pframe->mu[0])){
    ret_code=-1;
    fprintf(stderr,"WARNING # S[%10d]::T-COUPLING IS NOT FINITE\n",step);
  }
  if(!isfinite(pframe->mu[1])){
    ret_code=-1;
    fprintf(stderr,"WARNING # S[%10d]::P-COUPLING IS NOT FINITE\n",step);
  }

  if(!(isfinite(pframe->pc2))){
    ret_code=-1;
    fprintf(stderr,"WARNING # S[%10d]:: DEBROGLIE WAVELENGTH IS BAD\n",step);
  }
  if(!isfinite(pframe->T)){
    ret_code=-1;
    fprintf(stderr,"WARNING # S[%10d]::TEMP IS NOT FINITE\n",step);
  }
  if(!isfinite(pframe->p)){
    ret_code=-1;
    fprintf(stderr,"WARNING # S[%10d]::PRESSURE IS NOT FINITE\n",step);
  }

  for(ui=0;ui<pframe->nr_umols;ui++){
      for(mi=0;mi<pframe->nr_mols[pframe->mol_seq[ui]];mi++){
          for(ip=0;ip<pframe->mol[pframe->mol_seq[ui]]->nr_t;ip++){
              i=mi*pframe->mol[pframe->mol_seq[ui]]->nr_t+ip+pframe->n_start[ui];
              if(!(pframe->partyp[i]==pframe->mol[pframe->mol_seq[ui]]->t_seq[ip])){
                ret_code=-1;
                fprintf(stderr,"WARNING # S[%10d]::residue missmatch\n",step);
              }
          }
      }
  }

  for(i=0;i<pframe->nr_parts;i++){
    if(!isfinite(pframe->r[i][XX]) && !flag){
      ret_code=-1;
      fprintf(stderr,"WARNING # S[%10d]::X COORDINATE IS NOT FINITE\n",step);
    }
    if(!isfinite(pframe->r[i][YY]) && !flag){
      ret_code=-1;
      fprintf(stderr,"WARNING # S[%10d]::Y COORDINATE IS NOT FINITE\n",step);
    }
    if(!isfinite(pframe->r[i][ZZ]) && !flag){
      ret_code=-1;
      fprintf(stderr,"WARNING # S[%10d]::Z COORDINATE IS NOT FINITE\n",step);
    }
    if(!isfinite(pframe->v[i][XX]) && !flag){
      ret_code=-1;
      fprintf(stderr,"WARNING # S[%10d]::X VELOCITY IS NOT FINITE\n",step);
    }
    if(!isfinite(pframe->v[i][YY]) && !flag){
      ret_code=-1;
      fprintf(stderr,"WARNING # S[%10d]::Y VELOCITY IS NOT FINITE\n",step);
    }
    if(!isfinite(pframe->v[i][ZZ]) && !flag){
      ret_code=-1;
      fprintf(stderr,"WARNING # S[%10d]::Z VELOCITY IS NOT FINITE\n",step);
    }
    if(!isfinite(pframe->f[i][XX]) && !flag){
      ret_code=-1;
      fprintf(stderr,"WARNING # S[%10d]::X FORCE IS NOT FINITE\n",step);
    }
    if(!isfinite(pframe->f[i][YY]) && !flag){
      ret_code=-1;
      fprintf(stderr,"WARNING # S[%10d]::Y FORCE IS NOT FINITE\n",step);
    }
    if(!isfinite(pframe->f[i][ZZ]) && !flag){
      ret_code=-1;
      fprintf(stderr,"WARNING # S[%10d]::Z FORCE IS NOT FINITE\n",step);
    }
    if(ret_code==-1)
      flag=1;
  }

  if(ret_code==-1)
    fatal("BROKEN SIMULATION!");

  return(ret_code);
}

int gen_velocities(t_frame pframe)
{
  int seed=4711,i,j,k;
  double Ekin=0.0,temp=0.0,scale=0.0,ndf=1.0,unit=1.0,Ekin_ref=0.0;

  srand( (unsigned)time( NULL ) );
  vector s;
  unit         = AMU*1E6;

  s[XX]=0.0; s[YY]=0.0; s[ZZ]=0.0;

  for(i=0;i<pframe->nr_parts;i++){
    k               = pframe->partyp_key[pframe->partyp[i]];
    pframe->v[i][XX]=gauss()/pframe->M[pframe->nr_unique*k+k];
    pframe->v[i][YY]=gauss()/pframe->M[pframe->nr_unique*k+k];
    pframe->v[i][ZZ]=gauss()/pframe->M[pframe->nr_unique*k+k];
    s[XX]+=pframe->v[i][XX];     s[YY]+=pframe->v[i][YY];     s[ZZ]+=pframe->v[i][ZZ];
  }
  s[XX]/=((double)pframe->nr_parts);
  s[YY]/=((double)pframe->nr_parts);
  s[ZZ]/=((double)pframe->nr_parts);

  for(i=0;i<pframe->nr_parts;i++){
    pframe->v[i][XX]-=s[XX];    pframe->v[i][YY]-=s[YY];    pframe->v[i][ZZ]-=s[ZZ];
    k                = pframe->partyp_key[pframe->partyp[i]];
    Ekin += 0.5*iprod(pframe->v[i],pframe->v[i])*pframe->M[pframe->nr_unique*k+k]*unit*AVOGADRO*1E-3;
  }
  ndf          = 3.0*(pframe->nr_parts-1); 
  Ekin_ref     = pframe->Tref*ndf*BOLTZ*0.5;
  temp         = 2.0*Ekin/ndf/BOLTZ; 
  pframe->T    = isfinite(temp)?temp:(1.0);
  scale        = InvSqrt(pframe->T/pframe->Tref); 
  pframe->E[0] = Ekin;

  for(i=0;i<pframe->nr_parts;i++){
    pframe->v[i][XX]*=scale;    pframe->v[i][YY]*=scale;    pframe->v[i][ZZ]*=scale;
  }

  if(0){
    fprintf(stderr,"\n%f %f %f %f | %f %f \n",pframe->E[0],ndf,Ekin_ref,pframe->deBroglie,Ekin,pframe->T);
  }

  return 0;
}

int write_data(t_frame pframe, t_files fpkg, double t)
{ // fpkg->fil[5]
  fprintf(fpkg->fil[5],"% lf % lf % lf % lf % lf % lf % lf % lf     % lf % lf     % d     % lf % lf % lf   \n",t,pframe->E[1],pframe->E[2],pframe->E[0],pframe->T, pframe->box[XX][XX],pframe->p,pframe->vpres,pframe->kpres,pframe->pcorr,pframe->nr_parts,pframe->mu[0],pframe->mu[1],pframe->deBroglie);
  return 0;
}
