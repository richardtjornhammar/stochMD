#include <stdlib.h>
#include <math.h>

double rekin_noises(int n);
double rand1();
double gausdev();
double gdev(const int i);

double rekin(double k, double sigma, int ndeg, double taut){

  double factor,r;
  if(taut>0.1){
    factor=exp(-1.0/taut);
  } else{
    factor=0.0;
  }
  r = gausdev();
  return k + (1.0-factor)* (sigma*(rekin_noises(ndeg-1)+r*r)/ndeg-k)
            + 2.0*r*sqrt(k*sigma/ndeg*(1.0-factor)*factor);
}

double rekin_noises(int n){
/*
  sum n independent gaussians squared
*/
  double r;
  if(n==0) {
    return 0.0;
  } else if(n==1) {
    r=gausdev();
    return r*r;
  } else if(n%2==0) {
    return 2.0*gdev(n/2);
  } else {
    r=gausdev();
    return 2.0*gdev((n-1)/2)+r*r;
  }
}



double gdev(const int ia)
{
	int j;
	double am,e,s,v1,v2,x,y;

	if (ia < 1) {}; // FATAL ERROR
	if (ia < 6) {
		x=1.0;
		for (j=1;j<=ia;j++) x *= rand1();
		x = -log(x);
	} else {
		do {
			do {
				do {
					v1=rand1();
					v2=2.0*rand1()-1.0;
				} while (v1*v1+v2*v2 > 1.0);
				y=v2/v1;
				am=ia-1;
				s=sqrt(2.0*am+1.0);
				x=s*y+am;
			} while (x <= 0.0);
			e=(1.0+y*y)*exp(am*log(x/am)-s*y);
		} while (rand1() > e);
	}
	return x;
}


double gausdev()
{
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if (iset == 0) {
		do {
			v1=2.0*rand1()-1.0;
			v2=2.0*rand1()-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}


double rand1()
{
	return ((float)rand())/((float)RAND_MAX);
}
