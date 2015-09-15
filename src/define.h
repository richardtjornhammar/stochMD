#define SGN(X)    ((X)>=0?1.0:-1.0)
#define nm        0.000000001
#define k_b       1.3806503*nm*nm*nm*10000
#define e_q       1.60217646*nm*nm*0.1
#define epS       78.000000
#define Els       -2.837297
#define c_0       299792458.00
#define pi4       12.5663706143592
#define mu0       (pi4)*0.0000001
#define ep0       1/(c_0*c_0*mu0)
#define NA        6.02214179/nm/nm*100000
#define A2bohr    1.889725989
#define bohr2A    0.529177249

enum {
  etGRO=0, etPDB=1, etXYZ=2, etTOP=3, etITP=4, etDAT=5
};

enum {
  cmdMD=1,cmdCE=2,cmdMC=3
};

static struct ZNAME{
  char *name[120];
}atom={
  "LA","H ","He","Li","Be","B ","C ","N ","O ","F ","Ne","Na","Mg",
  "Al","Si","P ","S ","Cl","Ar","K ","Ca","Sc","Ti","V ","Cr","Mn",
  "Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr",
  "Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb",
  "Te","I ","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd",
  "Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W ","Re","Os","Ir",
  "Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",
  "Pa","U ","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
  "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg"};
