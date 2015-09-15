#define M_PI             3.14159265358979323846
#define ANGSTROM	 (1e-10)		  /* Old...	*/
#define KILO 		 (1e3)			  /* Thousand	*/
#define NANO		 (1e-9)			  /* A Number	*/
#define PICO		 (1e-12)		  /* A Number	*/
#define A2NM		 (ANGSTROM/NANO)	  /* NANO	        */
#define NM2A		 (NANO/ANGSTROM)	  /* 10.0		*/
#define RAD2DEG		 (180.0/M_PI)		  /* Conversion	*/
#define DEG2RAD		 (M_PI/180.0)		  /* id		*/
#define CAL2JOULE	 (4.184)		  /* id		*/
#define E_CHARGE         (1.60217733e-19)	  /* Coulomb	*/

#define AMU              (1.6605402e-27)          /* kg           */
#define BOLTZMANN	 (1.380658e-23)		  /* (J/K)	*/
#define AVOGADRO	 (6.0221367e23)		  /* ()		*/
#define RGAS             (BOLTZMANN*AVOGADRO)     /* (J/(mol K))  */
#define BOLTZ            (RGAS/KILO)              /* (kJ/(mol K)) */
#define FARADAY          (E_CHARGE*AVOGADRO)      /* (C/mol)      */
#define ELECTRONVOLT     (E_CHARGE*AVOGADRO/KILO) /* (kJ/mol)   */     
#define PLANCK1          (6.62606957e-34)         /* J s */
#define PLANCK           (6.6262e-34*AVOGADRO/(PICO*KILO)) /* (kJ/mol) ps */

#define EPSILON0 	 (5.72765E-4)		/* (e^2 / Na (kJ nm))     
						   == (e^2 mol/(kJ nm)) */
#define ELUNIT           (138.935437611)           //    kJ mol nm-1 e-2                                       
#define SPEED_OF_LIGHT   (2.9979245800E05)      /* nm/ps                */
#define SI_c             (299792458)            /* m/s */ 
#define SI_hc            (SI_c*PLANCK1)         /* hc */
#define ATOMICMASS_keV   (940000.0)             /* Atomic mass in keV   */
#define ELECTRONMASS_keV (512.0)                /* Electron mas in keV  */

#define FACEL		 (332.0636930*CAL2JOULE)/* (10 * (ONE_4PI_EPS0)) */
#define ONE_4PI_EPS0	 (FACEL*0.1)            /* 1/(4*pi*e0)*/
#define PRESFAC           (16.6054)             /* bar / pressure unity */
#define ENM2DEBYE         48.0321               /* Convert electron nm  *
						 * to debye             */
#define DEBYE2ENM         0.02081941
#define FIELDFAC          (FARADAY/KILO)

#define HARTREE2KJ        4.3597482e-21
#define BOHR2NM           0.0529177249
#define NM2BOHR           18.897259886
#define BOHR2AN           0.529177249 
#define HARTREE_BOHR2MD   (HARTREE2KJ*AVOGADRO/BOHR2NM)


