/* Structure */
#define LX 2               	/* lattice size in x-direction */
#define LY 2               	/* lattice size in y-direction */
#define LZ 256               	/* lattice size in z-direction */
#define RAND_PLUS 1         /*additional number of magnetic moments at each border (to secure exchange_interaction) */
#define SIZE_X  (LX+2*RAND_PLUS)    /*x-dimension for sx vectors */
#define SIZE_Y  (LY+2*RAND_PLUS)    /*y-dimension for sx vectors */
#define SIZE_Z  (LZ+2*RAND_PLUS)    /*z-dimension for sx vectors */
#define A  1.0      				/* lattice constant */

/*Flags */
#define KUBUS_SC
#define PBC_X
#define PBC_Y
#define MONO_Z
//#define PBC_Z

/* PARAMETER */
#define DT    (1.76E-3)       	/* time step */
#define NMAX  1000000            /* number of steps for output */
#define MU_S  1.0      		/* magn. saturartion magnetization mu_s */
#define GAMMA 1.              	/* gyroscopic ratio */
#define HX    (0.0*MU_S)      	/* magnetic field mu_s B in x-direction */
#define HY    (0.0*MU_S)      	/* magnetic field mu_s B in y-direction */
#define HZ    (0.0*MU_S)      	/* magnetic field mu_s B in z-direction */
#define J_1   1.0      		 /* 1st order exchange energy */
#define J_2   (0.1*J_1)     		 /* 1st order exchange energy */
#define J_3   (0.0*J_1)       		 /* 1st order exchange energy */
#define J_4   (0.0*J_1)     		 /* 1st order exchange energy */
#define J_5   (0.0*J_1)       		 /* 1st order exchange energy */
#define J_6   (0.0*J_1)       		 /* 1st order exchange energy */
#define J_7   (0.0*J_1)      		 /* 1st order exchange energy */
#define J_8   (0.0*J_1)      		 /* 1st order exchange energy */
#define J_9   (0.0*J_1)       		 /* 1st order exchange energy */
#define J_10  (0.0*J_1)     		 /* 1st order exchange energy */
#define J_IF  (0.0*J_1)   /* interface exchange energy*/
#define DX   (0.01*fabs(J_1))      		/* anistropy d V in x-direction (positive -> easy axis) */
#define DY   0.0      		/* anisotropy d V in y-direction */
#define DZ   (0.0*fabs(J_1))      		/* anisotropy d V in z-direction */
#define KB   1.0      		/* Boltzmann constant */
#define TEMP 0.0          /* Temperatur*/
#define ALPH  0.05          	/* damping constant */
#define OMEGA 0.5       /* precession frequency of the boundary conditions */
#define AMP  0.1		/* amplitude for initial alignment of the spin in function Kubus_angle */

/* Sonstiges */
#define PI 3.1415927
#define SEED 4711
