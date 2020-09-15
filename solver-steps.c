/* Function Prototypes */
double fx(double x, double y, double z, double hx, double hy, double hz);
double fy(double x, double y, double z,double hx, double hy, double hz);
double fz(double x, double y, double z,double hx, double hy, double hz);

double cx(double x, double y, double z,double ceta_x, double ceta_y, double ceta_z);
double cy(double x, double y, double z,double ceta_x, double ceta_y, double ceta_z);
double cz(double x, double y, double z,double ceta_x, double ceta_y, double ceta_z);

/* End Function Prototypes */

/*additional functions to calculate the crossproducts for the solver*/
/* DGL in x-direction */
double fx(double  x, double  y, double  z,double hx, double hy, double hz){
  double zz;
  zz = GAMMA/((1.0+ALPH*ALPH)*MU_S) * (- (y * hz - z* hy)-ALPH *(x * (x * hx + y * hy+z * hz)- hx));
  return zz;
}

/* DGL in y-direction */
double fy(double  x, double  y, double  z,double hx, double hy, double hz) {
  double zz;
  zz = GAMMA/((1.0+ALPH*ALPH)*MU_S) * (-(z * hx - x * hz)- ALPH * (y*(x*hx+y*hy+z*hz)-hy));
  return zz;
}

/* DGL in z-direction */
double fz(double  x, double  y, double  z,double hx, double hy, double hz) {
  double zz;
  zz = GAMMA/((1.0+ALPH*ALPH)*MU_S) * (- (x * hy - y * hx) - ALPH * (z*(x*hx+y*hy+z*hz)-hz));
  return zz;
}

/* DGL in x-direction (noise) */
double cx(double x, double y, double z, double ceta_x, double ceta_y, double ceta_z){
  double zz;
  zz = GAMMA/((1.0+ALPH*ALPH)*MU_S) *(- (y * ceta_z - z * ceta_y)- ALPH * (x*(x*ceta_x+y*ceta_y+z*ceta_z) - ceta_x));
  return zz;
}

/* DGL in y-direction (noise) */
double cy(double x, double y, double z,double ceta_x, double
ceta_y, double ceta_z) {
  double zz;
  zz = GAMMA/((1.0+ALPH*ALPH)*MU_S) * (-(z* ceta_x - x * ceta_z)- ALPH * (y*(x*ceta_x+y*ceta_y+z*ceta_z)-ceta_y));
  return zz;
}

/* DGL in z-direction (noise) */
double cz(double x, double y, double z, double ceta_x, double ceta_y, double ceta_z) {
  double zz;
  zz = GAMMA/((1.0+ALPH*ALPH)*MU_S) * (- (x * ceta_y - y * ceta_x)- ALPH * (z*(x*ceta_x+y*ceta_y+z*ceta_z)-ceta_z));
  return zz;

}
