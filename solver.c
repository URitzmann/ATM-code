
void heun(int spin, int n, int find_index[spin][3], double sx[SIZE_X][SIZE_Y][SIZE_Z], double sy[SIZE_X][SIZE_Y][SIZE_Z], double sz[SIZE_X][SIZE_Y][SIZE_Z],
double ceta_x[spin], double ceta_y[spin], double ceta_z[spin]);
double ex_sc(double scomp[SIZE_X][SIZE_Y][SIZE_Z], int i, int j, int k);



/* time integration */
void heun(int spin, int n, int find_index[spin][3], double sx[SIZE_X][SIZE_Y][SIZE_Z], double sy[SIZE_X][SIZE_Y][SIZE_Z], double sz[SIZE_X][SIZE_Y][SIZE_Z],
double ceta_x[spin], double ceta_y[spin], double ceta_z[spin]){
  static double sxs[SIZE_X][SIZE_Y][SIZE_Z]={0.0}, sys[SIZE_X][SIZE_Y][SIZE_Z]={0.0}, szs[SIZE_X][SIZE_Y][SIZE_Z]={0.0};
  static double sxn[SIZE_X][SIZE_Y][SIZE_Z]={0.0}, syn[SIZE_X][SIZE_Y][SIZE_Z]={0.0}, szn[SIZE_X][SIZE_Y][SIZE_Z]={0.0};
  double hx[spin],hy[spin],hz[spin];               /* effective field */
  double hx2,hy2,hz2;
  int i,j,k,l;
  double norm;								/* norm of the magnetic moment */
  double dt=DT;
  double q_z, ratio, amp1, amp2;

  const gsl_rng_type *rng_type;
  gsl_rng *rng;

  /* boundary conditions: open or periodic boundary condition of monochromatic spin wave excitation */
  /* default: open boundary conditions */
  #ifdef MONO_Z
   for(i=1; i<=LX;i++){
     for(j=1;j<=LY;j++){
       sx[i][j][0]=sqrt(1-AMP*AMP);
       sy[i][j][0]=AMP*cos(OMEGA*n*DT);
       sz[i][j][0]=AMP*sin(OMEGA*n*DT);
       sxs[i][j][0]=sqrt(1-AMP*AMP);
       sys[i][j][0]=AMP*cos(OMEGA*(n+1)*DT);
       szs[i][j][0]=AMP*sin(OMEGA*(n+1)*DT);
     }
   }
  #endif
  #ifdef PBC_Z
   for(i=1; i<=LX;i++){
     for(j=1;j<=LY;j++){
       sx[i][j][0]=sx[i][j][LZ];
       sy[i][j][0]=sy[i][j][LZ];
       sz[i][j][0]=sz[i][j][LZ];
       sx[i][j][LZ+1]=sx[i][j][1];
       sy[i][j][LZ+1]=sy[i][j][1];
       sz[i][j][LZ+1]=sz[i][j][1];
     }
   }
   #endif
   #ifdef PBC_Y
   for(i=1; i<=LX;i++){
     for(k=1;k<=LZ;k++){
       sx[i][0][k]=sx[i][LY][k];
       sy[i][0][k]=sy[i][LY][k];
       sz[i][0][k]=sz[i][LY][k];
       sx[i][LY+1][k]=sx[i][1][k];
       sy[i][LY+1][k]=sy[i][1][k];
       sz[i][LY+1][k]=sz[i][1][k];
     }
   }
   #endif
   #ifdef PBC_X
   for(j=1; j<=LY;j++){
     for(k=1;k<=LZ;k++){
       sx[0][j][k]=sx[LX][j][k];
       sy[0][j][k]=sy[LX][j][k];
       sz[0][j][k]=sz[LX][j][k];
       sx[LX+1][j][k]=sx[1][j][k];
       sy[LX+1][j][k]=sy[1][j][k];
       sz[LX+1][j][k]=sz[1][j][k];
     }
   }
   #endif

   for(l = 0; l < spin; l++){
       i = find_index[l][0]; j = find_index[l][1]; k = find_index[l][2];
       // effective field
       hx[l] = ex_sc(sx, i,j,k)+ 2*sx[i][j][k]*(DX);
       hy[l] = ex_sc(sy, i, j, k)+2*sy[i][j][k]*(DY);
       hz[l] = ex_sc(sz, i,j,k)+ 2*sz[i][j][k]*(DZ);

       /*Euler-step*/
       sxs[i][j][k] = sx[i][j][k] + fx(sx[i][j][k], sy[i][j][k], sz[i][j][k],hx[l], hy[l], hz[l]) * dt
                     + cx(sx[i][j][k], sy[i][j][k], sz[i][j][k], ceta_x[l], ceta_y[l], ceta_z[l]);
       sys[i][j][k] = sy[i][j][k] + fy(sx[i][j][k], sy[i][j][k], sz[i][j][k],hx[l], hy[l], hz[l]) * dt
                     + cy(sx[i][j][k], sy[i][j][k], sz[i][j][k], ceta_x[l], ceta_y[l], ceta_z[l]);
       szs[i][j][k] = sz[i][j][k] + fz(sx[i][j][k], sy[i][j][k], sz[i][j][k],hx[l], hy[l], hz[l]) * dt
                     + cz(sx[i][j][k], sy[i][j][k], sz[i][j][k], ceta_x[l], ceta_y[l], ceta_z[l]);
       norm = 1./sqrt(sxs[i][j][k]*sxs[i][j][k] + sys[i][j][k]*sys[i][j][k] + szs[i][j][k]*szs[i][j][k]);

       #ifdef KOMMENTAR
       //test if spin length is not wrong
       if ( norm > 1.01 || norm < 0.9999){
         fprintf(stderr,"Norm \n");
         fprintf(stderr,"Unnormalized Spin length is %f \n", norm);
         fprintf(stderr,"Unnormalized Spin length is %.16lf\t  %.16lf\t %.16lf \n", sxs[i][j][k],sys[i][j][k],szs[i][j][k]);
         fprintf(stderr,"Abort program!\n");
         exit(1);
       }
       #endif

       sxs[i][j][k] *= norm; sys[i][j][k] *= norm; szs[i][j][k] *= norm;
   }
   #ifdef PBC_Z
    for(i=1; i<=LX;i++){
      for(j=1;j<=LY;j++){
        sxs[i][j][0]=sxs[i][j][LZ];
        sys[i][j][0]=sys[i][j][LZ];
        szs[i][j][0]=szs[i][j][LZ];
        sxs[i][j][LZ+1]=sxs[i][j][1];
        sys[i][j][LZ+1]=sys[i][j][1];
        szs[i][j][LZ+1]=szs[i][j][1];
      }
    }
    #endif
    #ifdef PBC_Y
    for(i=1; i<=LX;i++){
      for(k=1;k<=LZ;k++){
        sxs[i][0][k]=sxs[i][LY][k];
        sys[i][0][k]=sys[i][LY][k];
        szs[i][0][k]=szs[i][LY][k];
        sxs[i][LY+1][k]=sxs[i][1][k];
        sys[i][LY+1][k]=sys[i][1][k];
        szs[i][LY+1][k]=szs[i][1][k];
      }
    }
    #endif

    #ifdef PBC_X
    for(j=1; j<=LY;j++){
      for(k=1;k<=LZ;k++){
        sxs[0][j][k]=sxs[LX][j][k];
        sys[0][j][k]=sys[LX][j][k];
        szs[0][j][k]=szs[LX][j][k];
        sxs[LX+1][j][k]=sxs[1][j][k];
        sys[LX+1][j][k]=sys[1][j][k];
        szs[LX+1][j][k]=szs[1][j][k];
      }
    }
   #endif
   //Heun-step
   for(l = 0; l < spin; l++){
       i = find_index[l][0]; j = find_index[l][1]; k = find_index[l][2];
       //effective field for (n+1)-step
       hx2 = ex_sc(sxs,i,j,k)+ 2*sxs[i][j][k]*(DX);
       hy2 = ex_sc(sys,i,j,k)+ 2*sys[i][j][k]*(DY);
       hz2 = ex_sc(szs,i,j,k)+ 2*szs[i][j][k]*(DZ);

       sxn[i][j][k] = sx[i][j][k] + dt*0.5*(fx(sx[i][j][k],  sy[i][j][k],  sz[i][j][k],hx[l],hy[l],hz[l]) +
                      fx(sxs[i][j][k], sys[i][j][k], szs[i][j][k],hx2,hy2,hz2))
                      + cx(sx[i][j][k], sy[i][j][k], sz[i][j][k], ceta_x[l], ceta_y[l], ceta_z[l]);
       syn[i][j][k] = sy[i][j][k] + dt*0.5*(fy(sx[i][j][k],  sy[i][j][k],  sz[i][j][k],hx[l],hy[l],hz[l]) +
                     fy(sxs[i][j][k], sys[i][j][k], szs[i][j][k],hx2,hy2,hz2))
                     + cy(sx[i][j][k], sy[i][j][k], sz[i][j][k], ceta_x[l], ceta_y[l], ceta_z[l]);
       szn[i][j][k] = sz[i][j][k] + dt*0.5*(fz(sx[i][j][k],  sy[i][j][k],  sz[i][j][k],hx[l],hy[l],hz[l]) +
                     fz(sxs[i][j][k], sys[i][j][k], szs[i][j][k],hx2,hy2,hz2))
                     + cz(sx[i][j][k], sy[i][j][k], sz[i][j][k], ceta_x[l], ceta_y[l], ceta_z[l]);
       norm = 1./sqrt(sxn[i][j][k]*sxn[i][j][k] + syn[i][j][k]*syn[i][j][k] + szn[i][j][k]*szn[i][j][k]);

 #ifdef KOMMENTAR
       //test if spin length is not wrong
       if ( norm > 1.01 || norm < 0.9999){
         fprintf(stderr,"Norm new spin\n");
         fprintf(stderr,"Unnormalized Spin length is %f \n", norm);
         fprintf(stderr,"Abort program!\n");
         exit(1);
       }
 #endif
   sxn[i][j][k] *= norm;
   syn[i][j][k] *= norm;
   szn[i][j][k] *= norm;
 }
 for(l = 0; l < spin; l++){
     i = find_index[l][0]; j = find_index[l][1]; k = find_index[l][2];
     sx[i][j][k] = sxn[i][j][k];
     sy[i][j][k] = syn[i][j][k];
     sz[i][j][k] = szn[i][j][k];
 }
}


/* exchange interaction*/
double ex_sc(double scomp[SIZE_X][SIZE_Y][SIZE_Z], int i, int j, int k){
	double zz;
	zz=J_1*(scomp[i+1][j][k]+scomp[i-1][j][k]) + J_1*(scomp[i][j+1][k]+scomp[i][j-1][k]) + J_1*(scomp[i][j][k+1]+scomp[i][j][k-1]);
	return zz;
}
