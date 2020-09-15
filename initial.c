/* collection of functions to initialize the system for atomistic spin dynamics simulation */
/* Version from april 2020*/
/******************** Function Prototypes ********************************/
int Kubus(double sx[SIZE_X][SIZE_Y][SIZE_Z], double sy[SIZE_X][SIZE_Y][SIZE_Z], double sz[SIZE_X][SIZE_Y][SIZE_Z], int find_index[LX*LY*LZ][3]);
int Kubus_angle (double sx[SIZE_X][SIZE_Y][SIZE_Z], double sy[SIZE_X][SIZE_Y][SIZE_Z], double sz[SIZE_X][SIZE_Y][SIZE_Z], int find_index[LX*LY*LZ][3]);

/***** Initialize cube with parallel alignment and open boundary conditions *****/
int Kubus(double sx[SIZE_X][SIZE_Y][SIZE_Z], double sy[SIZE_X][SIZE_Y][SIZE_Z], double sz[SIZE_X][SIZE_Y][SIZE_Z], int find_index[LX*LY*LZ][3]) {
  int i, j, k, s = 0;

  for(i = 0; i < SIZE_X; i++)  {
    for(j = 0; j < SIZE_Y; j++) {
      for(k = 0; k < SIZE_Z; k++) {
				if((i == 0) || (i == LX+1) || (j == 0) || (j == LY+1) || (k == 0) || (k == LZ+1)){
	  			sx[i][j][k] = sy[i][j][k] = sz[i][j][k] = 0.0;
				}
				else {
	  			sx[i][j][k] = 1.0;
	  			sy[i][j][k] = 0.0;
	  			sz[i][j][k] = 0.0;
     	  	find_index[s][0] = i;
	  			find_index[s][1] = j;
	  			find_index[s][2] = k;
	  			s++;
				}
      }
    }
  }
  return(s);
}

/* Initialize cube with parallel alignment and all spins are aligned under a certain angle (here in x-direction) with amplitude AMP */
int Kubus_angle (double sx[SIZE_X][SIZE_Y][SIZE_Z], double sy[SIZE_X][SIZE_Y][SIZE_Z], double sz[SIZE_X][SIZE_Y][SIZE_Z], int find_index[LX*LY*LZ][3]){
int i, j, k, s = 0;
double betrag=sqrt(1-AMP*AMP);

  for(i = 0; i < SIZE_X; i++)  {
    for(j = 0; j < SIZE_Y; j++) {
      for(k = 0; k < SIZE_Z; k++) {
			if((i == 0) || (i == SIZE_X-1) || (j == 0) || (j == SIZE_Y-1) || (k == 0) || (k == SIZE_Z-1)){
	  			sx[i][j][k] = sy[i][j][k] = sz[i][j][k] = 0.0;
			}
			else {
	  			sx[i][j][k] = betrag;
	 				sy[i][j][k] = AMP;
	  			sz[i][j][k] = 0.0;
     	  	find_index[s][0] = i;
	  			find_index[s][1] = j;
	  			find_index[s][2] = k;
	  			s++;
			}
      }
    }
  }
  return(s);
}
