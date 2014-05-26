/**============================================================================
                   U M 2 N e t C D F  V e r s i o n 2 . 0
                   --------------------------------------
 
       Main author: Mark Cheeseman
                    National Institute of Water & Atmospheric Research (Ltd)
                    Wellington, New Zealand
                    February 2014

   UM2NetCDF is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   UM2NetCDF is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License can be found in the main UM2NetCDF
   directory.  Alternatively, please see <http://www.gnu.org/licenses/>.
**============================================================================*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "field_def.h"


/***
 *** INTERP_DO_NOTHING
 ***
 *** Empty subroutine that does nothing to the input fata field.
 ***
 ***   Mark Cheeseman, NIWA
 ***   December 13, 2013
 ***/ 

double *interp_do_nothing(  double *val, int var_index ) {
       return val;
}


/***
 *** U_TO_P_POINT_INTERP_C_GRID 
 ***
 *** Subroutine that performs interpolation on a 2D array of points on 
 *** an Arakawa-C grid.  The data points are located on the U points
 *** and need to be translated onto the P points
 ***
 *** INPUT:
 ***      val -> data array with its points located on U-points of a C-grid
 ***       nx -> # of data points in X (longitude) direction
 ***       ny -> # of data points in Y (latitude) direction
 ***
 ***   Mark Cheeseman, NIWA
 ***   December 6, 2013
 ***/

double *u_to_p_point_interp_c_grid( double *val, int var_index ) {

       int     i, j, index[3], y_limit, ind[3];
       double *buf;

       j = stored_um_vars[var_index].nx*int_constants[6];
       buf = (double *) malloc( j*sizeof(double) );

   /** Find the appropriate Lat end limit for the interpolation calculations **/

       if ( int_constants[6]<=stored_um_vars[var_index].ny ) { y_limit=int_constants[6]; }
       else { y_limit = stored_um_vars[var_index].ny; }

   /** Take the average of the U points above & below the desired P-point **/

       for ( j=1; j<y_limit-1; j++ ) {
           ind[0] = stored_um_vars[var_index].nx*j;
           ind[1] = ind[0] - (int ) stored_um_vars[var_index].nx;
           ind[2] = ind[0] + (int ) stored_um_vars[var_index].nx;
           for ( i=0; i<stored_um_vars[var_index].nx; i++ ) {
               index[0] = i + ind[0];
               index[1] = i + ind[1];
               index[2] = i + ind[2];
               buf[index[0]] = 0.0;
               if ( (val[index[1]]>-1e-10)&&(val[index[1]]<1e10) ) { buf[index[0]] += 0.5*val[index[1]]; }
               if ( (val[index[2]]>-1e-10)&&(val[index[2]]<1e10) ) { buf[index[0]] += 0.5*val[index[2]]; }
           }
       }

       for ( i=0; i<stored_um_vars[var_index].nx; i++ ) { 
           index[0] = stored_um_vars[var_index].nx + i;
           buf[i] = buf[index[0]]; 
       }

       for ( j=y_limit-1; j<int_constants[6]; j++ ) {
       for ( i=0; i<stored_um_vars[var_index].nx; i++ ) {
           index[0] = i + stored_um_vars[var_index].nx*j;
           index[1] = i + stored_um_vars[var_index].nx*(j-1);
           buf[index[0]] = buf[index[1]]; 
       }
       }

   /** Deallocate the memory holding the uninterpolated data. Reassign pointer. **/

       free( val );
       return buf;
}


/***
 *** V_TO_P_POINT_INTERP_C_GRID 
 ***
 *** Subroutine that performs interpolation on a 2D array of points on 
 *** an Arakawa-C grid.  The data points are located on the V points
 *** and need to be translated onto the P points
 ***
 *** INPUT:
 ***      val -> data array with its points located on V-points of a C-grid
 ***       nx -> # of data points in X (longitude) direction
 ***       ny -> # of data points in Y (latitude) direction
 ***
 ***   Mark Cheeseman, NIWA
 ***   December 6, 2013
 ***/

double *v_to_p_point_interp_c_grid( double *val, int var_index ) {

       int     i, j, index[2], y_limit;
       double *buf;

       j = int_constants[6]*stored_um_vars[var_index].nx;
       buf = (double *) malloc( j*sizeof(double) );

   /** Find the appropriate Lat end limit for the interpolation calculations **/

       if ( int_constants[6]<=stored_um_vars[var_index].ny ) { y_limit=int_constants[6]; }
       else { y_limit = stored_um_vars[var_index].ny; }

 /** Take the average of the V points to the left & right of the desired P-point **/

       for ( j=0; j<y_limit; j++ ) {
       for ( i=1; i<stored_um_vars[var_index].nx-1; i++ ) {
           index[0] = j*stored_um_vars[var_index].nx + i;
           buf[index[0]] = 0.0;
           if ( (val[index[0]+1]>1.0e-10)&&(val[index[0]+1]<1.0e12) ) { buf[index[0]] += 0.5*val[index[0]+1]; }
           if ( (val[index[0]-1]>1.0e-10)&&(val[index[0]-1]<1.0e12) ) { buf[index[0]] += 0.5*val[index[0]-1]; }
       //    buf[index[0]] = 0.5*( val[index[0]+1] + val[index[0]-1] );
       }
       }

       for ( j=y_limit-1; j<int_constants[6]; j++ ) {
       for ( i=0; i<stored_um_vars[var_index].nx; i++ ) {
           index[0] = i + stored_um_vars[var_index].nx*j;
           index[1] = i + stored_um_vars[var_index].nx*(j-1);
           buf[index[0]] = buf[index[1]]; 
       }
       }

       for ( j=0; j<int_constants[6]; j++ ) {
           index[0] = j*stored_um_vars[var_index].nx;
           buf[index[0]] = buf[index[0]+1];
           index[1] = (j+1)*stored_um_vars[var_index].nx - 1;
           buf[index[1]] = buf[index[1]-1];
       }

       free( val );
       return buf;
}


/***
 *** B_TO_C_GRID_INTERP_U_POINTS 
 ***
 *** Subroutine that performs the interpolation of a 2D array of U-points on 
 *** an Arakawa-B grid to the corresponding P-points on an Arakawa-C grid.
 ***
 *** INPUT:
 ***      val -> data array with its points located on U-points of a B-grid
 ***       nx -> # of data points in X (longitude) direction of the original field
 ***       ny -> # of data points in Y (latitude) direction of the original field
 ***
 ***   Mark Cheeseman, NIWA
 ***   January 6, 2013
 ***/

double *b_to_c_grid_interp_u_points( double *val, int var_index ) {

     int     i, j, index[5], y_limit;
     double *buf;

     j = int_constants[6]*stored_um_vars[var_index].nx; 
     buf = (double *) calloc( j,sizeof(double) );

   /** Find the appropriate Lat end limit for the interpolation calculations **/

       if ( int_constants[6]<=stored_um_vars[var_index].ny ) { y_limit=int_constants[6]; }
       else { y_limit = stored_um_vars[var_index].ny; }

  /** Take the average of the 4 horizontal points surrounding the desired P-point location **/
  
     for ( j=1; j<y_limit-1; j++ ) {
     for ( i=1; i<stored_um_vars[var_index].nx-1; i++ ) {
         index[0] = j*stored_um_vars[var_index].nx + i; 
         index[1] = index[0] - 1; 
         index[2] = index[0] + 1; 
         index[3] = index[0] - stored_um_vars[var_index].nx; 
         index[4] = index[0] + stored_um_vars[var_index].nx; 
         buf[index[0]] = 0.0;
         if ( (val[index[0]]>1.0e-10)&&(val[index[0]]<1.0e12) ) { buf[index[0]] += 0.25*val[index[0]]; }
         if ( (val[index[1]]>1.0e-10)&&(val[index[1]]<1.0e12) ) { buf[index[0]] += 0.25*val[index[1]]; }
         if ( (val[index[2]]>1.0e-10)&&(val[index[2]]<1.0e12) ) { buf[index[0]] += 0.25*val[index[2]]; }
         if ( (val[index[3]]>1.0e-10)&&(val[index[3]]<1.0e12) ) { buf[index[0]] += 0.25*val[index[3]]; }
    //     buf[index[0]] = 0.25*( val[index[0]] + val[index[1]] + val[index[2]] + val[index[3]] ); 
     }
     }

  /** Fill the missing rows [0 and NY-1 -> INT_CONSTANTS(6)] **/

     for ( i=1; i<stored_um_vars[var_index].nx-1; i++ ) {
         index[0] = i + stored_um_vars[var_index].nx;
         buf[i] = buf[index[0]];
     } 

     for ( j=y_limit-1; j<int_constants[6]; j++ ) {
     for ( i=1; i<stored_um_vars[var_index].nx-1; i++ ) {
         index[0] = i + stored_um_vars[var_index].nx*j;
         index[1] = i + stored_um_vars[var_index].nx*(j-1);
         buf[index[0]] = buf[index[1]]; 
     }
     }

  /** Fill the missing columns [0 and NX-1] **/

     for ( j=0; j<int_constants[6]; j++ ) {
         index[0] = j*stored_um_vars[var_index].nx;
         index[1] = index[0] + 1;
         buf[index[0]] = buf[index[1]];
         index[2] = (j+1)*stored_um_vars[var_index].nx - 1;
         index[3] = index[2] - 1;
         buf[index[2]] = buf[index[3]];
     }

     free( val );
     return buf;
}

