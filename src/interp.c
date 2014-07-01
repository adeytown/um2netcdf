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
 *** No interpolation is done to the field.  Only the scale factor is applied. 
 ***
 *** INPUT:
 ***      val       --> data array ( double precision )
 ***      fval      --> data array ( single precision )
 ***      var_index --> index denoting the current variable in the global array
 ***                    of stored UM variables
 ***
 ***   Mark Cheeseman, NIWA
 ***   May 28, 2014
 ***/ 

void interp_do_nothing(  double *val, float *fval, int var_index ) {

       int n, cnt;

       cnt = (int )( stored_um_vars[var_index].nx*stored_um_vars[var_index].ny );

       for ( n=0; n<cnt; n++ ) 
           fval[n] = stored_um_vars[var_index].scale_factor*((float ) val[n]); 

       return;
}


/***
 *** V_TO_P_POINT_INTERP_C_GRID 
 ***
 *** Subroutine that performs interpolation on a 2D array of points on 
 *** an Arakawa-C grid.  The data points are located on the V points
 *** and need to be translated onto the P points
 ***
 *** INPUT:
 ***      val        --> data array with its points located on V-points of a C-grid
 ***      var_index --> index denoting the current variable in the global array
 ***                    of stored UM variables
 ***
 ***   Mark Cheeseman, NIWA
 ***   May 28, 2014
 ***/

void v_to_p_point_interp_c_grid( double *val, float *fval, int var_index ) {

       int    NY, i, j, index, index2, index3;
       double factor, tmp;

       NY = (int ) stored_um_vars[var_index].ny;
       if ( int_constants[6]<NY ) { NY = int_constants[6]; }
       factor = 0.5*((double )stored_um_vars[var_index].scale_factor);


       /** Interpolate in the Y direction from rows 1 to NY-2 **/
          for ( j=1; j<NY-1; j++ ) {
          for ( i=0; i<stored_um_vars[var_index].nx; i++ ) {
              index = i + j*stored_um_vars[var_index].nx;
              index2= index - stored_um_vars[var_index].nx;
              index3= index + stored_um_vars[var_index].nx;
              tmp = factor*( val[index2] + val[index3] );
              fval[index3] = (float ) tmp;
          }
          }

       /** Copy contents of Row 1 into Row 0 **/
          for ( i=0; i<stored_um_vars[var_index].nx; i++ ) {
              index = i + stored_um_vars[var_index].nx;
              fval[i] = fval[index];
          }

       /** Copy contents of Row NY-2 into Rows NY to INT_CONSTANTS[6]-1 **/
          for ( j=NY; j<int_constants[6]; j++ ) {
          for ( i=0; i<stored_um_vars[var_index].nx; i++ ) {
              index = i + (NY-1)*stored_um_vars[var_index].nx;
              index2= i + j*stored_um_vars[var_index].nx;
              fval[index2] = fval[index];
          }
          }

       return;
}


/***
 *** U_TO_P_POINT_INTERP_C_GRID 
 ***
 *** Subroutine that performs interpolation on a 2D array of points on 
 *** an Arakawa-C grid.  The data points are located on the V points
 *** and need to be translated onto the P points
 ***
 *** INPUT:
 ***      val        --> data array with its points located on U-points of a C-grid
 ***      var_index --> index denoting the current variable in the global array
 ***                    of stored UM variables
 ***
 ***   Mark Cheeseman, NIWA
 ***   May 29, 2014
 ***/

void u_to_p_point_interp_c_grid( double *val, float *fval, int var_index ) {

       int    i, j, index, index2, NY;
       double factor, tmp;

       NY = (int ) stored_um_vars[var_index].ny;
       if ( int_constants[6]<NY ) { NY = int_constants[6]; }
       factor = 0.5*((double )stored_um_vars[var_index].scale_factor);


       /** Interpolate in the X direction from rows 1 to NY-2 **/
          for ( j=0; j<NY; j++ ) {
          for ( i=1; i<stored_um_vars[var_index].nx-1; i++ ) {
              index = i + j*stored_um_vars[var_index].nx;
              tmp = factor*(val[index-1] + fval[index+1]);
              fval[index] = (float ) tmp;
          }
          }
 
          for ( j=0; j<NY; j++ ) {
              index = j*stored_um_vars[var_index].nx;
              fval[index] = fval[index+1];
              index += stored_um_vars[var_index].nx-1;
              fval[index] = fval[index-1];
          }

          for ( j=NY; j<int_constants[6]; j++ ) {
          for ( i=0; i<stored_um_vars[var_index].nx; i++ ) {
              index = i + (NY-1)*stored_um_vars[var_index].nx;
              index2= i + j*stored_um_vars[var_index].nx;
              fval[index2] = fval[index];
          }
          }

       return;
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

void b_to_c_grid_interp_u_points( double *val, float *fval, int var_index ) {

     int     NY, i, j, index[5];
     double  factor, tmp;

     NY = (int ) stored_um_vars[var_index].ny;
     if ( int_constants[6]<NY ) { NY = int_constants[6]; }
     factor = (double ) (0.25*stored_um_vars[var_index].scale_factor); 


     /** Take the average of the 4 horizontal points surrounding the desired P-point location **/
        for ( j=1; j<NY-1; j++ ) {
        for ( i=1; i<stored_um_vars[var_index].nx-1; i++ ) {
            index[0] = j*stored_um_vars[var_index].nx + i; 
            index[1] = index[0] - 1; 
            index[2] = index[0] + 1; 
            index[3] = index[0] - stored_um_vars[var_index].nx; 
            index[4] = index[0] + stored_um_vars[var_index].nx; 
            tmp = factor*( val[index[0]] + val[index[1]] + val[index[2]] + val[index[3]] ); 
            fval[index[0]] = (float ) tmp; 
        }
        }

     /** Fill the missing columns [0 and NX-1] **/
        for ( j=1; j<NY-1; j++ ) {
            index[0] = j*stored_um_vars[var_index].nx;
            index[1] = index[0] + 1;
            fval[index[0]] = fval[index[1]];
            index[0] += stored_um_vars[var_index].nx-1;
            index[1] = index[0] - 1;
            fval[index[0]] = fval[index[1]];
        }

     /** Copy contents of row 1 into row 0 **/
        for ( i=0; i<stored_um_vars[var_index].nx-1; i++ ) {
            index[0] = i + stored_um_vars[var_index].nx;
            fval[i] = fval[index[0]];
        } 

     /** Copy contents of row NY-2 into rows NY-1 to INT_CONSTANTS[6] **/
        for ( j=NY-1; j<int_constants[6]; j++ ) {
        for ( i=0; i<stored_um_vars[var_index].nx; i++ ) {
            index[0] = i + stored_um_vars[var_index].nx*(NY-2);
            index[1] = i + stored_um_vars[var_index].nx*j;
            fval[index[1]] = fval[index[0]]; 
        }
        }

     return;
}

