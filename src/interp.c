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

double *interp_do_nothing(  double *val, int nx, int ny ) {
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

double *u_to_p_point_interp_c_grid( double *val, int nx, int ny ) {

       int     i, j, index, ny_interp;
       double *interp_val;

  /** Is NY equal to the latitude dimension of the interpolated grid? **/
       if ( ny==int_constants[6] ) { ny_interp = ny; }
       else                        { ny_interp = ny+1; }

       interp_val = (double *)calloc( nx*ny_interp,sizeof(double) );

 /** Take the average of the U points above & below the desired P-point **/
       for ( i=0; i<nx; i++ ) {
       for ( j=ny-1; j>0; j-- ) {
           index = i + nx*j;
           interp_val[index] = 0.5*( val[index] + val[index-nx] ); 
       }
       }
       free( val );

       return interp_val;
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

double *v_to_p_point_interp_c_grid( double *val, int nx, int ny ) {

       int     i, j, index, ny_interp;
       double *interp_val;

  /** Is NY equal to the latitude dimension of the interpolated grid? **/
       if ( ny==int_constants[6] ) { ny_interp = ny; }
       else                        { ny_interp = ny+1; }

       interp_val = (double *)calloc( nx*ny_interp,sizeof(double) );

 /** Take the average of the V points to the left & right of the desired P-point **/
       for ( j=0; j<ny-1; j++ ) {
       for ( i=nx-1; i>1; i-- ) {
           index = i + nx*j;
           interp_val[index] = 0.5*( val[index] + val[index-1] ); 
       }
       }

 /** The top-most row must be replaced by the second top-most row. **/
       for ( i=0; i<nx; i++ ) {
           index = nx*(ny-1) + i;   
           interp_val[index] = val[index-nx];  
       }
       free( val );

       return interp_val;
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

double *b_to_c_grid_interp_u_points( double *val, int nx, int ny ) {

       int     i, j, index[4], ny_interp;
       double *interp_val;

  /** Is NY equal to the latitude dimension of the interpolated grid? **/
       if ( ny==int_constants[6] ) { ny_interp = ny; }
       else                        { ny_interp = ny+1; }

       interp_val = (double *)calloc( nx*ny_interp,sizeof(double) );

  /** Fill in the edges with un-interpolated values **/
       for ( i=0; i<nx; i++ ) {
           interp_val[i] = val[i+nx];
           interp_val[nx*ny-1-i] = val[nx*ny-1-i-nx];
       }
       for ( j=0; j<ny; j++ ) {
           interp_val[j*nx] = val[j*nx+1];
           interp_val[(j+1)*nx-1] = val[(j+1)*nx-2];
       }

  /** Take the average of the 4 horizontal points surrounding the desired P-point location **/
       for ( j=1; j<ny-1; j++ ) {
       for ( i=1; i<nx-1; i++ ) {
           index[0] = i + nx*j;
           index[1] = i + nx*j - 1;
           index[2] = i + nx*(j-1);
           index[3] = i + nx*(j-1) - 1;
           interp_val[index[0]] = 0.25*( val[index[0]] + val[index[1]] + val[index[2]] + val[index[3]]);
       }
       } 

       free( val );
       return interp_val;
}

