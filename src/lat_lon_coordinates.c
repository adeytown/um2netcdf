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
#include <netcdf.h>
#include <math.h>

/***
 *** CONSTRUCT ROTATED LAT LON ARRAYS
 ***
 *** Subroutine where the longitude and latitutde arrays are constructed for the
 *** rotated lon/lat grid used in the UM model.  
 ***
 *** INPUT: 
 ***      lon -> 2D longitude array
 ***      lat -> 2D latitude array
 ***
 ***   Mark Cheeseman, NIWA
 ***   May 1, 2014
 ***/

void construct_rotated_lat_lon_arrays( float *lon, float *lat ) {

     int    i, j;
     double tlat, tlon, degtorad, sock, cpart, t1, t2, longitude, latitude;
     double pseudolat, pseudolon, cos_pseudolat, sin_pseudolat, factor;


  /** Convert lon/lat position of rotated pole from degrees to radians**/
     degtorad=3.1415926535898/180.0;
     pseudolat = real_constants[4] * degtorad;
     pseudolon = real_constants[5] * degtorad;

  /** Take cosine/sine functions of the newly converted rotated pole lat/lon **/
     cos_pseudolat = cos( pseudolat );
     sin_pseudolat = sin( pseudolat );

     if ( pseudolon==0 ) { sock = 0.0; }
     else                { sock = pseudolon - 3.1415926535898; }

  /** Determine latitude values on 'unrotated' grid **/
     for ( i=0; i<int_constants[5]; i++) {

       /*** Find lon along model column j ***/
         tlon = ( real_constants[3] + ((double ) i)*real_constants[0] )*degtorad;
         tlon = cos( tlon );

         for ( j=0; j<int_constants[6]; j++) {
             tlat = ( real_constants[2] + ((double ) j)*real_constants[1] )*degtorad;
             cpart = tlon*cos(tlat);
             latitude = asin( cos_pseudolat*cpart + sin_pseudolat*sin(tlat) );
             lat[i+int_constants[5]*j] = (float ) (latitude / degtorad);
         } 
     }

  /** Determine longitude values on 'unrotated' grid **/
     for ( i=0; i<int_constants[6]; i++) {

         factor = 1.0;
         tlon = ( real_constants[3] + ((double ) i)*real_constants[0] )*degtorad;
         if ( (tlon>=0.0) && (tlon<=3.1415926535898) ) { factor = -1.0; }
         tlon = cos( tlon );

         for ( j=0; j<int_constants[5]; j++) {
 
             tlat = ( real_constants[2] + ((double ) j)*real_constants[1] )*degtorad;
             cpart = tlon*cos(tlat);
             latitude = asin( cos_pseudolat*cpart + sin_pseudolat*sin(tlat) );

             t1 = -cos_pseudolat*sin(tlat);
             t2 = sin_pseudolat*cpart; 
             if ( fabs( cos(latitude)+t1+t2 )<=1.0e-8 ) { longitude = 3.1415926535898; }
             else                                      { longitude = -acos((t1+t2)/cos(latitude)); }

             longitude = factor*longitude + sock;
             if ( longitude<0.0 ) { longitude += 2.0*3.1415926535898; }

             lon[i+int_constants[6]*j] = (float ) (longitude / degtorad);

         } 
     }
     return;
}


/***
 *** CONSTRUCT REG LAT LON ARRAYS
 ***
 *** Subroutine that constructs the 2D longitude and latitude arrays.  Note 
 *** that each is just a straight copy of values in the appropriate rlon/rlat
 *** 1D arrays.  I.E. no curvilinear coordinate aspects are used.
 ***
 *** INPUT: 
 ***      lon -> 2D longitude array
 ***      lat -> 2D latitude array
 ***
 ***   Mark Cheeseman, NIWA
 ***   May 1, 2014
 ***/

void construct_reg_lat_lon_arrays( float *lon, float *lat ) {

     int    i, j, ind;
     double val;

  /** Construct the "model" latitude values **/
     for ( j=0; j<int_constants[6]; j++ ) {
         val = real_constants[2] + real_constants[1]*((double ) j);
         for ( i=0; i<int_constants[5]; i++ ) {
             ind = j*int_constants[5] + i;
             lat[ind] = (float ) val;
         }
     }

  /** Construct the "model" longitude values **/
     for ( j=0; j<int_constants[5]; j++ ) { 
     for ( i=0; i<int_constants[6]; i++ ) {
         ind = j*int_constants[6] + i;
         val = real_constants[3] + real_constants[0]*((double ) i);
         lon[ind] = (float ) val;
     }
     }

     return;
}


/***
 *** CONSTRUCT LAT LON ARRAYS 
 ***
 *** Subroutine where the 2D longitude and latitude arrays are determined
 *** and written to hard disk.  These are the values from the model grid
 *** in curvilinear coordinates. 
 ***
 ***   INPUT: 
 ***           ncid -> file handler for the newly created NetCDF file 
 ***
 ***   Mark Cheeseman, NIWA
 ***   December 17, 2013
 ***/ 

void construct_lat_lon_arrays( int ncid ) {

     int    varID, ierr;
     float *lat, *lon;

     lat = (float *) malloc( int_constants[5]*int_constants[6]*sizeof(float) );
     lon = (float *) malloc( int_constants[5]*int_constants[6]*sizeof(float) );

  /** Check if a rotated coordinate system is being used **/
     if ( header[3]<99 ) { construct_reg_lat_lon_arrays( lon, lat ); }
     else                { construct_rotated_lat_lon_arrays( lon, lat ); }

  /** Output latitude values to hard disk **/
     ierr = nc_inq_varid( ncid, "longitude", &varID );
     ierr = nc_put_var_float( ncid, varID, lon );

  /** Output longitude values to hard disk **/
     ierr = nc_inq_varid( ncid, "latitude", &varID );
     ierr = nc_put_var_float( ncid, varID, lat );

     free( lat );
     free( lon );
     exit(1);

     return;
}
