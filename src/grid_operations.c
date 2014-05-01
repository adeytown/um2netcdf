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
     double tlat, tlon, degtorad, radtodeg, sock, cpart, t1, t2, longitude, latitude;
     double PI;

     PI = 3.1415926535898;

     degtorad = PI / 180.0;
     radtodeg = 180.0 / PI;
     if ( real_constants[4]==0 ) { sock = 0; }
     else                        { sock = real_constants[5] - PI; }

     for ( j=0; j<int_constants[6]; j++ ) {

         tlat = lat[j]*degtorad;
         t1   = -cos(real_constants[4])*sin(tlat);
         for ( i=0; i<int_constants[5]; i++ ) {

            tlon = lon[i]*degtorad;
            cpart = cos(tlon) * cos(tlat);
            t2    = sin(real_constants[4])*cpart;

            latitude = asin((cos(real_constants[4])*cpart)+(sin(real_constants[4])*sin(tlat)));
            if ( fabs(cos(latitude)+(t1+t2))<=1.0e-8 ) { longitude = PI; }
            else                                       { longitude = -acos((t1+t2)/cos(latitude)); }

            if ( tlon>=0.0 && tlon<=PI ) longitude = -1.0*longitude;
            longitude += sock;
            if ( longitude<0.0 ) { longitude += 2.0*PI; }

            lat[i+int_constants[5]*j] = (float )(latitude * radtodeg);
            lon[i+int_constants[5]*j] = (float )(longitude * radtodeg);

         }
     }

     return;
}

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
     for ( j=0; j<int_constants[6]; j++ ) {
     for ( i=0; i<int_constants[5]; i++ ) {
         val = real_constants[3] + real_constants[0]*((double ) i);
         ind = j*int_constants[5] + i;
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
     float *lat, *lon, val;

     lat = (float *) malloc( int_constants[5]*int_constants[6]*sizeof(float) );
     lon = (float *) malloc( int_constants[5]*int_constants[6]*sizeof(float) );

  /** Check if a rotated coordinate system is being used **/
     if ( header[3]<99 ) { construct_reg_lat_lon_arrays( lon, lat ); }
     else                { construct_rotated_lat_lon_arrays( lon, lat ); }

  /** Output latitude values to hard disk **/
     ierr = nc_inq_varid( ncid, "latitude", &varID );
     ierr = nc_put_var_float( ncid, varID, lon );

  /** Output longitude values to hard disk **/
     ierr = nc_inq_varid( ncid, "longitude", &varID );
     ierr = nc_put_var_float( ncid, varID, lat );

     free( lat );
     free( lon );

     return;
}
