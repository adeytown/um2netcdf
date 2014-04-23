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
 *** Subroutine where the true longitude and latitude values for each UM field datapoint 
 *** are determined.  It is assumed that the data is located on the P-points of an
 *** Arakawa-C grid.
 ***
 ***   Mark Cheeseman, NIWA
 ***   December 17, 2013
 ***/

void construct_rotated_lat_lon_arrays( int ncid, double *lon, double *lat ) {

     int    i, j;
     double tlat, tlon, degtorad, radtodeg, sock, cpart, t1, t2, longitude, latitude;
     double PI, *lat_rotated, *lon_rotated;

     PI = 3.1415926535898;
     lon_rotated = (double *) malloc( int_constants[5]*int_constants[6]*sizeof(double) );
     lat_rotated = (double *) malloc( int_constants[5]*int_constants[6]*sizeof(double) );

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

            lat_rotated[i+int_constants[5]*j] = latitude * radtodeg;
            lon_rotated[i+int_constants[5]*j] = longitude * radtodeg;

         }
     } 

  /** Output the lon/lat arrays into the NetCDF file **/
     i = nc_inq_varid( ncid, "longitude", &j );
     i = nc_put_var_double( ncid, j, lon_rotated );

     i = nc_inq_varid( ncid, "latitude", &j );
     i = nc_put_var_double( ncid, j, lat_rotated );

     free( lat_rotated );
     free( lon_rotated );
     return;
}

void construct_rotated_lat_lon_arrays_float( int ncid, float *flon, float *flat ) {

     int    i, j;
     float tlat, tlon, degtorad, radtodeg, sock, cpart, t1, t2, longitude, 
           latitude, *flon_rotated, *flat_rotated;
     const float PI = 3.1415926535898;

     flon_rotated = (float *) malloc( int_constants[5]*int_constants[6]*sizeof(float) );
     flat_rotated = (float *) malloc( int_constants[5]*int_constants[6]*sizeof(float) );

     degtorad = PI / 180.0;
     radtodeg = 180.0 / PI;
     if ( real_constants[4]==0 ) { sock = 0; }
     else                        { sock = (float) real_constants[5] - PI; }

     for ( j=0; j<int_constants[6]; j++ ) {
     for ( i=0; i<int_constants[5]; i++ ) {

   /** Convert from degrees to radians **/
         tlat = flat[j]*degtorad;
         tlon = flon[i]*degtorad;

         cpart = cos(tlon) * cos(tlat);
         t1    = (float ) (-cos(real_constants[4])*sin(tlat));
         t2    = (float ) (sin(real_constants[4])*cpart);

         latitude = (float ) (asin((cos(real_constants[4])*cpart)+(sin(real_constants[4])*sin(tlat))));
         if ( fabs(cos(latitude)+(t1+t2))<=1.0e-8 ) { longitude = PI; }
         else                                       { longitude = -acos((t1+t2)/cos(latitude)); }

         if ( tlon>=0.0 && tlon<=PI ) longitude = -1.0*longitude;
         longitude += sock;
         if ( longitude<0.0 ) { longitude += 2.0*PI; }

         flat_rotated[i+int_constants[5]*j] = latitude * radtodeg;
         flon_rotated[i+int_constants[5]*j] = longitude * radtodeg;

     }
     }

  /** Output the lon/lat arrays into the NetCDF file **/
     i = nc_inq_varid( ncid, "longitude", &j );
     i = nc_put_var_float( ncid, j, flon_rotated );

     i = nc_inq_varid( ncid, "latitude", &j );
     i = nc_put_var_float( ncid, j, flat_rotated );

     free( flat_rotated );
     free( flon_rotated );
     return;
}


/***
 *** CONSTRUCT LAT LON ARRAYS 
 ***
 *** Subroutine where the model longitude and latitude of the UM datapoints 
 *** are determined.
 ***
 ***   INPUT:  rflag -> denotes whether reduced precision is being used 
 ***
 ***   Mark Cheeseman, NIWA
 ***   December 17, 2013
 ***/ 

void construct_lat_lon_arrays( int ncid, int rflag, int iflag ) {

     int     varID, i, j, ind;
     double *buf;
     float  *fbuf;

     if ( rflag==1 ) {

        fbuf = (float *) malloc( int_constants[5]*int_constants[6]*sizeof(float) );
        
        for ( j=0; j<int_constants[6]; j++ ) {
        for ( i=0; i<int_constants[5]; i++ ) {
            ind = j*int_constants[5] + i;
            fbuf[ind] = (float ) (real_constants[2] + j*real_constants[1]); 
        }
        }
     
        i = nc_inq_varid( ncid, "latitude", &varID );
        i = nc_put_var_float( ncid, varID, fbuf );

        for ( j=0; j<int_constants[6]; j++ ) {
        for ( i=0; i<int_constants[5]; i++ ) {
            ind = j*int_constants[5] + i;
            fbuf[ind] = (float ) (real_constants[3] + j*real_constants[0]); 
        }
        }
     
        i = nc_inq_varid( ncid, "longitude", &varID );
        i = nc_put_var_float( ncid, varID, fbuf );

        free( fbuf );

     } else {

        buf = (double *) malloc( int_constants[5]*int_constants[6]*sizeof(double) );
        
        for ( j=0; j<int_constants[6]; j++ ) {
        for ( i=0; i<int_constants[5]; i++ ) {
            ind = j*int_constants[5] + i;
            buf[ind] = real_constants[2] + j*real_constants[1]; 
        }
        }
     
        i = nc_inq_varid( ncid, "latitude", &varID );
        i = nc_put_var_double( ncid, varID, buf );

        for ( j=0; j<int_constants[6]; j++ ) {
        for ( i=0; i<int_constants[5]; i++ ) {
            ind = j*int_constants[5] + i;
            buf[ind] = real_constants[3] + j*real_constants[0]; 
        }
        }
     
        i = nc_inq_varid( ncid, "longitude", &varID );
        i = nc_put_var_double( ncid, varID, buf );

        free( buf );
     }

     return;
}
