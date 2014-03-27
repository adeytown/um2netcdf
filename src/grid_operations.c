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

     int    n, var_ids[2];
     double *lon, *lat;
     float  *flon, *flat;

  /*** Determine the  IDs of the lat and lon NetCDF variables ***/
     n = nc_inq_varid( ncid, "lon", &var_ids[0] );
     n = nc_inq_varid( ncid, "lat", &var_ids[1] );

  /*** Construct the lon/lat values ***/
     lon = (double *) malloc( int_constants[5]*sizeof(double) );
     lat = (double *) malloc( int_constants[6]*sizeof(double) );

     for ( n=0; n<int_constants[6]; n++ )
         lat[n] = real_constants[2] + n*real_constants[1];
      
     for ( n=0; n<int_constants[5]; n++ )
         lon[n] = real_constants[3] + n*real_constants[0];

  /*** Write the values into the appropriate NetCDF variables ***/
     if ( rflag==0 ) {
        n = nc_put_var_double( ncid, var_ids[0], lon );
        n = nc_put_var_double( ncid, var_ids[1], lat );
        if ( header[3]>99 ) { construct_rotated_lat_lon_arrays( ncid, lon, lat ); }
     } else {
        flon = (float *) malloc( int_constants[5]*sizeof(float) );
        for ( n=0; n<int_constants[5]; n++) { flon[n] = (float ) lon[n]; }
        n = nc_put_var_float( ncid, var_ids[0], flon );

        flat = (float *) malloc( int_constants[6]*sizeof(float) );
        for ( n=0; n<int_constants[6]; n++) { flat[n] = (float ) lat[n]; }
        n = nc_put_var_float( ncid, var_ids[1], flat );
        if ( header[3]>99 ) { construct_rotated_lat_lon_arrays_float( ncid, flon, flat ); }
        free( flat );
        free( flon );
     }
     free( lon );

  /*** Fill the second Lat variable if interpolation is not being used ***/
     if ( iflag==0 ) {
        n = nc_inq_varid( ncid, "rlat2", &var_ids[1] );

        if ( rflag==0 ) { 
           lon = (double *) malloc( (int_constants[6]-1)*sizeof(double) );
           for ( n=0; n<int_constants[6]-1; n++ )
               lon[n] = lat[n];
           n = nc_put_var_double( ncid, var_ids[1], lon ); 
           free( lon );
        } else { 
           flat = (float *) malloc( (int_constants[6]-1)*sizeof(float) );
           for ( n=0; n<int_constants[6]-1; n++) { flat[n] = (float ) lat[n]; }
           n = nc_put_var_float( ncid, var_ids[1], flat ); 
           free( flat );
        }
     }
     
     free( lat );
     return;
}
