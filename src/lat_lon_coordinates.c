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
 *** CONSTRUCT LON ARRAY
 ***
 *** Subroutine that constructs the 2D array giving the 'true longitude' value on earth 
 *** for each model X,Y point.
 ***
 *** INPUT/OUTPUT:  lon -> ptr to 2D array that will hold the true lon values
 ***
 ***   Mark Cheeseman, NIWA
 ***   May 2, 2014
 ***/

void construct_lon_array( float *lon ) {

     int    i, j, ind;
     double tlat, tlon, tol, degtorad, sock, cpart, t1, t2, longitude, latitude;
     double pseudolat, pseudolon, cos_pseudolat, sin_pseudolat, factor;

     if ( header[3]<99 ) {

     /** For an unrotated coordinate system, don't do anything.  Just copy the **/
     /** lon values into the correct position in the 2D array.                 **/

        for ( j=0; j<int_constants[6]; j++ ) {
        for ( i=0; i<int_constants[5]; i++ ) {
            ind = j*int_constants[5] + i;
            tlon = real_constants[3] + real_constants[0]*((double ) i);
            lon[ind] = (float ) tlon;
        }
        }

     } else {

  /** Convert lon/lat position of rotated pole from degrees to radians**/
        degtorad=3.1415926535898/180.0;
        pseudolat = real_constants[4] * degtorad;
        pseudolon = real_constants[5] * degtorad;

        for ( j=0; j<int_constants[6]; j++ ) {
        for ( i=0; i<int_constants[5]; i++ ) {
             tlat = (real_constants[2] + real_constants[1]*((double ) j))*degtorad;
             tlon = (real_constants[3] + real_constants[0]*((double ) i))*degtorad;

             sock = pseudolon - 3.1415926535898;
             tol = pseudolon*pseudolon;
             if ( tol<1.0e-20 ) { sock = 0; }

             cpart = cos(tlon) * cos(tlat);

             latitude = asin((cos(pseudolat)*cpart)+(sin(pseudolat)*sin(tlat)));
             t1 = -cos(pseudolat)*sin(tlat);
             t2 = sin(pseudolat)*cpart;
 
             tol = (cos(latitude)+(t1+t2)) * (cos(latitude)+(t1+t2));
             if ( tol<=1.0e-16 ) { longitude = 3.1415926535898; }
             else               { longitude = -acos((t1+t2)/cos(latitude)); }
    
             if ( (tlon>-1.0e-20) && (tlon<3.1415926535898001) ) { longitude = -1.0*longitude; }
             longitude += sock;

             if ( longitude<0.0 ) { longitude += 2.0*3.1415926535898; }
             lon[i+int_constants[5]*j] = (float ) (longitude / degtorad);
        }
        }

     }
     return;
}


/***
 *** CONSTRUCT LAT ARRAY
 ***
 *** Subroutine that constructs the 2D array giving the 'true latitude' value on earth 
 *** for each model X,Y point.
 ***
 *** INPUT/OUTPUT:  lat -> ptr to 2D array that will hold the true lat values
 ***
 ***   Mark Cheeseman, NIWA
 ***   May 2, 2014
 ***/

void construct_lat_array( float *lat ) {

     int    i, j, ind;
     double tlat, tlon, degtorad, sock, cpart, latitude;
     double pseudolat, pseudolon, cos_pseudolat, sin_pseudolat;

     if ( header[3]<99 ) {

     /** For an unrotated coordinate system, don't do anything.  Just copy the **/
     /** lat values into the correct position in the 2D array.                 **/

        for ( j=0; j<int_constants[6]; j++ ) {
            tlat = real_constants[2] + real_constants[1]*((double ) j);
            for ( i=0; i<int_constants[5]; i++ ) {
                ind = j*int_constants[5] + i;
                lat[ind] = (float ) tlat;
            }
        }

     } else {

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

     }
     return;
}

/***
 *** CONSTRUCT LON BOUNDS ARRAY
 ***
 *** Subroutine that constructs the 3D array giving the 'true longitude' value on earth 
 *** for the 4 sides of each rectangular grid cell at each model X,Y point.
 ***
 *** INPUT/OUTPUT:  lon -> ptr to 3D array that will hold the true lon values
 ***
 ***   Mark Cheeseman, NIWA
 ***   May 2, 2014
 ***/

void construct_lon_bounds_array( float *lon ) {

     int    i, j, k, ind;
     double tlat, tlon, tol, degtorad, sock, cpart, t1, t2, longitude, latitude;
     double pseudolat, pseudolon, cos_pseudolat, sin_pseudolat, factor;

  /** For an unrotated coordinate system, just add/subtract one half of a grid **/
  /** cell width/height to find the proper cell bound values.                  **/

     if ( header[3]<99 ) {
        for ( j=0; j<int_constants[6]; j++ ) {
        for ( i=0; i<int_constants[5]; i++ ) {

        /*** Find longitude of center of grid cell at (i,j) ***/
            tlon = real_constants[3] + real_constants[0]*((double ) i);

        /*** Find & store longitude of sides 1 and 2 ***/
            t1 = tlon - 0.5*real_constants[0];
            ind = j*int_constants[5] + i;
            lon[ind] = (float ) t1;
            ind += int_constants[5]*int_constants[6];
            lon[ind] = (float ) t1;

        /*** Find & store longitude of sides 3 and 4 ***/
            t1 = tlon + 0.5*real_constants[0];
            ind += int_constants[5]*int_constants[6];
            lon[ind] = (float ) t1;
            ind += int_constants[5]*int_constants[6];
            lon[ind] = (float ) t1;
        }
        }

     } else {

  /** Convert lon/lat position of rotated pole from degrees to radians**/
        degtorad=3.1415926535898/180.0;
        pseudolat = real_constants[4] * degtorad;
        pseudolon = real_constants[5] * degtorad;

  /** Take cosine/sine functions of the newly converted rotated pole lat/lon **/
        cos_pseudolat = cos( pseudolat );
        sin_pseudolat = sin( pseudolat );

        if ( pseudolon==0 ) { sock = 0.0; }
        else                { sock = pseudolon - 3.1415926535898; }

        for ( i=0; i<int_constants[5]; i++) {
        for ( j=0; j<int_constants[6]; j++) {
        for ( k=1; k<5; k++ ) {

           /** Determine lat & lon for the grid cell side k**/
            switch( k ){
                  case 1:
                       tlon = ( real_constants[3] + ((double ) i - 0.5)*real_constants[0] )*degtorad;
                       tlat = ( real_constants[2] + ((double ) j - 0.5)*real_constants[1] )*degtorad;
                       break;
                  case 2:
                       tlon = ( real_constants[3] + ((double ) i - 0.5)*real_constants[0] )*degtorad;
                       tlat = ( real_constants[2] + ((double ) j + 0.5)*real_constants[1] )*degtorad;
                       break;
                  case 3:
                       tlon = ( real_constants[3] + ((double ) i + 0.5)*real_constants[0] )*degtorad;
                       tlat = ( real_constants[2] + ((double ) j + 0.5)*real_constants[1] )*degtorad;
                       break;
                  case 4:
                       tlon = ( real_constants[3] + ((double ) i + 0.5)*real_constants[0] )*degtorad;
                       tlat = ( real_constants[2] + ((double ) j - 0.5)*real_constants[1] )*degtorad;
                       break;
            }

            factor = 1.0;
            if ( (tlon>-1.0e-20) && (tlon<3.14159265358001) ) { factor = -1.0; }
            tlon = cos( tlon );

            cpart = tlon*cos(tlat);
            latitude = asin( cos_pseudolat*cpart + sin_pseudolat*sin(tlat) );

            t1 = -cos_pseudolat*sin(tlat);
            t2 = sin_pseudolat*cpart;
            tol =  (cos(latitude)+t1+t2) * (cos(latitude)+t1+t2);
            if ( tol<=1.0e-16 ) { longitude = 3.1415926535898; }
            else                { longitude = -acos((t1+t2)/cos(latitude)); }

            longitude = factor*longitude + sock;
            if ( longitude<0.0 ) { longitude += 2.0*3.1415926535898; }

            ind = (k-1)*int_constants[6]*int_constants[5] + int_constants[5]*j + i;
            lon[ind] = (float ) (longitude / degtorad);

        }
        }
        }

     }
     return;
}


/***
 *** CONSTRUCT LAT BOUNDS ARRAY
 ***
 *** Subroutine that constructs the 3D array giving the 'true latitude' value on earth 
 *** for the 4 sides of each rectangular grid cell at each model X,Y point.
 ***
 *** INPUT/OUTPUT:  lat -> ptr to 3D array that will hold the true lat values
 ***
 ***   Mark Cheeseman, NIWA
 ***   May 2, 2014
 ***/

void construct_lat_bounds_array( float *lat ) {

     int    i, j, k, ind;
     double tlat, tlon, degtorad, sock, cpart, t1, t2, latitude;
     double pseudolat, pseudolon, cos_pseudolat, sin_pseudolat;

  /** For an unrotated coordinate system, just add/subtract one half of a grid **/
  /** cell width/height to find the proper cell bound values.                  **/

     if ( header[3]<99 ) {
        for ( j=0; j<int_constants[6]; j++ ) {

        /*** Find latitude at center of grid cell at (i,j) ***/
            tlat = real_constants[2] + real_constants[1]*((double ) j);
            t1 = tlat - 0.5*real_constants[1];
            t2 = tlat + 0.5*real_constants[1];

            for ( i=0; i<int_constants[5]; i++ ) {
            
            /*** Store latitude of side 1 ***/
                ind = j*int_constants[5] + i;
                lat[ind] = (float ) t1;
            
            /*** Store latitude of side 2 ***/
                ind += int_constants[5]*int_constants[6];
                lat[ind] = (float ) t2;
            
            /*** Store latitude of side 3 ***/
                ind += int_constants[5]*int_constants[6];
                lat[ind] = (float ) t2;
            
            /*** Store latitude of side 4 ***/
                ind += int_constants[5]*int_constants[6];
                lat[ind] = (float ) t1;
            }
        }

     } else {

  /** Convert lon/lat position of rotated pole from degrees to radians**/
        degtorad=3.1415926535898/180.0;
        pseudolat = real_constants[4] * degtorad;
        pseudolon = real_constants[5] * degtorad;

  /** Take cosine/sine functions of the newly converted rotated pole lat/lon **/
        cos_pseudolat = cos( pseudolat );
        sin_pseudolat = sin( pseudolat );

        if ( pseudolon==0 ) { sock = 0.0; }
        else                { sock = pseudolon - 3.1415926535898; }

        for ( i=0; i<int_constants[5]; i++) {
        for ( j=0; j<int_constants[6]; j++) {
        for ( k=1; k<5; k++ ) {

           /** Determine lat & lon for the grid cell side k**/
            switch( k ){
                  case 1:
                       tlon = ( real_constants[3] + ((double ) i - 0.5)*real_constants[0] )*degtorad;
                       tlat = ( real_constants[2] + ((double ) j - 0.5)*real_constants[1] )*degtorad;
                       break;
                  case 2:
                       tlon = ( real_constants[3] + ((double ) i - 0.5)*real_constants[0] )*degtorad;
                       tlat = ( real_constants[2] + ((double ) j + 0.5)*real_constants[1] )*degtorad;
                       break;
                  case 3:
                       tlon = ( real_constants[3] + ((double ) i + 0.5)*real_constants[0] )*degtorad;
                       tlat = ( real_constants[2] + ((double ) j + 0.5)*real_constants[1] )*degtorad;
                       break;
                  case 4:
                       tlon = ( real_constants[3] + ((double ) i + 0.5)*real_constants[0] )*degtorad;
                       tlat = ( real_constants[2] + ((double ) j - 0.5)*real_constants[1] )*degtorad;
                       break;
            }

            tlon = cos( tlon );
            cpart = tlon*cos(tlat);
            latitude = asin( cos_pseudolat*cpart + sin_pseudolat*sin(tlat) );

            ind = (k-1)*int_constants[6]*int_constants[5] + int_constants[5]*j + i;
            lat[ind] = (float ) (latitude / degtorad);

        }
        }
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
     float *buf;

  /** Compute/output the values for the 2D longitude & latitude UM variables **/
     buf = (float *) malloc( int_constants[5]*int_constants[6]*sizeof(float) );
     construct_lon_array( buf );

     ierr = nc_inq_varid( ncid, "longitude", &varID );
     ierr = nc_put_var_float( ncid, varID, buf );

     construct_lat_array( buf );

     ierr = nc_inq_varid( ncid, "latitude", &varID );
     ierr = nc_put_var_float( ncid, varID, buf );
     free( buf );

  /** Compute/output the values for the 3D longitude & latitude cell bounds UM variables **/
     buf = (float *) malloc( 4*int_constants[5]*int_constants[6]*sizeof(float) );
     construct_lon_bounds_array( buf );

     ierr = nc_inq_varid( ncid, "longitude_cell_bnd", &varID );
     ierr = nc_put_var_float( ncid, varID, buf );

     construct_lat_bounds_array( buf );

     ierr = nc_inq_varid( ncid, "latitude_cell_bnd", &varID );
     ierr = nc_put_var_float( ncid, varID, buf );

     free( buf );

     return;
}
