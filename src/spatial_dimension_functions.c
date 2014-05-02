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
 *** SET_LON_LAT_DIMENSIONS
 ***
 *** Function that creates the horizontal dimensions (lon,lat) for the data  
 *** fields in the input UM fields file.  The dimensions are created
 *** in the new NetCDF file.
 ***
 *** INPUT:  ncid    -> file ID for the new NetCDF file
 ***         iflag   -> value of 1 indicates that interpolation is being used
 ***         rflag   -> value of 1 indicates that 32-bit output is to be used (64-bit otherwise) 
 ***
 ***    Mark Cheeseman, NIWA
 ***    December 19, 2013
 ***/

void set_lon_lat_dimensions( int ncid, int iflag, int rflag ) {

     int     ierr, *dimID, varID, dim_ids[4];
     nc_type vartype;
     float   tmp;

     ierr = nc_def_dim( ncid, "rlon", int_constants[5], &dim_ids[0] );
     ierr = nc_def_dim( ncid, "rlat", int_constants[6], &dim_ids[1] );

  /*** Y (Lat) dimension differs for some velocity variables which are defined  ***/
  /*** on an Arakawa-B grid. Not needed if interpolation is being used.         ***/
     if ( rflag==1 ) { vartype = NC_FLOAT; }
     else            { vartype = NC_DOUBLE; }

     dimID = (int *) malloc( sizeof(int) );
     if ( iflag==0 ) { 
        ierr = nc_def_dim( ncid, "rlat2", int_constants[6]-1, &dim_ids[2] ); 
        dimID[0] = dim_ids[2];
        ierr = nc_def_var( ncid, "rlat2", vartype, 1, dimID, &varID );
        ierr = nc_put_att_text( ncid, varID, "units", 7, "degrees" );
        ierr = nc_put_att_text( ncid, varID,  "axis",  1, "Y" );
        ierr = nc_put_att_text( ncid, varID, "standard_name", 13, "grid_latitude" );
        ierr = nc_put_att_text( ncid, varID, "long_name", 26, "latitude on a rotated grid" );
        ierr = nc_put_att_text( ncid, varID, "point_spacing", 4, "even" );
        ierr = nc_put_att_double( ncid, varID, "grid_north_pole_latitude",  NC_DOUBLE, 1, &real_constants[4] );
        ierr = nc_put_att_double( ncid, varID, "grid_north_pole_longitude", NC_DOUBLE, 1, &real_constants[5] );
     } 

  /*** Create NetCDF variables for each coordinate ***/
     dimID[0] = dim_ids[0];
     ierr = nc_def_var( ncid,  "rlon", NC_FLOAT, 1, dimID, &varID );
     ierr = nc_put_att_text( ncid, varID,         "units",  7, "degrees" );
     ierr = nc_put_att_text( ncid, varID,          "axis",  1, "X" );
     ierr = nc_put_att_text( ncid, varID, "standard_name", 14, "grid_longitude" );
     ierr = nc_put_att_text( ncid, varID, "point_spacing",  4, "even" );
     ierr = nc_put_att_text( ncid, varID,     "long_name", 27, "longitude on a rotated grid" );
     ierr = nc_put_att_double( ncid, varID,  "grid_north_pole_latitude", NC_DOUBLE, 1, &real_constants[4] );
     ierr = nc_put_att_double( ncid, varID, "grid_north_pole_longitude", NC_DOUBLE, 1, &real_constants[5] );

     dimID[0] = dim_ids[1];
     ierr = nc_def_var( ncid,  "rlat", NC_FLOAT, 1, dimID, &varID );
     ierr = nc_put_att_text( ncid, varID, "units", 7, "degrees" );
     ierr = nc_put_att_text( ncid, varID,  "axis",  1, "Y" );
     ierr = nc_put_att_text( ncid, varID, "standard_name", 13, "grid_latitude" );
     ierr = nc_put_att_text( ncid, varID, "long_name", 26, "latitude on a rotated grid" );
     ierr = nc_put_att_text( ncid, varID, "point_spacing", 4, "even" );
     ierr = nc_put_att_double( ncid, varID, "grid_north_pole_latitude",  NC_DOUBLE, 1, &real_constants[4] );
     ierr = nc_put_att_double( ncid, varID, "grid_north_pole_longitude", NC_DOUBLE, 1, &real_constants[5] );

  /*** Define dimensions for the lat/lon extents for cells on the coordinate grid ***/

     ierr = nc_def_dim( ncid, "lon_bnd", 4, &dim_ids[2] );
     ierr = nc_def_dim( ncid, "lat_bnd", 4, &dim_ids[3] );

  /*** Check if a rotated lat/lon coordinate system has been used in the UM ***/
  /*** fields file                                                          ***/

     if ( header[3]>99 ) {
        dimID[0] = dim_ids[1];

        ierr = nc_def_var( ncid, "rotated_pole", NC_CHAR, 1, dimID, &varID );
        ierr = nc_put_att_text( ncid, varID, "grid_mapping_name", 26, "rotated_latitude_longitude" );
        ierr = nc_put_att_double( ncid, varID, "grid_north_pole_latitude", NC_DOUBLE, 1, &real_constants[4] );
        ierr = nc_put_att_double( ncid, varID, "grid_north_pole_longitude", NC_DOUBLE, 1, &real_constants[5] );
     }
     free( dimID );

     dimID = (int *)malloc( 3*sizeof(int) );
     dimID[0] = dim_ids[2];
     dimID[1] = dim_ids[1];
     dimID[2] = dim_ids[0];

     ierr = nc_def_var( ncid, "longitude_cell_bnd", vartype, 3, dimID, &varID );
     ierr = nc_put_att_text( ncid, varID,  "long_name", 33, "longitude of cell bounds on earth" );
     ierr = nc_put_att_text(  ncid, varID,     "units", 13, "degrees_east" );

     dimID[0] = dim_ids[3];
     ierr = nc_def_var( ncid, "latitude_cell_bnd", vartype, 3, dimID, &varID );
     ierr = nc_put_att_text( ncid, varID, "long_name", 32, "latitude of cell bounds on earth" );
     ierr = nc_put_att_text( ncid, varID,     "units", 12, "degrees_north" );
     free( dimID );

     dimID = (int *)malloc( 2*sizeof(int) );
     dimID[0] = dim_ids[0];
     dimID[1] = dim_ids[1];

     ierr = nc_def_var( ncid,  "longitude", NC_FLOAT, 2, dimID, &varID );
     ierr = nc_put_att_text( ncid, varID, "standard_name", 9, "longitude" );
     ierr = nc_put_att_text( ncid, varID,     "long_name",18, "longitude on earth" );
     ierr = nc_put_att_text( ncid, varID,         "units",12, "degrees_east" );
     ierr = nc_put_att_text( ncid, varID,          "axis", 1, "X" );
     ierr = nc_put_att_text( ncid, varID,        "bounds",18, "longitude_cell_bnd" );

     ierr = nc_def_var( ncid,  "latitude", NC_FLOAT, 2, dimID, &varID );
     ierr = nc_put_att_text(  ncid, varID, "standard_name", 8, "latitude" );
     ierr = nc_put_att_text(  ncid, varID,     "long_name",17, "latitude on earth" );
     ierr = nc_put_att_text(  ncid, varID,         "units",13, "degrees_north" );
     ierr = nc_put_att_text(  ncid, varID,          "axis", 1, "Y" );
     ierr = nc_put_att_text( ncid, varID,         "bounds",17, "latitude_cell_bnd" );
     tmp = 90.0;
     ierr = nc_put_att_float( ncid, varID, "valid_max", NC_FLOAT, 1, &tmp );
     tmp = -90.0;
     ierr = nc_put_att_float( ncid, varID, "valid_min", NC_FLOAT, 1, &tmp );
      
     free( dimID );    

     return;
}


void set_horizontal_dimensions( int ncid ) {

     int    i, lon_id, lat_id;
     float *buf, dx;

    /*** Get ID of longitude & latitude vars ***/ 
     i = nc_inq_varid( ncid, "rlon", &lon_id );
     i = nc_inq_varid( ncid, "rlat", &lat_id );

    /*** Construct and fill an array to hold the longitude values ***/
     buf = malloc( int_constants[5]*sizeof(float) );

     buf[0] = real_constants[3];
     dx = (float ) real_constants[0];
     for ( i=1; i<int_constants[5]; i++ ) { buf[i] = buf[i-1] + dx; }

    /*** Write longitude values to hard disk ***/
     i = nc_enddef( ncid );
     i = nc_put_var( ncid, lon_id, buf );
     free( buf );

    /*** Construct and fill an array to hold the latitude values ***/
     buf = malloc( int_constants[6]*sizeof(float) );

     buf[0] = real_constants[2];
     dx = (float ) real_constants[1];
     for ( i=1; i<int_constants[6]; i++ ) { buf[i] = buf[i-1] + dx; }

    /*** Write longitude values to hard disk ***/
     i = nc_put_var( ncid, lat_id, buf );
     i = nc_redef( ncid );
     free( buf );
     return;
}
