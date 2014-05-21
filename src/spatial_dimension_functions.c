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

void construct_lon_array( int ny, float *lon );
void construct_lat_array( int ny, float *lat );
void construct_lon_bounds_array( int ny, float *lon );
void construct_lat_bounds_array( int ny, float *lat );

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

     int     n, ierr, varID, dim_1d[1], i, dim_3d[3], 
             dim_2d[2], lon_bnd_dimid, lat_bnd_dimid;
     float   tmp, *buf;
     char    latname[8], lat2name[12], lonname[13], lonbndname[23], coord_str[30],
             latbndname[22];

  /*** Define dimensions for the lat/lon extents for cells on the coordinate grid ***/

     ierr = nc_def_dim( ncid, "lon_bnd", 4, &lon_bnd_dimid );
     ierr = nc_def_dim( ncid, "lat_bnd", 4, &lat_bnd_dimid );

  /*** Create a longitudinal dimension, Set each stored UM variable's lon dimension **/

     ierr = nc_def_dim( ncid, "rlon", int_constants[5], &dim_1d[0] );
     dim_2d[1] = dim_1d[0];
     dim_3d[2] = dim_1d[0];
     for ( n=0; n<num_stored_um_fields; n++ )
         stored_um_vars[n].x_dim = (unsigned short int ) dim_1d[0];

  /*** Create & fill the 1D longitudinal NetCDF variable  **/

     ierr = nc_def_var( ncid, "rlon", NC_FLOAT, 1, dim_1d, &varID );
     ierr = nc_put_att_text( ncid, varID,         "units",  7, "degrees" );
     ierr = nc_put_att_text( ncid, varID,          "axis",  1, "X" );
     ierr = nc_put_att_text( ncid, varID, "standard_name", 14, "grid_longitude" );
     ierr = nc_put_att_text( ncid, varID, "point_spacing",  4, "even" );
     ierr = nc_put_att_text( ncid, varID,     "long_name", 27, "longitude on a rotated grid" );
     tmp = (float ) real_constants[4];
     ierr = nc_put_att_float( ncid, varID, "grid_north_pole_latitude",  NC_FLOAT, 1, &tmp );
     tmp = (float ) real_constants[5];
     ierr = nc_put_att_float( ncid, varID, "grid_north_pole_longitude", NC_FLOAT, 1, &tmp );

     buf = (float *) malloc( int_constants[5]*sizeof(float) );
     buf[0] = (float ) real_constants[3];
     tmp = (float ) real_constants[0];
     for ( n=1; n<int_constants[5]; n++ ) { buf[n] = buf[n-1] + tmp; }

     ierr = nc_enddef( ncid );
     ierr = nc_put_var( ncid, varID, buf );
     ierr = nc_redef( ncid );
     free( buf );

  /*** Create a latitudinal dimension, Set each stored UM variable's lon dimension **/
   
     if ( iflag==0 ) {
    
     for ( n=0; n<num_stored_um_fields; n++ ) {
         sprintf( latname, "rlat%hu", stored_um_vars[n].ny );
         ierr = nc_inq_dimid( ncid, latname, &dim_1d[0] );
         if ( ierr!=NC_NOERR ) { 
            ierr = nc_def_dim( ncid, latname, (int )stored_um_vars[n].ny, &dim_1d[0] );
            ierr = nc_def_var( ncid, latname, NC_FLOAT, 1, dim_1d, &varID );
            ierr = nc_put_att_text( ncid, varID, "units", 7, "degrees" );
            ierr = nc_put_att_text( ncid, varID,  "axis",  1, "Y" );
            ierr = nc_put_att_text( ncid, varID, "standard_name", 13, "grid_latitude" );
            ierr = nc_put_att_text( ncid, varID, "long_name", 26, "latitude on a rotated grid" );
            ierr = nc_put_att_text( ncid, varID, "point_spacing", 4, "even" );
            tmp = (float ) real_constants[4];
            ierr = nc_put_att_float( ncid, varID, "grid_north_pole_latitude",  NC_FLOAT, 1, &tmp );
            tmp = (float ) real_constants[5];
            ierr = nc_put_att_float( ncid, varID, "grid_north_pole_longitude", NC_FLOAT, 1, &tmp );
         
            buf = (float *) malloc( (int )stored_um_vars[n].ny*sizeof(float) );
            buf[0] = (float ) real_constants[2];
            tmp = (float ) real_constants[1];
            for ( i=1; i<stored_um_vars[n].ny; i++ ) { buf[i] = buf[i-1] + tmp; }

            ierr = nc_enddef( ncid );
            ierr = nc_put_var( ncid, varID, buf );
            ierr = nc_redef( ncid );
            free( buf );

         dim_2d[0] = dim_1d[0];
         sprintf( lonname, "longitude%hu", stored_um_vars[n].ny );
         ierr = nc_def_var( ncid, lonname, NC_FLOAT, 2, dim_2d, &varID );
         ierr = nc_put_att_text( ncid, varID, "standard_name", 9, "longitude" );
         ierr = nc_put_att_text( ncid, varID,     "long_name",18, "longitude on earth" );
         ierr = nc_put_att_text( ncid, varID,         "units",12, "degrees_east" );
         ierr = nc_put_att_text( ncid, varID,          "axis", 1, "X" );
         sprintf( lonbndname, "longitude_cell_bnd%hu", stored_um_vars[n].ny );
         if ( stored_um_vars[n].ny<100 ) { i=20; }
         else if ( (stored_um_vars[n].ny>99)&&(stored_um_vars[n].ny<1000) ) { i=21; }
         else if ( (stored_um_vars[n].ny>999)&&(stored_um_vars[n].ny<10000) ) { i=22; }
         else { i=23; }
         ierr = nc_put_att_text( ncid, varID, "bounds", i, lonbndname );
         sprintf( coord_str, "latitude%hu longitude%hu", stored_um_vars[n].ny, stored_um_vars[n].ny );
         if ( stored_um_vars[n].ny<100 ) { i=22; }
         else if ( (stored_um_vars[n].ny>99)&&(stored_um_vars[n].ny<1000) ) { i=24; }
         else if ( (stored_um_vars[n].ny>999)&&(stored_um_vars[n].ny<10000) ) { i=26; }
         else { i=28; }
         ierr = nc_put_att_text( ncid, varID, "coordinates", i, coord_str );

         i = (int ) stored_um_vars[n].ny;
         buf = (float *) malloc( i*int_constants[5]*sizeof(float) );
         construct_lon_array( i, buf );
         ierr = nc_enddef( ncid );
         ierr = nc_put_var( ncid, varID, buf );
         ierr = nc_redef( ncid );
 
         sprintf( lat2name, "latitude%hu", stored_um_vars[n].ny );
         ierr = nc_def_var( ncid, lat2name, NC_FLOAT, 2, dim_2d, &varID );
         ierr = nc_put_att_text(  ncid, varID, "standard_name", 8, "latitude" );
         ierr = nc_put_att_text(  ncid, varID,     "long_name",17, "latitude on earth" );
         ierr = nc_put_att_text(  ncid, varID,         "units",13, "degrees_north" );
         ierr = nc_put_att_text(  ncid, varID,          "axis", 1, "Y" );
         sprintf( latbndname, "latitude_cell_bnd%hu", stored_um_vars[n].ny );
         if ( stored_um_vars[n].ny<100 ) { i=19; }
         else if ( (stored_um_vars[n].ny>99)&&(stored_um_vars[n].ny<1000) ) { i=20; }
         else if ( (stored_um_vars[n].ny>999)&&(stored_um_vars[n].ny<10000) ) { i=21; }
         else { i=22; }
         ierr = nc_put_att_text( ncid, varID, "bounds", i, latbndname );
         sprintf( coord_str, "latitude%hu longitude%hu", stored_um_vars[n].ny, stored_um_vars[n].ny );
         if ( stored_um_vars[n].ny<100 ) { i=22; }
         else if ( (stored_um_vars[n].ny>99)&&(stored_um_vars[n].ny<1000) ) { i=24; }
         else if ( (stored_um_vars[n].ny>999)&&(stored_um_vars[n].ny<10000) ) { i=26; }
         else { i=28; }
         ierr = nc_put_att_text( ncid, varID, "coordinates", i, coord_str );
         tmp = 90.0;
         ierr = nc_put_att_float( ncid, varID, "valid_max", NC_FLOAT, 1, &tmp );
         tmp = -90.0;
         ierr = nc_put_att_float( ncid, varID, "valid_min", NC_FLOAT, 1, &tmp );

         i = (int ) stored_um_vars[n].ny;
         construct_lat_array( i, buf );
         ierr = nc_enddef( ncid );
         ierr = nc_put_var( ncid, varID, buf );
         ierr = nc_redef( ncid ); 
         free( buf );
        
         dim_3d[0] = lon_bnd_dimid;
         dim_3d[1] = dim_1d[0];
         ierr = nc_def_var( ncid, lonbndname, NC_FLOAT, 3, dim_3d, &varID );
         ierr = nc_put_att_text( ncid, varID, "long_name", 33, "longitude of cell bounds on earth" );
         ierr = nc_put_att_text(  ncid, varID,    "units", 12, "degrees_east" );

         i = (int ) stored_um_vars[n].ny;
         buf = (float *) malloc( 4*i*int_constants[5]*sizeof(float) );
         construct_lon_bounds_array( i, buf );
         ierr = nc_enddef( ncid );
         ierr = nc_put_var( ncid, varID, buf );
         ierr = nc_redef( ncid );
 
         dim_3d[0] = lat_bnd_dimid;
         ierr = nc_def_var( ncid, latbndname, NC_FLOAT, 3, dim_3d, &varID );
         ierr = nc_put_att_text( ncid, varID,"long_name", 32, "latitude of cell bounds on earth" );
         ierr = nc_put_att_text( ncid, varID,    "units", 13, "degrees_north" );

         construct_lat_bounds_array( i, buf );
         ierr = nc_enddef( ncid );
         ierr = nc_put_var( ncid, varID, buf );
         ierr = nc_redef( ncid );
         free( buf );
         }
         stored_um_vars[n].y_dim = (unsigned short int ) dim_1d[0]; 
     }

     } else {

       ierr = nc_def_dim( ncid, "rlat", (int ) int_constants[6], &dim_1d[0] );
       ierr = nc_def_var( ncid, "rlat", NC_FLOAT, 1, dim_1d, &varID );
       ierr = nc_put_att_text( ncid, varID, "units", 7, "degrees" );
       ierr = nc_put_att_text( ncid, varID,  "axis",  1, "Y" );
       ierr = nc_put_att_text( ncid, varID, "standard_name", 13, "grid_latitude" );
       ierr = nc_put_att_text( ncid, varID, "long_name", 26, "latitude on a rotated grid" );
       ierr = nc_put_att_text( ncid, varID, "point_spacing", 4, "even" );
       tmp = (float ) real_constants[4];
       ierr = nc_put_att_float( ncid, varID, "grid_north_pole_latitude",  NC_FLOAT, 1, &tmp );
       tmp = (float ) real_constants[5];
       ierr = nc_put_att_float( ncid, varID, "grid_north_pole_longitude", NC_FLOAT, 1, &tmp );
         
       buf = (float *) malloc( (int )int_constants[6]*sizeof(float) );
       buf[0] = (float ) real_constants[2];
       tmp = (float ) real_constants[1];
       for ( i=1; i<int_constants[6]; i++ ) { buf[i] = buf[i-1] + tmp; }

       ierr = nc_enddef( ncid );
       ierr = nc_put_var( ncid, varID, buf );
       ierr = nc_redef( ncid );
       free( buf );

       for ( n=0; n<num_stored_um_fields; n++ )
           stored_um_vars[n].y_dim = (unsigned short int ) dim_1d[0]; 

       dim_2d[0] = dim_1d[0];
       ierr = nc_def_var( ncid, "longitude", NC_FLOAT, 2, dim_2d, &varID );
       ierr = nc_put_att_text( ncid, varID, "standard_name", 9, "longitude" );
       ierr = nc_put_att_text( ncid, varID,     "long_name",18, "longitude on earth" );
       ierr = nc_put_att_text( ncid, varID,         "units",12, "degrees_east" );
       ierr = nc_put_att_text( ncid, varID,          "axis", 1, "X" );
       ierr = nc_put_att_text( ncid, varID, "bounds", 18, "longitude_cell_bnd" );
       ierr = nc_put_att_text( ncid, varID, "coordinates", 18, "latitude longitude" );

       i = (int ) int_constants[6];
       buf = (float *) malloc( i*int_constants[5]*sizeof(float) );
       construct_lon_array( i, buf );
       ierr = nc_enddef( ncid );
       ierr = nc_put_var( ncid, varID, buf );
       ierr = nc_redef( ncid );

       ierr = nc_def_var( ncid, "latitude", NC_FLOAT, 2, dim_2d, &varID );
       ierr = nc_put_att_text(  ncid, varID, "standard_name", 8, "latitude" );
       ierr = nc_put_att_text(  ncid, varID,     "long_name",17, "latitude on earth" );
       ierr = nc_put_att_text(  ncid, varID,         "units",13, "degrees_north" );
       ierr = nc_put_att_text(  ncid, varID,          "axis", 1, "Y" );
       ierr = nc_put_att_text( ncid, varID, "bounds", 17, "latitude_cell_bnd" );
       ierr = nc_put_att_text( ncid, varID, "coordinates", 18, "latitude longitude" );
       tmp = 90.0;
       ierr = nc_put_att_float( ncid, varID, "valid_max", NC_FLOAT, 1, &tmp );
       tmp = -90.0;
       ierr = nc_put_att_float( ncid, varID, "valid_min", NC_FLOAT, 1, &tmp );

       construct_lat_array( i, buf );
       ierr = nc_enddef( ncid );
       ierr = nc_put_var( ncid, varID, buf );
       ierr = nc_redef( ncid );
       free( buf );

       dim_3d[0] = lon_bnd_dimid;
       dim_3d[1] = dim_1d[0];
       ierr = nc_def_var( ncid, "longitude_cell_bnd", NC_FLOAT, 3, dim_3d, &varID );
       ierr = nc_put_att_text( ncid, varID, "long_name", 33, "longitude of cell bounds on earth" );
       ierr = nc_put_att_text(  ncid, varID,    "units", 12, "degrees_east" );

       buf = (float *) malloc( 4*i*int_constants[5]*sizeof(float) );
       construct_lon_bounds_array( i, buf );
       ierr = nc_enddef( ncid );
       ierr = nc_put_var( ncid, varID, buf );
       ierr = nc_redef( ncid );

       dim_3d[0] = lat_bnd_dimid;
       ierr = nc_def_var( ncid, "latitude_cell_bnd", NC_FLOAT, 3, dim_3d, &varID );
       ierr = nc_put_att_text( ncid, varID,"long_name", 32, "latitude of cell bounds on earth" );
       ierr = nc_put_att_text( ncid, varID,    "units", 13, "degrees_north" );

       construct_lat_bounds_array( i, buf );
       ierr = nc_enddef( ncid );
       ierr = nc_put_var( ncid, varID, buf );
       ierr = nc_redef( ncid );
       free( buf );
     }

  /*** If a rotated lon/lat grid is being used, create a rotated pole NetCDF variable ***/

     if ( header[3]>99 ) {
        ierr = nc_def_var( ncid, "rotated_pole", NC_CHAR, 1, dim_1d, &varID );
        ierr = nc_put_att_text( ncid, varID, "grid_mapping_name", 26, "rotated_latitude_longitude" );
        tmp = (float ) real_constants[4];
        ierr = nc_put_att_float( ncid, varID, "grid_north_pole_latitude",  NC_FLOAT, 1, &tmp );
        tmp = (float ) real_constants[5];
        ierr = nc_put_att_float( ncid, varID, "grid_north_pole_longitude", NC_FLOAT, 1, &tmp );
     }
     
/*==========================================================================================
   START OF SANITY CHECK
  ==========================================================================================*
     ierr = nc_enddef( ncid );
     for ( n=0; n<num_stored_um_fields; n++ )
         printf( "%hu %hu %hu\n", stored_um_vars[n].stash_code, stored_um_vars[n].x_dim,
                                  stored_um_vars[n].y_dim );
     exit(1);
 *==========================================================================================
   END OF SANITY CHECK
  ==========================================================================================*/

     return;
}

