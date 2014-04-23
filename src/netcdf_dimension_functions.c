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

void set_soil_levels( int ncid, int n, int id );
void set_pressure_levels( int ncid, int n, int id );
void set_altitude( int ncid, int n, int id );
void set_hybrid_levels( int ncid, int n, int id );
void set_surface_levels( int ncid, int n, int id );
void set_sea_surface_levels( int ncid, int n, int id );

/***
 *** ARE_TIME_BND_REQUIRED
 ***
 *** Function that determines if any of the variables in the UM field file are
 *** accummulation data fields.  If so, the time bounds for each accummulation are
 *** read in.
 ***
 ***  INPUT:  ncid         -> ID of the newly created NetCDF file
 ***          cnt          -> # of unique validity times defined in the NetCDF file
 ***          unique_times -> array of IDs for the validity times present in the file 
 ***
 ***    Mark Cheeseman, NIWA
 ***    January 3, 2014
 ***/

void are_time_bnd_required( int ncid, int cnt, int *unique_times ) {

    int    i, j, k, dim_ids[2], varID, flag;
    size_t count[2], offset[2];
    float  *buf;

 /** Do any of the present UM variables contain post-processing that requires **
  ** temporal boundaries?                                                     **/
    flag = 0;
    for ( i=0; i<num_stored_um_fields; i++ ) 
        if ( stored_um_fields[i].lbproc!=0 ) { flag=1; break; }

    if ( flag==1 ) {

       i = nc_def_dim( ncid, "nv", 2, &dim_ids[1] );
      
       buf = (float *) calloc( num_timesteps,sizeof(float) ); 
       count[0] = num_timesteps;
       count[1] = 1;
       offset[0] = 0;

       for ( i=0; i<num_stored_um_fields; i++ ) {
           if ( stored_um_fields[i].lbproc!=0 ) {
              dim_ids[0] = 0;
              for ( j=0; j<cnt; j++ ) 
                  if ( i==unique_times[j] ) { dim_ids[0] = j; }

              j = nc_def_var( ncid,  "time_bnd", NC_FLOAT, 2, dim_ids, &varID );
              j = nc_put_att_text( ncid, varID, "long_name", 37, "start & end times for the cell method" );
              j = nc_put_att_text( ncid, varID, "units", 5, "hours" );

              j = nc_enddef( ncid );

              offset[1] = 0;
              for ( k=0; k<num_timesteps; k++ ) { buf[k] = stored_um_fields[i].time_bnds[k][0]; }
              j = nc_put_vara_float( ncid, varID, offset, count, buf );

              offset[1] = 1;
              for ( k=0; k<num_timesteps; k++ ) { buf[k] = stored_um_fields[i].time_bnds[k][1]; }
              j = nc_put_vara_float( ncid, varID, offset, count, buf );

              j = nc_redef( ncid ); 
           }
       }
       free( buf );
    }

    return;
}

/***
 *** SET_TEMPORAL_DIMENSIONS
 ***
 *** Construct all time-related dimensions for the new NetCDF file.
 ***
 *** INPUT:   ncid -> ID of the newly created NetCDF file
 ***
 ***    Mark Cheeseman, NIWA
 ***    December 28, 2013
 ***/

 void set_temporal_dimensions( int ncid ) {

    int    *p, *dimID, i, j, k, flag, cnt, unique_times[30], dim_id[1];
    char   time_des[31], dim_name[6], calendar[9];
    double tdiff, tdiff2;

  /** Generate a list of the indices of the stored UM variables that possess **/
  /** unique validity date/times.                                            **/
    i = num_stored_um_fields-1;
    unique_times[0] = i;
    cnt = 1;

    for ( j=0; j<i; j++ ) {
        tdiff = difftime( mktime(&stored_um_fields[i].validity),
                          mktime(&stored_um_fields[j].validity) );
        if ( abs(tdiff)>0.01 ) {
           flag = 0; 
           for ( k=0; k<cnt; k++ ) {
               tdiff2 = difftime( mktime(&stored_um_fields[j].validity),
                                  mktime(&stored_um_fields[unique_times[k]].validity) );
               if ( abs(tdiff2)<0.01 ) { flag=1; break; }
           }
           if ( flag==0 ) {
              unique_times[cnt] = j;
              cnt++;
           }
        } 

    }

  /** Set the time dimension ID for each stored UM variable **/
    for ( i=0; i<num_stored_um_fields; i++ ) {
        flag = 0;
        for ( j=0; j<cnt; j++ ) {
            if ( i==unique_times[j] ) { flag=1; break; }  
        }
        if ( flag==0 ) { stored_um_fields[i].t_dim = 0; }
        else           { stored_um_fields[i].t_dim = j; } 
    }

  /** Create dimension & corresponding variable for each unique validity time **/
    dimID = (int *) malloc( cnt*sizeof(int) );

    if ( header[7]==1 ) { sprintf( calendar, "gregorian" ); }
    else                { sprintf( calendar, "360_day" ); }

    for ( i=0; i<cnt; i++ ) {

        sprintf( dim_name, "time%i", i );
        j = nc_def_dim( ncid, dim_name, num_timesteps, &dim_id[0] );

        j = nc_def_var( ncid, dim_name, NC_FLOAT, 1, dim_id, &k );
        j = nc_put_att_text( ncid, k, "calendar", 9, calendar ); 

        if ( stored_um_fields[unique_times[i]].validity.tm_mon<10 ) {
           if ( stored_um_fields[unique_times[i]].validity.tm_mday<10 ) {
              if ( stored_um_fields[unique_times[i]].validity.tm_hour<10 ) {
                 if ( stored_um_fields[unique_times[i]].validity.tm_min<10 ) {
                    sprintf( time_des, "hours since %d-0%d-0%d 0%d:0%d:%d", 
                             stored_um_fields[unique_times[i]].validity.tm_year,
                             stored_um_fields[unique_times[i]].validity.tm_mon,
                             stored_um_fields[unique_times[i]].validity.tm_mday,
                             stored_um_fields[unique_times[i]].validity.tm_hour,
                             stored_um_fields[unique_times[i]].validity.tm_min,
                             stored_um_fields[unique_times[i]].validity.tm_sec );
                 } else {
                    sprintf( time_des, "hours since %d-0%d-0%d 0%d:%d:%d", 
                             stored_um_fields[unique_times[i]].validity.tm_year,
                             stored_um_fields[unique_times[i]].validity.tm_mon,
                             stored_um_fields[unique_times[i]].validity.tm_mday,
                             stored_um_fields[unique_times[i]].validity.tm_hour,
                             stored_um_fields[unique_times[i]].validity.tm_min,
                             stored_um_fields[unique_times[i]].validity.tm_sec );
                 }
              } else {
                 if ( stored_um_fields[unique_times[i]].validity.tm_min<10 ) {
                    sprintf( time_des, "hours since %d-0%d-0%d %d:0%d:%d", 
                             stored_um_fields[unique_times[i]].validity.tm_year,
                             stored_um_fields[unique_times[i]].validity.tm_mon,
                             stored_um_fields[unique_times[i]].validity.tm_mday,
                             stored_um_fields[unique_times[i]].validity.tm_hour,
                             stored_um_fields[unique_times[i]].validity.tm_min,
                             stored_um_fields[unique_times[i]].validity.tm_sec );
                 } else {
                    sprintf( time_des, "hours since %d-0%d-0%d %d:%d:%d", 
                             stored_um_fields[unique_times[i]].validity.tm_year,
                             stored_um_fields[unique_times[i]].validity.tm_mon,
                             stored_um_fields[unique_times[i]].validity.tm_mday,
                             stored_um_fields[unique_times[i]].validity.tm_hour,
                             stored_um_fields[unique_times[i]].validity.tm_min,
                             stored_um_fields[unique_times[i]].validity.tm_sec );
                 }
              }
           } else {
              if ( stored_um_fields[unique_times[i]].validity.tm_hour<10 ) {
                 if ( stored_um_fields[unique_times[i]].validity.tm_min<10 ) {
                    sprintf( time_des, "hours since %d-0%d-%d 0%d:0%d:%d", 
                             stored_um_fields[unique_times[i]].validity.tm_year,
                             stored_um_fields[unique_times[i]].validity.tm_mon,
                             stored_um_fields[unique_times[i]].validity.tm_mday,
                             stored_um_fields[unique_times[i]].validity.tm_hour,
                             stored_um_fields[unique_times[i]].validity.tm_min,
                             stored_um_fields[unique_times[i]].validity.tm_sec );
                 } else {
                    sprintf( time_des, "hours since %d-0%d-%d 0%d:%d:%d", 
                             stored_um_fields[unique_times[i]].validity.tm_year,
                             stored_um_fields[unique_times[i]].validity.tm_mon,
                             stored_um_fields[unique_times[i]].validity.tm_mday,
                             stored_um_fields[unique_times[i]].validity.tm_hour,
                             stored_um_fields[unique_times[i]].validity.tm_min,
                             stored_um_fields[unique_times[i]].validity.tm_sec );
                 }
              } else {
                 if ( stored_um_fields[unique_times[i]].validity.tm_min<10 ) {
                    sprintf( time_des, "hours since %d-0%d-%d %d:0%d:%d", 
                             stored_um_fields[unique_times[i]].validity.tm_year,
                             stored_um_fields[unique_times[i]].validity.tm_mon,
                             stored_um_fields[unique_times[i]].validity.tm_mday,
                             stored_um_fields[unique_times[i]].validity.tm_hour,
                             stored_um_fields[unique_times[i]].validity.tm_min,
                             stored_um_fields[unique_times[i]].validity.tm_sec );
                 } else {
                    sprintf( time_des, "hours since %d-0%d-%d %d:%d:%d", 
                             stored_um_fields[unique_times[i]].validity.tm_year,
                             stored_um_fields[unique_times[i]].validity.tm_mon,
                             stored_um_fields[unique_times[i]].validity.tm_mday,
                             stored_um_fields[unique_times[i]].validity.tm_hour,
                             stored_um_fields[unique_times[i]].validity.tm_min,
                             stored_um_fields[unique_times[i]].validity.tm_sec );
                 }
              }
           }
        } else {
           if ( stored_um_fields[unique_times[i]].validity.tm_mday<10 ) {
              if ( stored_um_fields[unique_times[i]].validity.tm_hour<10 ) {
                 if ( stored_um_fields[unique_times[i]].validity.tm_min<10 ) {
                    sprintf( time_des, "hours since %d-%d-0%d 0%d:0%d:%d", 
                             stored_um_fields[unique_times[i]].validity.tm_year,
                             stored_um_fields[unique_times[i]].validity.tm_mon,
                             stored_um_fields[unique_times[i]].validity.tm_mday,
                             stored_um_fields[unique_times[i]].validity.tm_hour,
                             stored_um_fields[unique_times[i]].validity.tm_min,
                             stored_um_fields[unique_times[i]].validity.tm_sec );
                 } else {
                    sprintf( time_des, "hours since %d-%d-0%d 0%d:%d:%d", 
                             stored_um_fields[unique_times[i]].validity.tm_year,
                             stored_um_fields[unique_times[i]].validity.tm_mon,
                             stored_um_fields[unique_times[i]].validity.tm_mday,
                             stored_um_fields[unique_times[i]].validity.tm_hour,
                             stored_um_fields[unique_times[i]].validity.tm_min,
                             stored_um_fields[unique_times[i]].validity.tm_sec );
                 }
              } else {
                 if ( stored_um_fields[unique_times[i]].validity.tm_min<10 ) {
                    sprintf( time_des, "hours since %d-%d-0%d %d:0%d:%d", 
                             stored_um_fields[unique_times[i]].validity.tm_year,
                             stored_um_fields[unique_times[i]].validity.tm_mon,
                             stored_um_fields[unique_times[i]].validity.tm_mday,
                             stored_um_fields[unique_times[i]].validity.tm_hour,
                             stored_um_fields[unique_times[i]].validity.tm_min,
                             stored_um_fields[unique_times[i]].validity.tm_sec );
                 } else {
                    sprintf( time_des, "hours since %d-%d-0%d %d:%d:%d", 
                             stored_um_fields[unique_times[i]].validity.tm_year,
                             stored_um_fields[unique_times[i]].validity.tm_mon,
                             stored_um_fields[unique_times[i]].validity.tm_mday,
                             stored_um_fields[unique_times[i]].validity.tm_hour,
                             stored_um_fields[unique_times[i]].validity.tm_min,
                             stored_um_fields[unique_times[i]].validity.tm_sec );
                 }
              }
           } else {
              if ( stored_um_fields[unique_times[i]].validity.tm_hour<10 ) {
                 if ( stored_um_fields[unique_times[i]].validity.tm_min<10 ) {
                    sprintf( time_des, "hours since %d-%d-%d 0%d:0%d:%d", 
                             stored_um_fields[unique_times[i]].validity.tm_year,
                             stored_um_fields[unique_times[i]].validity.tm_mon,
                             stored_um_fields[unique_times[i]].validity.tm_mday,
                             stored_um_fields[unique_times[i]].validity.tm_hour,
                             stored_um_fields[unique_times[i]].validity.tm_min,
                             stored_um_fields[unique_times[i]].validity.tm_sec );
                 } else {
                    sprintf( time_des, "hours since %d-%d-%d 0%d:%d:%d", 
                             stored_um_fields[unique_times[i]].validity.tm_year,
                             stored_um_fields[unique_times[i]].validity.tm_mon,
                             stored_um_fields[unique_times[i]].validity.tm_mday,
                             stored_um_fields[unique_times[i]].validity.tm_hour,
                             stored_um_fields[unique_times[i]].validity.tm_min,
                             stored_um_fields[unique_times[i]].validity.tm_sec );
                 }
              } else {
                 if ( stored_um_fields[unique_times[i]].validity.tm_min<10 ) {
                    sprintf( time_des, "hours since %d-%d-%d %d:0%d:%d", 
                             stored_um_fields[unique_times[i]].validity.tm_year,
                             stored_um_fields[unique_times[i]].validity.tm_mon,
                             stored_um_fields[unique_times[i]].validity.tm_mday,
                             stored_um_fields[unique_times[i]].validity.tm_hour,
                             stored_um_fields[unique_times[i]].validity.tm_min,
                             stored_um_fields[unique_times[i]].validity.tm_sec );
                 } else {
                    sprintf( time_des, "hours since %d-%d-%d %d:%d:%d", 
                             stored_um_fields[unique_times[i]].validity.tm_year,
                             stored_um_fields[unique_times[i]].validity.tm_mon,
                             stored_um_fields[unique_times[i]].validity.tm_mday,
                             stored_um_fields[unique_times[i]].validity.tm_hour,
                             stored_um_fields[unique_times[i]].validity.tm_min,
                             stored_um_fields[unique_times[i]].validity.tm_sec );
                 }
              }
           }
        }

        j = nc_put_att_text( ncid, k,    "units", 31, time_des );
        j = nc_put_att_text( ncid, k,     "axis", 1, "T" );
        j = nc_put_att_text( ncid, k, "standard_name", 4, "time" ); 
        j = nc_put_att_text( ncid, k, "long_name", 43, "forecast period (end of reporting period)" ); 

    }
    free( dimID );

  /** Write the appropriate validity time offset values into the NetCDF file. **/
    j = nc_enddef( ncid );
    for ( i=0; i<cnt; i++ ) {  
        sprintf( dim_name, "time%i", i );
        j = nc_inq_varid( ncid, dim_name, &k );
        j = nc_put_var_float( ncid, k, stored_um_fields[unique_times[i]].time_offsets );
    }
    j = nc_redef( ncid );

    p = unique_times;
    are_time_bnd_required( ncid, cnt, p );

    return;
}

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

        free( dimID );

        dimID = (int *)malloc( 3*sizeof(int) );
        dimID[0] = dim_ids[2];
        dimID[1] = dim_ids[0];
        dimID[2] = dim_ids[1];

        ierr = nc_def_var( ncid, "longitude_cell_bnd", vartype, 3, dimID, &varID );
        ierr = nc_put_att_text( ncid, varID,  "long_name", 33, "longitude of cell bounds on earth" );
        ierr = nc_put_att_text(  ncid, varID,     "units", 13, "degrees_north" );

        dimID[0] = dim_ids[3];
        ierr = nc_def_var( ncid, "latitude_cell_bnd", vartype, 3, dimID, &varID );
        ierr = nc_put_att_text( ncid, varID, "long_name", 32, "latitude of cell bounds on earth" );
        ierr = nc_put_att_text( ncid, varID,     "units", 12, "degrees_east" );
        free( dimID );

        dimID = (int *)malloc( 2*sizeof(int) );
        dimID[0] = dim_ids[0];
        dimID[1] = dim_ids[1];

        ierr = nc_def_var( ncid,  "longitude", vartype, 2, dimID, &varID );
        ierr = nc_put_att_text( ncid, varID, "standard_name", 9, "longitude" );
        ierr = nc_put_att_text( ncid, varID,     "long_name",18, "longitude on earth" );
        ierr = nc_put_att_text( ncid, varID,         "units",12, "degrees_east" );
        ierr = nc_put_att_text( ncid, varID,          "axis", 1, "X" );
        ierr = nc_put_att_text( ncid, varID,        "bounds",18, "longitude_cell_bnd" );

        ierr = nc_def_var( ncid,  "latitude", vartype, 2, dimID, &varID );
        ierr = nc_put_att_text(  ncid, varID, "standard_name", 8, "latitude" );
        ierr = nc_put_att_text(  ncid, varID,     "long_name",17, "latitude on earth" );
        ierr = nc_put_att_text(  ncid, varID,         "units",13, "degrees_north" );
        ierr = nc_put_att_text(  ncid, varID,          "axis", 1, "Y" );
        ierr = nc_put_att_text( ncid, varID,         "bounds",17, "latitude_cell_bnd" );
        tmp = 90.0;
        ierr = nc_put_att_float( ncid, varID, "valid_max", NC_FLOAT, 1, &tmp );
        tmp = -90.0;
        ierr = nc_put_att_float( ncid, varID, "valid_min", NC_FLOAT, 1, &tmp );
     }
     free( dimID );    

     return;
}

/***
 *** SET_VERTICAL_DIMENSIONS
 ***
 *** Function that determines how many vertical dimensions are required to
 *** describe the datafields in the UM fields file.  The dimensions are created
 *** in the new NetCDF file.
 ***
 *** INPUT:  ncid    -> id of the newly created NetCDF file
 ***         dim_ids -> array containing the NetCDF dimension IDs
 ***         ndim    -> # of previously defined NetCDF dimensions
 ***         rflag   -> equal to 1 if 32-bit output precision requested by user
 ***
 ***    Mark Cheeseman, NIWA
 ***    December 19, 2013
 ***/

void set_vertical_dimensions( int ncid, int rflag ) {

     int    i, num_instances;

    /*** Check each stored UM variable its 'z-coordinate' ***/
     for ( i=0; i<num_stored_um_fields; i++ ) {

       /** determine # of 'z-levels' for this variable */
         num_instances = stored_um_fields[i].num_slices / num_timesteps;

       /** Create new dimension if variable is 4D **/
         if ( num_instances>1 ) {
            switch ( stored_um_fields[i].lbvc ) {
                 case 1:
                      set_altitude( ncid, num_instances, i );
                      break;
                 case 6:
                      set_soil_levels( ncid, num_instances, i );
                      break;
                 case 8:
                      set_pressure_levels( ncid, num_instances, i );
                      break;
                 case 65:
                      set_hybrid_levels( ncid, num_instances, i );
                      break;
                 case 128:
                      set_sea_surface_levels( ncid, num_instances, i );
                      break;
                 case 129:
                      set_surface_levels( ncid, num_instances, i );
                      break;
            }

         }

     }

    /*** Define dimensions to hold length of ETA arrays ***/
     i = nc_def_dim( ncid, "num_etaT_levels", (size_t ) header[110], 
                    &num_instances );
     i = nc_def_dim( ncid, "num_etaR_levels", (size_t ) (header[110]-1), 
                    &num_instances );

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
