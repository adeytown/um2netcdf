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
 *** ARE_TIME_BND_REQUIRED
 ***
 *** Function that determines if any of the variables in the UM field file are
 *** accummulation data fields.  If so, the time bounds for each accummulation
 *** are determined.
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

 /** Do any of the UM variables present contain post-processing that requires **
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

void create_time_dim( int ncid, int num_times, int *vars ) {

     int  *dimID, dim_id[2], i, j ,k;
     char time_des[40], dim_name[6], calendar[9], mth_str[3], day_str[3],
          hr_str[3], min_str[3], sec_str[3];


     dimID = (int *) malloc( num_times*sizeof(int) );

     if ( header[7]==1 ) { sprintf( calendar, "gregorian" ); }
     else                { sprintf( calendar, "360_day" ); }

     for ( i=0; i<num_times; i++ ) {

     /** Create the appropriate time dimension **/
        sprintf( dim_name, "time%i", i );
        j = nc_def_dim( ncid, dim_name, num_timesteps, &dim_id[0] );

     /** Construct timestamp for variable's validity time **/
        j = stored_um_fields[vars[i]].validity.tm_mon;
        if ( j<10 ) { snprintf( mth_str, sizeof mth_str, "0%d", j ); }
        else        { snprintf( mth_str, sizeof mth_str,  "%d", j ); }

        j = stored_um_fields[vars[i]].validity.tm_mday;
        if ( j<10 ) { snprintf( day_str, sizeof day_str, "0%d", j ); }
        else        { snprintf( day_str, sizeof day_str,  "%d", j ); }

        j = stored_um_fields[vars[i]].validity.tm_hour;
        if ( j<10 ) { snprintf( hr_str, sizeof hr_str, "0%d", j ); }
        else        { snprintf( hr_str, sizeof hr_str,  "%d", j ); }

        j = stored_um_fields[vars[i]].validity.tm_min;
        if ( j<10 ) { snprintf( min_str, sizeof min_str, "0%d", j ); }
        else        { snprintf( min_str, sizeof min_str,  "%d", j ); }

        j = stored_um_fields[vars[i]].validity.tm_sec;
        if ( j<10 ) { snprintf( sec_str, sizeof sec_str, "0%d", j ); }
        else        { snprintf( sec_str, sizeof sec_str,  "%d", j ); }

        sprintf( time_des, "hours since %d-%s-%s %s:%s:%s",
                 stored_um_fields[vars[i]].validity.tm_year,
                 mth_str, day_str, hr_str, min_str, sec_str );
        time_des[31] = '\0';

     /** Create variable to hold time offset values for time dimension **/
        j = nc_def_var( ncid, dim_name, NC_FLOAT, 1, dim_id, &k );
        j = nc_put_att_text( ncid, k,      "calendar",  9, calendar );
        j = nc_put_att_text( ncid, k,         "units", 31, time_des );
        j = nc_put_att_text( ncid, k,          "axis",  1, "T" );
        j = nc_put_att_text( ncid, k, "standard_name",  4, "time" );
        j = nc_put_att_text( ncid, k,     "long_name", 43, "forecast period (end of reporting period)" );

     /** Write time offset values that correspond to the newly created time dimension **/
        j = nc_enddef( ncid );
        j = nc_put_var_float( ncid, k, stored_um_fields[vars[i]].time_offsets );
        j = nc_redef( ncid );

    }
    free( dimID );
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

//    int    *p, *dimID, i, j, k, flag, cnt, unique_times[30], dim_id[1];
    int    i, cnt;
    int    num_validity_times, *field_ids;
//    char   time_des[40], dim_name[6], calendar[9], mth_str[3], day_str[3],
//           hr_str[3], min_str[3], sec_str[3];
    double tdiff; //, tdiff2;


  /*** How many unique validity exist in the input UM file? ***/
    num_validity_times = 1;
    cnt = 1;
    for ( i=1; i<num_stored_um_fields; i++ ) {
        tdiff = difftime( mktime(&stored_um_fields[0].validity),
                          mktime(&stored_um_fields[i].validity) );
        if ( tdiff>0.00001 ) { num_validity_times++; }
    }

    field_ids = (int *) calloc( num_validity_times,sizeof(int) );
    if ( num_validity_times==1 ) {

       for ( i=0; i<num_stored_um_fields; i++ ) { stored_um_fields[i].t_dim = 0; }
       create_time_dim( ncid, num_validity_times, field_ids );

    } else {

       for ( i=1; i<num_stored_um_fields; i++ ) {
           tdiff = difftime( mktime(&stored_um_fields[0].validity),
                          mktime(&stored_um_fields[i].validity) );
           if ( tdiff>0.00001 ) { field_ids[cnt] = i; cnt++; }
       }
       create_time_dim( ncid, num_validity_times, field_ids );

    }

    are_time_bnd_required( ncid, cnt, field_ids );
    free( field_ids );

    return;
}

