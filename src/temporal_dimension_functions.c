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
 *** SET_TIME_BND
 ***
 *** Function that determines if any of the variables in the UM field file are
 *** accummulation data fields.  If so, the time bounds for each accummulation
 *** are determined.
 ***
 ***  INPUT:  ncid -> ID of the newly created NetCDF file
 ***          dt   -> current timestep 
 ***
 ***    Mark Cheeseman, NIWA
 ***    May 1, 2014
 ***/

void set_time_bnd( int ncid, int var_index ) {

    int   i, ierr, dim_ids[2], varID;
    char  time_bnd_str[9], dim_name[6];
    float dt, tval[2];
    size_t index[2];

 /** Get the dimensions for the new time bounds variable **/

     sprintf( dim_name, "time%i", stored_um_vars[var_index].t_dim );
     ierr = nc_inq_dimid( ncid, dim_name, &dim_ids[0] );

     ierr = nc_inq_dimid( ncid,     "nv", &dim_ids[1] );
     if ( ierr!=NC_NOERR ) { ierr = nc_def_dim( ncid, "nv", 2, &dim_ids[1] ); }

 /** Create a time_bnd variable if an appropriate one is not present **/
 
     sprintf( time_bnd_str, "time_bnd%hu", stored_um_vars[var_index].t_dim );
     ierr = nc_inq_varid( ncid, time_bnd_str, &varID );
     if ( ierr!=NC_NOERR ) {
        ierr = nc_def_var( ncid, time_bnd_str, NC_FLOAT, 2, dim_ids, &varID ); 
        ierr = nc_put_att_text( ncid, varID, "long_name", 37, "start & end times for the cell method" );
        ierr = nc_put_att_text( ncid, varID,     "units",  5, "hours" );
     }
     ierr = nc_enddef( ncid );

 /** Determine the time values for the operation **/
     if ( stored_um_vars[var_index].nt==1 ) {
        tval[0] = 0.0;
        tval[1] = stored_um_vars[var_index].times[0];
        ierr = nc_put_var_float( ncid, varID, tval );
     } else {
        dt = 0.5*(stored_um_vars[var_index].times[1] - stored_um_vars[var_index].times[0]);
        for ( i=0; i<stored_um_vars[var_index].nt; i++ ) {
            index[0] = i;
            index[1] = 0;
            tval[0] = stored_um_vars[var_index].times[i] - dt;
            ierr = nc_put_var1_float( ncid, varID, index, &tval[0] );
            index[1] = 1;
            tval[1] = stored_um_vars[var_index].times[i] + dt;
            ierr = nc_put_var1_float( ncid, varID, index, &tval[1] );
        } 
     }   

     ierr = nc_redef( ncid );
     return;
}

void create_time_dim( int ncid, int var_index, int time_dim_cnt ) {

     int  dimID[1], j, ierr, varID;
     char time_des[40], dim_name[6], calendar[9], mth_str[3], day_str[3],
          hr_str[3], min_str[3], sec_str[3];

   /** Construct an appropriate name for the time dimension **/
     sprintf( dim_name, "time%i", time_dim_cnt );

   /** Define the dimension & a corresponding variable **/
     ierr = nc_def_dim( ncid, dim_name, (size_t ) stored_um_vars[var_index].nt, &dimID[0] );
     ierr = nc_def_var( ncid, dim_name, NC_FLOAT, 1, dimID, &varID );

   /** Add some attributes to the time variable **/

     if ( header[7]==1 ) { sprintf( calendar, "gregorian" ); }
     else                { sprintf( calendar, "360_day" ); }
     ierr = nc_put_att_text( ncid, varID, "calendar", 9, calendar );

//     j = stored_um_vars[var_index].slices[0][0].datatime.tm_mon;
     j = (int ) header[21];
     if ( j<10 ) { snprintf( mth_str, sizeof mth_str, "0%d", j ); }
     else        { snprintf( mth_str, sizeof mth_str,  "%d", j ); }

//     j = stored_um_vars[var_index].slices[0][0].datatime.tm_mday;
     j = (int ) header[22];
     if ( j<10 ) { snprintf( day_str, sizeof day_str, "0%d", j ); }
     else        { snprintf( day_str, sizeof day_str,  "%d", j ); }

//     j = stored_um_vars[var_index].slices[0][0].datatime.tm_hour;
     j = (int ) header[23];
     if ( j<10 ) { snprintf( hr_str, sizeof hr_str, "0%d", j ); }
     else        { snprintf( hr_str, sizeof hr_str,  "%d", j ); }

//     j = stored_um_vars[var_index].slices[0][0].datatime.tm_min;
     j = (int ) header[24];
     if ( j<10 ) { snprintf( min_str, sizeof min_str, "0%d", j ); }
     else        { snprintf( min_str, sizeof min_str,  "%d", j ); }

//     j = stored_um_vars[var_index].slices[0][0].datatime.tm_sec;
     j = (int ) header[25];
     if ( j<10 ) { snprintf( sec_str, sizeof sec_str, "0%d", j ); }
     else        { snprintf( sec_str, sizeof sec_str,  "%d", j ); }

     sprintf( time_des, "hours since %ld-%s-%s %s:%s:%s", header[20],
              mth_str, day_str, hr_str, min_str, sec_str );
//              stored_um_vars[var_index].slices[0][0].datatime.tm_year,
//              mth_str, day_str, hr_str, min_str, sec_str );
     time_des[31] = '\0';
     ierr = nc_put_att_text( ncid, varID, "units", 31, time_des );

     ierr = nc_put_att_text( ncid, varID,          "axis",  1, "T" );
     ierr = nc_put_att_text( ncid, varID, "standard_name",  4, "time" );
     ierr = nc_put_att_text( ncid, varID,     "long_name", 42, "forecast period (end of reporting period)" );

   /** Write time offset values that correspond to the newly created time dimension **/
     ierr = nc_enddef( ncid );
     ierr = nc_put_var_float( ncid, varID, stored_um_vars[var_index].times );
     ierr = nc_redef( ncid );

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

    float *sums, tol, flag;
    int    i, j, num_unique_times; //,k;

 /*
  * Start by determining the # of time dimensions required
  *-----------------------------------------------------------------------*/
    sums = (float *) calloc( num_stored_um_fields,sizeof(float) );
    for ( i=0; i<num_stored_um_fields; i++ ) {
    for ( j=0; j<stored_um_vars[i].nt; j++ ) {
        sums[i] += stored_um_vars[i].times[j]; 
    } 
    } 

    flag = -1.0;
    for ( i=0;   i<num_stored_um_fields; i++ ) {
        if ( sums[i]>-0.5 ) {
           for ( j=i+1; j<num_stored_um_fields; j++ ) {
               if ( sums[j]>-0.5 ) {
                  tol = (sums[i]-sums[j])*(sums[i]-sums[j]);
                  if ( tol<0.0001 ) { sums[j] = flag; }
               }
           }
           flag = flag-1.0;
        }
    } 
    
    num_unique_times = (int )(-1.0*(flag+1.0));

 /*
  * Create a NetCDF dimension for each required time dimension 
  *-----------------------------------------------------------------------*/
    j = 0;
    for ( i=0; i<num_stored_um_fields; i++ ) 
        if ( sums[i]>-0.5 ) { create_time_dim( ncid, i, j ); j++; }

 /*
  * Set t_dim value of each UM field to point to the appropriate NetCDF
  * time dimension for that variable.
  *-----------------------------------------------------------------------*/
    j = 0;
    for ( i=0; i<num_stored_um_fields; i++ ) {
        if ( sums[i]>-0.5 ) { stored_um_vars[i].t_dim = (unsigned short int) j; j++; }
        else                { stored_um_vars[i].t_dim = (unsigned short int) ( -1.0*(sums[i]+1.0) ); }
    }
    free( sums );

 /*
  * If any UM variable is some sort of temporal accummulation, we need to
  * set temporal bounds for that variable. 
  *-----------------------------------------------------------------------*/
    for ( i=0; i<num_stored_um_fields; i++ ) {
        if ( (stored_um_vars[i].lbproc==128)||(stored_um_vars[i].lbproc==4096)||(stored_um_vars[i].lbproc==8192) ) {
           set_time_bnd( ncid, i );
           j++;
        }
    }

/*==============================================================================
  START OF SANITY CHECK
 *==============================================================================
    for ( j=0; j<num_unique_times; j++ ) {
        k= 0;
        printf( "T_DIM = %d STASH_CODES = ", j );
        for ( i=0; i<num_stored_um_fields; i++ ) {
            if ( stored_um_vars[i].t_dim==(unsigned short int) j ) { 
               printf( "%hu ", stored_um_vars[i].stash_code ); 
               k++;
            }
        }
        printf( "NUM_VARS=%d\n", k );
    }
    exit(1);
 *==============================================================================
  END OF SANITY CHECK
 *==============================================================================*/

    return;
}

