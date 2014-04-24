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

/** Function prototypes **/

void set_vertical_dimensions( int ncid, int rflag );
void set_horizontal_dimensions( int ncid, int rflag );
void set_lon_lat_dimensions( int ncid, int iflag, int rflag );
void set_temporal_dimensions( int ncid );
void construct_lat_lon_arrays( int ncid, int rflag, int iflag );
int fill_variables( int ncid, FILE *fid, int iflag, int rflag );

/***
 *** CONSTRUCT_UM_VARIABLES
 ***
 *** Every stored UM variable is defined within the new NetCDF file. 
 ***
 ***  INPUT: ncid  -> ID of the newly created NetCDF file 
 ***         iflag -> equal to 1 if interpolation has been requested by user 
 ***
 ***   Mark Cheeseman, NIWA
 ***   December 29, 2013
 ***/

void construct_um_variables2( int ncid, int iflag ) {

     int     i, ierr, num_z_levels, *dim_ids, ndim, varID, loc, horz_dimids[3];
     size_t *chunksize; 

  /** Determine the IDs of lon & lat dimensions **/
     ierr = nc_inq_dimid( ncid,  "rlon", &horz_dimids[0] );
     ierr = nc_inq_dimid( ncid,  "rlat", &horz_dimids[1] );
     if ( iflag==0 ) { ierr = nc_inq_dimid( ncid, "rlat2", &horz_dimids[2] ); }

     for ( i=0; i<num_stored_um_fields; i++ ) {

  /** Is this a 3D (t,y,x) or 4D (t,z,y,x) variable? **/
         num_z_levels = stored_um_fields[i].num_slices / num_timesteps;
         if ( num_z_levels>1 ) { ndim = 4; }
         else                  { ndim = 3; }

  /** Set the dimensions describing the UM variable. **/
         dim_ids = (int *) malloc( ndim*sizeof(int) );

         dim_ids[0] = stored_um_fields[i].t_dim;
         dim_ids[ndim-1] = horz_dimids[0];
         dim_ids[ndim-2] = horz_dimids[1]; 

         if ( iflag==0 ) {  
            if ( stored_um_fields[i].ny==int_constants[6] ) { dim_ids[ndim-2] = horz_dimids[1]; }
            else                                            { dim_ids[ndim-2] = horz_dimids[2]; }
         }

         if ( ndim==4 ) { dim_ids[1] = stored_um_fields[i].z_dim; }

  /** Define the appropriate NetCDF variable **/
         ierr = nc_def_var( ncid, stored_um_fields[i].name, stored_um_fields[i].vartype, ndim, dim_ids, &varID );
         free( dim_ids );

 /** Set the chunking attribute for this variable **/
         chunksize = (size_t* ) malloc( ndim*sizeof(size_t) );
         chunksize[0] = 1;
         chunksize[ndim-1] = stored_um_fields[i].nx;
         if ( iflag==0 ) { chunksize[ndim-2] = stored_um_fields[i].ny; }
         else            { chunksize[ndim-2] = int_constants[6]; }
         if ( ndim==4 ) { chunksize[1] = 1; }

         ierr = nc_def_var_chunking( ncid, varID, NC_CHUNKED, chunksize );
         free( chunksize );

 /** Set the data compression attribute for this variable **/
         ierr = nc_def_var_deflate( ncid, varID, NC_NOSHUFFLE, 1, 4 );

 /*** Output details about the coordinate system used to describe field ***/
         if ( stored_um_fields[i].coordinates==101 ) {
            ierr = nc_put_att_text( ncid, varID, "grid_mapping", 12, "rotated_pole" );
         } 
         

 /*** Determine if any post-processing was performed on the data-field.  If so, ***/
 /*** indicate the operation performed.***/
        if ( stored_um_fields[i].lbproc>0 ) { 
           if ( (stored_um_fields[i].lbproc==8)||(stored_um_fields[i].lbproc==64) ) {  
              ierr = nc_put_att_text( ncid, varID, "cell_method", 11, "space:zonal" );
           }
           else if ( stored_um_fields[i].lbproc==128 )  { 
                if ( stored_um_fields[i].accum==0 ) { ierr = nc_put_att_text( ncid, varID, "cell_method", 9, "time:mean" ); }
                else                                { ierr = nc_put_att_text( ncid, varID, "cell_method", 8, "time:sum" ); }
           }
           else if ( stored_um_fields[i].lbproc==4096 ) { ierr = nc_put_att_text( ncid, varID, "cell_method", 8, "time:min" ); }
           else if ( stored_um_fields[i].lbproc==8192 ) { ierr = nc_put_att_text( ncid, varID, "cell_method", 8, "time:max" ); }
        }

/** Add the appropriate attributes **/
         loc = stored_um_fields[i].xml_index; 
         ierr = nc_put_att_int( ncid, varID,   "stash_model", NC_INT, 1, &
                                    um_vars[loc].model );
         ierr = nc_put_att_int( ncid, varID, "stash_section", NC_INT, 1, &
                                    um_vars[loc].section );
         ierr = nc_put_att_int( ncid, varID,    "stash_item", NC_INT, 1, &
                                    um_vars[loc].code );
         ierr = nc_put_att_float( ncid, varID, "valid_max", NC_FLOAT, 1, &
                                      um_vars[loc].validmax );
         ierr = nc_put_att_float( ncid, varID, "valid_min", NC_FLOAT, 1, &
                                      um_vars[loc].validmin );
         ierr = nc_put_att_text( ncid, varID, "long_name", 100, um_vars[loc].longname );
         ierr = nc_put_att_text( ncid, varID, "standard_name", 75, um_vars[loc].stdname );
         ierr = nc_put_att_text( ncid, varID, "units", 25, um_vars[loc].units );
         
     }

   /*** Create the ETA arrays -used to store the coefficients needed to determine ***
    *** model depth on hybrid levels.                                             ***/

     ndim = 1;
     dim_ids = (int *) malloc( ndim*sizeof(int) );

     ierr = nc_inq_dimid( ncid, "num_etaT_levels", &dim_ids[0] );
     ierr = nc_def_var( ncid, "eta_theta", NC_FLOAT, ndim, dim_ids, &varID );

     ierr = nc_inq_dimid( ncid, "num_etaR_levels", &dim_ids[0] );
     ierr = nc_def_var( ncid, "eta_rho",   NC_FLOAT, ndim, dim_ids, &varID );

     free( dim_ids );
     return;
}


/***
 *** CREATE_NETCDF_FILE 
 ***
 *** Subroutine that creates a NetCDF file with all the appropriate variables
 *** defined. 
 ***
 ***   Mark Cheeseman, NIWA
 ***   November 29, 2013
 ***/

int create_netcdf_file( char *um_file, int iflag, int rflag ) {
     
     int   ncid, ierr, pos;
     char  forecast_ref_time[18], netcdf_filename[50], *str, *dest;
     char  *nzlam = "nzlam";
     char  *nzcsm = "nzcsm";
     char  mth_str[2], day_str[2], min_str[2], hr_str[2], sec_str[2];
     FILE *fid;


 /**=========================================================================**
  ** STEP 0:  FILE CREATION                                                  **
  **=========================================================================**/ 

     fid = fopen( um_file, "r" );
     if ( fid==NULL ) { return 999; }

 /* 
  * 0a) Create an appropriate name for the NetCDF file 
  *---------------------------------------------------------------------------*/
     dest = strstr( um_file, ".um" );
     if ( dest!=NULL ) { 
        pos = dest - um_file;

        str = malloc( 1 + strlen(um_file) );
        if ( str ) { strncpy( str, um_file, pos ); }
        else       { return 999; }
        str[pos] = '\0';

        snprintf( netcdf_filename, sizeof netcdf_filename, "%s.nc", str ); 
        free( str );
     } else {
        snprintf( netcdf_filename, sizeof netcdf_filename, "%s.nc", um_file );
     }

     ierr = nc_set_chunk_cache( 129600000, 101, 0.75 );
     ierr = nc_create( netcdf_filename, NC_NETCDF4, &ncid ); 
     if ( ierr != NC_NOERR ) { return 999; }
 
 /**=========================================================================**
  ** STEP 1:  DIMENSIONS                                                     **
  **=========================================================================**/ 

 /* 
  * 1a) Temporal Dimensions
  *---------------------------------------------------------------------------*/
     set_temporal_dimensions( ncid );

 /* 
  * 1b) Horizontal (Lat/Lon) Dimensions
  *---------------------------------------------------------------------------*/
     set_lon_lat_dimensions( ncid, iflag, rflag );
     set_horizontal_dimensions( ncid, rflag );

 /* 
  * 1c) Vertical Dimensions
  *---------------------------------------------------------------------------*/
     set_vertical_dimensions( ncid, rflag );

 /**=========================================================================**
  ** STEP 2:  VARIABLES                                                      **
  **=========================================================================**/ 

     construct_um_variables2( ncid, iflag );

 /**=========================================================================**
  ** STEP 3:  GLOBAL ATTRIBUTES                                              **
  **=========================================================================**/ 

 /*
  * Construct a properly formatted forecast reference string 
  *--------------------------------------------------------------------------*/
     if ( header[28]>10 ) { sprintf( mth_str, "%ld", header[28] ); }
     else                 { sprintf( mth_str, "0%ld", header[28] ); }
     if ( header[29]>10 ) { sprintf( day_str, "%ld", header[29] ); }
     else                 { sprintf( day_str, "0%ld", header[29] ); }
     if ( header[30]>10 ) { sprintf(  hr_str, "%ld", header[30] ); }
     else                 { sprintf(  hr_str, "0%ld", header[30] ); }
     if ( header[31]>10 ) { sprintf( min_str, "%ld", header[31] ); }
     else                 { sprintf( min_str, "0%ld", header[31] ); }
     if ( header[32]>10 ) { sprintf( sec_str, "%ld", header[32] ); }
     else                 { sprintf( sec_str, "0%ld", header[32] ); }
     sprintf( forecast_ref_time, "%ld-%s-%s %s:%s:%s", header[27], mth_str, day_str, hr_str,
                                                       min_str, sec_str );
     ierr = nc_put_att_text( ncid, NC_GLOBAL, "forecast_reference_time", strlen(forecast_ref_time), 
                             forecast_ref_time ); 

 /*
  * Add other global attributes to output NetCDF file 
  *--------------------------------------------------------------------------*/
     ierr = nc_put_att_double( ncid, NC_GLOBAL, "grid_north_pole_latitude", NC_DOUBLE, 1, &real_constants[4] );
     ierr = nc_put_att_double( ncid, NC_GLOBAL, "grid_north_pole_longitude", NC_DOUBLE, 1, &real_constants[5] ); 
     ierr = nc_put_att_text( ncid, NC_GLOBAL, "history", 39, "UM fields file reformatted by um2netcdf" ); 
     ierr = nc_put_att_text( ncid, NC_GLOBAL, "input_uri", strlen(um_file), um_file ); 
     ierr = nc_put_att_text( ncid, NC_GLOBAL, "institution", 4, "NIWA" ); 
     ierr = nc_put_att_text( ncid, NC_GLOBAL, "conventions", 6, "CF-1.5" ); 
     ierr = nc_put_att_long( ncid, NC_GLOBAL, "um_version_number", NC_LONG, 1, &header[11] );
     ierr = nc_put_att_int( ncid, NC_GLOBAL, "met_office_ps", NC_INT, 1, 0 ); 
     ierr = nc_put_att_int( ncid, NC_GLOBAL,      "niwa_eps", NC_INT, 1, 0 ); 
     ierr = nc_put_att_int( ncid, NC_GLOBAL,       "rose_id", NC_INT, 1, 0 ); 

     if ( strstr( um_file, nzlam )!=NULL ) { 
        ierr = nc_put_att_text( ncid, NC_GLOBAL, "model_name", 5, "nzlam" ); 
        ierr = nc_put_att_text( ncid, NC_GLOBAL, "data_assimilation_method", 6, "pseudo" ); 
     }else {
        if ( strstr( um_file, nzcsm )!=NULL ) { 
           ierr = nc_put_att_text( ncid, NC_GLOBAL, "model_name", 5, "nzcsm" ); 
        } else { ierr = nc_put_att_text( ncid, NC_GLOBAL, "model_name", 7, "unknown" ); }
     }          
     if ( header[3]>100 ) { ierr = nc_put_att_text( ncid, NC_GLOBAL, "grid_mapping_name", 26, "rotated_latitude_longitude" ); } 
     if ( header[11]<=804 ) { ierr = nc_put_att_text( ncid, NC_GLOBAL, "dynamical_core", 12, "new_dynamics" ); }
     else                   { ierr = nc_put_att_text( ncid, NC_GLOBAL, "dynamical_core",  8, "end_game" ); }

 /*
  * Close the input UM fields file. 
  *--------------------------------------------------------------------------*/
     ierr = nc_enddef( ncid );
     fclose( fid );

     return ncid;
}


/***
 *** FILL_NETCDF_FILE
 ***
 *** Subroutine that transfers the raw UM field data into the newly created 
 *** NetCDF file under the appropriate variable.
 ***
 ***   Mark Cheeseman, NIWA
 ***   December 4, 2013
 ***/

int fill_netcdf_file( int ncid, char *filename, int iflag, int rflag ) {

    int   i; 
    FILE *fid;

 /*
  * Re-open the UM fields file. 
  --------------------------------------------------------------------------*/
     fid = fopen( filename, "r" );

 /*
  * Output the lon/lat data values 
  *-------------------------------------------------------------------------*/
     construct_lat_lon_arrays( ncid, rflag, iflag );

     printf( "-----------------------------------------------------------\n" );
     printf( "  Stash_Code        Variable_Name             Dimensions\n" );
     printf( "-----------------------------------------------------------\n" );

     i = fill_variables( ncid, fid, iflag, rflag );
     if ( i==-1 ) { 
        printf( "ERROR: write failed\n" );
        return i;
     }

 /*** Finish by closing the UM fields and NetCDF files ***/
     fclose( fid );
     free( um_vars );

     i = nc_close( ncid );

     return 1;
}
