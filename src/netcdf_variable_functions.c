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

#include <netcdf.h>
#include <stdlib.h>
#include <stdio.h>
#include "field_def.h"
#include <string.h>
#include <math.h>
#include <stdint.h>

/** Function prototypes **/

double *interp_do_nothing( double *val, int nx, int ny, double maxval, double minval );
double *u_to_p_point_interp_c_grid( double *val, int nx, int ny, double maxval, double minval );
double *v_to_p_point_interp_c_grid( double *val, int nx, int ny, double maxval, double minval );
double *b_to_c_grid_interp_u_points( double *val, int nx, int ny, double maxval, double minval );
void   endian_swap_4bytes( void *ptr, int nchunk );
void   wgdos_unpack( FILE *fh, unsigned short nx, unsigned short ny, double *buf,
                     double mdi );
void   ieee_usage_message();



void write_interpolated_fields( int ncid, FILE *fid, int rflag ) {

     int     n, i, j, k, ndim, cnt, loc, varid;
     size_t *count, *offset;
     double *buf, maxval, minval;
     float  *fbuf;
     char    name[45];

     for ( n=0; n<num_stored_um_fields; n++ ) {

         strcpy( name, stored_um_vars[n].name );
         i = nc_inq_varid( ncid, name, &varid );

       /** Determine # of dimensions for current UM variable **/
         ndim = 3;
         if ( stored_um_vars[n].nz>1 ) { ndim = 4; }

       /*** Allocate & set the sizes of a 2D slice in the UM variable ***/
       /*** Remember that COUNT = COUNT[NT,NZ,NY,NX]                  ***/

         count = (size_t *) malloc( ndim*sizeof(size_t) );
         count[ndim-2] = int_constants[6];
         count[ndim-1] = stored_um_vars[n].nx;
         count[0]      = 1;   // only 1 timeslice printed at a time

       /*** For 4D variables, we only want to output 1 depth level at a time ***/

         if ( ndim==4 ) { count[1] = 1; }

         offset = (size_t *) calloc( ndim,sizeof(size_t) );

       /** Initialize function pointer to proper interpolation function **/
         switch ( stored_um_vars[n].grid_type ) {
                 case 11: 
                        field_interpolation = &b_to_c_grid_interp_u_points;
                 case 18: 
                        field_interpolation = &u_to_p_point_interp_c_grid;
                 case 19: 
                        field_interpolation = &v_to_p_point_interp_c_grid;
                 default:
                        field_interpolation = &interp_do_nothing;
         }

       /** Create a buffer to hold the raw data values **/
         cnt = stored_um_vars[n].nx*stored_um_vars[n].ny;
         buf = (double *) calloc( cnt,sizeof(double) );

         loc = stored_um_vars[n].xml_index;
         minval = um_vars[loc].validmin;
         maxval = um_vars[loc].validmax;

       /*** Write the 2D data slices belonging to the UM field one at a time ***/
         if ( ndim==3 ) {

            offset[0] = 0;
            for ( k=0; k<stored_um_vars[n].nt; k++ ) {

                fseek( fid, stored_um_vars[n].slices[k][0].location*wordsize, SEEK_SET );
                fread( buf, wordsize, cnt, fid );
                endian_swap( buf, cnt );

                buf = field_interpolation( buf, stored_um_vars[n].nx, stored_um_vars[n].ny,
                                           maxval, minval );

                if ( rflag==0 ) {
                   i = nc_put_vara_double( ncid, varid, offset, count, buf );
                 } else {
                   fbuf = (float *) malloc( cnt*sizeof(float) );
                   for ( i=0; i<cnt; i++ ) { fbuf[i] = (float ) buf[i]; }
                   i = nc_put_vara_float( ncid, varid, offset, count, fbuf );
                   free( fbuf );
                }
                offset[0]++;
             }

         } else if ( ndim==4 ) {
 
            offset[0] = 0;
            for ( k=0; k<stored_um_vars[n].nt; k++ ) {
                offset[1] = 0;
                for ( j=0; j<stored_um_vars[n].nz; j++ ) {
                    fseek( fid, stored_um_vars[n].slices[k][j].location*wordsize, SEEK_SET );
                    fread( buf, wordsize, cnt, fid );
                    endian_swap( buf, cnt );

                    buf = field_interpolation( buf, stored_um_vars[n].nx, stored_um_vars[n].ny,
                                               maxval, minval );

                    if ( rflag==0 ) {
                       i = nc_put_vara_double( ncid, varid, offset, count, buf );
                    } else {
                       fbuf = (float *) malloc( cnt*sizeof(float) );
                       for ( i=0; i<cnt; i++ ) { fbuf[i] = (float ) buf[i]; }
                       i = nc_put_vara_float( ncid, varid, offset, count, fbuf );
                       free( fbuf );
                    }
                    offset[1]++;
                }
                offset[0]++;
             }

         }

         free( count );
         free( offset );
         free( buf );

     }  // End of FOR LOOP

     return;
}


void write_uninterpolated_fields( int ncid, FILE *fid, int rflag ) {

     int     n, i, j, k, ndim, cnt, loc, varid;
     size_t *count, *offset;
     double *buf, maxval, minval;
     float  *fbuf;
     char    name[45];

     for ( n=0; n<num_stored_um_fields; n++ ) {

         strcpy( name, stored_um_vars[n].name );
         i = nc_inq_varid( ncid, name, &varid );

       /** Determine # of dimensions for current UM variable **/
         ndim = 3;
         if ( stored_um_vars[n].nz>1 ) { ndim = 4; }

       /*** Allocate & set the sizes of a 2D slice in the UM variable ***/
       /*** Remember that COUNT = COUNT[NT,NZ,NY,NX]                  ***/

         count = (size_t *) malloc( ndim*sizeof(size_t) );
         count[ndim-2] = int_constants[6];
         count[ndim-1] = stored_um_vars[n].nx;
         count[0]      = 1;   // only 1 timeslice printed at a time

       /*** For 4D variables, we only want to output 1 depth level at a time ***/

         if ( ndim==4 ) { count[1] = 1; }

         offset = (size_t *) calloc( ndim,sizeof(size_t) );

       /** Create a buffer to hold the raw data values **/
         cnt = stored_um_vars[n].nx*stored_um_vars[n].ny;
         buf = (double *) calloc( cnt,sizeof(double) );

         loc = stored_um_vars[n].xml_index;
         minval = um_vars[loc].validmin;
         maxval = um_vars[loc].validmax;

       /*** Write the 2D data slices belonging to the UM field one at a time ***/
         if ( ndim==3 ) {

            offset[0] = 0;
            for ( k=0; k<stored_um_vars[n].nt; k++ ) {

                fseek( fid, stored_um_vars[n].slices[k][0].location*wordsize, SEEK_SET );
                fread( buf, wordsize, cnt, fid );
                endian_swap( buf, cnt );

                if ( rflag==0 ) {
                   i = nc_put_vara_double( ncid, varid, offset, count, buf );
                 } else {
                   fbuf = (float *) malloc( cnt*sizeof(float) );
                   for ( i=0; i<cnt; i++ ) { fbuf[i] = (float ) buf[i]; }
                   i = nc_put_vara_float( ncid, varid, offset, count, fbuf );
                   free( fbuf );
                }
                offset[0]++;
             }

         } else if ( ndim==4 ) {
 
            offset[0] = 0;
            for ( k=0; k<stored_um_vars[n].nt; k++ ) {
                offset[1] = 0;
                for ( j=0; j<stored_um_vars[n].nz; j++ ) {
                    fseek( fid, stored_um_vars[n].slices[k][j].location*wordsize, SEEK_SET );
                    fread( buf, wordsize, cnt, fid );
                    endian_swap( buf, cnt );

                    if ( rflag==0 ) {
                       i = nc_put_vara_double( ncid, varid, offset, count, buf );
                    } else {
                       fbuf = (float *) malloc( cnt*sizeof(float) );
                       for ( i=0; i<cnt; i++ ) { fbuf[i] = (float ) buf[i]; }
                       i = nc_put_vara_float( ncid, varid, offset, count, fbuf );
                       free( fbuf );
                    }
                    offset[1]++;
                }
                offset[0]++;
             }

         }

         free( count );
         free( offset );
         free( buf );

     }  // End of FOR LOOP

     return;
}

/***
 *** OUTPUT_UM_FIELDS 
 ***
 *** The input UM fields file is parsed for 2D dataslices belonging to
 *** each previously identified UM variable.  Appropriate dataslices
 *** are read, processed (if required) and written into the new NetCDF
 *** file.  
 ***
 ***  INPUT:  ncid -> ID of the newly created NetCDF file 
 ***           fid -> file pointer to the UM fields file
 ***         iflag -> denotes if interpolation is to be performed on the 
 ***                  input data fields (1=YES, 0=NO)
 ***         rflag -> denotes whether 32 or 64-bit output is desired
 ***                  (0->64 bit, 1->32 bit)
 ***          
 ***   Mark Cheeseman, NIWA
 ***   December 23, 2013
 ***/

int output_um_fields( int ncid, FILE *fid, int iflag, int rflag ) {

     int    ierr, n, varid;
     char   name[45];
     size_t dimlen;
     float  *buf;

  /*
   * Perform a check for WGDOS-packed fields (currently not supported)
   *-------------------------------------------------------------------*/
     for ( n=0; n<num_stored_um_fields; n++ ) {
         if ( stored_um_vars[n].slices[0][0].lbpack==1 ) {
            ieee_usage_message();
            return -1;
         }   
     }

  /*
   * List the UM variables to be written 
   *-------------------------------------------------------------------*/
     printf( "-----------------------------------------------------------\n" );
     printf( "  Stash_Code        Variable_Name             Dimensions\n" );
     printf( "-----------------------------------------------------------\n" );

     for ( n=0; n<num_stored_um_fields; n++ ) {
         strcpy( name, stored_um_vars[n].name ); 
         printf( "     %5d %25s     [%4d x %4d", stored_um_vars[n].stash_code, 
                                                 stored_um_vars[n].name, 
                                                 stored_um_vars[n].nx, 
                                                 stored_um_vars[n].ny );
         if ( stored_um_vars[n].nz>1 ) {  printf( " x %3d", stored_um_vars[n].nz ); }
         printf( "]\n" );
     }

  /*
   * Write the UM field to hard disk one 2D data slice at a time. 
   *-------------------------------------------------------------------*/
     if ( iflag==1 ) { write_interpolated_fields( ncid, fid, rflag ); }
     else            { write_uninterpolated_fields( ncid, fid, rflag ); }

    /*** Output the coefficients for the ETA arrays ***/

     dimlen = (size_t ) header[110];

     buf = (float *) malloc( dimlen*sizeof(float) );
     for ( n=0; n<dimlen; n++ ) { buf[n] = level_constants[0][n]; } 
     ierr = nc_inq_varid( ncid, "eta_theta", &varid );   
     ierr = nc_put_var_float( ncid, varid, buf );
     free( buf ); 

     buf = (float *) malloc( (dimlen-1)*sizeof(float) );
     for ( n=0; n<dimlen-1; n++ ) { buf[n] = level_constants[1][n]; } 
     ierr = nc_inq_varid( ncid, "eta_rho", &varid );   
     ierr = nc_put_var_float( ncid, varid, buf );

     return 1;
}

