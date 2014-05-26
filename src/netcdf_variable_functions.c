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

double *interp_do_nothing( double *val, int index );
double *u_to_p_point_interp_c_grid( double *val, int index );
double *v_to_p_point_interp_c_grid( double *val, int index );
double *b_to_c_grid_interp_u_points( double *val, int index );
float  *interp_do_nothing_float( double *val, int index );
float  *u_to_p_point_interp_c_grid_float( double *val, int index );
float  *v_to_p_point_interp_c_grid_float( double *val, int index );
float  *b_to_c_grid_interp_u_points_float( double *val, int index );
void endian_swap_4bytes( void *ptr, int nchunk );
void wgdos_unpack( FILE *fh, unsigned short nx, unsigned short ny, double *buf,
                   double mdi );
void ieee_usage_message();


/***
 *** WRITE_INTERPOLATED_FIELDS 
 ***
 *** Subroutine that reads in each UM variable by 1 2D data slice at a time.
 *** Each slice is interpolated until the P-grid (if required) and the resulting
 *** 2D slice is written into the appropriate NetCDF variable. 
 ***
 ***  INPUT:  ncid -> ID of the newly created NetCDF file 
 ***           fid -> file pointer to the UM fields file
 ***         rflag -> denotes whether 32 or 64-bit output is desired
 ***                  (0->64 bit, 1->32 bit)
 ***          
 ***   Mark Cheeseman, NIWA
 ***   May 19, 2014
 ***/

void write_interpolated_fields( int ncid, FILE *fid, int rflag ) {

     int     n, i, j=0, k, kk, ndim, cnt, varid;
     size_t *count, *offset;
     double *buf;
     float  *fbuf=NULL;
     char    name[45];

     for ( n=0; n<num_stored_um_fields; n++ ) {

         strcpy( name, stored_um_vars[n].name );
         i = nc_inq_varid( ncid, name, &varid );

       /** Determine # of dimensions for current UM variable **/
         ndim = 3;
         if ( stored_um_vars[n].nz>1 ) { ndim = 4; }

       /*** Allocate & set the sizes of a 2D slice in the UM variable ***/
       /*** Remember that COUNT = COUNT[NT,NZ,NY,NX]                  ***/

         offset = (size_t *) calloc( ndim,sizeof(size_t) );
         count = (size_t *) malloc( ndim*sizeof(size_t) );
         count[ndim-2] = int_constants[6];
         count[ndim-1] = stored_um_vars[n].nx;
         count[0]      = 1;   // only 1 timeslice printed at a time

       /*** Count the number of elements to be read for a single 2D data slice ***/

         cnt = stored_um_vars[n].nx*stored_um_vars[n].ny;

       /*** Initialize function pointer to proper interpolation function ***/

         switch ( stored_um_vars[n].grid_type ) {
                 case 11:
                        field_interpolation = &b_to_c_grid_interp_u_points; 
                        break;
                 case 18: 
                        field_interpolation = &u_to_p_point_interp_c_grid; 
                        break;
                 case 19: 
                        field_interpolation = &v_to_p_point_interp_c_grid; 
                        break;
                 default:
                        field_interpolation = &interp_do_nothing; 
                        if ( int_constants[6]!=stored_um_vars[n].ny ) {
                           printf( "ERROR: NY dimension of %s is not equal to %ld\n", name, int_constants[6] );
                           printf( "       Check the umgrid value in the XML stashfile for this field\n\n" );
                           exit(1);
                        }
                        break;
         }

       /*** For a 3D UM variable [NX,NY,NT], read in a single 2D data slice per data time   ***/
       /*** valid for this UM variable.  Then apply the appropriate interpolation and write ***/
       /*** the final field to hard disk.                                                   ***/

         if ( ndim==3 ) {

               offset[0] = 0;
               for ( k=0; k<stored_um_vars[n].nt; k++ ) {

               /* Read in a 2D data slice. Apply appropriate endian swap on the data */
                   buf = (double *) malloc( cnt*sizeof(double) );
                   fseek( fid, stored_um_vars[n].slices[k][0].location*wordsize, SEEK_SET );
                   fread( buf, wordsize, cnt, fid );
                   endian_swap( buf, cnt );

               /* Apply appropriate interpolation on values */
                   buf = field_interpolation( buf, n );

               /* Write interpolated 2D slice to hard disk */
                   if ( rflag==0 ) { 
                      i = nc_put_vara_double( ncid, varid, offset, count, buf ); 
                      free( buf );
                   } else {
                      i = int_constants[6]*((int ) stored_um_vars[n].nx);
                      fbuf = (float *) malloc( i*sizeof(float) );
                      for ( j=0; j<i; j++ ) { fbuf[j] = (float ) buf[j]; }
                      i = nc_put_vara_float( ncid, varid, offset, count, fbuf );
                      free( fbuf );
                      free( buf );
                   }

               /* Update the time counter for this UM variable */
                   offset[0]++;
               }

         } 

       /*** For a 4D UM variable [NX,NY,NZ,NT], read in a single 2D data slice per level and ***/
       /*** data time valid for this UM variable.  Then apply the appropriate interpolation  ***/
       /*** and write the final field to hard disk.                                          ***/

         else if ( ndim==4 ) {
              count[1] = 1; 
              for ( k=0; k<stored_um_vars[n].nt; k++ ) {
                  offset[1] = 0;
                  for ( j=0; j<stored_um_vars[n].nz; j++ ) {

               /* Read in a 2D data slice. Apply appropriate endian swap on the data */
                      buf = (double *) malloc( cnt*sizeof(double) );
                      fseek( fid, stored_um_vars[n].slices[k][j].location*wordsize, SEEK_SET );
                      fread( buf, wordsize, cnt, fid );
                      endian_swap( buf, cnt );

               /* Apply appropriate interpolation on values */
                      buf = field_interpolation( buf, n ); 

               /* Write interpolated 2D slice to hard disk */
                      if ( rflag==0 ) { 
                         i = nc_put_vara_double( ncid, varid, offset, count, buf ); 
                         free( buf );
                      } else {
                         i = int_constants[6]*((int ) stored_um_vars[n].nx);
                         fbuf = (float *) malloc( i*sizeof(float) );
                         for ( kk=0; kk<i; kk++ ) { fbuf[kk] = (float ) buf[kk]; }
                         i = nc_put_vara_float( ncid, varid, offset, count, fbuf );
                         free( buf );
                         free( fbuf );
                      }
               /* Update the level counter for this UM variable */
                      offset[1]++;
                  }
               /* Update the time counter for this UM variable */
                  offset[0]++;
              }
         } // End of NDIM if block 

         free( count );
         free( offset );

     }  // End of FOR LOOP

     return;
}


/***
 *** WRITE_UNINTERPOLATED_FIELDS 
 ***
 *** Subroutine that reads in each UM variable by 1 2D data slice at a time.
 *** The 2D data slice is then written into the appropriate NetCDF variable. 
 ***
 ***  INPUT:  ncid -> ID of the newly created NetCDF file 
 ***           fid -> file pointer to the UM fields file
 ***         rflag -> denotes whether 32 or 64-bit output is desired
 ***                  (0->64 bit, 1->32 bit)
 ***          
 ***   Mark Cheeseman, NIWA
 ***   May 19, 2014
 ***/

void write_uninterpolated_fields( int ncid, FILE *fid ) {

     int     n, i, j=0, k, cnt, varid;
     size_t offset_3d[3], count_3d[3], offset_4d[4], count_4d[4];
     double *buf;
     char    name[45];

     for ( n=0; n<3; n++ ) {
         offset_3d[n] = 0;
         offset_4d[n] = 0;
     }
     offset_4d[3] = 0;

     count_3d[0] = 1;   
     count_4d[0] = 1;   
     count_4d[1] = 1;   
 
     for ( n=0; n<num_stored_um_fields; n++ ) {

         strcpy( name, stored_um_vars[n].name );
         i = nc_inq_varid( ncid, name, &varid );

       /*** Allocate & set the sizes of a 2D slice in the UM variable ***/
       /*** Remember that COUNT = COUNT[NT,NZ,NY,NX]                  ***/

         count_3d[1] = (size_t ) stored_um_vars[n].ny;
         count_3d[2] = (size_t ) stored_um_vars[n].nx;
         count_4d[2] = (size_t ) stored_um_vars[n].ny;
         count_4d[3] = (size_t ) stored_um_vars[n].nx;

       /*** Count the number of elements to be read for a single 2D data slice ***/

         cnt = (int ) stored_um_vars[n].nx*stored_um_vars[n].ny;
         buf = (double *) malloc( cnt*sizeof(double) );

       /*** For a 3D UM variable [NX,NY,NT], read in a single 2D data slice per data time   ***/
       /*** valid for this UM variable.  Then apply the appropriate interpolation and write ***/
       /*** the final field to hard disk.                                                   ***/

         if ( stored_um_vars[n].nz<2 ) {
               offset_3d[0] = 0;
               for ( k=0; k<stored_um_vars[n].nt; k++ ) {
                   fseek( fid, stored_um_vars[n].slices[k][0].location*wordsize, SEEK_SET );
                   fread( buf, wordsize, cnt, fid );
                   endian_swap( buf, cnt );
                   i = nc_put_vara_double( ncid, varid, offset_3d, count_3d, buf );
                   offset_3d[0]++;
               }
         } 

       /*** For a 4D UM variable [NX,NY,NZ,NT], read in a single 2D data slice per level and ***/
       /*** data time valid for this UM variable.  Then apply the appropriate interpolation  ***/
       /*** and write the final field to hard disk.                                          ***/

         else {
               offset_4d[1] = 0;
               for ( j=0; j<stored_um_vars[n].nz; j++ ) {
                   offset_4d[0] = 0;
                   for ( k=0; k<stored_um_vars[n].nt; k++ ) {
                       fseek( fid, stored_um_vars[n].slices[k][j].location*wordsize, SEEK_SET );
                       fread( buf, wordsize, cnt, fid );
                       endian_swap( buf, cnt );
                       i = nc_put_vara_double( ncid, varid, offset_4d, count_4d, buf );
                       offset_4d[0]++;
                   }
                   offset_4d[1]++;
               }
         }
         free( buf );

     }  // End of FOR LOOP

     return;
}

void write_uninterpolated_fields_flt( int ncid, FILE *fid ) {

     int     n, i, j=0, k, kk, cnt, varid;
     size_t offset_3d[3], count_3d[3], offset_4d[4], count_4d[4];
     double *buf;
     float  *fbuf=NULL;
     char    name[45];

     for ( n=0; n<3; n++ ) {
         offset_3d[n] = 0;
         offset_4d[n] = 0;
     }
     offset_4d[3] = 0;

     count_3d[0] = 1;   
     count_4d[0] = 1;   
     count_4d[1] = 1;   
 
     for ( n=0; n<num_stored_um_fields; n++ ) {

         strcpy( name, stored_um_vars[n].name );
         i = nc_inq_varid( ncid, name, &varid );

       /*** Allocate & set the sizes of a 2D slice in the UM variable ***/
       /*** Remember that COUNT = COUNT[NT,NZ,NY,NX]                  ***/

         count_3d[1] = (size_t ) stored_um_vars[n].ny;
         count_3d[2] = (size_t ) stored_um_vars[n].nx;
         count_4d[2] = (size_t ) stored_um_vars[n].ny;
         count_4d[3] = (size_t ) stored_um_vars[n].nx;

       /*** Count the number of elements to be read for a single 2D data slice ***/

         cnt = (int ) stored_um_vars[n].nx*stored_um_vars[n].ny;
         buf = (double *) malloc( cnt*sizeof(double) );
         fbuf = (float *) malloc( cnt*sizeof(float) ); 

       /*** For a 3D UM variable [NX,NY,NT], read in a single 2D data slice per data time   ***/
       /*** valid for this UM variable.  Then apply the appropriate interpolation and write ***/
       /*** the final field to hard disk.                                                   ***/

         if ( stored_um_vars[n].nz<2 ) {
               offset_3d[0] = 0;
               for ( k=0; k<stored_um_vars[n].nt; k++ ) {
                   fseek( fid, stored_um_vars[n].slices[k][0].location*wordsize, SEEK_SET );
                   fread( buf, wordsize, cnt, fid );
                   endian_swap( buf, cnt );
                   for ( j=0; j<cnt; j++ ) { fbuf[j] = (float ) buf[j]; }
                   i = nc_put_vara_float( ncid, varid, offset_3d, count_3d, fbuf );
                   offset_3d[0]++;
               }
         } 

       /*** For a 4D UM variable [NX,NY,NZ,NT], read in a single 2D data slice per level and ***/
       /*** data time valid for this UM variable.  Then apply the appropriate interpolation  ***/
       /*** and write the final field to hard disk.                                          ***/

         else {
               offset_4d[1] = 0;
               for ( kk=0; kk<stored_um_vars[n].nz; kk++ ) {
                   offset_4d[0] = 0;
                   for ( k=0; k<stored_um_vars[n].nt; k++ ) {
                       fseek( fid, stored_um_vars[n].slices[k][kk].location*wordsize, SEEK_SET );
                       fread( buf, wordsize, cnt, fid );
                       endian_swap( buf, cnt );
                       for ( j=0; j<cnt; j++ ) { fbuf[j] = (float ) buf[j]; }
                       i = nc_put_vara_float( ncid, varid, offset_4d, count_4d, fbuf );
                       offset_4d[0]++;
                   }
                   offset_4d[1]++;
               }
         }

         free( buf );
         free( fbuf ); 

     }  // End of FOR LOOP

     return;
}


int output_um_fields( int ncid, FILE *fid, int iflag, int rflag ) {

     int    ierr, n, varid; //, i,j,k;
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
     printf( "--------------------------------------------------------------\n" );
     printf( "  Stash_Code        Variable_Name        Dimensions [NX,NY,NZ]\n" );
     printf( "--------------------------------------------------------------\n" );

     for ( n=0; n<num_stored_um_fields; n++ ) {
         strcpy( name, stored_um_vars[n].name ); 
         printf( "     %5d %25s", stored_um_vars[n].stash_code, stored_um_vars[n].name );
         if ( (iflag==1)&&(stored_um_vars[n].ny!=int_constants[6]) ) { printf( "*" ); }
         else { printf( " " ); }
         printf( "    [%4d x %4d", stored_um_vars[n].nx, stored_um_vars[n].ny ); 
         if ( stored_um_vars[n].nz>1 ) { printf( " x %d", stored_um_vars[n].nz ); } 
         printf( "]\n" );
     }
     printf( "--------------------------------------------------------------\n\n" );
     if ( iflag==1 ) {
        printf( "  * denotes that field is to be interpolated onto a P grid.\n\n" );
     }

  /*
   * Write the UM field to hard disk one 2D data slice at a time. 
   *-------------------------------------------------------------------*/
     if ( iflag==1 ) { write_interpolated_fields( ncid, fid, rflag ); }
     else            { 
        if ( rflag==1 ) { write_uninterpolated_fields_flt( ncid, fid ); }
        else            { write_uninterpolated_fields( ncid, fid ); }
     }

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

