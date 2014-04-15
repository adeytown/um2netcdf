/**============================================================================
                   U M 2 N e t C D F  V e r s i o n 2 . 0
                   --------------------------------------
 
       Main author: Mark Cheeseman
                    National Institute of Water & Atmospheric Research (Ltd)
                    Wellington, New Zealand
                    April 2014

       UM2NetCDF is free software: you can redistribute it and/or modify
       it under the terms of the GNU General Public License as published by
       the Free Software Foundation, either version 3 of the License, or
       any later version.
 
       UM2NetCDF is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

       A copy of the GNU General Public License can be found in the main 
       UM2NetCDF directory.  Alternatively, please see 
       <http://www.gnu.org/licenses/>.
**============================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "field_def.h"

/** Function Prototypes **/

double *interp_do_nothing( double *val, int nx, int ny );
double *u_to_p_point_interp_c_grid( double *val, int nx, int ny );
double *v_to_p_point_interp_c_grid( double *val, int nx, int ny );
double *b_to_c_grid_interp_u_points( double *val, int nx, int ny );
void   endian_swap_4bytes( void *ptr, int nchunk );


/***
 *** COPY_INPUT_FILE
 ***
 *** Subroutine that creates a copy of the input UM data file.
 ***
 ***   Mark Cheeseman, NIWA
 ***   April 14, 2014
 ***/

void copy_input_file( char *um_file, int iflag ) {

     int     pos, n, j, num_z_levels, cnt;
     char    ch, filename[50], *str, *dest;
     FILE   *fid_1, *fh;
     double *buf;

 /*
  * 1) Create an appropriate name for the output file
  *---------------------------------------------------------------------------*/
     dest = strstr( um_file, ".um" );
     pos = dest - um_file;

     str = malloc( 1 + strlen(um_file));
     strcpy( str, um_file ); 
     str[pos] = '\0'; 

     snprintf( filename, sizeof filename, "%s.um_new", str );
     free( str );

 /*
  * 2) Copy the contents of the input UM data file
  *---------------------------------------------------------------------------*/
     fid_1 = fopen( um_file, "r" );
     if ( fid_1==NULL ) {
        printf( "ERROR: Input UM file %s not found\n", um_file );
        exit(1);
     } 

     fh = fopen( filename, "w" );
     if ( fh==NULL ) {
        printf( "ERROR: could not open file %s\n", filename );
        exit(1);
     }

     while ( (ch=fgetc(fid_1)) != EOF )
           fputc( ch, fh );
 
     fclose( fid_1 );

     if ( iflag==0 ) { field_interpolation = &interp_do_nothing; }
     for ( n=0; n<num_stored_um_fields; n++ ) {

     /** Output the name and dimensions of the UM variable to be modified **/
         num_z_levels = stored_um_fields[n].num_slices/num_timesteps;
         printf( "     %5d     %s        [%4d x %4d", stored_um_fields[n].stash_code, stored_um_fields[n].name,
                                                      stored_um_fields[n].nx, stored_um_fields[n].ny );
         if ( num_z_levels>1 ) { printf( " x %2d", num_z_levels ); }
         printf( "]\n" );

     /** Initialize function pointer to proper interpolation function **/
         if ( iflag==1 ) { 
            if ( stored_um_fields[n].grid_type==11 ) {
               field_interpolation = &b_to_c_grid_interp_u_points;
            }
            else if ( stored_um_fields[n].grid_type==18 ) {
               field_interpolation = &u_to_p_point_interp_c_grid;
            }
            else if ( stored_um_fields[n].grid_type==19 ) {
               field_interpolation = &v_to_p_point_interp_c_grid;
            }
            else { field_interpolation = &interp_do_nothing; }
         }

     /** Create a buffer to hold the raw data values **/
         cnt = stored_um_fields[n].nx*stored_um_fields[n].ny;
         buf = (double *) calloc( cnt,sizeof(double) );

         for ( j=0; j<stored_um_fields[n].num_slices; j++ ) {
      
         /** Read in 2D slice of data for UM variable **/
             fseek( fh, stored_um_fields[n].slices[j].location*wordsize, SEEK_SET );
             if ( (stored_um_fields[n].slices[j].lbpack==0)||(stored_um_fields[n].slices[j].lbpack==3000) ) {
                fread( buf, wordsize, cnt, fh );
                endian_swap( buf, cnt );
             } else if ( stored_um_fields[n].slices[j].lbpack==1 ) {
                 printf( "WARNING: WGDOS packing not currently support in um2netcdf\n" );
                 continue;
             }

         /** Perform interpolation (if requested by user) **/
             buf = field_interpolation( buf, stored_um_fields[n].nx, stored_um_fields[n].ny );

 /**++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 *
 *   WHERE USER PLACES HIS/HER MODIFICATIONS TO FIELD VALUES
 * 
 **++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/   

         /** Place the modified 2D data slice back into the file **/
             fseek( fh, stored_um_fields[n].slices[j].location*wordsize, SEEK_SET );
             endian_swap( buf, cnt );
             fwrite( buf, wordsize, cnt, fh );

         }
         free( buf );
     }
     printf( "\n" );
     return;
}
