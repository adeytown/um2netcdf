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

double *interp_do_nothing( double *val, int nx, int ny );
double *u_to_p_point_interp_c_grid( double *val, int nx, int ny );
double *v_to_p_point_interp_c_grid( double *val, int nx, int ny );
double *b_to_c_grid_interp_u_points( double *val, int nx, int ny );
void   endian_swap_4bytes( void *ptr, int nchunk );
void   wgdos_unpack( FILE *fh, unsigned short nx, unsigned short ny, double *buf,
                     double mdi );

/***
 *** FILL_VARIABLES
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

int fill_variables( int ncid, FILE *fid, int iflag, int rflag ) {

     int      ierr, cnt, n, j, ndim, varid, i, num_z_levels;
     size_t   *offset, *count, dimlen;
     double   *buf;
     float    *fbuf;
     char     name[45];

     for ( n=0; n<num_stored_um_fields; n++ ) {

       //  strcpy( name,um_vars[stored_um_fields[n].xml_index].varname ); 
         strcpy( name, stored_um_fields[n].name ); 

       /** Determine ID & # of dimensions of variable **/
         ierr = nc_inq_varid( ncid, name, &varid );
         printf( "     %5d %25s     [%4d x %4d", stored_um_fields[n].stash_code, stored_um_fields[n].name, 
                                                   stored_um_fields[n].nx, stored_um_fields[n].ny ); 
 
         num_z_levels = stored_um_fields[n].num_slices/num_timesteps;
         if ( num_z_levels==1 ) { ndim = 3; } 
         else                   { ndim = 4; printf( " x %2d", num_z_levels ); }
         printf( "]\n" ); 

         offset = (size_t *) calloc( ndim,sizeof(size_t) );
         count = (size_t *) malloc( ndim*sizeof(size_t) );

       /** set the offsets for the upcoming NetCDF write **/
         count[0] = 1;
         if ( ndim==4 ) { count[1] = 1; }
         if ( iflag==0 ) { count[ndim-2] = stored_um_fields[n].ny; }
         else            { count[ndim-2] = int_constants[6]; }
         count[ndim-1] = stored_um_fields[n].nx; 

       /** Initialize function pointer to proper interpolation function **/
         if ( iflag==0 ) { field_interpolation = &interp_do_nothing; }
         else {
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
             
       /** Read in a 2D raw data array **/
             fseek( fid, stored_um_fields[n].slices[j].location*wordsize, SEEK_SET );
             if ( (stored_um_fields[n].slices[j].lbpack==0)||(stored_um_fields[n].slices[j].lbpack==3000) ) { 
                fread( buf, wordsize, cnt, fid ); 
                endian_swap( buf, cnt );
             } else if ( stored_um_fields[n].slices[j].lbpack==1 ) {
                wgdos_unpack( fid, stored_um_fields[n].nx, stored_um_fields[n].ny, buf, 
                              stored_um_fields[n].slices[j].mdi );
             }
             
       /** Perform interpolation if requested **/
             buf = field_interpolation( buf, stored_um_fields[n].nx, stored_um_fields[n].ny );            

       /** Write the raw data to disk **/
             if ( rflag==0 ) { 
                ierr = nc_put_vara_double( ncid, varid, offset, count, buf ); 
             } else {
                fbuf = (float *) malloc( cnt*sizeof(float) );
                for ( i=0; i<cnt; i++ ) { fbuf[i] = (float ) buf[i]; }
                ierr = nc_put_vara_float( ncid, varid, offset, count, fbuf ); 
                free( fbuf );   
             }  

       /** Update the offset counter **/
             if ( ndim==3 ) { offset[0]++;  }
             else {
                 offset[1]++;
                 if ( offset[1] == num_z_levels ) {
                    offset[0]++;
                    offset[1] = 0;
                 }   
             }

         }  
         free( buf ); 
         free( count );
         free( offset );  

     }
     printf( "\n" ); 

    /*** Output the coefficients for the ETA arrays ***/

     dimlen = (size_t ) header[110];

     fbuf = (float *) malloc( dimlen*sizeof(float) );
     for ( j=0; j<dimlen; j++ ) { fbuf[j] = level_constants[0][j]; } 
     ierr = nc_inq_varid( ncid, "eta_theta", &varid );   
     ierr = nc_put_var_float( ncid, varid, fbuf );
     free( fbuf ); 

     fbuf = (float *) malloc( (dimlen-1)*sizeof(float) );
     for ( j=0; j<dimlen-1; j++ ) { fbuf[j] = level_constants[1][j]; } 
     ierr = nc_inq_varid( ncid, "eta_rho", &varid );   
     ierr = nc_put_var_float( ncid, varid, fbuf );

     return 1;
}

