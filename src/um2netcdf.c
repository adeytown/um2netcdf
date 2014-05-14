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


/***
 *** UM2NetCDF
 ***
 *** Re-written implementation of the UM2NetCDF utility initially 
 *** developed at NIWA.
 ***
 ***    Mark Cheeseman, NIWA
 ***    Nov 28, 2013
 ***/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <netcdf.h>
#include "field_def.h"

/** Function prototypes **/

void status_check( int status, char *message );
void usage();
int read_stash_file( char *filename );
int check_um_file( char *filename, int rflag);
int create_netcdf_file( char *um_file, int iflag, int rflag );
int fill_netcdf_file( int ncid, char *filename, int iflag, int rflag );

int main( int argc, char *argv[] ) {

     int ncid, status, c, iflag, rflag, sflag, n, i;

 /*
  * Check if the user has included the correct number of commandline arguments
  *---------------------------------------------------------------------------*/ 
     if ( argc<3 ) {
        usage(); 
        status = 0;
        status_check( status, "ERROR: too few commandline arguments" );
     } 

     iflag = 0;
     rflag = 0;
     sflag = 0;
     num_stored_um_fields = 0;
     while ( (c = getopt(argc,argv,"irs")) != EOF ) {
           switch(c) {
               case 'i':
                       iflag = 1;
                       break;
               case 'r':
                       rflag = 1;
                       break;
               case 's':
                       sflag = 1;
                       break;
           }
     } 

 /*
  * If the user has requested specific stash codes (-s option), read them in. 
  *---------------------------------------------------------------------------*/ 
     if ( sflag==1 ) {
        num_stored_um_fields = argc - 4 - iflag - rflag;
        stored_um_vars = (new_um_variable *) malloc( num_stored_um_fields*sizeof(new_um_variable) );
        for ( n=0; n<num_stored_um_fields; n++ ) {
            stored_um_vars[n].stash_code = atoi(argv[2+iflag+rflag+n]);
            stored_um_vars[n].nt = 0;
            stored_um_vars[n].nz = 0;
        }
     }

 /*
  * Read in the variable definitions in the XML stash file 
  *---------------------------------------------------------------------------*/ 
     status = read_stash_file( argv[argc-1] );
     status_check( status, "ERROR: could not parse stash file" ); 

 /*
  * Determine the type of input UM file specified by user, its endianness
  * and word size. 
  *---------------------------------------------------------------------------*/ 
     status = check_um_file( argv[argc-2], rflag ); 
     status_check( status, "ERROR: could not determine UM filetype" );

     printf( "\n===========================================================\n" );
     printf( "                          UM2NetCDF\n" );
     printf( "===========================================================\n\n" );
     printf( "Input UM Fields File\n" );
     printf( "-----------------------------------------------------------\n" );
     printf( "   Filename  : %s\n", argv[argc-2] );
     printf( "   Wordsize  : %d\n\n", wordsize );
     printf( "UM Stash Codes->CF Metadata XML File\n" );
     printf( "-----------------------------------------------------------\n" );
     printf( "   Filename  : %s\n\n", argv[argc-1] );
     
 /*
  * Create a NetCDF file to hold the UM data 
  *---------------------------------------------------------------------------*/ 
     ncid = create_netcdf_file( argv[argc-2], iflag, rflag );  

 /*
  * Write the UM data into the NetCDF file 
  *---------------------------------------------------------------------------*/ 
     status = fill_netcdf_file( ncid, argv[argc-2], iflag, rflag );  
     status_check( status, "ERROR: data write to NetCDF file failed" );

 /*
  * Free allocated memory 
  *---------------------------------------------------------------------------*/ 
     for ( i=0; i<num_stored_um_fields; i++ ) {
         free( stored_um_vars[i].times );
         for ( n=0; n<stored_um_vars[i].nt; n++ )
             free( stored_um_vars[i].slices[n] );
         free( stored_um_vars[i].slices );
     }
 
     free( stored_um_vars );  

     return 1;
}
