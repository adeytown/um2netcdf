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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <netinet/in.h>
#include "field_def.h"

/***
 *** ENDIAN_SWAP_#BYTES 
 ***
 *** Subroutine that performs a byte-swap on #-byte long input words. This
 *** process is used to convert from big-endian to little-endian data.
 *** 
 ***  INPUT: ptr -> pointer to the data array whose contents are to be 
 ***                byte-swapped
 ***         N   -> # of elements of the input data array
 ***
 ***   Mark Cheeseman, NIWA
 ***   November 29, 2013
 ***/

void endian_swap_8bytes( void *ptr, int N ) {

     int      i;
     char *p, t;

     for ( i=0; i<N; i++ ) {
         p = (char*) ptr + 8*i;
         t=p[7]; p[7]=p[0]; p[0]=t;
         t=p[6]; p[6]=p[1]; p[1]=t;
         t=p[5]; p[5]=p[2]; p[2]=t;
         t=p[4]; p[4]=p[3]; p[3]=t; 
     }
     return;
}

void endian_swap_4bytes( void *ptr, int nchunk ) {

     int i;
     char *p, t;

     for ( i=0; i<nchunk; i++ ) {
         p = (char*) ptr + 4*i;
         t=p[3]; p[3]=p[0]; p[0]=t;
         t=p[2]; p[2]=p[1]; p[1]=t;
     }
     return;
}

void endian_swap_2bytes( void *ptr, int nchunk ) {

     int i;
     char *p, t;

     for ( i=0; i<nchunk; i++ ) {
         p = (char*) ptr + 2*i;
         t=p[1]; p[1]=p[0]; p[0]=t;
     }
     return;
}

void no_endian_swap( void *ptr, int nchunk ) { return; }


/***
 *** GET_FILE_ENDIANNESS_WORDSIZE 
 ***
 *** Function that determines the appropriate wordsize and endianness of the 
 *** user-supplied input UM fields file. 
 ***
 *** INPUT:  fh -> file handle/pointer of the input UM fields file
 ***
 ***   Mark Cheeseman, NIWA
 ***   January 16, 2014
 ***/

int get_file_endianness_wordsize( FILE *fh ) {

    int word_size;

    for ( word_size=4; word_size<9; word_size+=4 ) {

        fseek( fh, 0, SEEK_SET );
        fread( header, word_size, 256, fh );

        if ( (header[1]==1)&&(header[4]==3)&&(header[8]==3)&&(header[150]==64) ) {
           printf( "No swapping needed\n" );
           endian_swap = &no_endian_swap;
           return word_size;
        }
  
        if ( word_size==4 ) { endian_swap = &endian_swap_4bytes; }
        else                { endian_swap = &endian_swap_8bytes; }

        endian_swap( header, 256 );
        if ( (header[1]==1)&&(header[4]==3)&&(header[8]==3)&&(header[150]==64) ) {
//           printf( "Swapping needed\n" );
           return word_size;
        }
 
    }
    word_size = -1; 
    return wordsize;
}


/***
 *** CHECK_UM_FILE 
 ***
 *** Subroutine that opens the user-specified UM input file and determines 
 *** its type, its endianness and wordsize. 
 ***
 ***   Mark Cheeseman, NIWA
 ***   November 29, 2013
 ***/

int check_um_file( char *filename, int rflag ) {

     unsigned short *temp, num_lbproc;
     int    i, j, k, kk, ind, flag, cnt, nrec, *lbprocs;
     long   nn, **tmp, **lookup;
     size_t n;
     FILE   *fid;
     struct tm first, next;
     char   varname[60];

/**
 ** Attempt to open the file 
 **---------------------------------------------------------------------------*/
     fid = fopen( filename, "r" );
     if ( fid==NULL ) { return 0; }

/**
 ** Initialize the word size of the input UM fields file 
 **---------------------------------------------------------------------------*/
     wordsize = get_file_endianness_wordsize( fid );

/**
 ** Read in the REAL CONSTANTS array 
 **---------------------------------------------------------------------------*/
     fseek( fid, (header[104]-1)*wordsize, SEEK_SET );
     fread( real_constants, wordsize, 6, fid );
     endian_swap( real_constants,6 );
   
/**
 ** Read in the INTEGER CONSTANTS array 
 **---------------------------------------------------------------------------*/
     fseek( fid, (header[99]-1)*wordsize, SEEK_SET );
     fread( int_constants, wordsize, 46, fid );
     endian_swap( int_constants,46 );

/**
 ** Read in the LEVEL DEPENDENT CONSTANTS array 
 **---------------------------------------------------------------------------*/
     fseek( fid, (header[109]-1)*wordsize, SEEK_SET );
     level_constants = (double **) malloc( header[111]*sizeof(long *) );

     for ( nrec=0; nrec<header[111]; nrec++ ) {
         level_constants[nrec] = (double *) malloc( header[110]*sizeof(long) );
         n = fread( level_constants[nrec], wordsize, header[110], fid );
         endian_swap( level_constants[nrec],header[110] ); 
     }

/**
 ** Read in the LOOKUP table 
 **---------------------------------------------------------------------------*/
     fseek( fid, (header[149]-1)*wordsize, SEEK_SET );
     tmp = (long **) malloc( header[151]*sizeof(long *) );

     for ( nrec=0; nrec<header[151]; nrec++ ) {
         tmp[nrec] = (long *) malloc( header[150]*sizeof(long) );
         n = fread( tmp[nrec], wordsize, header[150], fid );
         endian_swap( tmp[nrec],header[150] ); 
     }

/**
 ** Determine how many read entries belong to valid UM data slices [NUM_UM_VARS] 
 **---------------------------------------------------------------------------*/
     cnt = 0;
     for ( nrec=0; nrec<header[151]; nrec++ ) {
         if ( tmp[nrec][28]!=-99 ) { 
            cnt++; 
         }
     }
     num_um_vars = cnt;

/**
 ** Gather all the valid lookup entries into a condensed array (called LOOKUP)) 
 **---------------------------------------------------------------------------*/
     lookup = (long **) malloc( cnt*sizeof(long *) );

     cnt = -1;
     for ( nrec=0; nrec<header[151]; nrec++ ) {
         if ( tmp[nrec][28]!=-99 ) {
            cnt++;
            lookup[cnt] = (long *) malloc( header[150]*sizeof(long) );
            for ( n=0; n<header[150]; n++ ) {
                lookup[cnt][n] = tmp[nrec][n]; 
            } 
         }
     }

     for ( nrec=0; nrec<header[151]; nrec++ )
           free( tmp[nrec] );
     free( tmp );

/**
 ** Determine # of timesteps covered in the input UM fields file 
 **---------------------------------------------------------------------------*/
     num_timesteps = 0;
     for ( nrec=0; nrec<num_um_vars; nrec++ ) {
         nn = (lookup[0][41]-lookup[nrec][41]) + (lookup[0][32]-lookup[nrec][32]);
         if ( nn==0 ) {  num_timesteps++; }
     } 

/**
 ** If the user has not requested specific stash codes, determine # of unique 
 ** UM variables found in the input UM fields file 
 **---------------------------------------------------------------------------*/
     if ( num_stored_um_fields==0) {
        temp = (unsigned short *) malloc( 250*sizeof(unsigned short) );

        flag = 0;
        for ( j=0; j<num_um_vars; j++ ) {

            if ( num_stored_um_fields>0 ) {
               flag = 0;
               for ( i=0; i<num_stored_um_fields; i++ )
                   if ( temp[i]==lookup[j][41] ) { flag=1; }
            }

            if ( flag==1 ) { continue; }  

            temp[num_stored_um_fields] = (unsigned short ) lookup[j][41];
            num_stored_um_fields++;
        }

/**
 ** Allocate UM variable data structure & copy list of unique stash codes to it 
 **---------------------------------------------------------------------------*/
        stored_um_fields = (um_variable* ) malloc( num_stored_um_fields*sizeof(um_variable) );
 
        for ( j=0; j<num_stored_um_fields; j++ ) {
            stored_um_fields[j].stash_code = temp[j];
            stored_um_fields[j].num_slices = 0;
        }
        free( temp );
     }

/**
 ** Determine # of 2D data slices belonging to each unique UM variable. 
 **---------------------------------------------------------------------------*/
     for ( j=0; j<num_stored_um_fields; j++ ) {
         flag = 0;
         for ( i=0; i<num_um_vars; i++ ) {
             if ( (unsigned short ) lookup[i][41]==stored_um_fields[j].stash_code ) {
                stored_um_fields[j].num_slices++; 
                flag = 1;
             }
         }
         if ( flag==0 ) { num_stored_um_fields--; }
     }

     if ( num_stored_um_fields<=0 ) {
        printf( "ERROR: supplied stash codes not found in input UM fields file\n\n" );
        exit(1);
     }

/**
 ** Record the ID of each 2d data slice belonging to each unique UM variable 
 **---------------------------------------------------------------------------*/
     for ( j=0; j<num_stored_um_fields; j++ ) {
         stored_um_fields[j].slices       = (um_dataslice * ) malloc( stored_um_fields[j].num_slices*sizeof(um_dataslice) );
         for ( i=0; i<stored_um_fields[j].num_slices; i++ ) { stored_um_fields[j].slices[i].id = 9999; }
     }

     for ( j=0; j<num_stored_um_fields; j++ ) {
     for ( i=0; i<num_um_vars; i++ ) {
         if ( (unsigned short ) lookup[i][41]==stored_um_fields[j].stash_code ) {
            for ( k=0; k<stored_um_fields[j].num_slices; k++ ) {
                if ( stored_um_fields[j].slices[k].id==9999 ) { 
                   stored_um_fields[j].slices[k].id       = (unsigned short ) i; 
                   stored_um_fields[j].slices[k].size     = lookup[i][14]; 
                   stored_um_fields[j].coordinates        = (unsigned short ) lookup[i][15]; 
                   stored_um_fields[j].slices[k].location = lookup[i][28]; 
                   stored_um_fields[j].slices[k].reclength= lookup[i][29]; 
                   stored_um_fields[j].slices[k].level    = (unsigned short ) lookup[i][32]; 
                   stored_um_fields[j].slices[k].lbproc   = lookup[i][24]; 
                   stored_um_fields[j].slices[k].lbpack   = (unsigned short ) lookup[i][20]; 
                   stored_um_fields[j].slices[k].mdi      = (double ) lookup[i][62]; 
                   break; 
                }
                if ( stored_um_fields[j].slices[k].id==(unsigned short ) i ) {  break; }
            }
         }
     } 
     } 

/**
 ** Fill the remaining attributes of the UM variable data structure   
 **---------------------------------------------------------------------------*/
     for ( j=0; j<num_stored_um_fields; j++ ) {
         i = stored_um_fields[j].slices[0].id;

         stored_um_fields[j].ny     = (unsigned short ) lookup[i][17];
         stored_um_fields[j].nx     = (unsigned short ) lookup[i][18];
         stored_um_fields[j].lbvc   = (unsigned short ) lookup[i][25];
         stored_um_fields[j].validity.tm_year = (int ) lookup[i][0];
         stored_um_fields[j].validity.tm_mon  = (int ) lookup[i][1];
         stored_um_fields[j].validity.tm_mday = (int ) lookup[i][2];
         stored_um_fields[j].validity.tm_hour = (int ) lookup[i][3];
         stored_um_fields[j].validity.tm_min  = (int ) lookup[i][4];
         stored_um_fields[j].validity.tm_sec  = (int ) lookup[i][5];
         stored_um_fields[j].lbproc           =        lookup[i][24]; 
         if ( rflag==0 ) {
            if ( lookup[i][38]==1 ) { stored_um_fields[j].vartype = NC_DOUBLE; }
            else                    { stored_um_fields[j].vartype = NC_LONG; }
         } else {
            if ( lookup[i][38]==1 ) { stored_um_fields[j].vartype = NC_FLOAT; }
            else                    { stored_um_fields[j].vartype = NC_INT; }
         }
     }

/**
 ** Locate the index of the XML file where each unique UM variable is given
 ** its CF-compliant description.   
 **---------------------------------------------------------------------------*/
     for ( i=0; i<num_stored_um_fields; i++ ) {
     for ( j=0; j<num_xml_vars; j++ ) {
         if ( stored_um_fields[i].stash_code == 1000*um_vars[j].section + um_vars[j].code ) {
            stored_um_fields[i].xml_index = j;
            strcpy( stored_um_fields[i].name, um_vars[j].varname );
            stored_um_fields[i].grid_type = um_vars[j].umgrid;
            stored_um_fields[i].accum = um_vars[j].accum;
            break;
         }
     }
     }

     for ( k=0; k<num_stored_um_fields; k++ ) {

/**
 ** How many unique LBPROC values does each UM variable possess? 
 **---------------------------------------------------------------------------*/
         lbprocs = (int*) malloc( stored_um_fields[k].num_slices*sizeof(int) );
          
         for ( i=0; i<stored_um_fields[k].num_slices; i++ ) { 
             lbprocs[i] = stored_um_fields[k].slices[i].lbproc;
         }

         for ( i=0;   i<stored_um_fields[k].num_slices; i++ ) { 
         for ( j=i+1; j<stored_um_fields[k].num_slices; j++ ) { 
             if ( lbprocs[i]==lbprocs[j] ) lbprocs[j] = -1;
         }
         }
          
         num_lbproc = 0;
         for ( i=0; i<stored_um_fields[k].num_slices; i++ ) { 
             if ( lbprocs[i]!=-1 ) { 
                num_lbproc++;
             }
         } 
         free( lbprocs );
    
   /*
    * If only 1 unique value of LBPROC exists, set LBPROC for the whole variable
    *------------------------------------------------------------------------*/ 
         if ( num_lbproc==1 ) { stored_um_fields[k].lbproc = stored_um_fields[k].slices[0].lbproc; }

   /*
    * If multiple unique values of LBPROC exist, spawn NUM_LBPROCS-1 additional
    * UM variables  
    *------------------------------------------------------------------------*/
         if ( num_lbproc>1 ) {

         /* Extend the current array of UM variables */
            num_lbproc--;
            cnt = num_stored_um_fields;
            num_stored_um_fields += num_lbproc;
            stored_um_fields = (um_variable* ) realloc( stored_um_fields, num_stored_um_fields*sizeof(um_variable) );

            for ( i=0; i<num_lbproc; i++ ) {

            /* Set whole variable attributes */
                stored_um_fields[i+cnt].stash_code  = stored_um_fields[k].stash_code;
                stored_um_fields[i+cnt].xml_index   = stored_um_fields[k].xml_index;
                stored_um_fields[i+cnt].grid_type   = stored_um_fields[k].grid_type;
                stored_um_fields[i+cnt].accum       = stored_um_fields[k].accum;
                stored_um_fields[i+cnt].vartype     = stored_um_fields[k].vartype;
                stored_um_fields[i+cnt].ny          = stored_um_fields[k].ny;
                stored_um_fields[i+cnt].nx          = stored_um_fields[k].nx;
                stored_um_fields[i+cnt].lbvc        = stored_um_fields[k].lbvc;
                stored_um_fields[i+cnt].coordinates = stored_um_fields[k].coordinates;
                stored_um_fields[i+cnt].validity.tm_year = stored_um_fields[k].validity.tm_year;
                stored_um_fields[i+cnt].validity.tm_mon  = stored_um_fields[k].validity.tm_mon;
                stored_um_fields[i+cnt].validity.tm_mday = stored_um_fields[k].validity.tm_mday;
                stored_um_fields[i+cnt].validity.tm_hour = stored_um_fields[k].validity.tm_hour;
                stored_um_fields[i+cnt].validity.tm_min  = stored_um_fields[k].validity.tm_min;
                stored_um_fields[i+cnt].validity.tm_sec  = stored_um_fields[k].validity.tm_sec;

            /* Set the name of the newly added UM variable */
                strcpy( stored_um_fields[i+cnt].name, stored_um_fields[k].name );

            /* Set the # of slices for this new UM variable */
                stored_um_fields[i+cnt].num_slices = stored_um_fields[k].num_slices / (num_lbproc+1);
                stored_um_fields[i+cnt].slices     = (um_dataslice * ) malloc( stored_um_fields[i+cnt].num_slices*sizeof(um_dataslice) );

            /* Set the attributes for each data slice of the new UM variable */
                for ( kk=0; kk<stored_um_fields[i+cnt].num_slices; kk++ ) {
                    ind = kk*num_timesteps + i + 1;
                    stored_um_fields[i+cnt].slices[kk].id       = stored_um_fields[k].slices[ind].id;
                    stored_um_fields[i+cnt].slices[kk].size     = stored_um_fields[k].slices[ind].size;
                    stored_um_fields[i+cnt].slices[kk].location = stored_um_fields[k].slices[ind].location;
                    stored_um_fields[i+cnt].slices[kk].reclength= stored_um_fields[k].slices[ind].reclength;
                    stored_um_fields[i+cnt].slices[kk].level    = stored_um_fields[k].slices[ind].level;
                    stored_um_fields[i+cnt].slices[kk].lbproc   = stored_um_fields[k].slices[ind].lbproc; 
                    stored_um_fields[i+cnt].slices[kk].lbpack   = stored_um_fields[k].slices[ind].lbpack; 
                    stored_um_fields[i+cnt].slices[kk].mdi      = stored_um_fields[k].slices[ind].mdi; 
                }
                stored_um_fields[i+cnt].lbproc  = stored_um_fields[i+cnt].slices[0].lbproc;

            /* Set the name of the newly added UM variable */
        /*        kk = stored_um_fields[i+cnt].slices[0].lbproc; 
                if ( kk==0 )    { strcpy( varname, "inst_" ); }
                if ( kk==64 )   { strcpy( varname, "mean_" ); }
                if ( kk==4096 ) { strcpy( varname, "min_" ); }
                if ( kk==8192 ) { strcpy( varname, "max_" ); }

                if ( kk==128 ) {
                   if (stored_um_fields[i+cnt].accum==0 ) { strcpy( varname, "mean_" ); }
                   else                                   { strcpy( varname, "sum_" ); }
                }
                strcat( varname, stored_um_fields[k].name );
                strcpy( stored_um_fields[i+cnt].name, varname ); */

            } 

            /* Rename the original UM variable that was split */
     //       stored_um_fields[k].lbproc = stored_um_fields[k].slices[0].lbproc;

      /*      kk = stored_um_fields[k].slices[0].lbproc; 
            if ( kk==0 )    { strcpy( varname, "inst_" ); }
            if ( kk==64 )   { strcpy( varname, "mean_" ); }
            if ( kk==4096 ) { strcpy( varname, "min_" ); }
            if ( kk==8192 ) { strcpy( varname, "max_" ); }

            if ( kk==128 ) {
               if (stored_um_fields[k].accum==0 ) { strcpy( varname, "mean_" ); }
               else                               { strcpy( varname, "sum_" ); }
            }
            strcat( varname, stored_um_fields[k].name );
            strcpy( stored_um_fields[k].name, varname ); */

            /* Resize the # of slices */
            stored_um_fields[k].num_slices = stored_um_fields[k].num_slices / (num_lbproc+1);
            for ( kk=0; kk<stored_um_fields[k].num_slices; kk++ ) {
                ind = kk*num_timesteps;
                stored_um_fields[k].slices[kk].id       = stored_um_fields[k].slices[ind].id;
                stored_um_fields[k].slices[kk].size     = stored_um_fields[k].slices[ind].size;
                stored_um_fields[k].slices[kk].location = stored_um_fields[k].slices[ind].location;
                stored_um_fields[k].slices[kk].level    = stored_um_fields[k].slices[ind].level;
                stored_um_fields[k].slices[kk].lbproc   = stored_um_fields[k].slices[ind].lbproc;
            }
            stored_um_fields[k].slices = (um_dataslice * ) realloc( stored_um_fields[k].slices, 
                                                                    stored_um_fields[k].num_slices*sizeof(um_dataslice) );

         }

     }

/**
 ** Add an appropriate prefix to the UM variable's name if post-processing was
 ** performed on the field. 
 **---------------------------------------------------------------------------*/
     for ( i=0; i<num_stored_um_fields; i++ ) {
         kk = stored_um_fields[i].lbproc;
         if ( kk!=0 ) { 
            if ( kk==64 )   { strcpy( varname, "mean_" ); }
            if ( kk==4096 ) { strcpy( varname, "min_"  ); }
            if ( kk==8192 ) { strcpy( varname, "max_"  ); }
            if ( kk==128 ) {
               if (stored_um_fields[i].accum==0 ) { strcpy( varname, "mean_" ); }
               else                               { strcpy( varname, "sum_" ); }
            }
            strcat( varname, stored_um_fields[i].name );
            strcpy( stored_um_fields[i].name, varname ); 
         }
     } 

/**
 ** Compute the time difference in hours for each successive timestep after 
 ** the first recorded one. 
 **---------------------------------------------------------------------------*/
     for ( i=0; i<num_stored_um_fields; i++ ) {
         stored_um_fields[i].time_offsets = (float * ) calloc( num_timesteps,sizeof(float) );
         for ( j=0; j<num_timesteps; j++ ) {

            cnt = j*stored_um_fields[i].num_slices/num_timesteps;
            k = stored_um_fields[i].slices[cnt].id;
            next.tm_year = lookup[k][0];
            next.tm_mon  = lookup[k][1];
            next.tm_mday = lookup[k][2];
            next.tm_hour = lookup[k][3];
            next.tm_min  = lookup[k][4];
            next.tm_sec  = lookup[k][5];

            stored_um_fields[i].time_offsets[j] = (float ) (difftime ( mktime(&next),
                                                                    mktime(&stored_um_fields[i].validity) ) 
                                                         / 3600.0);
         
        }
     }

     for ( i=0; i<num_stored_um_fields; i++ ) {
/**
 ** For UM variables that are a temporal mean or accummulation, compute the 
 ** duration of each mean/accummulatoin. 
 **---------------------------------------------------------------------------*/
         if ( (stored_um_fields[i].lbproc==32)||(stored_um_fields[i].lbproc==128)||
              (stored_um_fields[i].lbproc==4096)||(stored_um_fields[i].lbproc==8192) ) {

            stored_um_fields[i].time_bnds = (float **) malloc( num_timesteps*sizeof(float *) );
            for ( k=0; k<num_timesteps; k++ )
                stored_um_fields[i].time_bnds[k] = (float *) calloc( 2,sizeof(float) );

            for ( j=0; j<num_timesteps; j++ ) {

                cnt = j*stored_um_fields[i].num_slices/num_timesteps;
                k = stored_um_fields[i].slices[cnt].id;

                first.tm_year = lookup[k][0];
                first.tm_mon  = lookup[k][1];
                first.tm_mday = lookup[k][2];
                first.tm_hour = lookup[k][3];
                first.tm_min  = lookup[k][4];
                first.tm_sec  = lookup[k][5];

                next.tm_year = lookup[k][6];
                next.tm_mon  = lookup[k][7];
                next.tm_mday = lookup[k][8];
                next.tm_hour = lookup[k][9];
                next.tm_min  = lookup[k][10];
                next.tm_sec  = lookup[k][11];

                stored_um_fields[i].time_bnds[j][0] = (float ) (difftime ( mktime(&first),
                                               mktime(&stored_um_fields[i].validity) ) / 3600.0);
                stored_um_fields[i].time_bnds[j][1] = (float ) (difftime ( mktime(&next),
                                               mktime(&stored_um_fields[i].validity) ) / 3600.0);

            }
         }

/**
 ** For UM variables that are a spatial mean or accummulation, set the end 
 ** points for the accummulation operation.
 **---------------------------------------------------------------------------*/
         if ( (stored_um_fields[i].lbproc==8)||(stored_um_fields[i].lbproc==64) ) {
            stored_um_fields[i].space_bnds[0] = real_constants[3]; 
            stored_um_fields[i].space_bnds[1] = real_constants[3] + real_constants[0]*((float ) stored_um_fields[j].nx);
         } else if ( stored_um_fields[i].lbproc==16 ) {
            stored_um_fields[i].space_bnds[0] = real_constants[2]; 
            stored_um_fields[i].space_bnds[1] = real_constants[2] + real_constants[1]*((float ) stored_um_fields[j].ny);
         }  

     }

/**
 ** Free memory used by the LOOKUP 2D array
 **---------------------------------------------------------------------------*/
     for ( nrec=0; nrec<num_um_vars; nrec++ )
         free( lookup[nrec] );
     free( lookup );

     fclose( fid ); 
     return 1; 
}
