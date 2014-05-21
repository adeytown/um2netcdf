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
#include <omp.h>
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

#pragma omp parallel shared(N,ptr) private(i,p,t) 
{
#pragma omp for 
     for ( i=0; i<N; i++ ) {
         p = (char*) ptr + 8*i;
         t=p[7]; p[7]=p[0]; p[0]=t;
         t=p[6]; p[6]=p[1]; p[1]=t;
         t=p[5]; p[5]=p[2]; p[2]=t;
         t=p[4]; p[4]=p[3]; p[3]=t; 
     }
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

  /*
   * Check first if file was created on a big endian 8-byte word platform
   * like the IBM Power system
   *----------------------------------------------------------------------*/
    word_size = 8;
    fseek( fh, 0, SEEK_SET );
    fread( header, word_size, 256, fh );

    if ( (header[1]==1)&&(header[150]==64) ) {
       endian_swap = &no_endian_swap;
//       printf( "No swapping: %ld %ld\n", header[1], header[150] );
       return word_size;
    } else {
       endian_swap = &endian_swap_8bytes;
       endian_swap( header, 256 );
//       printf( "Big Endian swapping: %ld %ld\n", header[1], header[150] );
       if ( (header[1]==1)&&(header[150]==64) ) {
          return word_size;
       }
    }

  /*
   * Then check if file was created on a little endian 8-byte word platform
   * like a x86 system
   *----------------------------------------------------------------------*/
    word_size = 4;
    fseek( fh, 0, SEEK_SET );
    fread( header, word_size, 256, fh );
    printf( "%ld %ld\n", header[1], header[150] );

    if ( (header[1]==1)&&(header[150]==64) ) {
       endian_swap = &no_endian_swap;
//       printf( "No swapping: %ld %ld\n", header[1], header[150] );
       return word_size;
    } else {
       endian_swap = &endian_swap_4bytes;
       endian_swap( header, 256 );
//       printf( "Little Endian swapping: %ld %ld\n", header[1], header[150] );
       if ( (header[1]==1)&&(header[150]==64) ) {
          return word_size;
       }
    }

    word_size = -1;
    return word_size;

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

     unsigned short *temp, num_lbproc, *temp_id=NULL, ntmp;
     int    i, j, k, kk, ind, flag, cnt, nrec, *lbprocs;
     int    modified_num_stored_um_fields;
     long   **tmp, **lookup;
     size_t n;
     FILE   *fid;
     struct tm t1, t2;
     char   varname[60];
     double *tdiff;
     float tol;

/**
 ** Attempt to open the file 
 **---------------------------------------------------------------------------*/
     fid = fopen( filename, "r" );
     if ( fid==NULL ) { return 0; }

/**
 ** Initialize the word size of the input UM fields file 
 **---------------------------------------------------------------------------*/
     wordsize = get_file_endianness_wordsize( fid );
     if ( wordsize==-1 ) {
        printf( "ERROR: incorrect wordsize detected!\n" );
        exit(1);  
     }

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
     printf( "# of 2D slices found: %d\n", cnt );

/**
 ** Gather all the valid lookup entries into a condensed array (called LOOKUP)) 
 ** Only take the first 45 entires per UM field as there is a variable change
 ** from long to double at that point.
 **---------------------------------------------------------------------------*/
     lookup = (long **) malloc( cnt*sizeof(long *) );

     cnt = 0;
     for ( nrec=0; nrec<header[151]; nrec++ ) {
         if ( tmp[nrec][28]!=-99 ) {
            lookup[cnt] = (long *) malloc( 45*sizeof(long) );
            for ( n=0; n<45; n++ ) {
                lookup[cnt][n] = tmp[nrec][n]; 
            } 
            cnt++;
         }
     }
  //   num_um_vars = cnt+1;

     for ( nrec=0; nrec<header[151]; nrec++ )
           free( tmp[nrec] );
     free( tmp );

/**
 ** If the user has NOT requested specific stash codes, determine # of unique 
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

     /*** Allocate UM variable data structure ***/ 
        stored_um_vars = (new_um_variable *) malloc( num_stored_um_fields*sizeof(new_um_variable) );
        for ( j=0; j<num_stored_um_fields; j++ ) {
            stored_um_vars[j].stash_code = temp[j];
            stored_um_vars[j].nz         = 0;
            stored_um_vars[j].nt         = 0;
        }

        free( temp );
     }

/**
 ** Asign the 2D data slices found in the input UM data file to the appropriate
 ** UM data structure. 
 **---------------------------------------------------------------------------*/
     modified_num_stored_um_fields = num_stored_um_fields;
     for ( j=0; j<num_stored_um_fields; j++ ) {

     /** How many 2D data slices belong to this variable in total? **/
         for ( i=0; i<num_um_vars; i++ ) 
             if ( stored_um_vars[j].stash_code==(unsigned short )lookup[i][41] ) { stored_um_vars[j].nz++; }
 
     /** Which 2D data slices belong to this variable? **/
         temp_id = (unsigned short *) malloc( (int )stored_um_vars[j].nz*sizeof(unsigned short) );
         cnt = 0;
         for ( i=0; i<num_um_vars; i++ ) {
             if ( stored_um_vars[j].stash_code==(unsigned short )lookup[i][41] ) { 
                temp_id[cnt] = (unsigned short int ) i;
                cnt++;
            }
         }

     /** What is the offset from the reference forecast time for each 2D data slice  **/
     /** belonging to this variable?                                                 **/
         tdiff = (double *) malloc( (int )stored_um_vars[j].nz*sizeof(double) );
         for ( i=0; i<stored_um_vars[j].nz; i++ ) {
             t1.tm_year = (int ) lookup[temp_id[i]][0];
             t1.tm_mon  = (int ) lookup[temp_id[i]][1];
             t1.tm_mday = (int ) lookup[temp_id[i]][2];
             t1.tm_hour = (int ) lookup[temp_id[i]][3];
             t1.tm_min  = (int ) lookup[temp_id[i]][4];
             t1.tm_sec  = (int ) lookup[temp_id[i]][5];
             t2.tm_year = (int ) lookup[temp_id[i]][6];
             t2.tm_mon  = (int ) lookup[temp_id[i]][7];
             t2.tm_mday = (int ) lookup[temp_id[i]][8];
             t2.tm_hour = (int ) lookup[temp_id[i]][9];
             t2.tm_min  = (int ) lookup[temp_id[i]][10];
             t2.tm_sec  = (int ) lookup[temp_id[i]][11];
             tdiff[i] = difftime ( mktime(&t1), mktime(&t2) )/3600.0; 

         /*** Accumulated fields record only a relative time offset to a reference ***/
         /*** value stored elsewhere in the lookup array.                          ***/
             if ( tdiff[i]<0.0 ) { tdiff[i] = (float ) lookup[temp_id[i]][13] - tdiff[i]; }
         }

         for ( i=0; i<stored_um_vars[j].nz; i++ ) {
             tol = (tdiff[i]-999999.0)*(tdiff[i]-999999.0);
             if ( tol>0.0001 ) {
                for ( k=i+1; k<stored_um_vars[j].nz; k++ ) {
                    tol = (tdiff[k]-tdiff[i])*(tdiff[k]-tdiff[i]);
                    if ( tol<0.0001 ) { tdiff[k]=999999.0; }
                }
             } 
         }

     /** How many unique time offset values belong to this variable? **/
         for ( i=0; i<stored_um_vars[j].nz; i++ ) { 
             tol = (tdiff[i]-999999.0)*(tdiff[i]-999999.0);
             if ( tol>0.0001 ) { stored_um_vars[j].nt++; }
         }

     /** Store the unique time offset values belonging to this variable **/
         stored_um_vars[j].times = (float *) malloc( (int )stored_um_vars[j].nt*sizeof(float) );

         k = 0;
         for ( i=0; i<stored_um_vars[j].nz; i++ ) { 
             tol = (tdiff[i]-999999.0)*(tdiff[i]-999999.0);
             if ( tol>0.0001 ) { stored_um_vars[j].times[k]=tdiff[i]; k++; }
         }

     /** Does NT evenly divide into NZ? **/
         ntmp = stored_um_vars[j].nz;
         stored_um_vars[j].nz = stored_um_vars[j].nz/stored_um_vars[j].nt;

     /** If so, go ahead and allocate memory for the 2D data slices belonging to **/
     /** this UM variable.  Assign the proper IDs to the data slices.            **/
         if ( ntmp==stored_um_vars[j].nz*stored_um_vars[j].nt ) { 
            
         /*** Allocate memory for the 2D data slice array ***/
            stored_um_vars[j].slices = (um_dataslice **) malloc( stored_um_vars[j].nt*sizeof(um_dataslice *) );
            for ( i=0; i<stored_um_vars[j].nt; i++ )
                stored_um_vars[j].slices[i] = (um_dataslice *) malloc( stored_um_vars[j].nz*sizeof(um_dataslice) );

         /*** Assign ID of each 2D data slice that belongs to this UM variable ***/
            cnt = 0;
            for ( i=0; i<stored_um_vars[j].nt; i++ ) {
            for ( k=0; k<stored_um_vars[j].nz; k++ ) {
                stored_um_vars[j].slices[i][k].id = (unsigned short ) temp_id[cnt];
                cnt++;
            }
            }

            stored_um_vars[j].lbproc = (unsigned short int ) lookup[temp_id[0]][24];  

         }
 
     /** If not, this means that multiple processed fields have been amalgamated into the 1 **/
     /** UM variable.  They need to be separated into individual UM_VAR struct elements.    **/
         else {

         /*** How many unique LBPROC values exist for this UM variable? ***/
            lbprocs = (int *) malloc( ntmp*sizeof(int) );
            for ( k=0; k<ntmp; k++ )
                lbprocs[k] = (int ) lookup[temp_id[k]][24]; 

            for ( k=0; k<ntmp; k++ ) {
                if ( lbprocs[k]!=-1 ) {
                   for ( i=k+1; i<ntmp; i++ ) 
                      if ( lbprocs[k]==lbprocs[i] ) { lbprocs[i] = -1; } 
                }
            }

            num_lbproc = 0;
            for ( k=0; k<ntmp; k++ ) 
                if ( lbprocs[k]!=-1 ) { num_lbproc++; } 

            modified_num_stored_um_fields += (num_lbproc-1);
            stored_um_vars = (new_um_variable *) realloc( stored_um_vars, 
                                                         ((int ) modified_num_stored_um_fields)*sizeof(new_um_variable) ); 

         /*** Alter NT and NZ values to accommodate a multiply-defined UM variable ***/
            stored_um_vars[j].nt -= (num_lbproc-1);
            stored_um_vars[j].nz = ntmp/(num_lbproc*stored_um_vars[j].nt);
            for ( i=0; i<num_lbproc; i++ ) {
                cnt = modified_num_stored_um_fields - i;
                stored_um_vars[cnt].stash_code = stored_um_vars[j].stash_code;
                stored_um_vars[cnt].nt = stored_um_vars[j].nt;
                stored_um_vars[cnt].nz = stored_um_vars[j].nz;
                stored_um_vars[cnt].times = (float *) malloc( (int )stored_um_vars[j].nt*sizeof(float) );
                for ( k=0; k<stored_um_vars[j].nt; k++ )
                    stored_um_vars[cnt].times[k] = stored_um_vars[j].times[k+i];
            }

         /*** Allocate memory for the 2D data slice array ***/
            stored_um_vars[j].slices = (um_dataslice **) malloc( stored_um_vars[j].nt*sizeof(um_dataslice *) );
            for ( i=0; i<stored_um_vars[j].nt; i++ ) 
                stored_um_vars[j].slices[i] = (um_dataslice *) malloc( stored_um_vars[j].nz*sizeof(um_dataslice) );

            for ( k=0; k<num_lbproc; k++ ) {
                cnt = modified_num_stored_um_fields - k;
                stored_um_vars[cnt].slices = (um_dataslice **) malloc( stored_um_vars[j].nt*sizeof(um_dataslice *) );
                for ( i=0; i<stored_um_vars[j].nt; i++ ) 
                    stored_um_vars[cnt].slices[i] = (um_dataslice *) malloc( stored_um_vars[j].nz*sizeof(um_dataslice) );
            } 

         /*** Assign ID of each 2D data slice that belongs to this UM variable ***/
            cnt = 0;
            for ( i=0; i<stored_um_vars[j].nt; i++ ) {
            for ( k=0; k<stored_um_vars[j].nz; k++ ) {
                stored_um_vars[j].slices[i][k].id = (unsigned short ) temp_id[cnt];
                cnt++;
                ind = modified_num_stored_um_fields - num_lbproc + 1;
                for ( kk=0; kk<num_lbproc-1; kk++ ) {
                    stored_um_vars[ind].slices[i][k].id = (unsigned short ) temp_id[cnt];
                    cnt++;
                    ind++;
                }
            }
            }
            free( lbprocs );

          /*** Set the LBPROC value for the entire UM variable ***/
            stored_um_vars[j].lbproc = (unsigned short ) lookup[0][24];
            ind = modified_num_stored_um_fields - num_lbproc + 1;
            for ( k=0; k<num_lbproc-1; k++ ) {
                i = stored_um_vars[ind].slices[0][0].id;
                stored_um_vars[ind].lbproc = (unsigned short int ) lookup[i][24];
                ind++;
            }

         }

         free( temp_id );
         free( tdiff );
     }

     num_stored_um_fields = modified_num_stored_um_fields;

/**
 ** Add some necessary attributes to each UM variables (and its associated 2D
 ** data slices) 
 **---------------------------------------------------------------------------*/
     for ( j=0; j<num_stored_um_fields; j++ ) {
     for ( i=0; i<num_um_vars; i++ ) {

         if ( (unsigned short ) lookup[i][41]==stored_um_vars[j].stash_code ) {
            stored_um_vars[j].coordinates = (unsigned short ) lookup[i][15];
            stored_um_vars[j].ny          = (unsigned short ) lookup[i][17];
            stored_um_vars[j].nx          = (unsigned short ) lookup[i][18];
            stored_um_vars[j].lbvc        = (unsigned short ) lookup[i][25];
            if ( rflag==0 ) {
               if ( lookup[i][38]==1 ) { stored_um_vars[j].vartype = NC_DOUBLE; }
               else                    { stored_um_vars[j].vartype = NC_LONG; }
            } else {
               if ( lookup[i][38]==1 ) { stored_um_vars[j].vartype = NC_FLOAT; }
               else                    { stored_um_vars[j].vartype = NC_INT; }
            }

            for ( kk=0; kk<stored_um_vars[j].nt; kk++ ) {
            for ( k=0;   k<stored_um_vars[j].nz;  k++ ) {
                if ( stored_um_vars[j].slices[kk][k].id == (unsigned short int ) i ) {
                   stored_um_vars[j].slices[kk][k].size      = lookup[i][14];
                   stored_um_vars[j].slices[kk][k].location  = lookup[i][28];
                   stored_um_vars[j].slices[kk][k].reclength = lookup[i][29];
                   stored_um_vars[j].slices[kk][k].level     = (unsigned short ) lookup[i][32];
                   stored_um_vars[j].slices[kk][k].lbproc    = lookup[i][24];
                   stored_um_vars[j].slices[kk][k].lbpack    = (unsigned short ) lookup[i][20];
                   stored_um_vars[j].slices[kk][k].mdi       = (double ) lookup[i][62];
                   stored_um_vars[j].slices[kk][k].datatime.tm_year = (int ) lookup[i][0];
                   stored_um_vars[j].slices[kk][k].datatime.tm_mon  = (int ) lookup[i][1];
                   stored_um_vars[j].slices[kk][k].datatime.tm_mday = (int ) lookup[i][2];
                   stored_um_vars[j].slices[kk][k].datatime.tm_hour = (int ) lookup[i][3];
                   stored_um_vars[j].slices[kk][k].datatime.tm_min  = (int ) lookup[i][4];
                   stored_um_vars[j].slices[kk][k].datatime.tm_sec  = (int ) lookup[i][5];
                }
            }
            }
         }
     }
     }

/**
 ** Locate the index of the XML file where each unique UM variable is given
 ** its CF-compliant description.   
 **---------------------------------------------------------------------------*/
     for ( i=0; i<num_stored_um_fields; i++ ) {
     for ( j=0; j<num_xml_vars; j++ ) {
         if ( stored_um_vars[i].stash_code == 1000*um_vars[j].section + um_vars[j].code ) {
            stored_um_vars[i].xml_index   = j;
            strcpy( stored_um_vars[i].name, um_vars[j].varname );
            stored_um_vars[i].grid_type   = (unsigned short ) um_vars[j].umgrid;
            stored_um_vars[i].accum       = (unsigned short ) um_vars[j].accum;
            stored_um_vars[i].level_type  = (unsigned short ) um_vars[j].level_type;
            break;
         }
     }
     }

/**
 ** Add an appropriate prefix to the UM variable's name if post-processing was
 ** performed on the field. 
 **---------------------------------------------------------------------------*/
     for ( i=0; i<num_stored_um_fields; i++ ) {
         kk = stored_um_vars[i].lbproc;
         if ( kk!=0 ) { 
            if ( kk==64 )   { strcpy( varname, "mean_" ); }
            if ( kk==4096 ) { strcpy( varname, "min_"  ); }
            if ( kk==8192 ) { strcpy( varname, "max_"  ); }
            if ( kk==128 ) {
               if (stored_um_vars[i].accum==0 ) { strcpy( varname, "mean_" ); }
               else                             { strcpy( varname, "sum_" ); }
            }
            strcat( varname, stored_um_vars[i].name );
            strcpy( stored_um_vars[i].name, varname ); 
         }
     } 

     for ( i=0; i<num_stored_um_fields; i++ ) {
/**
 ** For UM variables that are a temporal mean or accummulation, compute the 
 ** duration of each mean/accummulation. 
 **---------------------------------------------------------------------------*/
         if ( (stored_um_vars[i].lbproc==32)||(stored_um_vars[i].lbproc==128)||
              (stored_um_vars[i].lbproc==4096)||(stored_um_vars[i].lbproc==8192) ) {

            stored_um_vars[i].time_bnds = (float **) malloc( stored_um_vars[i].nt*sizeof(float *) );
            for ( k=0; k<stored_um_vars[i].nt; k++ )
                stored_um_vars[i].time_bnds[k] = (float *) calloc( 2,sizeof(float) );

         }

/**
 ** For UM variables that are a spatial mean or accummulation, set the end 
 ** points for the accummulation operation.
 **---------------------------------------------------------------------------*/
         if ( (stored_um_vars[i].lbproc==8)||(stored_um_vars[i].lbproc==64) ) {
            stored_um_vars[i].space_bnds[0] = (float ) real_constants[3]; 
            stored_um_vars[i].space_bnds[1] = (float ) (real_constants[3] + real_constants[0]*stored_um_vars[i].nx);
         }  

     }


/*==============================================================================
 * START OF SANITY CHECK
 *============================================================================== * 
      for ( i=0; i<num_stored_um_fields; i++ ) {
         printf( "%s %hu [%hu, %hu %hu, %hu]\n", stored_um_vars[i].name, stored_um_vars[i].stash_code,
                                               stored_um_vars[i].nx, stored_um_vars[i].ny, stored_um_vars[i].nz,
                                               stored_um_vars[i].nt );
        printf( "TIMES: " );
        for ( j=0; j<stored_um_vars[i].nt; j++ ) { printf( "%f ", stored_um_vars[i].times[j] ); }
        printf( "\n" );
        printf( "SLICES: " );
        for ( j=0; j<stored_um_vars[i].nt; j++ ) { 
        for ( k=0; k<stored_um_vars[i].nz; k++ ) {
            printf( "%hu ", stored_um_vars[i].slices[j][k].id ); 
        }
        }
        printf( "\n" );
     }
     exit(1);
 *==============================================================================
 * END OF SANITY CHECK 
 *==============================================================================*/ 

/**
 ** Free memory used by the LOOKUP 2D array
 **---------------------------------------------------------------------------*/
     for ( nrec=0; nrec<num_um_vars; nrec++ )
         free( lookup[nrec] );
     free( lookup );

     fclose( fid );
     return 1; 
}
