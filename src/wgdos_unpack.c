#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "field_def.h"

/* Function Prototypes */

void endian_swap_4bytes( void *ptr, int nchunk );
int convert_float_ibm_to_ieee32(int ibm[], int ieee[], int* n);

#define   expon 0x7F000000
#define   sign  0x80000000
#define   tiss  0x00FFFFFF
#define   etis  0x007FFFFF
#define   nrm   0x00F00000

/***
 *** CONVERT IBM FLOAT TO IEEE FLOAT 
 ***
 *** Subroutine that converts a 32-bit IBM float into a 32-bit IEEE float. 
 ***
 ***  INPUT: ibm -> integer containing all bits in the IBM float
 ***
 ***   Mark Cheeseman, NIWA
 ***   April 2, 2014
 ***/
/*
double convert_ibm_float_to_ieee_float( int32_t ibm ) {

      int32_t ieee, ibs, ibe, ibt, it, k;
      union { int32_t i; float r; } u;
*/
    /* Grab sign bit */
//      ibs = ibm & sign;

    /* Determine the appropriate exponent & mantissa bits */
/*      ibt = ibm & tiss;
      if ( ibt==0 ) { ibe = 0; }
      else {

         ibe = ibm & expon;
         if ( (ibe != 0) && (ibt & nrm) == 0 ) {
            u.i = ibm;
            u.r = u.r + 0e0;
            ibe = u.i & expon;
            ibt = u.i & tiss;
         }
     
         it = ibt << 8;
         for ( k=0; (k < 5) && (it >= 0); k++ ) { it = it << 1; }
         
         if ( k < 4 ) {
            ibt = (it >> 8) & etis;
       //     ibe = (ibe >> 22) - 256 + 127 - k - 1;
            ibe = (ibe >> 22) - 130 - k;
            if (ibe < 0) { ibe = ibt = 0; }
        
            if (ibe >= 255) { 
 //              printf( "ERROR: IEEE exponent cannot exceed 8 bits\n" );
 //              exit(1);
               ibe = 255; ibt = 0; 
            }
            ibe = ibe << 23;
         }
      }
*/
    /* Construct the newly formed IEEE float */
//      ieee = ibs | ibe | ibt;
  
//      return (double ) ieee; 
//}


/***
 *** APPLY BITMAP
 ***
 *** Subroutine that applies the appropriate bitmap masking to the unpacked data
 *** array.
 ***
 ***  INPUT: buf -> character array containing the packed data
 ***         val -> bitmap value being applied (MDI, 0 or BASE_VAL)
 ***          nx -> X dimension of unpacked array
 ***
 ***  INPUT/OUTPUT: array       -> ptr to the unpacked data array
 ***                bit_index   -> current bit in the character buffer being evaluated
 ***                total_bmap  -> ptr to bitmap that masks all zero, MDI and min 
 ***                               values
 ***
 ***   Mark Cheeseman, NIWA
 ***   April 1, 2014
 ***/

void apply_bitmap( char *buf, int bit_index, double val, unsigned short nx, 
                   double *array, unsigned short *total_bmap ) {
  
     int            bit, start_byte; 
     unsigned short n, check, start_bit;

  /* Determine the starting byte in the char buffer (and start bit in that word) */
     start_byte = bit_index / 8; 
     start_bit  = bit_index % 8;

     bit = 0;
     while ( bit<nx ) {

           for ( n=start_bit; n<8; n++ ) {
               check = (buf[start_byte] >> n) & 1;
               if ( check==1 ) { 
                  array[bit] = val; 
                  total_bmap[bit] = 1;
               }
               bit++;
               if ( bit>=nx ) { break; }
           }
           start_bit = 0;
           start_byte++; 
     }
     bit_index += nx;

}


/***
 *** EXTRACT DATAVAL 
 ***
 *** Subroutine that extracts the packed data values. 
 ***
 ***  INPUT:  cbuf        -> char buffer holding the raw packed data 
 ***          bit_index   -> current starting bit in the char buffer 
 ***          array_index -> starting index of the unpacked data array 
 ***          nbit        -> # of bits each for each packed data value
 ***
 ***  OUTPUT: unpacked -> ptr to the 2D array of the unpacked data values
 ***          
 ***
 ***   Mark Cheeseman, NIWA
 ***   April 1, 2014
 ***/

double extract_dataval( char *cbuf, int bit_index, int16_t nbit ) { 

     int            ierr, baselen;
     int32_t        tmp;
     unsigned short start_bit, start_byte;
     char           *ptr = NULL;
     float          ftmp;

  /* Determine the starting byte in the char buffer (and start bit in that word) */
     start_byte = bit_index / 8;
     start_bit  = bit_index % 8; 

  /* Extract 4 bytes from the char buffer & convert them into a 32-bit int */
     ptr = cbuf;
     ptr += start_byte;
     memcpy( &tmp, &cbuf[start_byte], 4 );

  /* Bit-shift and mask out the unwanted bits */
     tmp = tmp>>start_bit;
     tmp = tmp & (int32_t )(pow(2,nbit)-1); 

  /* Construct the final 64-bit data value */
     baselen = 1;
     ierr = convert_float_ibm_to_ieee32( &tmp, (int*)&ftmp, &baselen ); 
     if ( ierr==-1 ) { printf( "ERROR: conversion failed" ); exit(1); }
     return (double ) ftmp;

}


/***
 *** WGDOS UNPACK 
 ***
 *** Subroutine that coordinates the unpacking of the 2D slice of an UM variable
 *** packed using the WGDOS packing scheme.
 ***
 ***  INPUT:  fid      -> file handle to the packed UM fieldsfile
 ***          nx       -> # of points in a data row
 ***          ny       -> # of rows in the 2D slice being processed
 ***          mdi      -> value used to denote a missing data point
 ***
 ***  OUTPUT: unpacked -> 2D array of the unpacked data values
 ***
 ***   Mark Cheeseman, NIWA
 ***   April 1, 2014
 ***/

void wgdos_unpack( FILE *fid, unsigned short nx, unsigned short ny, double *unpacked, double mdi ) {

     int      i, j, bit_index, baselen;
     int32_t  ibuf[3], ibuf2[2], k;
     int16_t  N;
     uint16_t nbit;
     char    *cbuf;
     float    ftmp;
     double   val, base_val, prec, *row;
     unsigned short *total_bitmap, min_bitmap_present, mdi_bitmap_present, zero_bitmap_present;

  /* Read & decode field header */
     fread( ibuf, 4, 3, fid );
     endian_swap_4bytes( ibuf, 3 );
     prec = pow( 2, ibuf[1] );

  /* Allocate memory for a unpacked row of data points */
     row = (double *) malloc( nx*sizeof(double) );

     for ( j=0; j<ny; j++ ) {

         fread( &ibuf2, 4, 2, fid );
         endian_swap_4bytes( &ibuf2, 2 );

     /* Decode base value */
         baselen = 1;
         k = ibuf2[0];
         i = convert_float_ibm_to_ieee32( &k, (int*)&ftmp, &baselen ); 
         if ( i==-1 ) { printf("ERROR: conversion failed\n"); exit(1); }
         base_val = (double )ftmp;

         for ( i=0; i<nx; i++ ) { row[i] = base_val; }

     /* Determine # of packed data points */
         N = ibuf2[1] & 65535;
         if ( N==0 ) { continue; }

     /* Determine which bitmaps are present in current row */
         mdi_bitmap_present  = (ibuf2[1]>>16) & 128;
         min_bitmap_present  = (ibuf2[1]>>16) & 64;
         zero_bitmap_present = (ibuf2[1]>>16) & 32;
         
     /* Determine how many bits used per packed data value in current row */
         nbit = (ibuf2[1]>>16) & 31;

     /* Read in character buffer that holds the packed bitmap(s) and data */
         cbuf = (char *) malloc( sizeof(int32_t)*N ); 
         fread( cbuf, 4, N, fid );

     /* Read in active bitmaps */
         total_bitmap = (unsigned short *) calloc( nx, sizeof(unsigned short) );

         bit_index = 0;
         if ( mdi_bitmap_present != 0 ) {
      //      printf( "MDI bitmap present\n" );
            apply_bitmap( cbuf, bit_index, mdi, nx, row, total_bitmap );
            bit_index += nx; 
         }
         if ( min_bitmap_present != 0 ) {
      //      printf( "MIN bitmap present\n" );
            apply_bitmap( cbuf, bit_index, base_val, nx, row, total_bitmap );
            bit_index += nx; 
         }
         if ( zero_bitmap_present != 0 ) {
     //       printf( "ZERO bitmap present\n" );
            apply_bitmap( cbuf, bit_index, 0.0, nx, row, total_bitmap );
            bit_index += nx; 
         }

        /* Extract the packed data values */
            for ( i=0; i<nx; i++ ) { 
                if ( total_bitmap[i]==0 ) {
                   val = extract_dataval( cbuf, bit_index, nbit ); 
                   row[i] = base_val + prec * val;
                   bit_index += 32;
                }
         //       printf( "%f ", row[i] );
            } 
           // printf( "\n" );

     /* Copy row values into main data array */
         for ( i=0; i<nx; i++ ) { unpacked[i + j*nx] = row[i]; } 

         free( total_bitmap ); 
         free( cbuf );
     }
     free( row );

}
