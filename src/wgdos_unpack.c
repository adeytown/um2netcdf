#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "field_def.h"

/* Function Prototypes */

void endian_swap_4bytes( void *ptr, int nchunk );

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

double convert_ibm_float_to_ieee_float( int32_t ibm ) {

      int32_t ieee, ibs, ibe, ibt, it, k;
      union { int32_t i; float r; } u;

    /* Grab sign bit */
      ibs = ibm & sign;

    /* Determine the appropriate exponent & mantissa bits */
      ibt = ibm & tiss;
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

    /* Construct the newly formed IEEE float */
      ieee = ibs | ibe | ibt;
  
      return (double ) ieee; 
}


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
 ***                array_index -> current starting element in unpacked data array 
 ***                total_bmap  -> ptr to bitmap that masks all zero, MDI and min 
 ***                               values
 ***
 ***   Mark Cheeseman, NIWA
 ***   April 1, 2014
 ***/

void apply_bitmap( char *buf, int bit_index, double val, unsigned short nx, 
                   double *array, int array_index, unsigned short *total_bmap ) {
  
     int            bit, start_byte; 
     unsigned short n, check, start_bit;

  /* Determine the starting byte in the char buffer (and start bit in that word) */
     start_byte = bit_index / 8; 
     start_bit  = bit_index % 8;

  /* Check the first byte separately as we may not be starting on a byte boundary */
     for ( n=start_bit; n<8; n++ ) {
         check = (buf[start_byte] >> n) & 1;
         if ( check==1 ) { array[array_index] = val; array_index++; }
         bit_index++;
     }
     start_byte++;

     bit = 0;
     while ( bit<nx ) {
           for ( n=0; n<8; n++ ) {
               check = (buf[start_byte] >> n) & 1;
               if ( check==1 ) { 
                  array[array_index] = val; 
                  array_index++; 
                  total_bmap[bit] = 1;
               }
               bit++;
               if ( bit>=nx ) { break; }
           }
           start_byte++; 
     }
     bit_index += bit;

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

double extract_dataval( char *cbuf, int bit_index, int16_t nbit, double base, double acc ) { 

     int            start_byte;
     int32_t        tmp;
     unsigned short start_bit, num_bytes;
     char           *ptr = NULL;
     double         val, dtmp;

  /* Determine # of bytes that NBIT bits occupy */
     num_bytes = nbit / 8;
     if ( nbit%8 > 0 ) { num_bytes++; }

  /* Determine the starting byte in the char buffer (and start bit in that word) */
     start_byte = bit_index / 8;
     start_bit  = bit_index % 8; 

  /* Extract 4 bytes from the char buffer & convert them into a 32-bit int */
     ptr = cbuf;
     ptr += start_byte;
     memcpy( &tmp, ptr, 4 );

  /* Bit-shift and mask out the unwanted bits */
     tmp = tmp>>start_bit;
     tmp = tmp & (int32_t )(pow(2,nbit)-1); 

  /* Construct the final 64-bit data value */
     dtmp = convert_ibm_float_to_ieee_float( tmp );
     val = base + acc * dtmp;

     return val;
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

     int     i, j, bitmap_flags[3], bit_index, array_index;
     int32_t ibuf[3], ibuf2[2];
     int16_t N, nbit;
     char    *cbuf;
     double  bitmap_val[3], base_val, prec;
     unsigned short *total_bitmap;

  /* Read & decode field header */
     fread( ibuf, 4, 3, fid );
     endian_swap_4bytes( ibuf, 3 );
     prec = pow( 2, ibuf[1] );

     bitmap_val[2] = mdi;
     bitmap_val[0] = 0.0;

     printf( "%d %f\n", ibuf[0], prec );

     for ( j=0; j<ny; j++ ) {

         fread( &ibuf2, 4, 2, fid );
         endian_swap_4bytes( &ibuf2, 2 );

     /* Decode base value */
         if ( ibuf2[0]==0 ) { base_val = 0.0; } 
         else { base_val = convert_ibm_float_to_ieee_float( ibuf2[0] ); }
         bitmap_val[1] = base_val;

     /* Determine which bitmaps are active in current row */
         for ( i=0; i<3; i++ ) { bitmap_flags[i] = (ibuf2[1]>>i) & 1; }

     /* Determine how many bits used per packed data value in current row */
         nbit = (ibuf2[1] & 248)>>3;

         N = (ibuf2[1] & 16776960)>>8;

         printf( "%f %d %d %d %d %d\n", base_val, bitmap_flags[0], bitmap_flags[1], bitmap_flags[2], nbit, N );
         if ( N>0 ) {
        /* Read in character buffer that holds the packed bitmap(s) and data */
            cbuf = (char *) malloc( sizeof(int32_t)*N ); 
            fread( cbuf, 4, N, fid );

        /* Read in active bitmaps */
             total_bitmap = (unsigned short *) calloc( nx, sizeof(unsigned short) );

             bit_index = 0;
             array_index = nx*j;
             for ( i=0; i<3; i++ ) 
                 if ( bitmap_flags[0]==1 ) { apply_bitmap( cbuf, bit_index, bitmap_val[i], 
                                                           nx, unpacked, array_index, total_bitmap ); } 

        /* Extract the packed data values */
            for ( i=0; i<nx; i++ ) { 
                if ( total_bitmap[i]==0 ) {
                   unpacked[array_index+i] = extract_dataval( cbuf, bit_index, nbit, 
                                                              base_val, prec );
                }
            }
   
            free( total_bitmap ); 
            free( cbuf );
         }
     }
         exit(1);

}
