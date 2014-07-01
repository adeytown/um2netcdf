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


#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "field_def.h"

/** Function prototypes **/

void endian_swap_4bytes( void *ptr, int nchunk );
int ibm2ieee(uint32_t ibm[], float ieee[], int n);


typedef union {
   uint32_t w;
   unsigned char b[4];
} W32;

typedef union {
   uint16_t w;
   unsigned char b[2];
} W16;

uint32_t byteswap32(unsigned char bytes[4])
{
   W32 in;
   memcpy(&in.b, bytes, 4);
#ifndef IBM_POWER
   unsigned char tmp = in.b[0];
   in.b[0] = in.b[3];
   in.b[3] = tmp;
   tmp = in.b[1];
   in.b[1] = in.b[2];
   in.b[2] = tmp;
#endif
   return in.w;
}

uint16_t byteswap16(unsigned char bytes[2])
{
   W16 in;
   memcpy(&in.b, bytes, 2);

#ifndef IBM_POWER
   unsigned char tmp = in.b[0];
   in.b[0] = in.b[1];
   in.b[1] = tmp;
#endif
   return in.w;
}

static inline
uint32_t getbits(unsigned char* bp, int pos, int nbits)
{
   assert(pos < 8 && nbits <= 32);

   uint32_t res = 0;
   int more = nbits;
   while (1) {
      int bits = ((8 - pos) < more) ? 8 - pos : more;
      uint8_t mask = ((1 << bits) - 1) << (8 - pos - bits);
      uint8_t val = (*bp & mask) >> (8 - pos - bits);
      more -= bits;
      res |= val << more;
      if (!more) break;
      bp++; pos = 0;
   }

   return res;
}


/***
 *** READ_BITMASKS
 ***
 *** Function that reads in the bitmap indicate where missing data points are
 *** located in the datafield. 
 ***
 ***  INPUT: bp -> pointer to the input data stream
 ***         nx -> # of columns present in the input 2D data slice
 ***
 ***  INPUT/OUTPUT: flag_array -> integer array that denotes if a MDI value
 ***                              should be written at the iTH position.
 ***
 ***  Function will return the number of 32 bit words read for the bitmap
 ***  extraction and decoding.  
 ***
 ***   Mark Cheeseman, NIWA
 ***   June 26, 2014
 ***/ 

int read_bitmasks( unsigned char* bp, uint16_t ncols, uint32_t *mdi_array, uint32_t *zero_array, 
                   bool a, bool c ) {

     int     n, start_bit, num_words_read; 
     uint8_t mask, val;

     start_bit = 0;
     num_words_read = 0;

 /*
  * Extract & read MDI bitmask (if present)
  *--------------------------------------------------------------------*/  
     if ( c ) {
        for ( n=0; n<ncols; n++ ) {
            mask = 1 << start_bit;
            val = *bp & mask;
            if ( val==0 ) { mdi_array[n] = 1; }
            start_bit++;
            if ( start_bit>7 ) { bp++; start_bit=0; num_words_read++; }
        }
     }
     if ( start_bit>0 ) { num_words_read++; }

 /*
  * Extract & read ZEROS bitmask (if present)
  *--------------------------------------------------------------------*/  
     if ( a ) {
        for ( n=0; n<ncols; n++ ) {
            mask = 1 << start_bit;
            val = *bp & mask;
            if ( val>0 ) { zero_array[n] = 1; }
            start_bit++;
            if ( start_bit>7 ) { bp++; start_bit=0; num_words_read++; }
        }
     }

 /* Round up to the nearest word */
     if ( start_bit>0 ) { num_words_read++; }

     return num_words_read;
}


/***
 *** WGDOS UNPACK 
 ***
 *** Subroutine that unpacks a 2D data slice that has undergone WGDOS packing &
 *** compression
 ***
 *** INPUT:   fh -> file handle to the input UM fields file
 ***         ind -> index to current UM variable being processed
 ***
 *** OUTPUT: unpacked_data -> pointer to the array of values for the unpacked 2D data 
 ***                          slice     
 ***
 ***   Mark Cheeseman, NIWA
 ***   June 26, 2014
 ***/

void wgdos_unpack( FILE *fh, double *unpacked_data, double mdi ) {

     int      i, j, nbits, pos, new_pos;
     uint16_t ncol, nrow, n;
     int32_t  prec;
     uint32_t ibm, val, *zero_mask, *mdi_mask;
   
     float  base; 
     double dval=0.0, scaling_factor; 
     bool   a, c;

     char          cba_nbit; 
     unsigned char buf[5258], *bp;
       
 
  /*
   * Read & decode field header
   *-------------------------------------------------------------------*/   
     bp = buf;
     fread(buf, 4, 5, fh);  
     bp += 4;
     prec = byteswap32(bp);
     scaling_factor = pow( 2.0, (double )prec );
     bp += 4;
     ncol = byteswap16(bp);
     bp += 2;
     nrow = byteswap16(bp);
     bp += 2;
#ifdef IBM_POWER
     n = ncol;
     ncol = nrow;
     nrow = n;
#endif

/*     printf( "scaling factor: %f\n", scaling_factor );
     printf( "# of rows: %d\n", nrow );
     printf( "# of columns: %d\n", ncol ); */

  /*
   * Allocate memory to hold zero and MDI masks for 1 unpacked row 
   *-------------------------------------------------------------------*/   
     mdi_mask  = (uint32_t *) calloc( ncol,sizeof(uint32_t) );
     zero_mask = (uint32_t *) calloc( ncol,sizeof(uint32_t) );

     for ( j=0; j<nrow; ++j) {

  /*
   * Read & decode a row's header
   *   BASE -> minimum value of all points in unpacked row
   *   NBITS-> # of bits used for each packed data point in this row
   *   N    -> # of 32 bit words used to hold all packed data points
   *           and bitmaps for this row
   *   C    -> boolean that denotes if zero bitmap is present
   *   A    -> boolean that denotes if missing value bitmap is present
   *-------------------------------------------------------------------*/   
         ibm = byteswap32(bp);
         ibm2ieee( &ibm, &base, 1 );
         bp += 4;

         cba_nbit = *(bp+1);          // want low byte only (big endian)
         c = cba_nbit & 0x80;
//         b = cba_nbit & 0x40;
         a = cba_nbit & 0x20;
         nbits = cba_nbit & 0x1F;
         bp += 2;
         n = byteswap16(bp);

  /*
   * Read in contents of packed row 
   *-------------------------------------------------------------------*/   
         if (fread(buf, 4, n+2, fh) < n+2) {      // row data + next row hdr
            puts("end of file reached...");
            break;
         }
      /*   printf( "base value: %f\n", base );
         printf( "nbit:       %d\n", nbits );
         if ( a ) { printf( "zeros bitmap present\n" ); }
         if ( b ) { printf( "min value bitmap present\n" ); }
         if ( c ) { printf( "mdi bitmap present\n" ); }
         exit(1);*/

     /*
      * Read & decode the masking bitmaps 
      *-------------------------------------------------------------------*/   
         bp = buf;
         if ( a | c ) { pos = read_bitmasks( bp, ncol, mdi_mask, zero_mask, a, c ); bp += pos; }   
//         if ( b )     { printf( "MINIUM VALUE mask needed\n" ); exit(1);         }

     /*
      * Decode the packed data points 
      *-------------------------------------------------------------------*/   
         pos = 0;
         for ( i=0; i<ncol; ++i) {

             if ( zero_mask[i]==1 ) { dval = 0.0; }
             if ( mdi_mask[i]==1 )  { dval = 1.0e-30; }
             if ( mdi_mask[i]+zero_mask[i]==0 ) {
                val = getbits( bp, pos, nbits );
                dval = (double ) (base + scaling_factor*val);
                new_pos = pos + nbits;
                bp += (new_pos) / 8;
                pos = (new_pos) % 8;
             }

             unpacked_data[i+j*ncol] = dval;

         }
         bp = buf + n*4;

     }
     free( zero_mask );
     free( mdi_mask );

     return;
}

