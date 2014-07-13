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

#define   expon 0x7F000000
#define   sign  0x80000000
#define   tiss  0x00FFFFFF
#define   etis  0x007FFFFF
#define   nrm   0x00F00000

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

typedef union {
   uint32_t w;
   unsigned char b[4];
} W32;

typedef union {
   uint16_t w;
   unsigned char b[2];
} W16;

#ifdef IBM_POWER
uint32_t byteswap32(unsigned char bytes[4]) {

         W32 in;
         memcpy( &in.b, bytes, 4 );
         return in.w;
}

uint16_t byteswap16(unsigned char bytes[2]) {
         
         W16 in;
         memcpy( &in.b, bytes, 2 );
         return in.w;
}
#else
uint32_t byteswap32(unsigned char bytes[4])
{
   W32 in;
   memcpy(&in.b, bytes, 4);

   unsigned char tmp = in.b[0];
   in.b[0] = in.b[3];
   in.b[3] = tmp;
   tmp = in.b[1];
   in.b[1] = in.b[2];
   in.b[2] = tmp;

   return in.w;
}

uint16_t byteswap16(unsigned char bytes[2])
{
   W16 in;
   memcpy(&in.b, bytes, 2);

   unsigned char tmp = in.b[0];
   in.b[0] = in.b[1];
   in.b[1] = tmp;

   return in.w;
}
#endif


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

void readBitmap(unsigned char* bp, int start, int cols, bool reverse, float value, float data[], bool bmap[])
{
   int           i, pos;
   unsigned char byte;

   byte = *bp;
   if (reverse) byte = ~byte;
   byte <<= start;
   pos = start;
   for ( i=0; i < cols; ++i) {
      if (byte & 0x80) {
         data[i] = value;
         bmap[i] = true;
      }
      if (pos < 7) {
         byte <<= 1;
         pos++;
      } else {
         byte = *++bp;
         if (reverse) byte = ~byte;
         pos = 0;
      }
   }
}

float ibm2ieee2(uint32_t ibm)
{
   int32_t ibe, it;
   uint32_t ibs, ibt;
   int k;
   union { uint32_t i; float r; } u, res;

   ibs = ibm & sign;
   ibe = ibm & expon ;
   ibt = ibm & tiss;
   if (ibt == 0) {
      ibe = 0 ;
   }
   else {
      if ( (ibe != 0) && (ibt & nrm) == 0 ) {
         u.i = ibm;
         u.r = u.r + 0e0 ;
         ibe = u.i & expon ;
         ibt = u.i & tiss ;
      }
      /* mantissa */
      it = ibt << 8;
      for (k = 0; (k < 5) && (it >= 0); k++ ) {
         it = it << 1;
      }
      if ( k < 4 ) {
         ibt = (it >> 8) & etis;
         ibe = (ibe >> 22) - 256 + 127 - k - 1;
         if (ibe < 0) {
            ibe = ibt = 0;
         }
         if (ibe >= 255) {
            ibe = 255; ibt = 0;
         }
         ibe = ibe << 23;
      }
   }
   res.i = ibs | ibe | ibt;
   return res.r;
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

     int            i, j, nbits, pos, new_pos;
     uint16_t       cols, rows, n;
     uint32_t       len;
     int32_t        prec;
     float          scale, base, *unpacked_row;
     char           cba_nbit;
     unsigned char  hdr[20], *buf, *bp;
     bool           a, b, c, use_bmaps, *bmap;

  /*
   * Read & decode field header
   *-------------------------------------------------------------------*/   
     bp = hdr;
     fread( hdr, 4, 5, fh );

     len = byteswap32(bp);
     bp += 4;
     prec = byteswap32(bp);
     scale = powf( 2.0, (float ) prec );
     bp += 4;
     cols = byteswap16(bp);
     bp += 2;
     rows = byteswap16(bp);
     bp += 2;

/*     printf( "scaling factor: %d %f\n", prec, scale );
     printf( "# of rows: %d\n", rows );
     printf( "# of columns: %d\n\n", cols ); */

  /*
   * Allocate memory to hold unpacked data points in 1 row 
   *-------------------------------------------------------------------*/   
     unpacked_row = (float *) malloc( cols*sizeof(float) );

  /*
   * Allocate memory to hold bitmap mask for 1 unpacked row 
   *-------------------------------------------------------------------*/   
     bmap  = (bool *) malloc( cols*sizeof(bool) );

  /*
   * Allocate memory to hold entire packed row (+ up to 3 bitmaps) 
   *-------------------------------------------------------------------*/   
     pos = ((3*cols+7) / 8) + (cols*4) + 8;
     buf = (unsigned char *) malloc( pos*sizeof(unsigned char) );

     for ( j=0; j<rows; ++j) {

  /*
   * Decode a row's header
   *   BASE -> minimum value of all points in unpacked row
   *   NBITS-> # of bits used for each packed data point in this row
   *   N    -> # of 32 bit words used to hold all packed data points
   *           and bitmaps for this row
   *   C    -> boolean that denotes if zero bitmap is present
   *   B    -> boolean that denotes if minimum value bitmap is present
   *   A    -> boolean that denotes if missing value bitmap is present
   *-------------------------------------------------------------------*/   
         base = ibm2ieee2( byteswap32(bp) );
         bp += 4;

         cba_nbit = *(bp+1);
         c = cba_nbit & 0x80;
         b = cba_nbit & 0x40;
         a = cba_nbit & 0x20;
         use_bmaps = (a || b || c);
         nbits = cba_nbit & 0x1F;
         bp += 2;
         n = byteswap16(bp);

     /*    printf( "base value: %f\n", base );
         printf( "nbit:       %d\n", nbits );
         if ( a ) { printf( "zeros bitmap present\n" ); }
         if ( b ) { printf( "min value bitmap present\n" ); }
         if ( c ) { printf( "mdi bitmap present\n\n" ); }*/

  /*
   * Set all values in unpacked row to BASE initially.  If nbits==0,
   * we leave the data points to this uniform value 
   *-------------------------------------------------------------------*/   
         for ( i=0; i<cols; i++ ) { unpacked_row[i] = base; }

  /*
   * Read in contents of packed row (data points + bitmaps) 
   *-------------------------------------------------------------------*/   
         if ( fread(buf, 4, n+2, fh)<n+2 ) { break; }
         bp = buf;
         pos = 0;

  /*
   * If required, extract the bitmap masks 
   *-------------------------------------------------------------------*/   
         if ( use_bmaps ) {
            memset( bmap, 0, cols*sizeof(bool) );

         /** Read in MISSING DATA VALUE bitmap (if pesent) **/
            if (a) { 
               readBitmap( bp, pos, cols, false, mdi, unpacked_row, bmap );
               bp += cols / 8;
               pos = cols % 8;
            }
         /** Read in MINIMUM DATA VALUE bitmap (if pesent) **/
            if (b) { 
               readBitmap(bp, pos, cols, false, base, unpacked_row, bmap);
               bp += (pos + cols) / 8;
               pos = (pos + cols) % 8;
            }
         /** Read in ZERO DATA VALUE bitmap (if pesent) **/
            if (c) {  
               readBitmap(bp, pos, cols, true, 0.0, unpacked_row, bmap);
               bp += (pos + cols) / 8;
               pos = (pos + cols) % 8;
            }
         /** Make sure data pointer is aligned with the next 32-bit word **/
            if ( pos || (bp - buf) % 4 ) {
               bp += 4 - ((bp - buf) % 4);
               pos = 0;
            }
         }

  /*
   * Extract the packed data points 
   *-------------------------------------------------------------------*/   
         if ( nbits>0 ) {
            for ( i=0; i<cols; ++i) {
                if ( !(use_bmaps && bmap[i]) ) { 
                   unpacked_row[i] = base + scale*getbits( bp, pos, nbits );
                   new_pos = pos + nbits;
                   bp += (new_pos) / 8;
                   pos = (new_pos) % 8;
                }

            }
         }

         for ( i=0; i<cols; ++i) 
             unpacked_data[i+j*cols] = (double ) unpacked_row[i];

         bp = buf + n*4;

     } // End of NROWS for loop

     free( bmap );
     free( unpacked_row );
     free( buf );

     return;
}

