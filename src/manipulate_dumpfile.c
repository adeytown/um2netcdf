#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

long header[256], **lookup;
void (*endian_swap)( void*, int );


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

void no_endian_swap( void *ptr, int nchunk ) { return; }

/***
 *** GET_FILE_ENDIANNESS_WORDSIZE
 ***
 *** Function that determines the appropriate wordsize and endianness of the
 *** user-supplied input UM dump file.
 ***
 *** INPUT:  fh -> file handle/pointer of the input UM dump file
 ***
 ***   Mark Cheeseman, NIWA
 ***   June 18, 2014
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
       return word_size;
    } else {
       endian_swap = &endian_swap_8bytes;
       endian_swap( header, 256 );
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
       return word_size;
    } else {
       endian_swap = &endian_swap_4bytes;
       endian_swap( header, 256 );
       if ( (header[1]==1)&&(header[150]==64) ) {
          return word_size;
       }
    }

    word_size = -1;
    return word_size;

}


/***
 *** USAGE
 ***
 *** Displays the command-line flags available
 ***
 ***   Mark Cheeseman, NIWA
 ***   June 18, 2014
 ***/

void usage() {
     printf( "\nUsage:  manip.x [ OPTIONS ] <dump-file> \n\n" );
     printf( "  where\n\n" );
     printf( "    dump-file  --> UM dumpfile containing the field that is to be altered.\n\n" );
     printf( "  The following options can be specified.  They MUST appear before the dump filename!\n\n" );
     printf( "    -h used to display this help message\n" );
     printf( "    -s used to specify the stash code of UM variable that is to be manipulated\n\n" );
}


/*** 
 *** MANIPULATE DUMPFILE
 ***
 *** C program that allows user to apply a specified function/process on a field that
 *** is selected by STASH CODE from the commandline.
 ***
 ***   Mark Cheeseman, NIWA
 ***   June 18, 2014
 ***/

int main( int argc, char *argv[] ) {

    int     c, wordsize, nrec, num_vals;
    long    stash_code;
    double *data_val;
    FILE *fid;

    if ( argc<3 ) {
        printf( "ERROR: too few commandline arguments\n\n" );
        usage();
        exit(1);
     }

     while ( (c = getopt(argc,argv,"hs:")) != EOF ) {
           switch(c) {
               case 'h':
                       usage();
                       exit(1);
               case 's':
                       optind--;
                       stash_code = atoi( argv[optind] );
          }
     }

     printf( "Dump file to be manipulated: %s\n", argv[argc-1] );
     printf( "STASH CODE of variable     : %ld\n\n", stash_code );

     fid = fopen( argv[argc-1], "r" );
     if ( fid==NULL ) { return 0; }

     wordsize = get_file_endianness_wordsize( fid );
     if ( wordsize==-1 ) {
        printf( "ERROR: incorrect wordsize detected!\n" );
        exit(1);
     }

/**
 ** Read in the LOOKUP table
 **---------------------------------------------------------------------------*/
     fseek( fid, (header[149]-1)*wordsize, SEEK_SET );
     lookup = (long **) malloc( header[151]*sizeof(long *) );

     for ( nrec=0; nrec<header[151]; nrec++ ) {
         lookup[nrec] = (long *) malloc( header[150]*sizeof(long) );
         c = fread( lookup[nrec], wordsize, header[150], fid );
         endian_swap( lookup[nrec],header[150] );
     }
 
/**
 ** Look for the specified variable 
 **---------------------------------------------------------------------------*/
     for ( nrec=0; nrec<header[151]; nrec++ ) {
         if ( lookup[nrec][41]==stash_code ) {

            num_vals = (int ) (lookup[nrec][17]*lookup[nrec][18]);
            data_val = (double *) malloc( num_vals*sizeof(double) ); 

            fseek( fid, lookup[nrec][28], SEEK_SET );
            fread( data_val, wordsize, num_vals, fid );
            endian_swap( data_val, num_vals );

            manipulating_function_( &lookup[nrec][17], &lookup[nrec][18], &data_val );

            endian_swap( data_val, num_vals );
            fseek( fid, lookup[nrec][28], SEEK_SET );
            fwrite( data_val, wordsize, num_vals, fid );

            free( data_val );
         }
     }
     fclose( fid );

     return 1;
}
