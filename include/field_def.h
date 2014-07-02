/**============================================================================
    This file is part of UM2NetCDF.

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


/** FIELD_DEF
 ** 
 ** Header file contianing constants, variables and structs relating to the
 ** layout of data in the input UM fieldsfile.
 **
 **   Mark Cheeseman, NIWA
 **   December 27, 2013
 **==========================================================================*/

#include <time.h>
#include <netcdf.h>
#include <stdint.h>

/*---------------------------------------------------------------------------*
 *  VARIABLES                                                                *
 *---------------------------------------------------------------------------*/

int num_unique_lats;

int num_xml_vars;         /* # of possible defined in the XML stash file */
int num_um_vars;          /* # of variables defined in the UM file */ 
int wordsize;             /* wordsize of the architecture on which the input UM fields file was created: 4 or 8 bytes */
int num_timesteps;        /* # of timesteps present in the UM fields file */
int num_stored_um_fields; /* # of UM variables found & processed in input UM fields file */

long header[256];        /* input UM fields file's header */
long int_constants[46];  /* input UM fields file's integer_constants array */

double **level_constants;   /* input UM fields file's 2d level-dependent constants array */
double   real_constants[6]; /* input UM fields files's real_constants array */

unsigned short int blacklist[10]; /* STASH CODES of UM variables to be avoided */
unsigned short int blacklist_cnt;

struct tm forecast_reference;

/*---------------------------------------------------------------------------*
 *  FUNCTION POINTERS                                                        *
 *---------------------------------------------------------------------------*/

void   (*field_interpolation)( double*, float*, int ); /* ptr to appropriate interpolation procedure */
void   (*endian_swap)( void*, int );                   /* ptr to appropriate endian swap procedure   */
void   (*endian_swap_4b)( void*, int );                /* ptr to appropriate 4 byte endian swap procedure   */
double (*ibm2ieee_convert)( uint32_t );                 /* ptr to appropriate IBM float to IEEE float function */

/*---------------------------------------------------------------------------*
 *  STRUCTS                                                                  *
 *---------------------------------------------------------------------------*/

/**
 ** UM_Field_Metadata - Struct that contains various pieces of metadata for a
 **                     UM variable.  This metadata is stored in the XML stash file.
 **/ 

typedef struct um_field_metadata {
       int code;           /* item code for the specific field */ 
       int model;          /* number of the UM model */
       int section;        /* section number where the field is located */
       int level_type;     /* denotes which type of depth level is used for field (used for Hybrid levels only) */
                           /*    level_type=1  -> points on model rho levels */
                           /*    level_type=2  -> points on model theta levels (default) */
       int umgrid;         /* does field need to be intepreted onto the Arakawa-C grid? */
                           /*    umgrid=11 -> intepolate U,V points on an Arakawa-B grid (middle of a Arakawa-C grid cell) */
                           /*    umgrid=18 -> intepolate U-points on North-South cell faces to lower left data point */
                           /*    umgrid=19 -> intepolate V-points on East-West cell faces to lower left data point */
                           /*      all other values of umgrid result in no interpolation */
       int accum;          /* is the variable the result of some sort of accummulation operation?  */
                           /*    accum=0  -> values are not accummulations */
                           /*    accum=1  -> values are sums */
                           /*    accum=2  -> values are maximums */
       float validmax;     /* max valid value for the field */
       float validmin;     /* min valid value for the field */
       float scale;        /* scale factor applied to files added to the output NetCDF file */
       char varname[45];   /* name of field in UM output file */ 
       char longname[100]; /* full descriptive name of field */
       char stdname[75];   /* CF-compliant name of field */
       char units[25];     /* units for field */
} um_field_metadata;

um_field_metadata *um_vars;

/**
 ** UM_dataslice - Struct that contains unique information about a 2D slice of an 
 **                UM variable stored in input UM fields file.
 **/ 

typedef struct um_dataslice {
        unsigned short id;
        unsigned short level;  /* Depth or z-level on which the data slice resides */
        unsigned short lbpack; /* Code used to denote the packing method used */ 
        long location;         /* Starting address of dataslice in UM fields file lookup[n][28] */
        long reclength;        /* Size of the stored dataslice in records */
        long size;             /* Size of the slice in words */
        long lbproc;           /* denotes whether post-processing has been performed on variable.  0 if not */
        double mdi;            /* value used to denote a missing data point */
        struct tm validity, datatime;
} um_dataslice;


/**
 ** UM_variable - Struct that contains essential attributes for a variable 
 **               stored in the input UM fields file.
 **/ 

typedef struct new_um_variable {
        char           name[55];   /* name of the UM variable */
        unsigned short stash_code;
        unsigned short xml_index;  /* location of field in the XML stash description file */ 
        unsigned short nt;
        unsigned short nz;
        unsigned short nx, ny;     /* # of points in the x (lon) and y (lat) directions */
        unsigned short lbvc;       /* indicates the vertical coordinate system used */
        unsigned short lbpack;     /* packing method used: 0 -> no packing, 1 -> WGDOS packing used */
        unsigned short accum;      /* indicates if field is a product of some sort of accummulation process */
        unsigned short coordinates;
        unsigned short grid_type;  /* used for interpolation; indicates the Arakawa grid used for the data */
        unsigned short lbproc;     /* denotes whether post-processing has been performed on variable.  0 if not */
        unsigned short t_dim ;     /* ID of the "time" dimension of variable */
        unsigned short z_dim ;     /* ID of the vertical coordinate of variable */
        unsigned short y_dim ;     /* ID of the latitudinal coordinate of variable */
        unsigned short x_dim ;     /* ID of the longitudinal coordinate of variable */
        unsigned short level_type; /* used to denote whether variable is on a RHO level (1) or a  */
                                   /* THETA level (2) */
        int lat_index;
        float *times;
        float **time_bnds;         /* difference between start & end times for accummulation during */
                                   /* each stored timestep of the variable */
        float space_bnds[2];       /* start & end points for either spatial accumulation done in the variable */
        float scale_factor;
        um_dataslice **slices;
        nc_type        vartype;    /* datatype of the UM variable (float/double/int/long) */     
} new_um_variable;

new_um_variable *stored_um_vars;


/**
 ** run_details - Struct that contains attributes about the UM run that produced  
 **               the input UM data file.
 **/ 

typedef struct run_details {
        char institution[25];  /* name of institution where field was generated */
        int  ps;               /* UK MetOffice PS version used in the run */
        int  eps;              /* EPS version number of the run */
        int  rose_id;          /* id number of the ROSE application used for the run */
        char model[15];        /* name of the specific UM model configuration used */
        char ref[100];         /* string listing location of input field reference documentation */
        char comment[100];     /* string containing additional info about run/config */
        char title[50];        /* title of the model run */
        char assim[50];        /* data assimilation method used in run */ 
} run_details;

run_details run_config;
