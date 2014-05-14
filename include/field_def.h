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

/*---------------------------------------------------------------------------*
 *  VARIABLES                                                                *
 *---------------------------------------------------------------------------*/

int num_xml_vars;         /* # of possible defined in the XML stash file */
int num_um_vars;          /* # of variables defined in the UM file */ 
int wordsize;             /* wordsize of the architecture on which the input UM fields file was created: 4 or 8 bytes */
int num_timesteps;        /* # of timesteps present in the UM fields file */
int num_stored_um_fields; /* # of UM variables found & processed in input UM fields file */

long header[256];        /* input UM fields file's header */
long int_constants[46];  /* input UM fields file's integer_constants array */

double **level_constants;   /* input UM fields file's 2d level-dependent constants array */
double   real_constants[6]; /* input UM fields files's real_constants array */


/*---------------------------------------------------------------------------*
 *  FUNCTION POINTERS                                                        *
 *---------------------------------------------------------------------------*/

double* (*field_interpolation)( double*, int, int, double, double );  /* ptr to appropriate interpolation procedure */
void    (*endian_swap)( void*, int );                 /* ptr to appropriate endian swap procedure   */


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
       char varname[45];   /* name of field in UM output file */ 
       char longname[100]; /* full descriptive name of field */
       char stdname[75];   /* CF-compliant name of field */
       char units[25];     /* units for field */
       char coord[25];     /* coordinate system used to describe field */
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
        long reclength;        /* Size of the storec dataslice in records */
        long size;             /* Size of the slice in words */
        long lbproc;           /* denotes whether post-processing has been performed on variable.  0 if not */
        double mdi;            /* value used to denote a missing data point */
        struct tm validity, datatime;
} um_dataslice;


/**
 ** UM_variable - Struct that contains essential attributes for a variable 
 **               stored in the input UM fields file.
 **/ 

typedef struct um_variable {
        char           name[55];   /* name of the UM variable */
        unsigned short num_slices; /* # of 2D slices beinging to this variable */
        unsigned short xml_index;  /* location of field in the XML stash description file */ 
        unsigned short grid_type;  /* used for interpolation; indicates the Arakawa grid used for the data */
        unsigned short nx, ny;     /* # of points in the x (lon) and y (lat) directions */
        unsigned short accum;      /* indicates if field is a product of some sort of accummulation process */
        long lbproc;     /* denotes whether post-processing has been performed on variable.  0 if not */
        unsigned short lbvc;       /* indicates the vertical coordinate system used */
        unsigned short stash_code;
        unsigned short coordinates;/* LBCODE that describes the coordinate system used to map field */
        int            z_dim ;        /* ID of the vertical coordinate of variable */
        int            t_dim ;        /* ID of the "time" dimension of variable */
        float          space_bnds[2]; /* start & end points for either spatial accumulation done in the variable */
        float          *time_offsets; /* array containing temporal offset of each timesetp of variable from */
                                      /* the initial validity date & time. [in hours] */
        float          **time_bnds;   /* difference between start & end times for accummulation during */
                                      /* each stored timestep of the variable */
        um_dataslice   *slices;       /* array of data slices belonging to this UM variable */
        nc_type        vartype;       /* datatype of the UM variable (float/double/int/long) */     
} um_variable;

um_variable *stored_um_fields;


typedef struct new_um_variable {
        char           name[55];   /* name of the UM variable */
        unsigned short stash_code;
        unsigned short xml_index;  /* location of field in the XML stash description file */ 
        unsigned short nt;
        unsigned short nz;
        unsigned short nx, ny;     /* # of points in the x (lon) and y (lat) directions */
        unsigned short lbvc;       /* indicates the vertical coordinate system used */
        unsigned short accum;      /* indicates if field is a product of some sort of accummulation process */
        unsigned short coordinates;
        unsigned short grid_type;  /* used for interpolation; indicates the Arakawa grid used for the data */
        unsigned short lbproc;     /* denotes whether post-processing has been performed on variable.  0 if not */
        unsigned short t_dim ;     /* ID of the "time" dimension of variable */
        unsigned short z_dim ;     /* ID of the vertical coordinate of variable */
        float *times;
        float **time_bnds;         /* difference between start & end times for accummulation during */
                                   /* each stored timestep of the variable */
        float space_bnds[2];       /* start & end points for either spatial accumulation done in the variable */
        um_dataslice **slices;
        nc_type        vartype;    /* datatype of the UM variable (float/double/int/long) */     
} new_um_variable;

new_um_variable *stored_um_vars;
