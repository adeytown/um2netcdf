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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "field_def.h"

/***
 *** SET_SOIL_LEVELS
 ***
 *** Set the depth of each model soil layer used in the output from the UM fields
 *** file.
 ***   INPUT:   ncid -> ID of the newly created NetCDF file
 ***               n -> number of soil levels to be filled
 ***              id -> index of the stored UM variable requiring this vertical
 ***                    dimension.
 ***
 ***   Mark Cheeseman, NIWA
 ***   December 29, 2012
 ***/

void set_soil_levels( int ncid, int n, int id ) {

     int    ierr, i, var_id, dim_id[1], nn;
     char   dim_name[12];
     float *buf;

  /** Create the appropriate NetCDF dimension **/
     sprintf( dim_name, "soil_level%d", n );
     ierr = nc_def_dim( ncid, dim_name, n, &dim_id[0] );

  /** Create a corresponding NetCDF variable **/
     ierr = nc_def_var( ncid, dim_name, NC_FLOAT, 1, dim_id, &var_id );

  /** Add some appropriate attributes **/
     ierr = nc_put_att_text( ncid, var_id,"standard_name", 5, "depth" );
     ierr = nc_put_att_text( ncid, var_id,        "units", 6, "meters" );
     ierr = nc_put_att_text( ncid, var_id,     "positive", 4, "down" );
     ierr = nc_put_att_text( ncid, var_id,         "axis", 1, "Z" );

  /** Fill the new NetCDF variable **/
     ierr = nc_enddef( ncid );

     buf = (float *) calloc( n,sizeof(float) );
     for ( i=0; i<n; i++ ) { 
         nn = (int ) (stored_um_fields[id].slices[i].level - 1);
         buf[i] = (float ) level_constants[3][nn]; 
     }

     ierr = nc_put_var_float( ncid, var_id, buf );
     free( buf );

  /** Need to set the ID of this vertical dimension for the stored UM variable **/
     stored_um_fields[id].z_dim = dim_id[0];

     return;
}


/***
 *** SET_PRESSURE_LEVELS
 ***
 *** Set the value of each model pressure layer used in the output from the UM fields
 *** file.
 ***   INPUT:   ncid -> ID of the newly created NetCDF file
 ***               n -> number of pressure levels to be filled
 ***              id -> index of the stored UM variable requiring this vertical
 ***                    dimension.
 ***
 ***   Mark Cheeseman, NIWA
 ***   December 29, 2012
 ***/

void set_pressure_levels( int ncid, int n, int id ) {

     int    i, ierr, var_id, dim_id[1];
     char   dim_name[10];
     float  pressure;
     size_t nc_index[1];

     sprintf( dim_name, "pressure%d", n );
     ierr = nc_inq_dimid( ncid, dim_name, &dim_id[0] ); 
     if ( ierr==NC_NOERR ) {
        stored_um_fields[id].z_dim = dim_id[0];
        return;
     }

  /** Create the appropriate NetCDF dimension **/
     ierr = nc_def_dim( ncid, dim_name, n, &dim_id[0] );

  /** Create a corresponding NetCDF variable **/
     ierr = nc_def_var( ncid, dim_name, NC_FLOAT, 1, dim_id, &var_id );

  /** Add some appropriate attributes **/
     ierr = nc_put_att_text( ncid, var_id,    "units", 2, "Pa" );
     ierr = nc_put_att_text( ncid, var_id,     "axis", 1, "Z" );
     ierr = nc_put_att_text( ncid, var_id,"standard_name", 12, "air_pressure" );
     ierr = nc_put_att_text( ncid, var_id, "positive", 4, "down" );
     ierr = nc_put_att_text( ncid, var_id, "long_name", 23, "standard level pressure" );

  /** Fill the new NetCDF variable **/
     ierr = nc_enddef( ncid );

     for ( i=0; i<n; i++ ) {
         pressure = (float ) stored_um_fields[id].slices[i].level; 
         nc_index[0] = (size_t) i;
         ierr = nc_put_var1_float( ncid, var_id, nc_index, &pressure );
     }

  /** Need to set the ID of this vertical dimension for the stored UM variable **/
     stored_um_fields[id].z_dim = dim_id[0];

     return;
}


/***
 *** SET_ALTITUDE
 ***
 *** Set the height of each model altitude layer used in the output from the UM fields
 *** file.
 ***   INPUT:   ncid -> ID of the newly created NetCDF file
 ***               n -> number of altitude levels to be filled
 ***              id -> index of the stored UM variable requiring this vertical
 ***                    dimension.
 ***
 ***   Mark Cheeseman, NIWA
 ***   December 29, 2012
 ***/

void set_altitude( int ncid, int n, int id ) {

     int    i, ierr, var_id, dim_id[1];
     char   dim_name[10];
     float  height;
     size_t nc_index[1];

  /** Get ID of the appropriate NetCDF dimension **/
     sprintf( dim_name, "altitude%d", n );
     ierr = nc_inq_dimid( ncid, dim_name, &dim_id[0] );
     if ( ierr==NC_NOERR ) {
        stored_um_fields[id].z_dim = dim_id[0];
        return;
     }

  /** Create the appropriate NetCDF dimension **/
     ierr = nc_def_dim( ncid, dim_name, n, &dim_id[0] );

  /** Create a corresponding NetCDF variable **/
     ierr = nc_def_var( ncid, dim_name, NC_FLOAT, 1, dim_id, &var_id );

  /** Add some appropriate attributes **/
     ierr = nc_put_att_text( ncid, var_id,    "units", 1, "m" );
     ierr = nc_put_att_text( ncid, var_id,     "axis", 1, "Z" );
     ierr = nc_put_att_text( ncid, var_id,"standard_name", 8, "altitude" );
     ierr = nc_put_att_text( ncid, var_id, "positive", 2, "up" );
     ierr = nc_put_att_text( ncid, var_id, "long_name", 22, "height above sea level" );

  /** Fill the new NetCDF variable **/
     ierr = nc_enddef( ncid );

     for ( i=0; i<n; i++ ) {
         height = (float ) stored_um_fields[id].slices[i].level;
         nc_index[0] = (size_t) i;
         ierr = nc_put_var1_float( ncid, var_id, nc_index, &height );
     }

  /** Need to set the ID of this vertical dimension for the stored UM variable **/
     stored_um_fields[id].z_dim = dim_id[0];

     return;
}


/***
 *** SET_HYBRID_LEVELS
 ***
 *** Set the height of each model hybrid layer used in the output from the UM fields
 *** file.
 ***   INPUT:   ncid -> ID of the newly created NetCDF file
 ***               n -> number of hybrid levels to be filled
 ***              id -> index of the stored UM variable requiring this vertical
 ***                    dimension.
 ***            grid -> type of Arakawa Grid on which variable is defined
 ***
 ***   Mark Cheeseman, NIWA
 ***   April 29, 2014
 ***/

void set_hybrid_levels( int ncid, int n, int id, unsigned short grid ) {

     int    i, ierr, var_id, dim_id[1], ind;
     char   dim_name[8];
     float  *height;

     sprintf( dim_name, "hybrid%d", n );
     ierr = nc_inq_dimid( ncid, dim_name, &dim_id[0] ); 
     if ( ierr==NC_NOERR ) {
        stored_um_fields[id].z_dim = dim_id[0];
        return;
     }

  /** Create the appropriate NetCDF dimension **/
     ierr = nc_def_dim( ncid, dim_name, n, &dim_id[0] );

  /** Create a corresponding NetCDF variable **/
     ierr = nc_def_var( ncid, dim_name, NC_FLOAT, 1, dim_id, &var_id );

  /** Add some appropriate attributes **/
     ierr = nc_put_att_text( ncid, var_id,    "units",  6, "meters" );
     ierr = nc_put_att_text( ncid, var_id, "positive",  2, "up" );
     ierr = nc_put_att_text( ncid, var_id,     "axis",  1, "Z" );
     ierr = nc_put_att_text( ncid, var_id,"long_name", 22, "height above sea level" );

  /** Fill the new NetCDF variable **/
     ierr = nc_enddef( ncid );

     height = (float *) malloc( n*sizeof(float) );

     ind = 4;
     if ( (grid==11)||(grid==18)||(grid==19) ) { ind = 6; } 
     for ( i=0; i<n; i++ ) 
         height[i] = (float ) level_constants[ind][stored_um_fields[id].slices[i].level]; 

     ierr = nc_put_var_float( ncid, var_id, height );
     free( height );

  /** Need to set the ID of this vertical dimension for the stored UM variable **/
     stored_um_fields[id].z_dim = dim_id[0];

     return;
}


/***
 *** SET_SURFACE_LEVELS
 ***
 *** Set the value of each model surface layer used in the output from the UM fields
 *** file.
 ***   INPUT:   ncid -> ID of the newly created NetCDF file
 ***               n -> number of surface levels to be filled
 ***              id -> index of the stored UM variable requiring this vertical
 ***                    dimension.
 ***
 ***   Mark Cheeseman, NIWA
 ***   December 29, 2012
 ***/

void set_surface_levels( int ncid, int n, int id ) {

     int    i, ierr, var_id, dim_id[1];
     char   dim_name[9];
     float *buf;

  /** Create the appropriate NetCDF dimension **/
     sprintf( dim_name, "surface%d", n );
     ierr = nc_def_dim( ncid, dim_name, n, &dim_id[0] );

  /** Create a corresponding NetCDF variable **/
     ierr = nc_def_var( ncid, dim_name, NC_FLOAT, 1, dim_id, &var_id );

  /** Add some appropriate attributes **/
     ierr = nc_put_att_text( ncid, var_id,     "axis", 1, "Z" );

  /** Fill the new NetCDF variable **/
     ierr = nc_enddef( ncid );

     buf = (float *) malloc( n*sizeof(float) );
     for ( i=0; i<n; i++ ) { buf[i] = (float ) i; }

     ierr = nc_put_var_float( ncid, var_id, buf );
     free( buf );

  /** Need to set the ID of this vertical dimension for the stored UM variable **/
     stored_um_fields[id].z_dim = dim_id[0];

     return;
}


/***
 *** SET_SEA_SURFACE_LEVELS
 ***
 *** Set the value of each model sea-surface layer used in the output from the UM fields
 *** file.
 ***   INPUT:   ncid -> ID of the newly created NetCDF file
 ***               n -> number of sea-surface levels to be filled
 ***              id -> index of the stored UM variable requiring this vertical
 ***                    dimension.
 ***
 ***   Mark Cheeseman, NIWA
 ***   December 29, 2012
 ***/

void set_sea_surface_levels( int ncid, int n, int id ) {

     int    i, ierr, var_id, dim_id[1];
     char   dim_name[13];
     float *buf;

  /** Create the appropriate NetCDF dimension **/
     sprintf( dim_name, "sea_surface%d", n );
     ierr = nc_def_dim( ncid, dim_name, n, &dim_id[0] );

  /** Create a corresponding NetCDF variable **/
     ierr = nc_def_var( ncid, dim_name, NC_FLOAT, 1, dim_id, &var_id );

  /** Add some appropriate attributes **/
     ierr = nc_put_att_text( ncid, var_id,    "units", 6, "meters" );
     ierr = nc_put_att_text( ncid, var_id, "positive", 2, "up" );
     ierr = nc_put_att_text( ncid, var_id,     "axis", 1, "Z" );

  /** Fill the new NetCDF variable **/
     ierr = nc_enddef( ncid );

     buf = (float *) malloc( n*sizeof(float) );
     for ( i=0; i<n; i++ ) { buf[i] = (float ) i; }

     ierr = nc_put_var_float( ncid, var_id, buf );
     free( buf );

  /** Need to set the ID of this vertical dimension for the stored UM variable **/
     stored_um_fields[id].z_dim = dim_id[0];

     return;
}

/***
 *** SET_VERTICAL_DIMENSIONS
 ***
 *** Function that determines how many vertical dimensions are required to
 *** describe the datafields in the UM fields file.  The dimensions are created
 *** in the new NetCDF file.
 ***
 *** INPUT:  ncid    -> id of the newly created NetCDF file
 ***         dim_ids -> array containing the NetCDF dimension IDs
 ***         ndim    -> # of previously defined NetCDF dimensions
 ***         rflag   -> equal to 1 if 32-bit output precision requested by user
 ***
 ***    Mark Cheeseman, NIWA
 ***    December 19, 2013
 ***/

void set_vertical_dimensions( int ncid, int rflag ) {

     int    i, num_instances;

    /*** Check each stored UM variable its 'z-coordinate' ***/
     for ( i=0; i<num_stored_um_fields; i++ ) {

       /** determine # of 'z-levels' for this variable */
         num_instances = stored_um_fields[i].num_slices / num_timesteps;

       /** Create new dimension if variable is 4D **/
         if ( num_instances>1 ) {
            switch ( stored_um_fields[i].lbvc ) {
                 case 1:
                      set_altitude( ncid, num_instances, i );
                      break;
                 case 6:
                      set_soil_levels( ncid, num_instances, i );
                      break;
                 case 8:
                      set_pressure_levels( ncid, num_instances, i );
                      break;
                 case 65:
                      set_hybrid_levels( ncid, num_instances, i, stored_um_fields[i].grid_type );
                      break;
                 case 128:
                      set_sea_surface_levels( ncid, num_instances, i );
                      break;
                 case 129:
                      set_surface_levels( ncid, num_instances, i );
                      break;
            }

         }

     }

    /*** Define dimensions to hold length of ETA arrays ***/
     i = nc_def_dim( ncid, "num_etaT_levels", (size_t ) header[110],
                    &num_instances );
     i = nc_def_dim( ncid, "num_etaR_levels", (size_t ) (header[110]-1),
                    &num_instances );

     return;
}
