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


/** FLAG_DEF
 ** 
 ** Header file containing global flag variables used for the configuration 
 ** of an UM2NetCDF run 
 **
 **   Mark Cheeseman, NIWA
 **   July 24, 2014
 **==========================================================================*/

int netcdf3_flag;  /* Flag variable denoting whether the output file should NOT
                      have any NetCDF4 features (HDF5 chunking and/or compression). */

