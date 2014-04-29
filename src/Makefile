##============================================================================
##                 U M 2 N e t C D F  V e r s i o n 2 . 0
##                 --------------------------------------
##
##    Main author: Mark Cheeseman
##                 National Institute of Water & Atmospheric Research (Ltd)
##                 Wellington, New Zealand
##                 February 2014
##
##    UM2NetCDF is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    any later version.
##
##    UM2NetCDF is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    A copy of the GNU General Public License can be found in the main UM2NetCDF
##    directory.  Alternatively, please see <http://www.gnu.org/licenses/>.
##============================================================================*/

##-----------------------------------------------------------------------------
## Setup environment for the specific build architecture
##-----------------------------------------------------------------------------

default: bld_message exe_build check 

bld_message:
	@echo " "
	@echo "***************************************************************"
	@echo "                      UM2NETCDF BUILD"
	@echo "***************************************************************"
	@echo " "
	@echo " PLATFORM: $(ARCH)"
	@echo " COMPILER: $(CC)"
	@echo " "
	@echo " Compiling..."
	@echo "---------------------------------------------------------------"
include ../config/make.inc.$(ARCH)_$(CC)

##-----------------------------------------------------------------------------
## Objects Listing
##-----------------------------------------------------------------------------

OBJS =	util.o stashfile_operations.o umfile_operations.o interp.o \
	vertical_dimensions.o grid_operations.o netcdf_dimension_functions.o \
	wgdos_unpack.o netcdf_variable_functions.o convert_float_ibm_to_ieee32.o \
	netcdf_functions.o um2netcdf.o 

##-----------------------------------------------------------------------------
## Define the name and location of the um2netcdf binary 
##-----------------------------------------------------------------------------

BINARY=/hpcf/working/cheesemanmp/um2netcdf/NZCSM_ED/um2netcdf.x

##-----------------------------------------------------------------------------
## Global Compile Rules
##-----------------------------------------------------------------------------

INCS2 = $(INCS) -I../include

%.o: %.c
	$(CC) $(OPT_FLAGS) $(CPPFLAGS) $(DEFS) $(INCS2) -c $<

exe_build: $(OBJS)
	@echo " "
	@echo " Linking..."
	@echo "---------------------------------------------------------------"
	$(CC) $(OPT_FLAGS) -o $(BINARY) $(OBJS) $(LIBS) -lm 

clean:
	@rm -f *.o $(BINARY)

check:
	@echo " "
	@echo " "
	@if [ -a $(BINARY) ]; then echo "SUCCESS!"; echo "Binary located at $(BINARY)"; fi;
	@echo " "

##-----------------------------------------------------------------------------
## Individual File Compile Rules
##-----------------------------------------------------------------------------

util.o:
umfile_operations.o:
interp.o:
stashfile_operations.o:  
grid_operations.o:  
vertical_dimensions.o:
convert_float_ibm_to_ieee32.o:
wgdos.o:
netcdf_dimension_functions.o: vertical_dimensions.o
netcdf_variable_functions.o: interp.o umfile_operations.o wgdos.o 
netcdf_functions.o: interp.o grid_operations.o netcdf_dimension_functions.o netcdf_variable_functions.o
um2netcdf.o: util.o stashfile_operations.o umfile_operations.o netcdf_functions.o
