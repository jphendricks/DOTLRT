SHELL=/bin/sh

###  Directory for netcdf library

## Agens
# NCDIR =  ${SSEC_NETCDF_DIR}

## CET
NCDIR = /usr/lib64/gfortran/modules

# libraries
LIBS = -L$(NCDIR) -lnetcdf -lnetcdff
LDFLAGS =

### Fortran compiler

## Sequential
FC = gfortran

## MPI
# FC = /usr/lib64/mpich/bin/mpif90

### Fortran Flags

## Debug
FFLAGS = -g -Wall -I$(NCDIR) -Wunused-parameter -ffree-line-length-none -fbounds-check

## Minimal Flags -ffpe-trap=invalid,zero,overflow 
# FFLAGS = -I$(NCDIR) -Wunused-parameter

### Source Files
SRC = \
	absh2o.f90 \
	absn2.f90 \
	gaussj.f90 \
	gaussq.f90 \
	handle_err.f90 \
	hydro_master_derivatives.f90 \
	jacobi.f90 \
	o2abs.f90 \
	rd.f90 \
	rf.f90 \
	InterpSR.f90 \
	calc_fbw_temp_weight_scat.f90 \
	calc_mon_temp_weight_scat.f90 \
	calc_passband_freq.f90 \
	calc_tot_ext.f90 \
	calcprofile.f90 \
	d3lec.f90 \
	do_tb_gvh94.f90 \
	get_instr_spec.f90 \
	fresnel_refl.f90 \
	get_data.f90 \
	kirchoff_ocean_refl.f90 \
	mrt.f90 \
	read_WRF_netcdf.f90 \
	read_namel.f90 \
	Dtb94.f90 \
	GeoJacobian.f90 \
	GetProfile.f90 \
	Hg_phmat.f90 \
	Jacobian.f90 \
	RT_Jacobian.f90 \
	Tb94.f90 \
	write_output.f90 \
	configure.f90 \
	construct_surf.f90 \
	core95.f90

SCAN = \
	scan_control.f90 \
	scan_read_input.f90 \
	scan_routines.f90 \
	scan_var_assign.f90 \
	scan_var_control.f90

MODS = \
        netcdf_utilities_mod.f90 \
	DOTLRT_output.f90 \
	DOTLRT_variables.f90 \
	scan_mod.f90\
	profiles.f90 

### List of .f90 files for program
scan = $(MODS) $(SRC) $(SCAN)
dot = $(MODS) $(SRC) 

### Variable that takes all .f90 files and replaces extension with .o
SCAN_OBJS = $(patsubst %.f90,%.o,$(scan))
DOT_OBJS = $(patsubst %.f90,%.o,$(dot))

### Make Targets
scan: $(SCAN_OBJS) scan.f90
	$(FC) $(FFLAGS) $(SCAN_OBJS) scan.f90 $(LIBS) -o scan

dot: $(DOT_OBJS) DOTLRT_main.f90
	$(FC) $(FFLAGS) $(DOT_OBJS) DOTLRT_main.f90 $(LIBS) -o dot

## Compile all .f90 files into .o files
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	-rm -f *.o *.out *.mod *.stb core 


