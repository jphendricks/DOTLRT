SHELL=/bin/sh
# 
# List of Source objects for each main calling sequence 
# (main first, subroutines second)

#test = DOTLRT_variables.o profiles.o DOTLRT_output.o scan_mod.o scan_control.o scan_routines.o scan_var_assign.o scan_var_control.o  scan_read_input.o absh2o.o absn2.o calc_fbw_temp_weight_scat.o calc_mon_temp_weight_scat.o calc_passband_freq.o calcprofile.o calc_tot_ext.o configure.o construct_surf.o core95.o d3lec.o do_tb_gvh94.o Dtb94.o fresnel_refl.o gaussj.o gaussq.o GeoJacobian.o get_instr_spec.o GetProfile.o Hg_phmat.o hydro_master_derivatives.o InterpSR.o Jacobian.o jacobi.o kirchoff_ocean_refl.o mrt.o o2abs.o rd.o rf.o RT_Jacobian.o Tb94.o handle_err.o read_WRF_netcdf.o get_data.o read_namel.o write_output.o scan.o

# this is where the mod files for netcdf live
NCDIR = $(HOME)/local


FC = gfortran

# DEBUG FLAGS
#FFLAGS = -g -Wall -Wextra -I$(NCDIR)/include -Wuninitialized -Wunused -ffree-line-length-none -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow 

# OPTIMIZED FLAGS
FFLAGS = -I$(NCDIR)/include

LIBS = -L$(NCDIR)/lib -lnetcdf -lnetcdff
LDFLAGS =

## Linux unique
#NETCDFINC = /lib64/gfortran/modules
#NETCDFLIB = /usr/lib64/
##PGPLOTLIB = /usr/lib64/
##PGPLOTINC = /usr/lib64/
#
#FC =		gfortran
##FFLAGS =	-Wall -Wextra -I$(NETCDFINC) -I$(PGPLOTINC) -fbounds-check
##LIBS =		-L/home/schaefek/LIB/ -L$(NETCDFLIB) -lnetcdf -lnetcdff -L$(PGPLOTLIB) -lpgplot 
#FFLAGS =	-Wall -Wextra -I$(NETCDFINC) -Wno-compare-reals -Wno-do-subscript -fbounds-check
#LIBS =		-L$(NETCDFLIB) -lnetcdf -lnetcdff 
#LDFLAGS =       	

# scan files
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
	InterpSR.f90

SCAN = \
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
	scan_control.f90 \
	scan_read_input.f90 \
	scan_routines.f90 \
	scan_var_assign.f90 \
	configure.f90 \
	construct_surf.f90 \
	core95.f90 \
	scan_var_control.f90

	#scan.f90 \

MODS = \
	DOTLRT_output.f90 \
	DOTLRT_variables.f90 \
	netcdf_utilities_mod.f90 \
	profiles.f90 \
	scan_mod.f90

PROGS = \
	write_netcdf_file.f90 \
	simple_read.f90 \
	simple_write.f90 \
	DOTLRT_main.f90

OTHER = \

test = $(SRC) $(MODS) $(SCAN)
mods = $(MODS)
src  = $(SRC)
OBJS = $(patsubst %.f90,%.o,$(test))
FILES = $(patsubst %.f90,%.o,$(test))

# make command list
src:
	$(FC) $(FFLAGS) $(LIBS) -c $(test)

print:
	@echo $(OBJS)
	@echo $(NCDIR)

scan: $(OBJS) scan.f90
	$(FC) $(FFLAGS) $(OBJS) scan.f90 $(LIBS) -o scan

dot: $(OBJS) DOTLRT_main.f90
	$(FC) $(FFLAGS) $(OBJS) DOTLRT_main.f90 $(LIBS) -o dot 

ind: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) simple_read.f90 $(LIBS) -o test_index_table 

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	-rm -f *.o *.out *.mod *.stb core scan dot test_index_table


