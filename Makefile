SHELL=/bin/sh
# 
# List of Source objects for each main calling sequence 
# (main first, subroutines second)

<<<<<<< HEAD
test= DOTLRT_variables.o profiles.o DOTLRT_output.o scan_mod.o scan.o scan_control.o scan_routines.o scan_var_assign.o scan_var_control.o  scan_read_input.o absh2o.o absn2.o calc_fbw_temp_weight_scat.o calc_mon_temp_weight_scat.o calc_passband_freq.o calcprofile.o calc_tot_ext.o configure.o construct_surf.o core95.o d3lec.o do_tb_gvh94.o Dtb94.o fresnel_refl.o gaussj.o gaussq.o GeoJacobian.o get_instr_spec.o GetProfile.o Hg_phmat.o hydro_master_derivatives.o InterpSR.o Jacobian.o jacobi.o kirchoff_ocean_refl.o mrt.o o2abs.o rd.o rf.o RT_Jacobian.o Tb94.o handle_err.o read_WRF_netcdf.o get_data.o read_namel.o write_output.o

# Tunable parameters
#
# FC		Name of the fortran compiling system to use
# FFLAGS	Flags to the compiler
#               -O optimize
#               -o output specified
#               -g use debugger
# LDFLAGS	Flags to the loader
# LIBS		List of libraries
# CMD		Name of the executable
#
# SGI unique

# Linux unique
NETCDFINC = /lib64/gfortran/modules
NETCDFLIB = /usr/lib64/
PGPLOTLIB = /usr/lib64/
PGPLOTINC = /usr/lib64/

FC =		gfortran
#FFLAGS =	-Wall -Wextra -I$(NETCDFINC) -I$(PGPLOTINC) -fbounds-check
#LIBS =		-L/home/schaefek/LIB/ -L$(NETCDFLIB) -lnetcdf -lnetcdff -L$(PGPLOTLIB) -lpgplot 
FFLAGS =	-Wall -Wextra -I$(NETCDFINC) -Wno-compare-reals -Wno-do-subscript -fbounds-check
LIBS =		-L$(NETCDFLIB) -lnetcdf -lnetcdff 
LDFLAGS =       	
=======
# netcdf crap /home/kevin/LIB/netcdf
#NETCDFINC = /home/kevin/LIB/netcdf/include
#NETCDFLIB = /home/kevin/LIB/netcdf/lib
#NETCDFbin = /home/kevin/LIB/netcdf/bin

NETCDF = ${SSEC_NETCDF_DIR}

NETCDFINC = ${NETCDF}/include
NETCDFLIB = ${NETCDF}/lib
NETCDFbin = ${NETCDF}/bin

FC =            mpiifort
#FC =		gfortran
#FFLAGS =	-Wall -Wextra -I$(NETCDFINC) -fbounds-check 	
FFLAGS =        -I$(NETCDFINC)
LIBS =		-L$(NETCDFLIB) -lnetcdf -lnetcdff
>>>>>>> DOTLRTv2.1-MPI
#
# end machine uniquescan_var_assign.o
CMD =		a.out
.f:
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $<

# scan files
#scan_plot.o: scan_plot.f90
#	$(FC) $(FFLAGS) -c scan_plot.f90
scan.o: scan.f90
	$(FC) $(FFLAGS) -c scan.f90
scan_routines.o: scan_routines.f90
	$(FC) $(FFLAGS) -c scan_routines.f90
scan_mod.o: scan_mod.f90
	$(FC) $(FFLAGS) -c scan_mod.f90
scan_read_input.o: scan_read_input.f90
	$(FC) $(FFLAGS) -c scan_read_input.f90
scan_control.o: scan_control.f90
	$(FC) $(FFLAGS) -c scan_control.f90
scan_var_control.o: scan_var_control.f90
	$(FC) $(FFLAGS) -c scan_var_control.f90
scan_var_assign.o: scan_var_assign.f90
	$(FC) $(FFLAGS) -c scan_var_assign.f90


# DOTLRT Files
DOTLRT_variables.o: DOTLRT_variables.f90
	$(FC) $(FFLAGS) -c DOTLRT_variables.f90
profiles.o: profiles.f90
	$(FC) $(FFLAGS) -c profiles.f90
DOTLRT_output.o: DOTLRT_output.f90
	$(FC) $(FFLAGS) -c DOTLRT_output.f90
DOTLRT_main.o: DOTLRT_main.f90
	$(FC) $(FFLAGS) -c DOTLRT_main.f90
absh2o.o: absh2o.f90
	$(FC) $(FFLAGS) -c absh2o.f90
absn2.o: absn2.f90
	$(FC) $(FFLAGS) -c absn2.f90
calc_fbw_temp_weight_scat.o: calc_fbw_temp_weight_scat.f90
	$(FC) $(FFLAGS) -c calc_fbw_temp_weight_scat.f90
calc_mon_temp_weight_scat.o: calc_mon_temp_weight_scat.f90
	$(FC) $(FFLAGS) -c calc_mon_temp_weight_scat.f90
calc_passband_freq.o: calc_passband_freq.f90
	$(FC) $(FFLAGS) -c calc_passband_freq.f90
calcprofile.o: calcprofile.f90
	$(FC) $(FFLAGS) -c calcprofile.f90
calc_tot_ext.o: calc_tot_ext.f90
	$(FC) $(FFLAGS) -c calc_tot_ext.f90
configure.o: configure.f90
	$(FC) $(FFLAGS) -c configure.f90
construct_surf.o: construct_surf.f90
	$(FC) $(FFLAGS) -c -ffree-line-length-none construct_surf.f90
core95.o: core95.f90
	$(FC) $(FFLAGS) -c core95.f90
d3lec.o: d3lec.f90
	$(FC) $(FFLAGS) -c d3lec.f90
do_tb_gvh94.o: do_tb_gvh94.f90
	$(FC) $(FFLAGS) -c do_tb_gvh94.f90
Dtb94.o: Dtb94.f90
	$(FC) $(FFLAGS) -c Dtb94.f90
fresnel_refl.o: fresnel_refl.f90
	$(FC) $(FFLAGS) -c fresnel_refl.f90
gaussj.o: gaussj.f90
	$(FC) $(FFLAGS) -c gaussj.f90
gaussq.o: gaussq.f90
	$(FC) $(FFLAGS) -c gaussq.f90
GeoJacobian.o: GeoJacobian.f90
	$(FC) $(FFLAGS) -c GeoJacobian.f90
get_instr_spec.o: get_instr_spec.f90
	$(FC) $(FFLAGS) -c get_instr_spec.f90
GetProfile.o: GetProfile.f90
	$(FC) $(FFLAGS) -c GetProfile.f90
Hg_phmat.o: Hg_phmat.f90
	$(FC) $(FFLAGS) -c Hg_phmat.f90
hydro_master_derivatives.o: hydro_master_derivatives.f90
	$(FC) $(FFLAGS) -c hydro_master_derivatives.f90
InterpSR.o: InterpSR.f90
	$(FC) $(FFLAGS) -c InterpSR.f90
Jacobian.o: Jacobian.f90
	$(FC) $(FFLAGS) -c Jacobian.f90
jacobi.o: jacobi.f90
	$(FC) $(FFLAGS) -c jacobi.f90
kirchoff_ocean_refl.o: kirchoff_ocean_refl.f90
	$(FC) $(FFLAGS) -c kirchoff_ocean_refl.f90
mrt.o: mrt.f90
	$(FC) $(FFLAGS) -c mrt.f90
o2abs.o: o2abs.f90
	$(FC) $(FFLAGS) -c o2abs.f90
read_namel.o: read_namel.f90
	$(FC) $(FFLAGS) -c read_namel.f90
rd.o: rd.f90
	$(FC) $(FFLAGS) -c rd.f90
rf.o: rf.f90
	$(FC) $(FFLAGS) -c rf.f90
RT_Jacobian.o: RT_Jacobian.f90
	$(FC) $(FFLAGS) -c RT_Jacobian.f90
Tb94.o: Tb94.f90
	$(FC) $(FFLAGS) -c Tb94.f90
handle_err.o: handle_err.f90
	$(FC) $(FFLAGS) -c handle_err.f90
read_WRF_netcdf.o: read_WRF_netcdf.f90
	$(FC) $(FFLAGS) -c read_WRF_netcdf.f90
get_data.o: get_data.f90
	$(FC) $(FFLAGS) -c get_data.f90
write_output.o: write_output.f90
	$(FC) $(FFLAGS) -c write_output.f90

# make command list

test: $(test)
	$(FC) $(FFLAGS) $(test) $(LIBS) -o test.out

# Cleanup stuff
clean:
	-rm -f *.o *.out *.mod *.stb core


