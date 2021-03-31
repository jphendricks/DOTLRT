!
!=============================================================
subroutine read_namel_dotlrt
!=============================================================
! reads in dotlrt control variables and inputs
!
! history:
!  10/19/20 Kevin Schaefer created routine from init_grid from SiBCASA
!-------------------------------------------------------------
  use dotlrt_variables
  implicit none

!
! Namelists
namelist /dotlrt_control/ & ! execution configuration control
  flag_print_full, flag_read_anc, ocean_mod, prof_src

namelist /dotlrt_input_files/ & ! input file names
  file_instr, file_var
  
namelist /dotlrt_rad_tran_param/ & ! input radiative transfer parameters
  nstream, nstream_surf, obs_height, obs_theta, nsub_freq 

namelist /dotlrt_grid_param/ & ! input grid indeces
  chan_strt, chan_stop, lon_strt, lon_stop, lat_strt, lat_stop, flag_reduce_nvar, new_nlev

namelist /dotlrt_output_control/ & ! output control parameters
  save_sing_prof, save_ilon, save_ilat, save_rad_file, save_jac_file, save_sing_rad, save_prof_file

! read in namel_dotlrt
    open(unit=2,file='namel_dotlrt',form='formatted')

    print *,'Read Configuration Parameters'
    read (2,dotlrt_control)

    print *,'Read Input File Names'
    read (2,dotlrt_input_files)

    print *,'Read Radiative Transfer Parameters'
    read (2,dotlrt_rad_tran_param)

    print *,'Read Grid Parameters'
    read (2,dotlrt_grid_param)

    print *,'Read Output Control Parameters'
    read (2,dotlrt_output_control)

end subroutine read_namel_dotlrt
