!
!===============================================================================
subroutine dielectric_const_water( frequency, temp_surf, salinity, surf_type )
!===============================================================================
! calculates the complex dielectric constant of water
! as a function of temperature, salinity, and frequency
! complex dielectric constant of water from Ray, Applied Optics Vol 11, No 8
!
! History:
!   4/1/1997 M. Klein Rewrote to Delphi Pascal
!   10/5/2020 Kevin Schaefer commented code
!  10/16/2020 Kevin Schaefer switched to variable module and replaced duplicate variables
!  12/1/2020  Kevin Schaefer corrected speed of light
!  12/1/2020  Kevin Schaefer removed negative sign for imaginary dielectric assignment
!  12/1/2020  Kevin Schaefer removed salinity conversion from mole fraction to parts per thousand (ppt)
!  12/1/2020  Kevin Schaefer corrected temperature conversion from Kelvin to centigrade
!  12/1/2020  Kevin Schaefer added seawater correction for normality
!  12/1/2020  Kevin Schaefer changed from temp_del to temp_surf relaxed wavelength and alpha
!------------------------------------------------------------------------------
Use dotlrt_variables
implicit none

! input variables
  real(8) frequency  ! (GHz) frequency < c/1000 microns
  real(8) temp_surf  ! (K) surface or skin temperature
  real(8) salinity   ! (ppt) ocean salinity
  integer surf_type  ! (-) surface type code (2=land, 1=lake, 0=ocean)

! internal variables
  real(8) freq_hz    ! (Hz) frequency in Hz
  real(8) temp_cent  ! (deg C) surface temperature
  real(8) temp_del   ! (deg C) delta temp from reference temperature
  real(8) wavelength ! (m) wavelength
  real(8) wave_relax ! (m) Relaxation wavelength 
  real(8) wave_relax_scale ! (-) Relaxation wavelength scaling factor
  real(8) e_infinity ! (-) high frequency dielectric constant
  real(8) e_static   ! (-) Static dielectric Constant
  real(8) e_static_scale ! (-) static dielectric constant scaling factor
  real(8) alpha      ! (-) spread parameter
  real(8) sigma      ! (1/m) Frequency independent conductivity
  real(8) sig_base   ! (1/m) reference independent conductivity
  real(8) sig_a      ! (1/C) constant conductivity salinity scaling factor
  real(8) sig_b      ! (1/C/ppt) slope conductivity salinity scaling factor
  real(8) sig_scale  ! (-) ocean conductivity salinity scaling factor
  real(8) norm       ! (mole charge/L) normality of seawater
  real(8) sin_a      ! (-) sin(alpha) term of dielectric calculation
  real(8) cos_a      ! (-) cos(alpha) term of dielectric calculation
  real(8) tau        ! (-) relaxation scale of dielectric calculation
  real(8) e_denom    ! (-) denominator of dielectric calculation
  real(8) e_real     ! (-) real component of dielectric constant
  real(8) e_imag     ! (-) imaginary component of dielectric constant

! unit conversions
  freq_hz = frequency*1.0e9       ! (Hz) frequency in Hz
  temp_cent = temp_surf - 273.15  ! (deg C) surface temperature
  wavelength = speed_lt/freq_hz   ! (m) wavelength

! (-) high frequency dielectric constant [Ray, 1972, eqn 7a]
  e_infinity = 5.27137d0 + 2.16474d-2*temp_cent - 1.31198d-3*temp_cent*temp_cent
    
! (deg C) delta temp from reference temperature [Ray, 1972, eqn 4]
  temp_del = temp_cent - 25.0d0
    
! (-) Static dielectric Constant [Ray, 1972, eqn 4]
  e_static = 78.54d0 * (1.0 - 4.579d-3*temp_del + 1.19d-5*temp_del**2 - 2.8d-8*temp_del**3 )
    
! (-) spread parameter [Ray, 1972, eqn 7b]
  alpha = 6.09265d-2 - 16.8129d0/temp_surf

! (m) Relaxation wavelength [Ray, 1972, eqn 7c] (converted from cm to m)
  wave_relax = 3.3836d-6*dexp(2513.98d0/temp_surf)

  if (surf_type == 1) then ! lake
    ! (1/m) Frequency independent conductivity for lake
    sigma = 1.25664d9   
  endif

  if (surf_type == 0) then ! ocean
    ! (mho/m) ionic conductivity of seawater [Klein and Swift, 1978, eqn 9-12; Stogryn, 1971, eqn 10]
    sig_base = 0.182521d0*salinity-1.46192d-3*salinity**2+2.09324d-5*salinity**3-1.28205d-7*salinity**4
    sig_a = 2.033d-2-1.266d-4*temp_del+2.464d-6*temp_del*temp_del
    sig_b = 1.849d-5+2.551d-7*temp_del+2.551d-8*temp_del*temp_del
    sig_scale = dexp(temp_del*(sig_a-salinity*sig_b))
    sigma = sig_base*sig_scale/e_free

    ! (mole charge/L) normality of seawater [Klein and Swift, 1978, eqn 19]
    norm = 0.9141d0*(1.707d-2*salinity + 1.205d-5*salinity**2 + 4.058d-9*salinity**3)
        
    ! static dielectric constant [Stogryn, 1971, eqn 4] for T <40C only
    e_static_scale = 1.0d0-0.2551d0*norm+5.151d-2*norm**2-6.889d-3*norm**3

    ! scale the static dielectric [Klein and Swift, 1978, eqn 13]
    e_static = e_static * e_static_scale
        
    ! scale the static dielectric [Klein and Swift, 1978, eqn 13]
    wave_relax_scale=1.0d0+0.1463d-2*temp_cent*norm-0.04896d0*norm-0.02967d0*norm**2+5.644d-3*norm**3
    wave_relax = wave_relax * wave_relax_scale
  endif

! Components of dielectric calculations [Ray, 1972, eqn 5-6]
  sin_a = dsin(alpha*90d0*pi/180d0)
  cos_a = dcos(alpha*90d0*pi/180d0)
  tau = (wave_relax/wavelength)**(1.0d0-alpha)
  e_denom = 1.0d0 + 2d0*tau*sin_a + tau*tau

! calculate dielectric [Ray, 1972, eqn 5-6]
  e_real = e_infinity+(e_static-e_infinity)*(1.0+tau*sin_a)/e_denom
  e_imag = (e_static-e_infinity)*tau*cos_a/e_denom + sigma*wavelength*1.88496d-11

! transfer to variable tree
  surf%diel = dcmplx(e_real,e_imag)

end subroutine dielectric_const_water
