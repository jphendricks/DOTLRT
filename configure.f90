! configure.f90
! This routine establishesw configuration parameters for DOTLRTv1
! 27 May 2005    Bob Weber
subroutine configure( inpheight,         & ! observation height (AGL, km)
                      inptheta,          & ! observation angle from nadir (degrees)
                      numsbfreqs,        & ! number of frequency points per sideband
                      instrspec,         &
                      nstreams,          & ! number of streams (total, up and down)
                      stream_angles,     & ! stream angles
                      nsurfangles,       & ! number of angles for surface reflectivity
                      surfinp,           & ! surface reflectivity (vertical, horizontal polarization)
                      atm_inp_num_levels )
  use variables
  implicit none
  integer numsbfreqs, atm_inp_num_levels
  real(8) instrspec(5), theta_obs
  real(8) inpheight, inptheta
  integer kpts, kind, num_stream_angles
  integer j, joc, nstreams, nsurfangles, next_angle, last_angle
  real(8) endpts(2), alpha, beta, b(32), t(32), law, naw
  real(8) stream_angles(nstreams)
  real(8) surfinp(nsurfangles,3)

  character*120 debugout

  inp_height = inpheight
  inp_theta = inptheta

  atm_inp%num_levels = atm_inp_num_levels
  nlr = atm_inp%num_levels
  m1 = 1
  nlr1 = nlr*m1

!  write(debugout,*) "nlevels=",atm_inp_num_levels
!  call mexPrintf(debugout//achar(10))
  
  ! instrument specification parameters input file:
  call get_instr_spec( instrspec )

  num_sb_freqs = numsbfreqs

  d_albedo = 1.0d-7  ! lower limit (hydrometeor scattering/absorption total)
                     !  for diagonalization of layer instrument parameters


  number_h2o_phases = 5    ! number of hydrometeor phases
                           ! liquid             : cloud liquid, rain
                           ! solid              : ice
                           ! mixed liquid/solid : snow, grauppel
  nvar = 2 + 2 * number_h2o_phases ! 12
                                   ! number of variables for radiation Jacobian
                                   ! T (temperature)
                                   ! Ka (absorption coefficient)
                                   ! Ks (scatter coefficient)      X number_h2o_phases
                                   ! g (scatter asymmetry)         X number_h2o_phases

  nang = nstreams ! number of stream angles, up and down (even)
  nangover2=nang/2
  if( 2*nangover2 /= nang ) then
    write(8,*) 'nang is odd integer, program STOPPED'
    stop
  end if
  ! stream angles and weights for Gaussian-Lobatto quadrature, including
  ! streams at 0 and 180 degrees
    kpts = 2
    kind = 1
    endpts(1) = -1.0D0
    endpts(2) = +1.0D0
    call gaussq(kind, nang, alpha, beta, kpts, endpts, b, t, cris_quad_wghts)
    t(1) = -1.0d0
    t(nang) = +1.0d0
    do j = 1, nang
      quad_angle_array(j) = (180.0d0/pi)*acos(t(nang-j+1))
    end do

    !write(debugout,*) "w = "
    !call mexPrintf(debugout//achar(10))
    !do j = 1, nang
    !   write(debugout,*) quad_angle_array(j), cris_quad_wghts(j)
    !   call mexPrintf(debugout//achar(10))
    !end do
    
  ! Henyey-Greenstein phase matrix
    call HG_phmat()
    do j = 1, nang
      stream_angles(j) = teta(j)
    end do

  ! inteprolate surface reflectivities to stream angles
    do j = 1, nangover2
      theta_obs = stream_angles(j)
      if( theta_obs < surfinp(1,1) )           theta_obs = surfinp(1,1)
      if( theta_obs > surfinp(nsurfangles,1) ) theta_obs = surfinp(nsurfangles,1)
      ! Determine angle boundaries
      next_angle = 2
      do while( ( surfinp(next_angle,1) < theta_obs   ) .and. &
                (          next_angle < nsurfangles ) )
           next_angle = next_angle + 1
      end do
      last_angle = next_angle - 1
      law = (surfinp(next_angle,1) -            theta_obs) &
          / (surfinp(next_angle,1) - surfinp(last_angle,1) )
      naw = 1.0d0 - law

      surf_inp%vr(j) = law * surfinp(last_angle,2) + naw * surfinp(next_angle,2)
      surf_inp%hr(j) = law * surfinp(last_angle,3) + naw * surfinp(next_angle,3)
    end do

end subroutine configure
