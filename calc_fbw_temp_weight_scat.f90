! calc_fbw_temp_weight_scat.f90
! 24 May 2005 Bob Weber
! Calculates full bandwidth brightness temperature profiles and Jacobians
!   by passband quadrature over monochromatic brightness temperature profiles and Jacobians
! Calls monochromatic routine implementing DOTLRTv1.0 technique

subroutine calc_fbw_temp_weight_scat( Tbo, tau, &
             Tb_obs, dTb_dT_obs, dTb_dp_obs, dTb_dq_obs, dTb_dw_obs, Tbo_streams, &
             dTb_dT_streams, dTb_dp_streams, dTb_dq_streams, dTb_dw_streams)
  use variables
  implicit none

  real(8) tb_pl_inp(0:atm_inp%num_levels), tb_mn_inp(0:atm_inp%num_levels)
  real(8) dtb_pl_inp(0:atm_inp%num_levels,nvar), dtb_mn_inp(0:atm_inp%num_levels,nvar)
  real(8) Tbo_streams_inp(nangover2)
  real(8) Tbo_streams(nangover2)

  integer jchan
  real(8) frequency

  integer alloc_err, last_level, next_level
  real(8) Tb_inp, tau(2)

  ! brightness temperatures passband quadrature
  real(8), DIMENSION(0:atm_inp%num_levels,2) :: cTb

  ! absorption/scattering coefficients passband quadrature
  real(8), dimension(atm_inp%num_levels,5,3) :: cKsa

  ! geophysical Jacobian passband quadrature
  real(8), dimension(atm_inp%num_levels) :: cdKab_dT, cdKab_dp, cdKab_dq
  real(8), dimension(atm_inp%num_levels,5) :: cdKsc_dT, cdg_dT, cdKab_dw, cdKsc_dw, cdg_dw

  ! total Jacobian passband quadrature
  real(8), dimension(atm_inp%num_levels,2) :: cdTb_dT, cdTb_dp, cdTb_dq
  real(8), dimension(atm_inp%num_levels,5,2) :: cdTb_dw

  ! total Jacobian passband quadrature over stream angles
  real(8), dimension(atm_inp%num_levels,nangover2,2) :: dTb_dT_streams, dTb_dp_streams, dTb_dq_streams
  real(8), dimension(atm_inp%num_levels,nangover2,5,2) :: dTb_dw_streams

  ! radiative transfer Jacobian passband quadrature
  real(8), dimension(atm_inp%num_levels,2) :: &
                                           cdTbdTr,      & ! temperature
                                           cdTbdKa,      & ! absorption
                                           cdTbdKsliq,   & ! scatter liquid
                                           cdTbdgliq,    & ! asymmetry liquid
                                           cdTbdKsrn,    & ! scatter rain
                                           cdTbdgrn,     & ! asymmetry rain
                                           cdTbdKsice,   & ! scatter ice
                                           cdTbdgice,    & ! asymmetry ice
                                           cdTbdKssnow,  & ! scatter snow
                                           cdTbdgsnow,   & ! asymmetry snow
                                           cdTbdKsgrpl,  & ! scatter graupel
                                           cdTbdggrpl      ! asymmetry graupel

  real(8), dimension(0:atm_inp%num_levels,2) :: Tb_obs
  real(8), dimension(atm_inp%num_levels,2) :: dTb_dT_obs
  real(8), dimension(atm_inp%num_levels,2) :: dTb_dp_obs
  real(8), dimension(atm_inp%num_levels,2) :: dTb_dq_obs
  real(8), dimension(atm_inp%num_levels,5,2) :: dTb_dw_obs
  real(8) Tbo
  real(8), dimension(atm_inp%num_levels) ::  cabs_total, cscat_cloud

  integer j, ilr, ilr1, k, i, hydrometeor_phase, jud, kpts, kind

  real(8) quadweight, norm, theta_obs

  DOUBLE PRECISION cos_observ_angle

  character*120 debugout

  ! Calculate passband quadrature frequencies
  ! num_freqs: number of freq points for passband quadrature
  jchan = 1
  call calc_passband_freq(jchan)
  ! Normalization of passband quadrature weights
  if( num_sb_freqs == 1 ) then
    norm = num_freqs
  else
    norm = num_freqs * ( 1.0d0 - 1.0d0 / dble(num_sb_freqs) )
  end if

  if(num_freqs .gt. 0) then
     ! Initialize passband quadrature start !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     Tbo = 0.0d0
     do j = 1, nangover2
        Tbo_streams(j) = 0.0d0
     end do
     
     do jud = 1, 2
        do ilr1 = 0, nlr1
           ! Brightness temperatures
           cTb(ilr1,jud) = 0.0d0
        end do
        do ilr1 = 1, nlr1
           ! Geophysical/Radiation Jacobian
           cdTb_dT(ilr1,jud)     = 0.0d0
           cdTb_dp(ilr1,jud)     = 0.0d0
           cdTb_dq(ilr1,jud)     = 0.0d0
           do hydrometeor_phase = 1, number_h2o_phases
              cdTb_dw(ilr1,hydrometeor_phase,jud) = 0.0d0
           end do
           ! Radiation Jacobian
           cdTbdTr(ilr1,jud)     = 0.0d0
           cdTbdKa(ilr1,jud)     = 0.0d0
           cdTbdKsliq(ilr1,jud)  = 0.0d0
           cdTbdgliq(ilr1,jud)   = 0.0d0
           cdTbdKsrn(ilr1,jud)   = 0.0d0
           cdTbdgrn(ilr1,jud)    = 0.0d0
           cdTbdKsice(ilr1,jud)  = 0.0d0
           cdTbdgice(ilr1,jud)   = 0.0d0
           cdTbdKssnow(ilr1,jud) = 0.0d0
           cdTbdgsnow(ilr1,jud)  = 0.0d0
           cdTbdKsgrpl(ilr1,jud) = 0.0d0
           cdTbdggrpl(ilr1,jud)  = 0.0d0
        end do ! ilr1
     end do ! jud

     do ilr1 = 1, nlr1
        ! Geophysical Jacobian
        cdKab_dT(ilr1) = 0.0d0
        cdKab_dp(ilr1) = 0.0d0
        cdKab_dq(ilr1) = 0.0d0
        do hydrometeor_phase = 1, number_h2o_phases
           cdKsc_dT(ilr1,hydrometeor_phase) = 0.0d0
           cdg_dT(ilr1,hydrometeor_phase)   = 0.0d0
           cdKab_dw(ilr1,hydrometeor_phase) = 0.0d0
           cdKsc_dw(ilr1,hydrometeor_phase) = 0.0d0
           cdg_dw(ilr1,hydrometeor_phase)   = 0.0d0
           ! absorption/scattering coefficients:
           cKsa(ilr1,hydrometeor_phase,1) = 0.0d0
           cKsa(ilr1,hydrometeor_phase,2) = 0.0d0
           cKsa(ilr1,hydrometeor_phase,3) = 0.0d0
        end do ! hydrometeor_phase
        cabs_total(ilr1) = 0.0d0
        cscat_cloud(ilr1) = 0.0d0
     end do ! ilr1

     ! Initialize Jacobian over stream angles, K. Zhang 12/07/15
     do ilr1 = 1, nlr1
        do j = 1, nangover2
           do jud = 1, 2
              dTb_dT_streams(ilr1, j, jud) = 0.0d0
              dTb_dp_streams(ilr1, j, jud) = 0.0d0
              dTb_dq_streams(ilr1, j, jud) = 0.0d0
              do hydrometeor_phase = 1, number_h2o_phases
                 dTb_dw_streams(ilr1, j, hydrometeor_phase, jud) = 0.0d0
              end do
           end do
        end do
     end do
     
    ! Initialize passband quadrature end   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! step through frequencies in channel - start !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do nfreq = 1, num_freqs
      frequency = passband_freq(nfreq)
      ! Passband quadrature over passband frequencies using trapezoid rule
      if( ( ( mod(nfreq,num_sb_freqs) == 0 ) .or. ( mod(nfreq,num_sb_freqs) == 1 ) ) &
        .and. ( num_sb_freqs > 1 ) ) then
        quadweight = 0.5d0
      else
        quadweight = 1.0d0
      end if
      quadweight = quadweight / norm

      !write(debugout,*) "quadweight = ",quadweight
      !call mexPrintf(debugout//achar(10))

      ! Compute weighting vector quadrature over passband frequencies using trapezoid rule
      call calc_mon_temp_weight_scat( Tb_inp, tau, tb_pl_inp, tb_mn_inp, dtb_pl_inp, dtb_mn_inp, Tbo_streams_inp )

      ! calculate Jacobian at streams angle, K. Zhang, 12/04/2015
      call GeoJacobian()
      call RT_Jacobian()
      call Jacobian()
      
      ! calculate Geophysical/Radiation Jacobian at streams angle, K. Zhang, 12/04/2015
      ! jud stands for upwelling or downwelling
      do ilr1 = 1, nlr1
         do j = 1, nangover2
            do jud = 1, 2
               dTb_dT_streams(ilr1,j,jud) = dTb_dT_streams(ilr1,j,jud) + quadweight * dTb_dT(ilr1,j,jud)
               dTb_dp_streams(ilr1,j,jud) = dTb_dp_streams(ilr1,j,jud) + quadweight * dTb_dp(ilr1,j,jud)
               dTb_dq_streams(ilr1,j,jud) = dTb_dq_streams(ilr1,j,jud) + quadweight * dTb_dq(ilr1,j,jud)
               do hydrometeor_phase = 1,number_h2o_phases
                  dTb_dw_streams(ilr1,j,hydrometeor_phase,jud) = dTb_dw_streams(ilr1,j,hydrometeor_phase,jud) &
                                                                  + quadweight * dTb_dw(ilr1,j,hydrometeor_phase,jud)
               end do
            end do
         end do
      end do

      ! calculate Jacobian at observation angle
      do ilr1 = 0, nlr1
        do k = 1, nvar
          dtb_pl(ilr1,1,k) = dtb_pl_inp(ilr1,k)
          dtb_mn(ilr1,1,k) = dtb_mn_inp(ilr1,k)
        end do
      end do
      call GeoJacobian()
      call RT_Jacobian()
      call Jacobian()

      ! passband quadrature of Tb streams at observation level
      do j = 1, nangover2
         Tbo_streams(j) = Tbo_streams(j) + quadweight * Tbo_streams_inp(j)
      end do

      ! passband quadrature
      Tbo = Tbo + quadweight * Tb_inp
        do ilr1 = 0, nlr1
          ! Brightness temperatures
          cTb(ilr1,1)    = cTb(ilr1,1)    + quadweight * tb_pl_inp(ilr1)
          cTb(ilr1,2)    = cTb(ilr1,2)    + quadweight * tb_mn_inp(ilr1)
        end do
        do ilr1 = 1, nlr1
          do jud = 1, 2
            ! Geophysical/Radiation Jacobian
            cdTb_dT(ilr1,jud)     = cdTb_dT(ilr1,jud)     + quadweight * dTb_dT(ilr1,1,jud)
            cdTb_dp(ilr1,jud)     = cdTb_dp(ilr1,jud)     + quadweight * dTb_dp(ilr1,1,jud)
            cdTb_dq(ilr1,jud)     = cdTb_dq(ilr1,jud)     + quadweight * dTb_dq(ilr1,1,jud)
            do hydrometeor_phase = 1, number_h2o_phases
              cdTb_dw(ilr1,hydrometeor_phase,jud) = cdTb_dw(ilr1,hydrometeor_phase,jud) &
                                                         + quadweight * dTb_dw(ilr1,1,hydrometeor_phase,jud)
            end do 
            ! Radiation Jacobian
            cdTbdTr(ilr1,jud)     = cdTbdTr(ilr1,jud)     + quadweight * dTbdTr(ilr1,1,jud)
            cdTbdKa(ilr1,jud)     = cdTbdKa(ilr1,jud)     + quadweight * dTbdKa(ilr1,1,jud)
            cdTbdKsliq(ilr1,jud)  = cdTbdKsliq(ilr1,jud)  + quadweight * dTbdKsliq(ilr1,1,jud)
            cdTbdgliq(ilr1,jud)   = cdTbdgliq(ilr1,jud)   + quadweight * dTbdgliq(ilr1,1,jud)
            cdTbdKsrn(ilr1,jud)   = cdTbdKsrn(ilr1,jud)   + quadweight * dTbdKsrn(ilr1,1,jud)
            cdTbdgrn(ilr1,jud)    = cdTbdgrn(ilr1,jud)    + quadweight * dTbdgrn(ilr1,1,jud)
            cdTbdKsice(ilr1,jud)  = cdTbdKsice(ilr1,jud)  + quadweight * dTbdKsice(ilr1,1,jud)
            cdTbdgice(ilr1,jud)   = cdTbdgice(ilr1,jud)   + quadweight * dTbdgice(ilr1,1,jud)
            cdTbdKssnow(ilr1,jud) = cdTbdKssnow(ilr1,jud) + quadweight * dTbdKssnow(ilr1,1,jud)
            cdTbdgsnow(ilr1,jud)  = cdTbdgsnow(ilr1,jud)  + quadweight * dTbdgsnow(ilr1,1,jud)
            cdTbdKsgrpl(ilr1,jud) = cdTbdKsgrpl(ilr1,jud) + quadweight * dTbdKsgrpl(ilr1,1,jud)
            cdTbdggrpl(ilr1,jud)  = cdTbdggrpl(ilr1,jud)  + quadweight * dTbdggrpl(ilr1,1,jud)
          end do ! jud
        end do ! ilr1
      do ilr1 = 1, nlr1
        ! Geophysical Jacobian
        cdKab_dT(ilr1) = cdKab_dT(ilr1) + quadweight * dKab_dT(ilr1)
        cdKab_dp(ilr1) = cdKab_dp(ilr1) + quadweight * dKab_dp(ilr1)
        cdKab_dq(ilr1) = cdKab_dq(ilr1) + quadweight * dKab_dq(ilr1)
        do hydrometeor_phase = 1, number_h2o_phases
          cdKsc_dT(ilr1,hydrometeor_phase) = cdKsc_dT(ilr1,hydrometeor_phase) &
                                           + quadweight * dKsc_dT(ilr1,hydrometeor_phase)
          cdg_dT(ilr1,hydrometeor_phase)   = cdg_dT(ilr1,hydrometeor_phase)   &
                                           + quadweight * dg_dT(ilr1,hydrometeor_phase)
          cdKab_dw(ilr1,hydrometeor_phase) = cdKab_dw(ilr1,hydrometeor_phase) &
                                           + quadweight * dKab_dw(ilr1,hydrometeor_phase)
          cdKsc_dw(ilr1,hydrometeor_phase) = cdKsc_dw(ilr1,hydrometeor_phase) &
                                           + quadweight * dKsc_dw(ilr1,hydrometeor_phase)
          cdg_dw(ilr1,hydrometeor_phase)   = cdg_dw(ilr1,hydrometeor_phase)   &
                                           + quadweight * dg_dw(ilr1,hydrometeor_phase)
        ! absorption/scattering coefficients:
          cKsa(ilr1,hydrometeor_phase,1) = cKsa(ilr1,hydrometeor_phase,1) &
                                         + quadweight * hydro_prof(ilr1,hydrometeor_phase)%cloudab
          cKsa(ilr1,hydrometeor_phase,2) = cKsa(ilr1,hydrometeor_phase,2) &
                                         + quadweight * hydro_prof(ilr1,hydrometeor_phase)%cloudsc
          cKsa(ilr1,hydrometeor_phase,3) = cKsa(ilr1,hydrometeor_phase,3) &
                                         + quadweight * hydro_prof(ilr1,hydrometeor_phase)%cloudg
        end do ! hydrometeor_phase
        cabs_total(ilr1) = cabs_total(ilr1) &
                         + quadweight * abs_total1(ilr1)
        cscat_cloud(ilr1) = cscat_cloud(ilr1) &
                          + quadweight * scat_cloud1(ilr1)
      end do ! ilr1
    end do ! nfreq

    ! step through frequencies in channel - end   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! brigthness temperature and Jacobian profiles for stream angles - start
    do ilr1 = 0, nlr1
       Tb_obs(ilr1,1) = cTb(ilr1,1)
       Tb_obs(ilr1,2) = cTb(ilr1,2)
    end do
    do ilr1 = 1, nlr1
       dTb_dT_obs(ilr1,1) = cdTb_dT(ilr1,1)
       dTb_dp_obs(ilr1,1) = cdTb_dp(ilr1,1)
       dTb_dq_obs(ilr1,1) = cdTb_dq(ilr1,1)
       dTb_dT_obs(ilr1,2) = cdTb_dT(ilr1,2)
       dTb_dp_obs(ilr1,2) = cdTb_dp(ilr1,2)
       dTb_dq_obs(ilr1,2) = cdTb_dq(ilr1,2)
       do hydrometeor_phase = 1, number_h2o_phases
          dTb_dw_obs(ilr1,hydrometeor_phase,1) = cdTb_dw(ilr1,hydrometeor_phase,1)
          dTb_dw_obs(ilr1,hydrometeor_phase,2) = cdTb_dw(ilr1,hydrometeor_phase,2)
       end do
    end do ! ilr1
    ! brigthness temperature and Jacobian profiles for stream angles - end
  end if ! num_freqs
end subroutine calc_fbw_temp_weight_scat
