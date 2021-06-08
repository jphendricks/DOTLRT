!=============================================================================
subroutine calc_fbw_temp_weight_scat( Tbo, tau, &
             Tb_obs, dTb_dT_obs, dTb_dp_obs, dTb_dq_obs, dTb_dw_obs, Tbo_streams, &
             dTb_dT_streams, dTb_dp_streams, dTb_dq_streams, dTb_dw_streams)
!=============================================================================
! Calculates full bandwidth brightness temperature profiles and Jacobians
! by passband quadrature over monochromatic brightness temperature profiles and Jacobians
! Calls monochromatic routine implementing DOTLRTv1.0 technique
!
! History:
!   5/24/2005 Bob Weber created routine
!   9/26/2020 Kevin Schaefer deleted unused variables
!-----------------------------------------------------------------------------
  use dotlrt_variables
  implicit none

  real(8) tb_pl_inp(0:nlev)
  real(8) tb_mn_inp(0:nlev)
  real(8) dtb_pl_inp(0:nlev,nvar), dtb_mn_inp(0:nlev,nvar)
  real(8) Tbo_streams_inp(nangover2)
  real(8) Tbo_streams(nangover2)

  integer ifreq  ! frequency index
  real(8) frequency
  real(8) Tb_inp, tau(2)

! brightness temperatures passband quadrature
  real(8), DIMENSION(0:nlev,2) :: cTb

! absorption/scattering coefficients passband quadrature
  real(8), dimension(nlev,5,3) :: cKsa

! geophysical Jacobian passband quadrature
  real(8), dimension(nlev) :: cdKab_dT
  real(8), dimension(nlev) :: cdKab_dp
  real(8), dimension(nlev) :: cdKab_dq
  real(8), dimension(nlev,5) :: cdKsc_dT
  real(8), dimension(nlev,5) :: cdg_dT
  real(8), dimension(nlev,5) :: cdKab_dw
  real(8), dimension(nlev,5) :: cdKsc_dw
  real(8), dimension(nlev,5) :: cdg_dw

! total Jacobian passband quadrature
  real(8), dimension(nlev,2) :: cdTb_dT
  real(8), dimension(nlev,2) :: cdTb_dp
  real(8), dimension(nlev,2) :: cdTb_dq
  real(8), dimension(nlev,5,2) :: cdTb_dw

! total Jacobian passband quadrature over stream angles
  real(8), dimension(nlev,nangover2,2) :: dTb_dT_streams
  real(8), dimension(nlev,nangover2,2) :: dTb_dp_streams
  real(8), dimension(nlev,nangover2,2) :: dTb_dq_streams
  real(8), dimension(nlev,nangover2,5,2) :: dTb_dw_streams

! radiative transfer Jacobian passband quadrature
  real(8), dimension(nlev,2) :: cdTbdTr     ! temperature
  real(8), dimension(nlev,2) :: cdTbdKa     ! absorption
  real(8), dimension(nlev,2) :: cdTbdKsliq  ! scatter liquid
  real(8), dimension(nlev,2) :: cdTbdgliq   ! asymmetry liquid
  real(8), dimension(nlev,2) :: cdTbdKsrn   ! scatter rain
  real(8), dimension(nlev,2) :: cdTbdgrn    ! asymmetry rain
  real(8), dimension(nlev,2) :: cdTbdKsice  ! scatter ice
  real(8), dimension(nlev,2) :: cdTbdgice   ! asymmetry ice
  real(8), dimension(nlev,2) :: cdTbdKssnow ! scatter snow
  real(8), dimension(nlev,2) :: cdTbdgsnow  ! asymmetry snow
  real(8), dimension(nlev,2) :: cdTbdKsgrpl ! scatter graupel
  real(8), dimension(nlev,2) :: cdTbdggrpl  ! asymmetry graupel

  real(8), dimension(0:nlev,2) :: Tb_obs
  real(8), dimension(nlev,2) :: dTb_dT_obs
  real(8), dimension(nlev,2) :: dTb_dp_obs
  real(8), dimension(nlev,2) :: dTb_dq_obs
  real(8), dimension(nlev,5,2) :: dTb_dw_obs
  real(8) Tbo
  real(8), dimension(nlev) :: cabs_total
  real(8), dimension(nlev) :: cscat_cloud

  integer j, k, jud
  integer iphase ! (-) hydrometeor phase index
  integer ilev ! (-) vertical level index

  real(8) quadweight, norm

! Calculate passband quadrature frequencies
! num_freqs: number of freq points for passband quadrature

  call calc_passband_freq()

! Normalization of passband quadrature weights
  if( nsub_freq == 1 ) then
    norm = num_freqs
  else
    norm = num_freqs * ( 1.0d0 - 1.0d0 / dble(nsub_freq) )
  end if

  if(num_freqs .gt. 0) then
     ! Initialize passband quadrature start
     Tbo = 0.0d0
     do j = 1, nangover2
        Tbo_streams(j) = 0.0d0
     end do
     
     do jud = 1, 2
        do ilev = 0, nlev
           ! Brightness temperatures
           cTb(ilev,jud) = 0.0d0
        end do
        do ilev = 1, nlev
           ! Geophysical/Radiation Jacobian
           cdTb_dT(ilev,jud)     = 0.0d0
           cdTb_dp(ilev,jud)     = 0.0d0
           cdTb_dq(ilev,jud)     = 0.0d0
           do iphase = 1, nphase
              cdTb_dw(ilev,iphase,jud) = 0.0d0
           end do
           ! Radiation Jacobian
           cdTbdTr(ilev,jud)     = 0.0d0
           cdTbdKa(ilev,jud)     = 0.0d0
           cdTbdKsliq(ilev,jud)  = 0.0d0
           cdTbdgliq(ilev,jud)   = 0.0d0
           cdTbdKsrn(ilev,jud)   = 0.0d0
           cdTbdgrn(ilev,jud)    = 0.0d0
           cdTbdKsice(ilev,jud)  = 0.0d0
           cdTbdgice(ilev,jud)   = 0.0d0
           cdTbdKssnow(ilev,jud) = 0.0d0
           cdTbdgsnow(ilev,jud)  = 0.0d0
           cdTbdKsgrpl(ilev,jud) = 0.0d0
           cdTbdggrpl(ilev,jud)  = 0.0d0
        end do ! ilev
     end do ! jud

     do ilev = 1, nlev
        ! Geophysical Jacobian
        cdKab_dT(ilev) = 0.0d0
        cdKab_dp(ilev) = 0.0d0
        cdKab_dq(ilev) = 0.0d0
        do iphase = 1, nphase
           cdKsc_dT(ilev,iphase) = 0.0d0
           cdg_dT(ilev,iphase)   = 0.0d0
           cdKab_dw(ilev,iphase) = 0.0d0
           cdKsc_dw(ilev,iphase) = 0.0d0
           cdg_dw(ilev,iphase)   = 0.0d0
           ! absorption/scattering coefficients:
           cKsa(ilev,iphase,1) = 0.0d0
           cKsa(ilev,iphase,2) = 0.0d0
           cKsa(ilev,iphase,3) = 0.0d0
        end do ! iphase
        cabs_total(ilev) = 0.0d0
        cscat_cloud(ilev) = 0.0d0
     end do ! ilev

     ! Initialize Jacobian over stream angles, K. Zhang 12/07/15
     do ilev = 1, nlev
        do j = 1, nangover2
           do jud = 1, 2
              dTb_dT_streams(ilev, j, jud) = 0.0d0
              dTb_dp_streams(ilev, j, jud) = 0.0d0
              dTb_dq_streams(ilev, j, jud) = 0.0d0
              do iphase = 1, nphase
                 dTb_dw_streams(ilev, j, iphase, jud) = 0.0d0
              end do
           end do
        end do
     end do
     
    ! Initialize passband quadrature end
    ! step through frequencies in channel - start
    do ifreq = 1, num_freqs
      frequency = passband_freq(ifreq)
      ! Passband quadrature over passband frequencies using trapezoid rule
      if( ( ( mod(ifreq,nsub_freq) == 0 ) .or. ( mod(ifreq,nsub_freq) == 1 ) ) &
        .and. ( nsub_freq > 1 ) ) then
        quadweight = 0.5d0
      else
        quadweight = 1.0d0
      end if
      quadweight = quadweight / norm

      ! Compute weighting vector quadrature over passband frequencies using trapezoid rule
      call calc_mon_temp_weight_scat( ifreq, Tb_inp, tau, tb_pl_inp, tb_mn_inp, dtb_pl_inp, dtb_mn_inp, Tbo_streams_inp )

      ! calculate Jacobian at streams angle, K. Zhang, 12/04/2015
      call GeoJacobian()
      call RT_Jacobian()
      call Jacobian()
      
      ! calculate Geophysical/Radiation Jacobian at streams angle, K. Zhang, 12/04/2015
      ! jud stands for upwelling or downwelling
      do ilev = 1, nlev
         do j = 1, nangover2
            do jud = 1, 2
               dTb_dT_streams(ilev,j,jud) = dTb_dT_streams(ilev,j,jud) + quadweight * dTb_dT(ilev,j,jud)
               dTb_dp_streams(ilev,j,jud) = dTb_dp_streams(ilev,j,jud) + quadweight * dTb_dp(ilev,j,jud)
               dTb_dq_streams(ilev,j,jud) = dTb_dq_streams(ilev,j,jud) + quadweight * dTb_dq(ilev,j,jud)
               do iphase = 1,nphase
                  dTb_dw_streams(ilev,j,iphase,jud) = dTb_dw_streams(ilev,j,iphase,jud) &
                                                                  + quadweight * dTb_dw(ilev,j,iphase,jud)
               end do
            end do
         end do
      end do

      ! calculate Jacobian at observation angle
      do ilev = 0, nlev
        do k = 1, nvar
          dtb_pl(ilev,1,k) = dtb_pl_inp(ilev,k)
          dtb_mn(ilev,1,k) = dtb_mn_inp(ilev,k)
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
        do ilev = 0, nlev
          ! Brightness temperatures
          cTb(ilev,1)    = cTb(ilev,1)    + quadweight * tb_pl_inp(ilev)
          cTb(ilev,2)    = cTb(ilev,2)    + quadweight * tb_mn_inp(ilev)
        end do
        do ilev = 1, nlev
          do jud = 1, 2
            ! Geophysical/Radiation Jacobian
            cdTb_dT(ilev,jud)     = cdTb_dT(ilev,jud)     + quadweight * dTb_dT(ilev,1,jud)
            cdTb_dp(ilev,jud)     = cdTb_dp(ilev,jud)     + quadweight * dTb_dp(ilev,1,jud)
            cdTb_dq(ilev,jud)     = cdTb_dq(ilev,jud)     + quadweight * dTb_dq(ilev,1,jud)
            do iphase = 1, nphase
              cdTb_dw(ilev,iphase,jud) = cdTb_dw(ilev,iphase,jud) &
                                                         + quadweight * dTb_dw(ilev,1,iphase,jud)
            end do 
            ! Radiation Jacobian
            cdTbdTr(ilev,jud)     = cdTbdTr(ilev,jud)     + quadweight * dTbdTr(ilev,1,jud)
            cdTbdKa(ilev,jud)     = cdTbdKa(ilev,jud)     + quadweight * dTbdKa(ilev,1,jud)
            cdTbdKsliq(ilev,jud)  = cdTbdKsliq(ilev,jud)  + quadweight * dTbdKsliq(ilev,1,jud)
            cdTbdgliq(ilev,jud)   = cdTbdgliq(ilev,jud)   + quadweight * dTbdgliq(ilev,1,jud)
            cdTbdKsrn(ilev,jud)   = cdTbdKsrn(ilev,jud)   + quadweight * dTbdKsrn(ilev,1,jud)
            cdTbdgrn(ilev,jud)    = cdTbdgrn(ilev,jud)    + quadweight * dTbdgrn(ilev,1,jud)
            cdTbdKsice(ilev,jud)  = cdTbdKsice(ilev,jud)  + quadweight * dTbdKsice(ilev,1,jud)
            cdTbdgice(ilev,jud)   = cdTbdgice(ilev,jud)   + quadweight * dTbdgice(ilev,1,jud)
            cdTbdKssnow(ilev,jud) = cdTbdKssnow(ilev,jud) + quadweight * dTbdKssnow(ilev,1,jud)
            cdTbdgsnow(ilev,jud)  = cdTbdgsnow(ilev,jud)  + quadweight * dTbdgsnow(ilev,1,jud)
            cdTbdKsgrpl(ilev,jud) = cdTbdKsgrpl(ilev,jud) + quadweight * dTbdKsgrpl(ilev,1,jud)
            cdTbdggrpl(ilev,jud)  = cdTbdggrpl(ilev,jud)  + quadweight * dTbdggrpl(ilev,1,jud)
          end do ! jud
        end do ! ilev
      do ilev = 1, nlev
        ! Geophysical Jacobian
        cdKab_dT(ilev) = cdKab_dT(ilev) + quadweight * dKab_dT(ilev)
        cdKab_dp(ilev) = cdKab_dp(ilev) + quadweight * dKab_dp(ilev)
        cdKab_dq(ilev) = cdKab_dq(ilev) + quadweight * dKab_dq(ilev)
        do iphase = 1, nphase
          cdKsc_dT(ilev,iphase) = cdKsc_dT(ilev,iphase) + quadweight * dKsc_dT(ilev,iphase)
          cdg_dT(ilev,iphase)   = cdg_dT(ilev,iphase)   + quadweight * dg_dT(ilev,iphase)
          cdKab_dw(ilev,iphase) = cdKab_dw(ilev,iphase) + quadweight * dKab_dw(ilev,iphase)
          cdKsc_dw(ilev,iphase) = cdKsc_dw(ilev,iphase) + quadweight * dKsc_dw(ilev,iphase)
          cdg_dw(ilev,iphase)   = cdg_dw(ilev,iphase)   + quadweight * dg_dw(ilev,iphase)

        ! absorption/scattering coefficients:
          cKsa(ilev,iphase,1) = cKsa(ilev,iphase,1) + quadweight * hydro_prof(ilev,iphase)%cloudab
          cKsa(ilev,iphase,2) = cKsa(ilev,iphase,2) + quadweight * hydro_prof(ilev,iphase)%cloudsc
          cKsa(ilev,iphase,3) = cKsa(ilev,iphase,3) + quadweight * hydro_prof(ilev,iphase)%cloudg
        end do ! iphase
        cabs_total(ilev) = cabs_total(ilev) + quadweight * abs_total1(ilev)
        cscat_cloud(ilev) = cscat_cloud(ilev) + quadweight * scat_cloud1(ilev)
      end do ! ilev
    end do ! ifreq

    ! step through frequencies in channel - end
    ! brigthness temperature and Jacobian profiles for stream angles - start
    do ilev = 0, nlev
       Tb_obs(ilev,1) = cTb(ilev,1)
       Tb_obs(ilev,2) = cTb(ilev,2)
    end do
    do ilev = 1, nlev
       dTb_dT_obs(ilev,1) = cdTb_dT(ilev,1)
       dTb_dp_obs(ilev,1) = cdTb_dp(ilev,1)
       dTb_dq_obs(ilev,1) = cdTb_dq(ilev,1)
       dTb_dT_obs(ilev,2) = cdTb_dT(ilev,2)
       dTb_dp_obs(ilev,2) = cdTb_dp(ilev,2)
       dTb_dq_obs(ilev,2) = cdTb_dq(ilev,2)
       do iphase = 1, nphase
          dTb_dw_obs(ilev,iphase,1) = cdTb_dw(ilev,iphase,1)
          dTb_dw_obs(ilev,iphase,2) = cdTb_dw(ilev,iphase,2)
       end do
    end do ! ilev
    ! brightness temperature and Jacobian profiles for stream angles - end
  end if ! num_freqs

end subroutine calc_fbw_temp_weight_scat
