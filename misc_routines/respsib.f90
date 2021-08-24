!
!=======================================================================
      subroutine respsib(len, nsib, nsoil, wopt, zm, &
       wsat, tgs, td, www, forcerestore,  &
       respfactor, Respg, soilscale, Moist, soilq10,  &
       test, test10)
!=======================================================================
! Calculates the rate of CO2 efflux from soils, according to the "R-star"
! approach of Denning et al (1996), adapted for use with the Bonan
! 6-layer soil thermodynamics module.  The temperature effect is a
! simple Q10, with a reference T of 298.15 (K).
!
! Modifications:
!  Scott Denning Changed soil Q10 value for soil resp from 2 to 2.4
!    following Raich and Schelsinger (1992, Tellus 44B, 81-89) (9/14/95)
!  Kevin Schaefer changed to implicit none (7/23/01)
!  Kevin Schaefer cleaned and commented the code (7/23/01)
!  Kevin Schaefer changed b to moistexp, zmstscale to moist (7/23/01)
!  Kevin Schaefer changed to moisture fraction in moistexp (7/23/01)
!  Kevin Schaefer absolute value for Rcon diagnostics (8/16/01)
!  Kevin Schaefer changed Respgtemp to JRespg, moisttemp to Jmoist (11/29/01)
!  Kevin Schaefer add canopy autotrophic resp to Rcon diagnostics (11/29/01)
!  Kevin Schaefer scaled Rcon diagnostics by Resp total (11/29/01)
!  Kevin Schaefer moved respT external so respsib works with phosib (1/21/02)
!  Kevin Schaefer Moved influence diagnostics to separate subroutine (1/21/02)
!
      implicit none
!
      integer  len   ! local number of SiB points &
      integer  nsib  ! total number of SiB points &
      integer  nsoil ! number of soil layers
!
! input/output variables
      real wopt(len)       ! (-) optimal soil moisture content for respiration
      real zm(len)         ! (-) skewness exponent for respiration exponent
      real www(len,3)      ! (-) Soil water fraction of saturation
      real wsat(len)       ! (-) determines the respiration rate at saturation
      real tgs(len)        ! (K) topsoil temperature
      real td(nsib,nsoil)  ! (K) Soil temperature
      real respfactor(nsib,nsoil+1) ! (mol/m2/s) mean unstressed respiration rate
      real Respg(len)      ! (mol/m2/s) ground respiration rate
      real soilscale(len,nsoil+1) ! (-) respiration scaling factor
      real Moist(len,2)    ! (-) soil moisture scaling factor for respiration
      real soilq10(len,nsoil+1) ! (-) soil temp scaling factor for respiration
      real MoistExp(len,2) ! (-) wetness exponent
!
! Local variables
      real woptzm          ! (-) local parameter for wetness exponent
      integer i ! SiB point index &
      integer l ! Soil layer index
      logical forcerestore
!
      real test, test10(10)   ! test variables
!
      if(.not.forcerestore) then
        do i = 1,len
!         optimal respiration wrt WWW constant
          woptzm=(wopt(i)/100.)**zm(i)
!
!         Start from the bottom and work upward
!
!         Deep soil layers 1 and 2 do not respire (no carbon below root zone)
          soilQ10(i,1)=0.
          soilscale(i,1)=0.
          soilQ10(i,2)=0.
          soilscale(i,2)=0.
!
!         root zone layers 3 to nsoil: use T=Td(3:nSoil) and W=WWW(2)
          MoistExp(i,2)=(((www(i,2))**zm(i)-woptzm)/(1.-woptzm))**2
          MoistExp(i,2)=min(MoistExp(i,2),10.)
          Moist(i,2)=0.8*wsat(i)**MoistExp(i,2)+.2
          do l = 3,nsoil
            soilQ10(i,L)=exp(0.087547*(td(i,L)-298.15))
            soilscale(i,L)=soilQ10(i,L)*Moist(i,2)
          enddo
!
!         Topsoil layer (nsoil+1): Use T=Tgs and W=WWW(1)
          MoistExp(i,1)=(((www(i,1))**zm(i)-woptzm)/(1.-woptzm))**2
          MoistExp(i,1)=min(MoistExp(i,1),10.)
          Moist(i,1)=0.8*wsat(i)**MoistExp(i,1)+.2
          soilQ10(i,nsoil+1)=exp(0.087547*(tgs(i)-298.15))
          soilscale(i,nsoil+1)=soilQ10(i,nsoil+1)*Moist(i,1)
!
!         Calculate soil resp flux to balance annual budget
          Respg(i)=respfactor(i,3)*soilscale(i,3)
          do L = 4, nSoil+1
            Respg(i)=Respg(i)+respfactor(i,L)*soilscale(i,L)
          enddo
        enddo
      else  ! (FORCERESTORE CASE ... only two soil T levels available)
        do i = 1,len
!         Moisture effect from TEM
          woptzm=(wopt(i)/100.)**zm(i)
          MoistExp(i,2)=(((www(i,2))**zm(i)-woptzm)/(1.-woptzm))**2
          MoistExp(i,2)=min(MoistExp(i,2),10.)
          Moist(i,1)=0.8*wsat(i)**MoistExp(i,2)+.2
!
!         Temperature effect is Q10 =2.4 from ref T of 25 C
          soilQ10(i,1)=exp(0.087547*(td(i,nsoil)-298.15))
          soilscale(i,1)=soilQ10(i,1)*Moist(i,1)
!
!         Dimensionalize soil resp flux to balance annual budget
          Respg(i)=respfactor(i,1)*soilscale(i,1)
        enddo
      endif
      
      test=Respg(1)
      test10(1)=Respg(1)
!
      return
      end
