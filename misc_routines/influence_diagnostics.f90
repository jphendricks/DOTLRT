!
!=======================================================================
      subroutine Influence_Diagnostics( &
      len, nsib, nsoil, wopt,Phc, &
      PhSat,bee, respfactor, Respg, Moist, &
      soilq10, Respc, RespT, scaleRsp, ScaleH2o, &
      RconWat,RconTem,RconLAI)
!=======================================================================
! Calculates the climate influence diagnostics for respiration
!
! Modifications:
!  Kevin Schaefer created subroutine (1/21/02)
!
      implicit none
!
      integer len   ! local number of SiB points
      integer nsib  ! total number of SiB points
      integer nsoil ! number of soil layers
!
! input/output variables
      real wopt(len)       ! (-) optimal soil moisture content for respiration
      real bee(len)        ! (-) soil wetness exponent
      real PhSat(len)      ! (m) Soil tension at saturation
      real Phc(len)        ! (m) 1/2 Critical leaf water potential limit
      real respfactor(nsib,nsoil+1) ! (mol/m2/s) mean unstressed respiration rate
      real Respg(len)      ! (mol/m2/s) ground respiration rate
      real Respc(len)      ! (mol/m2/s) canopy autotrophic respiration rate
      real RespT(len)      ! (mol/m2/s) Total respiration rate
      real Moist(len,2)    ! (-) soil moisture scaling factor for respiration
      real soilq10(len,nsoil+1) ! (-) soil temp scaling factor for respiration
      real RconWat(len)    ! (mol/m2/s) contribution of water to Respg
      real RconTem(len)    ! (mol/m2/s) contribution of Temperature to Respg
      real RconLAI(len)   ! (mol/m2/s) LAI effect on canopy respiration
      real ScaleRsp(len)  ! (-) autotrophic leaf respiration scaling factor
      real ScaleH2O(len)  ! (-) soil water scaling factor for photosynthesis
!
! Local variables
      real JRespg          ! (mol/m2/s) junk ground respiration rate
      real JScaleH2o       ! (-) junk H2o scaling factor autotrophic leaf resp
      integer i ! SiB point index
      integer l ! Soil layer index
!
      do i = 1,len
!       Total respiration
        RespT(i)=Respc(i)+respg(i)
!
!       ---------------------------------------------------
!       Contribution of temperature to respiration
!       ---------------------------------------------------
!       reference ground respiration rate  (T=298.15; soilQ10=1)
        JRespg=respfactor(i,3)*Moist(i,2) ! deepest layer
        do L = 4, nSoil ! rest of layers
          JRespg=JRespg+respfactor(i,L)*Moist(i,2)
        enddo
        JRespg=JRespg+respfactor(i,nSoil+1)*Moist(i,1) ! topsoil layer
!
!       temperature influence
        RconTem(i)=RespT(i)*abs(respc(i)/scaleRsp(i)+JRespg-RespT(i))
!
!       ---------------------------------------------------
!       Contribution of soil water to respiration (WWW=Wopt; Moist=1)
!       ---------------------------------------------------
!       reference autotrophic water scaling factor (WWW=Wopt)
        JScaleH2o=1./(1.+EXP(0.02*(Phc(i) &
       -PhSat(i)*Wopt(i)**(-Bee(i)))))
!
!       reference ground respiration rate  (WWW=Wopt; Moist=1)
        JRespg=respfactor(i,3)*soilQ10(i,3) ! deepest layer
        do L = 4, nSoil ! rest of layers
          JRespg=JRespg+respfactor(i,L)*soilQ10(i,L)
        enddo
        JRespg=JRespg+respfactor(i,nsoil+1)*soilQ10(i,nsoil+1) ! topsoil layer
!
!       soil water/precipitation influence
        RconWat(i)=RespT(i) &
       *abs(Respc(i)*JScaleH2o/ScaleH2o(i)-JRespg-RespT(i))
!
!       ---------------------------------------------------
!       Contribution of LAI to respiration (LAI=LAImax; aPARkk=aPARkkmax)
!       ---------------------------------------------------
!       only scaling done here, the rest done in phosib
        RconLAI(i)=RespT(i)*RconLAI(i)
      enddo
!
      return
      end
