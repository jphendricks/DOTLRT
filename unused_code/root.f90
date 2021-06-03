!
!=======================================================================
      SUBROUTINE RootFraction(len,nsoil,bionum,xroot,test,test10)
!=======================================================================
! Calculate the fraction of roots per layer given a standardized
! vertical root distribution
!
!  REFERENCES:
!  Denning et al., 1996a, Tellus 48B, 521-542
!  Denning at al., 1996b, Tellus 48B, 543-567 
!
! Modifications:
!  Kevin Schaefer converted this into a separate subroutine (1/18/01)
!  Kevin Schaefer converted into a SiB compatable subroutine (8/2/01)
!
! Input variables
      INTEGER bionum(len)
      INTEGER len
      INTEGER nsoil
!
! Output variable
      REAL xroot(len,nsoil+1)  ! map of root Fraction per soil layer
!
! Local variables
      REAL kRoot        ! Exp const for decrease root-zone organic matter w/ depth
      REAL kRootSiB(12) ! Lookup table for kRoot by biome  
      REAL dZ(nsoil+1)        ! Relative depth of each active soil layer
      REAL totalRoot    ! Sum of root-zone organic matter
      REAL ztop         ! top of soil layer
      REAL zbot         ! bottom of soil layer
      REAL temp(nsoil+1)      ! temporary storage variable for root fraction
      real test, test10(10)   ! test variables
!
      kRootSiB(1)=3.9
      kRootSiB(2)=3.9
      kRootSiB(3)=2.0
      kRootSiB(4)=5.5
      kRootSiB(5)=5.5
      kRootSiB(6)=2.0
      kRootSiB(7)=5.5
      kRootSiB(8)=2.0
      kRootSiB(9)=2.0
      kRootSiB(10)=5.5
      kRootSiB(11)=2.0
      kRootSiB(12)=5.5
!
! Calculate relative depth of each layer in root zone
      dZ(1)=0.          ! surface layer (no roots)
      dZ(2)=1./15.      ! bottom soil layer 6  (top of root zone)
      dZ(3)=2./15.      ! bottom soil layer 5  (root zone)
      dZ(4)=4./15.      ! bottom soil layer 4  (root zone)
      dZ(5)=8./15.      ! bottom soil layer 3  (bottom root zone)
      dZ(6)=0.          ! bottom soil layer 2  (abiotic recharge zone)
      dZ(7)=0.          ! bottom soil layer 1  (abiotic recharge zone)
!
! Scan through vegetation map
      xRoot=0.
      do i=1,len
!
!         Calculate total amount of roots
          kRoot=kRootSiB(bionum(i))
          totalRoot=(1.-exp(-kRoot))/kRoot
!
!         Calculate root fraction in each layer
          ztop=0.
          Do L=1,5
            zbot=ztop+dZ(L)
            temp(L)=(exp(-kRoot*ztop)-exp(-kRoot*zbot)) &
            /(kRoot*totalRoot)
            ztop=zbot
          End do
!
!         Switch order because soilscale level is opposite to root density
!         xroot(i,7)=0.
          xroot(i,6)=temp(2)
          xroot(i,5)=temp(3)
          xroot(i,4)=temp(4)
          xroot(i,3)=temp(5)
!         xroot(i,2)=0.
!         xroot(i,1)=0.
!
      end do
      
      test=dZ(1)
      test10(1)=dZ(1)
!
      RETURN
      END
!
!=======================================================================
      SUBROUTINE permacarb(d_active,test,test10)
!=======================================================================
! Calculate the fraction of roots per layer given a standardized
! vertical root distribution
!
!  REFERENCES:
!  Denning et al., 1996a, Tellus 48B, 521-542
!  Denning at al., 1996b, Tellus 48B, 543-567 
!
! Modifications:
!  Kevin Schaefer converted this into a separate subroutine (1/18/01)
!  Kevin Schaefer converted into a SiB compatable subroutine (8/2/01)
!
! Input variables
    real d_active                     ! active layer thickness
    real test, test10(10)   ! test variables
!
! Local variables
    real del_active ! (m) change in maximum active layer depth
    real del_carb   ! (mol/m2) change in soil carbon content
    real d_permc_min ! (m) top of permacarb layer
    real d_permc_max ! (m) bottom of permacarb layer
!  
! initialize variables
    del_active=0.
    del_carb=0.
    d_permc_min=1.5
    d_permc_max=3.
!
! test if permafrost carbon is thawed
    if(d_active>d_permc_min.and.d_permc_min<d_permc_max) then
!
! change in threshold depth
          del_active=d_active-d_permc_min
          if(d_active>d_permc_max) then
            del_active=d_permc_max-d_permc_min
          endif
          d_permc_min=min(d_active,d_permc_max)
    endif
!
! Output variable
      test = d_active
      test10(1)=d_active
      test10(2)=d_permc_min
      test10(3)=del_active
!
      RETURN
      END
