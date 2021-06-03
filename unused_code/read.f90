!=====================================================================
subroutine readsnowtable
!=====================================================================
! Purpose:
!   This subroutine reads in the snow type lookup table                                                     *
!
! Modifications:
!   Kevin Schaefer created routine (8/21/07)
! Note: this is not the exact version in mapper
!---------------------------------------------------------------------
!
use mapper_variables
!
implicit none
!
character*100 junk  ! junk variable to read character text from table
integer i  ! indeces
!
! read snow type table
    open(unit=55,file=trim(filename%snowtab),form='formatted', status='old')
    read(55,20) junk
    read(55,20) junk
    do i=1,nsnowclass
        read(55,*) snowtab(i)%snow_class,   &
            snowtab(i)%name,   &
            snowtab(i)%f_bot_max, &                    
            snowtab(i)%d_min,   &
            snowtab(i)%d_max,   &
            snowtab(i)%den_ref_top,   &
            snowtab(i)%den_ref_bot,   &
            snowtab(i)%den_bulk_obs,   &
            snowtab(i)%den_bulk_std  
    enddo
!
! generic formats
20    format (a100)
!
return
end subroutine readsnowtable
!
!=======================================================================
      subroutine ReadBioTable
!=======================================================================       
!     read in vegetation parameter lookup table                                
!-----------------------------------------------------------------------
!
      use Mapper_Variables
!      
      implicit none
!
      character*100 junk  ! junk variable to read character text from table
      integer iv        ! counter for biome type do-loop
      integer ivv       ! temporary variable for reading vegetation type

!                                                                               
      open(unit=55,file=FileName%BioTab,form='formatted', status='old')
!
      read(55,20) junk
      do iv=1,nv
         read(55,13) ivv, &
                     BioTab(iv)%z2, &
                     BioTab(iv)%z1, &
                     BioTab(iv)%fVcover, &
                     BioTab(iv)%Chil, &
                     BioTab(iv)%SoDep, &
                     BioTab(iv)%RootD, &
                     BioTab(iv)%Phi_half
!
!       assign biome type to BioTab
        BioTab(iv)%BioNum=ivv
      enddo

      read(55,20) junk
      read(55,20) junk
      do iv=1,nv
         read(55,14) ivv, &
                     BioTab(iv)%LTran(1,1), &
                     BioTab(iv)%LTran(2,1), &
                     BioTab(iv)%LTran(1,2), &
                     BioTab(iv)%LTran(2,2), &
                     BioTab(iv)%LRef(1,1), &
                     BioTab(iv)%LRef(2,1), &
                     BioTab(iv)%LRef(1,2), &
                     BioTab(iv)%LRef(2,2)
      enddo
       read(55,20) junk
      do iv=1,nv
        read(55,12) ivv, &
                     BioTab(iv)%vmax0,  &
                     BioTab(iv)%EffCon,   &
                     BioTab(iv)%gsSlope,  &
                     BioTab(iv)%gsMin, &
                     BioTab(iv)%Atheta,  &
                     BioTab(iv)%Btheta
      enddo
      read(55,20) junk
      do iv=1,nv
         read(55,14) ivv, &
                     BioTab(iv)%TRDA,    &
                     BioTab(iv)%TRDM,    &
                     BioTab(iv)%TROP,  &
                     BioTab(iv)%respcp,                       &
                     BioTab(iv)%SLTI,    &
                     BioTab(iv)%HLTI,    &
                     BioTab(iv)%SHTI,   &
                     BioTab(iv)%HHTI
      enddo
      read(55,20) junk
      do iv=1,nv
         read(55,14) ivv, &
                     BioTab(iv)%SoRef(1),    &
                     BioTab(iv)%SoRef(2)
      enddo
      close(55)
!
! assign other biome parameters
      do iv=1,nv
        BioTab(iv)%tfrost=283.15
      enddo
      BioTab(10)%tfrost=274.15
      
20    Format (a100)
12    format (i2,4x,e9.3,5(f8.3,1x))
13    format (i2,2x,6(f8.3,1x),f10.3)
14    format (i2,2x,8(f8.3,1x))
!
      return
      end
!
!======================================================================
      subroutine ReadMorphTable
!======================================================================
!     read in morphological and bounding ndvi values for each vegtype
!-----------------------------------------------------------------------
!
      use Mapper_Variables
!
      implicit none
!
      integer iv       ! biome type do-loop counter
      integer ivv      ! biome type number
!
      open(unit=7,file=FileName%MorphTab,form='formatted', status='old')
!
      do iv=1, nv
!
      read(7,*) ivv, &
       MorphTab(iv)%zc, &
       MorphTab(iv)%LWidth, &
       MorphTab(iv)%LLength, &
       MorphTab(iv)%LAImax, &
       MorphTab(iv)%stems, &
       MorphTab(iv)%NDVImax,  &
       MorphTab(iv)%NDVImin
!
! Convert maximum/minimum NDVI values to simple ratios 
       MorphTab(iv)%SRmax=(1.+MorphTab(iv)%NDVImax)/(1.-MorphTab(iv)%NDVImax)
       MorphTab(iv)%SRmin=(1.+MorphTab(iv)%NDVImin)/(1.-MorphTab(iv)%NDVImin)
!
      enddo
      close(7)
      return
      end
!
!======================================================================
      subroutine ReadSoilTable
!======================================================================
!     read in soil properties for each texture class
!----------------------------------------------------------------------
!
      use Mapper_Variables
!
      implicit none
!
      character*100 junk  ! junk variable to read character text from table
      integer s    ! soil type index
!
      open (unit=14,file=trim(Filename%SoilTab),form='formatted', status='old')
!
      read(14,20) junk
      do s=1,ns
        read(14,30) SoilTab(s)%SoilNum, &
        SoilTab(s)%BEE, &
        SoilTab(s)%PhiSat, &
        SoilTab(s)%SatCo, &
        SoilTab(s)%poros, &
        SoilTab(s)%Slope, &
        SoilTab(s)%Wopt, &
        SoilTab(s)%Skew, &
        SoilTab(s)%RespSat
      enddo
!
      close(unit=14)
!
20    Format (a100)
30    format (i2,1x,2(f7.3,1x),e9.3, 5(f7.3,1x))
!
      return
      end
!
!=======================================================================        
      subroutine ReadAeroTables
!=======================================================================
! This subroutine reads in interpolation tables of previously 
! calculated aerodynamic variables. 
!
      use Mapper_Variables 
!
      implicit none
!
! open interpolation tables
      open(unit=7,file=Filename%AeroVar, form='unformatted', status='old')
!
!   Read in interpolation tables
      read (7) LAIgrid
      read (7) fVCovergrid
      read (7) AeroVar%zo
      read (7) AeroVar%zp_disp
      read (7) AeroVar%RbC
      read (7) AeroVar%RdC
!
! close interpolation tables
      Close (7)
!
      return
      end
