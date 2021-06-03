!
!=======================================================================
      subroutine snow_test_setup
!=========================================================================
! reads input parameters for the scan program
!
! Modifications:
!  Kevin Schaefer split off from main program (3/34/04)
!--------------------------------------------------------------------------
!
      use Mapper_Variables
      use sibscan_Variables
!
      IMPLICIT NONE
!
integer m,n  ! indeces
!      
! new snow thickness criteria
    sib(1)%prog%dz_min(1)=0.02
    sib(1)%prog%dz_max(1)=0.03
    sib(1)%prog%dz_new(1)=sib(1)%prog%dz_max(1)+1.6*sib(1)%prog%dz_min(1)
    sib(1)%prog%snow_min_mass(1)=0.05
    do m=2,nsnow
      sib(1)%prog%dz_min(m)=1.6*sib(1)%prog%dz_min(m-1)
      sib(1)%prog%dz_max(m)=1.6*sib(1)%prog%dz_max(m-1)
      sib(1)%prog%dz_new(m)=sib(1)%prog%dz_max(m)+1.6*sib(1)%prog%dz_min(m)
      sib(1)%prog%snow_min_mass(m)=0.05
    enddo
    sib(1)%prog%dz_max(nsnow)=10000.
    sib(1)%prog%dz_new(nsnow)=10000.
!
! snowpack properties
    sib(1)%prog%snow_bulk=sib(1)%prog%snow_den(1)
    sib(1)%prog%dz(1+sib%prog%nsl)=sib(1)%prog%snow_depth
    sib(1)%prog%z_bot(1+sib%prog%nsl)=0.
    sib(1)%prog%z_top(1+sib%prog%nsl)=sib(1)%prog%snow_depth
    sib(1)%prog%node_z(1+sib%prog%nsl)=-sib(1)%prog%snow_depth/2.
    sib(1)%prog%www_ice(1+sib%prog%nsl)=sib(1)%prog%snow_den(1)*sib(1)%prog%snow_depth
!
! snow table variables
    n=snowclass
    sib%param%snow_class=snowclass
    sib(1)%param%f_bot_max=snowtab(n)%f_bot_max                   
    sib(1)%param%d_min=snowtab(n)%d_min
    sib(1)%param%d_max=snowtab(n)%d_max
    sib(1)%param%den_ref_top=snowtab(n)%den_ref_top
    sib(1)%param%den_ref_bot=snowtab(n)%den_ref_bot
    sib(1)%param%den_bulk_obs=snowtab(n)%den_bulk_obs
    sib(1)%param%den_bulk_std=snowtab(n)%den_bulk_std  
!
      return
      end
!
!======================================================================
subroutine test_snow_density(sib)
!======================================================================
! Description:
!   changes snow layer thickness based on standard density profiles
!
! Modifications:
!   Kevin Schaefer created routine (7/16/07)
!----------------------------------------------------------------------

use kinds
use sibtype

implicit none
!
! i/o variables
type(sib_t), intent(inout) :: sib
!
! local variables
integer(kind=int_kind) :: i,n  ! indeces
integer(kind=int_kind) :: num_snow   ! number of snow layers
real(kind=dbl_kind) :: den_slope     ! (1/m) slope of top layer density function
real(kind=dbl_kind) :: d_slope       ! (1/m) slope of bottom layer depth fraction function
real(kind=dbl_kind) :: d_half        ! (m) half point of bottom layer depth fraction
real(kind=dbl_kind) :: delta_den     ! (kg/m3) difference between top and bottom layer densities
real(kind=dbl_kind) :: snow_den_bot  ! (kg/m3) bottom layer density
real(kind=dbl_kind) :: snow_den_top  ! (kg/m3) top layer density
real(kind=dbl_kind) :: s_delta_den   ! (-) bulk density scaling factor
real(kind=dbl_kind) :: f_bot         ! (m) bottom layer depth fraction
real(kind=dbl_kind) :: d_bot         ! (m) bottom layer depth
real(kind=dbl_kind) :: dzsno(nsnow)  ! (m) snow layer thickness
real(kind=dbl_kind) :: dzsno_req(nsnow)  ! (m) snow layer thickness
real(kind=dbl_kind) :: z_botsno(nsnow)  ! (m) bottom of snow layer
real(kind=dbl_kind) :: z_midsno      ! (m) middle of snow layer
real(kind=dbl_kind) :: z_topsno(nsnow)  ! (m) top of snow layer
real(kind=dbl_kind) :: swice(nsnow)  ! (kg m-2) snow layer ice content
real(kind=dbl_kind) :: swliq(nsnow)  ! (kg m-2) snow layer water content
real(kind=dbl_kind) :: tsnow(nsnow)  ! (K) snow layer temperature
real(kind=dbl_kind) :: dsnow(nsnow)  ! (kg/m3) snow layer density
real(kind=dbl_kind) :: z_exceed      ! (m) excess layer thickness
real(kind=dbl_kind) :: move_frac     ! (-) fraction mass moved to another layer
real(kind=dbl_kind) :: stay_frac     ! (-) fraction mass staying in current layer
real(kind=dbl_kind) :: move_ice      ! (kg m-2) snow layer ice moved to another layer
real(kind=dbl_kind) :: move_liq      ! (kg m-2) snow layer water moved to another layer
real(kind=dbl_kind) :: mass_req(nsnow)      ! (kg m-2) snow layer required mass
real(kind=dbl_kind) :: mass_act(nsnow)      ! (kg m-2) snow layer actual mass
real(kind=dbl_kind) :: mass_tot      ! (kg m-2) total snow mass
real(kind=dbl_kind) :: del_mass      ! (kg m-2) total snow mass
real(kind=dbl_kind) :: sum           ! (kg m-2) total snow mass
real(kind=dbl_kind) :: fac_bot       ! (-) scaling fraction bottom layer
real(kind=dbl_kind) :: fac_top       ! (-) scaling fraction top layer
logical flag
!
! transfer snow data to local variables
    num_snow = abs(sib%prog%nsl)
    do n=1,num_snow
        dzsno(n) = sib%prog%dz(n+sib%prog%nsl)
        z_botsno(n) = -sib%prog%z_bot(n+sib%prog%nsl)
        z_topsno(n) = -sib%prog%z_top(n+sib%prog%nsl)
        swice(n) = sib%prog%www_ice(n+sib%prog%nsl)
        swliq(n) = sib%prog%www_liq(n+sib%prog%nsl)
        tsnow(n) = sib%prog%td(n+sib%prog%nsl)
    enddo
!
! bottom layer fraction
    d_slope=10.d0/(sib%param%d_max-sib%param%d_min)
    d_half=0.5d0*(sib%param%d_max+sib%param%d_min)
    f_bot=sib%param%f_bot_max/(1+exp(d_slope*(d_half-sib%prog%snow_depth)))
    if(f_bot<0.001) f_bot=0.
    d_bot=f_bot*sib%prog%snow_depth
!
! density difference between top and bottom layer
    delta_den=abs(sib%param%den_ref_top-sib%param%den_ref_bot)/sib%param%f_bot_max*f_bot
    s_delta_den=exp(-(sib%prog%snow_bulk-sib%param%den_bulk_obs)**2./sib%param%den_bulk_std**2.)
    delta_den=delta_den*s_delta_den
!
! prescribed density profiles
    i=sib%param%snow_class
    if(i==1.or.i==2.or.i==6) then ! constant density each layer
        snow_den_bot=sib%prog%snow_bulk-(1.-f_bot)*delta_den
        snow_den_top=sib%prog%snow_bulk+f_bot*delta_den
        num_snow = abs(sib%prog%nsl)
        do n=1,num_snow
            if(z_botsno(n)>d_bot) then
            dsnow(n)=snow_den_top
            elseif(z_topsno(n)<d_bot) then
              dsnow(n)=snow_den_bot
            else
              fac_bot=(d_bot-z_botsno(n))/dzsno(n)
              fac_top=1.-fac_bot
              dsnow(n)=fac_bot*snow_den_bot+fac_top*snow_den_top
            endif
        enddo
    elseif(i==3.or.i==4.or.i==5) then ! linear density top layer, constant density bottom layer
        snow_den_bot=sib%prog%snow_bulk+(1.-f_bot)*delta_den*.5
        den_slope=-delta_den/(sib%prog%snow_depth-d_bot)
        do n=1,num_snow
            if(z_botsno(n)>d_bot) then
              z_midsno=(z_botsno(n)+z_topsno(n))*.5
              dsnow(n)=den_slope*(z_midsno-d_bot)+snow_den_bot
            elseif(z_topsno(n)<d_bot) then
              dsnow(n)=snow_den_bot
            else
              z_midsno=(d_bot+z_topsno(n))*.5
              fac_bot=(d_bot-z_botsno(n))/dzsno(n)
              fac_top=1.-fac_bot
              dsnow(n)=fac_bot*snow_den_bot+fac_top*(den_slope*(z_midsno-d_bot)+snow_den_bot)
            endif
        enddo
    endif
    sib%diag%testvar1=delta_den
    sib%diag%testvar2=f_bot
    sib%diag%testvar3=delta_den
!
! calculate required mass per layer
    do n=1,num_snow
        mass_req(n)=dsnow(n)*dzsno(n)
        mass_act(n)=swice(n)+swliq(n)
    enddo
!
! check mass conservation
    mass_tot=0.d0
    sum=0.d0
    do n=1,num_snow
      mass_tot=mass_tot+swice(n)+swliq(n)
      del_mass=mass_req(n)-mass_act(n)
      sum=sum+del_mass
    enddo
!
! redistribute mass based on density
    do n=1,num_snow
      mass_act(n)=swice(n)+swliq(n)
      del_mass=mass_req(n)-mass_act(n)
      dzsno_req(n)=mass_act(n)/dsnow(n)

      if(num_snow > n  .AND.  del_mass < 0.d0 ) then  ! move snow down to lower layer
        z_exceed = 0.
        move_frac= abs(del_mass)/mass_act(n)
        stay_frac= 1.-move_frac
        move_ice = move_frac*swice(n)
        move_liq = move_frac*swliq(n)
        swice(n) = stay_frac*swice(n)
        swliq(n) = stay_frac*swliq(n)
        call clm_combo(dzsno(n+1),swliq(n+1),swice(n+1),tsnow(n+1),z_exceed,move_liq,move_ice,tsnow(n))
      endif

      if(num_snow > n  .AND.  del_mass > 0.d0 ) then  ! move snow up from below
        z_exceed = 0.
        move_frac= abs(del_mass)/mass_act(n+1)
        stay_frac= 1.-move_frac
        move_ice = move_frac*swice(n+1)
        move_liq = move_frac*swliq(n+1)
        swice(n+1) = stay_frac*swice(n+1)
        swliq(n+1) = stay_frac*swliq(n+1)
        call clm_combo(dzsno(n),swliq(n),swice(n),tsnow(n),z_exceed,move_liq,move_ice,tsnow(n+1))
      endif
    enddo
!
! check mass conservation again
    mass_tot=0.d0
    flag=.false.
    do n=1,num_snow
      mass_tot=mass_tot+swice(n)+swliq(n)
      mass_act(n)=swice(n)+swliq(n)
      del_mass=mass_req(n)-mass_act(n)
      if(abs(del_mass)>.001d0) flag=.true.
    enddo
    if(flag) stop 'snow mass not conserved'
!
! transfer back to variable tree
    do n=sib%prog%nsl+1,0
       sib%prog%www_ice(n)=swice(n-sib%prog%nsl)
       sib%prog%www_liq(n)=swliq(n-sib%prog%nsl)
       sib%prog%td(n)=tsnow(n-sib%prog%nsl)
       sib%prog%snow_den(n-sib%prog%nsl)=dsnow(n-sib%prog%nsl)
    enddo
!
end subroutine test_snow_density
!
!======================================================================
subroutine test_subdivide_snow(sib)
!======================================================================
! Description:
!   Subdivides snow layers if they exceed their prescribed
!   maximum thickness
!   Based on CLM subroutine CLM_SUBDIV
!
! Revision History:
!   15 September 1999: Yongjiu Dai, initial code
!   15 December  1999: Paul Houser and Jon Radakovich, F90 revision
!   30 January   2002: Ian Baker, SiB integration
!   12 July      2007: Kevin Schaefer moved thickness criteria to arrays
!   12 July      2007: Kevin Schaefer changed to generic do loop
!----------------------------------------------------------------------

use kinds
use sibtype

implicit none
!
! i/o variables
type(sib_t), intent(inout) :: sib
!
! local variables
integer(kind=int_kind) :: j,n  ! indeces
integer(kind=int_kind) :: num_snow   ! number of snow layers
real(kind=dbl_kind) :: dzsno(nsnow)  ! (m) snow layer thickness
real(kind=dbl_kind) :: swice(nsnow)  ! (kg m-2) snow layer ice content
real(kind=dbl_kind) :: swliq(nsnow)  ! (kg m-2) snow layer water content
real(kind=dbl_kind) :: tsnow(nsnow)  ! (K) snow layer temperature
real(kind=dbl_kind) :: z_exceed      ! (m) excess layer thickness
real(kind=dbl_kind) :: move_frac     ! (-) fraction mass moved to another layer
real(kind=dbl_kind) :: stay_frac     ! (-) fraction mass staying in current layer
real(kind=dbl_kind) :: move_ice      ! (kg m-2) snow layer ice moved to another layer
real(kind=dbl_kind) :: move_liq      ! (kg m-2) snow layer water moved to another layer
!
! transfer snow data to local variables
    num_snow = abs(sib%prog%nsl)
    do j=1,num_snow
        dzsno(j) = sib%prog%dz(j+sib%prog%nsl)
        swice(j) = sib%prog%www_ice(j+sib%prog%nsl)
        swliq(j) = sib%prog%www_liq(j+sib%prog%nsl)
        tsnow(j) = sib%prog%td(j+sib%prog%nsl)
    enddo
!
! loop through snow layers from top of snowpack to bottom
    do n=1,nsnow
!
! create new snow layer (if required)
      if(n == num_snow  .AND.  dzsno(n) > sib%prog%dz_new(n)) then ! new snow layer below current one
        num_snow = num_snow+1
        dzsno(n) = dzsno(n)/2.0
        swice(n) = swice(n)/2.0
        swliq(n) = swliq(n)/2.0
        dzsno(n+1) = dzsno(n)
        swice(n+1) = swice(n)
        swliq(n+1) = swliq(n)
        tsnow(n+1) = tsnow(n)
      endif
!
! redistribute snow (if required)
      if(num_snow > n  .AND.  dzsno(n) > sib%prog%dz_max(n) ) then  ! redistribute snow
        z_exceed = dzsno(n) - sib%prog%dz_max(n)
        move_frac= z_exceed/dzsno(n)
        stay_frac= 1.-move_frac
        move_ice = move_frac*swice(n)
        move_liq = move_frac*swliq(n)
        swice(n) = stay_frac*swice(n)
        swliq(n) = stay_frac*swliq(n)
        dzsno(n) = sib%prog%dz_max(n)
        call clm_combo(dzsno(n+1),swliq(n+1),swice(n+1),tsnow(n+1),z_exceed,move_liq,move_ice,tsnow(n))
      endif
    enddo
!
! transfer local data back to SiB variable tree
    sib%prog%nsl = -num_snow
    do j=sib%prog%nsl+1,0
        sib%prog%dz(j)      = dzsno(j-sib%prog%nsl)
        sib%prog%www_ice(j) = swice(j-sib%prog%nsl)
        sib%prog%www_liq(j) = swliq(j-sib%prog%nsl)
        sib%prog%td(j)      = tsnow(j-sib%prog%nsl)
    enddo
!
! recalculate snow layer geometry
    do j=0,sib%prog%nsl+1,-1
        sib%prog%z_bot(j)  = sib%prog%z_top(j+1)
        sib%prog%node_z(j) = sib%prog%z_bot(j) - 0.5 * sib%prog%dz(j)
        sib%prog%z_top(j)  = sib%prog%z_bot(j) - sib%prog%dz(j)
    enddo
!
end subroutine test_subdivide_snow
!
!----------------------------------------------------------------------
subroutine subdivide_snow(sib)
!----------------------------------------------------------------------
!
!   Based on CLM subroutine CLM_SUBDIV
!
!   Description
!   Subdivides snow layers if they exceed their prescribed 
!   maximum thickness
!
!   Revision History:
!   15 September 1999: Yongjiu Dai, initial code
!   15 December  1999: Paul Houser and Jon Radakovich, F90 revision
!   30 January   2002: Ian Baker, SiB integration
!----------------------------------------------------------------------

use kinds
use sibtype

implicit none

!----------------------------------------------------------------------

type(sib_t), intent(inout) :: sib

!----------------------------------------------------------------------  


!...local variables
integer(kind=int_kind) :: msno,j
real(kind=dbl_kind) :: dzsno(5)
real(kind=dbl_kind) :: swice(5)
real(kind=dbl_kind) :: swliq(5)
real(kind=dbl_kind) :: tsnow(5)
real(kind=dbl_kind) :: drr
real(kind=dbl_kind) :: propor
real(kind=dbl_kind) :: zwice
real(kind=dbl_kind) :: zwliq

    !if ( sib%prog%nsl == 0 ) return

    msno = abs(sib%prog%nsl)

    do j=1,msno
        dzsno(j) = sib%prog%dz(j+sib%prog%nsl)
        swice(j) = sib%prog%www_ice(j+sib%prog%nsl)
        swliq(j) = sib%prog%www_liq(j+sib%prog%nsl)
        tsnow(j) = sib%prog%td(j+sib%prog%nsl)
    enddo

    if(msno == 1) then
        if(dzsno(1) > 0.03) then
            msno = 2

            !...specify a new snow layer
            dzsno(1) = dzsno(1)/2.0
            swice(1) = swice(1)/2.0
            swliq(1) = swliq(1)/2.0

            dzsno(2) = dzsno(1)
            swice(2) = swice(1)
            swliq(2) = swliq(1)
            tsnow(2) = tsnow(1)
        endif

    endif   ! if msno == 1 condition


    if(msno > 1) then

        if(dzsno(1) > 0.02 ) then
            drr      = dzsno(1) - 0.02
            propor   = drr/dzsno(1)
            zwice    = propor*swice(1)
            zwliq    = propor*swliq(1)
            propor   = 0.02/dzsno(1)
            swice(1) = propor*swice(1)
            swliq(1) = propor*swliq(1)
            dzsno(1) = 0.02


            call clm_combo(dzsno(2),swliq(2),swice(2),tsnow(2),         &
                drr,zwliq,zwice,tsnow(1))


            if(msno <= 2  .AND. dzsno(2) > 0.07 ) then

                !...subdivide a new layer
                msno = 3
                dzsno(2) = dzsno(2)/2.0
                swice(2) = swice(2)/2.0
                swliq(2) = swliq(2)/2.0
                dzsno(3) = dzsno(2)
                swice(3) = swice(2)
                swliq(3) = swliq(2)
                tsnow(3) = tsnow(2)
            endif
        endif     ! if dzsno(1) > 0.02 condition
    endif       ! if msno > 1 condition


    if(msno > 2) then
        if(dzsno(2) > 0.05) then

            drr      = dzsno(2) - 0.05
            propor   = drr/dzsno(2)
            zwice    = propor*swice(2)
            zwliq    = propor*swliq(2)
            propor   = 0.05/dzsno(2)
            swice(2) = propor*swice(2)
            swliq(2) = propor*swliq(2)
            dzsno(2) = 0.05

            call clm_combo(dzsno(3),swliq(3),swice(3),tsnow(3),         &
                drr,zwliq,zwice,tsnow(2))



            if(msno <= 3  .AND.  dzsno(3) > 0.18) then

                !...subdivide a new layer
                msno = 4
                dzsno(3) = dzsno(3)/2.0
                swice(3) = swice(3)/2.0
                swliq(3) = swliq(3)/2.0
                dzsno(4) = dzsno(3)
                swice(4) = swice(3)
                swliq(4) = swliq(3)
                tsnow(4) = tsnow(3) 
            endif
        endif    ! if dzsno(2) > 0.05 condition
    endif      ! if msno > 2 condition


    if(msno > 3) then
        if(dzsno(3) > 0.11) then

            drr      = dzsno(3) - 0.11
            propor   = drr/dzsno(3)
            zwice    = propor*swice(3)
            zwliq    = propor*swliq(3)
            propor   = 0.11/dzsno(3)
            swice(3) = propor*swice(3)
            swliq(3) = propor*swliq(3)
            dzsno(3) = 0.11

            call clm_combo(dzsno(4),swliq(4),swice(4),tsnow(4),         &
                drr,zwliq,zwice,tsnow(3))


            if(msno <= 4  .AND.  dzsno(4) > 0.41) then

                !...subdivide a new layer
                msno = 5
                dzsno(4) = dzsno(4)/2.0
                swice(4) = swice(4)/2.0
                swliq(4) = swliq(4)/2.0
                dzsno(5) = dzsno(4)
                swice(5) = swice(4)
                swliq(5) = swliq(4)
                tsnow(5) = tsnow(4)
            endif
        endif    ! if dzsno(3) > 0.11 condition
    endif      ! if msno > 3 condition


    if(msno > 4) then
        if(dzsno(4) > 0.23) then
            drr      = dzsno(4) - 0.23
            propor   = drr/dzsno(4)
            zwice    = propor*swice(4)
            zwliq    = propor*swliq(4)
            propor   = 0.23/dzsno(4)
            swice(4) = propor*swice(4)
            swliq(4) = propor*swliq(4)
            dzsno(4) = 0.23

            call clm_combo(dzsno(5),swliq(5),swice(5),tsnow(5),         &
                drr,zwliq,zwice,tsnow(4))



        endif    ! if dzsno(4) > 0.23 condition
    endif      ! if msno > 4 condition

    sib%prog%nsl = -msno

    do j=sib%prog%nsl+1,0
        sib%prog%dz(j) = dzsno(j - sib%prog%nsl)
        sib%prog%www_ice(j) = swice(j - sib%prog%nsl)
        sib%prog%www_liq(j) = swliq(j - sib%prog%nsl)
        sib%prog%td(j)      = tsnow(j - sib%prog%nsl)
    enddo

    do j=0,sib%prog%nsl+1,-1

        sib%prog%node_z(j) = sib%prog%z_bot(j) - 0.5 * sib%prog%dz(j)
        sib%prog%z_bot(j-1) = sib%prog%node_z(j) - 0.5*sib%prog%dz(j)

    enddo


end subroutine subdivide_snow
