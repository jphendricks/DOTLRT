!
!=========================================================================
subroutine methane_core(sib)
!=========================================================================
! calculates methane fluxes from respiration flux
!
! Modifications:
!  Zhuxiao Li made subroutine (8/16/11)
!-------------------------------------------------------------------------

use kinds
use sibtype
use cfrax
use sib_const_module
use physical_parameters

implicit none
!
! input/output variables
type(sib_t), intent(inout) :: sib
!
! local variables
integer i,j,k,n  ! indeces
real (kind=dbl_kind) ms(nsoil)     ! (-) 
real (kind=dbl_kind) wetfrac     ! (-) 
!
! assign constants
   wetfrac=0.1_dbl_kind

!
! start methane calculations
   do i=1,nsoil
     j=sib%prog%satfrac(i)
!     sib%meth%ef_flood_lay(i)=0.19_dbl_kind/6._dbl_kind*sib%prog%satfrac(i)
!     sib%meth%ef_peat_lay(i)=0.005_dbl_kind/3._dbl_kind
   enddo
!
! methane calculation
!   sib%meth%flux=0._dbl_kind
 !  sib%meth%flux_lay(:)=0._dbl_kind
   do i=1,config%ndead_pool
     n=config%indx_dead(i)
     do k=1,nsoil
       ms(k) = 1.
!       sib%meth%flux=sib%meth%flux+sib%casa%resp_pool_lay(i,k)*sib%meth%ef_flood_lay(k)
!       sib%meth%flux_lay(k)=sib%meth%flux_lay(k)+sib%casa%resp_pool_lay(i,k)*ms(k)
     enddo
   enddo
!   sib%diag%testvar1=sib%meth%flux
!
! end methane calculations
!
end subroutine methane_core
!
!=========================================================================
subroutine methane_setup(sib)
!=========================================================================
! this sets up basic parameters for methane_core
! This routine will probably not be in SiBCASA
!
! Modifications:
!  Zhuxiao Li made subroutine (8/16/11)
!-------------------------------------------------------------------------

use kinds
use sibtype
use cfrax
use sib_const_module
use physical_parameters

implicit none
!
! input/output variables
type(sib_t), intent(inout) :: sib
!
! local variables
integer i,k,n  ! indeces
!
! identify dead pools
   npool=13
   config%ndead_pool=9
   n=0
   do i=5,13
     n=n+1
     config%indx_dead(n)=i
     sib%param%ncarb(n)=15
   enddo
   sib%param%ncarb(:)=15
!
! set input respiration rates (micromol/m2/s)
   sib%casa%resp_pool_lay(:,:)=0._dbl_kind
   do i=1,npool
     do k=1,sib%param%ncarb(i)
        sib%casa%resp_pool_lay(i,k) = .1_dbl_kind
     enddo
   enddo
!
end subroutine methane_setup
