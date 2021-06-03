!
!==========================================================================
subroutine get_pres_and_hght(ncdfID,outtim,we_ll,sn_ll,wepts,snpts, bt)
!==========================================================================
!  THIS SUBROUTINE CALCULATES THE FULL PRESSURE AND HEIGHT (THIS APPEARS TO
!  BE HEIGHT ABOVE 1000 OR 1013 MB; NEED TO CHECK) AT EACH GRID POINT.
!
! History:
!  01/20/2005 JASON OTKIN created routine SSEC/CIMSS/UW-MADISON
!  07/13/2007 JASON OTKIN MODIFIED TO IMPROVE EFFICIENCY
!  10/16/2020 Kevin Schaefer integrated routine into DOTLRT
!--------------------------------------------------------------------------

use netcdf
use WRF_MODEL_PARAMETERS
implicit none

include 'netcdf.inc'

integer,intent(in)			:: ncdfID, outtim
integer					:: we_ll,sn_ll,wepts,snpts,bt
real,allocatable,dimension(:,:,:)	:: pertp, basep
real,allocatable,dimension(:,:,:)	:: perth, baseh, fhght
integer					:: i,j,k
integer					:: status
integer					:: varID,vartype,natts
integer					:: vardims
integer					:: dims(4)
integer					:: vardimIDs(10)
character(len=15)				:: varnam

!-----------------------------------------------------------------------
!  RETRIEVE THE PERTURBATION AND BASE STATE PRESSURES
write(*,*) "Inside get_pres_and hght.F90"
allocate(pertp(wepts,snpts,bt))
allocate(basep(wepts,snpts,bt))

write(*,*) "Reading WRF Pressure Perturbation"
status = nf90_inq_varid(ncdfID, 'P', varID)
if(status /= nf90_noerr ) call handle_err(status)

status = nf90_get_var(ncdfID,varID,pertp,start=(/we_ll,sn_ll,1,outtim/), &
		      count=(/wepts,snpts,bt,1/) )
if(status /= nf90_noerr ) call handle_err(status)

write(*,*) "Reading WRF Base State Pressure"
status = nf90_inq_varid(ncdfID, 'PB', varID)
if(status /= nf90_noerr ) call handle_err(status)

status = nf90_get_var(ncdfID,varID,basep,start=(/we_ll,sn_ll,1,outtim/), &
		      count=(/wepts,snpts,bt,1/) )
if(status /= nf90_noerr ) call handle_err(status)

!-----------------------------------------------------------------------
!  CALCULATE THE TOTAL PRESSURE
write(*,*) "Compute WRF Total Pressure"
P = pertp + basep
		
deallocate(pertp,basep)

!------------------------------------------------------------------------------
!  RETRIEVE THE PERTURBATION AND BASE STATE GEOPOTENTIAL HEIGHT IN ORDER TO
!  CALCULATE THE GEOPOTENTIAL HEIGHT.  PERTURBATION AND BASE STATE GEOPOTENTIAL
!  HEIGHT ARE ON FULL ETA LEVELS SO MUST INTERPOLATE TO HALF LEVELS.

allocate(perth(wepts,snpts,bt+1))
allocate(baseh(wepts,snpts,bt+1))

status = nf90_inq_varid(ncdfID, 'PH', varID)
if(status /= nf90_noerr ) call handle_err(status)
!		
status = nf90_get_var(ncdfID,varID,perth,start=(/we_ll,sn_ll,1,outtim/), &
		      count=(/wepts,snpts,bt+1,1/) )
if(status /= nf90_noerr ) call handle_err(status)

status = nf90_inq_varid(ncdfID, 'PHB', varID)
if(status /= nf90_noerr ) call handle_err(status)

status = nf90_get_var(ncdfID,varID,baseh,start=(/we_ll,sn_ll,1,outtim/), &
		      count=(/wepts,snpts,bt+1,1/) )
if(status /= nf90_noerr ) call handle_err(status)

!------------------------------------------------------------------------------
! COMPUTE HEIGHT AT HALF LEVELS
!------------------------------------------------------------------------------

allocate(fhght(wepts,snpts,bt+1))
fhght = (perth+baseh)/9.81
ghgt=0.5*(fhght(:,:,1:bt)+fhght(:,:,2:bt+1))

deallocate(perth,baseh,fhght)

return
end subroutine get_pres_and_hght
