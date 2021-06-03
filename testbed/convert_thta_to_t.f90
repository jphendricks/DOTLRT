!---------------------------------------------------------------------------
!  THIS SUBROUTINE CONVERTS THE PERTURBATION POTENTIAL TEMPERATURE OUTPUT BY
!  THE WRF MODEL INTO THE ACTUAL TEMPERATURE TO BE OUTPUT INTO THE PSEUDO
!  MM5 FILE.  IT USES EQUATIONS MODIFIED FROM THE ripdp_wrf.f PROGRAM.

!  JASON OTKIN
!  SSEC/CIMSS/UW-MADISON
!  02-12-04
!
!  MODIFIED TO IMPROVE EFFICIENCY
!
!  JASON OTKIN
!  07/13/07
!  Modified on 27 May 2009 alim@ssec.wisc.edu
!---------------------------------------------------------------------------

	subroutine convert_thta_to_t(ncdfID, outtim, we_ll, sn_ll, wepts, &
					     snpts, bt)
	
	use netcdf
        use WRF_MODEL_PARAMETERS
	implicit none
	include 'netcdf.inc'

	integer,intent(in)			:: ncdfID, outtim
	
	integer					:: we_ll,sn_ll,wepts,snpts,bt
	
	integer					:: i,j,k
	
	integer					:: varID

	integer					:: status	
	real,allocatable,dimension(:,:,:)	:: pthta,qvp,totp,basep,pertp

	real						:: rgas,rgasmd,cp,cpmd,gamma,gammamd
	real						:: gammam
	real						:: bstatetemp

!---------------------------------------------------------------------------
!  DEFINE A FEW CONSTANTS
!---------------------------------------------------------------------------

	rgas=287.04  		!J/K/kg
      rgasmd=.608   		!rgas_moist=rgas*(1.+rgasmd*qvp)
      cp=1004.     		!J/K/kg  Note:not using Bolton's value of 1005.7
      cpmd=.887   		!cp_moist=cp*(1.+cpmd*qvp)
      gamma=rgas/cp
      gammamd=rgasmd-cpmd  	!gamma_moist=gamma*(1.+gammamd*qvp)

	bstatetemp=300.

!---------------------------------------------------------------------------
!  RETRIEVE THE PERTURBATION AND BASIC STATE PRESSURE, COMBINE, AND
!  CONVERT TO HPA
!---------------------------------------------------------------------------
		
	allocate(pertp(wepts,snpts,bt))
	allocate(basep(wepts,snpts,bt))
	allocate( totp(wepts,snpts,bt))
	
	status = nf90_inq_varid(ncdfID, 'P', varID)
	if(status /= nf90_noerr ) call handle_err(status)
	
	status = nf90_get_var(ncdfID,varID,pertp,start=(/we_ll,sn_ll,1,outtim/), &
			      count=(/wepts,snpts,bt,1/) )
	if(status /= nf90_noerr ) call handle_err(status)
	    
	status = nf90_inq_varid(ncdfID, 'PB', varID)
	if(status /= nf90_noerr ) call handle_err(status)
	
	status = nf90_get_var(ncdfID,varID,basep,start=(/we_ll,sn_ll,1,outtim/), &
			      count=(/wepts,snpts,bt,1/) )
	if(status /= nf90_noerr ) call handle_err(status)

	totp = .01 * (pertp+basep)

!---------------------------------------------------------------------------
!  RETRIEVE WATER VAPOR MIXING RATIO
!---------------------------------------------------------------------------

!	allocate(qvp(wepts,snpts,bt))

!	status = nf90_inq_varid(ncdfID, 'QVAPOR', varID)
!	if(status /= nf90_noerr ) call handle_err(status)

!	status = nf90_get_var(ncdfID,varID,qvp,start=(/we_ll,sn_ll,1,outtim/), &
!			      count=(/wepts,snpts,bt,1/) )
!	if(status /= nf90_noerr ) call handle_err(status)

!---------------------------------------------------------------------------
!  RETRIEVE PERTURBATION POTENTIAL TEMPERATURE
!---------------------------------------------------------------------------

	allocate(pthta(wepts,snpts,bt))

        write(*,*) "Reading WRF Potential Temperature Perturbation"
	status = nf90_inq_varid(ncdfID, 'T', varID)
	if(status /= nf90_noerr ) call handle_err(status)
	
	status = nf90_get_var(ncdfID,varID,pthta,start=(/we_ll,sn_ll,1,outtim/), &
			      count=(/wepts,snpts,bt,1/) )
	if(status /= nf90_noerr ) call handle_err(status)
	       
!---------------------------------------------------------------------------
!  CALCULATE TEMPERATURE
!---------------------------------------------------------------------------

        write(*,*) "Convert Potentual Temperature to Temperature"
	do k = 1,bt
	  do j = 1,snpts
	    do i = 1,wepts
	    	
!		gammam = gamma*(1.0+gammamd*qvp(i,j,k))
!		T(i,j,k) = (pthta(i,j,k)+bstatetemp)*(totp(i,j,k)/1000.)**gammam
                T(i,j,k) = (pthta(i,j,k)+bstatetemp)*(totp(i,j,k)/1000.)**gamma
	    enddo
	  enddo
	enddo	 
	
!	deallocate(pthta,qvp,totp,basep,pertp)
	deallocate(pthta,totp,basep,pertp)
	
	return
	end subroutine convert_thta_to_t
