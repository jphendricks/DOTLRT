! This file contains unused, outdated pre-fortran4 subroutines from DOTLRT



!====================================================================
subroutine read_stream_angles(filename, num_streams)
!====================================================================
! reads in single profile from WRF run
! Assumes csv format
! primarily for testing dotlrt
!
! History:
!  9/18/20  Kevin Schaefer created routine
!--------------------------------------------------------------------
  use variables
  implicit none
!
! input variables
  character*250 filename ! input filename
  integer num_streams    ! (-) number of streams
!
! internal variables  
  integer iang       ! (-) angle index
  real(8) temp(2)    ! (varies) temporary read variable
  character*250 junk ! junk read variable
!
! open profile file
  open(unit=20, file=trim(filename), form='formatted', status='old')
  read(20,*) junk
!
! read profile
  do iang = 1, num_streams
    read(20,*) temp
    stream_ang(iang)=temp(2)
    print*, iang, stream_ang(iang)
  end do
!
! close profile file
  close(unit=20)
!
end subroutine read_stream_angles


!====================================================================
subroutine getprofile( )
!====================================================================
!
! History:
!   9/26/2020 Kevin Schaefer deleted unused variables
!  10/16/2020 Kevin Schaefer deleted routine from MRT
!--------------------------------------------------------------------
  use variables
  implicit none
  integer ilev
  
  ! surface at first level: height = 0 (agl)
  surf_inp%surf_temp = atminp(1,3)

  do ilev = 1, atm_inp%num_levels
    atm_inp%prof(ilev)%height          = atminp(ilev+1,1) ! (km)
    atm_inp%prof(ilev)%pressure        = atminp(ilev+1,2) ! (mb)
    atm_inp%prof(ilev)%temperature     = atminp(ilev+1,3) ! (K)
    atm_inp%prof(ilev)%vapor_density   = atminp(ilev+1,4) ! (g/m^3)
    atm_inp%prof(ilev)%cloud_liq_dens  = atminp(ilev+1,5) ! (g/m^3)
    atm_inp%prof(ilev)%cloud_rn_dens   = atminp(ilev+1,6) ! (g/m^3)
    atm_inp%prof(ilev)%cloud_ice_dens  = atminp(ilev+1,7) ! (g/m^3)
    atm_inp%prof(ilev)%cloud_snow_dens = atminp(ilev+1,8) ! (g/m^3)
    atm_inp%prof(ilev)%cloud_grpl_dens = atminp(ilev+1,9) ! (g/m^3)
  end do
  atm_inp%inf= .true.
  atm_inp%fiveph= .true.

 end subroutine getprofile


!===================================================================
subroutine output( jchan, atm_data )
!===================================================================
! History:
!  24 May 2005    Bob Weber created routine
!  9/26/20 Kevin Schaefer deleted unused variabes
!-------------------------------------------------------------------
  use variables
  implicit none
  type(prof_data_type)    :: atm_data
  integer ilr1, hydrometeor_phase, jchan
  real(8) z, zout(0:max_num_levels)
  real(8) freqghz
  freqghz = instr_spec%chan(jchan)%lo_freq
  do i0 = 1, nangover2
    ! Product of geophysical and radiation Jacobians:
    z=0.d0
    do ilr1=1,nlr1
      z=z+h1(ilr1)
      zout(ilr1) = z
    end do
    jrec_tot(i0) = jrec_tot(i0) + 1
    write(tot_unit(i0),rec=jrec_tot(i0)) &
          sngl(atm_data%lat), sngl(atm_data%lng), nlr1, &
          sngl(freqghz), i0, sngl(dacos(cs(i0))*180.0d0/pi), &
        ( sngl(zout(ilr1)), sngl(dTb_dT(ilr1,i0,1)), sngl(dTb_dp(ilr1,i0,1)), &
          sngl(dTb_dq(ilr1,i0,1)), &
        ( sngl(dTb_dw(ilr1,i0,hydrometeor_phase,1)), &
                      hydrometeor_phase = 1, number_h2o_phases ), &
                    ilr1 = 1, nlr1 )
    z=0.d0
    zout(0) = z
    DO ilr1=1,nlr1
      z=z+altitude1(ilr1+1) - altitude1(ilr1)
      zout(ilr1) = z
    end do
              if( i0 == 1 ) then
                ! Temperatures:
                jrec_rt(3,i0) = jrec_rt(3,i0) + 1
                write( rt_unit(3,i0),rec=jrec_rt(3,i0)) sngl(atm_data%lat), sngl(atm_data%lng), &
                                            nlr1, i0, sngl(dacos(cs(i0))*180.0d0/pi), &
                                        sngl(zout(0)), sngl(surf_inp%surf_temp), &
                             ( sngl(zout(ilr1)), sngl(temperature1(ilr1)), ilr1 = 1, nlr1 )
                ! Gaseous absorption and combined hydrometeor ansorption/scattering coefficients:
                jrec_rt(4,i0) = jrec_rt(4,i0) + 1
                write( rt_unit(4,i0),rec=jrec_rt(4,i0)) sngl(atm_data%lat), sngl(atm_data%lng), &
                                               nlr1, i0, sngl(dacos(cs(i0))*180.0d0/pi), &
                                      (sngl(zout(ilr1)), sngl(al_gas1(ilr1)), &
                        sngl(scat_cloud1(ilr1)), sngl(abs_cloud1(ilr1)), sngl(abs_O2_1(ilr1)), &
                        sngl(abs_H2O_1(ilr1)), ilr1 = 1, nlr1 )
              end if
              z=0.d0
              DO ilr1=0, nlr1
                IF(ilr1 /= 0) z=z+h1(ilr1)
                zout(ilr1) = z
              end do
              ! Brightness temperatures:
              jrec_rt(1,i0) = jrec_rt(1,i0) + 1
              write( rt_unit(1,i0),rec=jrec_rt(1,i0)) sngl(atm_data%lat), sngl(atm_data%lng), &
                                            nlr1, i0, sngl(dacos(cs(i0))*180.0d0/pi), sngl(freqghz), &
                                           ( sngl(zout(ilr1)), sngl(tb_pl(ilr1,i0)), &
                                           sngl(tb_mn(ilr1,i0)), ilr1 = 0, nlr1 )
              z=0.d0
              DO ilr1=1,nlr1
                z=z+h1(ilr1)
                zout(ilr1) = z
              end do
              ! Radiation Jacobian:
              jrec_rt(2,i0) = jrec_rt(2,i0) + 1
                write( rt_unit(2,i0),rec=jrec_rt(2,i0)) sngl(atm_data%lat), sngl(atm_data%lng), &
                                              (nlr1-1), i0, sngl(dacos(cs(i0))*180.0d0/pi), sngl(freqghz), &
                                            (sngl(zout(ilr1)),                                    &
                        sngl(dTbdKa(ilr1,i0,1)), sngl(dTbdKa(ilr1,i0,2)), &
                        sngl(dTbdTr(ilr1,i0,1)), sngl(dTbdTr(ilr1,i0,2)), &
                        sngl(dTbdKsliq(ilr1,i0,1)), sngl(dTbdKsliq(ilr1,i0,2)), &
                        sngl(dTbdgliq(ilr1,i0,1)), sngl(dTbdgliq(ilr1,i0,2)), &
                        sngl(dTbdKsrn(ilr1,i0,1)), sngl(dTbdKsrn(ilr1,i0,2)), &
                        sngl(dTbdgrn(ilr1,i0,1)), sngl(dTbdgrn(ilr1,i0,2)), &
                        sngl(dTbdKsice(ilr1,i0,1)), sngl(dTbdKsice(ilr1,i0,2)), &
                        sngl(dTbdgice(ilr1,i0,1)), sngl(dTbdgice(ilr1,i0,2)), &
                        sngl(dTbdKssnow(ilr1,i0,1)), sngl(dTbdKssnow(ilr1,i0,2)), &
                        sngl(dTbdgsnow(ilr1,i0,1)), sngl(dTbdgsnow(ilr1,i0,2)), &
                        sngl(dTbdKsgrpl(ilr1,i0,1)), sngl(dTbdKsgrpl(ilr1,i0,2)), &
                        sngl(dTbdggrpl(ilr1,i0,1)), sngl(dTbdggrpl(ilr1,i0,2)), &
                                                       ilr1 = 1, nlr1)
            end do ! i0
            ! Geophysical Jacobian:
            jrec_geo = jrec_geo + 1
            write(geo_unit,rec=jrec_geo) sngl(atm_data%lat), sngl(atm_data%lng), nlr1, sngl(freqghz), &
                                           ( sngl(zout(ilr1)), &
                                             sngl(dKab_dT(ilr1)), &
                                             sngl(dKab_dp(ilr1)), &
                                             sngl(dKab_dq(ilr1)), &
                                           ( sngl(dKsc_dT(ilr1,hydrometeor_phase)), &
                                             sngl(  dg_dT(ilr1,hydrometeor_phase)), &
                                              hydrometeor_phase = 1, number_h2o_phases ), &
                                           ( sngl(dKab_dw(ilr1,hydrometeor_phase)), &
                                             sngl(dKsc_dw(ilr1,hydrometeor_phase)), &
                                             sngl(  dg_dw(ilr1,hydrometeor_phase)), &
                                              hydrometeor_phase = 1, number_h2o_phases ), ilr1 = 1, nlr1 )
            ! Hydrometeor absorption and scattering coefficients and asymmetry factor for all phases:
            z=0.d0
            do ilr1=1,nlr1
              z=z+h1(ilr1)
              zout(ilr1) = z
            end do
            jrec_Ksa = jrec_Ksa + 1
  write( Ksa_unit,rec=jrec_Ksa) sngl(atm_data%lat), sngl(atm_data%lng), nlr1, sngl(freqghz), &
        ( sngl(      zout(ilr1)), &
        ( sngl(hydro_prof(ilr1,hydrometeor_phase)%cloudab), &
          sngl(hydro_prof(ilr1,hydrometeor_phase)%cloudsc), &
          sngl(hydro_prof(ilr1,hydrometeor_phase)%cloudg),  &
              hydrometeor_phase = 1, number_h2o_phases ), &
              ilr1 = 1, nlr1)
end subroutine output


!**********************************************************************
!**********************************************************************
      subroutine seteru (messg, nmessg, nerr, iopt)
      common /cseter/ iunflo
      integer messg(1)
      data iunflo / 0 /
!
      if (iopt.ne.0) call seterr (messg, nmessg, nerr, iopt)
      if (iopt.ne.0) return
!
      if (iunflo.le.0) return
      call seterr (messg, nmessg, nerr, 1)
!
      return
      end

!**********************************************************************
!**********************************************************************
      subroutine seterr (messg, nmessg, nerr, iopt)
!
!  this version modified by w. fullerton to dump if iopt = 1 and
!  not recovering.
!  seterr sets lerror = nerr, optionally prints the message and dumps
!  according to the following rules...
!
!    if iopt = 1 and recovering      - just remember the error.
!    if iopt = 1 and not recovering  - print, dump and stop.
!    if iopt = 2                     - print, dump and stop.
!
!  input
!
!    messg  - the error message.
!    nmessg - the length of the message, in characters.
!    nerr   - the error number. must have nerr non-zero.
!    iopt   - the option. must have iopt=1 or 2.
!
!  error states -
!
!    1 - message length not positive.
!    2 - cannot have nerr=0.
!    3 - an unrecovered error followed by another error.
!    4 - bad value for iopt.
!
!  only the first 72 characters of the message are printed.
!
!  the error handler calls a subroutine named fdump to produce a
!  symbolic dump. to complete the package, a dummy version of fdump
!  is supplied, but it should be replaced by a locally written version
!  which at least gives a trace-back.
!
      integer messg(1)
      external i1mach, i8save
!
!  the unit for error messages.
!
      iwunit=i1mach(4)
!
      if (nmessg.ge.1) go to 10
!
!  a message of non-positive length is fatal.
!
        write(iwunit,9000)
 9000   format(52h1error    1 in seterr - message length not positive.)
        go to 60
!
!  nw is the number of words the message occupies.
!
 10   nw=(min0(nmessg,72)-1)/i1mach(6)+1
!
      if (nerr.ne.0) go to 20
!
!  cannot turn the error state off using seterr.
!
        write(iwunit,9001)
 9001   format(42h1error    2 in seterr - cannot have nerr=0// &
               34h the current error message follows///)
        call e9rint(messg,nw,nerr,.true.)
        itemp=i8save(1,1,.true.)
        go to 50
!
!  set lerror and test for a previous unrecovered error.
!
 20   if (i8save(1,nerr,.true.).eq.0) go to 30
!
        write(iwunit,9002)
 9002   format(23h1error    3 in seterr -, &
               48h an unrecovered error followed by another error.// &
               48h the previous and current error messages follow.///)
        call eprint
        call e9rint(messg,nw,nerr,.true.)
        go to 50
!
!  save this message in case it is not recovered from properly.
!
 30   call e9rint(messg,nw,nerr,.true.)
!
      if (iopt.eq.1 .or. iopt.eq.2) go to 40
!
!  must have iopt = 1 or 2.
!
        write(iwunit,9003)
 9003   format(42h1error    4 in seterr - bad value for iopt// &
               34h the current error message follows///)
        go to 50
!
!  test for recovery.
!
 40   if (iopt.eq.2) go to 50
!
      if (i8save(2,0,.false.).eq.1) return
!
 50   call eprint
 60   call fdump
      stop
!
      end

!**********************************************************************
!**********************************************************************
      subroutine e9rint(messg,nw,nerr,save)
!
!  this routine stores the current error message or prints the old one,
!  if any, depending on whether or not save = .true. .
!
      integer messg(nw)
      logical save
      external i1mach, i8save
!
!  messgp stores at least the first 72 characters of the previous
!  message. its length is machine dependent and must be at least
!
!       1 + 71/(the number of characters stored per integer word).
!
      integer messgp(36),fmt(14),ccplus
!
!  start with no previous message.
!
      data messgp(1)/1h1/, nwp/0/, nerrp/0/
!
!  set up the format for printing the error message.
!  the format is simply (a1,14x,72axx) where xx=i1mach(6) is the
!  number of characters stored per integer word.
!
      data ccplus  / 1h+ /
!
      data fmt( 1) / 1h( /
      data fmt( 2) / 1ha /
      data fmt( 3) / 1h1 /
      data fmt( 4) / 1h, /
      data fmt( 5) / 1h1 /
      data fmt( 6) / 1h4 /
      data fmt( 7) / 1hx /
      data fmt( 8) / 1h, /
      data fmt( 9) / 1h7 /
      data fmt(10) / 1h2 /
      data fmt(11) / 1ha /
      data fmt(12) / 1hx /
      data fmt(13) / 1hx /
      data fmt(14) / 1h) /
!
      if (.not.save) go to 20
!
!  save the message.
!
        nwp=nw
        nerrp=nerr
        do 10 i=1,nw
 10     messgp(i)=messg(i)
!
        go to 30
!
 20   if (i8save(1,0,.false.).eq.0) go to 30
!
!  print the message.
!
        iwunit=i1mach(4)
        write(iwunit,9000) nerrp
 9000   format(7h error ,i4,4h in )
!
        call s88fmt(2,i1mach(6),fmt(12))
        write(iwunit,fmt) ccplus,(messgp(i),i=1,nwp)
!
 30   return
!
      end

!**********************************************************************
!**********************************************************************
      subroutine eprint
!
!  this subroutine prints the last error message, if any.
!
      integer messg(1)
!
      call e9rint(messg,1,1,.false.)
      return
!
      end

!======================================================================
      subroutine s88fmt( n, w, ifmt )
!======================================================================
! Replaces ifmt(1-n) with the characters corresponding to the n
! least significant digits of w.
!
      integer n,w,ifmt(n)
!
      integer nt,wt,digits(10)
!
      data digits( 1) / 1h0 /
      data digits( 2) / 1h1 /
      data digits( 3) / 1h2 /
      data digits( 4) / 1h3 /
      data digits( 5) / 1h4 /
      data digits( 6) / 1h5 /
      data digits( 7) / 1h6 /
      data digits( 8) / 1h7 /
      data digits( 9) / 1h8 /
      data digits(10) / 1h9 /
!
      nt = n
      wt = w
!
 10   if (nt .le. 0) return
        idigit = mod( wt, 10 )
        ifmt(nt) = digits(idigit+1)
        wt = wt/10
        nt = nt - 1
        go to 10
!
      end

!**********************************************************************
!**********************************************************************
!DECK FDUMP
      SUBROUTINE FDUMP
!***BEGIN PROLOGUE  FDUMP
!***PURPOSE  Symbolic dump (should be locally written).
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3
!***TYPE      ALL (FDUMP-A)
!***KEYWORDS  ERROR, XERMSG
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!        ***Note*** Machine Dependent Routine
!        FDUMP is intended to be replaced by a locally written
!        version which produces a symbolic dump.  Failing this,
!        it should be replaced by a version which prints the
!        subprogram nesting list.  Note that this dump must be
!        printed on each of up to five files, as indicated by the
!         XGETUA routine.  See XSETUA and XGETUA for details.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  FDUMP
!***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END

!***********************************************************************
!***********************************************************************






!**********************************************************************
!**********************************************************************
      subroutine d9gaml (xmin, xmax)
! june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
!
! calculate the minimum and maximum legal bounds for x in gamma(x).
! xmin and xmax are not the only bounds, but they are the only non-
! trivial ones to calculate.
!
!             output arguments --
! xmin   dble prec minimum legal value of x in gamma(x).  any smaller
!        value of x might result in underflow.
! xmax   dble prec maximum legal value of x in gamma(x).  any larger
!        value of x might cause overflow.
!
      double precision xmin, xmax, alnbig, alnsml, xln, xold, d1mach, &
        dlog
      external d1mach, dlog
!
      alnsml = dlog(d1mach(1))
      xmin = -alnsml
      do 10 i=1,10
        xold = xmin
        xln = dlog(xmin)
        xmin = xmin - xmin*((xmin+0.5d0)*xln - xmin - 0.2258d0 + &
          alnsml) / (xmin*xln+0.5d0)
        if (dabs(xmin-xold).lt.0.005d0) go to 20
 10   continue
 20   xmin = -xmin + 0.01d0
!
      alnbig = dlog (d1mach(2))
      xmax = alnbig
      do 30 i=1,10
        xold = xmax
        xln = dlog(xmax)
        xmax = xmax - xmax*((xmax-0.5d0)*xln - xmax + 0.9189d0 - &
        alnbig) / (xmax*xln-0.5d0)
        if (dabs(xmax-xold).lt.0.005d0) go to 40
 30   continue
 40   xmax = xmax - 0.01d0
      xmin = dmax1 (xmin, -xmax+1.d0)
!
      return
      end

!**********************************************************************
!**********************************************************************
      double precision function d9lgmc (x)
! august 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!
! compute the log gamma correction factor for x .ge. 10. so that
! dlog (dgamma(x)) = dlog(dsqrt(2*pi)) + (x-.5)*dlog(x) - x + d9lgmc(x)
!
      double precision x, algmcs(15), xbig, xmax, dcsevl, d1mach, &
        dexp, dlog, dsqrt
      external d1mach, dcsevl, dexp, dlog, dsqrt, initds
!
! series for algm       on the interval  0.          to  1.00000e-02
!                                        with weighted error   1.28e-31
!                                         log weighted error  30.89
!                               significant figures required  29.81
!                                    decimal places required  31.48
!
      data algmcs(  1) /  .1666389480451863247205729650822d+0      /
      data algmcs(  2) / -.1384948176067563840732986059135d-4      /
      data algmcs(  3) /  .9810825646924729426157171547487d-8      /
      data algmcs(  4) / -.1809129475572494194263306266719d-10     /
      data algmcs(  5) /  .6221098041892605227126015543416d-13     /
      data algmcs(  6) / -.3399615005417721944303330599666d-15     /
      data algmcs(  7) /  .2683181998482698748957538846666d-17     /
      data algmcs(  8) / -.2868042435334643284144622399999d-19     /
      data algmcs(  9) /  .3962837061046434803679306666666d-21     /
      data algmcs( 10) / -.6831888753985766870111999999999d-23     /
      data algmcs( 11) /  .1429227355942498147573333333333d-24     /
      data algmcs( 12) / -.3547598158101070547199999999999d-26     /
      data algmcs( 13) /  .1025680058010470912000000000000d-27     /
      data algmcs( 14) / -.3401102254316748799999999999999d-29     /
      data algmcs( 15) /  .1276642195630062933333333333333d-30     /
!
      data nalgm, xbig, xmax / 0, 2*0.d0 /
!
      if (nalgm.ne.0) go to 10
      nalgm = initds (algmcs, 15, sngl(d1mach(3)) )
      xbig = 1.0d0/dsqrt(d1mach(3))
      xmax = dexp (dmin1(dlog(d1mach(2)/12.d0), -dlog(12.d0*d1mach(1))))
!
 10   continue
      if (x.ge.xmax) go to 20
!
      d9lgmc = 1.d0/(12.d0*x)
      if (x.lt.xbig) d9lgmc = dcsevl (2.0d0*(10.d0/x)**2-1.d0, algmcs, &
        nalgm) / x
      return
!
 20   d9lgmc = 0.d0
      return
!
      end

!**********************************************************************
!**********************************************************************
      double precision function d9pak (y, n)
! december 1979 edition. w. fullerton, c3, los alamos scientific lab.
!
! pack a base 2 exponent into floating point number x.  this routine is
! almost the inverse of d9upak.  it is not exactly the inverse, because
! dabs(x) need not be between 0.5 and 1.0.  if both d9pak and 2.d0**n
! were known to be in range we could compute
!                d9pak = x * 2.0d0**n
!
      double precision y, aln2b, aln210, d1mach
      external d1mach, i1mach
      data nmin, nmax / 2*0 /
      data aln210 / 3.321928094887362347870319429489d0 /
!
      if (nmin.ne.0) go to 10
      aln2b = 1.0d0
      if (i1mach(10).ne.2) aln2b = d1mach(5)*aln210
      nmin = aln2b*dble(float(i1mach(15)))
      nmax = aln2b*dble(float(i1mach(16)))
!
 10   call d9upak (y, d9pak, ny)
!
      nsum = n + ny
      if (nsum.lt.nmin) go to 40
!
      if (nsum.eq.0) return
      if (nsum.gt.0) go to 30
!
 20   d9pak = 0.5d0*d9pak
      nsum = nsum + 1
      if (nsum.ne.0) go to 20
      return
!
 30   d9pak = 2.0d0*d9pak
      nsum = nsum - 1
      if (nsum.ne.0) go to 30
      return
!
 40   continue 
      d9pak = 0.0d0
      return
!
      end

!**********************************************************************
!**********************************************************************
      subroutine d9upak (x, y, n)
! august 1980 portable edition.  w. fullerton, los alamos scientific lab
!
! unpack floating point number x so that x = y * 2.0**n, where
! 0.5 .le. abs(y) .lt. 1.0 .
!
      double precision x, y, absx
!
      absx = dabs(x)
      n = 0
      y = 0.0d0
      if (x.eq.0.0d0) return
!
 10   if (absx.ge.0.5d0) go to 20
      n = n - 1
      absx = absx*2.0d0
      go to 10
!
 20   if (absx.lt.1.0d0) go to 30
      n = n + 1
      absx = absx*0.5d0
      go to 20
!
 30   y = dsign (absx, x)
      return
!
      end

!**********************************************************************
!**********************************************************************
      double precision function dcsevl (x, a, n)
!
! evaluate the n-term chebyshev series a at x.  adapted from
! r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973).
!
!             input arguments --
! x      dble prec value at which the series is to be evaluated.
! a      dble prec array of n terms of a chebyshev series.  in eval-
!        uating a, only half the first coef is summed.
! n      number of terms in array a.
!
      double precision a(n), x, twox, b0, b1, b2
!
      twox = 2.0d0*x
      b1 = 0.d0
      b0 = 0.d0
      do 10 i=1,n
        b2 = b1
        b1 = b0
        ni = n - i + 1
        b0 = twox*b1 - b2 + a(ni)
 10   continue
!
      dcsevl = 0.5d0 * (b0-b2)
!
      return
      end

!**********************************************************************
!**********************************************************************
      double precision function dexp (x)
! may 1980 edition.   w. fullerton, c3, los alamos scientific lab.
      double precision x, expcs(14), twon16(17), aln216, f, xint, xmax,&
        xmin, y,  d1mach, dint, d9pak, dcsevl, dlog
      external d1mach, d9pak, dcsevl, dint, dlog, initds
!
! series for exp        on the interval -1.00000e+00 to  1.00000e+00
!                                        with weighted error   2.30e-34
!                                         log weighted error  33.64
!                               significant figures required  32.28
!                                    decimal places required  34.21
!
      data expcs(  1) / .866569493314985712733404647266231d-1    /
      data expcs(  2) / .938494869299839561896336579701203d-3    /
      data expcs(  3) / .677603970998168264074353014653601d-5    /
      data expcs(  4) / .366931200393805927801891250687610d-7    /
      data expcs(  5) / .158959053617461844641928517821508d-9    /
      data expcs(  6) / .573859878630206601252990815262106d-12   /
      data expcs(  7) / .177574448591421511802306980226000d-14   /
      data expcs(  8) / .480799166842372422675950244533333d-17   /
      data expcs(  9) / .115716376881828572809260000000000d-19   /
      data expcs( 10) / .250650610255497719932458666666666d-22   /
      data expcs( 11) / .493571708140495828480000000000000d-25   /
      data expcs( 12) / .890929572740634240000000000000000d-28   /
      data expcs( 13) / .148448062907997866666666666666666d-30   /
      data expcs( 14) / .229678916630186666666666666666666d-33   /
!
! twon16(i) is 2.0**((i-1)/16) - 1.0
      data twon16(  1) / 0.0d0                                   /
      data twon16(  2) / .44273782427413840321966478739929d-1    /
      data twon16(  3) / .90507732665257659207010655760707d-1    /
      data twon16(  4) / .13878863475669165370383028384151d0     /
      data twon16(  5) / .18920711500272106671749997056047d0     /
      data twon16(  6) / .24185781207348404859367746872659d0     /
      data twon16(  7) / .29683955465100966593375411779245d0     /
      data twon16(  8) / .35425554693689272829801474014070d0     /
      data twon16(  9) / .41421356237309504880168872420969d0     /
      data twon16( 10) / .47682614593949931138690748037404d0     /
      data twon16( 11) / .54221082540794082361229186209073d0     /
      data twon16( 12) / .61049033194925430817952066735740d0     /
      data twon16( 13) / .68179283050742908606225095246642d0     /
      data twon16( 14) / .75625216037329948311216061937531d0     /
      data twon16( 15) / .83400808640934246348708318958828d0     /
      data twon16( 16) / .91520656139714729387261127029583d0     /
      data twon16( 17) / 1.d0                                    /
!
! aln216 is 16.0/alog(2.) - 23.0
      data aln216 / .83120654223414517758794896030274d-1    /
      data nterms, xmin, xmax /0, 2*0.d0 /
!
      if (nterms.ne.0) go to 10
      nterms = initds (expcs, 14, 0.1*sngl(d1mach(3)))
      xmin = dlog (d1mach(1)) + .01d0
      xmax = dlog (d1mach(2)) - 0.001d0
!
 10   if (x.lt.xmin) go to 20
!
      xint = dint (x)
      y = x - xint
!
      y = 23.d0*y + x*aln216
      n = y
      f = y - dble(float(n))
      n = 23.d0*xint + dble(float(n))
      n16 = n/16
      if (n.lt.0) n16 = n16 - 1
      ndx = n - 16*n16 + 1
!
      dexp = 1.0d0 + (twon16(ndx) + f*(1.0d0 + twon16(ndx)) * &
        dcsevl (f, expcs, nterms) )
!
      dexp = d9pak (dexp, n16)
      return
!
 20   continue 
      dexp = 0.d0
      return
!
      end

!**********************************************************************
!**********************************************************************
      double precision function dgamma (x)
! jan 1984 edition.  w. fullerton, c3, los alamos scientific lab.
! jan 1994 wpp@ips.id.ethz.ch, ehg@research.att.com   declare xsml
      double precision x, gamcs(42), dxrel, pi, sinpiy, sq2pil, xmax, &
        xmin, y, d9lgmc, dcsevl, d1mach, dexp, dint, dlog, &
        dsin, dsqrt, xsml
      external d1mach, d9lgmc, dcsevl, dexp, dint, dlog, dsin, dsqrt, &
        initds
!
! series for gam        on the interval  0.          to  1.00000e+00
!                                        with weighted error   5.79e-32
!                                         log weighted error  31.24
!                               significant figures required  30.00
!                                    decimal places required  32.05
!
      data gamcs(  1) / .8571195590989331421920062399942d-2      /
      data gamcs(  2) / .4415381324841006757191315771652d-2      /
      data gamcs(  3) / .5685043681599363378632664588789d-1      /
      data gamcs(  4) / -.4219835396418560501012500186624d-2     /
      data gamcs(  5) / .1326808181212460220584006796352d-2      /
      data gamcs(  6) / -.1893024529798880432523947023886d-3     /
      data gamcs(  7) / .3606925327441245256578082217225d-4      /
      data gamcs(  8) / -.6056761904460864218485548290365d-5     /
      data gamcs(  9) / .1055829546302283344731823509093d-5      /
      data gamcs( 10) / -.1811967365542384048291855891166d-6     /
      data gamcs( 11) / .3117724964715322277790254593169d-7      /
      data gamcs( 12) / -.5354219639019687140874081024347d-8     /
      data gamcs( 13) / .9193275519859588946887786825940d-9      /
      data gamcs( 14) / -.1577941280288339761767423273953d-9     /
      data gamcs( 15) / .2707980622934954543266540433089d-10     /
      data gamcs( 16) / -.4646818653825730144081661058933d-11    /
      data gamcs( 17) / .7973350192007419656460767175359d-12     /
      data gamcs( 18) / -.1368078209830916025799499172309d-12    /
      data gamcs( 19) / .2347319486563800657233471771688d-13     /
      data gamcs( 20) / -.4027432614949066932766570534699d-14    /
      data gamcs( 21) / .6910051747372100912138336975257d-15     /
      data gamcs( 22) / -.1185584500221992907052387126192d-15    /
      data gamcs( 23) / .2034148542496373955201026051932d-16     /
      data gamcs( 24) / -.3490054341717405849274012949108d-17    /
      data gamcs( 25) / .5987993856485305567135051066026d-18     /
      data gamcs( 26) / -.1027378057872228074490069778431d-18    /
      data gamcs( 27) / .1762702816060529824942759660748d-19     /
      data gamcs( 28) / -.3024320653735306260958772112042d-20    /
      data gamcs( 29) / .5188914660218397839717833550506d-21     /
      data gamcs( 30) / -.8902770842456576692449251601066d-22    /
      data gamcs( 31) / .1527474068493342602274596891306d-22     /
      data gamcs( 32) / -.2620731256187362900257328332799d-23    /
      data gamcs( 33) / .4496464047830538670331046570666d-24     /
      data gamcs( 34) / -.7714712731336877911703901525333d-25    /
      data gamcs( 35) / .1323635453126044036486572714666d-25     /
      data gamcs( 36) / -.2270999412942928816702313813333d-26    /
      data gamcs( 37) / .3896418998003991449320816639999d-27     /
      data gamcs( 38) / -.6685198115125953327792127999999d-28    /
      data gamcs( 39) / .1146998663140024384347613866666d-28     /
      data gamcs( 40) / -.1967938586345134677295103999999d-29    /
      data gamcs( 41) / .3376448816585338090334890666666d-30     /
      data gamcs( 42) / -.5793070335782135784625493333333d-31    /
!
      data pi / 3.14159265358979323846264338327950d0 /
! sq2pil is 0.5*alog(2*pi) = alog(sqrt(2*pi))
      data sq2pil / 0.91893853320467274178032973640562d0 /
      data ngam, xmin, xmax, xsml, dxrel / 0, 4*0.d0 /
!
      if (ngam.ne.0) go to 10
      ngam = initds (gamcs, 42, 0.1*sngl(d1mach(3)) )
!
      call d9gaml (xmin, xmax)
      xsml = dexp (dmax1 (dlog(d1mach(1)), -dlog(d1mach(2)))+0.01d0)
      dxrel = dsqrt (d1mach(4))
!
 10   y = dabs(x)
      if (y.gt.10.d0) go to 50
!
! compute gamma(x) for -xbnd .le. x .le. xbnd.  reduce interval and find
! gamma(1+y) for 0.0 .le. y .lt. 1.0 first of all.
!
      n = x
      if (x.lt.0.d0) n = n - 1
      y = x - dble(float(n))
      n = n - 1
      dgamma = 0.9375d0 + dcsevl (2.d0*y-1.d0, gamcs, ngam)
      if (n.eq.0) return
!
      if (n.gt.0) go to 30
!
! compute gamma(x) for x .lt. 1.0
!
      n = -n
!
      do 20 i=1,n
        dgamma = dgamma/(x+dble(float(i-1)) )
 20   continue
      return
!
! gamma(x) for x .ge. 2.0 and x .le. 10.0
!
 30   do 40 i=1,n
        dgamma = (y+dble(float(i))) * dgamma
 40   continue
      return
!
! gamma(x) for dabs(x) .gt. 10.0.  recall y = dabs(x).
!
 50   continue 
!
      dgamma = 0.d0
      if (x.lt.xmin) return
!
      dgamma = dexp ((y-0.5d0)*dlog(y) - y + sq2pil + d9lgmc(y) )
      if (x.gt.0.d0) return
!
      sinpiy = dsin (pi*y)
!
      dgamma = -pi/(y*sinpiy*dgamma)
!
      return
      end

!**********************************************************************
!**********************************************************************
      double precision function dint (x)
! december 1983 edition. w. fullerton, c3, los alamos scientific lab.
!
! dint is the double precision equivalent of aint.  this portable
! version is quite efficient when the argument is reasonably small (a
! common case), and so no faster machine-dependent version is needed.
!
      double precision x, xscl, scale, xbig, xmax, part, d1mach, &
        dlog
      external d1mach, dlog, i1mach, r1mach
      data npart, scale, xbig, xmax / 0, 3*0.0d0 /
!
      if (npart.ne.0) go to 10
      ibase = i1mach(10)
      xmax = 1.0d0/d1mach (4)
      xbig = amin1 (float (i1mach(9)), 1.0/r1mach(4))
      scale = ibase**int(dlog(xbig)/dlog(dble(float(ibase)))-0.5d0)
      npart = dlog(xmax)/dlog(scale) + 1.0d0
!
 10   if (x.lt.(-xbig) .or. x.gt.xbig) go to 20
!
      dint = int(sngl(x))
      return
!
 20   xscl = dabs(x)
      if (xscl.gt.xmax) go to 50
!
      do 30 i=1,npart
        xscl = xscl/scale
 30   continue
!
      dint = 0.0d0
      do 40 i=1,npart
        xscl = xscl*scale
        ipart = xscl
        part = ipart
        xscl = xscl - part
        dint = dint*scale + part
 40   continue
!
      if (x.lt.0.0d0) dint = -dint
      return
!
 50   continue 
      dint = x
      return
!
      end

!**********************************************************************
!**********************************************************************
      double precision function dlog (x)
! june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
      double precision x, alncs(11), center(4), alncen(5), aln2, y, t, &
        t2, xn,  dcsevl, d1mach
      external d1mach, dcsevl, initds
!
! series for aln        on the interval  0.          to  3.46021e-03
!                                        with weighted error   4.15e-32
!                                         log weighted error  31.38
!                               significant figures required  31.21
!                                    decimal places required  31.90
!
      data alncs(  1) / .13347199877973881561689386047187d+1     /
      data alncs(  2) / .69375628328411286281372438354225d-3     /
      data alncs(  3) / .42934039020450834506559210803662d-6     /
      data alncs(  4) / .28933847795432594580466440387587d-9     /
      data alncs(  5) / .20512517530340580901741813447726d-12    /
      data alncs(  6) / .15039717055497386574615153319999d-15    /
      data alncs(  7) / .11294540695636464284521613333333d-18    /
      data alncs(  8) / .86355788671171868881946666666666d-22    /
      data alncs(  9) / .66952990534350370613333333333333d-25    /
      data alncs( 10) / .52491557448151466666666666666666d-28    /
      data alncs( 11) / .41530540680362666666666666666666d-31    /
!
      data center(1) / 1.0d0 /
      data center(2) / 1.25d0 /
      data center(3) / 1.50d0 /
      data center(4) / 1.75d0 /
!
      data alncen(  1) / 0.0d0                                         /
      data alncen(  2) / .22314355131420975576629509030983d0     /
      data alncen(  3) / .40546510810816438197801311546434d0     /
      data alncen(  4) / .55961578793542268627088850052682d0     /
      data alncen(  5) / .69314718055994530941723212145817d0     /
!
! aln2 = alog(2.0) - 0.625
      data aln2 / 0.06814718055994530941723212145818d0 /
      data nterms / 0 /
!
      if (nterms.eq.0) nterms = initds (alncs, 11, 28.9*sngl(d1mach(3)))
!
      call d9upak (x, y, n)
!
      xn = n - 1
      y = 2.0d0*y
      ntrval = 4.0d0*y - 2.5d0
!
      if (ntrval.eq.5) t = ((y-1.0d0)-1.0d0) / (y+2.0d0)
      if (ntrval.lt.5) t = (y-center(ntrval)) / (y+center(ntrval))
      t2 = t*t
      dlog = 0.625d0*xn + (aln2*xn + alncen(ntrval) + 2.0d0*t + &
        t*t2*dcsevl(578.d0*t2-1.0d0, alncs, nterms) )
!
      return
      end

!**********************************************************************
!**********************************************************************
      double precision function dsin (x)
! august 1980 edition.  w. fullerton, los alamos scientific lab.
!
! this routine is based on the algorithm of cody and waite in
! argonne tm-321, software manual working note number 1
!
      double precision x, sincs(15), pihi, pilo, pirec, pi2rec, xsml, &
        xwarn, xmax, y, xn, sgn, f, dint, dcsevl, d1mach, dsqrt
      external d1mach, dcsevl, dint, dsqrt, initds
!
! series for sin    on the interval  0.00000e+00 to  2.46740e+00
!                                        with weighted error   2.56e-34
!                                         log weighted error  33.59
!                               significant figures required  33.01
!                                    decimal places required  34.18
!
      data sincs(  1) / -0.374991154955873175839919279977323464d0/
      data sincs(  2) / -0.181603155237250201863830316158004754d0/
      data sincs(  3) /  0.005804709274598633559427341722857921d0/
      data sincs(  4) / -0.000086954311779340757113212316353178d0/
      data sincs(  5) /  0.000000754370148088851481006839927030d0/
      data sincs(  6) / -0.000000004267129665055961107126829906d0/
      data sincs(  7) /  0.000000000016980422945488168181824792d0/
      data sincs(  8) / -0.000000000000050120578889961870929524d0/
      data sincs(  9) /  0.000000000000000114101026680010675628d0/
      data sincs( 10) / -0.000000000000000000206437504424783134d0/
      data sincs( 11) /  0.000000000000000000000303969595918706d0/
      data sincs( 12) / -0.000000000000000000000000371357734157d0/
      data sincs( 13) /  0.000000000000000000000000000382486123d0/
      data sincs( 14) / -0.000000000000000000000000000000336623d0/
      data sincs( 15) /  0.000000000000000000000000000000000256d0/
!
! pihi + pilo = pi.  pihi is exactly representable on all machines
! with at least 8 bits of precision.  whether it is exactly
! represented depends on the compiler.  this routine is more
! accurate if it is exactly represented.
      data pihi   / 3.140625d0 /
      data pilo   / 9.6765358979323846264338327950288d-4/
      data pirec  / 0.31830988618379067153776752674503d0 /
      data pi2rec / 0.63661977236758134307553505349006d0 /
      data ntsn, xsml, xwarn, xmax / 0, 3*0.0d0 /
!
      if (ntsn.ne.0) go to 10
      ntsn = initds (sincs, 15, 0.1*sngl(d1mach(3)))
!
      xsml = dsqrt (2.0d0*d1mach(3))
      xmax = 1.0d0/d1mach(4)
      xwarn = dsqrt (xmax)
!
 10   y = dabs (x)
!
      dsin = x
      if (y.lt.xsml) return
!
      xn = dint (y*pirec+0.5d0)
      n2 = dmod (xn, 2.0d0) + 0.5d0
      sgn = x
      if (n2.ne.0) sgn = -sgn
      f = (y-xn*pihi) - xn*pilo
!
      dsin = f + f*dcsevl(2.0d0*(f*pi2rec)**2-1.0d0, sincs, ntsn)
      if (sgn.lt.0.0d0) dsin = -dsin
      if (dabs(dsin).gt.1.0d0) dsin = dsign (1.0d0, dsin)
!
      return
      end

!**********************************************************************
!**********************************************************************
      double precision function dsqrt (x)
! june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
      double precision x, sqrt2(3), y,  d9pak, d1mach
      external alog, d1mach, d9pak
      data sqrt2(1) / 0.70710678118654752440084436210485d0 /
      data sqrt2(2) / 1.0d0 /
      data sqrt2(3) / 1.41421356237309504880168872420970d0 /
!
      data niter / 0 /
!
      if (niter.eq.0) niter = 1.443*alog(-0.104*alog(0.1*sngl(d1mach(3)))) + 1.0
!
      if (x.le.0.d0) go to 20
!
      call d9upak (x, y, n)
      ixpnt = n/2
      irem = n - 2*ixpnt + 2
!
! the approximation below has accuracy of 4.16 digits.
      z = y
      dsqrt = .261599e0 + z*(1.114292e0 + z*(-.516888e0 + z*.141067e0))
!
      do 10 iter=1,niter
        dsqrt = dsqrt + 0.5d0*(y - dsqrt*dsqrt) / dsqrt
 10   continue
!
      dsqrt = d9pak (sqrt2(irem)*dsqrt, ixpnt)
      return
!
 20   continue
      dsqrt = 0.0d0
      return
!
      end

!**********************************************************************
!**********************************************************************
      integer function i8save(isw,ivalue,set)
!
!  if (isw = 1) i8save returns the current error number and
!               sets it to ivalue if set = .true. .
!
!  if (isw = 2) i8save returns the current recovery switch and
!               sets it to ivalue if set = .true. .
!
      logical set
!
      integer iparam(2)
!  iparam(1) is the error number and iparam(2) is the recovery switch.
!
!  start execution error free and with recovery turned off.
!
      data iparam(1) /0/,  iparam(2) /2/
!
      i8save=iparam(isw)
      if (set) iparam(isw)=ivalue
!
      return
!
      end

!**********************************************************************
!**********************************************************************
      function initds (dos, nos, eta)
! june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
!
! initialize the double precision orthogonal series dos so that initds
! is the number of terms needed to insure the error is no larger than
! eta.  ordinarily eta will be chosen to be one-tenth machine precision.
!
!             input arguments --
! dos    dble prec array of nos coefficients in an orthogonal series.
! nos    number of coefficients in dos.
! eta    requested accuracy of series.
!
      double precision dos(nos)
!
      err = 0.
      do 10 ii=1,nos
        i = nos + 1 - ii
        err = err + abs(sngl(dos(i)))
        if (err.gt.eta) go to 20
 10   continue
!
 20   continue 
      initds = i
!
      return
      end
      
!======================================================================
      function alog (x)
!======================================================================
! series for aln on the interval 0. to 3.46021d-03
!                                        with weighted error   1.50e-16
!                                         log weighted error  15.82
!                               significant figures required  15.65
!                                    decimal places required  16.21
! History:
!   june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
!----------------------------------------------------------------------
      dimension alncs(6), center(4), alncen(5)
      external csevl, inits, r1mach
!
      data alncs( 1) /   1.3347199877973882e0   /
      data alncs( 2) /    .000693756283284112e0 /
      data alncs( 3) /    .000000429340390204e0 /
      data alncs( 4) /    .000000000289338477e0 /
      data alncs( 5) /    .000000000000205125e0 /
      data alncs( 6) /    .000000000000000150e0 /
!
      data center(1) / 1.0 /
      data center(2) / 1.25 /
      data center(3) / 1.50 /
      data center(4) / 1.75 /
!
      data alncen(  1) / 0.0e0                                    /
      data alncen(  2) /  .223143551314209755e0                   /
      data alncen(  3) /  .405465108108164381e0                   /
      data alncen(  4) /  .559615787935422686e0                   /
      data alncen(  5) /  .693147180559945309e0                   /
!
      data aln2 / 0.068147180559945309e0 /
      data nterms / 0 /
!
      if (nterms.eq.0) nterms = inits (alncs, 6, 28.9*r1mach(3))
!
      call r9upak (x, y, n)
!
      xn = n - 1
      y = 2.0*y
      ntrval = 4.0*y - 2.5
      if (ntrval.eq.5) t = ((y-1.0)-1.0) / (y+2.0)
      if (ntrval.lt.5) t = (y-center(ntrval))/(y+center(ntrval))
      t2 = t*t
!
      alog = 0.625*xn + (aln2*xn + alncen(ntrval) + 2.0*t + &
        t*t2*csevl(578.0*t2-1.0, alncs, nterms) )
!
      return
      end
