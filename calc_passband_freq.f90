!===========================================================
subroutine calc_passband_freq()
!===========================================================
! This rountine:
!     1 - calculates passband quadrature frequency values
!         from channel specification
!     2 - chan_spec is global in variables_unit
! History:
!     1 - Original provided in PASCAL by Marian Klein and Albin Gasiewski
!         NOAA ETL/ET1 Microwave System Development Division
!     2 - Converted June 2003 from PASCAL to FORTRAN by Bill Otto
!         william.d.otto@noaa.gov NOAA/ETL SET
!         FORTRAN 90 Portland Group Compiler on Red Hat Linux
!     3 - Modified 31 July 2003 by Bob Weber
!         FORTRAN 90 COMPAQ Compiler on Windows 2000 / XP
!      9/26/2020 Kevin Schaefer deleted unused variables 
!----------------------------------------------------------
  use dotlrt_variables
    
  implicit none

  real(8) lo_freq   ! (GHz) LO frequency
  real(8) if1_freq  ! (GHz) IF1 frequency
  real(8) if2_freq  ! (GHz) IF2 frequency
  real(8) bandwidth ! (GHz) filter bandwidth
  real(8) dtrms     ! (K) observation noise in degrees Kelvin
  real(8) freq_incr ! (GHz) delta frequency in passband integral
  real(8) freq      ! (GHz) frequency in passband integral
  integer ifrq      ! (-) frequency index

    lo_freq   = channel%lo_freq
    if1_freq  = channel%if1_freq
    if2_freq  = channel%if2_freq
    bandwidth = channel%bandwidth
    dtrms     = channel%dtrms

    ! if invalid arguments, then
    num_freqs = 0
    if( (nsub_freq .le. 0) .or. (lo_freq .le. 0.0d0) .or.  &
        (if1_freq .lt. 0.0d0) .or. (if2_freq .lt. 0.0d0) .or. &
        (bandwidth .lt. 0.0d0) ) return
    if( if1_freq .eq. 0.0d0) then
        ! Single sideband or single conversion superhet receiver with
        ! only low pass IF filtering : one contigous passband
        if( ( (lo_freq-bandwidth/2.0d0) .gt. 0.0d0) .and. &
            ( (lo_freq+bandwidth/2.0d0) .lt. max_freq) ) then
            if( nsub_freq .gt. 1 ) then
                if( nsub_freq .gt. max_nfreq ) then
                    freq_incr = bandwidth/dble(max_nfreq-1)
                    freq = lo_freq-bandwidth/2.0d0
                    do ifrq = 1, max_nfreq 
                        passband_freq(ifrq) = freq
                        freq = freq+freq_incr
                    end do
                    num_freqs = max_nfreq
                    ! WARNING: Requested number of quadrature frequencies
                    !          not available. Maximum number being used
                else
                    freq_incr = bandwidth/dble(nsub_freq-1)
                    freq = lo_freq-bandwidth/2.0d0
                    do ifrq = 1, nsub_freq
                        passband_freq(ifrq) = freq
                        freq = freq+freq_incr
                    end do
                    num_freqs = nsub_freq
                end if
            else
                passband_freq(1) = lo_freq
                num_freqs = 1
            end if
        else
            num_freqs = 0
            ! ERROR: Channel frequency out of range
        end if
    else ! if1_freq
        if( if2_freq .eq. 0.0d0 ) then
            ! Double sideband or double conversion superhet receiver with only
            ! low pass 2nd IF filtering one or two contigous passbands
            if( ((lo_freq-if1_freq-bandwidth/2.0d0) .gt. 0.0d0) .and.    &
                ((lo_freq+if1_freq+bandwidth/2.0d0) .lt. max_freq) .and. &
                (if1_freq .gt. bandwidth/2.0d0) .and.                    &
                (lo_freq .gt. if1_freq) ) then
                if( nsub_freq .gt. 1 ) then
                    if( nsub_freq .gt. idint(dble(max_nfreq)/2.0d0)) THEN
                        freq_incr = bandwidth/(idint(dble(max_nfreq)/2.0d0)-1.0d0)
                        freq = lo_freq-if1_freq-bandwidth/2.0d0
                        do ifrq = 1, idint(dble(max_nfreq)/2.0d0) 
                            passband_freq(ifrq) = freq
                            freq = freq+freq_incr
                        end do
                        freq = lo_freq+if1_freq-bandwidth/2.0d0
                        do ifrq = idint(dble(max_nfreq)/2.0d0)+1, &
                               2*idint(dble(max_nfreq)/2.0d0)
                            passband_freq(ifrq) = freq
                            freq = freq+freq_incr
                        end do
                        num_freqs = 2*idint(dble(max_nfreq)/2.0d0)
                        ! WARNING: Requested number of quadrature frequencies not available.
                        !          Maximum number being used
                    else
                        freq_incr = bandwidth/dble(nsub_freq-1)
                        freq = lo_freq-if1_freq-bandwidth/2.0d0
                        do ifrq = 1, nsub_freq
                            passband_freq(ifrq) = freq
                            freq = freq+freq_incr
                        end do
                        freq = lo_freq+if1_freq-bandwidth/2.0d0
                        do ifrq = nsub_freq+1, 2*nsub_freq
                            passband_freq(ifrq) = freq
                            freq = freq+freq_incr
                        end do
                        num_freqs = 2*nsub_freq
                    end if
                else
                    passband_freq(1) = lo_freq-if1_freq
                    passband_freq(2) = lo_freq+if1_freq
                    num_freqs = 2
                end if
            else
                num_freqs = 0
                ! ERROR: Channel frequency out of range
            end if
        else ! if2_freq
            ! Double conversion superhet receiver with passband filtering in 2nd IF
            ! one, two, or four contigous passbands
            if( ((lo_freq-if1_freq-if2_freq-bandwidth/2.0d0) .gt. 0.0d0) .and.       &
                ((lo_freq+if1_freq+if2_freq+bandwidth/2.0d0) .lt. max_freq) .and.    &
                (if2_freq .gt. bandwidth/2.0d0) .and. (if1_freq .gt. if2_freq) .and. &
                (lo_freq .gt. if1_freq) ) then
                if( nsub_freq .gt. 1 ) then
                    if( nsub_freq .gt. idint(dble(max_nfreq)/4.0d0) ) then
                        freq_incr = bandwidth/(idint(dble(max_nfreq)/4.0d0)-1)
                        freq = lo_freq-if1_freq-if2_freq-bandwidth/2.0d0
                        do ifrq = 1, idint(dble(max_nfreq)/4.0d0)
                            passband_freq(ifrq) = freq
                            freq = freq+freq_incr
                        end do
                        freq = lo_freq-if1_freq+if2_freq-bandwidth/2.0d0
                        do ifrq = idint(dble(max_nfreq)/4.0d0)+1, &
                               2*idint(dble(max_nfreq)/4.0d0)
                            passband_freq(ifrq) = freq
                            freq = freq+freq_incr
                        end do
                        freq = lo_freq+if1_freq-if2_freq-bandwidth/2.0
                        do ifrq = 2*idint(dble(max_nfreq)/4.0d0)+1, &
                               3*idint(dble(max_nfreq)/4.0d0)
                            passband_freq(ifrq) = freq
                            freq = freq+freq_incr
                        end do
                        freq = lo_freq+if1_freq+if2_freq-bandwidth/2.0d0
                        do ifrq = 3*idint(dble(max_nfreq)/4.0d0)+1, &
                               4*idint(dble(max_nfreq)/4.0d0)
                            passband_freq(ifrq) = freq
                            freq = freq+freq_incr
                        end do
                        num_freqs = 4*idint(dble(max_nfreq)/4.0d0)
                        ! WARNING: Requested number of quadrature frequencies not available.
                        !          Maximum number being used
                    else
                        freq_incr = bandwidth/dble(nsub_freq-1)
                        freq = lo_freq-if1_freq-if2_freq-bandwidth/2.0d0
                        do ifrq = 1, nsub_freq
                            passband_freq(ifrq) = freq
                            freq = freq+freq_incr
                        end do
                        freq = lo_freq-if1_freq+if2_freq-bandwidth/2.0d0
                        do ifrq = nsub_freq+1, 2*nsub_freq
                            passband_freq(ifrq) = freq
                            freq = freq+freq_incr
                        end do
                        freq = lo_freq+if1_freq-if2_freq-bandwidth/2.0d0
                        do ifrq = 2*nsub_freq+1, 3*nsub_freq
                            passband_freq(ifrq) = freq
                            freq = freq+freq_incr
                        end do
                        freq = lo_freq+if1_freq+if2_freq-bandwidth/2.0
                        do ifrq = 3*nsub_freq+1, 4*nsub_freq
                            passband_freq(ifrq) = freq
                            freq = freq+freq_incr
                        end do
                        num_freqs = 4*nsub_freq
                    end if
                else
                    passband_freq(1) = lo_freq-if1_freq-if2_freq
                    passband_freq(2) = lo_freq-if1_freq+if2_freq
                    passband_freq(3) = lo_freq+if1_freq-if2_freq
                    passband_freq(4) = lo_freq+if1_freq+if2_freq
                    num_freqs = 4
                end if
            else
                num_freqs = 0
                ! ERROR: Channel frequency out of range
            end if
        end if ! if2_freq
    end if ! if1_freq
    return
end subroutine calc_passband_freq
