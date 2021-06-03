!=========================================================================
     program scan
!=========================================================================
! calculates output for various subroutines by scanning through ranges of input variables
!
! Modifications:
! 2/1/2000  Kevin Schaefer created scan
! 7/1/2001  Kevin Schaefer added generic scanning
! 1/17/2021 Kevin Schaefer split off from scan and restructures for MRT
! 1/30/2021 Kevin Schaefer got rid of biome loop and if_scan
! 4/4/2021  Kevin Schaefer moved write code to separate routine
!--------------------------------------------------------------------------

  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output

  IMPLICIT NONE

  integer ix   ! x value index
  integer iy   ! y value index
  integer ilin ! line index
  integer iseg ! segment index

! variables for timing
  integer clock_rate   ! conversion between count and clock time
  integer clock_start  ! time count at stat
  integer clock_stop   ! time count at end
  real(8) e_time       ! elapsed time

! allocate sib variable tree
  call clock

! read inputs
  call scan_read_input

! start measuring execution time
  call system_clock(clock_start) ! start counting

! set start values X variable
  call Plot2D3D('StartX  ')
  do ix=1,maxi
!
! set start values Y variables
    if (flagDim==3) call Plot2D3D('StartY  ')
    do iy=1,maxj

! call scan control
      call scan_control()

! assign plotting variables
      call Plot2D3D('AssXval ')
      X(1,ix)=Xval
      if(Flagdim==2) then
        do ilin=1,nline
          X(ilin,ix)=Xval
          Y(ilin,ix)=test10(ilin)
        enddo
        print*, trim(Label(xvar)), value(xvar), test10(1), Tbo_obs(1,:)
      endif
      if(Flagdim==3) then
        zval(1,ix,iy)=test
        zval(2,ix,iy)=testvar1
        x3(ix)=value(xvar)
        y3(iy)=value(yvar)
        print*, trim(Label(xvar)), value(xvar), trim(Label(yvar)), value(yvar), test
      endif

! increment Y variable
      if (flagDim==3) call Plot2D3D('IncrY   ')
    enddo ! end y-loop

! Increment X variable
    call Plot2D3D('IncrX   ')
  enddo ! end x-loop

! write Output
  call write_output

! print execution times
  if(print_ex) then
    do iseg = 1, nseg
      print*, seg_name(iseg), time_seg(iseg), num_call(iseg)
    enddo
  endif
  call system_clock(clock_stop,clock_rate) ! stop counting
  e_time = real(clock_stop-clock_start)/real(clock_rate)
  print*, 'Overall Execution time (s):', e_time

!--------------------------------------------------
! INTERNAL SUBROUTINES
!---------------------------------------------------
    contains
!
!---------------------------------------------------
    subroutine clock
!---------------------------------------------------
! Call system clock again to determine how long code took to run. 
! The values of Tstart and Tend are in internal system time units. 
! The variable specifier count_rate provides the conversion from 
! internal system time units to seconds. 
!value(xvar), trim(Label(yvar)), value(yvar)
    call system_clock(count_rate=clock2)
    call system_clock(count=Tend)
    if(Tend >= Tstart) then
               Time=float(Tend-Tstart)/float(clock2)
    else
               Time=0.0
    endif
!
    end subroutine clock

    end program scan
!
!=======================================================================
    subroutine output_manipulation
!=========================================================================
! manipulates output from scan subroutine
!
! Modifications:
!  Kevin Schaefer made routine (6/8/14))
!--------------------------------------------------------------------------
!
    use scan_Variables
!
  IMPLICIT NONE
!
    real(8) var1, var2 ! local temporary variables
    integer iman ! manipulation index
    integer ilin ! line number index
    integer iptx  ! x point number index
    integer ipty  ! y point number index
    integer ibio ! biome index
!
! output manipulations
    do iman=1,nout_man
!
! 2-D value Integral 
      if(out_man(iman)%doit.and.flagDim==2.and.out_man(iman)%typ=='sum') then
        print*, '2-D Integral sum'
        var2=(xmax-xmin)/real(numpts)
     do ilin=1,nline
       var1=0.
       do iptx=1,Numpts
         var1=var1+Y(ilin,iptx)
       enddo
       var1=var1*var2
       print*, trim(legend(ilin)),': ',var1
     enddo
      endif
!
! 3-D value Integral nv, MinBiome,MaxBiome
      if(out_man(iman)%doit.and.flagDim==3.and.out_man(iman)%typ=='sum') then
        print*, '3-D Integral sum'
          var1=0.
          do ipty=1,Numpts
            do iptx=1,Numpts
              var1=var1+zval(ibio,iptx,ipty)
            enddo
          enddo
          var1=var1*dx*dy
          print*, 'biome ',ibio,' sum: ',var1
      endif
!
! end manipulation loop
    enddo
!
    return
    end
!
!=======================================================================
    subroutine write_output
!=========================================================================
! writes output from scan program
!
! Modifications:
!  04/04/2021 Kevin Schaefer created routine
!--------------------------------------------------------------------------
  use scan_Variables
  implicit none

! Internal variables
  integer ix   ! x value index
  integer iy   ! y value index
  integer ipt  ! point index
  integer i3d  ! 3d index
  Character*100 form       ! writing format
  Character*500 write_line ! writing line
  Character*20 var_char    ! character version of variable

!--------------------------------------------------------------------------
! save line data as formatted text file
!--------------------------------------------------------------------------
  if(save_2d==1) then
    open(87, file=trim(plotfile)//'_'//trim(scantype)//'.dat', form='formatted')

! header
    write(form,'(i2)') nline
    form='(a15,1x,'//trim(form)//'(a15,1x))'
    Write(87,form) labx,legend(1:nline)

! main body
    write(form,'(i2)') nline
    form='(e15.8,1x,'//trim(form)//'(e15.8,1x))'
    do ipt=1,Numpts
      Write(87,form) x(1,ipt),Y(1:nline,ipt)
    enddo
    close(unit=87)
  endif

!--------------------------------------------------------------------------
! Save line data as CSV file
!--------------------------------------------------------------------------
  if(save_2d==2) then
    open(87, file=trim(plotfile)//'_'//trim(scantype)//'.csv', form='formatted')
	
! header
    write_line=trim(labx)
    do iy=1,nline
      write_line=trim(adjustl(write_line))//','//trim(adjustl(legend(iy)))
    enddo
    Write(87,*) adjustl(write_line)

! main body
    form='(e15.8,1x)'
    do ipt=1,Numpts
      write(var_char, '(e15.8)') x(1,ipt)
      write_line=trim(adjustl(var_char))
      do iy=1,nline
        write(var_char, '(e15.8)') Y(iy,ipt)
        write_line=trim(adjustl(write_line))//','//trim(adjustl(var_char))
      enddo
      Write(87,*) adjustl(write_line)
    enddo

    close(unit=87)
  endif

!--------------------------------------------------------------------------
! Save contour data as CSV file
!--------------------------------------------------------------------------
  if(flagDim==3 .and. save_3d==2) then
    do i3d = 1,2
      if(i3d == 1) open(87, file=trim(plotfile)//'_'//trim(scantype)//'_3d_01.csv', form='formatted')
      if(i3d == 2) open(87, file=trim(plotfile)//'_'//trim(scantype)//'_3d_02.csv', form='formatted')

! header
      write(var_char,'(i2)') Numpts
      write_line=trim(adjustl(var_char))
      form='(e15.8,1x)'
      do ix=1,Numpts
        write(var_char, '(e15.8)') x3(ix)
        write_line=trim(adjustl(write_line))//','//trim(adjustl(var_char))
      enddo
      Write(87,*) adjustl(write_line)

! main body
      form='(e15.8,1x)'
      do iy=1,Numpts
        write(var_char, '(e15.8)') y3(iy)
        write_line=trim(adjustl(var_char))
        do ix=1,Numpts
          write(var_char, '(e15.8)') zval(i3d,ix,iy)
          write_line=trim(adjustl(write_line))//','//trim(adjustl(var_char))
        enddo
        Write(87,*) adjustl(write_line)
      enddo

      close(unit=87)
    enddo
  endif

  return
  end
