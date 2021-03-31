! mex_matlab_mrt.f90
! Matlab gateway routine to mrt
subroutine mexFunction(nlhs, plhs, nrhs, prhs)
  implicit none
!---------------------------------------------------------
  integer*4 plhs(*), prhs(*)
  integer*4 mxCreateDoubleMatrix, mxCreateNumericArray
  integer*4 y1_pr, y2_pr, y3_pr, y4_pr, y5_pr, y6_pr, y7_pr, &
            y8_pr, y9_pr, y0_pr
  integer*4 ndim, dims(8), classid, ComplexFlag
  integer*4 mxClassIDFromClassName
  character*8 classname
  integer*4 z1_pr, z2_pr, z3_pr, z4_pr, z5_pr, z6_pr, z7_pr, &
            z8_pr, z9_pr, z0_pr
!---------------------------------------------------------
  integer(4) nlhs, nrhs
  integer, parameter :: max_num_levels = 1024
  integer, parameter :: max_number_streams = 32
  integer, parameter :: max_nsurfangles = 32

  integer(4) mxGetM, mxGetN, mxIsNumeric, mxGetPr
  integer(4) nang, nsb, nsurfangles, nlevels, i, j, k, jud, jpol, alloc_err
  real(8) num_streams, num_surf_angles, num_sb_freqs, nz
  real(8) inp_height, inp_theta
  real(8) instr_spec(5)
  real(8) atm_inp(max_num_levels*9), surf_inp(max_nsurfangles*3)
  real(8), dimension((max_num_levels+1)*2) :: Tb
  real(8), dimension(max_num_levels*2) :: dTb_dT, dTb_dp, dTb_dq
  real(8), dimension(max_num_levels*5*2) :: dTb_dw
  real(8), dimension(2) :: Tb_inp
  real(8), dimension(2) :: tau
  real(8), dimension(max_number_streams) :: stream_angles
  real(8), dimension(:,:), allocatable :: atminp, surfinp
  real(8), dimension(:,:), allocatable :: aTb, adTb_dT, adTb_dp, adTb_dq
  real(8), dimension(:,:,:), allocatable :: adTb_dw


!     Check for proper number of arguments. 
      if (nrhs .ne. 9) then
         call mexErrMsgTxt('9 input required.')
      elseif (nlhs .ne. 8) then
         call mexErrMsgTxt('8 outputs required.')
      endif

      z1_pr = mxGetPr(prhs( 1))
      z2_pr = mxGetPr(prhs( 2))
      z3_pr = mxGetPr(prhs( 3))
      z4_pr = mxGetPr(prhs( 4))
      z5_pr = mxGetPr(prhs( 5))
      z6_pr = mxGetPr(prhs( 6))
      z7_pr = mxGetPr(prhs( 7))
      z8_pr = mxGetPr(prhs( 8))
      z9_pr = mxGetPr(prhs( 9))

!     Load the input into a Fortran array.
      ! observation parameters
      call mxCopyPtrToReal8(z1_pr, inp_height, 1)
      call mxCopyPtrToReal8(z2_pr, inp_theta, 1)
      ! instrument parameters
      call mxCopyPtrToReal8(z3_pr, num_sb_freqs, 1)
      nsb                 = num_sb_freqs
      call mxCopyPtrToReal8(z4_pr, instr_spec, 5)
      ! surface parameters
      call mxCopyPtrToReal8(z5_pr, num_streams, 1)
      nang              = num_streams
      call mxCopyPtrToReal8(z6_pr, num_surf_angles, 1)
      nsurfangles              = num_surf_angles
      call mxCopyPtrToReal8(z7_pr, surf_inp, nsurfangles*3)
      ! atmosphere parameters
      call mxCopyPtrToReal8(z8_pr, nz, 1)
      nlevels = nz - 1    ! first level at surface
      call mxCopyPtrToReal8(z9_pr, atm_inp, (nlevels+1)*9)

  allocate(atminp(nlevels+1,9), surfinp(nsurfangles,3), &
           aTb(0:nlevels,2), adTb_dT(nlevels,2),  &
           adTb_dp(nlevels,2), adTb_dq(nlevels,2), adTb_dw(nlevels,5,2) )
  i = 0
  do k = 1, 9
    do j = 1, nlevels+1
      i = i + 1
      atminp(j,k) = atm_inp(i)
    end do
  end do
  i = 0
  do k = 1, 3
    do j = 1, nsurfangles
      i = i + 1
      surfinp(j,k) = surf_inp(i)
    end do
  end do

!     call fortran routine to read crtm binary data for input
      call mrt( inp_height,                     & ! observation height (AGL, km)
                inp_theta,                      & ! observation angle from nadir (degrees)
                nsb,                            & ! number of frequency points per sideband
                instr_spec,                     &
                nang,                           & ! number of streams (total, up and down)
                stream_angles,                  &
                nsurfangles,                    & ! number of angles for surface reflectivity
                surfinp,                        & ! surface reflectivity (vertical, horizontal polarization)
                                                  ! for specified angles of incidence
                nlevels,                        &
                atminp,                         &
                aTb, adTb_dT, adTb_dp,          &
                adTb_dq, adTb_dw, Tb_inp, tau )

  i = 0
  do jpol = 1, 2
    do k = 1, 5
        do j = 1, nlevels
          i = i + 1
          dTb_dw(i) = adTb_dw(j,k,jpol)
        end do
    end do
  end do

  i = 0
  do jpol = 1, 2
      do j = 1, nlevels
        i = i + 1
        dTb_dT(i) = adTb_dT(j,jpol)
        dTb_dp(i) = adTb_dp(j,jpol)
        dTb_dq(i) = adTb_dq(j,jpol)
      end do
  end do

  i = 0
  do jpol = 1, 2
      do j = 0, nlevels
        i = i + 1
        Tb(i) = aTb(j,jpol)
      end do
  end do

!     Create matrices for the return arguments.
!      classname = 'int16'
      classname = 'double'
      classid = mxClassIDFromClassName(classname)
      ComplexFlag = 0

      ndim = 1
      dims(1) = 2
      plhs(1) = mxCreateNumericArray(ndim, dims, classid, ComplexFlag)

      ndim = 1
      dims(1) = 2
      plhs(2) = mxCreateNumericArray(ndim, dims, classid, ComplexFlag)

      ndim = 2
      dims(1) = nlevels + 1
      dims(2) = 2
      plhs(3) = mxCreateNumericArray(ndim, dims, classid, ComplexFlag)

      do j = 4, 6
        ndim = 2
        dims(1) = nlevels
        dims(2) = 2
        plhs(j) = mxCreateNumericArray(ndim, dims, classid, ComplexFlag)
      end do

      ndim = 3
      dims(1) = nlevels
      dims(2) = 5
      dims(3) = 2
      plhs(7) = mxCreateNumericArray(ndim, dims, classid, ComplexFlag)

      ndim = 1
      dims(1) = nang
      plhs(8) = mxCreateNumericArray(ndim, dims, classid, ComplexFlag)


      y1_pr = mxGetPr(plhs( 1))
      y2_pr = mxGetPr(plhs( 2))
      y3_pr = mxGetPr(plhs( 3))
      y4_pr = mxGetPr(plhs( 4))
      y5_pr = mxGetPr(plhs( 5))
      y6_pr = mxGetPr(plhs( 6))
      y7_pr = mxGetPr(plhs( 7))
      y8_pr = mxGetPr(plhs( 8))

!     Load the data into MatLab arrays.
      call mxCopyReal8ToPtr(Tb_inp,            y1_pr, 2)
      call mxCopyReal8ToPtr(tau,               y2_pr, 2)
      call mxCopyReal8ToPtr(Tb,                y3_pr, (nlevels+1)*2)
      call mxCopyReal8ToPtr(dTb_dT,            y4_pr, nlevels*2)
      call mxCopyReal8ToPtr(dTb_dp,            y5_pr, nlevels*2)
      call mxCopyReal8ToPtr(dTb_dq,            y6_pr, nlevels*2)
      call mxCopyReal8ToPtr(dTb_dw,            y7_pr, nlevels*5*2)
      call mxCopyReal8ToPtr(stream_angles,     y8_pr, nang)

    deallocate( atminp, surfinp, aTb, adTb_dT,  &
           adTb_dp, adTb_dq, adTb_dw, stat = alloc_err)
  return
  end