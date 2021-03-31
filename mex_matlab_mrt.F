! mex_matlab_mrt.f90
! Matlab gateway routine to mrt
#include "fintrf.h"

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      implicit none
      !---------------------------------------------------------
      mwPointer plhs(*), prhs(*)
      mwPointer mxCreateDoubleMatrix, mxCreateNumericArray
      mwPointer y1_pr, y2_pr, y3_pr, y4_pr, y5_pr, y6_pr, y7_pr, y8_pr
      mwPointer y9_pr, y10_pr, y11_pr, y12_pr, y13_pr, y14_pr
      mwPointer y15_pr, y16_pr, y17_pr
      mwSize    ndim, dims(8)
      integer*4 classid
      integer*4 ComplexFlag
      integer*4 mxClassIDFromClassName
      character*8 classname
      mwPointer z1_pr, z2_pr, z3_pr, z4_pr, z5_pr, z6_pr, z7_pr, z8_pr
      mwPointer z9_pr, z0_pr
      !---------------------------------------------------------
      integer  nlhs, nrhs
      integer, parameter :: max_num_levels = 1024
      integer, parameter :: max_number_streams = 32
      integer, parameter :: max_nsurfangles = 32

      mwPointer mxGetM, mxGetN, mxGetPr
      integer mxIsNumeric
      integer*4 nang, nsb, nsurfangles, nlevels, i, j, k, jud, jpol,jj
      integer*4 alloc_err
      real*8 num_streams, num_surf_angles, num_sb_freqs, nz
      real*8 inp_height(1), inp_theta(1)
      real*8 instr_spec(5)
      real*8 atm_inp(max_num_levels*9), surf_inp(max_nsurfangles*3)
      real*8, dimension((max_num_levels+1)*2) :: Tb
      real*8, dimension(max_num_levels*2) :: dTb_dT, dTb_dp, dTb_dq
      real*8, dimension(max_num_levels*5*2) :: dTb_dw
      real*8, dimension(2) :: Tb_inp
      real*8, dimension(2) :: tau
      real*8, dimension(max_number_streams) :: stream_angles
      real*8, dimension(max_number_streams*2) :: Tbo_streams
      real*8, dimension(max_num_levels*max_number_streams*2) :: 
     + dTb_dT_streams, dTb_dp_streams, dTb_dq_streams
      real*8, dimension(max_num_levels*max_number_streams*2) :: 
     + dTb_dhc_streams, dTb_dhr_streams, dTb_dhi_streams
      real*8, dimension(max_num_levels*max_number_streams*2) ::
     + dTb_dhs_streams, dTb_dhg_streams
      real*8, dimension(:,:), allocatable :: atminp, surfinp
      real*8, dimension(:,:), allocatable :: aTb, adTb_dT, adTb_dp
      real*8, dimension(:,:), allocatable :: adTb_dq
      real*8, dimension(:,:,:), allocatable :: adTb_dw
      real*8, dimension(:,:), allocatable :: Tbo_streams_inp
      real*8, dimension(:,:,:), allocatable :: dTb_dT_streams_inp
      real*8, dimension(:,:,:), allocatable :: dTb_dp_streams_inp
      real*8, dimension(:,:,:), allocatable :: dTb_dq_streams_inp
      real*8, dimension(:,:,:,:), allocatable :: dTb_dw_streams_inp

      character*120 debugout

      ! Check for proper number of arguments. 
      if (nrhs .ne. 9) then
         call mexErrMsgTxt('9 input required.')
      elseif (nlhs .ne. 17) then
         call mexErrMsgTxt('17 outputs required.')
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

      ! Load the input into a Fortran array.
      ! observation parameters
      call mxCopyPtrToReal8(z1_pr, inp_height, 1)
      !write(debugout,*) "height = ",inp_height
      !call mexPrintf(debugout//achar(10))

      call mxCopyPtrToReal8(z2_pr, inp_theta, 1)
      !write(debugout,*) "theta = ",inp_theta
      !call mexPrintf(debugout//achar(10))

      ! instrument parameters
      call mxCopyPtrToReal8(z3_pr, num_sb_freqs, 1)
      nsb = num_sb_freqs
      !write(debugout,*) "num_sb =",num_sb_freqs
      !call mexPrintf(debugout//achar(10))

      call mxCopyPtrToReal8(z4_pr, instr_spec, 5)
      !write(debugout,*) "instr_spec(1) = ",instr_spec(1)
      !call mexPrintf(debugout//achar(10))
      !write(debugout,*) "instr_spec(2) = ",instr_spec(2)
      !call mexPrintf(debugout//achar(10))
      !write(debugout,*) "instr_spec(3) = ",instr_spec(3)
      !call mexPrintf(debugout//achar(10))
      !write(debugout,*) "instr_spec(4) = ",instr_spec(4)
      !call mexPrintf(debugout//achar(10))
      !write(debugout,*) "instr_spec(5) = ",instr_spec(5)
      !call mexPrintf(debugout//achar(10))

      ! surface parameters
      call mxCopyPtrToReal8(z5_pr, num_streams, 1)
      nang = num_streams
      !write(debugout,*) "num streams = ",num_streams
      !call mexPrintf(debugout//achar(10))

      call mxCopyPtrToReal8(z6_pr, num_surf_angles, 1)
      nsurfangles = num_surf_angles
      !write(debugout,*) "num_surf_angles = ",nsurfangles
      !call mexPrintf(debugout//achar(10))

      call mxCopyPtrToReal8(z7_pr, surf_inp, nsurfangles*3)
      !write(debugout,*) "surf_inp angle = ",surf_inp(1), surf_inp(8)
      !call mexPrintf(debugout//achar(10))
      !write(debugout,*) "surf_inp Rv = ",surf_inp(9), surf_inp(16)
      !call mexPrintf(debugout//achar(10))
      !write(debugout,*) "surf_inp Rh = ",surf_inp(17), surf_inp(24)
      !call mexPrintf(debugout//achar(10))

      ! atmosphere parameters
      call mxCopyPtrToReal8(z8_pr, nz, 1)
      nlevels = nz - 1    ! first level at surface
      !write(debugout,*) "nlevels = ",nz
      !call mexPrintf(debugout//achar(10))

      call mxCopyPtrToReal8(z9_pr, atm_inp, (nlevels+1)*9)

      !write(debugout,*) "atm_inp height = ",atm_inp(1), atm_inp(61)
      !call mexPrintf(debugout//achar(10))
      !write(debugout,*) "atm_inp pres = ",atm_inp(62), atm_inp(122)
      !call mexPrintf(debugout//achar(10))
      !write(debugout,*) "atm_inp vd = ",atm_inp(184), atm_inp(244)
      !call mexPrintf(debugout//achar(10))
      !write(debugout,*) "atm_inp density = ",atm_inp(488), atm_inp(548)
      !call mexPrintf(debugout//achar(10))

      allocate(atminp(nlevels+1,9))
      allocate(surfinp(nsurfangles,3))
      allocate(aTb(0:nlevels,2))
      allocate(adTb_dT(nlevels,2))
      allocate(adTb_dp(nlevels,2))
      allocate(adTb_dq(nlevels,2))
      allocate(adTb_dw(nlevels,5,2))
      allocate(Tbo_streams_inp(nang/2,2))
      allocate(dTb_dT_streams_inp(nlevels,nang/2,2))
      allocate(dTb_dp_streams_inp(nlevels,nang/2,2))
      allocate(dTb_dq_streams_inp(nlevels,nang/2,2))
      allocate(dTb_dw_streams_inp(nlevels,nang/2,5,2))

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

      ! call fortran routine to read crtm binary data for input
      ! inp_height,      observation height (AGL, km)
      ! inp_theta,       observation angle from nadir (degrees)
      ! nsb,             number of frequency points per sideband
      ! instr_spec,      
      ! nang,            number of streams (total, up and down)
      ! stream_angles,   
      ! nsurfangles,     number of angles for surface reflectivity
      ! surfinp,         surface reflectivity (vertical, horizontal polarization)
      !                  for specified angles of incidence

      call mrt( inp_height, inp_theta, nsb, instr_spec, nang, 
     +          stream_angles, nsurfangles, surfinp, nlevels, 
     +          atminp, aTb, adTb_dT, adTb_dp,
     +          adTb_dq, adTb_dw, Tb_inp, tau, Tbo_streams_inp,
     +          dTb_dT_streams_inp, dTb_dp_streams_inp,
     +          dTb_dq_streams_inp, dTb_dw_streams_inp)

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

      i = 0
      do jpol = 1,2
         do j = 1,nang/2
            i = i + 1
            Tbo_streams(i) = Tbo_streams_inp(j,jpol)
         end do
      end do

      i = 0
      do jpol = 1, 2
         do jj = 1, nang/2
            do j = 1, nlevels
               i = i + 1
               dTb_dT_streams(i) = dTb_dT_streams_inp(j,jj,jpol)
               dTb_dp_streams(i) = dTb_dp_streams_inp(j,jj,jpol)
               dTb_dq_streams(i) = dTb_dq_streams_inp(j,jj,jpol)
            end do
         end do
      end do

      i = 0
      do jpol = 1, 2
         do jj = 1, nang/2
            do j = 1, nlevels
               i = i + 1
               dTb_dhc_streams(i) = dTb_dw_streams_inp(j,jj,1,jpol)
               dTb_dhr_streams(i) = dTb_dw_streams_inp(j,jj,2,jpol)
               dTb_dhi_streams(i) = dTb_dw_streams_inp(j,jj,3,jpol)
               dTb_dhs_streams(i) = dTb_dw_streams_inp(j,jj,4,jpol)
               dTb_dhg_streams(i) = dTb_dw_streams_inp(j,jj,5,jpol)
            end do
         end do
      end do
      
      ! Create matrices for the return arguments.
      ! classname = 'int16'
      classname = 'double'
      classid = mxClassIDFromClassName(classname)
      ComplexFlag = 0

      ndim = 1
      dims(1) = 2
      plhs(1) = mxCreateNumericArray(ndim, dims, 
     +                               classid, 
     +                               ComplexFlag)

      ndim = 1
      dims(1) = 2
      plhs(2) = mxCreateNumericArray(ndim, dims, 
     +                               classid, 
     +                               ComplexFlag)

      ndim = 2
      dims(1) = nlevels + 1
      dims(2) = 2
      plhs(3) = mxCreateNumericArray(ndim, dims, 
     +                               classid, 
     +                               ComplexFlag)

      do j = 4, 6
        ndim = 2
        dims(1) = nlevels
        dims(2) = 2
        plhs(j) = mxCreateNumericArray(ndim, dims, 
     +                                 classid, 
     +                                 ComplexFlag)
      end do

      ndim = 3
      dims(1) = nlevels
      dims(2) = 5
      dims(3) = 2
      plhs(7) = mxCreateNumericArray(ndim, dims, 
     +                               classid, 
     +                               ComplexFlag)

      ndim = 1
      dims(1) = nang
      plhs(8) = mxCreateNumericArray(ndim, dims, 
     +                               classid, 
     +                               ComplexFlag)

      ndim = 2
      dims(1) = nang/2
      dims(2) = 2
      plhs(9) = mxCreateNumericArray(ndim,dims,
     +                               classid,
     +                               ComplexFlag)

      ndim = 3
      dims(1) = nlevels
      dims(2) = nang/2
      dims(3) = 2
      plhs(10) = mxCreateNumericArray(ndim,dims,
     +                                classid,
     +                                ComplexFlag)

      ndim = 3
      dims(1) = nlevels
      dims(2) = nang/2
      dims(3) = 2
      plhs(11) = mxCreateNumericArray(ndim,dims,
     +                                classid,
     +                                ComplexFlag)

      ndim = 3
      dims(1) = nlevels
      dims(2) = nang/2
      dims(3) = 2
      plhs(12) = mxCreateNumericArray(ndim,dims,
     +                                classid,
     +                                ComplexFlag)

      ndim = 3
      dims(1) = nlevels
      dims(2) = nang/2
      dims(3) = 2
      plhs(13) = mxCreateNumericArray(ndim,dims,
     +                                classid,
     +                                ComplexFlag)

      ndim = 3
      dims(1) = nlevels
      dims(2) = nang/2
      dims(3) = 2
      plhs(14) = mxCreateNumericArray(ndim,dims,
     +                                classid,
     +                                ComplexFlag)

      ndim = 3
      dims(1) = nlevels
      dims(2) = nang/2
      dims(3) = 2
      plhs(15) = mxCreateNumericArray(ndim,dims,
     +                                classid,
     +                                ComplexFlag)

      ndim = 3
      dims(1) = nlevels
      dims(2) = nang/2
      dims(3) = 2
      plhs(16) = mxCreateNumericArray(ndim,dims,
     +                                classid,
     +                                ComplexFlag)

      ndim = 3
      dims(1) = nlevels
      dims(2) = nang/2
      dims(3) = 2
      plhs(17) = mxCreateNumericArray(ndim,dims,
     +                                classid,
     +                                ComplexFlag)

      y1_pr = mxGetPr(plhs(1))
      y2_pr = mxGetPr(plhs(2))
      y3_pr = mxGetPr(plhs(3))
      y4_pr = mxGetPr(plhs(4))
      y5_pr = mxGetPr(plhs(5))
      y6_pr = mxGetPr(plhs(6))
      y7_pr = mxGetPr(plhs(7))
      y8_pr = mxGetPr(plhs(8))
      y9_pr = mxGetPr(plhs(9))
      y10_pr = mxGetPr(plhs(10))
      y11_pr = mxGetPr(plhs(11))
      y12_pr = mxGetPr(plhs(12))
      y13_pr = mxGetPr(plhs(13))
      y14_pr = mxGetPr(plhs(14))
      y15_pr = mxGetPr(plhs(15))
      y16_pr = mxGetPr(plhs(16))
      y17_pr = mxGetPr(plhs(17))

      ! Load the data into MatLab arrays.
      call mxCopyReal8ToPtr(Tb_inp,        y1_pr, 2)
      call mxCopyReal8ToPtr(tau,           y2_pr, 2)
      call mxCopyReal8ToPtr(Tb,            y3_pr, (nlevels+1)*2)
      call mxCopyReal8ToPtr(dTb_dT,        y4_pr, nlevels*2)
      call mxCopyReal8ToPtr(dTb_dp,        y5_pr, nlevels*2)
      call mxCopyReal8ToPtr(dTb_dq,        y6_pr, nlevels*2)
      call mxCopyReal8ToPtr(dTb_dw,        y7_pr, nlevels*5*2)
      call mxCopyReal8ToPtr(stream_angles, y8_pr, nang)
      call mxCopyReal8ToPtr(Tbo_streams,   y9_pr, nang)
      call mxCopyReal8ToPtr(dTb_dT_streams,y10_pr,nlevels*nang)
      call mxCopyReal8ToPtr(dTb_dp_streams,y11_pr,nlevels*nang)
      call mxCopyReal8ToPtr(dTb_dq_streams,y12_pr,nlevels*nang)
      call mxCopyReal8ToPtr(dTb_dhc_streams,y13_pr,nlevels*nang)
      call mxCopyReal8ToPtr(dTb_dhr_streams,y14_pr,nlevels*nang)
      call mxCopyReal8ToPtr(dTb_dhi_streams,y15_pr,nlevels*nang)
      call mxCopyReal8ToPtr(dTb_dhs_streams,y16_pr,nlevels*nang)
      call mxCopyReal8ToPtr(dTb_dhg_streams,y17_pr,nlevels*nang)

      deallocate( atminp, surfinp, aTb, adTb_dT,
     +            adTb_dp, adTb_dq, adTb_dw, Tbo_streams_inp,
     +            dTb_dT_streams_inp, dTb_dp_streams_inp,
     +            dTb_dq_streams_inp, dTb_dw_streams_inp,
     +            stat = alloc_err)
      return
      end
