!
!=======================================================================
      subroutine Histogram (NumPts,nBins,X, Xmin,Xmax, &
       DevFlag,LabX,LabY,Title,nVec, minP, maxP, Plotfile)
!=======================================================================       
! generic subroutine to create a histogram from unbinned data
!-----------------------------------------------------------------------
!
      implicit none
!
! Inputs
      integer NumPts      ! Number of points per vector
      integer nBins       ! Number of bins to be plotted
      real X(nVec,NumPts) ! X values for histogram
      real Xmin           ! minimum X value in Domain
      real Xmax           ! Maximum X value in Domain
      Integer DevFlag     ! flag for plotting to screen or PS file
      character*25 LabX   ! Label for x-axis
      character*25 LabY   ! Label for y-axis
      character*11 Plotfile  ! variable name
      integer nVec        ! number of vectors in input array
      Integer minP        ! minimum plot number for plot subset option
      Integer maxP        ! maximum plot number for plot subset option
!
! internal variables
      Character*30 device ! name of plotting device (screen or file)
      integer IER         ! error code for correct plotting setup
      integer PGBEG       ! internal function to check plotting setup
      Character*3 nPlot   ! plot number used to construct plot title
      Character*50 Title  ! title for individual plot
      integer np          ! total number of possible plots
      integer k           ! plot index
!
! assign target device
      if(DevFlag.eq.1) device='/xserve'
      if(DevFlag.eq.2) device=trim(Plotfile)//'/cps'
!
! open plot
      np=maxP-minP+1
      if(np.gt.6.and.np.le.9) IER=PGBEG(0,device,3,3)
      if(np.gt.3.and.np.le.6) IER=PGBEG(0,device,3,2)
      if(np.le.3) IER=PGBEG(0,device,np,1)
      IF (IER.NE.1) STOP
!
      do k=MinP,MaxP
!
! set plot domain
        CALL PGENV(Xmin,Xmax, 0.,30.,0,2)
!
! set labels
        write(nPlot,'(i3)') k
        Title=trim(title)//' band '//nPlot
        CALL PGLAB(LabX, LabY, Title)
!
! plot first set of Y values (solid, thinner)
        CALL PGHIST(NumPts,X(k,:),minval(X),maxval(X), nBins,1)
!
      enddo
!
! End plotting
      CALL PGEND
!
      return                                                                    
      end
