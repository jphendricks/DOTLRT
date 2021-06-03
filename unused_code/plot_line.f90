!
!=======================================================================
      subroutine Line (Numpts,X,Y,Numpts2,X2,Y2, &
       Xmin,Xmax,Ymin,Ymax, &
       DevFlag,LabX,LabY,title,nVec, minP, maxP, Plotfile)
!=======================================================================       
! plots 1-2 lines per window
!-----------------------------------------------------------------------
!
      implicit none
!
! input variables
      integer numpts      ! Number of points to be plotted for defined grid
      real X(Numpts)      ! X values of points to be plotted
      real Y(nVec,Numpts) ! Y values of points to be plotted
      integer numpts2     ! Number of points for second set to be plotted
      real X2(Numpts2)    ! Second set of X values to be plotted
      real Y2(nVec,Numpts2) ! second set of Y values to be plotted
      real Xmin           ! minimum X value in Domain
      real Ymin           ! Minimum Y value in Domain
      real Xmax           ! Maximum X value in Domain
      real Ymax           ! Maximum Y value in Domain
      Integer DevFlag     ! flag for plotting to screen or PS file
      character*25 LabX   ! Label for x-axis
      character*25 LabY   ! Label for y-axis
      Character*100 Title  ! title for individual plot
      integer nVec        ! total number of vectors in array
      Integer minP        ! minimum plot number for plot subset option
      Integer maxP        ! maximum plot number for plot subset option
      Character *11 Plotfile  ! output file name
!
! internal variables
      integer IER         ! error code for correct plotting setup
      integer PGBEG       ! internal function to check plotting setup
      Character*30 device ! name of plotting device (screen or file)
      Character*3 nPlot   ! plot number used to construct plot title

      integer k           ! biome index
      integer np          ! total number of plots to make
!
! assign target device
      if(DevFlag.eq.1) device='/xserve'
      if(DevFlag.eq.2) device=trim(Plotfile)//'/.eps'
!
! open plot
      np=MaxP-MinP+1
      if(np.gt.9.and.np.le.12) IER=PGBEG(0,device,4,3)
      if(np.gt.6.and.np.le.9) IER=PGBEG(0,device,3,3)
      if(np.gt.3.and.np.le.6) IER=PGBEG(0,device,3,2)
      if(np.le.3) IER=PGBEG(0,device,np,1)
      IF (IER.NE.1) STOP
!
      do k=MinP,MaxP
!
! set plot domain
        CALL PGSLW(2)
        CALL PGENV(Xmin,Xmax, Ymin,Ymax,0,2)
!
! set labels
        if(MaxP-MinP+1.gt.1) then
          write(nPlot,'(i3.0)') k
          Title=trim(title)//' '//nPlot
        endif
        CALL PGLAB(LabX, LabY, Title)
!
! plot first set of Y values (solid, thinner)
        CALL PGLINE(numpts,X,Y(k,:))
!
! Plot second set of Y Values (solid, fatter)
        CALL PGSLW(4)
        CALL PGLINE(numpts2,X2,Y2(k,:))
!
! reset line width for plotting borders and labels
        CALL PGSLW(2)
!
      enddo
!
! End plotting
      CALL PGEND
!
      return                                                                    
      end
