!
!=======================================================================
      subroutine MultiLine (Numpts,X,Y,Numpts2,X2,Y2, &
       Xmin,Xmax,Ymin,Ymax, &
       DevFlag,LabX,LabY,title,nVec, minP, maxP,  &
       Plotfile, legflag, legend)
!=======================================================================       
! plots one window, multiple lines
!-----------------------------------------------------------------------
!
      implicit none
!
! input variables
      integer numpts      ! Number of points to be plotted for defined grid
      real X(nVec,Numpts) ! X values of points to be plotted
      real Y(nVec,Numpts) ! Y values of points to be plotted
      integer numpts2     ! Number of points for second set to be plotted
      real X2(nVec,Numpts2) ! Second set of X values to be plotted
      real Y2(nVec,Numpts2) ! second set of Y values to be plotted
      real Xmin           ! minimum X value in Domain
      real Ymin           ! Minimum Y value in Domain
      real Xmax           ! Maximum X value in Domain
      real Ymax           ! Maximum Y value in Domain
      Integer DevFlag     ! flag for plotting to screen or PS file
      character*25 LabX   ! Label for x-axis
      character*25 LabY   ! Label for y-axis
      integer nVec        ! total number of vectors in array
      Integer minP        ! minimum plot number for plot subset option
      Integer maxP        ! maximum plot number for plot subset option
      Character*100 Plotfile  ! output file name
      integer LegFlag     ! flag to show legends
      Character*25 legend(nVec)  ! legend for each line
!
! internal variables

      integer IER         ! error code for correct plotting setup
      integer PGBEG       ! internal function to check plotting setup
      Character*100 device ! name of plotting device (screen or file)
      Character*100 Title  ! title for individual plot
      integer k           ! biome index
      integer np          ! total number of plots to make
      integer color       ! color number for plots
      real Xvalue         ! x value for labels
      real Yvalue         ! y value for labels
      real DY             ! incremental Y value between labels
!
! assign target device
      if(DevFlag==1) device='/XWINDOW'
      if(DevFlag==2) device=trim(Plotfile)//'.eps/cps'
      if(DevFlag==3) device=trim(Plotfile)//'.gif/gif'
      if(DevFlag==4) device=trim(Plotfile)//'.png/png'
      print*, 'ok1', DevFlag, device
!
! open plot
      np=MaxP-MinP+1
      IER=PGBEG(0,device,1,1)
      IF (IER.NE.1) STOP
!
! reverse background colors for gif format
      if(DevFlag==3.or.DevFlag==4) then
        call pgscr(0,1.,1.,1.)
        call pgscr(1,0.,0.,0.)
      endif
!
! set line width
      CALL PGSLW(4)
!
! set plot domain
      CALL PGENV(Xmin,Xmax, Ymin,Ymax,0,2)
!
! set labels
      Title=trim(title)
      CALL PGLAB(LabX, LabY, Title)
!
! set color value
      color=1
!
! set line label location
      Xvalue=(Xmax-Xmin)*0.75+Xmin
      Yvalue=Ymax
      DY=(Ymax-Ymin)*.05
!
! scan through lines
      do k=MinP,MaxP
!
! print min and max values
        print*, '  max=',maxval(Y(k,:)),'min=',minval(Y(k,:))
!
! set line color (start with black)
        color=color+1
        CALL PGSCI(color)
!
! print label for line
        Yvalue=yValue-Dy
        if(legflag==1) CALL PGPTXT(Xvalue,Yvalue,0.,0.,legend(k))
!
! plot first line
        CALL PGLINE(numpts,X(k,:),Y(k,:))
!
! Plot second line (solid, fatter)
        CALL PGSLW(15)
        CALL PGLINE(numpts2,X2(k,:),Y2(k,:))
!
! reset line width
        CALL PGSLW(4)
!
      enddo
!
! close plot
      CALL PGEND
!
      return                                                                    
      end
