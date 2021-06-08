
! plot output
!  if (PlotFlag==1) then
!    if (flagDim==2) then
!      if(Yscale==1) then
!        Ymin=minval(Y)
!        Ymax=maxval(Y)
!      endif
!      call output_manipulation
!      call MultiLine (Numpts,X,Y,Numpts,X,Y, Xmin,Xmax,Ymin,Ymax, DevFlag,LabX,LabY,title,nchan, 1,nline, Plotfile, legFlag, Legend)
!    else if (flagDim==3) then
!      call output_manipulation
!      call Contour_color (zval,numpts,numpts, PlotMin(Xvar),PlotMax(Xvar),PlotMin(yvar),PlotMax(yvar), Dx, Dy,  &
!      DevFlag, mincont, maxcont, nc, nchannel, Minchan,Maxchan,  &
!      LabX,LabY,Title,Plotfile, barformat)
!    endif
!  end if ! plot data
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
!
!=======================================================================
      subroutine Contour_color (Z,XNumpts,YNumpts, &
       Xmin, xmax, Ymin, ymax, Dx, Dy, &
       DevFlag, valmin, valmax, numcont, nArray, minP, maxP,  &
       LabX,LabY,Titlebase,Plotfile, barformat)
!=======================================================================       
! generic subroutine to create a black and white contour plot
!-----------------------------------------------------------------------
!
      implicit none
!
      integer IER         ! error code for correct plotting setup
      integer PGBEG       ! internal function to check plotting setup
      Character*100 device ! name of plotting device (screen or file)
      character*25 LabX   ! Label for x-axis
      character*25 LabY   ! Label for y-axis
      Character*100 Plotfile ! variable name
      Character*100 Titlebase  ! title for individual plot
      Character*100 Title  ! title for individual plot
      Character*2 text    ! temporary text variable for creating plot titles
      Integer DevFlag     ! flag for plotting to screen or PS file
      real Valmax   ! max value for color contouring
      real Valmin   ! Min value for color contouring
      real Dval     ! delta value between color contours
      integer numcont          ! Number of Contours in plot
      integer narray      ! total number of arrays in z
      integer np          ! total number of plots
      Integer minP        ! minimum plot number for plot subset option
      Integer maxP        ! maximum plot number for plot subset option

      integer i,j,k,n  !indexes
      real Xmin           ! minimum X value in Domain
      real Ymin           ! Minimum Y value in Domain
      real Xmax           ! Maximum X value in Domain
      real Ymax           ! Maximum Y value in Domain
      real Dx             ! grid Increment in x direction
      real Dy             ! grid Increment in y direction
      integer Xnumpts     ! Number of grid points x-axis
      integer Ynumpts     ! Number of grid points y-axis
      real Z(nArray,Xnumpts,Ynumpts) ! values to be plotted
      real xbox(4) ! array lat pts for pixel box
      real ybox(4) ! array lon pts for pixel box
      real xval
      real yval
! 
! Color bar variables
      integer color ! color index
      integer maxcolor ! max color index
      real dxBar   ! lon width for color bar boxes
      real dyBar   ! lon width for color bar boxes
      integer del   ! delta color bar index for labels
      integer delmax   ! max color bar boxs before skipping labels
      Character*50 CoLab ! Label for color bar values
      character*8 barformat ! format for color bar labels
!
! assign target device
      if(DevFlag==1) device='/xserve'
      if(DevFlag==2) device=trim(Plotfile)//'.eps/cps'
      if(DevFlag==3) device=trim(Plotfile)//'.gif/gif'
      if(DevFlag==4) device=trim(Plotfile)//'.png/png'
!
! open plot 
      np=maxP-minP+1
      if(np.gt.9.and.np.le.12) IER=PGBEG(0,device,4,3)
      if(np.gt.6.and.np.le.9) IER=PGBEG(0,device,3,3)
      if(np.gt.3.and.np.le.6) IER=PGBEG(0,device,3,2)
      if(np.le.3) IER=PGBEG(0,device,np,1)
      IF (IER.NE.1) STOP
!
! reverse background colors for gif format
      if(DevFlag==3.or.DevFlag==4) then
        call pgscr(0,1.,1.,1.)
        call pgscr(1,0.,0.,0.)
      endif
!
! set contour interval
      Dval=(Valmax-Valmin)/real(numcont)
!
! se maximum color index
      maxcolor=19+numcont
!
! set color bar
      call ColorBar(numcont)
!
! create contour plots
      do k=MinP,MaxP
!
! print min and max values
      print*, 'min=',minval(z(k,:,:)),'max=',maxval(z(k,:,:))
!
! set up plot domain
        CALL PGENV(Xmin,Xmax, Ymin,Ymax,0,1)
!
! draw colorbar
!       turn off viewport clipping
        call PGSCLP(0)
!
!       set width of colorbar boxes
        dxbar=(xmax-xmin)/real(numcont)
        dybar=(ymax-ymin)*.02
!
!       set color bar index at minimum
        color=20
!
!       plot box for each color index
        do i=1,numcont
          call pgsci(color)
          xval=xmin+(real(i)-0.5)*dxbar
          yval=Ymax+2*dybar
          xBox(1)=xval-.5*dxbar
          yBox(1)=yval-.5*dybar
          xBox(2)=xval-.5*dxbar
          yBox(2)=yval+.5*dybar
          xBox(3)=xval+.5*dxbar
          yBox(3)=yval+.5*dybar
          xBox(4)=xval+.5*dxbar
          yBox(4)=yval-.5*dybar
          call pgrect(xBox(1), xBox(3),yBox(1), yBox(2))
          color=color+1
        enddo
!
!       set label increment
        delmax=10
        del=1
        if(numcont.gt.20) del=2
        if(numcont.gt.30) del=4
!
!       draw color box labels
        call pgsci(1)
        call pgsch(.75)
        call PGSLW(2)
        do i=1,numcont,del
          xval=xmin+(real(i)-1.)*dxbar
          yval=Ymax+.5*dyBar
          if(barformat(2:2)=='i') then
            n=Valmin+dVal*real(i-1)
            write(colab, barformat) n
          else
            write(colab, barformat) Valmin+dVal*real(i-1)
          endif
          CALL PGPTXT(xval,yval,0.,0.,trim(colab))
        enddo
        call pgsch(1.)
        call PGSCLP(1)
!
! assign axis labels
        call PGNUMB(k, 0,1, text, 30)
        Title=trim(titlebase)
        CALL PGLAB(LabX, LabY,trim(Title))
!
        do i=1,Xnumpts
          xval=xmin+(real(i)-0.5)*dx
          do j=1,Ynumpts
            yval=Ymin+(real(j)-0.5)*dy
! set pixel color
            color=int((z(k,i,j)-Valmin)/Dval)+20
            if (color.gt.maxcolor) color=maxcolor
            if (color.lt.20) color=20
            call pgsci(color)
!
! create box to plot as pixel
            xBox(1)=xval-.5*dx
            yBox(1)=yval-.5*dy
            xBox(2)=xval-.5*dx
            yBox(2)=yval+.5*dy
            xBox(3)=xval+.5*dx
            yBox(3)=yval+.5*dy
            xBox(4)=xval+.5*dx
            yBox(4)=yval-.5*dy
!
! plot pixel
            if(z(k,i,j)/=0.) call pgpoly(4, xbox, ybox)
          Enddo
        Enddo
!
      Enddo
!
! End plotting
      CALL PGEND
      return                                                                    
      end
!
!=======================================================================
      subroutine ColorBar(ncontours)
!=======================================================================       
! creates a color bar for by assigning color indeces.  
! Assumes pgplot has already opened plot and defined viewport
! Given number of contours, divide the color bar into 4 quadrants
! assign color at each quadrant boundary: 0.0, 0.25, 0.5, 0.75, 1.0
! (4 quadrants, 5 boundaries)
! scale Red/Green/Blue intensities between colors at each quadrant boundary
! to change color bar, change colors at quadrant boundaries
! minimum number of contours is 4
! 
! Modifications:
!  Kevin Schaefer created routine (1/10/02)
!  Kevin Schaefer made ncontours an input variable (1/30/04)
!  Kevin Schaefer eliminated link to stats_variabbles module (1/30/04) 
!
      implicit none
!
! internal variables
      integer i        ! color bar index
      integer ilow(4)  ! color bar index lower limit
      integer ihi(4)   ! color bar index upper limit
      integer j        ! segment index
      real Red         ! red color intensity
      real Green       ! green color intensity
      real Blue        ! blue color intensity
      real DRed        ! delta red color intensity
      real DGreen      ! delta green color intensity
      real DBlue       ! delta blue color intensity
      real QRed(5)     ! red color intensity each quadrant boundary
      real QGreen(5)   ! green color intensity each quadrant boundary
      real QBlue(5)    ! blue color intensity each quadrant boundary
      integer color    ! color index
      integer quad(4)  ! contour levels defining the 4 color quadrants
      integer ncontours ! number of contour levels
!
! check for minimum number of contours
      if(ncontours.lt.4) print*, 'Error: need at least 4 contours'
!
! set color bar quadrants
      quad(1)=ncontours/4
      quad(2)=ncontours/2
      quad(3)=ncontours*3/4
      quad(4)=ncontours
!
! set lower bar index limits for each segment
      ilow(1)=1
      ilow(2)=quad(1)+1
      ilow(3)=quad(2)+1
      ilow(4)=quad(3)+1
!
! set upper bar index limits for each segment
      ihi(1)=quad(1)
      ihi(2)=quad(2)
      ihi(3)=quad(3)
      ihi(4)=quad(4)
!
! set colors at each segment boundary (4 segments, 5 boundaries)
!              R    G    B
!   purple     1    0    1
!   Blue       0    0    1
!   Green      0    1    0
!   Yellow     1    1    0
!   Red        1    0    0
!   white      1    1    1
!   black      0    0    0
!   steelgrn   0    .5   .5
!   grey       .5   .5   .5
!   tourquise  0    1   1
!   brown      .5   .5   0
!  start (0.0)
      Qred(1)=1.
      QGreen(1)=0.
      QBlue(1)=1.
!  first quad (0.25)
      Qred(2)=.25
      QGreen(2)=0.25
      QBlue(2)=.75
!  middle (0.5)
      Qred(3)=0.1
      QGreen(3)=.81
      QBlue(3)=.1
!  third quad (0.75)
      Qred(4)=1.
      QGreen(4)=1.
      QBlue(4)=0.
!  end (1.0)
      Qred(5)=1.
      QGreen(5)=0.
      QBlue(5)=0.
!
! set lowest color index
      color=20
!
! scan each quadrant, interpolate betrween intensities at quad boundaries
      do j=1,4
!       set boundary color
        Red=QRed(j)
        Green=QGreen(j)
        Blue=QBlue(j)
!
!       delta color between boundary colors
        DRed=(QRed(j+1)-QRed(j))/real(ihi(j)-ilow(j)+1)
        DGreen=(QGreen(j+1)-QGreen(j))/real(ihi(j)-ilow(j)+1)
        DBlue=(QBlue(j+1)-QBlue(j))/real(ihi(j)-ilow(j)+1)
!
!       increment colors between boundaries
        do i=ilow(j),ihi(j)
!         set color index
          call pgscr(color,Red,Green,Blue)
!
!         increment color
          Red=Red+DRed
          green=green+DGreen
          Blue=Blue+DBlue
!
!         increment color index
          color=color+1
        enddo
      enddo
!
      return                                                                    
      end
!
