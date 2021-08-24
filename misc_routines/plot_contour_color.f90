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
        Title=trim(titlebase)//' BioNum '//text
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
