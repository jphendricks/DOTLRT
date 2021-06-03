!
!=======================================================================
      subroutine Contour (Z,XNumpts,YNumpts,Xmin, Ymin, Dx, Dy, &
       DevFlag, Scale, mincont,maxcont,base, NC, &
       nArray, minP, maxP, &
       LabX,LabY,Title,Plotfile)
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
      Character*100 Title  ! title for individual plot
      Character*2 text    ! temporary text variable for creating plot titles
      Integer DevFlag     ! flag for plotting to screen or PS file
      Integer Scale       ! Contour interval scaling factor for plotting
      real base           ! base forContour levels
      integer NC          ! Number of Contours in plot
      integer narray      ! total number of arrays in z
      integer np          ! total number of plots
      Integer minP        ! minimum plot number for plot subset option
      Integer maxP        ! maximum plot number for plot subset option
      real Cont(NC)       ! Contour levels
      real DCont          ! increment between contour levels
      real MinCont        ! minimum contour level
      real MaxCont        ! minimum contour level
      character*30 Label  ! char version of contour level for labeling contour

      integer k           ! biome index
      integer l           ! contour level index
      real Grid(6)        ! grid definition variables for contour plotting
      real Xmin           ! minimum X value in Domain
      real Ymin           ! Minimum Y value in Domain
      real Xmax           ! Maximum X value in Domain
      real Ymax           ! Maximum Y value in Domain
      real Dx             ! grid Increment in x direction
      real Dy             ! grid Increment in y direction
      integer Xnumpts     ! Number of grid points x-axis
      integer Ynumpts     ! Number of grid points y-axis
      real Z(nArray,Xnumpts, Ynumpts) ! values to be plotted
      integer intbase     ! integer value for converting integer to character
!
! assign target device
      if(DevFlag==1) device='/xserve'
      if(DevFlag==2) device=trim(Plotfile)//'/cps'
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
! create contour plots
      do k=MinP,MaxP
        print*, '   veg type=', k &
      , ' max=',Maxval(Z(k,:,:)), ' min=',Minval(Z(k,:,:))
!
! automatic contour levels
      if(Scale==0) Then
        DCont=(Maxval(Z(k,:,:))-Minval(Z(k,:,:)))/real(NC)
        MinCont=MinVal(Z(k,:,:))
        do l=1,NC
          cont(l)=MinCont+DCont*real(l)
        enddo
      EndIf
!
! manual contout levels, automatic interval
      DCont=(maxcont-mincont)/real(NC)
      if(Scale==1) Then
        do l=1,NC
          cont(l)=MinCont+dcont*real(l-1)
        enddo
      EndIf
!
! manual contours, specify interval
      if(Scale==2) Then
        dCont=base
        do l=1,NC
          cont(l)=mincont+dcont*real(l-1)
        enddo
      EndIf
!
! power scale
      if(Scale==3) Then
        do l=1,NC
          cont(l)=base**real(l)
        enddo
      EndIf
!
! calculate maximum x and y values
      Xmax=Xmin+Dx*(Xnumpts-1)
      Ymax=Ymin+DY*(Ynumpts-1)
!
! set up plot domain
      CALL PGENV(Xmin,Xmax, Ymin,Ymax,0,1)
!
! set up grid calculations for contours (assume uniform grid of points)
      Grid(1)=Xmin-Dx
      Grid(2)=Dx
      Grid(3)=0
      Grid(4)=Ymin-Dy
      Grid(5)=0
      Grid(6)=DY
!
! plot contours
      CALL PGcont(Z(k,:,:),Xnumpts,Ynumpts, &
       1,Xnumpts,1,Ynumpts,cont, NC,Grid)
!
! create plot contour labels
      do l=1,NC
        intbase=aint(cont(l)*1000)
        call PGNUMB(intbase, -3,1, Label, 30)
        CALL PGconl(Z(k,:,:),Xnumpts,Ynumpts, &
       1,Xnumpts,1,Ynumpts,cont(l),Grid,Label, 10, 5)
      enddo
!
! assign axis labels
        call PGNUMB(k, 0,1, text, 30)
        if (np.gt.1) Title=trim(title)//' BioNum '//text
        CALL PGLAB(LabX, LabY,trim(Title))
!
      Enddo
!
! End plotting
      CALL PGEND
      return                                                                    
      end
