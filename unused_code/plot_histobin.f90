!
!=======================================================================
      subroutine HistogramBin (nbins,X, minBin, MaxBin, &
       DevFlag,LabX,LabY,title,nVec, minP, maxP, Plotfile)
!=======================================================================       
! generic subroutine to plot a histogram of previously binned data
!-----------------------------------------------------------------------
!
      implicit none
!
! Input variables
      integer nbins       ! number of bins
      real X(nbins,nVec)  ! X values for histogram
      real minBin         ! minimum value for smallest bin
      real maxBin         ! maximum values for largest bin
      Integer DevFlag     ! flag for plotting to screen or PS file
      character*25 LabX   ! Label for x-axis
      character*25 LabY   ! Label for y-axis
      Character*100 title   ! variable used to construct plot title
      integer nVec        ! number of vectors in input array
      Integer minP        ! minimum plot number for plot subset option
      Integer maxP        ! maximum plot number for plot subset option
      Character *11 Plotfile  ! output file name
!
!Internal Variables
      integer IER         ! error code for correct plotting setup
      integer PGBEG       ! internal function to check plotting setup
      Character*30 device ! name of plotting device (screen or file)
      Character*2 nPlot   ! plot number used to construct plot title
      integer np          ! total number of possible plots
      integer k           ! plot index
      integer m           ! bin index
      real Dbin           ! bin width
      real XBins(nbins)   ! left edge values for each bin
      real X2(nVec,nbins)   ! X values for histogram
      real Total          ! total number of points
      real Mean(nVec)       ! total number of points
      Character*7 Tmean   ! mean character version to construct plot title
      real Variance(nVec)   ! total number of points
      Character*7 TVar    ! variance character version to construct plot title
      real Value          ! temp variable
!
! Calculate bin sizes
      Dbin=(MaxBin-MinBin)/real(nBins)
      do m=1,nbins
        Xbins(m)=MinBin+Dbin*(m-1)
      enddo
!
! Normalize the histograms by dividing by the total; calculate mean
      mean=0.
      X2=0.
      do k=MinP,MaxP
        Total=0.
        do m=1,nbins
          Total=Total+X(m,k)
          mean(k)=mean(k)+(Xbins(m)+Dbin/2.)*X(m,k)
        enddo
        if(Total.eq.0.) Total=1.
        X2(k,:)=X(k,:)/Total
        mean(k)=mean(k)/Total
      enddo
!
! Calculate variance
      Variance=0.
      do k=MinP,MaxP
        do m=1,nbins
          Value=(Xbins(m)+Dbin/2.-mean(k))**2.*X2(m,k)
          Variance(k)=Variance(k)+Value
        enddo
      enddo
!
      do k=MinP,MaxP
        print*,'k=', k,' mean=', mean(k), ' Variance=', Variance(k)
      enddo
!
! assign target device
      if(DevFlag.eq.1) device='/xserve'
      if(DevFlag.eq.2) device=trim(Plotfile)//'/cps'
!
! open plot
      np=maxp-minP+1
      if(np.gt.9.and.np.le.12) IER=PGBEG(0,device,4,3)
      if(np.gt.6.and.np.le.9) IER=PGBEG(0,device,3,3)
      if(np.gt.3.and.np.le.6) IER=PGBEG(0,device,3,2)
      if(np.le.3) IER=PGBEG(0,device,np,1)
      IF (IER.NE.1) STOP
!
      do k=MinP,MaxP
!
! set plot domain
        CALL PGENV(minBin,maxBin, 0.,maxval(X2),0,2)
!
! set labels
        write(nPlot,'(i2)') k
        write(Tmean,'(f7.4)') mean(k)
        write(Tvar,'(f7.4)') Variance(k)
        Title=trim(title)//' mean='//Tmean//' var='// &
       Tvar//' biomum '//nPlot
        CALL PGLAB(LabX, LabY, Title)
!
! plot first set of Y values (solid, thinner)
        CALL PGbin(nbins,Xbins,X2(:,k), .false.)
!
      enddo
!
! End plotting
      CALL PGEND
!
      return                                                                    
      end
