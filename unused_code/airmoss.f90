subroutine permaf_soil_param(SOM,sand,clay,FC,Db,phi,b,psi_s,k_s)

    ! A new parametrization for permafrost soils considering both organic matter content
    ! and degree of decomposition

    ! Author:           Richard H. Chen
    ! Email:            chenrh@usc.edu
    ! Date:             02/22/2019

    implicit none
    
    ! Input variables
    real(kind=8) :: SOM     ! (g/g)     mass fraction of soil organic matter (SOM)
    real(kind=8) :: sand    ! (%)       sand fraction
    real(kind=8) :: clay    ! (%)       clay fraction

    ! Output variables
    real(kind=8) :: FC       ! (g/g)     mass fraction of fibrous materials in organic matter
    real(kind=8) :: Db       ! (g/cm3)   bulk density
    real(kind=8) :: phi      ! (m3/m3)   porosity
    real(kind=8) :: b        ! (-)       shape parameter of Campbell soil-water retention curve
    real(kind=8) :: psi_s    ! (cm)      soil matric potential at saturation
    real(kind=8) :: k_s      ! (m/s)     hydraulic conductivity at saturation
    real(kind=8) :: theta_wp ! (m3/m3)  volumetric water content at wilting point

    ! Local variables
    real(kind=8) :: Db_m    ! (g/cm3)   bulk density of pure mineral soil
    real(kind=8) :: phi_m   ! (m3/m3)   porosity of pure mineral soil
    real(kind=8) :: r       ! (-)       decreasing rate of bulk density from pure mineral to pure organic
    real(kind=8) :: b_m     ! (-)       b of pure mineral soil
    real(kind=8) :: psi_sm  ! (cm)      psi_s of pure mineral soil
    real(kind=8) :: k_sm    ! (cm)      k_s of pure mineral soil
    real(kind=8) :: v_m     ! (m3/m3)   volumetric fraction of mineral soil
    real(kind=8) :: v_s     ! (m3/m3)   volumetric fraction of sapric soil
    real(kind=8) :: v_f     ! (m3/m3)   volumetric fraction of fibric soil
    real(kind=8) :: fm_m    ! (m3/m3)   sub-phase mass fraction of mineral soil
    real(kind=8) :: fm_s    ! (m3/m3)   sub-phase mass fraction of sapric soil
    real(kind=8) :: fm_f    ! (m3/m3)   sub-phase mass fraction of fibric soil
    real(kind=8) :: fv_m    ! (m3/m3)   sub-phase volumetric fraction of mineral soil
    real(kind=8) :: fv_s    ! (m3/m3)   sub-phase volumetric fraction of sapric soil
    real(kind=8) :: fv_f    ! (m3/m3)   sub-phase volumetric fraction of fibric soil

    real(kind=8), parameter :: rho_m  = 2.65d0           ! (g/cm3) specific density of mineral particles
    real(kind=8), parameter :: rho_s  = 1.80d0           ! (g/cm3) specific density of fibric particles
    real(kind=8), parameter :: rho_f  = 0.60d0           ! (g/cm3) specific density of sapric particles
    real(kind=8), parameter :: Db_o   = 0.5d0            ! (g/cm3) bulk density of pure organic soil
    real(kind=8), parameter :: SOM_s  = 0.45d0           ! (g/g)   SOM fraction of sapric peat in Letts et al. 2000
    real(kind=8), parameter :: SOM_f  = 1.d0             ! (g/g)   SOM fraction of fibric peat in Letts et al. 2000
    real(kind=8), parameter :: b_s    = 26.754d0         ! (-)     b of pure sapric soil
    real(kind=8), parameter :: b_f    = 0.035d0          ! (-)     b of pure fibric soil
    real(kind=8), parameter :: psi_ss = -1.01d0          ! (cm)    psi_s of pure sapric soil
    real(kind=8), parameter :: psi_sf = -1.03d0          ! (cm)    psi_s of pure fibric soil
    real(kind=8), parameter :: k_ss   = 7.9644d-12       ! (-)     k_s of pure sapric soil
    real(kind=8), parameter :: k_sf   = 1.9195d-3        ! (-)     k_s of pure fibric soil
    real(kind=8), parameter :: psi_wp = -15000.d0        ! (cm)    soil matric potential at wilting point
  
    ! Pedotransfer functions of mineral soils in Cosby et al. 1984
    b_m = 2.91d0 + 0.159d0*clay
    psi_sm = -10.d0**(1.88d0-0.0131d0*sand)
    k_sm = 2.54d-2/3600.d0*10.d0**(-0.884d0+0.0153d0*sand)
    phi_m = 0.489d0 - 0.00126d0*sand
    Db_m = (1.d0 - phi_m)*rho_m

    r = -log(Db_o / Db_m)

    Db = Db_m*exp(-r*SOM)
    FC = 0.9887d0*exp(-5.512d0*Db)      ! Exponential fit from Boelter 1969

    ! Volumetric fractions of solid materials
    v_m = Db/rho_m*(1.d0-SOM)
    v_s = Db/rho_s*SOM*(1.d0-FC)
    v_f = Db/rho_f*SOM*FC

    ! Porosity
    phi = 1.d0 - v_m - v_s - v_f

    ! Sub-phase mass fractions of solid materials
    fm_m = 1.d0-SOM
    fm_s = SOM*(1.d0-FC)
    fm_f = SOM*FC

    ! Sub-phase volumetric fractions of solid materials
    fv_m = fm_m*Db/(1.d0-phi)/rho_m
    fv_s = fm_s*Db/(1.d0-phi)/rho_s
    fv_f = fm_f*Db/(1.d0-phi)/rho_f

    b = fv_m*b_m + fv_s*b_s + fv_f*b_f

    if(SOM >= SOM_s) then
        psi_s = psi_ss + (SOM-SOM_s)/(SOM_f-SOM_s)*(psi_sf-psi_ss)
    else
        psi_s = (psi_sm-psi_ss)/SOM_s**2*(SOM-SOM_s)**2 + psi_ss
    endif

    k_s = k_sm**fv_m*k_ss**fv_s*k_sf**fv_f

    theta_wp = phi*(psi_wp/psi_s)**(-1.d0/b)
    
end subroutine permaf_soil_param
!
!=======================================================================
subroutine dielectric_constant_exp(depth1, depth2, sat_top, watertab, sand, clay,mass_org, di_org_exp, test, test10)
!=======================================================================
! Calculate volumetric water content, porosity, and dielectric constant
! as a function of depth using the ReSALT exponential organic matter model
!
! Modifications:
!  Kevin Schaefer made subroutine from dielectric_constant (2/16/19)
!  Kevin Schaefer deleted mineral soil and organic layer model (2/16/19)
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
! Input variables
real(kind=8) :: depth1       ! (m) depth from surface to top of layer
real(kind=8) :: depth2       ! (m) depth from surface to bottom of layer
real(kind=8) :: sat_top      ! (-) saturation fraction at top of soil column (zero depth)
real(kind=8) :: watertab     ! (m) depth to water table
real         :: sand         ! (%) sand fraction of soil texture
real         :: clay         ! (%) clay fraction of soil texture
real         :: mass_org     ! (kg/m2) mass of organic matter in top 1 m of soil
!
! Output variables
real(kind=8) :: di_org_exp   ! (-) dielectric constant mix of organic and mineral soil
real :: test, test10(10)     ! test variables
!
! Local variables
real(kind=8) :: sat_slop     ! (-) saturation fraction slope
real(kind=8) :: FC           ! (g/g) mass fraction of fibrous materials in organic matter
real(kind=8) :: Db           ! (g/cm3) bulk density
real(kind=8) :: b_shape      ! (-) shape parameter of Campbell soil-water retention curve
real(kind=8) :: psi_s        ! (cm) soil matric potential at saturation
real(kind=8) :: k_s          ! (m/s)  hydraulic conductivity at saturation
real(kind=8) :: depth        ! (m) depth from surfaceairmoss
real(kind=8) :: satfrac      ! (-) saturation fraction
real(kind=8) :: loc_sand     ! (-) sand fraction of soil texture
real(kind=8) :: loc_clay     ! (-) sand fraction of soil texture
real(kind=8) :: loc_mass_org ! (kg/m2) mass of organic matter in top 1 m of soil
real(kind=8) :: kRoot        ! (1/m) Exp const for decrease root-zone organic matter w/ depth
real(kind=8) :: root_depth   ! (m) maximum rooting depth
real(kind=8) :: rho_om_max   ! (kg/m3) maximum organic matter density
real(kind=8) :: vwc_org_exp  ! (-) vwc mix of organic and mineral soil
real(kind=8) :: vwc_perc     ! (%) vwc in percent
real(kind=8) :: poros        ! (-) porosity for organic mineral soil mix
real(kind=8) :: poros_om     ! (-) soil porosity for pure organic soil
real(kind=8) :: poros_min    ! (-) mineral soil porosity
real(kind=8) :: orgfrac      ! (-) organic soil fraction
real(kind=8) :: a, b, c, d   ! (-) dielectric curve fit coefficients
integer iz           ! z index
integer numz         ! number z slices
real(kind=8) :: dz   ! (m) delta z
real(kind=8) :: sum  ! (-) total
!
! calculate saturation fraction slope
    sat_slop=(1.d0-sat_top)/watertab
!
! set vertical organic matter distribution parameters
    loc_sand=sand/100.
    loc_clay=clay/100.
    loc_mass_org=mass_org
    kroot=5.5
    rho_om_max=140.
    poros_om=.9
    root_depth=0.7
!
! set dielectric constant curve fit coefficients
    a=-0.0001
    b=0.0208
    c=-0.2
    d=4.
!
! calculate mineral soil porosity
    poros_min=0.489-0.00126*sand
!
! integrate over depths
    sum=0
    numz=100
    dz=(depth2-depth1)/real(numz)
    depth=depth1
    do iz=1,numz
!
! organic soil fraction
      orgfrac=kRoot*mass_org/(1.-exp(-kRoot*root_depth))*exp(-kRoot*(depth))/rho_om_max
      if(orgfrac>1.) orgfrac=1.
      if(orgfrac<0.) orgfrac=0.
      if(depth>root_depth) orgfrac=0.
!
! soil porosity
      poros=(1.-orgfrac)*poros_min+orgfrac*poros_om
      if(poros>1.) poros=1.
      call permaf_soil_param(orgfrac,sand,clay,FC,Db,poros,b_shape,psi_s,k_s)
!
! Saturation fraction
      satfrac=sat_top+depth*sat_slop
      satfrac=max(satfrac, 0.d0)
      satfrac=min(satfrac, 1.d0)
!
! volume of water and pore space
      vwc_org_exp=poros*satfrac
      vwc_perc=vwc_org_exp*100.
      sum=sum+(a*vwc_perc**3.+b*vwc_perc**2.+c*vwc_perc+d)
      depth=depth+dz
    enddo
    di_org_exp=sum/real(numz)
!
! assign test variables
    test=di_org_exp
    test10(1)=di_org_exp
!
RETURN
END
!
!--------------------------------------------------------------------------------
subroutine airmoss_fwd_model(nlay,freq,thetainc,ls,sh,dlay,diel,sigma_vv,sigma_hh, test, test10)
!--------------------------------------------------------------------------------
! calculates the bistatic scattering coefficients of a 
! two layer dielectric structure with two rough boundaries using the 
! method in A .tabatabaeenejad and M. moghaddam, "bistatic 
! scattering from three-dimensional layered rough surfaces," ieee 
! transactions on geoscience and remote sensing, vol. 44, no. 8, pp. 
! 2102-2114, aug. 2006. 
!
! code by Alireza Tabatabaeenejad
! radiation laboratory
! department of electrical engineering and computer science
! university of michigan, ann arbor
!
! note:
! - top, middle, and bottom layers are numbered 0, 1 & 2, respectively
! - subscripts 'h' and 'v' refer to h-pol and v-pol, respectively
!
! history
!  1/26/19: Schaefer cleaned up code
!  2/15/19  Schaefer moved angle calculations inside routine
!  2/15/19  Schaefer changed routine name to airmoss_fwd_model
!--------------------------------------------------------------------------------
implicit none
!
! define input variables
integer,         intent(in)  :: nlay       ! (-) number of soil layers
real(kind=8),    intent(in)  :: freq       ! (hz) frequency
real(kind=8),    intent(in)  :: thetainc   ! (deg) scattering angle
real(kind=8),    intent(in)  :: ls(nlay)   ! (m) correlation lengths of layer boundaries
real(kind=8),    intent(in)  :: sh(nlay)   ! (m) standard deviation of layer boundaries
real(kind=8),    intent(in)  :: dlay(nlay-1) ! (m) average z coordinate of layer boundaries
complex(kind=8), intent(in)  :: diel(nlay) ! (-) relative dielectric constant of layers.
real test, test10(10)     ! test variables
!
! define output variables
real(kind=8),    intent(out) :: sigma_vv ! (-) Scattering coefficients of the layered structure for VV polarization
real(kind=8),    intent(out) :: sigma_hh ! (-) Scattering coefficients of the layered structure for HH polarization
!
! Define local variables
real(kind=8)    :: thi           ! (deg) Incidence angle. >90 and is measured from positive z axis in spherical coordinate system.
real(kind=8)    :: phi           ! (deg) Incidence azimuth angle
real(kind=8)    :: ths           ! (deg) Scattered angle
real(kind=8)    :: phs           ! (deg) Scattered azimuth angle
real(kind=8)    :: wvl           ! (m) Free-space wavelength
real(kind=8)    :: k0            ! Incidence region wavenumber
real(kind=8)    :: z0            ! Incidence region intrinsic impedance
complex(kind=8) :: k(nlay)       ! (?) Wavenumber per layer
complex(kind=8) :: z(nlay)       ! Intrinsic impedance per layer
real(kind=8)    :: kxi           ! x-componenet of incidence wavenumber vector
real(kind=8)    :: kxs           ! x-componenet of scattering wavenumber vector
real(kind=8)    :: kyi           ! y-componenet of incidence wavenumber vector
real(kind=8)    :: kys           ! y-componenet of scattering wavenumber vector
real(kind=8)    :: aem           ! V-pol components of incident electric field intensity
real(kind=8)    :: ahm           ! H-pol components of incident electric field intensity
complex(kind=8) :: x0(4*nlay)    ! Vector containing zeroth-order unknown coefficients
complex(kind=8) :: alpha_h(nlay)
complex(kind=8) :: alpha_v(nlay)
real(kind=8)    :: spec_den(nlay) ! Spectral density of the rough boundary per layer
complex(kind=8) :: a0(4*nlay,4*nlay)
complex(kind=8) :: y0(4*nlay,4*nlay)
complex(kind=8) :: an(4*nlay,4*nlay)
complex(kind=8) :: yn(4*nlay,4*nlay)
complex(kind=8) :: b1(4*nlay)
real(kind=8)    :: b0(4*nlay)
real(kind=8)    :: bisct(2,2)    ! (-) bistatic scattering coeficients
                                 !     bisct(1,1): sigma_vv
                                 !     bisct(1,2): sigma_vh
                                 !     bisct(2,1): sigma_hv
                                 !     bisct(2,2): sigma_hh
real(kind=8)    :: sigma_hv      ! (-) Scattering coefficients of the layered structure for HV polarization
real(kind=8)    :: sigma_vh      ! (-) Scattering coefficients of the layered structure for VH polarization
complex(kind=8) :: sums
integer :: ilay ! layer index
integer :: ipol ! polarization index
real(kind=8), parameter :: pi = 3.141592653589793d0 ! value of pi
real(kind=8), parameter :: pi180 = pi/180.d0 ! value of pi
real(kind=8), parameter :: c0 = 299792458.d0        ! (m/s) free space speed of light
real(kind=8), parameter :: mu0 = 4.d0*pi*1.d-7      ! (?) Free-space permeability
real(kind=8), parameter :: eps0 = 1.d0/(c0**2*mu0)  ! (?) Free-space permittivity
!
! backscatter angles in radians
    thi = (180.d0 - thetainc)*pi180
    ths = thetainc*pi180
    phi = 0.d0
    phs = pi
!
! free space wavelength and wave number
    wvl=c0/freq
    k0=2.d0*pi/wvl
!
! free space impedence
    z0=sqrt(mu0/eps0)
!
! wave number and impedence per layer
    do ilay=1,nlay
        k(ilay)=k0*sqrt(diel(ilay))
        z(ilay)=z0/sqrt(diel(ilay))
    enddo
!
! x and y components of waves
    kxi=k0*sin(thi)*cos(phi)
    kyi=k0*sin(thi)*sin(phi)
    kxs=k0*sin(ths)*cos(phs)
    kys=k0*sin(ths)*sin(phs)
!
! spectral densities of rough boundaries.
    do ilay=1,nlay
!      spec_den(ilay)=exp((-(kxs-kxi)**2-(kys-kyi)**2)*ls(ilay)**2/4.d0)*sh(ilay)**2*ls(ilay)**2/4.d0/pi       ! gaussian
      spec_den(ilay)=(1.d0+((kxs-kxi)**2+(kys-kyi)**2)*ls(ilay)**2)**(-1.5d0)*sh(ilay)**2*ls(ilay)**2/2.d0/pi ! exponential
    enddo
!
! calculate coefficient matrices incident wave
    !call calculate_coefficients(nlay,dlay,k0,k,kxi,kyi,z0,z,a0)
    call EQU_COEFF(nlay,dlay,k0,k,kxi,kyi,z0,z,a0)

    call inverse_matrix(a0,4*nlay,4*nlay,y0)
!
! calculate coefficient matrices scattered wave
    !call calculate_coefficients(nlay,dlay,k0,k,kxs,kys,z0,z,an)
    call EQU_COEFF(nlay,dlay,k0,k,kxs,kys,z0,z,an)
    call inverse_matrix(an,4*nlay,4*nlay,yn)
!
! calculate sigma_vv and sigma_hv
    aem=0.d0
    ahm=1.d0

    call zero_order_unknowns(nlay,k0,kxi,kyi,z0,aem,ahm,b0)
    
    sums=(0.d0,0.d0)
    do ilay=1,4*nlay
        do ipol=1,4*nlay
            sums=sums+y0(ilay,ipol)*b0(ipol)
        enddo
        x0(ilay)=sums
        sums=(0.d0,0.d0)
    enddo
    
    call e_rgt1(nlay,aem,ahm,dlay,k0,k,kxi,kxs,kyi,kys,z0,z,x0,b1)

    do ilay=1,nlay
        alpha_h(ilay)=(0.d0,0.d0)
        alpha_v(ilay)=(0.d0,0.d0)
        do ipol=1,4
            alpha_h(ilay)=alpha_h(ilay)+yn(1,4*(ilay-1)+ipol)*b1(4*(ilay-1)+ipol)
            alpha_v(ilay)=alpha_v(ilay)+yn(2,4*(ilay-1)+ipol)*b1(4*(ilay-1)+ipol)
        enddo
    enddo

    bisct(1,1)=0.d0
    bisct(2,1)=0.d0
    do ilay=1,nlay
        bisct(1,1)=bisct(1,1)+abs(alpha_v(ilay))**2*spec_den(ilay)
        bisct(2,1)=bisct(2,1)+abs(alpha_h(ilay))**2*spec_den(ilay)
    enddo
    bisct(1,1)=4*pi*(k0**2)*(cos(ths)**2)*bisct(1,1)
    bisct(2,1)=4*pi*(k0**2)*(cos(ths)**2)*bisct(2,1)
    sigma_vv=bisct(1,1)
    sigma_hv=bisct(2,1)

! calculate sigma_hh and sigma_vh
    aem=1.d0
    ahm=0.d0

    call zero_order_unknowns(nlay,k0,kxi,kyi,z0,aem,ahm,b0)
    
    sums=(0.d0,0.d0)
    do ilay=1,4*nlay
        do ipol=1,4*nlay
            sums=sums+y0(ilay,ipol)*b0(ipol)
        enddo
        x0(ilay)=sums
        sums=(0.d0,0.d0)
    enddo

    call e_rgt1(nlay,aem,ahm,dlay,k0,k,kxi,kxs,kyi,kys,z0,z,x0,b1)

    do ilay=1,nlay
        alpha_h(ilay)=(0.d0,0.d0)
        alpha_v(ilay)=(0.d0,0.d0)
        do ipol=1,4
            alpha_h(ilay)=alpha_h(ilay)+yn(1,4*(ilay-1)+ipol)*b1(4*(ilay-1)+ipol)
            alpha_v(ilay)=alpha_v(ilay)+yn(2,4*(ilay-1)+ipol)*b1(4*(ilay-1)+ipol)
        enddo
    enddo

    bisct(1,2)=0.d0
    bisct(2,2)=0.d0
    do ilay=1,nlay
        bisct(1,2)=bisct(1,2)+abs(alpha_v(ilay))**2*spec_den(ilay)
        bisct(2,2)=bisct(2,2)+abs(alpha_h(ilay))**2*spec_den(ilay)
    enddo
    bisct(1,2)=4*pi*(k0**2)*(cos(ths)**2)*bisct(1,2)
    bisct(2,2)=4*pi*(k0**2)*(cos(ths)**2)*bisct(2,2)
    sigma_vh=bisct(1,2)
    sigma_hh=bisct(2,2)

    test=sigma_vh
    test10(1)=sigma_vh
    test10(2)=sigma_hh
return  
end subroutine airmoss_fwd_model
!
!--------------------------------------------------------------------------------
subroutine zero_order_unknowns(nlay,k0,kx,ky,z0,aem,ahm,b0)
!--------------------------------------------------------------------------------
! this subroutine calculates the zeroth-order unknowns based on the 
! formulation presented in the paper, and saves them in vector x0.  
!
implicit none
!
! define input variables
integer,      intent(in)  :: nlay  ! (-) number of soil layers
real(kind=8), intent(in)  :: k0  ! 
real(kind=8), intent(in)  :: z0  ! 
real(kind=8), intent(in)  :: kx  ! 
real(kind=8), intent(in)  :: ky  ! 
real(kind=8), intent(in)  :: aem ! V-pol components of incident electric field intensity
real(kind=8), intent(in)  :: ahm ! H-pol components of incident electric field intensity
!
! define output variables
real(kind=8),    intent(out) :: b0(4*nlay)
!
! define local variables
real(kind=8) kro
real(kind=8) kz0
integer ilay       ! index
!
    kro=sqrt(kx**2+ky**2)
    kz0=sqrt(k0**2-kro**2)

    b0(1)=-(kx/kro*aem-ky*kz0/k0/kro*ahm)
    b0(2)=-(ky/kro*aem+kx*kz0/k0/kro*ahm)
    b0(3)=-(ky*kz0/k0/kro*aem+kx/kro*ahm)/z0
    b0(4)=+(kx*kz0/k0/kro*aem-ky/kro*ahm)/z0
    do ilay=5,4*nlay
      b0(ilay)=0.d0
    end do

return
end subroutine zero_order_unknowns
!
!--------------------------------------------------------------------------------
subroutine inverse_matrix(a,n,np,y)
!--------------------------------------------------------------------------------
! inverts matrix a and returns the result to y using subroutine gausj from
! numerical recipes in fortran 77.
!
implicit none
!
! define input variables
integer,         intent(in)  :: n
integer,         intent(in)  :: np
complex(kind=8), intent(in)  :: a(np,np)
!
! define output variables
complex(kind=8), intent(out) :: y(np,np)
!
! define local variables
integer, parameter :: m = 1
integer, parameter :: mp = 1
complex(kind=8) :: b(np,mp)
integer ilay ! row index
integer jj ! column index
!
! Transfer to local matrix
    do ilay=1,np
      do jj=1,np
        y(ilay,jj)=a(ilay,jj)
      end do
    end do
!
! invert matrix
    call gaussj(y,n,np,b,m,mp)
!
return
end subroutine inverse_matrix
!
!--------------------------------------------------------------------------------
subroutine gaussj(a,n,np,b,m,mp)
!--------------------------------------------------------------------------------
! this subroutine inverts matrix a and returns the result to y. 
! source: numerical recipes in fortran 77.
implicit none

integer,         intent(in)  :: m,mp,n,np
complex(kind=8)              :: a(np,np)
complex(kind=8), intent(out) :: b(np,mp)
integer, parameter :: nmax = 10000
integer         :: i, icol, irow, j, k, l, ll, indxc(nmax), indxr(nmax), ipiv(nmax)
real(kind=8)    :: big
complex(kind=8) :: dum,pivinv

        do j=1,n
            ipiv(j)=0
        enddo
        do 22 i=1,n
            big=0.
            do 13 j=1,n
                if (ipiv(j).ne.1) then
                    do 12 k=1,n
                    if (ipiv(k).eq.0) then
                    if (abs(a(j,k)).ge.big)then
                    big=abs(a(j,k))
                    irow=j
                    icol=k
                    endif
                    else if (ipiv(k).gt.1) then
                    print *, 'gaussj: singular matrix'
                    stop
                    endif
    12              enddo
                endif
    13      enddo
            if (icol.eq.0) then
                print *,'gaussj: icol not set!'
                stop
            endif
            ipiv(icol)=ipiv(icol)+1
            if (irow.ne.icol) then
                do 14 l=1,n
                    dum=a(irow,l)
                    a(irow,l)=a(icol,l)
                    a(icol,l)=dum
    14          enddo
                do 15 l=1,m
                    dum=b(irow,l)
                    b(irow,l)=b(icol,l)
                    b(icol,l)=dum
    15          enddo
            endif
            indxr(i)=irow
            indxc(i)=icol
            if (a(icol,icol).eq.0.) stop 'gaussj: singular matrix'
            pivinv=1./a(icol,icol)
            a(icol,icol)=1.
                do 16 l=1,n
                    a(icol,l)=a(icol,l)*pivinv
    16          enddo
                do 17 l=1,m
                    b(icol,l)=b(icol,l)*pivinv
    17          enddo
                do 21 ll=1,n
                    if (ll.ne.icol) then
                        dum=a(ll,icol)
                        a(ll,icol)=0.
                        do 18 l=1,n
                        a(ll,l)=a(ll,l)-a(icol,l)*dum
    18                  enddo
                        do 19 l=1,m
                        b(ll,l)=b(ll,l)-b(icol,l)*dum
    19                  enddo
                    endif
    21          enddo
    22  enddo
        do 24 l=n,1,-1
            if (indxr(l).ne.indxc(l)) then
                do 23 k=1,n
                    dum=a(k,indxr(l))
                    a(k,indxr(l))=a(k,indxc(l))
                    a(k,indxc(l))=dum
    23          enddo
            endif
    24  enddo

return
end subroutine gaussj
!
!--------------------------------------------------------------------------------
subroutine e_rgt1_new(nlay,aem,ahm,dlay,k0,k,kxi,kxs,kyi,kys,z0,z,x0,b)
!--------------------------------------------------------------------------------
! this subroutine calculates the elements of the excitation vector in 
! a1x1=b1 based on the formulation in the paper, and saves them in vector b.
!
implicit none
!
! define input variables
integer,         intent(in)  :: nlay
real(kind=8),    intent(in)  :: aem
real(kind=8),    intent(in)  :: ahm
real(kind=8),    intent(in)  :: dlay(nlay-1)
real(kind=8),    intent(in)  :: k0
real(kind=8),    intent(in)  :: kxi
real(kind=8),    intent(in)  :: kxs
real(kind=8),    intent(in)  :: kyi
real(kind=8),    intent(in)  :: kys
real(kind=8),    intent(in)  :: z0
complex(kind=8), intent(in)  :: x0(4*nlay) ! Vector containing zeroth-order unknown coefficients
complex(kind=8), intent(in)  :: k(nlay)    ! (?) Wavenumber per layer
complex(kind=8), intent(in)  :: z(nlay)    ! Intrinsic impedance per layer
!
! define output variable
complex(kind=8), intent(out) :: b(4*nlay)
!
! define local variables
real(kind=8)    :: kroi ! Transverse component of incidence wavenumber vector
real(kind=8)    :: kzi0
complex(kind=8) :: i
complex(kind=8) :: kzi(nlay)
complex(kind=8) :: r1
complex(kind=8) :: r2
complex(kind=8) :: r3
complex(kind=8) :: r4
complex(kind=8) :: r5
complex(kind=8) :: r6
complex(kind=8) :: r7
complex(kind=8) :: r8
integer         :: ilay
complex(kind=8) :: exp_01
complex(kind=8) :: exp_02
complex(kind=8) :: exp_03
complex(kind=8) :: exp_04
complex(kind=8) :: exp_05
complex(kind=8) :: exp_06
complex(kind=8) :: exp_07


        i=(0.d0,1.d0)

        kroi=sqrt(kxi**2+kyi**2)
        kzi0=sqrt(k0**2-kroi**2)

        do ilay=1,nlay
            kzi(ilay)=sqrt(k(ilay)**2-kroi**2)
        end do

        r1=kzi0*kxi/kroi*x0(1)  
        r2=(kzi0*kzi0*kyi/k0/kroi-(kys-kyi)*kroi/k0)*x0(2)
        r3=-(kzi(1)*kxi/kroi)*x0(3)
        r4=((kys-kyi)*kroi/k(1)-kzi(1)**2*kyi/k(1)/kroi)*x0(4)
        r5=(kzi(1)*kxi/kroi)*x0(5)
        r6=((kys-kyi)*kroi/k(1)-kyi*kzi(1)**2/k(1)/kroi)*x0(6)
        r7=-(kzi0*kxi/kroi)*aem
        r8=(kyi*kzi0*kzi0/k0/kroi-(kys-kyi)*kroi/k0)*ahm
        b(1)=-i*(r1+r2+r3+r4+r5+r6+r7+r8)

        r1=(kzi0*kyi/kroi)*x0(1)
        r2=(-kzi0*kzi0*kxi/k0/kroi+(kxs-kxi)*kroi/k0)*x0(2)
        r3=-(kzi(1)*kyi/kroi)*x0(3)
        r4=-((kxs-kxi)*kroi/k(1)-kzi(1)*kzi(1)*kxi/k(1)/kroi)*x0(4)
        r5=(kzi(1)*kyi/kroi)*x0(5)
        r6=-((kxs-kxi)*kroi/k(1)-kxi*kzi(1)*kzi(1)/k(1)/kroi)*x0(6)
        r7=-(kzi0*kyi/kroi)*aem
        r8=(-kxi*kzi0*kzi0/k0/kroi+(kxs-kxi)*kroi/k0)*ahm
        b(2)=-i*(r1+r2+r3+r4+r5+r6+r7+r8)

        r1=(-(kys-kyi)*kroi/k0+kzi0*kzi0*kyi/k0/kroi)*x0(1)/z0
        r2=-(kxi*kzi0/kroi)*x0(2)/z0
        r3=((kys-kyi)*kroi/k(1)-kyi*kzi(1)*kzi(1)/k(1)/kroi)*x0(3)/z(1)
        r4=(kxi*kzi(1)/kroi)*x0(4)/z(1)
        r5=((kys-kyi)*kroi/k(1)-kyi*kzi(1)*kzi(1)/k(1)/kroi)*x0(5)/z(1)
        r6=-(kxi*kzi(1)/kroi)*x0(6)/z(1)
        r7=(kyi*kzi0*kzi0/k0/kroi-(kys-kyi)*kroi/k0)*aem/z0
        r8=(kxi*kzi0/kroi)*ahm/z0
        b(3)=i*(r1+r2+r3+r4+r5+r6+r7+r8)

        r1=((kxs-kxi)*kroi/k0-kzi0*kzi0*kxi/k0/kroi)*x0(1)/z0
        r2=-(kyi*kzi0/kroi)*x0(2)/z0
        r3=(-(kxs-kxi)*kroi/k(1)+kxi*kzi(1)*kzi(1)/k(1)/kroi)*x0(3)/z(1)
        r4=(kyi*kzi(1)/kroi)*x0(4)/z(1)
        r5=(-(kxs-kxi)*kroi/k(1)+kxi*kzi(1)*kzi(1)/k(1)/kroi)*x0(5)/z(1)
        r6=-(kyi*kzi(1)/kroi)*x0(6)/z(1)
        r7=(-kxi*kzi0*kzi0/k0/kroi+(kxs-kxi)*kroi/k0)*aem/z0
        r8=(kyi*kzi0/kroi)*ahm/z0
        b(4)=i*(r1+r2+r3+r4+r5+r6+r7+r8)

        do ilay=2,nlay-1
            exp_01=exp(-i*kzi(ilay-1)*dlay(ilay-1))*x0(4*ilay-5)
            exp_02=exp(-i*kzi(ilay-1)*dlay(ilay-1))*x0(4*ilay-4)
            exp_03=exp(i*kzi(ilay-1)*dlay(ilay-1))*x0(4*ilay-3)
            exp_04=exp(i*kzi(ilay-1)*dlay(ilay-1))*x0(4*ilay-2)
            exp_05=exp(-i*kzi(ilay)*dlay(ilay-1))*x0(4*ilay)
            exp_06=exp(i*kzi(ilay)*dlay(ilay-1))*x0(4*ilay+1)
            exp_07=exp(i*kzi(ilay)*dlay(ilay-1))*x0(4*ilay+2)
            
            r1=-(kxi*kzi(ilay-1)/kroi)*exp_01
            r2=((kys-kyi)*kroi/k(ilay-1)-kyi*kzi(ilay-1)*kzi(ilay-1)/k(ilay-1)/kroi)*exp_02
            r3=(kxi*kzi(ilay-1)/kroi)*exp_03
            r4=((kys-kyi)*kroi/k(ilay-1)-kyi*kzi(ilay-1)*kzi(ilay-1)/k(ilay-1)/kroi)*exp_04
            r5=(kxi*kzi(ilay)/kroi)*exp(-i*kzi(ilay)*dlay(ilay-1))*x0(4*ilay-1)
            r6=-((kys-kyi)*kroi/k(ilay)-kyi*kzi(ilay)*kzi(ilay)/k(ilay)/kroi)*exp_05
            r7=-(kxi*kzi(ilay)/kroi)*exp_06
            r8=(kyi*kzi(ilay)*kzi(ilay)/k(ilay)/kroi-(kys-kyi)*kroi/k(ilay))*exp_07
            b(4*ilay-3)=i*(r1+r2+r3+r4+r5+r6+r7+r8)

            r1=-(kyi*kzi(ilay-1)/kroi)*exp_01
            r2=(-(kxs-kxi)*kroi/k(ilay-1)+kxi*kzi(ilay-1)*kzi(ilay-1)/k(ilay-1)/kroi)*exp_02
            r3=(kyi*kzi(ilay-1)/kroi)*exp_03
            r4=(-(kxs-kxi)*kroi/k(ilay-1)+kxi*kzi(ilay-1)*kzi(ilay-1)/k(ilay-1)/kroi)*exp_04
            r5=(kyi*kzi(ilay)/kroi)*exp(-i*kzi(ilay)*dlay(ilay-1))*x0(4*ilay-1)
            r6=-(-(kxs-kxi)*kroi/k(ilay)+kxi*kzi(ilay)*kzi(ilay)/k(ilay)/kroi)*exp_05
            r7=-(kyi*kzi(ilay)/kroi)*exp_06
            r8=(-kxi*kzi(ilay)*kzi(ilay)/k(ilay)/kroi+(kxs-kxi)*kroi/k(ilay))*exp_07
            b(4*ilay-2)=i*(r1+r2+r3+r4+r5+r6+r7+r8)

            r1=(-(kys-kyi)*kroi/k(ilay-1)+kyi*kzi(ilay-1)*kzi(ilay-1)/k(ilay-1)/kroi)*exp_01/z(ilay-1)
            r2=-(kxi*kzi(ilay-1)/kroi)*exp_02/z(ilay-1)
            r3=(-(kys-kyi)*kroi/k(ilay-1)+kyi*kzi(ilay-1)*kzi(ilay-1)/k(ilay-1)/kroi)*exp_03/z(ilay-1)
            r4=(kxi*kzi(ilay-1)/kroi)*exp_04/z(ilay-1)
            r5=-(-(kys-kyi)*kroi/k(ilay)+kyi*kzi(ilay)*kzi(ilay)/k(ilay)/kroi) &
             *exp(-i*kzi(ilay)*dlay(ilay-1))*x0(4*ilay-1)/z(ilay)
            r6=(kxi*kzi(ilay)/kroi)*exp_05/z(ilay)
            r7=(-kyi*kzi(ilay)*kzi(ilay)/k(ilay)/kroi+(kys-kyi)*kroi/k(ilay))*exp_06/z(ilay)
            r8=-(kxi*kzi(ilay)/kroi)*exp_07/z(ilay)
            b(4*ilay-1)=i*(r1+r2+r3+r4+r5+r6+r7+r8)

            r1=((kxs-kxi)*kroi/k(ilay-1)-kxi*kzi(ilay-1)*kzi(ilay-1)/k(ilay-1)/kroi)*exp_01/z(ilay-1)
            r2=-(kyi*kzi(ilay-1)/kroi)*exp_02/z(ilay-1)
            r3=((kxs-kxi)*kroi/k(ilay-1)-kxi*kzi(ilay-1)*kzi(ilay-1)/k(ilay-1)/kroi)*exp_03/z(ilay-1)
            r4=(kyi*kzi(ilay-1)/kroi)*exp_04/z(ilay-1)
            r5=-((kxs-kxi)*kroi/k(ilay)-kxi*kzi(ilay)*kzi(ilay)/k(ilay)/kroi) &
             *exp(-i*kzi(ilay)*dlay(ilay-1))*x0(4*ilay-1)/z(ilay)
            r6=(kyi*kzi(ilay)/kroi)*exp_05/z(ilay)
            r7=(kxi*kzi(ilay)*kzi(ilay)/k(ilay)/kroi-(kxs-kxi)*kroi/k(ilay))*exp_06/z(ilay)
            r8=-(kyi*kzi(ilay)/kroi)*exp_07/z(ilay)
            b(4*ilay)=i*(r1+r2+r3+r4+r5+r6+r7+r8)
            
        end do

        r1=-(kxi*kzi(nlay-1)/kroi)*exp(-i*kzi(nlay-1)*dlay(nlay-1))*x0(4*nlay-5)
        r2=((kys-kyi)*kroi/k(nlay-1)-kyi*kzi(nlay-1)*kzi(nlay-1)/k(nlay-1)/kroi) &
           *exp(-i*kzi(nlay-1)*dlay(nlay-1))*x0(4*nlay-4)
        r3=(kxi*kzi(nlay-1)/kroi)*exp(i*kzi(nlay-1)*dlay(nlay-1))*x0(4*nlay-3)
        r4=((kys-kyi)*kroi/k(nlay-1)-kyi*kzi(nlay-1)*kzi(nlay-1)/k(nlay-1)/kroi) &
           *exp(i*kzi(nlay-1)*dlay(nlay-1))*x0(4*nlay-2)
        r5=-(kxi*kzi(nlay)/kroi)*exp(i*kzi(nlay)*dlay(nlay-1))*x0(4*nlay-1)
        r6=(kyi*kzi(nlay)*kzi(nlay)/k(nlay)/kroi-(kys-kyi)*kroi/k(nlay)) &
           *exp(i*kzi(nlay)*dlay(nlay-1))*x0(4*nlay)
        b(4*nlay-3)=i*(r1+r2+r3+r4+r5+r6)

        r1=-(kyi*kzi(nlay-1)/kroi)*exp(-i*kzi(nlay-1)*dlay(nlay-1))*x0(4*nlay-5)
        r2=(-(kxs-kxi)*kroi/k(nlay-1)+kxi*kzi(nlay-1)*kzi(nlay-1)/k(nlay-1)/kroi) &
           *exp(-i*kzi(nlay-1)*dlay(nlay-1))*x0(4*nlay-4)
        r3=(kyi*kzi(nlay-1)/kroi)*exp(i*kzi(nlay-1)*dlay(nlay-1))*x0(4*nlay-3)
        r4=(-(kxs-kxi)*kroi/k(nlay-1)+kxi*kzi(nlay-1)*kzi(nlay-1)/k(nlay-1)/kroi) &
           *exp(i*kzi(nlay-1)*dlay(nlay-1))*x0(4*nlay-2)
        r5=-(kyi*kzi(nlay)/kroi)*exp(i*kzi(nlay)*dlay(nlay-1))*x0(4*nlay-1)
        r6=(-kxi*kzi(nlay)*kzi(nlay)/k(nlay)/kroi+(kxs-kxi)*kroi/k(nlay)) &
           *exp(i*kzi(nlay)*dlay(nlay-1))*x0(4*nlay)
        b(4*nlay-2)=i*(r1+r2+r3+r4+r5+r6)

        r1=(-(kys-kyi)*kroi/k(nlay-1)+kyi*kzi(nlay-1)*kzi(nlay-1)/k(nlay-1)/kroi) &
           *exp(-i*kzi(nlay-1)*dlay(nlay-1))*x0(4*nlay-5)/z(nlay-1)
        r2=-(kxi*kzi(nlay-1)/kroi)*exp(-i*kzi(nlay-1)*dlay(nlay-1))*x0(4*nlay-4)/z(nlay-1)
        r3=(-(kys-kyi)*kroi/k(nlay-1)+kyi*kzi(nlay-1)*kzi(nlay-1)/k(nlay-1)/kroi) &
           *exp(i*kzi(nlay-1)*dlay(nlay-1))*x0(4*nlay-3)/z(nlay-1)
        r4=(kxi*kzi(nlay-1)/kroi)*exp(i*kzi(nlay-1)*dlay(nlay-1))*x0(4*nlay-2)/z(nlay-1)
        r5=(-kyi*kzi(nlay)*kzi(nlay)/k(nlay)/kroi+(kys-kyi)*kroi/k(nlay)) &
           *exp(i*kzi(nlay)*dlay(nlay-1))*x0(4*nlay-1)/z(nlay)
        r6=-(kxi*kzi(nlay)/kroi)*exp(i*kzi(nlay)*dlay(nlay-1))*x0(4*nlay)/z(nlay)
        b(4*nlay-1)=i*(r1+r2+r3+r4+r5+r6)

        r1=((kxs-kxi)*kroi/k(nlay-1)-kxi*kzi(nlay-1)*kzi(nlay-1)/k(nlay-1)/kroi) &
           *exp(-i*kzi(nlay-1)*dlay(nlay-1))*x0(4*nlay-5)/z(nlay-1)
        r2=-(kyi*kzi(nlay-1)/kroi)*exp(-i*kzi(nlay-1)*dlay(nlay-1))*x0(4*nlay-4)/z(nlay-1)
        r3=((kxs-kxi)*kroi/k(nlay-1)-kxi*kzi(nlay-1)*kzi(nlay-1)/k(nlay-1)/kroi) &
           *exp(i*kzi(nlay-1)*dlay(nlay-1))*x0(4*nlay-3)/z(nlay-1)
        r4=(kyi*kzi(nlay-1)/kroi)*exp(i*kzi(nlay-1)*dlay(nlay-1))*x0(4*nlay-2)/z(nlay-1)
        r5=(kxi*kzi(nlay)*kzi(nlay)/k(nlay)/kroi-(kxs-kxi)*kroi/k(nlay))*exp(i*kzi(nlay)*dlay(nlay-1))*x0(4*nlay-1)/z(nlay)
        r6=-(kyi*kzi(nlay)/kroi)*exp(i*kzi(nlay)*dlay(nlay-1))*x0(4*nlay)/z(nlay)
        b(4*nlay)=i*(r1+r2+r3+r4+r5+r6)

return
end subroutine e_rgt1_new

    SUBROUTINE E_RGT1(NL,AEM,AHM,DL,K0,K,KXI,KXS,KYI,KYS,Z0,Z,X0,B)

        IMPLICIT NONE

        INTEGER,         INTENT(IN)  :: NL
        REAL(KIND=8),    INTENT(IN)  :: AEM, AHM, DL(NL-1), K0, KXI, KXS, KYI, KYS, Z0
        COMPLEX(KIND=8), INTENT(IN)  :: K(NL), Z(NL), X0(4*NL)

        COMPLEX(KIND=8), INTENT(OUT) :: B(4*NL)

        REAL(KIND=8)    :: KROI, KZI0
        COMPLEX(KIND=8) :: I, KZI(NL)
        COMPLEX(KIND=8) :: R1, R2, R3, R4, R5, R6, R7, R8

        INTEGER :: II

        I=(0.D0,1.D0)

        KROI=SQRT(KXI**2+KYI**2)
        KZI0=SQRT(K0**2-KROI**2)

        DO II=1,NL
            KZI(II)=SQRT(K(II)**2-KROI**2)
        END DO

        R1=KZI0*KXI/KROI*X0(1)  
        R2=(KZI0*KZI0*KYI/K0/KROI-(KYS-KYI)*KROI/K0)*X0(2)
        R3=-(KZI(1)*KXI/KROI)*X0(3)
        R4=((KYS-KYI)*KROI/K(1)-KZI(1)**2*KYI/K(1)/KROI)*X0(4)
        R5=(KZI(1)*KXI/KROI)*X0(5)
        R6=((KYS-KYI)*KROI/K(1)-KYI*KZI(1)**2/K(1)/KROI)*X0(6)
        R7=-(KZI0*KXI/KROI)*AEM
        R8=(KYI*KZI0*KZI0/K0/KROI-(KYS-KYI)*KROI/K0)*AHM
        B(1)=-I*(R1+R2+R3+R4+R5+R6+R7+R8)

        R1=(KZI0*KYI/KROI)*X0(1)
        R2=(-KZI0*KZI0*KXI/K0/KROI+(KXS-KXI)*KROI/K0)*X0(2)
        R3=-(KZI(1)*KYI/KROI)*X0(3)
        R4=-((KXS-KXI)*KROI/K(1)-KZI(1)*KZI(1)*KXI/K(1)/KROI)*X0(4)
        R5=(KZI(1)*KYI/KROI)*X0(5)
        R6=-((KXS-KXI)*KROI/K(1)-KXI*KZI(1)*KZI(1)/K(1)/KROI)*X0(6)
        R7=-(KZI0*KYI/KROI)*AEM
        R8=(-KXI*KZI0*KZI0/K0/KROI+(KXS-KXI)*KROI/K0)*AHM
        B(2)=-I*(R1+R2+R3+R4+R5+R6+R7+R8)

        R1=(-(KYS-KYI)*KROI/K0+KZI0*KZI0*KYI/K0/KROI)*X0(1)/Z0
        R2=-(KXI*KZI0/KROI)*X0(2)/Z0
        R3=((KYS-KYI)*KROI/K(1)-KYI*KZI(1)*KZI(1)/K(1)/KROI)*X0(3)/Z(1)
        R4=(KXI*KZI(1)/KROI)*X0(4)/Z(1)
        R5=((KYS-KYI)*KROI/K(1)-KYI*KZI(1)*KZI(1)/K(1)/KROI)*X0(5)/Z(1)
        R6=-(KXI*KZI(1)/KROI)*X0(6)/Z(1)
        R7=(KYI*KZI0*KZI0/K0/KROI-(KYS-KYI)*KROI/K0)*AEM/Z0
        R8=(KXI*KZI0/KROI)*AHM/Z0
        B(3)=I*(R1+R2+R3+R4+R5+R6+R7+R8)

        R1=((KXS-KXI)*KROI/K0-KZI0*KZI0*KXI/K0/KROI)*X0(1)/Z0
        R2=-(KYI*KZI0/KROI)*X0(2)/Z0
        R3=(-(KXS-KXI)*KROI/K(1)+KXI*KZI(1)*KZI(1)/K(1)/KROI) &
        *X0(3)/Z(1)
        R4=(KYI*KZI(1)/KROI)*X0(4)/Z(1)
        R5=(-(KXS-KXI)*KROI/K(1)+KXI*KZI(1)*KZI(1)/K(1)/KROI) &
        *X0(5)/Z(1)
        R6=-(KYI*KZI(1)/KROI)*X0(6)/Z(1)
        R7=(-KXI*KZI0*KZI0/K0/KROI+(KXS-KXI)*KROI/K0)*AEM/Z0
        R8=(KYI*KZI0/KROI)*AHM/Z0
        B(4)=I*(R1+R2+R3+R4+R5+R6+R7+R8)


        DO II=2,NL-1
            R1=-(KXI*KZI(II-1)/KROI)*EXP(-I*KZI(II-1)*DL(II-1)) &
            *X0(4*II-5)
            R2=((KYS-KYI)*KROI/K(II-1)-KYI*KZI(II-1)*KZI(II-1)/K(II-1) &
            /KROI)*EXP(-I*KZI(II-1)*DL(II-1))*X0(4*II-4)
            R3=(KXI*KZI(II-1)/KROI)*EXP(I*KZI(II-1)*DL(II-1))*X0(4*II-3)
            R4=((KYS-KYI)*KROI/K(II-1)-KYI*KZI(II-1)*KZI(II-1)/K(II-1) &
            /KROI)*EXP(I*KZI(II-1)*DL(II-1))*X0(4*II-2)
            R5=(KXI*KZI(II)/KROI)*EXP(-I*KZI(II)*DL(II-1))*X0(4*II-1)
            R6=-((KYS-KYI)*KROI/K(II)-KYI*KZI(II)*KZI(II)/K(II)/KROI)* &
            EXP(-I*KZI(II)*DL(II-1))*X0(4*II)
            R7=-(KXI*KZI(II)/KROI)*EXP(I*KZI(II)*DL(II-1))*X0(4*II+1)
            R8=(KYI*KZI(II)*KZI(II)/K(II)/KROI-(KYS-KYI)*KROI/K(II))* &
            EXP(I*KZI(II)*DL(II-1))*X0(4*II+2)
            B(4*II-3)=I*(R1+R2+R3+R4+R5+R6+R7+R8)

            R1=-(KYI*KZI(II-1)/KROI)*EXP(-I*KZI(II-1)*DL(II-1)) &
            *X0(4*II-5)
            R2=(-(KXS-KXI)*KROI/K(II-1)+KXI*KZI(II-1)*KZI(II-1)/K(II-1) &
            /KROI)*EXP(-I*KZI(II-1)*DL(II-1))*X0(4*II-4)
            R3=(KYI*KZI(II-1)/KROI)*EXP(I*KZI(II-1)*DL(II-1))*X0(4*II-3)
            R4=(-(KXS-KXI)*KROI/K(II-1)+KXI*KZI(II-1)*KZI(II-1)/K(II-1) &
            /KROI)*EXP(I*KZI(II-1)*DL(II-1))*X0(4*II-2)
            R5=(KYI*KZI(II)/KROI)*EXP(-I*KZI(II)*DL(II-1))*X0(4*II-1)
            R6=-(-(KXS-KXI)*KROI/K(II)+KXI*KZI(II)*KZI(II)/K(II)/KROI)* &
            EXP(-I*KZI(II)*DL(II-1))*X0(4*II)
            R7=-(KYI*KZI(II)/KROI)*EXP(I*KZI(II)*DL(II-1))*X0(4*II+1)
            R8=(-KXI*KZI(II)*KZI(II)/K(II)/KROI+(KXS-KXI)*KROI/K(II))* &
            EXP(I*KZI(II)*DL(II-1))*X0(4*II+2)
            B(4*II-2)=I*(R1+R2+R3+R4+R5+R6+R7+R8)

            R1=(-(KYS-KYI)*KROI/K(II-1)+KYI*KZI(II-1)*KZI(II-1)/K(II-1) &
            /KROI)*EXP(-I*KZI(II-1)*DL(II-1))*X0(4*II-5)/Z(II-1)
            R2=-(KXI*KZI(II-1)/KROI)*EXP(-I*KZI(II-1)*DL(II-1)) &
            *X0(4*II-4)/Z(II-1)
            R3=(-(KYS-KYI)*KROI/K(II-1)+KYI*KZI(II-1)*KZI(II-1)/K(II-1) &
            /KROI)*EXP(I*KZI(II-1)*DL(II-1))*X0(4*II-3)/Z(II-1)
            R4=(KXI*KZI(II-1)/KROI)*EXP(I*KZI(II-1)*DL(II-1)) &
            *X0(4*II-2)/Z(II-1)
            R5=-(-(KYS-KYI)*KROI/K(II)+KYI*KZI(II)*KZI(II)/K(II)/KROI)* &
            EXP(-I*KZI(II)*DL(II-1))*X0(4*II-1)/Z(II)
            R6=(KXI*KZI(II)/KROI)*EXP(-I*KZI(II)*DL(II-1))*X0(4*II) &
            /Z(II)
            R7=(-KYI*KZI(II)*KZI(II)/K(II)/KROI+(KYS-KYI)*KROI/K(II))* &
            EXP(I*KZI(II)*DL(II-1))*X0(4*II+1)/Z(II)
            R8=-(KXI*KZI(II)/KROI)*EXP(I*KZI(II)*DL(II-1))*X0(4*II+2) &
            /Z(II)
            B(4*II-1)=I*(R1+R2+R3+R4+R5+R6+R7+R8)

            R1=((KXS-KXI)*KROI/K(II-1)-KXI*KZI(II-1)*KZI(II-1)/K(II-1) &
            /KROI)*EXP(-I*KZI(II-1)*DL(II-1))*X0(4*II-5)/Z(II-1)
            R2=-(KYI*KZI(II-1)/KROI)*EXP(-I*KZI(II-1)*DL(II-1)) &
            *X0(4*II-4)/Z(II-1)
            R3=((KXS-KXI)*KROI/K(II-1)-KXI*KZI(II-1)*KZI(II-1)/K(II-1) &
            /KROI)*EXP(I*KZI(II-1)*DL(II-1))*X0(4*II-3)/Z(II-1)
            R4=(KYI*KZI(II-1)/KROI)*EXP(I*KZI(II-1)*DL(II-1)) &
            *X0(4*II-2)/Z(II-1)
            R5=-((KXS-KXI)*KROI/K(II)-KXI*KZI(II)*KZI(II)/K(II)/KROI)* &
            EXP(-I*KZI(II)*DL(II-1))*X0(4*II-1)/Z(II)
            R6=(KYI*KZI(II)/KROI)*EXP(-I*KZI(II)*DL(II-1))*X0(4*II) &
            /Z(II)
            R7=(KXI*KZI(II)*KZI(II)/K(II)/KROI-(KXS-KXI)*KROI/K(II))* &
            EXP(I*KZI(II)*DL(II-1))*X0(4*II+1)/Z(II)
            R8=-(KYI*KZI(II)/KROI)*EXP(I*KZI(II)*DL(II-1))*X0(4*II+2) &
            /Z(II)
            B(4*II)=I*(R1+R2+R3+R4+R5+R6+R7+R8)
            
        END DO

        R1=-(KXI*KZI(NL-1)/KROI)*EXP(-I*KZI(NL-1)*DL(NL-1)) &
        *X0(4*NL-5)
        R2=((KYS-KYI)*KROI/K(NL-1)-KYI*KZI(NL-1)*KZI(NL-1)/K(NL-1) &
        /KROI)*EXP(-I*KZI(NL-1)*DL(NL-1))*X0(4*NL-4)
        R3=(KXI*KZI(NL-1)/KROI)*EXP(I*KZI(NL-1)*DL(NL-1))*X0(4*NL-3)
        R4=((KYS-KYI)*KROI/K(NL-1)-KYI*KZI(NL-1)*KZI(NL-1)/K(NL-1) &
        /KROI)*EXP(I*KZI(NL-1)*DL(NL-1))*X0(4*NL-2)
        R5=-(KXI*KZI(NL)/KROI)*EXP(I*KZI(NL)*DL(NL-1))*X0(4*NL-1)
        R6=(KYI*KZI(NL)*KZI(NL)/K(NL)/KROI-(KYS-KYI)*KROI/K(NL))* &
        EXP(I*KZI(NL)*DL(NL-1))*X0(4*NL)
        B(4*NL-3)=I*(R1+R2+R3+R4+R5+R6)

        R1=-(KYI*KZI(NL-1)/KROI)*EXP(-I*KZI(NL-1)*DL(NL-1)) &
        *X0(4*NL-5)
        R2=(-(KXS-KXI)*KROI/K(NL-1)+KXI*KZI(NL-1)*KZI(NL-1)/K(NL-1) &
        /KROI)*EXP(-I*KZI(NL-1)*DL(NL-1))*X0(4*NL-4)
        R3=(KYI*KZI(NL-1)/KROI)*EXP(I*KZI(NL-1)*DL(NL-1))*X0(4*NL-3)
        R4=(-(KXS-KXI)*KROI/K(NL-1)+KXI*KZI(NL-1)*KZI(NL-1)/K(NL-1) &
        /KROI)*EXP(I*KZI(NL-1)*DL(NL-1))*X0(4*NL-2)
        R5=-(KYI*KZI(NL)/KROI)*EXP(I*KZI(NL)*DL(NL-1))*X0(4*NL-1)
        R6=(-KXI*KZI(NL)*KZI(NL)/K(NL)/KROI+(KXS-KXI)*KROI/K(NL))* &
        EXP(I*KZI(NL)*DL(NL-1))*X0(4*NL)
        B(4*NL-2)=I*(R1+R2+R3+R4+R5+R6)

        R1=(-(KYS-KYI)*KROI/K(NL-1)+KYI*KZI(NL-1)*KZI(NL-1)/K(NL-1) &
        /KROI)*EXP(-I*KZI(NL-1)*DL(NL-1))*X0(4*NL-5)/Z(NL-1)
        R2=-(KXI*KZI(NL-1)/KROI)*EXP(-I*KZI(NL-1)*DL(NL-1))* &
        X0(4*NL-4)/Z(NL-1)
        R3=(-(KYS-KYI)*KROI/K(NL-1)+KYI*KZI(NL-1)*KZI(NL-1)/K(NL-1) &
        /KROI)*EXP(I*KZI(NL-1)*DL(NL-1))*X0(4*NL-3)/Z(NL-1)
        R4=(KXI*KZI(NL-1)/KROI)*EXP(I*KZI(NL-1)*DL(NL-1))*X0(4*NL-2) &
        /Z(NL-1)
        R5=(-KYI*KZI(NL)*KZI(NL)/K(NL)/KROI+(KYS-KYI)*KROI/K(NL))* &
        EXP(I*KZI(NL)*DL(NL-1))*X0(4*NL-1)/Z(NL)
        R6=-(KXI*KZI(NL)/KROI)*EXP(I*KZI(NL)*DL(NL-1))*X0(4*NL)/Z(NL)
        B(4*NL-1)=I*(R1+R2+R3+R4+R5+R6)

        R1=((KXS-KXI)*KROI/K(NL-1)-KXI*KZI(NL-1)*KZI(NL-1)/K(NL-1) &
        /KROI)*EXP(-I*KZI(NL-1)*DL(NL-1))*X0(4*NL-5)/Z(NL-1)
        R2=-(KYI*KZI(NL-1)/KROI)*EXP(-I*KZI(NL-1)*DL(NL-1))* &
        X0(4*NL-4)/Z(NL-1)
        R3=((KXS-KXI)*KROI/K(NL-1)-KXI*KZI(NL-1)*KZI(NL-1)/K(NL-1) &
        /KROI)*EXP(I*KZI(NL-1)*DL(NL-1))*X0(4*NL-3)/Z(NL-1)
        R4=(KYI*KZI(NL-1)/KROI)*EXP(I*KZI(NL-1)*DL(NL-1))*X0(4*NL-2) &
        /Z(NL-1)
        R5=(KXI*KZI(NL)*KZI(NL)/K(NL)/KROI-(KXS-KXI)*KROI/K(NL))* &
        EXP(I*KZI(NL)*DL(NL-1))*X0(4*NL-1)/Z(NL)
        R6=-(KYI*KZI(NL)/KROI)*EXP(I*KZI(NL)*DL(NL-1))*X0(4*NL)/Z(NL)
        B(4*NL)=I*(R1+R2+R3+R4+R5+R6)


        RETURN
    END
!
!--------------------------------------------------------------------------------
subroutine calculate_coefficients(nlay,dlay,k0,k,kx,ky,z0,z,a)
!--------------------------------------------------------------------------------
! calculates the elements of the coefficient matrix in ax=b based on formulation
! in the paper and saves them in matrix a.

implicit none
!
! define input variables
integer,         intent(in)  :: nlay
real(kind=8),    intent(in)  :: dlay(nlay-1)
real(kind=8),    intent(in)  :: k0
real(kind=8),    intent(in)  :: kx
real(kind=8),    intent(in)  :: ky
real(kind=8),    intent(in)  :: z0
complex(kind=8), intent(in)  :: k(nlay) ! (?) Wavenumber per layer
complex(kind=8), intent(in)  :: z(nlay) ! Intrinsic impedance per layer
!
! define output variable
complex(kind=8), intent(out) :: a(4*nlay,4*nlay)
!
! define local variables
real(kind=8)    :: kro
complex(kind=8) :: i
complex(kind=8) :: kz0
complex(kind=8) :: kz(nlay)
complex(kind=8) :: exp_01
complex(kind=8) :: exp_02
complex(kind=8) :: exp_03
complex(kind=8) :: exp_04
integer :: ilay
integer :: jj

        i=(0.d0,1.d0)

        kro=sqrt(kx**2+ky**2)
        kz0=sqrt(k0**2-kro**2)
        do ilay=1,nlay
            kz(ilay)=sqrt(k(ilay)**2-kro**2)
        end do

        a(1,1)=kx/kro 
        a(1,2)=kz0*ky/k0/kro 
        a(1,3)=-kx/kro 
        a(1,4)=-kz(1)*ky/k(1)/kro 
        a(1,5)=-kx/kro
        a(1,6)=ky*kz(1)/k(1)/kro
        
        a(2,1)=ky/kro 
        a(2,2)=-kz0*kx/k0/kro 
        a(2,3)=-ky/kro 
        a(2,4)=kz(1)*kx/k(1)/kro 
        a(2,5)=-ky/kro 
        a(2,6)=-kx*kz(1)/k(1)/kro
        
        a(3,1)=-kz0*ky/k0/kro/z0
        a(3,2)=kx/kro/z0
        a(3,3)=kz(1)*ky/k(1)/kro/z(1)
        a(3,4)=-kx/kro/z(1)
        a(3,5)=-ky*kz(1)/k(1)/kro/z(1)
        a(3,6)=-kx/kro/z(1)
        
        a(4,1)=kz0*kx/k0/kro/z0
        a(4,2)=ky/kro/z0
        a(4,3)=-kz(1)*kx/k(1)/kro/z(1)
        a(4,4)=-ky/kro/z(1)
        a(4,5)=kx*kz(1)/k(1)/kro/z(1)
        a(4,6)=-ky/kro/z(1)
        
        do ilay=1,4
            do jj=7,4*nlay
                a(ilay,jj)=(0.d0,0.d0)
            enddo
        end do


        do ilay=2,nlay-1

            do jj=1,4*ilay-6
                a(4*ilay-3,jj)=(0.d0,0.d0)
                a(4*ilay-2,jj)=(0.d0,0.d0)
                a(4*ilay-1,jj)=(0.d0,0.d0)
                a(4*ilay,jj)=(0.d0,0.d0)
            enddo
            exp_01=exp(-i*kz(ilay-1)*dlay(ilay-1))
            exp_02=exp(i*kz(ilay-1)*dlay(ilay-1))
            exp_03=exp(-i*kz(ilay)*dlay(ilay-1))
            exp_04=exp(i*kz(ilay)*dlay(ilay-1))

            a(4*ilay-3,4*ilay-5)=kx/kro*exp_01
            a(4*ilay-3,4*ilay-4)=kz(ilay-1)*ky/k(ilay-1)/kro*exp_01
            a(4*ilay-3,4*ilay-3)=kx/kro*exp_02
            a(4*ilay-3,4*ilay-2)=-ky*kz(ilay-1)/k(ilay-1)/kro*exp_02
            a(4*ilay-3,4*ilay-1)=-kx/kro*exp_03
            a(4*ilay-3,4*ilay)=-kz(ilay)*ky/k(ilay)/kro*exp_03
            a(4*ilay-3,4*ilay+1)=-kx/kro*exp_04
            a(4*ilay-3,4*ilay+2)=ky*kz(ilay)/k(ilay)/kro*exp_04

            a(4*ilay-2,4*ilay-5)=ky/kro*exp_01
            a(4*ilay-2,4*ilay-4)=-kz(ilay-1)*kx/k(ilay-1)/kro*exp_01
            a(4*ilay-2,4*ilay-3)=ky/kro*exp_02
            a(4*ilay-2,4*ilay-2)=kx*kz(ilay-1)/k(ilay-1)/kro*exp_02
            a(4*ilay-2,4*ilay-1)=-ky/kro*exp_03
            a(4*ilay-2,4*ilay)=kz(ilay)*kx/k(ilay)/kro*exp_03
            a(4*ilay-2,4*ilay+1)=-ky/kro*exp_04
            a(4*ilay-2,4*ilay+2)=-kx*kz(ilay)/k(ilay)/kro*exp_04
                
            a(4*ilay-1,4*ilay-5)=-ky*kz(ilay-1)/k(ilay-1)/kro/z(ilay-1)*exp_01
            a(4*ilay-1,4*ilay-4)=kx/kro/z(ilay-1)*exp_01
            a(4*ilay-1,4*ilay-3)=ky*kz(ilay-1)/k(ilay-1)/kro/z(ilay-1)*exp_02
            a(4*ilay-1,4*ilay-2)=kx/kro/z(ilay-1)*exp_02
            a(4*ilay-1,4*ilay-1)=ky*kz(ilay)/k(ilay)/kro/z(ilay)*exp_03
            a(4*ilay-1,4*ilay)=-kx/kro/z(ilay)*exp_03
            a(4*ilay-1,4*ilay+1)=-ky*kz(ilay)/k(ilay)/kro/z(ilay)*exp_04
            a(4*ilay-1,4*ilay+2)=-kx/kro/z(ilay)*exp_04      

            a(4*ilay,4*ilay-5)=kx*kz(ilay-1)/k(ilay-1)/kro/z(ilay-1)*exp_01
            a(4*ilay,4*ilay-4)=ky/kro/z(ilay-1)*exp_01
            a(4*ilay,4*ilay-3)=-kx*kz(ilay-1)/k(ilay-1)/kro/z(ilay-1)*exp_02
            a(4*ilay,4*ilay-2)=ky/kro/z(ilay-1)*exp_02
            a(4*ilay,4*ilay-1)=-kx*kz(ilay)/k(ilay)/kro/z(ilay)*exp_03
            a(4*ilay,4*ilay)=-ky/kro/z(ilay)*exp_03
            a(4*ilay,4*ilay+1)=kx*kz(ilay)/k(ilay)/kro/z(ilay)*exp_04
            a(4*ilay,4*ilay+2)=-ky/kro/z(ilay)*exp_04

            do jj=4*ilay+3,4*nlay
                a(4*ilay-3,jj)=(0.d0,0.d0)
                a(4*ilay-2,jj)=(0.d0,0.d0)
                a(4*ilay-1,jj)=(0.d0,0.d0)
                a(4*ilay,jj)=(0.d0,0.d0)
            enddo
        enddo

        do ilay=4*nlay-3,4*nlay
            do jj=1,4*nlay-6
                a(ilay,jj)=(0.d0,0.d0)
            enddo
        enddo
        
        exp_01=exp(-i*kz(nlay-1)*dlay(nlay-1))
        exp_02=exp(i*kz(nlay-1)*dlay(nlay-1))
        exp_03=exp(i*kz(nlay)*dlay(nlay-1))

        a(4*nlay-3,4*nlay-5)=kx/kro*exp_01
        a(4*nlay-3,4*nlay-4)=kz(nlay-1)*ky/k(nlay-1)/kro*exp_01
        a(4*nlay-3,4*nlay-3)=kx/kro*exp_02
        a(4*nlay-3,4*nlay-2)=-ky*kz(nlay-1)/k(nlay-1)/kro*exp_02
        a(4*nlay-3,4*nlay-1)=-kx/kro*exp_03
        a(4*nlay-3,4*nlay)=ky*kz(nlay)/k(nlay)/kro*exp_03

        a(4*nlay-2,4*nlay-5)=ky/kro*exp_01
        a(4*nlay-2,4*nlay-4)=-kz(nlay-1)*kx/k(nlay-1)/kro*exp_01
        a(4*nlay-2,4*nlay-3)=ky/kro*exp_02
        a(4*nlay-2,4*nlay-2)=kx*kz(nlay-1)/k(nlay-1)/kro*exp_02
        a(4*nlay-2,4*nlay-1)=-ky/kro*exp_03
        a(4*nlay-2,4*nlay)=-kx*kz(nlay)/k(nlay)/kro*exp_03
        
        a(4*nlay-1,4*nlay-5)=-ky*kz(nlay-1)/k(nlay-1)/kro/z(nlay-1)*exp_01
        a(4*nlay-1,4*nlay-4)=kx/kro/z(nlay-1)*exp_01
        a(4*nlay-1,4*nlay-3)=ky*kz(nlay-1)/k(nlay-1)/kro/z(nlay-1)* exp_02
        a(4*nlay-1,4*nlay-2)=kx/kro/z(nlay-1)*exp_02
        a(4*nlay-1,4*nlay-1)=-ky*kz(nlay)/k(nlay)/kro/z(nlay)*exp_03
        a(4*nlay-1,4*nlay)=-kx/kro/z(nlay)*exp_03

        a(4*nlay,4*nlay-5)=kx*kz(nlay-1)/k(nlay-1)/kro/z(nlay-1)*exp_01
        a(4*nlay,4*nlay-4)=ky/kro/z(nlay-1)*exp_01
        a(4*nlay,4*nlay-3)=-kx*kz(nlay-1)/k(nlay-1)/kro/z(nlay-1)*exp_02
        a(4*nlay,4*nlay-2)=ky/kro/z(nlay-1)*exp_02
        a(4*nlay,4*nlay-1)=kx*kz(nlay)/k(nlay)/kro/z(nlay)*exp_03
        a(4*nlay,4*nlay)=-ky/kro/z(nlay)*exp_03

        return
    end
!
!--------------------------------------------------------------------------------
subroutine nlayerspm(nl,freq,thetainc,ls,sh,dlay,diel,sigma_vv,sigma_hh, test, test10)
!--------------------------------------------------------------------------------
implicit none
! define input variables
integer,         intent(in)  :: nl       ! (-) number of soil layers
real(kind=8),    intent(in)  :: freq     ! (hz) frequency
real(kind=8),    intent(in)  :: thetainc ! (deg) scattering angle
real(kind=8),    intent(in)  :: ls(nl)   ! (?) correlation lengths of layer boundaries
real(kind=8),    intent(in)  :: sh(nl)   ! (?) standard deviation of layer boundaries
real(kind=8),    intent(in)  :: dlay(nl-1) ! (m) average z coordinate of layer boundaries
complex(kind=8), intent(in)  :: diel(nl) ! (-) relative dielectric constant of layers.
real test, test10(10)     ! test variables
!
! define output variables
real(kind=8),    intent(out) :: sigma_vv ! (-) Scattering coefficients of the layered structure for VV polarization
real(kind=8),    intent(out) :: sigma_hh ! (-) Scattering coefficients of the layered structure for HH polarization
!
    call airmoss_fwd_model(nl,freq,thetainc,ls,sh,dlay,diel,sigma_vv,sigma_hh, test, test10)

return
end subroutine nlayerspm

    SUBROUTINE EQU_COEFF(NL,DL,K0,K,KX,KY,Z0,Z,A)

        IMPLICIT NONE
        
        INTEGER,         INTENT(IN)  :: NL
        REAL(KIND=8),    INTENT(IN)  :: DL(NL-1), K0, KX, KY, Z0
        COMPLEX(KIND=8), INTENT(IN)  :: K(NL), Z(NL)

        COMPLEX(KIND=8), INTENT(OUT) :: A(4*NL,4*NL)

        REAL(KIND=8)    :: KRO
        COMPLEX(KIND=8) :: I, KZ0, KZ(NL)

        INTEGER :: II, JJ

        I=(0.D0,1.D0)

        KRO=SQRT(KX**2+KY**2)
        KZ0=SQRT(K0**2-KRO**2)
        DO II=1,NL
            KZ(II)=SQRT(K(II)**2-KRO**2)
        END DO

        A(1,1)=KX/KRO 
        A(1,2)=KZ0*KY/K0/KRO 
        A(1,3)=-KX/KRO 
        A(1,4)=-KZ(1)*KY/K(1)/KRO 
        A(1,5)=-KX/KRO
        A(1,6)=KY*KZ(1)/K(1)/KRO
        
        A(2,1)=KY/KRO 
        A(2,2)=-KZ0*KX/K0/KRO 
        A(2,3)=-KY/KRO 
        A(2,4)=KZ(1)*KX/K(1)/KRO 
        A(2,5)=-KY/KRO 
        A(2,6)=-KX*KZ(1)/K(1)/KRO
        
        A(3,1)=-KZ0*KY/K0/KRO/Z0
        A(3,2)=KX/KRO/Z0
        A(3,3)=KZ(1)*KY/K(1)/KRO/Z(1)
        A(3,4)=-KX/KRO/Z(1)
        A(3,5)=-KY*KZ(1)/K(1)/KRO/Z(1)
        A(3,6)=-KX/KRO/Z(1)
        
        A(4,1)=KZ0*KX/K0/KRO/Z0
        A(4,2)=KY/KRO/Z0
        A(4,3)=-KZ(1)*KX/K(1)/KRO/Z(1)
        A(4,4)=-KY/KRO/Z(1)
        A(4,5)=KX*KZ(1)/K(1)/KRO/Z(1)
        A(4,6)=-KY/KRO/Z(1)
        
        DO II=1,4
            DO JJ=7,4*NL
                A(II,JJ)=(0.D0,0.D0)
            ENDDO
        END DO


        DO II=2,NL-1

            DO JJ=1,4*II-6
                A(4*II-3,JJ)=(0.D0,0.D0)
                A(4*II-2,JJ)=(0.D0,0.D0)
                A(4*II-1,JJ)=(0.D0,0.D0)
                A(4*II,JJ)=(0.D0,0.D0)
            ENDDO

            A(4*II-3,4*II-5)=KX/KRO*EXP(-I*KZ(II-1)*DL(II-1))
            A(4*II-3,4*II-4)=KZ(II-1)*KY/K(II-1)/KRO*EXP(-I*KZ(II-1) &
            *DL(II-1))
            A(4*II-3,4*II-3)=KX/KRO*EXP(I*KZ(II-1)*DL(II-1))
            A(4*II-3,4*II-2)=-KY*KZ(II-1)/K(II-1)/KRO*EXP(I*KZ(II-1) &
            *DL(II-1))
            A(4*II-3,4*II-1)=-KX/KRO*EXP(-I*KZ(II)*DL(II-1))
            A(4*II-3,4*II)=-KZ(II)*KY/K(II)/KRO*EXP(-I*KZ(II)*DL(II-1))
            A(4*II-3,4*II+1)=-KX/KRO*EXP(I*KZ(II)*DL(II-1))
            A(4*II-3,4*II+2)=KY*KZ(II)/K(II)/KRO*EXP(I*KZ(II)*DL(II-1))

            A(4*II-2,4*II-5)=KY/KRO*EXP(-I*KZ(II-1)*DL(II-1))
            A(4*II-2,4*II-4)=-KZ(II-1)*KX/K(II-1)/KRO*EXP(-I*KZ(II-1) &
            *DL(II-1))
            A(4*II-2,4*II-3)=KY/KRO*EXP(I*KZ(II-1)*DL(II-1))
            A(4*II-2,4*II-2)=KX*KZ(II-1)/K(II-1)/KRO*EXP(I*KZ(II-1) &
            *DL(II-1))
            A(4*II-2,4*II-1)=-KY/KRO*EXP(-I*KZ(II)*DL(II-1))
            A(4*II-2,4*II)=KZ(II)*KX/K(II)/KRO*EXP(-I*KZ(II)*DL(II-1))
            A(4*II-2,4*II+1)=-KY/KRO*EXP(I*KZ(II)*DL(II-1))
            A(4*II-2,4*II+2)=-KX*KZ(II)/K(II)/KRO*EXP(I*KZ(II)*DL(II-1))
                
            A(4*II-1,4*II-5)=-KY*KZ(II-1)/K(II-1)/KRO/Z(II-1)* &
            EXP(-I*KZ(II-1)*DL(II-1))
            A(4*II-1,4*II-4)=KX/KRO/Z(II-1)*EXP(-I*KZ(II-1)*DL(II-1))
            A(4*II-1,4*II-3)=KY*KZ(II-1)/K(II-1)/KRO/Z(II-1)* &
            EXP(I*KZ(II-1)*DL(II-1))
            A(4*II-1,4*II-2)=KX/KRO/Z(II-1)*EXP(I*KZ(II-1)*DL(II-1))
            A(4*II-1,4*II-1)=KY*KZ(II)/K(II)/KRO/Z(II)* &
            EXP(-I*KZ(II)*DL(II-1))
            A(4*II-1,4*II)=-KX/KRO/Z(II)*EXP(-I*KZ(II)*DL(II-1))
            A(4*II-1,4*II+1)=-KY*KZ(II)/K(II)/KRO/Z(II)* &
            EXP(I*KZ(II)*DL(II-1))
            A(4*II-1,4*II+2)=-KX/KRO/Z(II)*EXP(I*KZ(II)*DL(II-1))      

            A(4*II,4*II-5)=KX*KZ(II-1)/K(II-1)/KRO/Z(II-1)* &
            EXP(-I*KZ(II-1)*DL(II-1))
            A(4*II,4*II-4)=KY/KRO/Z(II-1)*EXP(-I*KZ(II-1)*DL(II-1))
            A(4*II,4*II-3)=-KX*KZ(II-1)/K(II-1)/KRO/Z(II-1)* &
            EXP(I*KZ(II-1)*DL(II-1))
            A(4*II,4*II-2)=KY/KRO/Z(II-1)*EXP(I*KZ(II-1)*DL(II-1))
            A(4*II,4*II-1)=-KX*KZ(II)/K(II)/KRO/Z(II)* &
            EXP(-I*KZ(II)*DL(II-1))
            A(4*II,4*II)=-KY/KRO/Z(II)*EXP(-I*KZ(II)*DL(II-1))
            A(4*II,4*II+1)=KX*KZ(II)/K(II)/KRO/Z(II)* &
            EXP(I*KZ(II)*DL(II-1))
            A(4*II,4*II+2)=-KY/KRO/Z(II)*EXP(I*KZ(II)*DL(II-1))

            DO JJ=4*II+3,4*NL
                A(4*II-3,JJ)=(0.D0,0.D0)
                A(4*II-2,JJ)=(0.D0,0.D0)
                A(4*II-1,JJ)=(0.D0,0.D0)
                A(4*II,JJ)=(0.D0,0.D0)
            ENDDO
        ENDDO

        DO II=4*NL-3,4*NL
            DO JJ=1,4*NL-6
                A(II,JJ)=(0.D0,0.D0)
            ENDDO
        ENDDO

        A(4*NL-3,4*NL-5)=KX/KRO*EXP(-I*KZ(NL-1)*DL(NL-1))
        A(4*NL-3,4*NL-4)=KZ(NL-1)*KY/K(NL-1)/KRO*EXP(-I*KZ(NL-1) &
        *DL(NL-1))
        A(4*NL-3,4*NL-3)=KX/KRO*EXP(I*KZ(NL-1)*DL(NL-1))
        A(4*NL-3,4*NL-2)=-KY*KZ(NL-1)/K(NL-1)/KRO*EXP(I*KZ(NL-1) &
        *DL(NL-1))
        A(4*NL-3,4*NL-1)=-KX/KRO*EXP(I*KZ(NL)*DL(NL-1))
        A(4*NL-3,4*NL)=KY*KZ(NL)/K(NL)/KRO*EXP(I*KZ(NL)*DL(NL-1))

        A(4*NL-2,4*NL-5)=KY/KRO*EXP(-I*KZ(NL-1)*DL(NL-1))
        A(4*NL-2,4*NL-4)=-KZ(NL-1)*KX/K(NL-1)/KRO*EXP(-I*KZ(NL-1) &
        *DL(NL-1))
        A(4*NL-2,4*NL-3)=KY/KRO*EXP(I*KZ(NL-1)*DL(NL-1))
        A(4*NL-2,4*NL-2)=KX*KZ(NL-1)/K(NL-1)/KRO*EXP(I*KZ(NL-1) &
        *DL(NL-1))
        A(4*NL-2,4*NL-1)=-KY/KRO*EXP(I*KZ(NL)*DL(NL-1))
        A(4*NL-2,4*NL)=-KX*KZ(NL)/K(NL)/KRO*EXP(I*KZ(NL)*DL(NL-1))
        
        A(4*NL-1,4*NL-5)=-KY*KZ(NL-1)/K(NL-1)/KRO/Z(NL-1)* &
        EXP(-I*KZ(NL-1)*DL(NL-1))
        A(4*NL-1,4*NL-4)=KX/KRO/Z(NL-1)*EXP(-I*KZ(NL-1)*DL(NL-1))
        A(4*NL-1,4*NL-3)=KY*KZ(NL-1)/K(NL-1)/KRO/Z(NL-1)* &
        EXP(I*KZ(NL-1)*DL(NL-1))
        A(4*NL-1,4*NL-2)=KX/KRO/Z(NL-1)*EXP(I*KZ(NL-1)*DL(NL-1))
        A(4*NL-1,4*NL-1)=-KY*KZ(NL)/K(NL)/KRO/Z(NL)* &
        EXP(I*KZ(NL)*DL(NL-1))
        A(4*NL-1,4*NL)=-KX/KRO/Z(NL)*EXP(I*KZ(NL)*DL(NL-1))

        A(4*NL,4*NL-5)=KX*KZ(NL-1)/K(NL-1)/KRO/Z(NL-1)* &
        EXP(-I*KZ(NL-1)*DL(NL-1))
        A(4*NL,4*NL-4)=KY/KRO/Z(NL-1)*EXP(-I*KZ(NL-1)*DL(NL-1))
        A(4*NL,4*NL-3)=-KX*KZ(NL-1)/K(NL-1)/KRO/Z(NL-1)* &
        EXP(I*KZ(NL-1)*DL(NL-1))
        A(4*NL,4*NL-2)=KY/KRO/Z(NL-1)*EXP(I*KZ(NL-1)*DL(NL-1))
        A(4*NL,4*NL-1)=KX*KZ(NL)/K(NL)/KRO/Z(NL)* &
        EXP(I*KZ(NL)*DL(NL-1))
        A(4*NL,4*NL)=-KY/KRO/Z(NL)*EXP(I*KZ(NL)*DL(NL-1))

        RETURN
    END
