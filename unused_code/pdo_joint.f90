!
!=======================================================================
      SUBROUTINE saturation_fraction(depth_top, depth_bot, vwc, &
      sandfrac, kroot, mass_org, root_depth, org_depth, rho_om_max, poros_om, &
      satfrac, test, test10)
!=======================================================================
! Calculate depth-averaged saturation fraction for a soil layer
! given volumetric water content
!
! Modifications:
!  Kevin Schaefer made subroutine (10/23/17)
!  Kevin Schaefer added primary output variable (10/23/17)
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
! Input variables
    real depth_top  ! (m) depth from surface to top of soil layer
    real depth_bot  ! (m) depth from surface to bottom of soil layer
    real vwc        ! (-) volumetric water content
!
! input constants
    real sandfrac   ! (%) sand fraction of soil texture
    real kRoot      ! (1/m) Exp const for decrease root-zone organic matter w/ depth
    real mass_org   ! (kg/m2) mass of organic matter in top 1 m of soil
    real root_depth ! (m) maximum rooting depth
    real org_depth  ! (m) thickness of organic layer
    real rho_om_max ! (kg/m3) maximum organic matter density
    real poros_om   ! (-) soil porosity for pure organic soil
!
! Output variables
    real satfrac    ! (-) saturation fraction of soil porosity using mixed model
    real test, test10(10) ! test variables
!
! Local variables
    integer numsteps  ! (-) number of integration steps
    real satfrac_mix  ! (-) saturation fraction mix of organic and mineral soil
    real satfrac_min  ! (-) saturation fraction pure mineral soil
    real satfrac_lay  ! (-) saturation fraction organic soil layer
    integer idep      ! (-) depth index
    real dz           ! (m) incremental active layer thickness
    real depth_loc    ! (m) local depth in integral
    real poros        ! (-) porosity for organic mineral soil mix
    real poros_min    ! (-) mineral soil porosity
    real apor_mix     ! (-) ave porosity mix of organic and mineral soil
    real apor_min     ! (-) ave porosity pure mineral soil
    real apor_lay     ! (-) ave porosity organic soil layer
    real orgfrac      ! (-) organic soil fraction
!
! error check
    if(depth_top>depth_bot) then
      print*, 'error: depth_top>depth_bot'
      print*, depth_top,depth_bot
      return
    endif
!
! calculate mineral soil porosity
    poros_min=0.489-0.00126*Sandfrac
!
!-----------------------------------------------------------------------
! saturation fraction for pure mineral soil
!-----------------------------------------------------------------------
    apor_min=poros_min
    satfrac_min=vwc/apor_min
    satfrac_min=min(1.0,satfrac_min)
!
!-----------------------------------------------------------------------
! saturation fraction for mix of organic and mineral soil
!-----------------------------------------------------------------------
! initialize
    numsteps=1000
    dz=(depth_bot-depth_top)/real(numsteps)
    depth_loc=depth_top
    satfrac_mix=0.
    apor_mix=0.
!
! numerical integration to calculate total pore space for mix of organic and mineral soil 
    do idep=1,numsteps
!
! organic soil fraction
      orgfrac=kRoot*mass_org/(1.-exp(-kRoot*root_depth))*exp(-kRoot*(depth_loc+0.5*dz))/rho_om_max
      if(orgfrac>1.) orgfrac=1.
      if(orgfrac<0.) orgfrac=0.
      if(depth_loc>root_depth) orgfrac=0.
!
! soil porosity
      poros=(1.-orgfrac)*poros_min+orgfrac*poros_om
      if(poros>1.) poros=1.
!
! increment total pore space
      apor_mix=apor_mix+poros*dz
!
! increment depth
      depth_loc=depth_loc+dz
    enddo
!
! convert total pore space to average pore space
    apor_mix=apor_mix/(depth_bot-depth_top)
!
! calculate saturation fraction
    satfrac_mix=vwc/apor_mix
    satfrac_mix=min(1.0,satfrac_mix)
!
!-----------------------------------------------------------------------
! saturation fraction for organic soil layer
!-----------------------------------------------------------------------
! organic soil fraction
    orgfrac=mass_org/org_depth/rho_om_max
    if(orgfrac>1.) orgfrac=1.
    if(orgfrac<0.) orgfrac=0.
!
! soil porosity
    poros=(1.-orgfrac)*poros_min+orgfrac*poros_om
    if(poros>1.) poros=1.
!
! Calculate average total pore space
    if(depth_top<org_depth.and.depth_bot<=org_depth) then !  all in organic layer
      apor_lay=poros

    elseif(depth_top<=org_depth.and.depth_bot>=org_depth) then ! part in and part out of organic layer
        apor_lay=poros*(org_depth-depth_top)/(depth_top-depth_bot)
        apor_lay=apor_lay+poros_min*(depth_bot-org_depth)/(depth_top-depth_bot)

    elseif(depth_top>org_depth) then ! all below organic layer
      apor_lay=poros_min
    endif

    satfrac_lay=vwc/apor_lay
    satfrac_lay=min(1.0,satfrac_lay)
!
! assign test variables
    test=satfrac_mix
    test10(1)=satfrac_mix
    test10(2)=satfrac_lay
    test10(3)=satfrac_min
!
! assign primary output variable
    satfrac=satfrac_mix

    RETURN
    END
!
!=======================================================================
    SUBROUTINE alt_from_deform_joint(depth_top,depth_bot,vwc,deform,&
    sandfrac, kroot, mass_org, root_depth, org_depth, rho_om_max, poros_om,&
    satfrac, alt_wat, alt_org_exp, alt_org_lay, alt_min, test, test10)
!=======================================================================
! Joint ReSALT/Airmoss retrieval: Calculate the active layer thickness
! given surface deformation and volumetric water content
!
! Modifications:
!  Kevin Schaefer made subroutine (11/5/17)
!-----------------------------------------------------------------------
!
    IMPLICIT NONE
!
! Input variables
    real depth_top  ! (m) depth from surface to top of soil layer
    real depth_bot  ! (m) depth from surface to bottom of soil layer
    real vwc        ! (-) volumetric water content
    real deform     ! (cm) surface deformation
!
! input constants
    real sandfrac   ! (%) sand fraction of soil texture
    real kRoot      ! (1/m) Exp const for decrease root-zone organic matter w/ depth
    real mass_org   ! (kg/m2) mass of organic matter in top 1 m of soil
    real root_depth ! (m) maximum rooting depth
    real org_depth  ! (m) thickness of organic layer
    real rho_om_max ! (kg/m3) maximum organic matter density
    real poros_om   ! (-) soil porosity for pure organic soil
!
! Output variables
    real satfrac     ! (-) saturation fraction of soil porosity
    real alt_wat     ! (cm) active layer thickness pure water
    real alt_org_exp ! (cm) active layer thickness exponential mix of organic and mineral soil
    real alt_org_lay ! (cm) active layer thickness organic soil layer
    real alt_min     ! (cm) active layer thickness pure mineral soil
    real test, test10(10) ! (TBD) dummy variables for testing algorithm
!
! Local variables
    integer iter    ! (-) iteration index
    integer n_iter  ! (-) number iterations
    real dz         ! (m) incremental active layer thickness
    real depth      ! (m) current depth
    real def_test   ! (cm) current deformation
    real del_def    ! (cm) change in deformation
    real def_max    ! (cm) deformation for layer of organic soil
    real rho_water  ! (kg/m3) density of water
    real rho_ice    ! (kg/m3) density of ice
    real poros      ! (-) porosity for organic mineral soil mix
    real poros_min  ! (-) mineral soil porosity
    real orgfrac    ! (-) organic soil fraction
    real fsat       ! (-) soil air space saturation scaling factor
    real sat_max    ! (-) soil saturation fraction where fsat is a maximum
    real sat_zero   ! (-) soil saturation fraction where fsat is zero
    real sat_eff    ! (-) effective soil saturation fraction accounting for expansion into pore air
!
! execution control variables
    logical :: flg_fsat=.false. ! (-) apply saturation fraction scaling factor
!
!-----------------------------------------------------------------------
! Assign constants
!-----------------------------------------------------------------------
! assign ice and water densities
    rho_water=1000.
    rho_ice=917.
!
! vertical integration delta soil depth
! tests indicate that this value produces minimal error
    dz=0.01
!
! calculate mineral soil porosity
   poros_min=0.489-0.00126*Sandfrac
!
!-----------------------------------------------------------------------
! saturation fraction
!-----------------------------------------------------------------------
! Convert depth average vwc into saturation fraction
      call saturation_fraction(depth_top, depth_bot, vwc,&
      sandfrac, kroot, mass_org, root_depth, org_depth, rho_om_max, poros_om,&
      satfrac, test, test10)
!
! soil saturation scaling factor to account for air space
! transfer saturation fraction to local variable 
! to allow data assimilation of satfrac without stomping on it
! bottom stop saturation since retrieval does not work for low saturation
   if(flg_fsat) then 
!
! scale saturation using Degesse et al. [2010] curve fit for 11.1% clay
! curve fit has artifact peak at sat_max, which is less than 1.0
     sat_max=0.9629210       
     if(satfrac>=sat_max) then
       fsat=1.
     else
       fsat=-19.63107*satfrac**2.+37.8063409*satfrac-17.2022565
!
! bottom stop fsat to prevent unrealistically large ALT near sat_zero
       fsat=max(fsat, 0.1)
     endif
!
! set lower limit on saturation based on zero deformation from Degesse et al. [2010] curve
     sat_zero=0.7372225
     sat_eff=max(satfrac,sat_zero)
!
! scale saturation to account for expansion into pore air space
     sat_eff=sat_eff*fsat
   else
!
! set lower limit of 0.2 on saturation (arbitrarily chosen)
     sat_zero=0.2
     sat_eff=max(satfrac,sat_zero)
   endif
!
! make sure saturation does not exceed 1.0
    sat_eff=min(sat_eff,1.0)
!
!-----------------------------------------------------------------------
! alt_wat: ALT for pure water column
!-----------------------------------------------------------------------
    alt_wat=rho_ice/(rho_water-rho_ice)*deform
!
!-----------------------------------------------------------------------
! alt_min: ALT for pure mineral soil
!-----------------------------------------------------------------------
    alt_min=deform/sat_eff/poros_min*rho_ice/(rho_water-rho_ice)
!
!-----------------------------------------------------------------------
! alt_org_exp: exponential decrease in organic matter
!-----------------------------------------------------------------------
! numerical integration to calculate ALT for mix of organic and mineral soil 
    depth=0.
    def_test=0.
    n_iter=0
    do iter=1,10000
!
! count total number iterations
      n_iter=n_iter+1
!
! organic soil fraction
      orgfrac=kRoot*mass_org/(1.-exp(-kRoot*root_depth))*exp(-kRoot*(depth+0.5*dz))/rho_om_max
      if(orgfrac>1.) orgfrac=1.
      if(orgfrac<0.) orgfrac=0.
      if(depth>root_depth) orgfrac=0.
!
! soil porosity
      poros=(1.-orgfrac)*poros_min+orgfrac*poros_om
      if(poros>1.) poros=1.
!
! change in deformation
      del_def=poros*sat_eff*(rho_water-rho_ice)/rho_ice*dz*100.
!
! check on current deformation and depth
      if(def_test+del_def>=deform) then ! deformation too big
        dz=0.5*dz
      else ! continue with integration
        def_test=def_test+del_def
        depth=depth+dz
      endif
!
! convergence test (deformation within 0.0001 cm)
! tests indicate this is the best compromise between 
! accuracy of estimated ALT and the number of iterations
      if(abs(def_test-deform)<.0001) exit
    enddo
!
! save ALT for mixed soil, converting from meters to cm
    alt_org_exp=depth*100.
!
!-----------------------------------------------------------------------
! alt_org_lay: organic soil layer
!-----------------------------------------------------------------------
! organic soil fraction
    orgfrac=mass_org/org_depth/rho_om_max
    orgfrac=max(orgfrac,0.)
    orgfrac=min(orgfrac,1.)
!
! soil porosity
    poros=(1.-orgfrac)*poros_min+orgfrac*poros_om
    if(poros>1.) poros=1.
!
! max deformation due to organic layer
    def_max=org_depth*sat_eff*poros*(rho_water-rho_ice)/rho_ice*100.
    if(deform<=def_max) then
      alt_org_lay=deform/sat_eff/poros*rho_ice/(rho_water-rho_ice)
    else
      alt_org_lay=def_max/sat_eff/poros*rho_ice/(rho_water-rho_ice)
      alt_org_lay=alt_org_lay+(deform-def_max)/sat_eff/poros_min*rho_ice/(rho_water-rho_ice)
    endif
!
! assign test variables
    test=alt_org_exp
    test10(1)=alt_wat
    test10(2)=alt_org_exp
    test10(3)=alt_org_lay
    test10(4)=alt_min

    RETURN
    END
!
