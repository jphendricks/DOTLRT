! from scan_read_inputs
!
!
!--------------------------------------------------------------------------
! Read input data and tables
!--------------------------------------------------------------------------
! read file containing morphology characteristics for each biome type
      k=0
      do i=1,numfiles
        if(files(i)%type=='morph_tab') k=i
      enddo
      if(k==0) then
        print*, 'Error: no morph_tab specified'
        stop
      endif
      Filename%MorphTab=files(k)%path
      allocate(MorphTab(nv))
      call ReadMorphTable
!
! read look-up table of biome dependant variables
      k=0
      do i=1,numfiles
        if(files(i)%type=='biome_tab') k=i
      enddo
      if(k==0) then
        print*, 'Error: no biome_tab specified'
        stop
      endif
      Filename%BioTab=files(k)%path
      allocate(BioTab(nv))
      call ReadBioTable
!
! read look-up table of soil dependant variables
      k=0
      do i=1,numfiles
        if(files(i)%type=='soil_tab') k=i
      enddo
      if(k==0) then
        print*, 'Error: no soil_tab specified'
        stop
      endif
      Filename%SoilTab=files(k)%path
      allocate(SoilTab(ns))
      call ReadSoilTable
!
! Read aerodynamic parameter interpolation tables
      If (scantype=='rasite') then
        k=0
        do i=1,numfiles
          if(files(i)%type=='aero_tab') k=i
        enddo
        if(k==0) then
          print*, 'Error: no aero_tab specified'
          stop
        endif
        open(unit=7,file=trim(files(k)%path),form='unformatted')
        read (7) LAIgrid
        read (7) fvcovergrid
        read (7) NewAeroTab
        Close (7)
      endif
!
! read snow class table
      k=0
      do i=1,numfiles
        if(files(i)%type=='snow_tab') k=i
      enddo
      if(k==0) then
        print*, 'Error: no snow_tab specified'
        stop
      endif
      Filename%snowtab=files(k)%path
      allocate(snowtab(nsnowclass))
      call readsnowtable


$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! From scan_control
!
! airmoss forward model
    If (scantype=='airmoss') then
      frequency=430.d6
      LAMBDA0 = light_double/frequency
      K0 = 2.D0*PI_double/LAMBDA0
      len_lay = 3.D0/K0
      std_lay(1) = 0.3D0/K0
      std_lay(2:nlay) = 0.1D0/K0
      profile=.true.
      if(profile) then
        do ilay=1, nlay-1
          call dielectric_constant_exp(DLayer(ilay), DLayer(ilay)+thickness,sat_top, watertab,&
               sand, clay,m_om, di_org_exp, test, test10)
          diel(ilay)=cmplx(di_org_exp, ltan*di_org_exp)
        enddo
        diel(nlay)=diel(nlay-1)
      endif
!      print*, sat_top, watertab
      
      call airmoss_fwd_model(nlay,frequency,thetainc,len_lay,std_lay,dlayer,diel,sigma_vv,sigma_hh, test, test10)
      sigma_vv = 10.d0*log10(sigma_vv)
      sigma_hh = 10.d0*log10(sigma_hh)
      test10(1)=sigma_vv
      test10(2)=sigma_hh
!      test=sigma_vv
      cost_hh=(sigma_hh_obs-sigma_hh)**2
      cost_vv=(sigma_vv_obs-sigma_vv)**2
      cost=cost_hh+cost_vv
!      test10(1)=real(diel(1))
!      test10(2)=real(diel(nlay))
!      test10(1)=cost_vv
!      test10(2)=cost_hh
!      test10(3)=cost
      test=cost
    endif
!
! photosynthesis freeze factor
    if(scantype=='frz_fac') call phosib_temp_factor(bionum)
!
! leaf growth and decay scaling factors
    if(scantype=='leaf_grow') call leaf_temp_factor(bionum)
!
! SiB 2-stream radiative transfer model
    If (scantype=='Rada2') call rada2(sib,sib_loc,test,test10)
!
! sibdrv diffuse/direct parameterization
    If (scantype=='Raddrv') call raddrv(1,vdcsib1,cosz, radvbc_sib,radvdc_sib,radnbc_sib,radndc_sib, test,test10)
!
! SiB3 phosib
    If (scantype=='phosib') then
!
! split radiation into direct and diffuse
      call raddrv(1,sib%drvr%sw_dwn,cosz,sib%drvr%radvbc, sib%drvr%radvdc,sib%drvr%radnbc,sib%drvr%radndc,test,test10)
!
! set some constants
      sib(1)%diag%cas_cap_co2  = max(4.0_dbl_kind,sib(1)%param%z2)
      sib(1)%diag%ra=1.
      sib(1)%diag%rb=1.
      sib(1)%diag%resp_grnd=1.e-6
      sib(1)%param%effcon(:)=effcon
      sib(1)%param%vmax0(:)=vmax0
!
! call phosib
      call phosib(sib,sib_loc)
      test10(1)=sib(1)%diag%assim(6)*1.e6
      test=sib(1)%diag%assim(6)*1.e6
!
! call phosib
      sib(1)%param%effcon(:)=effcon*0.5
      call phosib(sib,sib_loc)
      test10(2)=sib(1)%diag%assim(6)*1.e6
    endif
! 
! methane subroutine
    If(scantype=='methane') then
      call methane_setup(sib)
      call methane_core(sib)
      test10(1)=sib(1)%diag%testvar1
      test10(2)=sib(1)%prog%td(1)
      test10(3)=3.
      test=sib(1)%diag%testvar1
    endif
!
! prognostic fpar and leaf growth scaling factor
    if(scantype=='prog_fpar') then
      sib(1)%casa%lai=lai
      sib(1)%param%lai_max=MorphTab(bionum)%LAImax
      call prognostic_fpar(sib,generic)
      test10(1)=sib(1)%casa%fpar
      test10(2)=sib(1)%casa%fpar_slp
      test10(3)=sib(1)%diag%testvar1
      test10(4)=sib(1)%diag%testvar2
      test10(5)=sib(1)%diag%testvar3
      test=sib(1)%diag%testvar3
    endif
!
! SiB aerodynamics
    If (scantype=='SiBx') then
      call aero (LAI,  &
       fVcover, &
       BioTab(bionum)%chil, &
       BioTab(bionum)%z2, &
       BioTab(bionum)%z1, &
       MorphTab(bionum)%zc, &
       MorphTab(bionum)%LWidth, &
       MorphTab(bionum)%LLength, &
       NewAeroVar%zo, &
       NewAeroVar%zp_disp, &
       NewAeroVar%RbC, &
       NewAeroVar%RdC, &
       NewAeroVar%G2, &
       NewAeroVar%G3, &
       NewAeroVar%CORB1, &
       NewAeroVar%CORB2, &
       NewAeroVar%HA, &
       sigma, &
       test)
!
      if(scantype=='sibx_table') NewAeroTab(bionum,ix,iy)=NewAeroVar      
      test10(1)=NewAeroVar%zo
      test=NewAeroVar%zo
    endif
!
! Aerodynamic resistances
    If(scantype=='rasite') then
      Call NewAeroInterpolate ( &
      LAI, &
      fVCover,  &
      LAIgrid, &
      fvcovergrid, &
      NewAeroTab(bionum,:,:), &
      NewAeroVar)
!
      call Resistance(sib(1)%drvr%tm,Ta,Tc,Tgs,rhoair,um,Hflux,LAI, &
      NewAeroVar%zo, &
      NewAeroVar%zp_disp, &
      NewAeroVar%G3, &
      NewAeroVar%CORB2, &
      NewAeroVar%Rbc, &
      NewAeroVar%Rdc, &
      BioTab(bionum)%z2, &
      zwind,zmet,u2,Ra,Rb,Rd,drag,ustar,test,test10) 
    endif
!
! Vertical root distribution
    If (scantype=='root') then
      call RootFraction(len,nsoil,bionum,xroot,test,test10)
      kroot(1)=3.9
      kroot(2)=3.9
      kroot(3)=2.0
      kroot(4)=5.5
      kroot(5)=5.5
      kroot(6)=2.0
      kroot(7)=5.5
      kroot(8)=2.0
      kroot(9)=2.0
      kroot(10)=5.5
      kroot(11)=2.0
      kroot(12)=5.5
      test10(1)=(1.-exp(-kRoot(1)*generic))/kRoot(1)
      test10(2)=(1.-exp(-kRoot(3)*generic))/kRoot(3)
      test10(3)=(1.-exp(-kRoot(7)*generic))/kRoot(7)

      test10(1)=(1.-(1.-exp(-kRoot(1)*generic))/(1.-exp(-kRoot(1)*100.)))*100.
      test10(2)=(1.-(1.-exp(-kRoot(3)*generic))/(1.-exp(-kRoot(3)*100.)))*100.
      test10(3)=(1.-(1.-exp(-kRoot(7)*generic))/(1.-exp(-kRoot(7)*100.)))*100.
    endif
!
! SiB3 temperature respiration scaling factor
    If (scantype=='resp_scale') call respiration_scale_factors(bionum)
!
! liquid water fraction as function of temperature (test Code)
    If (scantype=='water_frac') then
      call liquid_water_fraction_test
!
! Assign variables calculated by SiBCASA
      sib(1)%casa%org_frac(:)=f_om
      sib(1)%prog%td(:)=td(1)
!
! subroutine that will be called in SiBCASA
      sib(1)%param%poros= 0.489-0.00126*Sand
      sib(1)%param%hfrzi=273.25

      sib(1)%param%sandfrac=sand
      sib(1)%param%clayfrac=clay

      sib(1)%param%tstar_clay=0.01
      sib(1)%param%tstar_silt=0.01
      sib(1)%param%tstar_sand=0.01
      sib(1)%param%tstar_org =0.01

      sib(1)%param%beta_clay=-0.3
      sib(1)%param%beta_silt=-0.5
      sib(1)%param%beta_sand=-0.9
      sib(1)%param%beta_org =-1.

      call liquid_water_fraction(sib,test,test10)
    endif
!
! SiB3 temperature canopy respiration scaling factor
    If (scantype=='resp_can') then
      wssp     = 0.2
      pot_fc   = -15.0   ! field capacity (J/kg)
      pot_wp   = -1500.0 ! wilt point (J/kg)
      phsat    = -10.*10**(1.88-0.0131*Sand)/1000.
      poros    = 0.489-0.00126*Sand
      bee      = 2.91+0.159*clay      
      fieldcap = poros*((pot_fc/9.8)/phsat)**(-1.0/bee)
      wilt     = poros*((pot_wp/9.8)/phsat)**(-1.0/bee)
      paw      = max((WWW(1)-wilt),0.0)
      paw_max  = (fieldcap - wilt)
      pawfrac  = paw/paw_max
      pawfrac  = MAX(0.0, pawfrac)
      pawfrac  = MIN(pawfrac, 1.0)
      rstfac(2)= ((1+wssp) * pawfrac)/(wssp + pawfrac)
!      rstfac(2)= MAX(rstfac(2), 0.1_dbl_kind)
!
! canopy scaling factor
      scatp=sqrt(1.-BioTab(bionum)%LTran(1,1)-BioTab(bionum)%Lref(1,1))
      aparkk=fpar/scatp/gmudmu
!
! canopy autotrophic respiration, with patch to prevent underflow if temp is too cool...
      var1 = 0.1 * ( tc - BioTab(bionum)%Trop )
      var2 = BioTab(bionum)%respcp * BioTab(bionum)%vmax0 * rstfac(2)
      if(tc >= BioTab(bionum)%trdm )then
        respc = var2 * 2.0**var1/( 1. + EXP( BioTab(bionum)%trda*(tc-BioTab(bionum)%trdm )))
      else
        respc = var2 * 2.0**var1
      endif
      respc=respc*aparkk
      test10(1)=respc*1.e6
      test10(2)=aparkk
      test=respc*1.e6
    endif
!
! CASA coefficients
    If (scantype=='casa') then
!      test10(1)=0.25+0.75*sand  ! oxygen
!      test10(2)=1.-0.19277*sand  ! esoilmicslow
!      test10(1)=1.-3.*clay  ! f_slow_arm
!      if(test10(1)<0.) test10(1)=0.
!      test10(2)=1.+10.6667*clay  ! f_soilmic_arm
!
! wood growth temperature response
      scale_t=1./(1.+exp(BioTab(bionum)%slti*(BioTab(bionum)%hlti-tc))) /(1.+EXP(BioTab(bionum)%shti*(tc-BioTab(bionum)%hhti)))
!
! wood growth temperature response with frost inhibition
      scale_t=1./(1.+exp(BioTab(bionum)%slti*(BioTab(bionum)%hlti-tc)))/(1.+EXP(BioTab(bionum)%shti*(tc-BioTab(bionum)%hhti)))
      scale_f=1./(1.+EXP(1.3*(278.-tc)))
!
! wood growth moisture response function
      scale_f=1./(1.+exp(50.*(.4-WWW(1))))
      test10(1)=scale_f
      test10(2)=scale_f
    endif
!
! CASA sapwood respiration
    If (scantype=='sapwood') then
!
! constants
      leafsap_ratio=0.5
      cn_ratio=124.  ! assuming 124:1 cn ratio
      e_store2wood=0.5
      k_sapwood=1./365./24./3600.
!
! wood growth temperature response with frost inhibition
      scale_t=1./(1.+exp(BioTab(bionum)%slti*(BioTab(bionum)%hlti-tc))) &
        /(1.+EXP(BioTab(bionum)%shti*(tc-BioTab(bionum)%hhti)))
      scale_f=1./(1.+EXP(1.3*(278.-tc)))
      scale_t=scale_t*scale_f
!
! sapwood
      sapwood=lai*rho_wood*BioTab(bionum)%z2*100./leafsap_ratio/mwc
!      if(sapwood>wood) sapwood=wood
      f_sapwood=sapwood/wood
      
      resp_sap=(1.-e_store2wood)*f_sapwood*wood*k_sapwood*scale_t/cn_ratio/e_store2wood
      resp_wood=(1.-e_store2wood)*1.1574e-8*scale_t*80.
    
      test10(1)=sapwood
      test10(2)=sapwood/cn_ratio*scale_t*k_sapwood
      test10(3)=resp_wood*1.e6*.5
      test=resp_sap*1.e6
    endif
!
! assim soil moisture scaling factor
    If(scantype=='gpp_moist') then
!
! old version
      phsat=-10.*10**(1.88-0.0131*Sand)/1000.
      bee=2.91+0.159*Clay
      rstfac(2)=1./(1.+EXP(0.02*(BioTab(bionum)%phi_half-phsat*WWW(2)**(-Bee))))
      rstfac(2)=MAX(0.0001, rstfac(2))
      rstfac(2)=MIN(1., rstfac(2))
      temp=rstfac(2)
!
! new plant available water version
      wssp     = 0.2
      pot_fc   = -15.0   ! field capacity (J/kg)
      pot_wp   = -1500.0 ! wilt point (J/kg)
      phsat    = -10.*10**(1.88-0.0131*Sand)/1000.
      poros    = 0.489-0.00126*Sand
      bee      = 2.91+0.159*clay
      fieldcap = poros*((pot_fc/9.8)/phsat)**(-1.0/bee)
      wilt     = poros*((pot_wp/9.8)/phsat)**(-1.0/bee)
      paw      = max((WWW(1)-wilt),0.0)
      paw_max  = (fieldcap - wilt)
      pawfrac  = paw/paw_max
      pawfrac  = MAX(0.0, pawfrac)
      pawfrac  = MIN(pawfrac, 1.0)
      rstfac(2)= ((1+wssp) * pawfrac)/(wssp + pawfrac)
      rstfac(2)= MAX(rstfac(2), 0.1)

      test10(1)=rstfac(2)
      test10(2)=temp
      test=rstfac(2)
    endif
!
! drought stress factor
    If(scantype=='drought') then
      test10(1)=((1+0.0) * WWW(1))/ (0.0 + WWW(1))
      test10(2)=((1+0.05) * WWW(1))/ (0.05 + WWW(1))
      test10(3)=((1+0.1) * WWW(1))/ (0.1 + WWW(1))
      test10(4)=((1+0.15) * WWW(1))/ (0.15 + WWW(1))
      test10(5)=((1+0.2) * WWW(1))/ (0.2 + WWW(1))
      test=((1+generic) * WWW(1))/ (generic + WWW(1))
    endif
!
! soil thermodynamic properties
    If (scantype=='thermo') call soil_thermo_properties(ix,iy)
!
! soil organic layer dynamics
    If (scantype=='org_lay') call soil_organic_layer(bionum)
!
! Active layer porosity profiles
    if(scantype=='ALT_por') then
!
! constants
      kroot(10)=5.5
      rho_om_max=140.
      poros_om=.9
      root_depth=0.7
      org_depth=0.18
!
! mineral porosity
      test10(1)=0.489-0.00126*Sand
!
! mixed soil organic fraction (generic is depth)
      var1=kroot(10)*m_om*exp(-kroot(10)*depth)/(1.-exp(-kroot(10)*root_depth))/rho_om_max
      if(var1>1.) var1=1.
      if(var1<0.) var1=0.
      if(depth>root_depth) var1=0.
!
! mixed soil porosity
      test10(2)=(1.-var1)*test10(1)+var1*poros_om
      if(test10(2)>1.) test10(2)=1.
!
! organic layer porosity
      if(depth<=org_depth) then
        var1=0.489-0.00126*Sand
        var2=m_om/org_depth/rho_om_max
        if(var2>1.) var1=1.
        if(var2<0.) var1=0.
        test10(3)=(1.-var2)*var1+var2*poros_om
         if(test10(3)>1.) test10(3)=1.
      else
        test10(3)=0.489-0.00126*Sand
      endif
    endif
!
! ALT with a water table
    if(scantype=='ALT_water') then
      kroot(10)=5.5
      rho_om_max=140.
      poros_om=.9
      root_depth=0.7
      org_depth=0.18

      call alt_from_deform_water_table(deform, satfrac,wat_tab,sand,kroot(10),m_om,&
                    root_depth, org_depth, rho_om_max, poros_om,var1, var2, var3, var4)
      test=var2
      test10(1)=var1
      test10(2)=var2
      test10(3)=var3
      test10(4)=var4
    endif
!
! ALT for riverbeds
    if(scantype=='ALT_river') then
      kroot(10)=5.5
      rho_om_max=140.
      poros_om=.9
      root_depth=0.7
      org_depth=0.18
      call river_alt_from_deform(deform, generic,satfrac,sand,kroot(10),m_om,&
                    root_depth, rho_om_max, poros_om,test,test10)
      ref_alt=var2
      uncert_cum=0.
    endif
!
! Depth average VWC
    if(scantype=='ALT_VWC') then
!
! constants
      kroot(10)=5.5
      rho_om_max=140.
      poros_om=.9
      root_depth=0.7
      org_depth=0.18

      call volumetric_water_content(depth, satfrac,sand,kroot(10),m_om,&
           root_depth, org_depth, rho_om_max, poros_om,test,test10)
    endif
!
! Depth average saturation fraction
    if(scantype=='satfrac') then
!
! constants
      kroot(10)=5.5
      rho_om_max=140.
      poros_om=.9
      root_depth=0.7
      org_depth=0.18

      call saturation_fraction(depth_top, depth_bot, vwc,&
           sand, kroot(10), m_om, root_depth, org_depth, rho_om_max, poros_om,&
           junkvar,test,test10)
    endif
!
! joint ReSALT/Airmoss retrieval: ALT from deformation and vwc
   if(scantype=='ALT_joint') then
      kroot(10)=5.5
      rho_om_max=140.
      poros_om=.9
      root_depth=0.7
      org_depth=0.18

      call alt_from_deform_joint(depth_top,depth_bot,vwc,deform,&
           sand, kroot(10), m_om, root_depth, org_depth, rho_om_max, poros_om,&
           junkvar, var1, var2, var3, var4, test, test10)
   endif
!
! Active layer/surface deformation subroutine
    if(scantype=='ALT_def') then
!
! constants
      kroot(10)=5.5
      rho_om_max=140.
      poros_om=.9
      root_depth=0.7
      org_depth=0.18

      call alt_from_deform_v2(deform, satfrac,sand,kroot(10),m_om,&
        root_depth, org_depth, rho_om_max, poros_om,var1, var2, var3, var4,test,test10)
      ref_alt=var2
      uncert_cum=0.
!
! parameter uncertainties
      uncert_param(1)=0.5   ! deformation
      uncert_param(2)=0.1   ! satfrac
      uncert_param(3)=10.   ! m_om
      uncert_param(4)=.05   ! poros_om
      uncert_param(5)=10.   ! rho_om_max
      uncert_param(6)=14.7  ! fsand
      uncert_param(7)=0.5   ! kroot
      uncert_param(8)=0.1   ! droot
!
! pick uncertainty to save
!  n=0 nothing thank you
!  n=1 absolute uncertainty
!  n=2 uncertainty relative to ALT
!  n=3 cummulative uncertainty (total uncertainty all variables is test10(9))
!  n=3 relative contribution uncertainty
              n=0
              if(n>0) then
!
! 1) deformation uncertainty
            slope=ref_alt
            junkvar=deform-uncert_param(1)
            if(junkvar<0.0) junkvar=0.
            call alt_from_deform_v2(junkvar, satfrac,sand,kroot(10),m_om,&
                    root_depth, org_depth, rho_om_max, poros_om,var1, var2, var3, var4)
            slope=(slope-var2)/(deform-junkvar)
            uncert_alt(1)=slope*uncert_param(1)
            uncert_cum(1)=uncert_alt(1)*uncert_alt(1)      
!
! 2) satfrac uncertainty
            slope=ref_alt
            junkvar=satfrac-uncert_param(2)
            if(junkvar<0.0) junkvar=satfrac+uncert_param(2)
            call alt_from_deform_v2(deform, junkvar,sand,kroot(10),m_om,&
                    root_depth, org_depth, rho_om_max, poros_om,var1, var2, var3, var4)
            slope=(slope-var2)/(satfrac-junkvar)
            uncert_alt(2)=slope*uncert_param(2)
            uncert_cum(2)=uncert_cum(1)+uncert_alt(2)*uncert_alt(2)
!
! 3) organic matter uncertainty
            slope=ref_alt
            junkvar=m_om-uncert_param(3)
            if(junkvar<0.0) junkvar=m_om+uncert_param(1)
            call alt_from_deform_v2(deform, satfrac,sand,kroot(10),junkvar,&
                    root_depth, org_depth, rho_om_max, poros_om,var1, var2, var3, var4)
            slope=(slope-var2)/(m_om-junkvar)
            uncert_alt(3)=slope*uncert_param(3)
            uncert_cum(3)=uncert_cum(2)+uncert_alt(3)*uncert_alt(3)
!
! 4) porosity organic matter uncertainty
            slope=ref_alt
            junkvar=poros_om-uncert_param(4)
            if(junkvar<0.0) junkvar=poros_om+uncert_param(4)
            call alt_from_deform_v2(deform, satfrac,sand,kroot(10),m_om,&
                    root_depth, org_depth, rho_om_max, junkvar,var1, var2, var3, var4)
            slope=(slope-var2)/(poros_om-junkvar)
            uncert_alt(4)=slope*uncert_param(4)
            uncert_cum(4)=uncert_cum(3)+uncert_alt(4)*uncert_alt(4)
!
! 5) density organic matter max uncertainty
            slope=ref_alt
            junkvar=rho_om_max-uncert_param(5)
            if(junkvar<0.0) junkvar=rho_om_max+uncert_param(5)
            call alt_from_deform_v2(deform, satfrac,sand,kroot(10),m_om,&
                    root_depth, org_depth, junkvar, poros_om,var1, var2, var3, var4)
            slope=(slope-var2)/(rho_om_max-junkvar)
            uncert_alt(5)=slope*uncert_param(5)
            uncert_cum(5)=uncert_cum(4)+uncert_alt(5)*uncert_alt(5)
!
! 6) sand fraction uncertainty
            slope=ref_alt
            junkvar=sand-uncert_param(6)
            if(junkvar<0.0) junkvar=sand+uncert_param(6)
            call alt_from_deform_v2(deform, satfrac,junkvar,kroot(10),m_om,&
                    root_depth, org_depth, rho_om_max, poros_om,var1, var2, var3, var4)
            slope=(slope-var2)/(sand-junkvar)
            uncert_alt(6)=slope*uncert_param(6)
            uncert_cum(6)=uncert_cum(5)+uncert_alt(6)*uncert_alt(6)
!
! 7) kroot uncertainty
            slope=ref_alt
            junkvar=kroot(10)-uncert_param(7)
            if(junkvar<0.0) junkvar=kroot(10)+uncert_param(7)
            call alt_from_deform_v2(deform, satfrac,sand,junkvar,m_om,&
                    root_depth, org_depth, rho_om_max, poros_om,var1, var2, var3, var4)
            slope=(slope-var2)/(kroot(10)-junkvar)
            uncert_alt(7)=slope*uncert_param(7)
            uncert_cum(7)=uncert_cum(6)+uncert_alt(7)*uncert_alt(7)
!
! 8) rooting depth uncertainty
            slope=ref_alt
            junkvar=root_depth-uncert_param(8)
            if(junkvar<0.0) junkvar=root_depth+uncert_param(8)
            call alt_from_deform_v2(deform, satfrac,sand,kroot(10),m_om,&
                    junkvar, org_depth, rho_om_max, poros_om,var1, var2, var3, var4)
            slope=(slope-var2)/(root_depth-junkvar)
            uncert_alt(8)=slope*uncert_param(8)
            uncert_cum(8)=uncert_cum(7)+uncert_alt(8)*uncert_alt(8)
!
! absolute uncertainty (total uncertainty all variables is test10(10))
            if(n==1) then
              test10(1)=ref_alt
              do n=1,8
                test10(n+1)=abs(uncert_alt(n))
                enddo
              test10(10)=sqrt(uncert_cum(8))
            endif
!
! uncertainty relative to ALT
            if(n==2) then
              test10(1)=ref_alt
              if(ref_alt>0.) then
                do n=1,8
                  test10(n+1)=abs(uncert_alt(n))/ref_alt*100.
                  test10(n+1)=min(test10(n+1),100.)
                enddo
                test10(10)=sqrt(uncert_cum(8))/ref_alt*100.
                test10(10)=min(test10(10),100.)
              else
                test10(10)=100.
              endif
            endif
!
! cummulative uncertainty (total uncertainty all variables is test10(9))
            if(n==3) then
              test10(1)=ref_alt
              do n=1,8
                test10(n+1)=sqrt(uncert_cum(n))
              enddo
              test10(10)=sqrt(uncert_cum(8))
            endif
!
! relative contribution uncertainty
            if(n==4) then
              test10(1)=ref_alt
              test10(2)=sqrt(uncert_cum(1))/sqrt(uncert_cum(8))*100.
              do n=2,8
                test10(n+1)=(sqrt(uncert_cum(n))-sqrt(uncert_cum(n-1)))/sqrt(uncert_cum(8))*100.
              enddo
              test10(10)=100.
            endif
          endif
    endif
!
! permafrost thaw
    if(scantype=='perma_thaw') call permacarb(generic, test,test10)
!
! snow thermodynamic properties
    If (scantype=='snow') then
            call snow_test_setup
    
              ! compaction: thermal aging (desctructive metamorphism)
              comp_dm= -sno_c3*exp(-sno_c4*(tice - snow_temp))
              if(sib(1)%prog%snow_bulk > sno_dm) comp_dm = comp_dm*exp(-46.e-3*(sib(1)%prog%snow_bulk-sno_dm))
              if(m_liq > 0.0) comp_dm = comp_dm*sno_c5

              comp_dm= -sno_c3*exp(-sno_c4*(tice - snow_temp))
              comp_dm= comp_dm/(1+exp(generic*(sib(1)%prog%snow_bulk-100.)))
              if(m_liq > 0.0) comp_dm = comp_dm*sno_c5

              ! compaction: weight
              comp_burden = -(m_liq+m_ice)*exp(-0.08*(tice - snow_temp)-sno_c2*sib(1)%prog%snow_bulk)/sno_eta0
              comp_burden=comp_burden*86400.      

              ! snowfall density as a function of air temperature (original)
              if(sib(1)%drvr%tm > tice + 2.0) then
                snow_den = 189.0
              elseif(sib(1)%drvr%tm > tice - 15.0) then
                snow_den = 50.0 + 1.7 * (sib(1)%drvr%tm - tice + 15.0)**1.5
              else
                snow_den = 50.0
              endif
              
              ! snowfall density as a function of air temperature (original with discontinuity removed)
              if(sib(1)%drvr%tm > tice + 2.0) then 
                snow_den = 169.1577
              elseif(sib(1)%drvr%tm > tice - 15.0) then
                snow_den = 50.0 + 1.7 * (sib(1)%drvr%tm - tice + 15.0)**1.5
              else
                snow_den = 50.0
              endif
              
              ! snowfall density as a function of air temperature (final with 2.5 C threshold)
              snow_den = 50.0
              if(sib(1)%drvr%tm > tice - 15.0) snow_den = snow_den + 1.7 * (sib(1)%drvr%tm - tice + 15.0)**1.5 
              delta_den=(400.-snow_den)/(1+exp(1.5*(3.5-wind)))
              delta_den=(400.-50.)*(wind-2.)/3.
              delta_den=max(delta_den,0.)
              snow_den=snow_den+delta_den
              snow_den=min(snow_den,400.)
              if(sib(1)%drvr%tm > tice + 2.5) snow_den = 0.
              test10(1)=snow_den
              test=snow_den
              snow_den=50+140./(1+exp(.3*(270.-sib(1)%drvr%tm)))

! subdivide snow
              !  call test_subdivide_snow(sib)
              !  call test_snow_density(sib)
              snow_den=0.
              do n=abs(sib(1)%prog%nsl),1,-1
                if(snow_z>=-sib(1)%prog%z_bot(n+sib(1)%prog%nsl)) snow_den=sib(1)%prog%snow_den(n)
              enddo
              if(snow_z>sib(1)%prog%snow_depth) snow_den=0.
!
! combine snow
         !     sib%prog%nsl=-nsnow
         !     sib(1)%prog%snow_depth=sib(1)%prog%snow_mass/sib(1)%prog%snow_den
         !     sib(1)%prog%dz(1+sib%prog%nsl:0)=sib(1)%prog%snow_depth/real(nsnow)
         !     sib(1)%prog%www_ice(1+sib%prog%nsl:0)=sib(1)%prog%snow_mass/real(nsnow)
         !       call test_combine_snow(sib)
            !
            ! snow table
            n=snowclass
            !
            ! bottom layer fraction
            d_slope=10./(snowtab(n)%d_max-snowtab(n)%d_min)
            d_half=0.5*(snowtab(n)%d_max+snowtab(n)%d_min)
            f_bot=snowtab(n)%f_bot_max/(1+exp(d_slope*(d_half-sib(1)%prog%snow_depth)))
          !  if(f_bot<0.01) f_bot=0.
            d_bot=f_bot*sib(1)%prog%snow_depth
            !
            ! delta density function
            delta_den=abs(snowtab(n)%den_ref_top-snowtab(n)%den_ref_bot)*f_bot/snowtab(n)%f_bot_max
            s_delta_den=exp(-(sib(1)%prog%snow_bulk-snowtab(n)%den_bulk_obs)**2./snowtab(n)%den_bulk_std**2.)
            delta_den=delta_den*s_delta_den
           ! if(delta_den<.001) delta_den=0.
            !
            ! standard snow density profile
            if(n==1.or.n==2.or.n==6) then
              snow_den_bot=sib(1)%prog%snow_bulk-(1-f_bot)*delta_den
              snow_den_top=sib(1)%prog%snow_bulk+f_bot*delta_den
                if(zbot>d_bot) then
                snow_den=snow_den_top
              elseif(ztop<d_bot) then
                snow_den=snow_den_bot
              else
              fac_bot=(d_bot-zbot)/lay_dz
              fac_top=1.-fac_bot
                snow_den=fac_bot*snow_den_bot+fac_top*snow_den_top
              endif
              elseif(n==3.or.n==4.or.n==5) then ! linear density top layer, constant density bottom layer
                snow_den_bot=sib(1)%prog%snow_bulk+(1.-f_bot)*delta_den*.5
              den_slope=-delta_den/(sib(1)%prog%snow_depth-d_bot)
                if(zbot>d_bot) then
                z_midsno=(zbot+ztop)*.5
                snow_den=den_slope*(z_midsno-d_bot)+snow_den_bot
              elseif(ztop<d_bot) then
                snow_den=snow_den_bot
              else
                z_midsno=(d_bot+ztop)*.5
              fac_bot=(d_bot-zbot)/lay_dz
              fac_top=1.-fac_bot
                snow_den=fac_bot*snow_den_bot+fac_top*(den_slope*(z_midsno-d_bot)+snow_den_bot)
              endif
              endif
              lay_mass=snow_den*lay_dz
            if(snow_z<Stop(57)) tot_mass=tot_mass+lay_mass
            if(snow_z>sib(1)%prog%snow_depth) snow_den=0.

            ! snow thermal conductivity
            snow_tcon = tcon_air + (7.75E-5*sib(1)%prog%snow_bulk + &
              1.105E-6*sib(1)%prog%snow_bulk*sib(1)%prog%snow_bulk)*(tcon_ice-tcon_air) 
            !
            ! ice fraction of precipitation (old)
            if(sib(1)%drvr%tm <= tice) then
                fliq = 0.0
              elseif(sib(1)%drvr%tm < tice + 2.0) then
                fliq = 0.2*sib(1)%drvr%tm-54.64
              else
                fliq = 0.40
              endif
            if(sib(1)%drvr%tm>=tice+2.5) fliq=1.
              fliq = max(0.0,fliq)
            !
            ! ice fraction of precipitation (new)
            fliq = 0.4*sib(1)%drvr%tm-109.26
            if(sib(1)%drvr%tm>=tice+2.5) fliq=1.
              fliq = max(0.0,fliq)
            !
            ! water loss due to drainage (current)
            lay_dz=0.02
            sib(1)%prog%www_ice=m_ice
            sib(1)%prog%www_liq=m_liq
            
            sib(1)%prog%vol_ice(1)=sib(1)%prog%www_ice(1)/(lay_dz*rho_ice)
            sib(1)%prog%vol_ice(1)=min(sib(1)%prog%vol_ice(1),1._dbl_kind)
            sib(1)%prog%vol_ice=sib(1)%prog%vol_ice(1)
            
              sib(1)%diag%eff_poros = 1.0 - sib(1)%prog%vol_ice(1)
            
              sib(1)%prog%vol_liq(1)=m_liq/lay_dz/rho_water
            if(sib(1)%prog%vol_liq(1)>sib(1)%diag%eff_poros(1)) then
            v_run=sib(1)%prog%vol_liq(1)-sib(1)%diag%eff_poros(1)
              sib(1)%prog%vol_liq(1)=sib(1)%diag%eff_poros(1)
            endif
            sib(1)%prog%vol_liq=sib(1)%prog%vol_liq(1)
            runoff=v_run*1000.0*lay_dz
            
            v_air=sib(1)%diag%eff_poros(1)-sib(1)%prog%vol_liq(1)
            
              if(sib(1)%diag%eff_poros(1)<0.05 .or. sib(1)%diag%eff_poros(2)< 0.05) then
                 qout = 0.0
              else
                 qout = max(0.0_dbl_kind,(sib(1)%prog%vol_liq(1)-ssi*sib(1)%diag%eff_poros(1)))
                 qout = min(qout,sib(1)%diag%eff_poros(2)-sib(1)%prog%vol_liq(2))
              endif
              qout=qout*1000.0*lay_dz
            qout=min(qout, sib(1)%prog%www_liq(1))

            lay_dz=0.02
              sib(1)%prog%vol_ice(1) = m_ice/(lay_dz*rho_ice)
              sib(1)%diag%eff_poros(1) = 1.0 - sib(1)%prog%vol_ice(1)
              sib(1)%prog%vol_liq(1)=m_liq/(lay_dz*rho_water)
            
            qout=max(0.0_dbl_kind,(sib(1)%prog%vol_liq(1)- ssi*sib(1)%diag%eff_poros(1))*lay_dz)
              qout=min(qout,(1.0_dbl_kind-sib(1)%prog%vol_ice(1)- sib(1)%prog%vol_liq(1))*lay_dz)
            
            ! snow hydraulic conductivity
            lay_dz=0.02
            sib(1)%prog%vol_ice(1) = m_ice/(lay_dz*rho_ice)
            if(sib(1)%prog%vol_ice(1)>1.) then
              m_ice=lay_dz*rho_ice
            endif
              sib(1)%diag%eff_poros(1) = 1.0 - sib(1)%prog%vol_ice(1)

            sib(1)%prog%vol_liq(1) = min(sib(1)%diag%eff_poros(1),m_liq/(lay_dz*rho_water))
            if(sib(1)%prog%vol_liq(1)>=sib(1)%diag%eff_poros(1)) then
              m_liq=sib(1)%diag%eff_poros(1)*lay_dz*rho_water
            endif
           ! snow_den=(m_liq+m_ice)/lay_dz
            snow_hcon=8.e-23*10.**(15.*generic+0.029*snow_den)
            snow_dif=4.167e-11*10.**(10.*generic+0.008*snow_den)
            thetaslope=generic*sib(1)%prog%vol_liq(1)/lay_dz
            qout=(snow_dif*thetaslope+snow_hcon)*600./lay_dz
            qout=min(qout, m_liq)

              ! runoff fraction of liquid water
            frun=1./(1+exp(.05*(500.-snow_den)))
            frun=0.
            if(snow_den>400.) frun=0.005*snow_den-2.
            frun=min(frun, 1.)
    endif
!
! leaf area index
    If (scantype=='lai') call TestLAI (fPAR,fPARm,fPARmax,fVCover, MorphTab(bionum)%stems,MorphTab(bionum)%LAImax,&
              Green,LAI, test, test10)
!
! Mixed Layer CO2
    If(scantype=='co2_ppm') then ! calculate global average co2 concentration
      if(tau<2008.) then ! use curve fit of global observed co2
        var1=280.+0.27*exp(0.019325329*(tau-1700.))
      else ! use IPCC A1B projection to stabilize at 688 ppm by 2100
        var1=3.3887*(tau-2008.)+384.85
      endif
      var1=min(var1,688.536)
      
      call calculate_co2(tau,lat,var2,var3, test,test10)
      test10(1)=var2-var3
      test10(2)=var2
      test=var2
    endif
!
! calculating LAI from leaf biomass
    If(scantype=='pool_lai') call LAI_from_leaf_pool()
!
! scaling of specific leaf area
    If(scantype=='sla') call sla_scaling()
!
! leaf to canopy scaling factor
    If(scantype=='park') then
      scatp=sqrt(1.-BioTab(bionum)%LTran(1,1)-BioTab(bionum)%Lref(1,1))
      aparkk=.95/scatp/gmudmu
      test10(1)=aparkk
    endif
!
! respiration
    If (scantype=='Resp') then
              call respsib( &
              len, nsib, nsoil, Soiltab(bionum)%wopt,Soiltab(bionum)%skew, &
              Soiltab(bionum)%RespSat,tgs, td, www,forcerestore, &
              respfactor, respg, soilscale, Moist, soilq10, &
              Respc, scaleRsp, ScaleH2o, RconWat,RconTem, &
              test, test10)
    endif
!
! Kevins modified version of phosib
    If (scantype=='kphosib') then
              call kevinphosib (test,test10,Tc, Psur, &
              Biotab(bionum)%vmax0, BioTab(bionum)%Trop, &
              BioTab(bionum)%SHTI, BioTab(bionum)%HHTI, &
              BioTab(bionum)%SLTI, BioTab(bionum)%HLTI, &
              BioTab(bionum)%TRDA, BioTab(bionum)%TRDM, &
              BioTab(bionum)%Atheta, BioTab(bionum)%Btheta, &
              BioTab(bionum)%effcon, BioTab(bionum)%respcp, &
              BioTab(bionum)%LTran, BioTab(bionum)%LRef, &
              SoilTab(bionum)%Bee, SoilTab(bionum)%phisat,BioTab(bionum)%phi_half, &
              fPAR,Green,gmudmu, &
              C3,C4,  &
              po2m, pco2i, WWW, PAR, Tice)
    endif
!
! absorbed fraction of par (fpar)
    If(scantype=='fPAR') then
            var1=fPARmax
            var2=fParmin
            call testfPAR (ndvi,MorphTab(bionum)%NDVImin,MorphTab(bionum)%NDVImax,MorphTab(bionum)%SRmin,MorphTab(bionum)%SRmax,&
              var1,var2,fPAR,test,test10)
            call TestLAI (fPAR,fPAR,var1,fVCover, MorphTab(bionum)%stems,MorphTab(bionum)%LAImax,Green,LAI,test,test10)
              test10(1)=LAI

            call testfPAR (ndvi,MorphTab(bionum)%NDVImin,MorphTab(bionum)%NDVImax,MorphTab(bionum)%SRmin,MorphTab(bionum)%SRmax,&
              var1,var2,fPAR,test,test10)
              call TestLAI (fPAR,fPAR,var1,fVCover, MorphTab(bionum)%stems,3.,Green,LAI,test,test10)
              test10(2)=LAI

            call testfPAR (ndvi,MorphTab(bionum)%NDVImin,MorphTab(bionum)%NDVImax,MorphTab(bionum)%SRmin,MorphTab(bionum)%SRmax,&
              var1,var2,fPAR,test,test10)
            call TestLAI (fPAR,fPAR,var1,fVCover, MorphTab(bionum)%stems,generic,Green,LAI,test,test10)
              test10(3)=LAI
            test=lai
     
            call aparnew (LAI, &
                    Green, &
                    BioTab(bionum)%LTran, &
                    BioTab(bionum)%LRef, &
                    gmudmu, &
                    BioTab(bionum)%fVCover, &
                    fPAR, &
                    fPARmax, &
                    fPARmin)
!              test10(1)=fpar

            fvcover=test10(3)/fPARmax
            call TestLAI (fpar,fpar,fPARmax,fVCover, MorphTab(bionum)%stems,MorphTab(bionum)%LAImax, Green,LAI, test, test10)
    endif
!
! gmudmu
    If (scantype=='gmudmu') then
              call Oldgmuder (Lat, DOY, biotab(bionum)%ChiL, gmuold, &
              test, test10)
              test=gmuold
              call gmuder (Lat, DOY, biotab(bionum)%ChiL, gmudmu)
    endif

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!=======================================================================
      subroutine scan_non_std_routines
!=========================================================================
! controls all calls to all non-standard (non0generic) scans
!
! Modifications:
!  Kevin Schaefer split off from main program (5/27/13)
!--------------------------------------------------------------------------
!
    use Mapper_Variables
    use sibscan_Variables
!
    IMPLICIT NONE
!
integer i,j, k
!
!--------------------------------------------------------------------------
! scan though SiB declination calculations
!--------------------------------------------------------------------------
      If (scantype=='sun_dec') then
      print *, 'scan dec calculations using ', numpts, ' points'
!
      allocate(X(nv,numpts))
      allocate(Y(nv,numpts))
      allocate(X1(nv,numpts*numpts))
      allocate(Y1(nv,numpts*numpts))
      allocate(Z(nv,numpts,numpts))
      allocate(X3(numpts))
!
! Calculate the solar declination
      eqnx    = 80
      RECCN=1./(1.-ECCN*ECCN)**1.5
      H=2.*PI/DAYPYR
!
      Dx=H
      open(unit=44, file='dectest',form='formatted')
!
      angle=0.-real(eqnx)*pi/180.
      do j=1,numpts
!
!     Calculate the right ascension of the Earth
          T1=RECCN*(1.-ECCN*COS(angle-PERHLR))**2*H
          T2=RECCN*(1.-ECCN*COS(angle+T1*.5-PERHLR))**2*H
          T3=RECCN*(1.-ECCN*COS(angle+T2*.5-PERHLR))**2*H
          T4=RECCN*(1.-ECCN*COS(angle+T3-PERHLR))**2*H
          angle=angle+(T1+2.*(T2+T3)+T4)/6.
!
!       Calculate the sine and cosine of Solar declination
        SIND=sin(decmax*(pi/180.)) * SIN(angle)
        COSD=SQRT(1.-SIND*SIND)

        dec=asin(SIND)*180/Pi
        Y(1,j)=dec
        dec=decmax*sin(2*PI/DAYPYR*real(j-eqnx))
        Y(2,j)=dec
!        Y(1,j)=asin(SIND)*180./PI-decmax*sin(2*PI/DAYPYR*real(j-eqnx))
!        Y(1,j)=(SIND-sin(dec*pi/180.))/sin(decmax/PI*180)*100

        X(1,j)=real(j)
        X(2,j)=real(j)
!
! end biome loop
      enddo
!
      i=0
      do j=1,numpts,2
        i=i+1
        read(44,*) X(3,i), Y(3,i)
        Y(4,i)=Y(3,i)-Y(1,j)
        Y(5,i)=Y(3,i)-Y(2,j)
        X(4,i)=X(3,i)
        X(5,i)=X(3,i)
      enddo
!
      If (PlotFlag==1) then
        LabX='Day of Year'
        LabY='dec err (deg)'
        Title='solar declination'
!        call Line (i,X,Y,1.,365.,-.5,2.,DevFlag,
!     &  LabX,LabY,Title,nv,4,5)
      End If
!
      End If
!
!--------------------------------------------------------------------------
! scan 3-D stable Ustar error calculations
!--------------------------------------------------------------------------
      If (scantype=='Stable3D') then
      print *, 'scan 3-D stable Ustar error ', numpts, ' points'
!
      allocate(X(nv,numpts))
      allocate(Y(nv,numpts))
      allocate(X1(nv,numpts*numpts))
      allocate(Y1(nv,numpts*numpts))
      allocate(Z(nv,numpts,numpts))
      allocate(X3(numpts))
!
      sib(1)%drvr%tm = 301.0
      zp_disp=26
      z0=2.5
      um = 5.0
      Hflux = 300.0
      rhoair = 1.225
      zwind = 100.0
      zmet = 100.0
      vkc = 0.41
      cpair = 1004.67
      gamma = 16.0
      g = 9.8
      argz=dlog((zwind-zp_disp)/z0)
      Dy=5./real(numpts-1)
      Dx=100./real(numpts-1)
      Dz=5.
      ustar=vkc*um/dlog((zwind-zp_disp)/z0)
!
      do k=1,nv
        uest=0.001
        do j=1,numpts
          Hflux=-100.
          do i=1,numpts
            Obhukov=-uest**3*rhoair*cpair*sib(1)%drvr%tm/(vkc*g*Hflux)
            psione=-4.7*(zwind-zp_disp)/Obhukov
            psione=amax1(-4.7,psione)
            Error=um-uest*(argz-psione)/vkc
            z(k,i,j)=psione
            Hflux=HFlux+Dx
          enddo
          uest=uest+Dy
        enddo
        um=um+Dz
      enddo
!
      If (PlotFlag==1) then
        LabX='Hflux (W/m**2)'
        LabY='ustar (m/s)'
        Title='psione'
        print*, title
        Plotfile=trim(Title)
        call Contour (Z,numpts,Numpts,-100., 0.001, Dx, Dy,  &
        DevFlag, Scale, mincont,maxcont,base,NC,nv,MinBiome,MaxBiome, &
        LabX,LabY,Title, Plotfile)
      End If
!
      End If
!
!--------------------------------------------------------------------------
! scan 3-D unstable Ustar error calculations
!--------------------------------------------------------------------------
      If (scantype=='Unstable3D') then
      print *, 'scan 3-D unstable Ustar error ', numpts, ' points'
!
      allocate(X(nv,numpts))
      allocate(Y(nv,numpts))
      allocate(X1(nv,numpts*numpts))
      allocate(Y1(nv,numpts*numpts))
      allocate(Z(nv,numpts,numpts))
      allocate(X3(numpts))
!
      sib(1)%drvr%tm = 301.0
      zp_disp=26
      z0=2.5
      um = 5.0
      Hflux = 300.0
      rhoair = 1.225
      zwind = 100.0
      zmet = 100.0
      vkc = 0.41
      cpair = 1004.67
      gamma = 16.0
      g = 9.8
      argz=dlog((zwind-zp_disp)/z0)
      Dy=5./real(numpts-1)
      Dx=100./real(numpts-1)
      Dz=5.
      ustar=vkc*um/dlog((zwind-zp_disp)/z0)
!
      do k=1,nv
        uest=0.001
        do j=1,numpts
          Hflux=0.001
          do i=1,numpts
            Obhukov=-uest**3*rhoair*cpair*sib(1)%drvr%tm/(vkc*g*Hflux)
            xCon=(1.-gamma*(zwind-zp_disp)/Obhukov)**0.25
            psione=2.*alog((1.+xCon)/2.)+alog((1.+xCon*xCon)/2.) &
                -2.*atan(xCon)+pi/2.
            psione=amin1(argz*0.75,psione)
            Error=um-uest*(argz-psione)/vkc
            z(k,i,j)=argz-psione
            Hflux=HFlux+Dx
          enddo
          uest=uest+Dy
        enddo
        um=um+Dz
      enddo
!
      If (PlotFlag==1) then
        LabX='Hflux (W/m**2)'
        LabY='ustar (m/s)'
        Title='slope'
        print*, title
        Plotfile=trim(Title)
        call Contour (Z,numpts,Numpts,0.001, 0.001, Dx, Dy,  &
        DevFlag, Scale,mincont,maxcont,base,NC,nv,MinBiome,MaxBiome, &
        LabX,LabY,Title, Plotfile)
      End If
!
      End If
!
!--------------------------------------------------------------------------
! scan 2-D unstable Ustar error calculations
!--------------------------------------------------------------------------
      If (scantype=='Unstable2D') then
      print *, 'scan 2-D unstable Ustar error ', numpts, ' points'
!
      allocate(X3(numpts))
      allocate(Y(nv,numpts))
!
      sib(1)%drvr%tm = 301.0
      zp_disp=26
      z0=2.5
      Hflux = 100.
      rhoair = 1.225
      zwind = 100.0
      zmet = 100.0
      vkc = 0.41
      cpair = 1004.67
      gamma = 16.0
      g = 9.8
      argz=dlog((zwind-zp_disp)/z0)
      Dx=5./real(numpts-1)
      Dz=5.
!
! plot values for neutral ustar
      ustar=vkc*um/dlog((zwind-zp_disp)/z0)
!
      um = 0.001
      do k=1,nv
        uest=0.001
        do j=1,numpts
          Obhukov=-uest**3*rhoair*cpair*sib(1)%drvr%tm/(vkc*g*Hflux)
          xCon=(1.-gamma*(zwind-zp_disp)/Obhukov)**0.25
          psione=2.*alog((1.+xCon)/2.)+alog((1.+xCon*xCon)/2.) &
                -2.*atan(xCon)+pi/2.
          psione=amin1(argz*0.75,psione)
          Error=um-Uest*(argz-psione)/vkc
!
          y(k,j)=argz-psione
          x3(j)=uest
!
          uest=uest+Dx
        enddo
        um=um+Dz
      enddo
!
      If (PlotFlag==1) then
        Laby='Obhukov (m)'
        LabX='ustar (m/s)'
        Title='Unstable'
        title=trim(Title)
        Plotfile='Error'
        call MultiLine (Numpts,X3,Y,Numpts,X3,Y, &
        0.,5.,minval(Y),maxval(Y), &
        DevFlag,LabX,LabY,title,nv, 1, nv, Plotfile, legFlag, Legend)
      End If
!
      End If
!
!--------------------------------------------------------------------------
! scan 2-D stable Ustar error calculations
!--------------------------------------------------------------------------
      If (scantype=='Stable2D') then
      print *, 'scan 2-D table Ustar error ', numpts, ' points'
!
      allocate(X3(numpts))
      allocate(Y(nv,numpts))
!
      sib(1)%drvr%tm = 301.0
      zp_disp=26
      z0=2.5
      um = 0.001
      Hflux = -20.
      rhoair = 1.225
      zwind = 100.0
      zmet = 100.0
      vkc = 0.41
      cpair = 1004.67
      gamma = 16.0
      g = 9.8
      argz=dlog((zwind-zp_disp)/z0)
      Dx=1./real(numpts-1)
      Dz=5.
      ustar=vkc*um/dlog((zwind-zp_disp)/z0)
!
      do k=1,nv
        uest=0.001
        do j=1,numpts
          Obhukov=-uest**3*rhoair*cpair*sib(1)%drvr%tm/(vkc*g*Hflux)
          psione=-4.7*(zwind-zp_disp)/Obhukov
          psione=amax1(-4.7,psione)
          Error=um-Uest*(argz-psione)/vkc
!
          y(k,j)=psione
          x3(j)=uest
!
          uest=uest+Dx
        enddo
        um=um+Dz
      enddo
!
      If (PlotFlag==1) then
        Laby='Psione'
        LabX='ustar (m/s)'
        Title='Stable Error'
        title=trim(Title)
        Plotfile='Error'
        call MultiLine (Numpts,X3,Y,Numpts,X3,Y, &
        0.,1.,-5.,maxval(Y), &
        DevFlag,LabX,LabY,title,nv, 1, nv, Plotfile, legFlag, Legend)
      End If
!
      End If
!
!--------------------------------------------------------------------------
! scan 2-D ra1 error calculations
!--------------------------------------------------------------------------
      If (scantype=='RaError2D') then
      print *, 'scan 2-D ra1 error ', numpts, ' points'
!
      allocate(X3(numpts))
      allocate(Y(nv,numpts))
!
      sib(1)%drvr%tm = 301.0
      z2=BioTab(1)%z2
      ha=23.5
      Corb1=32.
      Corb2=202.
      U2=5.
      Dx=7./real(numpts-1)
      Dz=50.
!
      Hflux = -100.
      do k=1,nv
        ra1=-4.
        do j=1,numpts
          coef3=corb1*Hflux/sib(1)%drvr%tm/(z2-ha)
          Error=coef3*ra1**3+(u2*ra1)**2-corb2
!
          y(k,j)=error
          x3(j)=ra1
!
          ra1=ra1+Dx
        enddo
        Hflux=Hflux+Dz
      enddo
!
      If (PlotFlag==1) then
        Laby='Error (s/m)'
        LabX='ra1 (s/m)'
        Title='Ra1 Error'
        title=trim(Title)
        Plotfile='Error'
        call MultiLine (Numpts,X3,Y,Numpts,X3,Y, &
        minval(X3),maxval(X3),minval(Y),maxval(Y), &
        DevFlag,LabX,LabY,title,nv, 1, nv, Plotfile, legFlag, Legend)
      End If
!
      End If
!
!--------------------------------------------------------------------------
! scan 3-D ra1 error calculations
!--------------------------------------------------------------------------
      If (scantype=='RaError3D') then
      print *, 'scan 3-D ra1 error ', numpts, ' points'
!
      allocate(Z(nv,numpts,numpts))
!
      sib(1)%drvr%tm = 301.0
      z2=BioTab(1)%z2
      ha=23.5
      Corb1=32.
      Corb2=202.
      U2=5.
      Dx=800./real(numpts-1)
      Dy=40./real(numpts-1)
      Dz=1.
!
      U2=0.001
      do k=1,nv
        ra1=-20.
        do j=1,numpts
          Hflux=-400.
          do i=1,numpts
            coef3=corb1*Hflux/sib(1)%drvr%tm/(z2-ha)
            Error=coef3*ra1**3+(u2*ra1)**2-corb2
            z(k,i,j)=error
            Hflux=Hflux+Dx
          enddo
          ra1=ra1+Dy
        enddo
        U2=U2+Dz
      enddo
!
      If (PlotFlag==1) then
        LabX='Hflux (W/m**2)'
        LabY='ra1 (s/m)'
        Title='ra1 Error'
        Plotfile='Error'
        call Contour (Z,numpts,Numpts,-400.,-20., Dx, Dy,  &
        DevFlag, Scale,mincont,maxcont,base,NC,nv,MinBiome,MaxBiome, &
        LabX,LabY,Title, Plotfile)
      End If
!
      End If
!
!--------------------------------------------------------------------------
! scan Obhukov length
!--------------------------------------------------------------------------
      If (scantype=='Obhukov3D') then
      print *, 'scan 3-D stable rasite error ', numpts, ' points'
!
      allocate(Z(nv,numpts,numpts))
      z=0.
!
      zp_disp=26
      zwind = 100.0
      sib(1)%drvr%tm = 301.0
      rhoair = 1.225
      vkc = 0.41
      cpair = 1004.67
      gamma = 16.0
      g = 9.8
      Dy=5./real(numpts-1)
      Dx=450./real(numpts-1)
!
      Uest=0.001
      do j=1,numpts
        Hflux=-100.001
        do i=1,numpts
          Obhukov=-uest**3*rhoair*cpair*sib(1)%drvr%tm/(vkc*g*Hflux)
          zeta=(zwind-zp_disp)/Obhukov
          if (Hflux.gt.0.) Rich=zeta
          if (Hflux.lt.0.) Rich=zeta*(0.74+4.7*zeta)/(1.+4.7*zeta)**2.
          z(1,i,j)=abs(Rich)
          Hflux=HFlux+Dx
        enddo
        uest=uest+Dy
      enddo
!
      If (PlotFlag==1) then
        LabX='Hflux (W/m**2)'
        LabY='ustar (m/s)'
        Title='Richardson Number'
        Plotfile='Obk'
        call Contour (Z,numpts,Numpts,-100., 0., Dx, Dy,  &
        DevFlag, Scale,mincont,maxcont,base,NC,nv,MinBiome,MaxBiome, &
        LabX,LabY,Title, Plotfile)
      End If
!
      End If
!
!--------------------------------------------------------------------------
! scan 2_D U2
!--------------------------------------------------------------------------
      If (scantype=='U22D') then
      print *, 'scan U2 ', numpts, ' points'
!
      allocate(X3(numpts))
      allocate(Y(nv,numpts))
!
      sib(1)%drvr%tm = 301.0
      zp_disp=26
      z0=2.5
      Hflux = 100.
      rhoair = 1.225
      zwind = 100.0
      zmet = 100.0
      vkc = 0.41
      cpair = 1004.67
      gamma = 16.0
      g = 9.8
      G3=0.79
      corb2=200.
      zt=z2+11.875*z0
      argz=dlog((zwind-zp_disp)/z0)
      Dx=2./real(numpts-1)
      Dz=10.
!
      Hflux = 0.1
      do k=1,nv
        um=0.001
        do j=1,numpts

          ustarN=vkc*um/argz
          U2n=ustarn/vkc*alog((BioTab(1)%z2-zp_disp)/z0)

          Obhukov=-ustarN**3*rhoair*cpair*sib(1)%drvr%tm/(vkc*g*Hflux)
          xCon=(1.-gamma*(zwind-zp_disp)/Obhukov)**0.25
          psione=2.*alog((1.+xCon)/2.)+alog((1.+xCon*xCon)/2.) &
                -2.*atan(xCon)+pi/2.
!          psione=amin1(argz*0.75,psione)
          ustar=Um/(argz-psione)*vkc
!
          Obhukov=-ustar**3*rhoair*cpair*sib(1)%drvr%tm/(vkc*g*Hflux)
          xCon=(1.-gamma*(BioTab(1)%z2-zp_disp)/Obhukov)**0.25
          psione=2.*alog((1.+xCon)/2.)+alog((1.+xCon*xCon)/2.) &
                -2.*atan(xCon)+pi/2.
          psione=amin1(argz*0.75,psione)

          U2=Um*(alog((BioTab(1)%z2-zp_disp)/z0)-psione)/(argz-psione)
          u2=min(U2,Um)
          u2=max(u2,u2n)
!
          ra1=sqrt(corb2)/u2
          ra=raf+ra1
!
          y(k,j)=u2
          x3(j)=Um
!
          um=um+Dx
        enddo
        Hflux=Hflux+Dz
      enddo
!
      If (PlotFlag==1) then
        Laby='U2 (m/s)'
        LabX='Um (m/s)'
        Title='Hflux 0-120 W/m**2'
        title=trim(Title)
        Plotfile='Error'
        call MultiLine (Numpts,X3,Y,Numpts,X3,Y, &
        0.,maxval(X3),minval(Y),maxval(Y), &
        DevFlag,LabX,LabY,title,nv,  &
        minbiome, maxbiome, Plotfile, legFlag, Legend)
      End If
!
      End If
!
!--------------------------------------------------------------------------
! scan Turbulent-laminar flow transition
!--------------------------------------------------------------------------
      If (scantype=='HTL2D') then
      print *, 'scan Turb-laminar transition ', numpts, ' points'
!
      allocate(X3(numpts))
      allocate(Y(nv,numpts))
      X3=0.
      Y=0.
!
      sib(1)%drvr%tm = 301.0
      zp_disp=26
      z0=2.5
      um = 0.001
      rhoair = 1.225
      zwind = 100.0
      zmet = 100.0
      vkc = 0.41
      cpair = 1004.67
      gamma = 16.0
      g = 9.8
      argz=dlog((zwind-zp_disp)/z0)
      Dx=50./real(numpts-1)
!
      um=0.001
      do j=1,numpts
        hLT=-0.95*sib(1)%drvr%tm*rhoair*cpair/(2.0*4.7*g*(zwind-zp_disp))* &
              (2.0*um/3.0)**3*(vkc/argz)**2
        hTL=5.0*hLT
        ustar=vkc*um/argz
        y(1,j)=HLT
        y(2,j)=HTL
        x3(j)=ustar
        um=um+Dx
      enddo
!
      If (PlotFlag==1) then
        Laby='Ustar neutral (m/s)'
        LabX='Hflux (W/m**2)'
        Title='Turb-laminar trans flux'
        title=trim(Title)
        Plotfile='HTL'
        call MultiLine (Numpts,y(1,:),X3,Numpts,Y(2,:),X3, &
        -100.,0.,0.,5., &
        DevFlag,LabX,LabY,title,1, 1, 1, Plotfile, legFlag, Legend)
      End If
!
      End If
!
!--------------------------------------------------------------------------
! scan velocity from TKE
!--------------------------------------------------------------------------
      If (scantype=='TKE') then
      print *, 'scan TKE ', numpts, ' points'
!
      allocate(X3(numpts))
      allocate(Y(nv,numpts))
      X3=0.
      Y=0.
!
      sib(1)%drvr%tm = 301.0
      rhoair = 1.225
      cpair = 1004.67
      gamma = 16.0
      g = 9.8

      Dx=450./real(numpts-1)
!
      Hflux=0.001
      do j=1,numpts
        ustd=sqrt(2.*g*Hflux*600./3./sib(1)%drvr%tm/cpair/rhoair)
        y(1,j)=Ustd
        x3(j)=Hflux
        Hflux=Hflux+Dx
      enddo
!
      If (PlotFlag==1) then
        Laby='Ustd (m/s)'
        LabX='Hflux (W/m**2/s)'
        Title='Turbulent Velocity'
        title=trim(Title)
        Plotfile='TKE'
        call MultiLine (Numpts,X3,Y,Numpts,X3,Y, &
        0.,450.,minval(Y),maxval(Y), &
        DevFlag,LabX,LabY,title,nv, 1, 1, Plotfile, legFlag, Legend)
      End If
!
      End If
!
!--------------------------------------------------------------------------
! scan neutral friction velocity
!--------------------------------------------------------------------------
      If (scantype=='Ustar2D') then
      print *, 'scan nuetral friction velocity ', numpts, ' points'
!
      allocate(X3(numpts))
      allocate(Y(nv,numpts))
      X3=0.
      Y=0.
!
      corb2=200.
      z2=BioTab(1)%z2
      sib(1)%drvr%tm = 301.0
      zp_disp=26
      z0=2.5
      zwind = 100.0
      vkc = 0.41
      zt=z2+11.875*z0
      logZw=dlog((zwind-zp_disp)/z0)
      logZtZ2=alog((zt-zp_disp)/(z2-zp_disp))
      logZwZt=dlog((zwind-zp_disp)/(zt-zp_disp))
      Dx=5./real(numpts-1)
!
      um=0.001
      do j=1,numpts
        ustarN=vkc*um/logZw
        raf=(G3*logZtZ2+logZwZt)/(vkc*ustarN)
        U2N=Um-raf*ustarN**2
        y(1,j)=U2N
        y(2,j)=U2N+0.1
        y(3,j)=Um
        y(4,j)=sqrt(corb2)/U2N
        y(5,j)=sqrt(corb2)/Um
        x3(j)=um
        um=um+Dx
      enddo
!
      If (PlotFlag==1) then
        Laby='U* (m/s)'
        LabX='Um (m/s)'
        Title='Nuetral Friction Velocity'
        title=trim(Title)
        Plotfile='ustar'
        call MultiLine (Numpts,X3,Y,Numpts,X3,Y, &
        0.,maxval(X3),0.,50., &
        DevFlag,LabX,LabY,title,nv, 1, 5, Plotfile, legFlag, Legend)
      End If
!
      End If

    return
    end
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$





!
!--------------------------------------------------------------------
    function psi1(Dz) result(psione)
!--------------------------------------------------------------------
! Calculates the stability correction factor to velocity for momentum
! for the unstable case (Hflux>0) (Paulson, 1970)
!
    real psione ! output value of psitwo
    real Dz     ! input delta z from zero plane displacment
    real x   ! dummy variable
!
    x=(1-gamma*Dz/Obhukov)**.25
    psione=2.*alog((1.+x)/2.)+alog((1.+x*x)/2.)-2.*atan(x)+pi/2.
!
    end function psi1
!
!--------------------------------------------------------------------
    function psi2(Dz) result(psitwo)
!--------------------------------------------------------------------
! Calculates the stability correction factor to velocity for heat
! for the unstable case (Hflux>0) (Paulson, 1970)
!
    real psitwo ! output value of psitwo
    real Dz     ! input delta z from zero plane displacment
    real x      ! dummy variable
!
    x=(1-gamma*Dz/Obhukov)**.25
    psitwo=2*alog((1+x*x)/2.)
!
    end function psi2
