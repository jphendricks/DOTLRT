PROCEDURE calc_atm_rad(       frequency  : Double               ;
                       const  inp_height : Double               ;
                       const  inp_theta  : Double               ;
                       var    atm_inp    : atm_inp_type         ;
                       const  surf_inp   : surf_inp_type        ;
                       var    tb         : stokes_par_type     );

        {Specific intensity in watts/meter**2/steradian/hertz at for given height and frequency
        {for two orthogonal polarizations and surf_inp.num_angles propagation angles}
        {Propagation angle measured WRT zenith in deg, frequency in GHz, (geometric) height in kilometers above surface level}

        {Assumes planar stratified model with specularly reflecting surface, absorbing (but non-scattering) atmosphere, isotropic background}
        {Based on quadrature weights for integration WRT differental transmission elements}
        {Opacity increments are calculated from differential pressure elements where the hydrostatic equation is used
         to relate pressure differentials to height differentials. A first order correction for gravitational weakening
         at high altitudes is employed in the hydrostatic pressure-to-height differential relationship.}
        {Extinction coefficient used to calculate opacity increments.}

        {CONST
            {rad_earth = 6.3568E3     ;  {Radius of the earth (km)}
            {pi = 3.141592654         ;  {pumpkin}

        VAR
            u             : Double;   {secant of observation angle}
            i_spec_up     : stokes_par_type;
            tb_down       : stokes_par_type;
            i_spec_down   : stokes_par_type;
            theta         : Double;   {angle of propagation measured wrt zenith in deg}
            height        : Double;   {height of level in km at height}
            pressure      : Double;   {pressure of level in mb at height}
            temperature   : Double;   {temperature of level in K at height}
            vapor_density : Double;   {water vapor density of level in g/m**3 at height}
            lp,lq,lk0,la0 : Double;   {hydrometeor liquid size distribution parameters}
            ip,iq,ik0,ia0 : Double;   {hydrometeor liquid size distribution parameters}
            abs_o2        : Double;   {oxygen and nitrogen absorption in nepers/km at height}
            abs_h2o       : Double;   {water vapor absorption in nepers/km at height}
            abs_cloud_liq : Double;   {cloud liquid absorption in nepers/km at height}
            abs_cloud_ice : Double;   {cloud frozen absorption in nepers/km at height}
            scat_cloud_liq: Double;   {cloud liquid scattering in nepers/km at height}
            scat_cloud_ice: Double;   {cloud frozen scattering in nepers/km at height}
            g_cloud_liq   : Double;   {cloud liquid particle scattering asymmetry at height}
            g_cloud_ice   : Double;   {cloud frozen particle scattering asymmetry at height}
            ext_tot       : Double;   {total atmospheric extinction in nepers/km at height}
            bb_spec_int   : Double;   {specific black body intensity for either polarization at height in W/m**2/Hz/ster}
            num_levels    : level_type; {level counter}
            next_level    : level_type; {next level of atmospheric input data}
            last_level    : level_type; {last level of atmospheric input data}
            level,levelp1 : level_type; {level counter for quadrature}
            accum         : Double;
            vr            : Double;   {vertical surface reflectivity}
            hr            : Double;   {horizontal surface reflectivity}
            bb_surf_si    : Double;   {black body surface specific intensity}
            vesi          : Double;   {vertical emission specific intensity}
            hesi          : Double;   {horizontal emission specific intensity}
            cb_int        : Double;   {cosmic background specific intensity}
            llw,nlw       : Double;   {level interpolation weights}
            law,naw       : Double;   {angle interpolation weights}
            next_angle    : angle_type; {next angle of surface reflectivity data}
            last_angle    : angle_type; {last angle of surface reflectivity data}
            tauip1        : Double;   {atmospheric opacity}
            taui          : Double;   {atmospheric opacity}
            ti,tip1       : Double;   {atmospheric transmittance}
            i_spec        : stokes_par_type;

        BEGIN
            num_levels:=atm_inp.num_levels;
            {Calculate absorption and black-body spectral intensity for all levels}
            calc_tot_ext(atm_inp,frequency);
         (* calc_bb_spec_int(atm_inp,frequency); *)
            {Check for valid height}
            height:=inp_height;
            IF inp_height < 0 THEN height:=0;
            IF inp_height > atm_inp.prof[num_levels].height-0.0001 then height:=atm_inp.prof[num_levels].height-0.0001;
            {Determine level boundaries}
            next_level:=2;
            WHILE atm_inp.prof[next_level].height < height DO next_level:=next_level+1;
            last_level:=next_level-1;
            {Interpolate over atmospheric information at 'height'}
            llw:=(atm_inp.prof[next_level].height-height)/(atm_inp.prof[next_level].height-atm_inp.prof[last_level].height);
            nlw:=1.0-llw;
            temperature:=llw*atm_inp.prof[last_level].temperature + nlw*atm_inp.prof[next_level].temperature;
            pressure:=exp(llw*ln(atm_inp.prof[last_level].pressure) + nlw*ln(atm_inp.prof[next_level].pressure));
            vapor_density:=llw*atm_inp.prof[last_level].vapor_density + nlw*atm_inp.prof[next_level].vapor_density;
            lp:=0;
            lq:=0;
            lk0:=0;
            la0:=0;
            IF ((atm_inp.prof[last_level].cloud_liq_k0>0) AND (atm_inp.prof[next_level].cloud_liq_k0>0)) THEN
                BEGIN
                    lp:=llw*atm_inp.prof[last_level].cloud_liq_p + nlw*atm_inp.prof[next_level].cloud_liq_p;
                    lq:=llw*atm_inp.prof[last_level].cloud_liq_q + nlw*atm_inp.prof[next_level].cloud_liq_q;
                    lk0:=llw*atm_inp.prof[last_level].cloud_liq_k0 + nlw*atm_inp.prof[next_level].cloud_liq_k0;
                    la0:=llw*atm_inp.prof[last_level].cloud_liq_a0 + nlw*atm_inp.prof[next_level].cloud_liq_a0;
                END
            ELSE
                BEGIN
                    IF atm_inp.prof[last_level].cloud_liq_k0>0 THEN
                        BEGIN
                            lp:=atm_inp.prof[last_level].cloud_liq_p;
                            lq:=atm_inp.prof[last_level].cloud_liq_q;
                            lk0:=llw*atm_inp.prof[last_level].cloud_liq_k0;
                            la0:=atm_inp.prof[last_level].cloud_liq_a0;
                        END;
                    IF atm_inp.prof[next_level].cloud_liq_k0>0 THEN
                        BEGIN
                            lp:=atm_inp.prof[next_level].cloud_liq_p;
                            lq:=atm_inp.prof[next_level].cloud_liq_q;
                            lk0:=nlw*atm_inp.prof[next_level].cloud_liq_k0;
                            la0:=atm_inp.prof[next_level].cloud_liq_a0;
                        END;
                END;
            ip:=0;
            iq:=0;
            ik0:=0;
            ia0:=0;
            IF ((atm_inp.prof[last_level].cloud_ice_k0>0) AND (atm_inp.prof[next_level].cloud_ice_k0>0)) THEN
                BEGIN
                    ip:=llw*atm_inp.prof[last_level].cloud_ice_p + nlw*atm_inp.prof[next_level].cloud_ice_p;
                    iq:=llw*atm_inp.prof[last_level].cloud_ice_q + nlw*atm_inp.prof[next_level].cloud_ice_q;
                    ik0:=llw*atm_inp.prof[last_level].cloud_ice_k0 + nlw*atm_inp.prof[next_level].cloud_ice_k0;
                    ia0:=llw*atm_inp.prof[last_level].cloud_ice_a0 + nlw*atm_inp.prof[next_level].cloud_ice_a0;
                END
            ELSE
                BEGIN
                    IF atm_inp.prof[last_level].cloud_ice_k0>0 THEN
                        BEGIN
                            ip:=atm_inp.prof[last_level].cloud_ice_p;
                            iq:=atm_inp.prof[last_level].cloud_ice_q;
                            ik0:=llw*atm_inp.prof[last_level].cloud_ice_k0;
                            ia0:=atm_inp.prof[last_level].cloud_ice_a0;
                        END;
                    IF atm_inp.prof[next_level].cloud_ice_k0>0 THEN
                        BEGIN
                            ip:=atm_inp.prof[next_level].cloud_ice_p;
                            iq:=atm_inp.prof[next_level].cloud_ice_q;
                            ik0:=nlw*atm_inp.prof[next_level].cloud_ice_k0;
                            ia0:=atm_inp.prof[next_level].cloud_ice_a0;
                        END;
                END;
            {This is providing the option to switch between MRT and Rozenkranz absorption model}
            {m.klein January 1998}
            if mrt_checked then
              begin
                abs_o2:=o2abs(temperature,pressure,vapor_density,frequency) + absn2(temperature,pressure,frequency);
                abs_h2o:=abh2o(temperature,pressure,vapor_density,frequency);
              end
            else
              begin
                o2a_r(temperature,pressure,vapor_density,frequency,abs_o2);
                n2r(temperature,pressure,frequency,ext_tot);
                abs_o2:= abs_o2+ext_tot;
                h2o_r(temperature,pressure,vapor_density,frequency,abs_h2o);
              end;
            hydrometeor_ext(temperature,frequency,lp,lq,lk0,la0,abs_cloud_liq,scat_cloud_liq,g_cloud_liq,false);
            hydrometeor_ext(temperature,frequency,ip,iq,ik0,ia0,abs_cloud_ice,scat_cloud_ice,g_cloud_ice,true);
            ext_tot:=abs_o2+abs_h2o+abs_cloud_liq+scat_cloud_liq+abs_cloud_ice+scat_cloud_ice;
            bb_spec_int:=temperature; (* plank_int(frequency,temperature); *)

            {Check for valid angle}
            theta:=ABS(inp_theta);
            WHILE theta > 180 DO theta:=theta-180;
            {Check for direction of propagation}
            IF theta = 90 THEN
                BEGIN {Theta = 90 deg : Side propagating}
                    i_spec.vert:=bb_spec_int;
                    i_spec.horz:=bb_spec_int;
                END  {Side propagating}
            ELSE
                BEGIN
                    u:=0.5/cos(theta*pi/180);  {abs averaging weight of 1/2 multiplied here for computational efficiency}
                    u:=u*0.029276;  {Multiply by R/(Smx*g)=8.314/(28.96*9.80616). Value of g used for 45 degrees}
                    {Correct gravitational acceleration to latitude if known, else used default acceleration at 45 degrees}
                    IF ((atm_inp.latitude >= -90) AND (atm_inp.latitude <= 90)) THEN
                        u:=u*1.0026384/(1+0.0052885*sqr(sin(atm_inp.latitude*pi/180))-0.0000059*sqr(sin(atm_inp.latitude*pi/90)));
                    IF theta > 90 THEN
                        BEGIN {Theta > 90 deg : Downward propagating}
                            tauip1:= u * (atm_inp.prof[next_level].ext_tot*atm_inp.prof[next_level].temperature/atm_inp.prof[next_level].pressure +
                                          ext_tot*temperature/pressure) *
                                         (atm_inp.prof[next_level].pressure-pressure) *
                                         (1+(atm_inp.prof[next_level].height+height)/rad_earth);
                            ti:=exp(-tauip1);
                            accum:= (atm_inp.prof[next_level].bb_spec_int+bb_spec_int) * (1.0-ti);
                            FOR level:=next_level TO num_levels-1 DO
                                BEGIN
                                    levelp1:=level+1;
                                    tauip1:= tauip1 + u * (atm_inp.prof[levelp1].ext_tot*atm_inp.prof[levelp1].temperature/atm_inp.prof[levelp1].pressure +
                                                           atm_inp.prof[level].ext_tot*atm_inp.prof[level].temperature/atm_inp.prof[level].pressure) *
                                                          (atm_inp.prof[levelp1].pressure-atm_inp.prof[level].pressure) *
                                                          (1+(atm_inp.prof[levelp1].height+atm_inp.prof[level].height)/rad_earth);
                                    tip1:=exp(-tauip1);
                                    accum:= accum + (atm_inp.prof[level].bb_spec_int+atm_inp.prof[levelp1].bb_spec_int) * (ti-tip1);
                                    ti:=tip1;
                                END;
                            cb_int:=2.73; (* plank_int(frequency,2.73); *)
                            i_spec.vert:=accum/2 + cb_int*ti;
                            i_spec.horz:=i_spec.vert;
                        END   {Downward propagating}
                    ELSE
                        BEGIN {Theta < 90 deg : Upward propagating}
                            {Calculate downwelling radiation at ground level}                                          
                            calc_atm_rad(frequency,0,180-theta,atm_inp,surf_inp,tb_down);
                            i_spec_down.vert:=tb_down.vert; (* plank_int(frequency,tb_down.vert); *)
                            i_spec_down.horz:=tb_down.horz; (* plank_int(frequency,tb_down.horz); *)
                            {Determine angle boundaries}
                            next_angle:=2;
                            WHILE ((surf_inp.theta[next_angle] < theta) AND (next_angle < surf_inp.num_angles)) DO next_angle:=next_angle+1;
                            last_angle:=next_angle-1;
                            law:=(surf_inp.theta[next_angle]-theta)/(surf_inp.theta[next_angle]-surf_inp.theta[last_angle]);
                            naw:=1.0-law;
                            {Interpolate over surface information at 'theta'}
                            vr:=law*surf_inp.vr[last_angle] + naw*surf_inp.vr[next_angle];
                            hr:=law*surf_inp.hr[last_angle] + naw*surf_inp.hr[next_angle];
                            bb_surf_si:=surf_inp.surf_temp; (* plank_int(frequency,surf_inp.surf_temp); *)
                            vesi:=(1-vr)*bb_surf_si;
                            hesi:=(1-hr)*bb_surf_si;
                            {Calculate upwelling contribution }
                            taui:= u * (atm_inp.prof[last_level].ext_tot*atm_inp.prof[last_level].temperature/atm_inp.prof[last_level].pressure +
                                        ext_tot*temperature/pressure) *
                                       (atm_inp.prof[last_level].pressure-pressure) *
                                       (1+(atm_inp.prof[last_level].height+height)/rad_earth);
                            tip1:=exp(-taui);
                            ti:=1.0;
                            accum:= (atm_inp.prof[last_level].bb_spec_int+bb_spec_int) * (ti-tip1);
                            FOR level:=last_level-1 DOWNTO 1 DO
                                BEGIN
                                    levelp1:=level+1;
                                    taui:= taui + u * (atm_inp.prof[level].ext_tot*atm_inp.prof[level].temperature/atm_inp.prof[level].pressure +
                                                       atm_inp.prof[levelp1].ext_tot*atm_inp.prof[levelp1].temperature/atm_inp.prof[levelp1].pressure) * 
                                                      (atm_inp.prof[level].pressure-atm_inp.prof[levelp1].pressure) *
                                                      (1+(atm_inp.prof[level].height+atm_inp.prof[levelp1].height)/rad_earth);
                                    ti:=exp(-taui);
                                    accum:= accum + (atm_inp.prof[level].bb_spec_int+atm_inp.prof[levelp1].bb_spec_int) * (tip1-ti);
                                    tip1:=ti;
                                END;
                            i_spec_up.vert:=accum/2;
                            i_spec_up.horz:=i_spec_up.vert;
                            i_spec.vert:=i_spec_up.vert + (i_spec_down.vert*vr+vesi)*ti;
                            i_spec.horz:=i_spec_up.horz + (i_spec_down.horz*hr+hesi)*ti;
                        END;  {Upward propagating}
                END;  {Theta <> 90}
            {Convert specific intensities into equivalent brightness temperatures}
            tb:=i_spec; {R-J approx/constant scaling factors neglected throughout}
        END;  {calc_atm_rad}