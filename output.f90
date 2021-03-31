! output.f90
! DOTLRTv1.0
! 24 May 2005    Bob Weber
subroutine output( jchan, atm_data, Tbo, Tb_obs, dTb_dT_obs, dTb_dp_obs, dTb_dq_obs, dTb_dw_obs )
  use variables
  implicit none
  type(prof_data_type)    :: atm_data
  integer j, ilr, ilr1, k, i, hydrometeor_phase, jchan
  real(8) z, zout(0:max_num_levels)
  real(8), dimension(0:max_num_levels,max_num_quad_angles) :: Tb_obs
  real(8), dimension(max_num_levels,max_num_quad_angles) :: dTb_dT_obs
  real(8), dimension(max_num_levels,max_num_quad_angles) :: dTb_dp_obs
  real(8), dimension(max_num_levels,max_num_quad_angles) :: dTb_dq_obs
  real(8), dimension(max_num_levels,max_num_quad_angles,5) :: dTb_dw_obs
  real(8) Tbo, freqghz
  freqghz = instr_spec%chan(jchan)%lo_freq
  do i0 = 1, nangover2
    ! Product of geophysical and radiation Jacobians:
    z=0.d0
    do ilr1=1,nlr1
      z=z+h1(ilr1)
      zout(ilr1) = z
    end do
    jrec_tot(i0) = jrec_tot(i0) + 1
    write(tot_unit(i0),rec=jrec_tot(i0)) &
          sngl(atm_data%lat), sngl(atm_data%lng), nlr1, &
          sngl(freqghz), i0, sngl(dacos(cs(i0))*180.0d0/pi), &
        ( sngl(zout(ilr1)), sngl(dTb_dT(ilr1,i0,1)), sngl(dTb_dp(ilr1,i0,1)), &
          sngl(dTb_dq(ilr1,i0,1)), &
        ( sngl(dTb_dw(ilr1,i0,hydrometeor_phase,1)), &
                      hydrometeor_phase = 1, number_h2o_phases ), &
                    ilr1 = 1, nlr1 )
    z=0.d0
    zout(0) = z
    DO ilr1=1,nlr1
      z=z+altitude1(ilr1+1) - altitude1(ilr1)
      zout(ilr1) = z
    end do
              if( i0 == 1 ) then
                ! Temperatures:
                jrec_rt(3,i0) = jrec_rt(3,i0) + 1
                write( rt_unit(3,i0),rec=jrec_rt(3,i0)) sngl(atm_data%lat), sngl(atm_data%lng), &
                                            nlr1, i0, sngl(dacos(cs(i0))*180.0d0/pi), &
                                        sngl(zout(0)), sngl(surf_inp%surf_temp), &
                             ( sngl(zout(ilr1)), sngl(temperature1(ilr1)), ilr1 = 1, nlr1 )
                ! Gaseous absorption and combined hydrometeor ansorption/scattering coefficients:
                jrec_rt(4,i0) = jrec_rt(4,i0) + 1
                write( rt_unit(4,i0),rec=jrec_rt(4,i0)) sngl(atm_data%lat), sngl(atm_data%lng), &
                                               nlr1, i0, sngl(dacos(cs(i0))*180.0d0/pi), &
                                      (sngl(zout(ilr1)), sngl(al_gas1(ilr1)), &
                        sngl(scat_cloud1(ilr1)), sngl(abs_cloud1(ilr1)), sngl(abs_O2_1(ilr1)), &
                        sngl(abs_H2O_1(ilr1)), ilr1 = 1, nlr1 )
              end if
              z=0.d0
              DO ilr1=0, nlr1
                IF(ilr1 /= 0) z=z+h1(ilr1)
                zout(ilr1) = z
              end do
              ! Brightness temperatures:
              jrec_rt(1,i0) = jrec_rt(1,i0) + 1
              write( rt_unit(1,i0),rec=jrec_rt(1,i0)) sngl(atm_data%lat), sngl(atm_data%lng), &
                                            nlr1, i0, sngl(dacos(cs(i0))*180.0d0/pi), sngl(freqghz), &
                                           ( sngl(zout(ilr1)), sngl(tb_pl(ilr1,i0)), &
                                           sngl(tb_mn(ilr1,i0)), ilr1 = 0, nlr1 )
              z=0.d0
              DO ilr1=1,nlr1
                z=z+h1(ilr1)
                zout(ilr1) = z
              end do
              ! Radiation Jacobian:
              jrec_rt(2,i0) = jrec_rt(2,i0) + 1
                write( rt_unit(2,i0),rec=jrec_rt(2,i0)) sngl(atm_data%lat), sngl(atm_data%lng), &
                                              (nlr1-1), i0, sngl(dacos(cs(i0))*180.0d0/pi), sngl(freqghz), &
                                            (sngl(zout(ilr1)),                                    &
                        sngl(dTbdKa(ilr1,i0,1)), sngl(dTbdKa(ilr1,i0,2)), &
                        sngl(dTbdTr(ilr1,i0,1)), sngl(dTbdTr(ilr1,i0,2)), &
                        sngl(dTbdKsliq(ilr1,i0,1)), sngl(dTbdKsliq(ilr1,i0,2)), &
                        sngl(dTbdgliq(ilr1,i0,1)), sngl(dTbdgliq(ilr1,i0,2)), &
                        sngl(dTbdKsrn(ilr1,i0,1)), sngl(dTbdKsrn(ilr1,i0,2)), &
                        sngl(dTbdgrn(ilr1,i0,1)), sngl(dTbdgrn(ilr1,i0,2)), &
                        sngl(dTbdKsice(ilr1,i0,1)), sngl(dTbdKsice(ilr1,i0,2)), &
                        sngl(dTbdgice(ilr1,i0,1)), sngl(dTbdgice(ilr1,i0,2)), &
                        sngl(dTbdKssnow(ilr1,i0,1)), sngl(dTbdKssnow(ilr1,i0,2)), &
                        sngl(dTbdgsnow(ilr1,i0,1)), sngl(dTbdgsnow(ilr1,i0,2)), &
                        sngl(dTbdKsgrpl(ilr1,i0,1)), sngl(dTbdKsgrpl(ilr1,i0,2)), &
                        sngl(dTbdggrpl(ilr1,i0,1)), sngl(dTbdggrpl(ilr1,i0,2)), &
                                                       ilr1 = 1, nlr1)
            end do ! i0
            ! Geophysical Jacobian:
            jrec_geo = jrec_geo + 1
            write(geo_unit,rec=jrec_geo) sngl(atm_data%lat), sngl(atm_data%lng), nlr1, sngl(freqghz), &
                                           ( sngl(zout(ilr1)), &
                                             sngl(dKab_dT(ilr1)), &
                                             sngl(dKab_dp(ilr1)), &
                                             sngl(dKab_dq(ilr1)), &
                                           ( sngl(dKsc_dT(ilr1,hydrometeor_phase)), &
                                             sngl(  dg_dT(ilr1,hydrometeor_phase)), &
                                              hydrometeor_phase = 1, number_h2o_phases ), &
                                           ( sngl(dKab_dw(ilr1,hydrometeor_phase)), &
                                             sngl(dKsc_dw(ilr1,hydrometeor_phase)), &
                                             sngl(  dg_dw(ilr1,hydrometeor_phase)), &
                                              hydrometeor_phase = 1, number_h2o_phases ), ilr1 = 1, nlr1 )
            ! Hydrometeor absorption and scattering coefficients and asymmetry factor for all phases:
            z=0.d0
            do ilr1=1,nlr1
              z=z+h1(ilr1)
              zout(ilr1) = z
            end do
            jrec_Ksa = jrec_Ksa + 1
  write( Ksa_unit,rec=jrec_Ksa) sngl(atm_data%lat), sngl(atm_data%lng), nlr1, sngl(freqghz), &
        ( sngl(      zout(ilr1)), &
        ( sngl(hydro_prof(ilr1,hydrometeor_phase)%cloudab), &
          sngl(hydro_prof(ilr1,hydrometeor_phase)%cloudsc), &
          sngl(hydro_prof(ilr1,hydrometeor_phase)%cloudg),  &
              hydrometeor_phase = 1, number_h2o_phases ), &
              ilr1 = 1, nlr1)
end subroutine output
