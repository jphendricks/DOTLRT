  ! {Calculate passband quadrature frequencies}
  ! calc_passband_freq(passband_freq,num_freqs,chan_spec,num_sb_freqs);

  if( num_freqs > 0 ) then

    ! {Initialize weighting vector}
    temp_weight.surf:=0
    temp_weight.cos_back:=0
    do level = 1, atm_inp.num_levels
      temp_weight.atm[level] = 0
    end do

    ! {Compute weighting vector quadrature over passband frequencies using trapezoid rule}
    do nfreq = 1, num_freqs

      if( ( surf_inp.freq <> passband_freq[nfreq] ) .and. ( surf_inp.surf_type = 'O') ) then
        construct_surf(surf_inp,surf_inp.surf_subtype,surf_inp.surf_type,surf_inp.param1,passband_freq[nfreq],
                       surf_inp.param2, surf_inp.param3,surf_inp.num_angles)
      end if

      calc_mon_temp_weight_scat_iterative( passband_freq[nfreq], inp_height, inp_theta, atm_inp, surf_inp,
                                           frac_error, num_quad_angles, num_rpm_terms, mon_temp_weight_v, mon_temp_weight_h)

      ! polarization:
      v_pol = inp_pol * inp_pol
      h_pol = 1 - v_pol
      mon_temp_weight.surf = v_pol*mon_temp_weight_v.surf + h_pol*mon_temp_weight_h.surf
      mon_temp_weight.cos_back = v_pol*mon_temp_weight_v.cos_back + h_pol*mon_temp_weight_h.cos_back

      do level = 1, atm_inp.num_levels
        mon_temp_weight.atm[level] = v_pol*mon_temp_weight_v.atm[level] + h_pol*mon_temp_weight_h.atm[level]
      end do

      ! The mod operator returns the remainder obtained by dividing its operands:
      !     x mod y = x - int(x / y) * y
      ! if ((((nfreq MOD num_sb_freqs)=0) OR ((nfreq MOD num_sb_freqs)=1)) AND (num_sb_freqs > 1 )) THEN
      if( ( ( mod(nfreq,num_sb_freqs) == 0 ) .or. ( mod(nfreq,num_sb_freqs) == 1 ) ) &
         .and. ( num_sb_freqs > 1 ) ) )then
        temp_weight.surf = temp_weight.surf + 0.5 * mon_temp_weight.surf
        temp_weight.cos_back:=temp_weight.cos_back+0.5*mon_temp_weight.cos_back;
        do level = 1, atm_inp.num_levels
          temp_weight.atm[level] = temp_weight.atm[level] + 0.5 * mon_temp_weight.atm[level]
        end do
      else
        temp_weight.surf = temp_weight.surf + mon_temp_weight.surf
        temp_weight.cos_back = temp_weight.cos_back + mon_temp_weight.cos_back
        do level = 1, atm_inp.num_levels
          temp_weight.atm[level]:=temp_weight.atm[level]+mon_temp_weight.atm[level]
        end do
      end if

    end do ! nfreq = 1, num_freqs
    if( num_freqs > 1 ) then
      ! {Normalize weighting vector passband quadrature weights}
      if( num_sb_freqs = 1 ) THEN
        norm = num_freqs
      else
        norm = num_freqs * ( 1 - 1 / num_sb_freqs )
      end if
      temp_weight.surf = temp_weight.surf / norm
      temp_weight.cos_back = temp_weight.cos_back / norm
      do level = 1, atm_inp.num_levels
        temp_weight.atm[level] = temp_weight.atm[level] / norm
      end do
    end if
    ! {Record misc information}
    temp_weight.obs_lo_freq = chan_spec.lo_freq
    temp_weight.obs_if1_freq = chan_spec.if1_freq
    temp_weight.obs_if2_freq = chan_spec.if2_freq
    temp_weight.obs_bandwidth = chan_spec.bandwidth
    temp_weight.obs_height = inp_height
    temp_weight.obs_theta = inp_theta
    temp_weight.obs_pol = inp_pol
    temp_weight.inf = true
  end if ! num_freqs > 0