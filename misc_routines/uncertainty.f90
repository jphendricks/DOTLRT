!
! Active layer/surface deformation subroutine
            if(scantype=='ALT_def') then

! constants
              kroot(10)=5.5
	      rho_om_max=140.
	      poros_om=.9
	      root_depth=1.
              org_depth=0.18

	      call alt_from_deform(generic, satfrac,sand,kroot(10),m_om,&
                    root_depth, org_depth, rho_om_max, poros_om,var1, var2, var3, var4)
	      test=var2
	      test10(1)=var1
	      test10(2)=var2
	      test10(3)=var3
	      test10(4)=var4
	      ref_alt=var2
	      uncert_cum=0.
!
! deformation uncertainty
	      slope=ref_alt
	      junkvar=0.9*generic
	      if(junkvar<0.1) junkvar=0.1
	      call alt_from_deform(junkvar, satfrac,sand,kroot(10),m_om,&
                    root_depth, org_depth, rho_om_max, poros_om,var1, var2, var3, var4)
              slope=(slope-var2)/(generic-junkvar)
	      uncert(1)=slope*0.5
	      uncert_cum(1)=uncert(1)*uncert(1)	      
!
! satfrac uncertainty
              slope=ref_alt
	      junkvar=0.9*satfrac
	      call alt_from_deform(generic, junkvar,sand,kroot(10),m_om,&
                    root_depth, org_depth, rho_om_max, poros_om,var1, var2, var3, var4)
              slope=(slope-var2)/(satfrac-junkvar)
	      uncert(2)=slope*0.2
	      uncert_cum(2)=uncert_cum(1)+uncert(2)*uncert(2)
!
! organic matter uncertainty
              slope=ref_alt
	      junkvar=0.9*m_om
	      call alt_from_deform(generic, satfrac,sand,kroot(10),junkvar,&
                    root_depth, org_depth, rho_om_max, poros_om,var1, var2, var3, var4)
              slope=(slope-var2)/(m_om-junkvar)
	      uncert(3)=slope*5.
	      uncert_cum(3)=uncert_cum(2)+uncert(3)*uncert(3)
!
! porosity organic matter uncertainty
              slope=ref_alt
	      junkvar=0.9*poros_om
	      call alt_from_deform(generic, satfrac,sand,kroot(10),m_om,&
                    root_depth, org_depth, rho_om_max, junkvar,var1, var2, var3, var4)
              slope=(slope-var2)/(poros_om-junkvar)
	      uncert(4)=slope*.05
	      uncert_cum(4)=uncert_cum(3)+uncert(4)*uncert(4)
!
! density organic matter max uncertainty
              slope=ref_alt
	      junkvar=0.9*rho_om_max
	      call alt_from_deform(generic, satfrac,sand,kroot(10),m_om,&
                    root_depth, org_depth, junkvar, poros_om,var1, var2, var3, var4)
              slope=(slope-var2)/(rho_om_max-junkvar)
	      uncert(5)=slope*10.
	      uncert_cum(5)=uncert_cum(4)+uncert(5)*uncert(5)
!
! sand fraction uncertainty
              slope=ref_alt
	      junkvar=0.9*sand
	      call alt_from_deform(generic, satfrac,junkvar,kroot(10),m_om,&
                    root_depth, org_depth, rho_om_max, poros_om,var1, var2, var3, var4)
              slope=(slope-var2)/(sand-junkvar)
	      uncert(6)=slope*5.
	      uncert_cum(6)=uncert_cum(5)+uncert(6)*uncert(6)
!
! kroot uncertainty
              slope=ref_alt
	      junkvar=0.9*kroot(10)
	      call alt_from_deform(generic, satfrac,sand,junkvar,m_om,&
                    root_depth, org_depth, rho_om_max, poros_om,var1, var2, var3, var4)
              slope=(slope-var2)/(kroot(10)-junkvar)
	      uncert(7)=slope*.1
	      uncert_cum(7)=uncert_cum(6)+uncert(7)*uncert(7)
!
! rooting depth uncertainty
              slope=ref_alt
	      junkvar=0.9*root_depth
	      call alt_from_deform(generic, satfrac,sand,kroot(10),m_om,&
                    junkvar, org_depth, rho_om_max, poros_om,var1, var2, var3, var4)
              slope=(slope-var2)/(root_depth-junkvar)
	      uncert(8)=slope*.1
	      uncert_cum(8)=uncert_cum(7)+uncert(8)*uncert(8)
!
! absolute uncertainty
	      !test10(1)=ref_alt
	      !do n=1,8
	      !  test10(n+1)=abs(uncert(n))
              !enddo
!
! uncertainty relative to ALT
	      !test10(1)=ref_alt
	      !do n=1,8
	      !  test10(n+1)=abs(uncert(n))/ref_alt*100.
	      !  test10(n+1)=min(test10(n+1),100.)
              !enddo
!
! cummulative uncertainty (total uncertainty all variables is test10(8))
	      test10(1)=ref_alt
	      do n=1,8
	        test10(n+1)=sqrt(uncert_cum(n))
              enddo
!
! relative contribution uncertainty
	      !test10(1)=ref_alt
	      !test10(2)=sqrt(uncert_cum(1))/sqrt(uncert_cum(8))*100.
	      !do n=2,8
	      !  test10(n+1)=(sqrt(uncert_cum(n))-sqrt(uncert_cum(n-1)))/sqrt(uncert_cum(8))*100.
              !enddo
	    endif
