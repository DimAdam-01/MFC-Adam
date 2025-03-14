#:def Hardcoded2DVariables()

    real(wp) :: eps
    real(wp) :: r, rmax, gam, umax, p0
    real(wp) :: rhoH, rhoL, pRef, pInt, h, lam, wl, amp, intH, intL, alph
    real(wp) :: x1c,y1c,x2c,y2c,Rvortex,r1c,r2c,cvortex,u1c,u2c,v1c,v2c,U,pres! This are coordinates for vortices
    real(wp) :: y1,y2,y3,y4,y5,y6,y7,y8

    eps = 1e-9_wp

#:enddef

#:def Hardcoded2D()

    select case (patch_icpp(patch_id)%hcid) ! 2D_hardcoded_ic example case

    case (200)
        if (y_cc(j) <= (-x_cc(i)**3 + 1)**(1._wp/3._wp)) then
            ! Volume Fractions
            q_prim_vf(advxb)%sf(i, j, 0) = eps
            q_prim_vf(advxe)%sf(i, j, 0) = 1._wp - eps
            ! Denssities
            q_prim_vf(contxb)%sf(i, j, 0) = eps*1000._wp
            q_prim_vf(contxe)%sf(i, j, 0) = (1._wp - eps)*1._wp
            ! Pressure
            q_prim_vf(E_idx)%sf(i, j, 0) = 1000._wp
        end if
    case (202) ! Gresho vortex (Gouasmi et al 2022 JCP)
        r = ((x_cc(i) - 0.5_wp)**2 + (y_cc(j) - 0.5_wp)**2)**0.5_wp
        rmax = 0.2_wp

        gam = 1._wp + 1._wp/fluid_pp(1)%gamma
        umax = 2*pi*rmax*patch_icpp(patch_id)%vel(2)
        p0 = umax**2*(1._wp/(gam*patch_icpp(patch_id)%vel(2)**2) - 0.5_wp)

        if (r < rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -(y_cc(j) - 0.5_wp)*umax/rmax
            q_prim_vf(momxe)%sf(i, j, 0) = (x_cc(i) - 0.5_wp)*umax/rmax
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2/2._wp)
        else if (r < 2*rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -((y_cc(j) - 0.5_wp)/r)*umax*(2._wp - r/rmax)
            q_prim_vf(momxe)%sf(i, j, 0) = ((x_cc(i) - 0.5_wp)/r)*umax*(2._wp - r/rmax)
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2/2._wp + 4*(1 - (r/rmax) + log(r/rmax)))
        else
            q_prim_vf(momxb)%sf(i, j, 0) = 0._wp
            q_prim_vf(momxe)%sf(i, j, 0) = 0._wp
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*(-2 + 4*log(2._wp))
        end if
    case (203) ! Gresho vortex (Gouasmi et al 2022 JCP) with density correction
        r = ((x_cc(i) - 0.5_wp)**2._wp + (y_cc(j) - 0.5_wp)**2)**0.5_wp
        rmax = 0.2_wp

        gam = 1._wp + 1._wp/fluid_pp(1)%gamma
        umax = 2*pi*rmax*patch_icpp(patch_id)%vel(2)
        p0 = umax**2*(1._wp/(gam*patch_icpp(patch_id)%vel(2)**2) - 0.5_wp)

        if (r < rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -(y_cc(j) - 0.5_wp)*umax/rmax
            q_prim_vf(momxe)%sf(i, j, 0) = (x_cc(i) - 0.5_wp)*umax/rmax
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2._wp/2._wp)
        else if (r < 2*rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -((y_cc(j) - 0.5_wp)/r)*umax*(2._wp - r/rmax)
            q_prim_vf(momxe)%sf(i, j, 0) = ((x_cc(i) - 0.5_wp)/r)*umax*(2._wp - r/rmax)
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2/2._wp + 4._wp*(1._wp - (r/rmax) + log(r/rmax)))
        else
            q_prim_vf(momxb)%sf(i, j, 0) = 0._wp
            q_prim_vf(momxe)%sf(i, j, 0) = 0._wp
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2._wp*(-2._wp + 4*log(2._wp))
        end if

        q_prim_vf(contxb)%sf(i, j, 0) = q_prim_vf(E_idx)%sf(i, j, 0)**(1._wp/gam)
    case (204) ! Rayleigh-Taylor instability
        rhoH = 3._wp
        rhoL = 1._wp
        pRef = 1.e5_wp
        pInt = pRef
        h = 0.7_wp
        lam = 0.2_wp
        wl = 2._wp*pi/lam
        amp = 0.05_wp/wl

        intH = amp*sin(2._wp*pi*x_cc(i)/lam - pi/2._wp) + h

        alph = 0.5_wp*(1._wp + tanh((y_cc(j) - intH)/2.5e-3_wp))

        if (alph < eps) alph = eps
        if (alph > 1._wp - eps) alph = 1._wp - eps

        if (y_cc(j) > intH) then
            q_prim_vf(advxb)%sf(i, j, 0) = alph
            q_prim_vf(advxe)%sf(i, j, 0) = 1._wp - alph
            q_prim_vf(contxb)%sf(i, j, 0) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, 0) = (1._wp - alph)*rhoL
            q_prim_vf(E_idx)%sf(i, j, 0) = pref + rhoH*9.81_wp*(1.2_wp - y_cc(j))
        else
            q_prim_vf(advxb)%sf(i, j, 0) = alph
            q_prim_vf(advxe)%sf(i, j, 0) = 1._wp - alph
            q_prim_vf(contxb)%sf(i, j, 0) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, 0) = (1._wp - alph)*rhoL
            pInt = pref + rhoH*9.81_wp*(1.2_wp - intH)
            q_prim_vf(E_idx)%sf(i, j, 0) = pInt + rhoL*9.81_wp*(intH - y_cc(j))
        end if

    case (205) ! 2D lung wave interaction problem
        h = 0.0_wp           !non dim origin y
        lam = 1.0_wp         !non dim lambda
        amp = patch_icpp(patch_id)%a(2)         !to be changed later!       !non dim amplitude

        intH = amp*sin(2*pi*x_cc(i)/lam - pi/2) + h

        if (y_cc(j) > intH) then
            q_prim_vf(contxb)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(1)
            q_prim_vf(contxe)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(2)
            q_prim_vf(E_idx)%sf(i, j, 0) = patch_icpp(1)%pres
            q_prim_vf(advxb)%sf(i, j, 0) = patch_icpp(1)%alpha(1)
            q_prim_vf(advxe)%sf(i, j, 0) = patch_icpp(1)%alpha(2)
        end if

    case (206) ! 2D lung wave interaction problem - horizontal domain
        h = 0.0_wp           !non dim origin y
        lam = 1.0_wp         !non dim lambda
        amp = patch_icpp(patch_id)%a(2)

        intL = amp*sin(2*pi*y_cc(j)/lam - pi/2) + h

        if (x_cc(i) > intL) then        !this is the liquid
            q_prim_vf(contxb)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(1)
            q_prim_vf(contxe)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(2)
            q_prim_vf(E_idx)%sf(i, j, 0) = patch_icpp(1)%pres
            q_prim_vf(advxb)%sf(i, j, 0) = patch_icpp(1)%alpha(1)
            q_prim_vf(advxe)%sf(i, j, 0) = patch_icpp(1)%alpha(2)
        end if

    case (207) !2D Boundary Lido for multicomponent transport problems

      x1c=0.015_wp
      y1c=0.0125_wp
      x2c=0.015_wp
      y2c=0.0175_wp
      r1c=(x_cc(i)-x1c)**(2.0_wp)+(y_cc(j)-y1c)**(2.0_wp)
      r2c=(x_cc(i)-x2c)**(2.0_wp)+(y_cc(j)-y2c)**(2.0_wp)
      rvortex=0.0030_wp 
      cvortex=8.0_wp*10.0_wp**(4.0_wp)
      U=25.0_wp

      y1 = 0.2_wp - x_cc(i)/0.03_wp * 0.2_wp   ! H2
      y2 = 0.2_wp - x_cc(i)/0.03_wp * 0.2_wp   ! H
      y3 = 0.4_wp - x_cc(i)/0.03_wp * 0.4_wp   ! O2
      y4 =      x_cc(i)/0.03_wp * 0.2_wp       ! H2O
      y5 = 0.2_wp - x_cc(i)/0.03_wp * 0.2_wp   ! CH4
      y6 =      x_cc(i)/0.03_wp * 0.2_wp       ! CO
      y7 =      x_cc(i)/0.03_wp * 0.2_wp       ! CO2
      y8 =      x_cc(i)/0.03_wp * 0.4_wp       ! N2



      pres = 0.96325_wp * 10.0_wp**(5.0_wp) + 10.0_wp**(4.0_wp) * &
       ( exp( -r1c/(Rvortex)**2 ) + exp( -r2c/(Rvortex)**2 ) )
      u1c=-cvortex*((y_cc(j)-y1c))*exp(-r1c/(rvortex)**(2.0_wp)/2.0_wp)
      v1c=cvortex*((x_cc(j)-x1c))*exp(-r1c/(rvortex)**(2.0_wp)/2.0_wp)

      u2c=-cvortex*((y_cc(j)-y2c))*exp(-r2c/(rvortex)**(2.0_wp)/2.0_wp)
      v2c=cvortex*((x_cc(j)-x2c))*exp(-r2c/(rvortex)**(2.0_wp)/2.0_wp)


       q_prim_vf(1)%sf(i,j,0)= pres / ( 300.0_wp * gas_constant * ( &
           y1/molecular_weights(1) + y2/molecular_weights(2) + &
           y3/molecular_weights(3) + y4/molecular_weights(4) + &
           y5/molecular_weights(5) + y6/molecular_weights(6) + &
           y7/molecular_weights(7) + y8/molecular_weights(8) ) )


      q_prim_vf(2)%sf(i,j,0)=u+u1c+u2c
      q_prim_vf(3)%sf(i,j,0)=v1c+v2c



      q_prim_vf(4)%sf(i,j,0)=pres

      q_prim_vf(5)%sf(i,j,0)=1.0_wp
      q_prim_vf(6)%sf(i,j,0)=0.2_wp-x_cc(i)/0.03_wp*0.2_wp !H_2
      q_prim_vf(7)%sf(i,j,0)=0.2_wp-x_cc(i)/0.03_wp*0.2_wp !H
      q_prim_vf(8)%sf(i,j,0)=0.4_wp-x_cc(i)/0.03_wp*0.4_wp ! O2
      q_prim_vf(9)%sf(i,j,0)= x_cc(i)/0.03_wp*0.2_wp  !H2O
      q_prim_vf(10)%sf(i,j,0)=0.2_wp-x_cc(i)/0.03_wp*0.2_wp !CH4
      q_prim_vf(11)%sf(i,j,0)=x_cc(i)/0.03_wp*0.2_wp !CO
      q_prim_vf(12)%sf(i,j,0)=x_cc(i)/0.03_wp*0.2_wp !CO2
      q_prim_vf(13)%sf(i,j,0)=x_cc(i)/0.03_wp*0.4_wp !N2

   
    print *, x_cc(i)


    case default
        if (proc_rank == 0) then
            call s_int_to_str(patch_id, iStr)
            call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
        end if

    end select

#:enddef
