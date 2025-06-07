#:def Hardcoded2DVariables()
    real(wp) :: eps
    real(wp) :: r, rmax, gam, umax, p0,delta_omega
    real(wp) :: rhoH, rhoL, pRef, pInt, h, lam, wl, amp, intH, intL, alph
    real(wp) :: x1c,y1c,x2c,y2c,Rvortex,r1c,r2c,cvortex,u1c,u2c,v1c,v2c,u,pres! This are coordinates for vortices
    real(wp) :: y1,y2,y3,y4,y5,y6,y7,y8,Lx,Ly,Sl,lambda,temp_val
    real(wp) :: YN2_O,YN2_F,YH2_F,YH2_O,YO2_O,YO2_F,YH2O_F,YH2O_O, YH_F,YH_O,YO_F,YO_O
    real(wp) ::  YOH_F,YOH_O,YHO2_F,YHO2_O,YH2O2_F,YH2O2_O,YAR_O,YAR_F
      real(wp), parameter :: U_fuel = 972.0_wp      ! m/s
      real(wp), parameter :: U_oxidizer = 1633.0_wp
      real(wp), parameter :: beta_deg = 33.0_wp
      real(wp), parameter :: beta     = beta_deg * acos(-1.0_wp)/180.0_wp
      real(wp), parameter :: T_F =  545.0_wp
      real(wp), parameter :: T_O = 1475.0_wp
      real(wp), parameter :: rho_F = 0.354187_wp   ! fuel‐stream density (H₂/N₂ @ 545 K)
      real(wp), parameter :: rho_O = 0.202889_wp   ! oxidizer‐stream density (vitiated‐air @ 1475 K)

    integer  :: Nx1,Nx2,Ny1,Ny2,Ny3,Ny4,Ny5,Ny6,mpi_size,ierr,global_idx
    integer  :: mpi_rank,lol,idx,stdout,i_min,i_max,local_x_min,local_x_max
integer :: local_start, local_end, local_n

  integer, parameter :: nFiles = 15
  integer, parameter :: xRows  = 1201
  integer, parameter :: Nrows  = xRows

  ! Variables for file reading and domain data
  real(wp)                ::  domain_start, domain_end, x_step
  real(wp), dimension(Nrows)            :: x_coords
  real(wp), dimension(Nrows, nFiles)      :: stored_values
  character(len=100), dimension(nFiles)   :: fileNames
  character(len=20)       :: file_num_str
  character(len=6), parameter :: zeros_default = "045360"
  character(len=*), parameter :: init_dir = "/home/pain/MFC-Adam/examples/1D_reactive_shocktube/D/"
  integer                 :: f, iter, ios, unit
  logical                 :: files_loaded

  files_loaded=.false.

  eps = 1e-9_wp

    do f = 1, nFiles
        ! Convert file number to string with proper formatting
        if (f < 10) then
            write(file_num_str, '(I1)') f  ! Single digit
        else
            write(file_num_str, '(I2)') f  ! Double digit
            ! For more than 99 files, you might need to adjust this format
        end if
        
        ! Create the filename with the pattern "prim.X.00.000000.dat"
        if (f == 15) then
            fileNames(f) = trim(init_dir) // "T.15.00." // zeros_default // ".dat"
        else
            fileNames(f) = trim(init_dir) // "prim." // trim(file_num_str) // ".00." // zeros_default // ".dat"
        end if

    end do



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

      x1c=0.1_wp
      y1c=0.1_wp
 !     x2c=0.015_wp
  !    y2c=0.0175_wp
      r1c=(x_cc(i)-x1c)**(2.0_wp)+(y_cc(j)-y1c)**(2.0_wp)
    !  r2c=(x_cc(i)-x2c)**(2.0_wp)+(y_cc(j)-y2c)**(2.0_wp)
      rvortex=0.015_wp 
      cvortex=-0.06_wp
      u=1234.81_wp

      y1 = 0.2_wp - x_cc(i)/0.2_wp * 0.2_wp   ! H2
      y2 = 0.2_wp - x_cc(i)/0.2_wp * 0.2_wp   ! H
      y3 = 0.4_wp - x_cc(i)/0.2_wp * 0.4_wp   ! O2
      y4 =      x_cc(i)/0.2_wp * 0.2_wp       ! H2O
      y5 = 0.2_wp - x_cc(i)/0.2_wp * 0.2_wp   ! CH4
      y6 =      x_cc(i)/0.2_wp * 0.2_wp       ! CO
      y7 =      x_cc(i)/0.2_wp * 0.2_wp       ! CO2
      y8 =      x_cc(i)/0.2_wp * 0.4_wp       ! N2



  
      u1c=cvortex*((y_cc(j)-y1c))*exp(-r1c/(rvortex**2.0_wp)/2.0_wp) 
      v1c=-cvortex*((x_cc(i)-x1c))*exp(-r1c/(rvortex**2.0_wp)/2.0_wp) 

     ! u2c=-cvortex*((y_cc(j)-y2c))*exp(-r2c/(rvortex)**(2.0_wp)/2.0_wp)
     ! v2c=cvortex*((x_cc(j)-x2c))*exp(-r2c/(rvortex)**(2.0_wp)/2.0_wp)

    if (1 .eq. 207) then 
    pres = 0.96325_wp * 10.0_wp**(5.0_wp) * (1.0_wp+sqrt(u1c**2.0_wp+v1c**2.0_wp))

       q_prim_vf(1)%sf(i,j,0)= pres / ( 300.0_wp * gas_constant )!* ( &
     !      y1/molecular_weights(1) + y2/molecular_weights(2) + &
     !      y3/molecular_weights(3) + y4/molecular_weights(4) + &
     !      y5/molecular_weights(5) + y6/molecular_weights(6) + &
     !      y7/molecular_weights(7) + y8/molecular_weights(8) ) )


      q_prim_vf(2)%sf(i,j,0)=u+u1c
      q_prim_vf(3)%sf(i,j,0)=v1c



      q_prim_vf(4)%sf(i,j,0)=pres

      q_prim_vf(5)%sf(i,j,0)=1.0_wp
      q_prim_vf(6)%sf(i,j,0)=0.2_wp-x_cc(i)/0.2_wp*0.2_wp !H_2
      q_prim_vf(7)%sf(i,j,0)=0.2_wp-x_cc(i)/0.2_wp*0.2_wp !H
      q_prim_vf(8)%sf(i,j,0)=0.4_wp-x_cc(i)/0.2_wp*0.4_wp ! O2
      q_prim_vf(9)%sf(i,j,0)= x_cc(i)/0.2_wp*0.2_wp  !H2O
      q_prim_vf(10)%sf(i,j,0)=0.2_wp-x_cc(i)/0.2_wp*0.2_wp !CH4
      q_prim_vf(11)%sf(i,j,0)=x_cc(i)/0.2_wp*0.2_wp !CO
      q_prim_vf(12)%sf(i,j,0)=x_cc(i)/0.2_wp*0.2_wp !CO2
      q_prim_vf(13)%sf(i,j,0)=x_cc(i)/0.2_wp*0.4_wp !N2
    end if

   
 !@!   print *, pres

    case (208) !2D Boundary Lido for multicomponent transport problems

      x1c=0.1_wp
      y1c=0.1_wp
 !     x2c=0.015_wp
  !    y2c=0.0175_wp
      r1c=(x_cc(i)-x1c)**(2.0_wp)+(y_cc(j)-y1c)**(2.0_wp)
    !  r2c=(x_cc(i)-x2c)**(2.0_wp)+(y_cc(j)-y2c)**(2.0_wp)
      rvortex=0.015_wp 
      cvortex=-0.06_wp
      u=450.81_wp

    

  
      u1c=cvortex*((y_cc(j)-y1c))/rvortex**(2.0_wp)*exp(-r1c/(rvortex**2.0_wp)/2.0_wp) 
      v1c=-cvortex*((x_cc(i)-x1c))/rvortex**(2.0_wp)*exp(-r1c/(rvortex**2.0_wp)/2.0_wp) 

     ! u2c=-cvortex*((y_cc(j)-y2c))*exp(-r2c/(rvortex)**(2.0_wp)/2.0_wp)
     ! v2c=cvortex*((x_cc(j)-x2c))*exp(-r2c/(rvortex)**(2.0_wp)/2.0_wp)
     pres = 0.96325_wp * 10.0_wp**(5.0_wp) * (1.0_wp+sqrt(u1c**2.0_wp+v1c**2.0_wp))

     q_prim_vf(1)%sf(i,j,0)= pres / ( 300.0_wp * 287.0_wp)
     q_prim_vf(2)%sf(i,j,0)=u+u1c
     q_prim_vf(3)%sf(i,j,0)=v1c
     q_prim_vf(4)%sf(i,j,0)=pres
     q_prim_vf(5)%sf(i,j,0)=1.0_wp

    case (209)
   
        ! Determine the rank’s local domain.
        ! (Assuming local_x_min and local_x_max are known from the simulation domain and rank offset)
        ! For example, they might be computed as:
        ! local_x_min = simulation_x(1) + rank_offset
        ! local_x_max = local_x_min + (n_local_cells - 1) * delta_x
        !
        ! Compute indices within the file corresponding to the rank’s domain:
                
        ! Now read the necessary portion of each file for this rank
        if (.not. files_loaded) then
            ! Print status message
            print *, "Loading all data files..."
            
            do f = 1, nFiles
              ! Open the file for reading
                open(newunit=unit, file=trim(fileNames(f)), status='old', action='read', iostat=ios)
                if (ios /= 0) then
                    print *, "Error opening file: ", trim(fileNames(f))
                    cycle  ! Skip this file on error
                endif
                
                ! Read all rows at once into memory
                do iter = 1, nRows
                    read(unit, *, iostat=ios) x_coords(iter), stored_values(iter, f)
                    if (ios /= 0) then
                        print *, "Error reading file ", trim(fileNames(f)), " at row ", iter
                        exit  ! Exit loop on error
                    endif
                end do
                close(unit)
            end do
            
            ! Calculate domain information for mapping
            domain_start = x_coords(1)
            domain_end = x_coords(nRows)
            x_step = (domain_end - domain_start) / (nRows - 1)
            
            print *, "All data files loaded. Domain range:", domain_start, "to", domain_end
            files_loaded = .true.
        endif            
    ! Now, for each cell in the simulation, check if its coordinate falls within the file domain
    if (x_cc(i).ge. domain_start-x_step/2.0_wp.and. x_cc(i) .le. domain_end) then
        ! Calculate the index in the file data array corresponding to x_cc(i)
        idx = nint((x_cc(i) - domain_start) / x_step) + 1
        
        ! Ensure idx is within valid bounds
        if (idx < 1) idx = 1
        if (idx > Nrows) idx = Nrows

        ! Assign values from stored data for each file (with your small adjustment for f > 2)
        do f = 1, nFiles - 1
            if (f > 2) then
                lol = 1
            else 
                lol = 0
            end if
            q_prim_vf(f+lol)%sf(i, j, 0) = stored_values(idx, f)
        enddo

        ! Set element 3 explicitly to zero
        q_prim_vf(3)%sf(i, j, 0) = 0.0_wp
    endif

    ! Create rectangular pockets with unburnt fuel
if (x_cc(i) < 0.05200_wp .and. x_cc(i) > 0.0470_wp) then
    ! First pocket
    if (y_cc(j) > 0.06_wp/4.0_wp .and. y_cc(j) < 0.06_wp/4.0_wp+0.06_wp/9.0_wp) then
        do f = 1, sys_size-1
            if (f /= 2 ) then
                ! Apply the same lol offset
                if (f > 2) then
                    lol = 1
                else 
                    lol = 0
                end if
                
                ! Special handling for pressure (index 4 in 2D system)
                    q_prim_vf(f+lol)%sf(i, j, 0) = stored_values(nRows-5, f)
            end if
        end do
        
    ! Second pocket
    else if (y_cc(j) > 0.06_wp/2.0_wp-0.06_wp/18.0_wp .and. y_cc(j) < 0.06_wp/2.0_wp+0.06_wp/18.0_wp) then
        do f = 1, sys_size-1
            if (f /= 2 ) then
                ! Apply the same lol offset
                if (f > 2) then
                    lol = 1
                else 
                    lol = 0
                end if
                
                ! Special handling for pressure
                    q_prim_vf(f+lol)%sf(i, j, 0) = stored_values(nRows-5, f)
            end if
        end do
        
    ! Third pocket
    else if (y_cc(j) < 3.0_wp*0.06_wp/4.0_wp .and. y_cc(j) > 3.0_wp*0.06_wp/4.0_wp-0.06_wp/9.0_wp) then
        do f = 1, sys_size-1
            if (f /= 2  )then
                ! Apply the same lol offset
                if (f > 2) then
                    lol = 1
                else 
                    lol = 0
                end if
                
                ! Special handling for pressure
                    q_prim_vf(f+lol)%sf(i, j, 0) = stored_values(nRows-5, f)
            end if
        end do
    end if
end if                        !  q_T_sf%sf(i, j, 0) = stored_values(i+1, nFiles)


   case (210) !2D Boundary Lido for multicomponent transport problems

        if (.not. files_loaded) then
            ! Print status message
            print *, "Loading all data files..."
            
            do f = 1, nFiles
              ! Open the file for reading
                open(newunit=unit, file=trim(fileNames(f)), status='old', action='read', iostat=ios)
                if (ios /= 0) then
                    print *, "Error opening file: ", trim(fileNames(f))
                    cycle  ! Skip this file on error
                endif
                
                ! Read all rows at once into memory
                do iter = 1, nRows
                    read(unit, *, iostat=ios) x_coords(iter), stored_values(iter, f)
                    if (ios /= 0) then
                        print *, "Error reading file ", trim(fileNames(f)), " at row ", iter
                        exit  ! Exit loop on error
                    endif
                end do
                close(unit)
            end do
            
            ! Calculate domain information for mapping
            domain_start = x_coords(1)
            domain_end = x_coords(nRows)
            x_step = (domain_end - domain_start) / (nRows - 1)
            
            print *, "All data files loaded. Domain range:", domain_start, "to", domain_end
            files_loaded = .true.

        endif
        
        ! Simple direct mapping - find the closest index without interpolation
        idx = nint((x_cc(i) - domain_start) / x_step) + 1
        
        ! Boundary protection
        ! if (idx < 1) idx = 1
        ! if (idx > nRows) idx = nRows
        
        ! Assign values directly from stored data for each file
 do f = 1, nFiles-1
    if (f > 2) then
        lol = 1
    else 
        lol = 0
    end if
    q_prim_vf(f+lol)%sf(i, j, 0) = stored_values(i+1, f)
end do

! Set element 3 to zero (as requested)
!q_prim_vf(3)%sf(i, j, 0) = 0.0_wp


      x1c=0.0027_wp
      y1c=0.005_wp
      x2c=0.0027_wp
      y2c=0.003_wp
      r1c=(x_cc(i)-x1c)**(2.0_wp)+(y_cc(j)-y1c)**(2.0_wp)
      r2c=(x_cc(i)-x2c)**(2.0_wp)+(y_cc(j)-y2c)**(2.0_wp)
      rvortex=0.0005_wp 
      cvortex=6000.0_wp

      u1c=-cvortex*((y_cc(j)-y1c))*exp(-r1c/(2.0_wp*rvortex**2.0_wp)) 
      v1c=cvortex*((x_cc(i)-x1c))*exp(-r1c/(2.0_wp*rvortex**2.0_wp))

      u2c=cvortex*((y_cc(j)-y2c))*exp(-r2c/(2.0_wp*rvortex**2.0_wp))
      v2c=-cvortex*((x_cc(i)-x2c))*exp(-r2c/(2.0_wp*rvortex**2.0_wp))
      q_prim_vf(2)%sf(i,j,0)=q_prim_vf(2)%sf(i,j,0)+u1c+u2c
      q_prim_vf(3)%sf(i,j,0)=v1c+v2c
  

  case (211) !2D Boundary Lido for multicomponent transport problems

        if (.not. files_loaded) then
            ! Print status message
            print *, "Loading all data files..."

            do f = 1, nFiles
              ! Open the file for reading
                open(newunit=unit, file=trim(fileNames(f)), status='old', action='read', iostat=ios)
                if (ios /= 0) then
                    print *, "Error opening file: ", trim(fileNames(f))
                    cycle  ! Skip this file on error
                endif

                ! Read all rows at once into memory
                do iter = 1, nRows
                    read(unit, *, iostat=ios) x_coords(iter), stored_values(iter, f)
                    if (ios /= 0) then
                        print *, "Error reading file ", trim(fileNames(f)), " at row ", iter
                        exit  ! Exit loop on error
                    endif
                end do
                close(unit)
            end do

            ! Calculate domain information for mapping
            domain_start = x_coords(1)
            domain_end = x_coords(nRows)
            x_step = (domain_end - domain_start) / (nRows - 1)

            print *, "All data files loaded. Domain range:", domain_start, "to", domain_end
            files_loaded = .true.

        endif

        ! Simple direct mapping - find the closest index without interpolation
        idx = nint((x_cc(i) - domain_start) / x_step) + 1

        ! Boundary protection
        ! if (idx < 1) idx = 1
        ! if (idx > nRows) idx = nRows

        ! Assign values directly from stored data for each file
 do f = 1, nFiles-1
    if (f > 2) then
        lol = 1
    else
        lol = 0
    end if
    q_prim_vf(f+lol)%sf(i, j, 0) = stored_values(i+1, f)
end do

! Set element 3 to zero (as requested)
    q_prim_vf(3)%sf(i, j, 0) = 0.0_wp
    ! print *, y_cc(0)
     q_prim_vf(2)%sf(i,j,0)=q_prim_vf(2)%sf(i,j,0)+(0.1_wp*0.00001*366.11_wp)*sin(2.0_wp*pi*(y_cc(j)*6.0_wp/(0.0292888_wp/2.0_wp)))


   case(212)

      delta_omega = 1.0_wp
      YH2_F   = 0.05_wp          ! Hydrogen in fuel stream
      YO2_F   = 0.0_wp           ! Oxygen in fuel stream  
      YH2O_F  = 0.0_wp           ! Water vapor in fuel stream
      YN2_F   = 0.95_wp          ! Nitrogen in fuel stream (1.0 - 0.05 - 0.0 - 0.0)
      YH_F    = 0.0_wp           ! H radical in fuel stream
      YO_F    = 0.0_wp           ! O radical in fuel stream
      YOH_F   = 0.0_wp           ! OH radical in fuel stream
      YHO2_F  = 0.0_wp           ! HO2 radical in fuel stream
      YH2O2_F = 0.0_wp           ! H2O2 in fuel stream
      YAR_F = 0.0_wp           ! H2O2 in fuel stream
      ! Oxidizer stream (O) mass fractions  
      YH2_O   = 0.0_wp           ! Hydrogen in oxidizer stream
      YO2_O   = 0.278_wp         ! Oxygen in oxidizer stream
      YH2O_O  = 0.17_wp          ! Water vapor in oxidizer stream
      YN2_O   = 0.552_wp         ! Nitrogen in oxidizer stream (1.0 - 0.278 - 0.17 - radicals)
      YH_O    = 5.60e-7_wp       ! H radical in oxidizer stream
      YO_O    = 1.55e-4_wp       ! O radical in oxidizer stream
      YOH_O   = 1.83e-3_wp       ! OH radical in oxidizer stream
      YHO2_O  = 5.10e-6_wp       ! HO2 radical in oxidizer stream
      YH2O2_O = 2.50e-7_wp       ! H2O2 in oxidizer stream
      YAR_F = 0.0_wp           ! H2O2 in fuel stream

    if (y_cc(j) .ge. x_cc(i)*tan(beta)-60.0_wp) then

      q_prim_vf(1)%sf(i,j,0) = 0.5_wp * ((rho_F + rho_O) + &
                                        (rho_F - rho_O) * &
                                        tanh(2.0_wp * y_cc(j) / delta_omega))
        

        ! x-velocity component (u1) with hyperbolic tangent profile
        q_prim_vf(2)%sf(i,j,0) = 0.5_wp * ((U_fuel + U_oxidizer) + &
                                           (U_fuel - U_oxidizer) * &
                                           tanh(2.0_wp * y_cc(j) / delta_omega))
        
        ! Alternative simplified form:
        ! q_prim_vf(2)%sf(0,j,0) = 1303.5_wp - 330.5_wp * tanh(2.0_wp * y_cc(j) / delta_omega)
        
        ! y-velocity component (u2) - initially zero
        q_prim_vf(3)%sf(i,j,0) = 0.0_wp
        
        ! z-velocity component (u3) - initially zero  
        q_prim_vf(4)%sf(i,j,0) = 94232.25

        q_prim_vf(6)%sf(i,j,0) = 0.5_wp * ((YH2_F + YH2_O) + &
                                        (YH2_F - YH2_O) * &
                                        tanh(2.0_wp * y_cc(j) / delta_omega))
        
        ! Y_{H} (index 7)  
          q_prim_vf(7)%sf(i,j,0) = 0.5_wp * ((YH_F + YH_O) + &
                                          (YH_F - YH_O) * &
                                          tanh(2.0_wp * y_cc(j) / delta_omega))
        
        ! Y_{O} (index 8)
        q_prim_vf(8)%sf(i,j,0) = 0.5_wp * ((YO_F + YO_O) + &
                                          (YO_F - YO_O) * &
                                          tanh(2.0_wp * y_cc(j) / delta_omega))
        
        ! Y_{O2} (index 9)
        q_prim_vf(9)%sf(i,j,0) = 0.5_wp * ((YO2_F + YO2_O) + &
                                          (YO2_F - YO2_O) * &
                                          tanh(2.0_wp * y_cc(j) / delta_omega))
        
        ! Y_{OH} (index 10)
        q_prim_vf(10)%sf(i,j,0) = 0.5_wp * ((YOH_F + YOH_O) + &
                                            (YOH_F - YOH_O) * &
                                            tanh(2.0_wp * y_cc(j) / delta_omega))
        
        ! Y_{H2O} (index 11)
        q_prim_vf(11)%sf(i,j,0) = 0.5_wp * ((YH2O_F + YH2O_O) + &
                                            (YH2O_F - YH2O_O) * &
                                            tanh(2.0_wp * y_cc(j) / delta_omega))
        
        ! Y_{HO2} (index 12)
        q_prim_vf(12)%sf(i,j,0) = 0.5_wp * ((YHO2_F + YHO2_O) + &
                                            (YHO2_F - YHO2_O) * &
                                            tanh(2.0_wp * y_cc(j) / delta_omega))
        
        ! Y_{H2O2} (index 13)
        q_prim_vf(13)%sf(i,j,0) = 0.5_wp * ((YH2O2_F + YH2O2_O) + &
                                            (YH2O2_F - YH2O2_O) * &
                                            tanh(2.0_wp * y_cc(j) / delta_omega))
        
        ! Y_{AR} (index 14) - No Argon
        q_prim_vf(14)%sf(i,j,0) = 0.0_wp
        
        ! Y_{N2} (index 15) - Calculate as remainder to ensure sum = 1
        q_prim_vf(15)%sf(i,j,0) = 1.0_wp - (q_prim_vf(6)%sf(0,j,0) + &   ! H2
                                            q_prim_vf(7)%sf(0,j,0) + &   ! H
                                            q_prim_vf(8)%sf(0,j,0) + &   ! O
                                            q_prim_vf(9)%sf(0,j,0) + &   ! O2
                                            q_prim_vf(10)%sf(0,j,0) + &  ! OH
                                            q_prim_vf(11)%sf(0,j,0) + &  ! H2O
                                            q_prim_vf(12)%sf(0,j,0) + &  ! HO2
                                            q_prim_vf(13)%sf(0,j,0) + &  ! H2O2
                                            q_prim_vf(14)%sf(0,j,0))     ! AR

        
            
      else 
      q_prim_vf(1)%sf(i,j,0)=0.260771_wp
        q_prim_vf(2)%sf(i,j,0) = 1526.4_wp  ! u1 (streamwise velocity)
        q_prim_vf(3)%sf(i,j,0) = 165.7_wp 
       q_prim_vf(4)%sf(i,j,0) = 129951.0_wp
       q_prim_vf(6)%sf(i,j,0) = YH2_O      ! H2
    q_prim_vf(7)%sf(i,j,0) = YH_O       ! H
    q_prim_vf(8)%sf(i,j,0) = YO_O       ! O
    q_prim_vf(9)%sf(i,j,0) = YO2_O      ! O2
    q_prim_vf(10)%sf(i,j,0) = YOH_O     ! OH
    q_prim_vf(11)%sf(i,j,0) = YH2O_O    ! H2O
    q_prim_vf(12)%sf(i,j,0) = YHO2_O    ! HO2
    q_prim_vf(13)%sf(i,j,0) = YH2O2_O   ! H2O2
    q_prim_vf(14)%sf(i,j,0) = YAR_O     ! AR
    
    ! Y_{N2} - Calculate as remainder
    q_prim_vf(15)%sf(i,j,0) = 1.0_wp - (YH2_O + YH_O + YO_O + YO2_O + &
                                        YOH_O + YH2O_O + YHO2_O + YH2O2_O + YAR_O)
        ! For other i locations, you might want to copy the profile or set differently

    endif

        q_prim_vf(5)%sf(i,j,0) = 1.0_wp  ! volume fraction is always 1
        case default
        if (proc_rank == 0) then
            call s_int_to_str(patch_id, iStr)
            call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
        end if

    end select

#:enddef
