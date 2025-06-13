#:def Hardcoded2DVariables()

    real(wp) :: eps,P_ref
    real(wp) :: r, rmax, gam, umax, p0
    real(wp) :: rhoH, rhoL, pRef, pInt, h, lam, wl, amp, intH, intL, alph
    real(wp) :: factor
    integer  :: Nx1,Nx2,Ny1,Ny2,Ny3,Ny4,Ny5,Ny6,mpi_size,ierr,global_idx
    integer  :: mpi_rank,lol,idx,stdout,i_min,i_max,local_x_min,local_x_max,k
integer :: local_start, local_end, local_n
    real(wp) :: YN2_O,YN2_F,YH2_F,YH2_O,YO2_O,YO2_F,YH2O_F,YH2O_O, YH_F,YH_O,YO_F,YO_O
    real(wp) ::  YOH_F,YOH_O,YHO2_F,YHO2_O,YH2O2_F,YH2O2_O,YAR_O,YAR_F, temp, sum_INVMW
      real(wp), parameter :: U_fuel = 724.0741_wp      ! m/s
      real(wp), parameter :: U_oxidizer = 2257.1586_wp
       real(wp), parameter :: delta_omega = 0.0001_wp
      real(wp), parameter :: beta_deg = 33.0_wp ! Angle of Oblique
      real(wp), parameter :: beta     = beta_deg * acos(-1.0_wp)/180.0_wp ! Angle in rads
      real(wp), parameter :: temp_fu =  300.0_wp
      real(wp), parameter :: temp_ox = 1500.0_wp
      real(wp), parameter :: rho_F = 0.609939148553946_wp   ! fuel‐stream density (H₂/N₂)
      real(wp), parameter :: rho_O = 0.2343940698879096_wp   ! oxidizer‐stream density 
  integer, parameter :: nFiles = 14
  integer, parameter :: xRows  = 1281
  integer, parameter :: Nrows  = xRows

  ! Variables for file reading and domain data
  real(wp)                ::  domain_start, domain_end, x_step
  real(wp), dimension(Nrows)            :: x_coords
  real(wp), dimension(Nrows, nFiles)      :: stored_values
  character(len=100), dimension(nFiles)   :: fileNames
  character(len=20)       :: file_num_str
  character(len=6), parameter :: zeros_default = "000000"
  character(len=*), parameter :: init_dir = "/scratch/bbsc/dadam/MFC-Adam/examples/1D_reactive_shocktube/"
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
         do f = 1, nFiles
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
            if (x_cc(i) .ge.  0.010378) then
             q_prim_vf(2)%sf(i,j,0)=q_prim_vf(2)%sf(i,j,0)+(0.1_wp*94.67*10**(-6))*sin(2.0_wp*pi*(y_cc(j)*6.0_wp/(0.015147_wp/2.0_wp)))
          end if
    case (250) ! MHD Orszag-Tang vortex
        ! gamma = 5/3
        !   rho = 25/(36*pi)
        !     p = 5/(12*pi)
        !     v = (-sin(2*pi*y), sin(2*pi*x), 0)
        !     B = (-sin(2*pi*y)/sqrt(4*pi), sin(4*pi*x)/sqrt(4*pi), 0)

        q_prim_vf(momxb)%sf(i, j, 0) = -sin(2._wp*pi*y_cc(j))
        q_prim_vf(momxb + 1)%sf(i, j, 0) = sin(2._wp*pi*x_cc(i))

        q_prim_vf(B_idx%beg)%sf(i, j, 0) = -sin(2._wp*pi*y_cc(j))/sqrt(4._wp*pi)
        q_prim_vf(B_idx%beg + 1)%sf(i, j, 0) = sin(4._wp*pi*x_cc(i))/sqrt(4._wp*pi)

    case (251) ! RMHD Cylindrical Blast Wave [Mignone, 2006: Section 4.3.1]

        if (x_cc(i)**2 + y_cc(j)**2 < 0.08_wp**2) then
            q_prim_vf(contxb)%sf(i, j, 0) = 0.01
            q_prim_vf(E_idx)%sf(i, j, 0) = 1.0
        elseif (x_cc(i)**2 + y_cc(j)**2 <= 1._wp**2) then
            ! Linear interpolation between r=0.08 and r=1.0
            factor = (1.0_wp - sqrt(x_cc(i)**2 + y_cc(j)**2))/(1.0_wp - 0.08_wp)
            q_prim_vf(contxb)%sf(i, j, 0) = 0.01_wp*factor + 1.e-4_wp*(1.0_wp - factor)
            q_prim_vf(E_idx)%sf(i, j, 0) = 1.0_wp*factor + 3.e-5_wp*(1.0_wp - factor)
        else
            q_prim_vf(contxb)%sf(i, j, 0) = 1.e-4_wp
            q_prim_vf(E_idx)%sf(i, j, 0) = 3.e-5_wp
        end if
  case (270)
        YH2_F   =  0.067133_wp          ! Hydrogen in fuel stream
      YO2_F   = 0.0_wp           ! Oxygen in fuel stream  
      YH2O_F  = 0.0_wp           ! Water vapor in fuel stream
      YN2_F   = 0.932867_wp          ! Nitrogen in fuel stream (1.0 - 0.05 - 0.0 - 0.0)
      YH_F    = 0.0_wp           ! H radical in fuel stream
      YO_F    = 0.0_wp           ! O radical in fuel stream
      YOH_F   = 0.0_wp           ! OH radical in fuel stream
      YHO2_F  = 0.0_wp           ! HO2 radical in fuel stream
      YH2O2_F = 0.0_wp           ! H2O2 in fuel stream
      ! Oxidizer stream (O) mass fractions  
      YH2_O   = 0.0_wp           ! Hydrogen in oxidizer stream
      YO2_O   = 0.232909_wp         ! Oxygen in oxidizer stream
      YH2O_O  = 0.0_wp          ! Water vapor in oxidizer stream
      YN2_O   =  0.767091_wp         ! Nitrogen in oxidizer stream (1.0 - 0.278 - 0.17 - radicals)
      YH_O    = 0.0_wp       ! H radical in oxidizer stream
      YO_O    = 0.0_wp       ! O radical in oxidizer stream
      YOH_O   = 0.0_wp       ! OH radical in oxidizer stream
      YHO2_O  = 0.0_wp       ! HO2 radical in oxidizer stream
      YH2O2_O = 0.0_wp       ! H2O2 in oxidizer stream
       P_ref=101325.0_wp

        

        ! x-velocity component (u1) with hyperbolic tangent profile
        q_prim_vf(2)%sf(i,j,0) = 0.5_wp * ((u_fuel + u_oxidizer) + &
                                           (u_oxidizer - u_fuel) * &
                                           tanh(2.0_wp * y_cc(j) / delta_omega))
        
        ! Alternative simplified form:
        ! q_prim_vf(2)%sf(0,j,0) = 1303.5_wp - 330.5_wp * tanh(2.0_wp * y_cc(j) / delta_omega)
        
        ! y-velocity component (u2) - initially zero
        q_prim_vf(3)%sf(i,j,0) = 0.0_wp
        
        ! z-velocity component (u3) - initially zero  
        q_prim_vf(4)%sf(i,j,0) = P_ref

        ! Volume fraction advidx always 1
        q_prim_vf(5)%sf(i,j,0) = 1.0_wp


         ! Y_{H2} (index 6)  
        q_prim_vf(6)%sf(i,j,0) = 0.5_wp * ((YH2_F + YH2_O) - &
                                        (YH2_F - YH2_O) * &
                                        tanh(2.0_wp * y_cc(j) / delta_omega))
        
        ! Y_{H} (index 7)  
          q_prim_vf(7)%sf(i,j,0) = 0.5_wp * ((YH_F + YH_O) - &
                                          (YH_F - YH_O) * &
                                          tanh(2.0_wp * y_cc(j) / delta_omega))
        
        ! Y_{O2} (index 8)
        q_prim_vf(8)%sf(i,j,0) = 0.5_wp * ((YO2_F + YO2_O) - &
                                          (YO2_F - YO2_O) * &
                                          tanh(2.0_wp * y_cc(j) / delta_omega))
        ! Y_{OH} (index 9)
        q_prim_vf(9)%sf(i,j,0) = 0.5_wp * ((YOH_F + YOH_O) - &
                                            (YOH_F - YOH_O) * &
                                            tanh(2.0_wp * y_cc(j) / delta_omega))
         
        ! Y_{O} (index 10)
        q_prim_vf(10)%sf(i,j,0) = 0.5_wp * ((YO_F + YO_O) -&
                                          (YO_F - YO_O) * &
                                          tanh(2.0_wp * y_cc(j) / delta_omega))
        
       
        ! Y_{H2O} (index 11)
        q_prim_vf(11)%sf(i,j,0) = 0.5_wp * ((YH2O_F + YH2O_O) - &
                                            (YH2O_F - YH2O_O) * &
                                            tanh(2.0_wp * y_cc(j) / delta_omega))
        
        ! Y_{HO2} (index 12)
        q_prim_vf(12)%sf(i,j,0) = 0.5_wp * ((YHO2_F + YHO2_O) - &
                                            (YHO2_F - YHO2_O) * &
                                            tanh(2.0_wp * y_cc(j) / delta_omega))
        
        ! Y_{H2O2} (index 13)
        q_prim_vf(13)%sf(i,j,0) = 0.5_wp * ((YH2O2_F + YH2O2_O) - &
                                            (YH2O2_F - YH2O2_O) * &
                                            tanh(2.0_wp * y_cc(j) / delta_omega))
        
       
        ! Y_{N2} (index 14) - Calculate as remainder to ensure sum = 1
        q_prim_vf(14)%sf(i,j,0) = 1.0_wp - (q_prim_vf(6)%sf(0,j,0) + &   ! H2
                                            q_prim_vf(7)%sf(0,j,0) + &   ! H
                                            q_prim_vf(8)%sf(0,j,0) + &   ! O2
                                            q_prim_vf(9)%sf(0,j,0) + &   ! OH
                                            q_prim_vf(10)%sf(0,j,0) + &  ! O
                                            q_prim_vf(11)%sf(0,j,0) + &  ! H2O
                                            q_prim_vf(12)%sf(0,j,0) + &  ! HO2
                                            q_prim_vf(13)%sf(0,j,0) )   ! H2O2
 
      sum_INVMW = 0.0_wp                              
     do k = 1, num_species
      ! species slot = 5 + k  (since slot 6→k=1, 7→k=2, …, 14→k=9)
      sum_INVMW = sum_INVMW +  q_prim_vf(chemxb - 1 + k)%sf(i,j,0)/molecular_weights(k)
    end do
    temp = temp_fu + 0.5_wp * (temp_ox - temp_fu) * &
           (1.0_wp + tanh(2.0_wp * y_cc(j) / delta_omega))
   
    q_prim_vf(1)%sf(i,j,0)= P_ref/(gas_constant*temp)/sum_INVMW

        case default
        if (proc_rank == 0) then
            call s_int_to_str(patch_id, iStr)
            call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
        end if

    end select

#:enddef
