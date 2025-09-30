!!>
!! @file   m_chemistry.f90
!! @brief  Contains module m_chemistry
!! @author Henry Le Berre <hberre3@gatech.edu>

#:include 'macros.fpp'
#:include 'case.fpp'

module m_chemistry

    use m_thermochem, only: &
        num_species, molecular_weights, get_temperature, get_net_production_rates, &
        get_mole_fractions, get_species_binary_mass_diffusivities, &
        get_species_mass_diffusivities_mixavg, gas_constant, get_mixture_molecular_weight, &
        get_mixture_energy_mass, get_mixture_thermal_conductivity_mixavg, get_species_enthalpies_rt, &
        get_mixture_viscosity_mixavg, get_mixture_specific_heat_cp_mass

    use m_global_parameters

    implicit none

    type(int_bounds_info) :: isc1, isc2, isc3
    $:GPU_DECLARE(create='[isc1, isc2, isc3]')
    integer, dimension(3) :: offsets
    $:GPU_DECLARE(create='[offsets]')

contains

    subroutine compute_viscosity_and_inversion(T_L, Ys_L, T_R, Ys_R, Re_L, Re_R)

        $:GPU_ROUTINE(function_name='compute_viscosity_and_inversion',parallelism='[seq]', &
            & cray_inline=True)

        real(wp), intent(inout) :: T_L, T_R, Re_L, Re_R
        real(wp), dimension(num_species), intent(inout) :: Ys_R, Ys_L
        real(wp) :: muref, Tref, S

        muref = 1.716_wp*10.0_wp**(-5.0_wp)
        Tref = 273.15_wp
        S= 110.14_wp

        Re_L = 1.0_wp/(muref*(T_L/Tref)**(1.5_wp)*(Tref+S)/(T_L+S))
        Re_R = 1.0_wp/(muref*(T_R/Tref)**(1.5_wp)*(Tref+S)/(T_R+S))

    end subroutine compute_viscosity_and_inversion
    subroutine s_compute_q_T_sf(q_T_sf, q_cons_vf, bounds)

        ! Initialize the temperature field at the start of the simulation to
        ! reasonable values. Temperature is computed the regular way using the
        ! conservative variables.

        type(scalar_field), intent(inout) :: q_T_sf
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        type(int_bounds_info), dimension(1:3), intent(in) :: bounds

        integer :: x, y, z, eqn
        real(wp) :: energy
        real(wp), dimension(num_species) :: Ys

        do z = bounds(3)%beg, bounds(3)%end
            do y = bounds(2)%beg, bounds(2)%end
                do x = bounds(1)%beg, bounds(1)%end
                    $:GPU_LOOP(parallelism='[seq]')
                    do eqn = chemxb, chemxe
                        Ys(eqn - chemxb + 1) = &
                            q_cons_vf(eqn)%sf(x, y, z)/q_cons_vf(contxb)%sf(x, y, z)
                    end do

                    ! e = E - 1/2*|u|^2
                    ! cons. E_idx     = \rho E
                    ! cons. contxb    = \rho         (1-fluid model)
                    ! cons. momxb + i = \rho u_i
                    energy = q_cons_vf(E_idx)%sf(x, y, z)/q_cons_vf(contxb)%sf(x, y, z)
                    $:GPU_LOOP(parallelism='[seq]')
                    do eqn = momxb, momxe
                        energy = energy - &
                                 0.5_wp*(q_cons_vf(eqn)%sf(x, y, z)/q_cons_vf(contxb)%sf(x, y, z))**2._wp
                    end do

                    call get_temperature(energy, dflt_T_guess, Ys, .true., q_T_sf%sf(x, y, z))
                end do
            end do
        end do

    end subroutine s_compute_q_T_sf

    subroutine s_compute_T_from_primitives(q_T_sf, q_prim_vf, bounds)

        type(scalar_field), intent(inout) :: q_T_sf
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(int_bounds_info), dimension(1:3), intent(in) :: bounds

        integer :: x, y, z, i
        real(wp), dimension(num_species) :: Ys
        real(wp) :: mix_mol_weight

        do z = bounds(3)%beg, bounds(3)%end
            do y = bounds(2)%beg, bounds(2)%end
                do x = bounds(1)%beg, bounds(1)%end
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = chemxb, chemxe
                        Ys(i - chemxb + 1) = q_prim_vf(i)%sf(x, y, z)
                    end do

                    call get_mixture_molecular_weight(Ys, mix_mol_weight)
                    q_T_sf%sf(x, y, z) = q_prim_vf(E_idx)%sf(x, y, z)*mix_mol_weight/(gas_constant*q_prim_vf(1)%sf(x, y, z))
                end do
            end do
        end do

    end subroutine s_compute_T_from_primitives

    subroutine s_compute_chemistry_reaction_flux(rhs_vf, q_cons_qp, q_T_sf, q_prim_qp, bounds)

        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(scalar_field), intent(inout) :: q_T_sf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_qp, q_prim_qp
        type(int_bounds_info), dimension(1:3), intent(in) :: bounds

        integer :: x, y, z
        integer :: eqn
        real(wp) :: T
        real(wp) :: rho, omega_m
        real(wp), dimension(num_species) :: Ys
        real(wp), dimension(num_species) :: omega

        $:GPU_PARALLEL_LOOP(collapse=3, private='[Ys, omega]')
        do z = bounds(3)%beg, bounds(3)%end
            do y = bounds(2)%beg, bounds(2)%end
                do x = bounds(1)%beg, bounds(1)%end

                    $:GPU_LOOP(parallelism='[seq]')
                    do eqn = chemxb, chemxe
                        Ys(eqn - chemxb + 1) = q_prim_qp(eqn)%sf(x, y, z)
                    end do

                    rho = q_cons_qp(contxe)%sf(x, y, z)
                    T = q_T_sf%sf(x, y, z)

                    call get_net_production_rates(rho, T, Ys, omega)

                    $:GPU_LOOP(parallelism='[seq]')
                    do eqn = chemxb, chemxe

                        omega_m = molecular_weights(eqn - chemxb + 1)*omega(eqn - chemxb + 1)

                        rhs_vf(eqn)%sf(x, y, z) = rhs_vf(eqn)%sf(x, y, z) + omega_m

                    end do

                end do
            end do
        end do

    end subroutine s_compute_chemistry_reaction_flux
    subroutine s_compute_chemistry_diffusion_flux(idir, q_prim_qp, flux_src_vf, irx, iry, irz)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_qp
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_src_vf
        type(int_bounds_info), intent(in) :: irx, iry, irz

        integer, intent(in) :: idir

        ! >>> Declarations simplified for a standard 2nd-order scheme
        real(wp), dimension(num_species) :: Ys_L, Ys_R, Ys_cell
        real(wp) :: mass_diff_L, mass_diff_R, mass_diff_Cell
        real(wp), dimension(num_species) :: dYk_dxi, h_l, h_r, h_k
        real(wp), dimension(num_species) :: Mass_Diffu_Flux
        real(wp) :: Mass_Diffu_Energy
        real(wp) :: MW_L, MW_R, Rgas_L, Rgas_R
        real(wp) :: T_L, T_R, P_L, P_R, rho_L, rho_R, rho_cell, rho_Vic
        real(wp) :: lambda_L, lambda_R, lambda_Cell, dT_dxi
        real(wp) :: mu_L, mu_R
        real(wp) :: Cp_L, Cp_R
        real(wp) :: muref, Tref, S, Prandtl, grid_spacing

        integer :: x, y, z, i, n, eqn
        integer, dimension(3) :: offsets
        ! >>> END Simplified declarations

        ! >>> User-defined physical model constants
        muref = 1.716_wp*10.0_wp**(-5.0_wp)
        Tref = 273.15_wp
        S = 110.14_wp
        Prandtl = 0.739_wp

        isc1 = irx; isc2 = iry; isc3 = irz
        $:GPU_UPDATE(device='[isc1,isc2,isc3]')

        if (chemistry) then
            offsets = 0
            offsets(idir) = 1

            $:GPU_PARALLEL_LOOP(collapse=3,  private='[Ys_L, Ys_R, Ys_cell, &
            & h_l, h_r, h_k, &
            & dYk_dxi, Mass_Diffu_Flux]', copyin='[offsets]')
            do z = isc3%beg, isc3%end
                do y = isc2%beg, isc2%end
                    do x = isc1%beg, isc1%end
                        select case (idir)
                        case (1); grid_spacing = x_cc(x + 1) - x_cc(x)
                        case (2); grid_spacing = y_cc(y + 1) - y_cc(y)
                        case (3); grid_spacing = z_cc(z + 1) - z_cc(z)
                        end select

                        ! >>> Extract Ys and primitive variables for a 2-point stencil (L and R)
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = chemxb, chemxe
                            Ys_L(i - chemxb + 1)  = q_prim_qp(i)%sf(x, y, z)
                            Ys_R(i - chemxb + 1)  = q_prim_qp(i)%sf(x + offsets(1), y + offsets(2), z + offsets(3))
                        end do

                        P_L  = q_prim_qp(E_idx)%sf(x, y, z)
                        P_R  = q_prim_qp(E_idx)%sf(x + offsets(1), y + offsets(2), z + offsets(3))

                        rho_L  = q_prim_qp(1)%sf(x, y, z)
                        rho_R  = q_prim_qp(1)%sf(x + offsets(1), y + offsets(2), z + offsets(3))

                        ! >>> Calculate thermodynamic properties for the 2 stencil points
                        call get_mixture_molecular_weight(Ys_L, MW_L); Rgas_L = gas_constant/MW_L; T_L = P_L/(rho_L*Rgas_L)
                        call get_mixture_molecular_weight(Ys_R, MW_R); Rgas_R = gas_constant/MW_R; T_R = P_R/(rho_R*Rgas_R)

                        ! >>> Calculate transport properties based on your model for the 2 points
                        mu_L = muref * (T_L / Tref)**(1.5_wp) * (Tref + S) / (T_L + S)
                        mu_R = muref * (T_R / Tref)**(1.5_wp) * (Tref + S) / (T_R + S)

                        call get_mixture_specific_heat_cp_mass(T_L, Ys_L, Cp_L); lambda_L = mu_L * Cp_L / Prandtl
                        call get_mixture_specific_heat_cp_mass(T_R, Ys_R, Cp_R); lambda_R = mu_R * Cp_R / Prandtl

                        mass_diff_L = mu_L / (rho_L * Prandtl)
                        mass_diff_R = mu_R / (rho_R * Prandtl)

                        ! >>> Get enthalpies and convert to mass basis for the 2 points
                        call get_species_enthalpies_rt(T_L, h_l)
                        call get_species_enthalpies_rt(T_R, h_r)
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = chemxb, chemxe
                           h_l(i-chemxb+1) = h_l(i-chemxb+1) * gas_constant * T_L / molecular_weights(i-chemxb+1)
                           h_r(i-chemxb+1) = h_r(i-chemxb+1) * gas_constant * T_R / molecular_weights(i-chemxb+1)
                        end do

                        ! >>> Using robust 2nd-order averages for all cell face properties
                        rho_cell       = 0.5_wp * (rho_L + rho_R)
                        lambda_Cell    = 0.5_wp * (lambda_L + lambda_R)
                        mass_diff_Cell = 0.5_wp * (mass_diff_L + mass_diff_R)

                        ! >>> Using robust 2nd-order derivatives for all gradients
                        dT_dxi = (T_R - T_L) / grid_spacing

                        Mass_Diffu_Energy = 0.0_wp
                        rho_Vic = 0.0_wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do eqn = chemxb, chemxe
                            Ys_cell(eqn-chemxb+1) = 0.5_wp * (Ys_L(eqn-chemxb+1) + Ys_R(eqn-chemxb+1))
                            h_k(eqn-chemxb+1) = 0.5_wp * (h_l(eqn-chemxb+1) + h_r(eqn-chemxb+1))
                            dYk_dxi(eqn-chemxb+1) = (Ys_R(eqn-chemxb+1) - Ys_L(eqn-chemxb+1)) / grid_spacing

                            ! >>> RESTORED: Sign removed as requested by user.
                            Mass_Diffu_Flux(eqn - chemxb + 1) = rho_cell * mass_diff_Cell * dYk_dxi(eqn - chemxb + 1)
                            rho_Vic = rho_Vic + Mass_Diffu_Flux(eqn - chemxb + 1)
                        end do

                        ! >>> This loop correctly applies the mass conservation correction and calculates enthalpy transport
                        $:GPU_LOOP(parallelism='[seq]')
                        do eqn = chemxb, chemxe
                            Mass_Diffu_Flux(eqn - chemxb + 1) = Mass_Diffu_Flux(eqn - chemxb + 1) - Ys_cell(eqn-chemxb+1) * rho_Vic
                            Mass_Diffu_Energy = Mass_Diffu_Energy + h_k(eqn - chemxb + 1) * Mass_Diffu_Flux(eqn - chemxb + 1)
                        end do

                        ! >>> RESTORED: Sign removed and original formulation restored as requested by user.
                        Mass_Diffu_Energy = lambda_Cell*dT_dxi + Mass_Diffu_Energy

                        ! >>> Update the right-hand-side flux arrays. The subtraction is kept as requested.
                        ! This implies the flux term is for the right face of the cell (i+1/2).
                        flux_src_vf(E_idx)%sf(x, y, z) = flux_src_vf(E_idx)%sf(x, y, z) - Mass_Diffu_Energy
                        $:GPU_LOOP(parallelism='[seq]')
                        do eqn = chemxb, chemxe
                            flux_src_vf(eqn)%sf(x, y, z) = flux_src_vf(eqn)%sf(x, y, z) - Mass_diffu_Flux(eqn - chemxb + 1)
                        end do
                    end do
                end do
            end do
        end if

    end subroutine s_compute_chemistry_diffusion_flux 
  end module m_chemistry
