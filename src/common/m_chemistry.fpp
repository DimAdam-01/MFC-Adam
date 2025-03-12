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
        get_mixture_energy_mass,get_mixture_thermal_conductivity_mixavg,get_species_enthalpies_rt

    use m_global_parameters

    implicit none

contains

    subroutine s_compute_q_T_sf(q_T_sf, q_cons_vf, bounds)

        ! Initialize the temperature field at the start of the simulation to
        ! reasonable values. Temperature is computed the regular way using the
        ! conservative variables.

        type(scalar_field), intent(inout) :: q_T_sf
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        type(int_bounds_info), dimension(1:3), intent(in) :: bounds

        integer :: x, y, z, eqn
        real(wp) :: energy, mean_molecular_weight
        real(wp), dimension(num_species) :: Ys

        do z = bounds(3)%beg, bounds(3)%end
            do y = bounds(2)%beg, bounds(2)%end
                do x = bounds(1)%beg, bounds(1)%end
                    !$acc loop seq
                    do eqn = chemxb, chemxe
                        Ys(eqn - chemxb + 1) = &
                            q_cons_vf(eqn)%sf(x, y, z)/q_cons_vf(contxb)%sf(x, y, z)
                    end do

                    ! e = E - 1/2*|u|^2
                    ! cons. E_idx     = \rho E
                    ! cons. contxb    = \rho         (1-fluid model)
                    ! cons. momxb + i = \rho u_i
                    energy = q_cons_vf(E_idx)%sf(x, y, z)/q_cons_vf(contxb)%sf(x, y, z)
                    !$acc loop seq
                    do eqn = momxb, momxe
                        energy = energy - &
                                 0.5_wp*(q_cons_vf(eqn)%sf(x, y, z)/q_cons_vf(contxb)%sf(x, y, z))**2._wp
                    end do

                    call get_temperature(energy, dflt_T_guess, Ys, .true., q_T_sf%sf(x, y, z))
                    
                end do
             !   print *, q_cons_vf(E_idx)%sf(bounds(1)%end+1, 0, 0)
            end do
        end do

    end subroutine s_compute_q_T_sf

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

        !$acc parallel loop collapse(3) gang vector default(present) &
        !$acc private(Ys, omega)
        do z = bounds(3)%beg, bounds(3)%end
            do y = bounds(2)%beg, bounds(2)%end
                do x = bounds(1)%beg, bounds(1)%end

                    !$acc loop seq
                    do eqn = chemxb, chemxe
                        Ys(eqn - chemxb + 1) = q_prim_qp(eqn)%sf(x, y, z)
                    end do

                    rho = q_prim_qp(contxe)%sf(x, y, z)
                    T = q_T_sf%sf(x, y, z)

                    call get_net_production_rates(rho, T, Ys, omega)

                    !$acc loop seq
                    do eqn = chemxb, chemxe

                        omega_m = molecular_weights(eqn - chemxb + 1)*omega(eqn - chemxb + 1)

                        rhs_vf(eqn)%sf(x, y, z) = rhs_vf(eqn)%sf(x, y, z) + omega_m

                    end do

                end do
            end do
        end do

    end subroutine s_compute_chemistry_reaction_flux



  subroutine s_compute_chemistry_diffusion_flux(idir, q_prim_qp, flux_src_vf, q_T_sf, bounds)

    
      type (scalar_field),  dimension(sys_size), intent(in) ::  q_prim_qp
      type (scalar_field),  dimension(sys_size), intent(inout) :: flux_src_vf
      type(scalar_field), intent(inout) :: q_T_sf
      type(int_bounds_info), dimension(1:3), intent(in) :: bounds

      integer, intent(in) :: idir

      real(wp), dimension(num_species) :: Xs_L,Xs_R,Xs_cell,Ys_L,Ys_R,Ys_cell
      real(wp), dimension(num_species) :: mass_diffusivities_mixavg1, mass_diffusivities_mixavg2
      real(wp), dimension(num_species) :: mass_diffusivities_mixavg_Cell, dXk_dx_n,h_k,h_l,h_r
      real(wp), dimension(num_species) :: Mass_Diffu_Flux
      real(wp) :: Mass_Diffu_Energy
      real(wp) ::  MW_L,MW_R,MW_cell,Rgas_L,Rgas_R,T_L,T_R,P_L,P_R,rho_L,rho_R,rho_cell,rho_Vic,E_L,E_R
      real(wp) :: Tref,c1,b,lamda_L,lamda_R,Pr,mu_L,mu_R,Cp,Diffu,lamda_Cell,dT_dx_n,dx,dy,dz

      integer :: x,y,z,i,n,eqn


      if (chemistry) then

          if (idir == 1) then
              !$acc parallel loop collapse(3) gang vector default(present)
              do z = bounds(3)%beg, bounds(3)%end
                  do y = bounds(2)%beg, bounds(2)%end
                      do x = bounds(1)%beg, bounds(1)%end

                          dx=x_cc(x+1) - x_cc(x)

                          !$acc loop seq
                          do i = chemxb,chemxe
                              Ys_L(i-chemxb+1) = q_prim_qp(i)%sf(x,y,z)
                              Ys_R(i-chemxb+1) = q_prim_qp(i)%sf(x+1,y,z)
                          end do

                          Ys_cell=0.5_wp*(Ys_L+Ys_R)
                          

                          call get_mixture_molecular_weight(Ys_L,MW_L) 
                          call get_mixture_molecular_weight(Ys_R,MW_R) 

                          MW_cell=0.5_wp*(MW_L+MW_R)

                          call get_mole_fractions(MW_L,Ys_L,Xs_L)
                          call get_mole_fractions(MW_R,Ys_R,Xs_R)

                          Rgas_L = gas_constant/MW_L
                          Rgas_R = gas_constant/MW_R

                          P_L =q_prim_qp(E_idx)%sf(x,y,z)
                          P_R =q_prim_qp(E_idx)%sf(x+1,y,z)

                          rho_L =q_prim_qp(1)%sf(x,y,z)
                          rho_R =q_prim_qp(1)%sf(x+1,y,z)

                          T_L=  P_L/rho_L/Rgas_L
                          T_R = P_R/rho_R/Rgas_R

                          Xs_cell=0.5_wp*(Xs_L+Xs_R)


                          rho_cell=0.5_wp*(rho_L+rho_R)

                          dT_dx_n=(T_R-T_L)/dx

                        !  print *, rho_R,x

                          call get_species_enthalpies_rt(T_L,h_l)
                          call get_species_enthalpies_rt(T_R,h_r)

                          h_l=h_l*gas_constant*T_L/molecular_weights(:)
                          h_r=h_r*gas_constant*T_R/molecular_weights(:)

                          h_k=0.5_wp*(h_l+h_r)

                          call get_species_mass_diffusivities_mixavg(&
                            P_L, T_L, Ys_L, mass_diffusivities_mixavg1)
                          call get_species_mass_diffusivities_mixavg(&
                            P_R, T_R, Ys_R, mass_diffusivities_mixavg2)

                          call get_mixture_thermal_conductivity_mixavg(&
                            T_L, Ys_L,lamda_L)

                          call get_mixture_thermal_conductivity_mixavg(&
                            T_R,Ys_R,lamda_R )

                          !$acc loop seq
                          do i = chemxb,chemxe
                              mass_diffusivities_mixavg_Cell(i-chemxb+1)=(mass_diffusivities_mixavg2(i-chemxb+1) &
                                  +mass_diffusivities_mixavg1(i-chemxb+1))/2.0_wp*(1.0_wp-Xs_cell(i-chemxb+1))/(1.0_wp-Ys_cell(i-chemxb+1))
                          end do

                          dXk_dx_n=(Xs_R(:)-Xs_L(:))/dx

                          lamda_Cell=0.5_wp*(lamda_R+lamda_L)

                          dT_dx_n=(T_R-T_L)/dx

                          Mass_Diffu_Flux(:)=rho_cell*mass_diffusivities_mixavg_Cell(:)*molecular_weights(:)/MW_cell*dXk_dx_n(:)

                          rho_Vic=0.0_wp
                          Mass_Diffu_Energy=0.0_wp
                          
                          !$acc loop seq
                          do eqn=chemxb,chemxe
                              rho_Vic=rho_Vic+Mass_Diffu_Flux(eqn-chemxb+1)
                              Mass_Diffu_Energy=Mass_Diffu_Energy+h_k(eqn-chemxb+1)*Mass_Diffu_Flux(eqn-chemxb+1)
                          end do

                          !$acc loop seq
                          do eqn=chemxb,chemxe
                              Mass_Diffu_Energy=Mass_Diffu_Energy-h_k(eqn-chemxb+1)*Ys_cell(eqn-chemxb+1)*rho_Vic
                          end do

                          Mass_Diffu_Flux=Mass_Diffu_Flux -rho_Vic*Ys_cell
                          Mass_Diffu_Energy=lamda_Cell*dT_dx_n+Mass_Diffu_Energy


                          flux_src_vf(E_idx)%sf(x, y, z) = -Mass_Diffu_Energy

                          !$acc loop seq
                          do eqn=chemxb, chemxe
                              flux_src_vf(eqn)%sf(x,y,z)=0.0_wp
                              flux_src_vf(eqn)%sf(x,y,z)= & 
                              flux_src_vf(eqn)%sf(x,y,z) - &
                              Mass_diffu_Flux(eqn-chemxb+1)
                          end do

                         ! print *, q_prim_qp(E_idx)%sf(bounds(1)%end+1,y,z), bounds(1)%end+1

                      end do
                  end do
              end do

          elseif (idir == 2) then
              !$acc parallel loop collapse(3) gang vector default(present)
              do z = bounds(3)%beg, bounds(3)%end
                  do y = bounds(2)%beg, bounds(2)%end
                      do x = bounds(1)%beg, bounds(1)%end

                          dy=y_cc(y+1) - y_cc(y)

                          !$acc loop seq
                          do i = chemxb,chemxe
                              Ys_L(i-chemxb+1) = q_prim_qp(i)%sf(x,y,z)
                              Ys_R(i-chemxb+1) = q_prim_qp(i)%sf(x,y+1,z)
                          end do

                          Ys_cell=0.5_wp*(Ys_L+Ys_R)

                          call get_mixture_molecular_weight(Ys_L,MW_L) 
                          call get_mixture_molecular_weight(Ys_R,MW_R) 

                          MW_cell=0.5_wp*(MW_L+MW_R)

                          call get_mole_fractions(MW_L,Ys_L,Xs_L)
                          call get_mole_fractions(MW_R,Ys_R,Xs_R)

                          Rgas_L = gas_constant/MW_L
                          Rgas_R = gas_constant/MW_R

                          P_L =q_prim_qp(E_idx)%sf(x,y,z)
                          P_R =q_prim_qp(E_idx)%sf(x,y+1,z)

                          rho_L =q_prim_qp(1)%sf(x,y,z)
                          rho_R =q_prim_qp(1)%sf(x,y+1,z)

                          T_L=  P_L/rho_L/Rgas_L
                          T_R = P_R/rho_R/Rgas_R

                          Xs_cell=0.5_wp*(Xs_L+Xs_R)

                          rho_cell=0.5_wp*(rho_L+rho_R)

                          dT_dx_n=(T_R-T_L)/dy

                          call get_species_enthalpies_rt(T_L,h_l)
                          call get_species_enthalpies_rt(T_R,h_r)

                          h_l=h_l*gas_constant*T_L/molecular_weights(:)
                          h_r=h_r*gas_constant*T_R/molecular_weights(:)

                          h_k=0.5_wp*(h_l+h_r)

                          call get_species_mass_diffusivities_mixavg(&
                            P_L, T_L, Ys_L, mass_diffusivities_mixavg1)
                          call get_species_mass_diffusivities_mixavg(&
                            P_R, T_R, Ys_R, mass_diffusivities_mixavg2)

                          call get_mixture_thermal_conductivity_mixavg(&
                            T_L, Ys_L,lamda_L)

                          call get_mixture_thermal_conductivity_mixavg(&
                            T_R,Ys_R,lamda_R )

                          !$acc loop seq
                          do i = chemxb,chemxe
                              mass_diffusivities_mixavg_Cell(i-chemxb+1)=(mass_diffusivities_mixavg2(i-chemxb+1) &
                                  +mass_diffusivities_mixavg1(i-chemxb+1))/2.0_wp*(1.0_wp-Xs_cell(i-chemxb+1))/(1.0_wp-Ys_cell(i-chemxb+1))
                          end do

                          dXk_dx_n=(Xs_R(:)-Xs_L(:))/dx

                          lamda_Cell=0.5_wp*(lamda_R+lamda_L)

                          dT_dx_n=(T_R-T_L)/dy

                          Mass_Diffu_Flux(:)=rho_cell*mass_diffusivities_mixavg_Cell(:)*molecular_weights(:)/MW_cell*dXk_dx_n(:)

                          rho_Vic=0.0_wp
                          Mass_Diffu_Energy=0.0_wp

                          !$acc loop seq
                          do eqn=chemxb,chemxe
                              rho_Vic=rho_Vic+Mass_Diffu_Flux(eqn-chemxb+1)
                              Mass_Diffu_Energy=Mass_Diffu_Energy+h_k(eqn-chemxb+1)*Mass_Diffu_Flux(eqn-chemxb+1)
                          end do

                          !$acc loop seq
                          do eqn=chemxb,chemxe
                              Mass_Diffu_Energy=Mass_Diffu_Energy-h_k(eqn-chemxb+1)*Ys_cell(eqn-chemxb+1)*rho_Vic
                          end do

                          Mass_Diffu_Flux=Mass_Diffu_Flux -rho_Vic*Ys_cell
                          Mass_Diffu_Energy=lamda_Cell*dT_dx_n+Mass_Diffu_Energy


                          flux_src_vf(E_idx)%sf(x, y, z) = flux_src_vf(E_idx)%sf(x, y, z)-Mass_Diffu_Energy

                          !$acc loop seq
                          do eqn=chemxb, chemxe
                              flux_src_vf(eqn)%sf(x,y,z)=0.0_wp
                              flux_src_vf(eqn)%sf(x,y,z)= & 
                              flux_src_vf(eqn)%sf(x,y,z) - &
                              Mass_diffu_Flux(eqn-chemxb+1)
                          end do
                      end do
                  end do
              end do

          elseif (idir == 3) then
              !$acc parallel loop collapse(3) gang vector default(present)
              do z = bounds(3)%beg, bounds(3)%end
                  do y = bounds(2)%beg, bounds(2)%end
                      do x = bounds(1)%beg, bounds(1)%end

                          dz=z_cc(z+1) - z_cc(z)

                          !$acc loop seq
                          do i = chemxb,chemxe
                              Ys_L(i-chemxb+1) = q_prim_qp(i)%sf(x,y,z)
                              Ys_R(i-chemxb+1) = q_prim_qp(i)%sf(x,y,z+1)
                          end do

                          Ys_cell=0.5_wp*(Ys_L+Ys_R)

                          call get_mixture_molecular_weight(Ys_L,MW_L) 
                          call get_mixture_molecular_weight(Ys_R,MW_R) 

                          MW_cell=0.5_wp*(MW_L+MW_R)

                          call get_mole_fractions(MW_L,Ys_L,Xs_L)
                          call get_mole_fractions(MW_R,Ys_R,Xs_R)

                          Rgas_L = gas_constant/MW_L
                          Rgas_R = gas_constant/MW_R

                          P_L =q_prim_qp(E_idx)%sf(x,y,z)
                          P_R =q_prim_qp(E_idx)%sf(x,y,z+1)

                          rho_L =q_prim_qp(1)%sf(x,y,z)
                          rho_R =q_prim_qp(1)%sf(x,y,z+1)

                          T_L=  P_L/rho_L/Rgas_L
                          T_R = P_R/rho_R/Rgas_R

                          Xs_cell=0.5_wp*(Xs_L+Xs_R)

                          rho_cell=0.5_wp*(rho_L+rho_R)

                          dT_dx_n=(T_R-T_L)/dz

                          call get_species_enthalpies_rt(T_L,h_l)
                          call get_species_enthalpies_rt(T_R,h_r)

                          h_l=h_l*gas_constant*T_L/molecular_weights(:)
                          h_r=h_r*gas_constant*T_R/molecular_weights(:)

                          h_k=0.5_wp*(h_l+h_r)

                          call get_species_mass_diffusivities_mixavg(&
                            P_L, T_L, Ys_L, mass_diffusivities_mixavg1)
                          call get_species_mass_diffusivities_mixavg(&
                            P_R, T_R, Ys_R, mass_diffusivities_mixavg2)

                          call get_mixture_thermal_conductivity_mixavg(&
                            T_L, Ys_L,lamda_L)

                          call get_mixture_thermal_conductivity_mixavg(&
                            T_R,Ys_R,lamda_R )

                          !$acc loop seq
                          do i = chemxb,chemxe
                              mass_diffusivities_mixavg_Cell(i-chemxb+1)=(mass_diffusivities_mixavg2(i-chemxb+1) &
                                  +mass_diffusivities_mixavg1(i-chemxb+1))/2.0_wp*(1.0_wp-Xs_cell(i-chemxb+1))/(1.0_wp-Ys_cell(i-chemxb+1))
                          end do

                          dXk_dx_n=(Xs_R(:)-Xs_L(:))/dx

                          lamda_Cell=0.5_wp*(lamda_R+lamda_L)

                          dT_dx_n=(T_R-T_L)/dz

                          Mass_Diffu_Flux(:)=rho_cell*mass_diffusivities_mixavg_Cell(:)*molecular_weights(:)/MW_cell*dXk_dx_n(:)

                          rho_Vic=0.0_wp
                          Mass_Diffu_Energy=0.0_wp
                          
                          !$acc loop seq
                          do eqn=chemxb,chemxe
                              rho_Vic=rho_Vic+Mass_Diffu_Flux(eqn-chemxb+1)
                              Mass_Diffu_Energy=Mass_Diffu_Energy+h_k(eqn-chemxb+1)*Mass_Diffu_Flux(eqn-chemxb+1)
                          end do

                          !$acc loop seq
                          do eqn=chemxb,chemxe
                              Mass_Diffu_Energy=Mass_Diffu_Energy-h_k(eqn-chemxb+1)*Ys_cell(eqn-chemxb+1)*rho_Vic
                          end do

                          Mass_Diffu_Flux=Mass_Diffu_Flux -rho_Vic*Ys_cell
                          Mass_Diffu_Energy=lamda_Cell*dT_dx_n+Mass_Diffu_Energy


                          flux_src_vf(E_idx)%sf(x, y, z) = flux_src_vf(E_idx)%sf(x, y, z)-Mass_Diffu_Energy

                          !$acc loop seq
                          do eqn=chemxb, chemxe
                              flux_src_vf(eqn)%sf(x,y,z)=0.0_wp
                              flux_src_vf(eqn)%sf(x,y,z)= & 
                              flux_src_vf(eqn)%sf(x,y,z) - &
                              Mass_diffu_Flux(eqn-chemxb+1)
                          end do
                      end do
                  end do
              end do
          end if

      end if


  end subroutine

end module m_chemistry
