!>
!! @file m_perturbation.fpp
!! @brief Contains module m_perturbation

!> @brief This module contains subroutines that compute perturbations to the
!!              initial mean flow fields.
module m_perturbation

    use m_derived_types         ! Definitions of the derived types

    use m_global_parameters     ! Global parameters for the code

    use m_mpi_proxy              !< Message passing interface (MPI) module proxy

    use m_boundary_common   ! Boundary conditions module

    use m_helper

    use ieee_arithmetic

    implicit none

    real(wp), allocatable, dimension(:, :, :, :) :: q_prim_temp

contains

    impure subroutine s_initialize_perturbation_module()

        if (elliptic_smoothing) then
            allocate (q_prim_temp(0:m, 0:n, 0:p, 1:sys_size))
        end if

    end subroutine s_initialize_perturbation_module

    impure subroutine s_perturb_sphere(q_prim_vf)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer :: i, j, k, l !< generic loop operators

        real(wp) :: perturb_alpha

        real(wp) :: rand_real
        call random_seed()

        do k = 0, p
            do j = 0, n
                do i = 0, m
                    call random_number(rand_real)

                    perturb_alpha = q_prim_vf(E_idx + perturb_sph_fluid)%sf(i, j, k)

                    ! Perturb partial density fields to match perturbed volume fraction fields
                    !    IF ((perturb_alpha >= 25e-2_wp) .AND. (perturb_alpha <= 75e-2_wp)) THEN
                    if ((.not. f_approx_equal(perturb_alpha, 0._wp)) .and. (.not. f_approx_equal(perturb_alpha, 1._wp))) then

                        ! Derive new partial densities
                        do l = 1, num_fluids
                            q_prim_vf(l)%sf(i, j, k) = q_prim_vf(E_idx + l)%sf(i, j, k)*fluid_rho(l)
                        end do

                    end if
                end do
            end do
        end do

    end subroutine s_perturb_sphere

  impure subroutine E_vonK_sub(k, Ek_out)
    real(wp), intent(in)  :: k

    real(wp), intent(out) :: Ek_out
    real(wp)              :: kappa, Lambda, qq

    Lambda = 0.96_wp/1000.0_wp/3.0_wp

    qq = 3.0_wp/2.0_wp*25.0_wp

    Ek_out = 2.0_wp/3.0_wp*Lambda*(1.606_wp*(k*Lambda)**(4.0_wp))/(1.359+(k*Lambda)**(2.0_wp))**(17.0_wp/6.0_wp)*qq
  end subroutine E_vonK_sub

    impure subroutine s_perturb_surrounding_flow(q_prim_vf)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer :: i, j, k !<  generic loop iterators

        real(wp) :: perturb_alpha
        real(wp) :: rand_real
        call random_seed()

        ! Perturb partial density or velocity of surrounding flow by some random small amount of noise
        
    end subroutine s_perturb_surrounding_flow

    impure subroutine s_elliptic_smoothing(q_prim_vf, bc_type)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(integer_field), dimension(1:num_dims, -1:1), intent(in) :: bc_type
        integer :: i, j, k, l, q

        do q = 1, elliptic_smoothing_iters

            ! Communication of buffer regions and apply boundary conditions
            call s_populate_variables_buffers(bc_type, q_prim_vf, pb%sf, mv%sf)

            ! Perform smoothing and store in temp array
            if (n == 0) then
                do j = 0, m
                    do i = 1, sys_size
                        q_prim_temp(j, 0, 0, i) = (1._wp/4._wp)* &
                                                  (q_prim_vf(i)%sf(j + 1, 0, 0) + q_prim_vf(i)%sf(j - 1, 0, 0) + &
                                                   2._wp*q_prim_vf(i)%sf(j, 0, 0))
                    end do
                end do
            else if (p == 0) then
                do k = 0, n
                    do j = 0, m
                        do i = 1, sys_size
                            q_prim_temp(j, k, 0, i) = (1._wp/8._wp)* &
                                                      (q_prim_vf(i)%sf(j + 1, k, 0) + q_prim_vf(i)%sf(j - 1, k, 0) + &
                                                       q_prim_vf(i)%sf(j, k + 1, 0) + q_prim_vf(i)%sf(j, k - 1, 0) + &
                                                       4._wp*q_prim_vf(i)%sf(j, k, 0))
                        end do
                    end do
                end do
            else
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do i = 1, sys_size
                                q_prim_temp(j, k, l, i) = (1._wp/12._wp)* &
                                                          (q_prim_vf(i)%sf(j + 1, k, l) + q_prim_vf(i)%sf(j - 1, k, l) + &
                                                           q_prim_vf(i)%sf(j, k + 1, l) + q_prim_vf(i)%sf(j, k - 1, l) + &
                                                           q_prim_vf(i)%sf(j, k, l + 1) + q_prim_vf(i)%sf(j, k, l - 1) + &
                                                           6._wp*q_prim_vf(i)%sf(j, k, l))
                            end do
                        end do
                    end do
                end do
            end if

            ! Copy smoothed data back to array of scalar fields
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        do i = 1, sys_size
                            q_prim_vf(i)%sf(j, k, l) = q_prim_temp(j, k, l, i)
                        end do
                    end do
                end do
            end do
        end do

    end subroutine s_elliptic_smoothing

    !>  This subroutine computes velocity perturbations for a temporal mixing
        !!              layer with a hyperbolic tangent mean streamwise velocity
        !!              profile, using an inverter version of the spectrum-based
        !!              synthetic turbulence generation method proposed by
        !!              Guo et al. (2023, JFM).
subroutine s_perturb_mixlayer(q_prim_vf)
    type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
    real(wp), dimension(100) :: k, Ek,dk
    real(wp), dimension(3, 3) :: Rij, Lmat
   real(wp), dimension(3) :: velfluc, sig_tmp, sig, khat, xi
    real(wp) ::  arg, Eint, qi, psi, phase,phi
    integer  :: i, j, r, ierr, seed,l
    real(wp) :: Lx, Ly, dx, dy, kmin, kmax

     mixlayer_perturb_nk = 100
    ! --- box lengths and spacings (2D) ---
    Lx = x_cc(m) - x_cc(0)
    Ly = y_cc(n) - y_cc(0)
    dx = x_cc(1) - x_cc(0)
    dy = y_cc(1) - y_cc(0)

    print *, dx
    ! --- spectral bounds (grid-safe) ---
     kmin = 2*pi/Lx
    kmax = 10*pi / min(dx, dy)
! Logarithmic spacing concentrates points at low k


    ! --- energy spectrum & its integral ---
    Eint = 0._wp
  do i = 1, mixlayer_perturb_nk
        k(i) = kmin * (kmax/kmin)**((i-1.0_wp)/(mixlayer_perturb_nk-1.0_wp))
    end do

    do i = 1, mixlayer_perturb_nk-1
        dk(i) = k(i+1)-k(i)
        print *, k(i), i
        call E_vonK_sub(k(i), Ek(i))
         print *, Ek(i)
        Eint = Eint + Ek(i)*dk(i)
    end do

print *, Eint
   

    ! --- no anisotropy mapping for now ---
    Rij = 0._wp
    Lmat = 0._wp
    Lmat(1,1) = 1._wp;  Lmat(2,2) = 1._wp;  Lmat(3,3) = 1._wp


    ! --- build perturbations (2D: r loops over y) ---
    do r = 0, n

        do i = 1, mixlayer_perturb_nk-1
                if (proc_rank == 0) then
                    call s_generate_random_perturbation(khat, xi, phi, i, y_cc(r))
                end if

#ifdef MFC_MPI
                call MPI_BCAST(khat, 3, mpi_p, 0, MPI_COMM_WORLD, ierr)
                call MPI_BCAST(xi, 3, mpi_p, 0, MPI_COMM_WORLD, ierr)
                call MPI_BCAST(phi, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
#endif

                ! Compute mode direction by two-time cross product
                sig_tmp = f_cross(xi, khat)
                sig_tmp = sig_tmp/sqrt(sum(sig_tmp**2._wp))
                sig = f_cross(khat, sig_tmp)
                
            ! correct per-mode amplitude: qi = sqrt(E(k_i) * dk)
            qi = sqrt( max(Ek(i)*dk(i), 0._wp) )

            ! add perturbation (2D: l=0 only)
                           do l = 0, p
            do j = 0, m
                        arg = k(i)*(khat(1)*x_cc(j) + khat(2)*y_cc(r) + khat(3)*z_cc(l)) + 2._wp*pi*phi
                        velfluc = 2._wp*qi*sig*cos(arg)
                        velfluc = matmul(Lmat, velfluc)
                        q_prim_vf(momxb)%sf(j, r, l) = q_prim_vf(momxb)%sf(j, r, l) + velfluc(1)
                        q_prim_vf(momxb + 1)%sf(j, r, l) = q_prim_vf(momxb + 1)%sf(j, r, l) + velfluc(2)
                        q_prim_vf(momxb + 2)%sf(j, r, l) = q_prim_vf(momxb + 2)%sf(j, r, l) + velfluc(3)
            end do
            end do
        end do
    end do
    end subroutine s_perturb_mixlayer


    subroutine s_generate_random_perturbation(khat, xi, phi, ik, yloc)
        integer, intent(in) :: ik
        real(wp), intent(in) :: yloc
        real(wp), dimension(3), intent(out) :: khat, xi
        real(wp), intent(out) :: phi
        real(wp) :: theta, eta
        integer :: seed, kfac, yfac

        kfac = ik*amplifier
        yfac = nint((sin(yloc) + 1._wp)*amplifier)
        seed = nint(0.5_wp*modmul(kfac) + 0.5_wp*modmul(yfac))

        call s_prng(theta, seed)
        call s_prng(eta, seed)
        khat = f_unit_vector(theta, eta)

        call s_prng(theta, seed)
        call s_prng(eta, seed)
        xi = f_unit_vector(theta, eta)

        call s_prng(phi, seed)

    end subroutine s_generate_random_perturbation

    ! Generate a random unit vector (spherical distribution)
    function f_unit_vector(theta, eta) result(vec)
        real(wp), intent(in) :: theta, eta
        real(wp) :: zeta, xi
        real(wp), dimension(3) :: vec

        xi = 2._wp*pi*theta
        zeta = acos(2._wp*eta - 1._wp)
        vec(1) = sin(zeta)*cos(xi)
        vec(2) = sin(zeta)*sin(xi)
        vec(3) = cos(zeta)

    end function f_unit_vector

    !>  This function generates a pseudo-random number between 0 and 1 based on
    !!  linear congruential generator.
    subroutine s_prng(var, seed)
        integer, intent(inout) :: seed
        real(wp), intent(out) :: var
        integer :: i

        seed = mod(modmul(seed), modulus)
        var = seed/real(modulus, wp)

    end subroutine s_prng

    function modmul(a) result(val)
        integer, intent(in) :: a
        integer :: val
        real(wp) :: x, y

        x = (multiplier/real(modulus, wp))*a + (increment/real(modulus, wp))
        y = nint((x - floor(x))*decimal_trim)/decimal_trim
        val = nint(y*modulus)

    end function modmul

    impure subroutine s_finalize_perturbation_module()

        if (elliptic_smoothing) then
            deallocate (q_prim_temp)
        end if

    end subroutine s_finalize_perturbation_module

end module m_perturbation
