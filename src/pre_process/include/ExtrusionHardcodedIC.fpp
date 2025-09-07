!===============================================================================
! Read 3D prim.* files: x y z value  (k fastest, then j, then i)
! Detect Nx, Ny, Nz; assign all variables (no extrusion)
!===============================================================================

#:def HardcodedDimensionsExtrusion()
    integer :: xRows, yRows, zRows, nRows
    integer :: iix, iiy, iiz, max_files
    integer :: f, iter, ios, ios2, unit, unit2, idx, idy, idz
    integer :: index_x, index_y, index_z, jump, line_count, ycount
    real(wp) :: x_len, x_step, y_len, y_step, z_step, z_len, delta_z
    real(wp) :: dummy_x, dummy_y, dummy_z, x0, y0
    integer :: global_offset_x, global_offset_y, global_offset_z
    real(wp) :: delta_x, delta_y
    character(len=100), dimension(sys_size) :: fileNames
    character(len=200) :: errmsg
    real(wp), allocatable :: stored_values(:, :, :, :)          ! (x,y,z,var)
    real(wp), allocatable :: x_coords(:), y_coords(:), z_coords(:)
    logical :: files_loaded = .false.
    real(wp) :: domain_xstart, domain_xend, domain_ystart, domain_yend
    character(len=*), parameter :: init_dir = "/u/dadam/MFC-Adam/examples/2D_Diffusion_Instabillities/" ! For example /home/MFC/examples/1D_Shock/D/
    character(len=20) :: file_num_str
    character(len=20) :: zeros_part
    character(len=6),  parameter :: zeros_default = "000000"
#:enddef

#:def HardcodedReadValues()

        if (.not. files_loaded) then
        max_files = sys_size

        ! --- Try prim.1 with a few timestep paddings; keep the first that opens ---
        write (file_num_str, '(I0)') 1

        zeros_part = zeros_default
        open (newunit=unit2, file=trim(init_dir)//"prim."//trim(file_num_str)//".00."//trim(zeros_part)//".dat", &
              status='old', action='read', form='formatted', access='sequential', iostat=ios2, iomsg=errmsg)

        if (ios2 /= 0) then
            zeros_part = "000001"
            open (newunit=unit2, file=trim(init_dir)//"prim."//trim(file_num_str)//".00."//trim(zeros_part)//".dat", &
                  status='old', action='read', form='formatted', access='sequential', iostat=ios2, iomsg=errmsg)
        end if
        if (ios2 /= 0) then
            zeros_part = "000000"
            open (newunit=unit2, file=trim(init_dir)//"prim."//trim(file_num_str)//".00."//trim(zeros_part)//".dat", &
                  status='old', action='read', form='formatted', access='sequential', iostat=ios2, iomsg=errmsg)
        end if

        if (ios2 /= 0) then
            call s_mpi_abort( "OPEN failed for: "// &
                trim(init_dir)//"prim."//trim(file_num_str)//".00."//trim(zeros_part)//".dat ; iomsg="//trim(errmsg) )
        end if

        ! --- Build full list with the detected zeros_part ---
        do f = 1, max_files
            write (file_num_str, '(I0)') f
            fileNames(f) = trim(init_dir)//"prim."//trim(file_num_str)//".00."//trim(zeros_part)//".dat"
        end do

        select case (num_dims)
        case (1, 2)
            ! Count lines
            line_count = 0
            do
                read (unit2, *, iostat=ios2) dummy_x, dummy_y
                if (ios2 /= 0) exit
                line_count = line_count + 1
            end do
            close (unit2)

            xRows = line_count
            yRows = 1
            zRows = 1
            index_x = 0
            if (num_dims == 2) index_x = i
            @:ALLOCATE (x_coords(xRows), stored_values(xRows, 1, 1, sys_size))

            ! Read data from all files
            do f = 1, max_files
                open (newunit=unit, file=trim(fileNames(f)), status='old', action='read', iostat=ios)
                if (ios /= 0) call s_mpi_abort("Error opening file: "//trim(fileNames(f)))
                do iter = 1, xRows
                    read (unit, *, iostat=ios) x_coords(iter), stored_values(iter, 1, 1, f)
                    if (ios /= 0) call s_mpi_abort("Error reading file: "//trim(fileNames(f)))
                end do
                close (unit)
            end do

            ! Offsets
            domain_xstart = x_coords(1)
            x_step = x_cc(1) - x_cc(0)
            delta_x = merge(x_cc(0) - domain_xstart + x_step/2.0_wp, &
                            x_cc(index_x) - domain_xstart + x_step/2.0_wp, num_dims == 1)
            global_offset_x = nint(abs(delta_x)/x_step)

        case (3)
            ! 3D: read true 3D files (x y z value), z fastest then y then x
            xRows = 128
            yRows = 256
            zRows = 64
            nRows = xRows * yRows * zRows

            @:ALLOCATE (x_coords(nRows), y_coords(nRows), z_coords(nRows))
            @:ALLOCATE (stored_values(xRows, yRows, zRows, sys_size))
            index_x = i
            index_y = j
            index_z = k

            do f = 1, max_files
                open (newunit=unit, file=trim(fileNames(f)), status='old', action='read', iostat=ios)
                if (ios /= 0) then
                    if (f == 1) call s_mpi_abort("Error opening file: "//trim(fileNames(f)))
                    cycle
                end if

                iter = 0
                do iix = 1, xRows
                    do iiy = 1, yRows
                        do iiz = 1, zRows
                            iter = iter + 1
                            if (f == 1) then
                                read (unit, *, iostat=ios) x_coords(iter), y_coords(iter), z_coords(iter), stored_values(iix, iiy, iiz, f)
                            else
                                read (unit, *, iostat=ios) dummy_x, dummy_y, dummy_z, stored_values(iix, iiy, iiz, f)
                            end if
                            if (ios /= 0) call s_mpi_abort("Error reading data")
                        end do
                    end do
                end do
                close (unit)
            end do

            ! Steps from solver grid (assumed equal to file spacing)
            x_step = x_cc(1) - x_cc(0)
            y_step = y_cc(1) - y_cc(0)
            z_step = z_cc(1) - z_cc(0)

            ! Offsets from first file point (no 0.5 shift)
            delta_x = x_cc(index_x) - x_coords(1)
            delta_y = y_cc(index_y) - y_coords(1)
            delta_z = z_cc(index_z) - z_coords(1)

            global_offset_x = nint(delta_x / x_step)
            global_offset_y = nint(delta_y / y_step)
            global_offset_z = nint(delta_z / z_step)
        end select

        files_loaded = .true.
    end if

    ! Data assignment
    select case (num_dims)
    case (1)
        idx = i + 1 + global_offset_x
        do f = 1, sys_size
            q_prim_vf(f)%sf(i, 0, 0) = stored_values(idx, 1, 1, f)
        end do

    case (2)
        idx = i + 1 + global_offset_x - index_x
        do f = 1, sys_size - 1
            jump = merge(1, 0, f >= momxe)
            q_prim_vf(f + jump)%sf(i, j, 0) = stored_values(idx, 1, 1, f)
        end do
        q_prim_vf(momxe)%sf(i, j, 0) = 0.0_wp

    case (3)
        idx = i + 1 + global_offset_x - index_x
        idy = j + 1 + global_offset_y - index_y
        idz = k + 1 + global_offset_z - index_z
        do f = 1, sys_size
            q_prim_vf(f)%sf(i, j, k) = stored_values(idx, idy, idz, f)
        end do
    end select
#:enddef

#:def HardcodedDellacation()
    if (allocated(stored_values)) then
        @:DEALLOCATE (stored_values)
    end if
    if (allocated(x_coords)) then
        @:DEALLOCATE (x_coords)
    end if
    if (allocated(y_coords)) then
        @:DEALLOCATE (y_coords)
    end if
    if (allocated(z_coords)) then
        @:DEALLOCATE (z_coords)
    end if
#:enddef

