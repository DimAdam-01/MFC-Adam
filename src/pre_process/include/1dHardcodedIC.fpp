#:def Hardcoded1DVariables()


    integer, parameter :: nFiles = 15   ! Can be changed to any number
    integer, parameter :: nRows  = 2001
    integer :: f, iter, ios, unit, idx
    real(8) :: x_len, x_step
    character(len=100), dimension(nFiles) :: fileNames
    ! Arrays to store all data from files - read once, use many times
    real(kind(0d0)), dimension(nRows, nFiles) :: stored_values
    real(kind(0d0)), dimension(nRows) :: x_coords
    logical :: files_loaded = .false.
    real(kind(0d0)) :: domain_start, domain_end
    
    character(len=*), parameter :: init_dir = "/home/pain/ChemMFC/examples/Initial/"
    character(len=20) :: file_num_str     ! For storing the file number as a string
    character(len=20) :: zeros_part       ! For the trailing zeros part
    character(len=6), parameter :: zeros_default = "000000"  ! Default zeros (can be changed)
    
    ! Generate file names dynamically in a loop
    do f = 1, nFiles
        ! Convert file number to string with proper formatting
        if (f < 10) then
            write(file_num_str, '(I1)') f  ! Single digit
        else
            write(file_num_str, '(I2)') f  ! Double digit
            ! For more than 99 files, you might need to adjust this format
        end if
        
        ! Create the filename with the pattern "prim.X.00.000000.dat"
        fileNames(f) = trim(init_dir) // "prim." // trim(file_num_str) // ".00." // zeros_default // ".dat"
    end do

#:enddef

#:def Hardcoded1D()

    select case (patch_icpp(patch_id)%hcid)
    case (100)
       q_prim_vf(momxb)%sf(i,0,0)=8*exp(-200*((x_cc(i)-0.004/2)/0.004)**2)
       q_prim_vf(E_idx)%sf(i,0,0)=1.01325*10**5+8.0d0*0.3d0*800.74d0*exp(-200*((x_cc(i)-0.004/2)/0.004)**2)
       ! q_prim_vf(1)%sf(i,0,0)=(1.01325*10**5+2.2d0*0.227d0*800.74d0*exp(-400*((x_cc(i)-0.002/2)/0.002)**2))/(8.314*1000*((0.4+x_cc(i)/0.002*0.2)/2.016+(0.6-x_cc(i)/0.002*0.2)/32) *(300+30*exp(-400*((x_cc(i)-0.002/2)/0.002)**2)))
      ! q_prim_vf(15)%sf(i,0,0)=300+50*exp(-300*((x_cc(i)-0.002/2)/0.002)**2)
    !  q_prim_vf(1)%sf(i,0,0)=(1.01325*10**5)/(8.314*1000*((0.4+x_cc(i)/0.002*0.2)/2.016+(0.6-x_cc(i)/0.002*0.2)/32) *(300+50*exp(-400*((x_cc(i)-0.002/2)/0.002)**2)))
       q_prim_vf(1)%sf(i,0,0)=0.08+0.08/800.74*(8*exp(-200*((x_cc(i)-0.004/2)/0.004)**2))
     !  q_prim_vf(5)%sf(i,0,0)=0.4+x_cc(i)/0.004*0.2
           q_prim_vf(7)%sf(i,0,0)=1.0
   !  q_prim_vf(8)%sf(i,0,0)=0.6-x_cc(i)/0.004*0.2
    ! q_prim_vf(10)%sf(i,0,0)=0.1
    !q_prim_vf(18)%sf(i,0,0)=0.1
    !q_prim_vf(19)%sf(i,0,0)=0.1
    !q_prim_vf(20)%sf(i,0,0)=0.1
    !q_prim_vf(52)%sf(i,0,0)=x_cc(i)/0.0054*0.4
    !q_prim_vf(5)%sf(i,0,0)=0.4
    ! q_prim_vf(9)%sf(i,0,0)=0.6-x_cc(i)/0.04*0.2
    !   q_prim_vf(9)%sf(i,0,0)=1.0d0
    !   q_prim_vf(15)%sf(i,0,0)=300.0d0
    !q_prim_vf(1)%sf(1,0,0)=(1.01325*10**5+8.0d0*0.227d0*800.74d0*exp(-400*((x_cc(i)-0.002/2)/0.002)**2))/(8314/17.008*300)

    ! ! (0.4+x_cc(i)/0.002*0.2)/2.016+(0.6-x_cc(i)/0.002*0.2)/32)!
    ! Put your variable assignments here
    ! (0.4+x_cc(i)/0.002*0.2)/2.016+(0.6-x_cc(i)/0.002*0.2)/32)


    case(101)
    q_prim_vf(momxb)%sf(i,0,0)=8*exp(-600*((x_cc(i)-0.0035/2)/0.002)**2)
       q_prim_vf(E_idx)%sf(i,0,0)=1.01325*10**5+2.2d0*0.227d0*818.74d0*exp(-600*((x_cc(i)-0.0035/2)/0.002)**2)
       q_prim_vf(1)%sf(i,0,0)=0.30!(1.01325*10**5+1100*exp(-400*((x_cc(i)-0.002/2)/0.002)**2))/(8.314*1000*((0.4+x_cc(i)/0.002*0.2)/2.016+(0.6-x_cc(i)/0.002*0.2)/32) *(300+30*exp(-400*((x_cc(i)-0.002/2)/0.002)**2)))
      ! q_prim_vf(1)%sf(i,0,0)=0.227+0.227/818.74*(8*exp(-400*((x_cc(i)-0.002/2)/0.002)**2))
      !
      !
      !
    case(102)
      

      !call execute_command_line("pwd", wait=.true.)

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
            q_prim_vf(f)%sf(i, 0, 0) = stored_values(idx+1, f)
        end do

    case default
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
    end select


#:enddef
