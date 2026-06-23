program calculate_mean_phi
    implicit none
    
    ! Use double precision for accurate accumulation
    integer, parameter :: dp = kind(1.0d0) 
    integer :: iostatus, in_unit, out_unit
    integer :: count
    real(dp) :: time_val, phi_val, sum_phi, mean_phi
    character(len=256) :: line

    ! Initialize variables
    count = 0
    sum_phi = 0.0_dp

    ! Open the input file (Assuming your data is saved in 'phi_data.txt')
    open(newunit=in_unit, file='ReadData/surfaceFieldValue.dat', status='old', action='read', iostat=iostatus)
    if (iostatus /= 0) then
        print *, "Error: Could not open 'phi_data.txt'. Please check the file name."
        stop
    end if

    ! Loop through the file line by line
    do
        ! Read the entire line as a string first
        read(in_unit, '(A)', iostat=iostatus) line
        
        ! Check for End of File (iostatus < 0)
        if (iostatus < 0) exit 
        if (iostatus > 0) then
            print *, "Error encountered while reading the file."
            stop
        end if

        ! Remove leading whitespace
        line = adjustl(line)
        
        ! Skip empty lines or header lines starting with '#'
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle

        ! Parse the time and phi values from the valid data lines
        read(line, *, iostat=iostatus) time_val, phi_val
        if (iostatus == 0) then
            sum_phi = sum_phi + phi_val
            count = count + 1
        end if
    end do

    close(in_unit)

    ! Calculate the mean and write to output file
    if (count > 0) then
        mean_phi = sum_phi / real(count, dp)
        
        print *, "Successfully processed ", count, " data points."
        print *, "Calculated Mean sum(phi) = ", mean_phi

        ! Open output file and write the result
        open(newunit=out_unit, file='phi_ref.txt', status='replace', action='write', iostat=iostatus)
        if (iostatus == 0) then
            write(out_unit, '(ES16.8)') mean_phi  ! Write using scientific notation
            close(out_unit)
            print *, "Mean successfully saved to 'phi_ref.txt'."
        else
            print *, "Error: Could not create 'phi_ref.txt' for writing."
        end if
    else
        print *, "Warning: No valid numeric data found in the file."
    end if

end program calculate_mean_phi