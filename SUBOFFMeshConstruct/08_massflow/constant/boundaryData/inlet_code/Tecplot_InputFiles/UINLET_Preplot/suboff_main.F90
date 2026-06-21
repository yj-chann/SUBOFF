program generate_ps1
    implicit none
    integer :: i
    integer :: max_files = 99
    integer, parameter :: unit_out = 10
    character(len=100) :: filename1, filename2

    ! Open the output PowerShell script file
    open(unit=unit_out, file='run_preplot.ps1', status='replace', action='write')

    write(unit_out, '(A)') 'Write-Host "Starting preplot batch job..."'
    write(unit_out, '(A)') ''

    ! Loop through the file indices (0 to 18)
    do i = 0, max_files
        ! Construct the filename string dynamically
        write(filename1, '(A,I0,A)') '..\UINLET\UINLET_Animate_', i, '.plt'
        write(filename2, '(A,I0,A)') '.\UINLET_Animate_', i, '.plt'
        ! Write the preplot command to the .ps1 file
        ! Syntax assumption: preplot <input_file> <output_file>
        write(unit_out, '(A,A,A,A)') 'preplot ', trim(filename1), ' ',filename2
    end do

    ! Add a completion message to the PowerShell script
    write(unit_out, '(A)') ''
    write(unit_out, '(A)') 'Write-Host "Preplotting complete!"'

    ! Close the file unit
    close(unit_out)

    print *, 'Successfully generated run_preplot.ps1'

end program generate_ps1