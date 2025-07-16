    !     Code for Anisotropies in the Microwave Background
    !     by Antony Lewis (http://cosmologist.info/) and Anthony Challinor
    !     This is a sample driver routine that reads
    !     in one set of parameters and produdes the corresponding output.

    program driver
    use CAMB
    implicit none
    character(len=:), allocatable :: InputFile

    integer :: start_time, end_time, clock_rate, clock_max
    real :: elapsed_time, start, finish

    InputFile = ''
    if (GetParamCount() /= 0)  InputFile = GetParam(1)
    if (InputFile == '') error stop 'No parameter input file'

    call system_clock(count_rate=clock_rate, count_max=clock_max)
    if (clock_rate == 0) then
      print *, "Error: System clock rate is zero. Cannot measure time."
      stop
    end if
    call system_clock(count=start_time)
    call cpu_time(start)

    call CAMB_CommandLineRun(InputFile)
    deallocate(InputFile) ! Just so no memory leaks in valgrind

    call cpu_time(finish)
    call system_clock(count=end_time)
    if (end_time < start_time) then
      write(*,*) 'end_time < start_time', end_time, start_time, clock_rate, clock_max
      elapsed_time = (real(clock_max - start_time) + real(end_time) + 1.0) / real(clock_rate)
    else
      elapsed_time = real(end_time - start_time) / real(clock_rate)
    end if
    write(*,'(a,1x,f12.5)') 'Total     Time taken :', elapsed_time
    write(*,'(a,1x,f12.5)') 'Total CPU Time taken :', finish-start


    end program driver

