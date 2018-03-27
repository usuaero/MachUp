!gfortran -fdefault-real-8 math.f90 dataset.f90 airfoil.f90 section.f90 wing.f90 plane.f90 special_functions.f90 view.f90 main.f90 -o AF3D.out
program main
    use view_m
    use special_functions_m
    use loads_m
    use plane_m
    implicit none
    type(plane_t) :: myplane
    type(json_value),pointer :: json_run, json_command
    character(100) :: filename
    integer :: i,nrun_types,dorun

    real :: time1,time2

    call cpu_time(time1)
    write(*,*) '-----------------------------------------------'
    write(*,*) '|                                             |'
    write(*,*) '|                         BBBBB               |'
    write(*,*) '|                       BB   BB               |'
    write(*,*) '|                     BB     BB               |'
    write(*,*) '|                    BB      BB               |'
    write(*,*) '|                  BB        BB               |'
    write(*,*) '|                 BB         BB               |'
    write(*,*) '|               BB           BB               |'
    write(*,*) '|              BB        BBBBBB               |'
    write(*,*) '|                                             |'
    write(*,*) '|                 MachUp 1.0                  |'
    write(*,*) '|                                             |'
    write(*,*) '|        (c) USU Aero Lab, LLC, 2016          |'
    write(*,*) '|                                             |'
    write(*,*) '|          This software comes with           |'
    write(*,*) '| ABSOLUTELY NO WARRANTY EXPRESSED OR IMPLIED |'
    write(*,*) '|                                             |'
    write(*,*) '|           Submit bug reports to:            |'
    write(*,*) '|           doug.hunsaker@usu.edu             |'
    write(*,*) '-----------------------------------------------'
    write(*,*)
    write(*,*)

    call get_command_argument(1, filename)
    myplane%master_filename = filename

    call plane_set_defaults(myplane)
    call plane_load_json(myplane)
    call plane_init_setup(myplane)
    myplane%verbose = 1

    call myplane%json%get('run', json_run)
    nrun_types = json_value_count(json_run)
!    write(*,*) 'Number of commands to run : ',nrun_types

    do i=1,nrun_types
        call json_value_get(json_run,i,json_command)
        run_type = trim(json_command%name)
        dorun = 1
        call json_get(json_command,'run',dorun,json_found)
        if(json_failed()) dorun = 1; !automatically run command if .run sub-command doesn't exist
        call json_clear_exceptions()

        if(dorun .eq. 1) then
            write(*,*)
            write(*,*) 'Running command : ',run_type

            myplane%verbose = 1

            select case (run_type)
                case ('stl')
                    call view_stl(myplane)
                case ('panair')
                    call view_panair(myplane, json_command)
                case ('plot')
                    call view_plotmtv(myplane)
                case ('forces')
                    call plane_run_current(myplane)
    !               call plane_distributed_loads(myplane)
    !               call view_plotmtv(myplane)
                case ('distributions')
                    call plane_run_current(myplane)
                    call sf_distributions(myplane,json_command)
!                case ('save')
!                    myplane%verbose = 1
!                    call plane_save(myplane)

                !special functions
                case ('derivatives')
                    myplane%verbose = 0
                    call sf_start_json_file(myplane,json_command)
                    call sf_derivs_stability(myplane)
                    call sf_derivs_control(myplane)
                    call sf_derivs_damping(myplane)
                    call sf_end_json_file(myplane)
                case ('damping')
                    myplane%verbose = 0
                    call sf_start_json_file(myplane,json_command)
                    call sf_derivs_damping(myplane)
                    call sf_end_json_file(myplane)
                case ('aerocenter')
                    myplane%verbose = 0
                    call sf_aerocenter(myplane,json_command)
                case ('stallonset')
                    myplane%verbose = 0
                    call sf_stallonset(myplane,json_command)
                case ('pointloads')
                    myplane%verbose = 1
                    call plane_run_current(myplane)
                    call loads_point(myplane,json_command)
                case ('spanloads')
                    myplane%verbose = 1
                    call plane_run_current(myplane)
                    call loads_per_span(myplane,json_command)
                case ('pitchtrim')
                    myplane%verbose = 0
                    call sf_pitch_trim(myplane,json_command)
                case ('targetcl')
                    myplane%verbose = 0
                    call sf_target_CL(myplane,json_command)
                case ('target')
                    call sf_target(myplane)
                case ('report')
                    call sf_report(myplane,json_command)
                case default
                    write(*,*) 'Command not recognized.'
            end select
        end if
    end do


    call plane_deallocate(myplane)

    call cpu_time(time2)
    write(*,*) 'CPU time total (sec): ',time2-time1

!    filename = trim(filename)//'_out.json'
!    call plane_write_json_file(myplane,filename)

end program main
!---notes---
!you have to connect to a wing with a lower id than yourself. So for loads calcs, you can sum
!by starting with the last wing and working backwards, adding the forces and moments
!to the parent wing
!To Do---------------------------------------------------------------
!apply Mach effects to airfoil data read from a 2D coeff file rather than the slope info.
!add flaps
