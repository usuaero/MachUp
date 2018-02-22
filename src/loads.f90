module loads_m
#ifdef dnad
    use dnadmod
#define real type(dual)
#endif
    use plane_m
    implicit none

contains

!-----------------------------------------------------------------------------------------------------------
subroutine loads_point(t,json_command)
    type(plane_t) :: t
    type(json_value),intent(in),pointer :: json_command
    character(len=:),allocatable :: cval
    character(100) :: filename
    integer :: iwing,isec,ierror
    type(section_t),pointer :: si
    real :: P(3),percent,chord
    120 Format(A15, 100ES25.13)


    !Get filename if specified
    call json_get(json_command,'filename', cval,json_found);
    if(json_failed() .or. (trim(cval).eq.'')) then !No filename specified
        call json_clear_exceptions()
        filename = trim(adjustl(t%master_filename))//'_pointloads.txt'
    else
        filename = trim(cval)
    end if

    open(unit = 10, File = filename, action = 'write', iostat = ierror)
    write(10,*) 'WingName            ControlPoint(x)          ControlPoint(y)          ControlPoint(z)          &
           &SpanwiseCoordinate       Chord                    &
           &F(x)/q                   F(y)/q                   F(z)/q                   &
           &M(x)/q                   M(y)/q                   M(z)/q'

    do iwing=1,t%nwings
        do isec=1,t%wings(iwing)%nSec
            si => t%wings(iwing)%sec(isec)
            P(:) = si%PC(:)
            percent = si%percent_c
            chord = si%chord_c
            write(10,120) t%wings(iwing)%name,P(:),percent*t%wings(iwing)%span,chord,si%F(:),si%M(:)
        end do
        write(10,*)
    end do
    close(10)
    write(*,'(A)') 'Point load results written to: '//trim(filename)

end subroutine loads_point

!-----------------------------------------------------------------------------------------------------------
subroutine loads_per_span(t,json_command)
    type(plane_t) :: t
    type(json_value),intent(in),pointer :: json_command
    character(len=:),allocatable :: cval
    character(100) :: filename
    integer :: iwing,isec,ierror
    type(section_t),pointer :: si
    type(dataset_t) :: dist
    real :: ans(7),P(3),percent,chord
    real :: zero
    120 Format(A15, 100ES25.13)

    zero = 0.0

    !Get filename if specified
    call json_get(json_command,'filename', cval,json_found);
    if(json_failed() .or. (trim(cval).eq.'')) then !No filename specified
        call json_clear_exceptions()
        filename = trim(adjustl(t%master_filename))//'_spanloads.txt'
    else
        filename = trim(cval)
    end if

    open(unit = 10, File = filename, action = 'write', iostat = ierror)
    write(10,*) 'WingName            X                        Y                        Z                        &
           &SpanwiseCoordinate       Chord                    &
           &CF(x)                    CF(y)                    CF(z)                    &
           &CM(x)                    CM(y)                    CM(z)'

    do iwing=1,t%nwings
        call loads_setup_cubic(t%wings(iwing),dist)

        !Load at root
        call ds_cubic_interpolate(dist, zero, 0, ans(:))

        si => t%wings(iwing)%sec(1)
        P(:) = si%P1(:)
        percent = si%percent_1
        chord = si%chord_1
        if(t%wings(iwing)%side.eq.'left') then
            P(:) = si%P2(:)
            percent = si%percent_2
            chord = si%chord_2
        end if

        write(10,120) t%wings(iwing)%name, P(:), percent*t%wings(iwing)%span,chord,ans(2:4)/chord,ans(5:7)/chord**2

        !Solve for all other cell vertex points
        do isec=1,t%wings(iwing)%nSec
            si => t%wings(iwing)%sec(isec)
            P(:) = si%P2(:)
            percent = si%percent_2
            chord = si%chord_2
            if(t%wings(iwing)%side.eq.'left') then
                P(:) = si%P1(:)
                percent = si%percent_1
                chord = si%chord_1
            end if

            call ds_cubic_interpolate(dist,percent*t%wings(iwing)%span,0,ans(:))

            write(10,120) t%wings(iwing)%name, P(:), percent*t%wings(iwing)%span,chord,ans(2:4)/chord,ans(5:7)/chord**2
        end do
        call ds_deallocate(dist)
        write(10,*)
    end do
    close(10)
    write(*,'(A)') 'Span load results written to: '//trim(filename)

end subroutine loads_per_span

!-----------------------------------------------------------------------------------------------------------
subroutine loads_setup_cubic(t,dist)
    type(wing_t) :: t
    type(dataset_t) :: dist
    type(section_t),pointer :: si
    integer :: isec
    real :: secSpan,rawdata(t%nSec,7)
    real :: zero

    zero = 0.0

    do isec=1,t%nSec
        si => t%sec(isec)
        secSpan = abs(si%percent_2-si%percent_1)*t%span
        rawdata(isec,1) = si%percent_c*t%span
        rawdata(isec,2:4) = si%F(:)/secSpan
        rawdata(isec,5:7) = si%M(:)/secSpan
    end do

    call ds_create_from_data(dist,t%nSec,7,rawdata)
    call ds_cubic_setup(dist, 1, 2, zero, 2, zero)

end subroutine loads_setup_cubic

!-----------------------------------------------------------------------------------------------------------

end module loads_m
