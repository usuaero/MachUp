module special_functions_m
#ifdef dnad
    use dnadmod
#define real type(dual)
#endif
    use plane_m
    implicit none

contains

!-----------------------------------------------------------------------------------------------------------
subroutine sf_distributions(t,json_command)
    type(plane_t) :: t
    type(json_value),intent(in),pointer :: json_command
    character(:), allocatable :: cval
    character(100) :: filename,output_type
    integer :: iwing,isec,ierror
    type(section_t),pointer :: si
    120 Format(A15, 100ES25.13)

    !Get filename if specified
    filename = trim(adjustl(t%master_filename))//'_distributions'
    call myjson_get(json_command, 'filename', cval, filename)
    filename = trim(adjustl(cval))

    !Get file type
    output_type = 'text'
    call myjson_get(json_command, 'output', cval, output_type);
    output_type = trim(adjustl(cval))

    if(trim(adjustl(output_type)).eq.'json') then
        filename = trim(adjustl(filename))//'.json'
        open(unit = 10, File = filename, action = 'write', iostat = ierror)
        write(10,*) '{'
        write(10,'(A84)',advance="no") '   "data" : "WingName             ControlPoint(x)          ControlPoint(y)          '
        write(10,'(A75)',advance="no") 'ControlPoint(z)          Chord                    Twist(deg)               '
        write(10,'(A75)',advance="no") 'Sweep(deg)               Dihedral(deg)            Area                     '
        write(10,'(A75)',advance="no") 'Alpha(deg)               Section_CL               Section_CD_parasitic     '
        write(10,'(A75)',advance="no") 'Section_Cm               CL(Ref)                  Section_alpha_L0(deg)    '
        do iwing=1,t%nwings
            write(10,'(A7)',advance="no") 'newline'
            do isec=1,t%wings(iwing)%nSec
                si => t%wings(iwing)%sec(isec)
                write(10,120,advance="no") t%wings(iwing)%name,si%PC(:),&
                              & si%chord_c,si%twist*180.0/pi,si%sweep*180.0/pi,si%dihedral*180.0/pi,&
                              & si%ds,&
                              & si%alpha*180.0/pi,sec_CL(si),sec_CD(si),sec_Cm(si),sec_CL(si)*si%chord_c/t%long_r,&
                              & sec_alpha_L0(si)*180.0/pi
                write(10,'(A7)',advance="no") 'newline'
            end do
        end do
        write(10,'(A)') '"'
        write(10,*) '}'
        close(10)
    else
        filename = trim(adjustl(filename))//'.txt'

        open(unit = 10, File = filename, action = 'write', iostat = ierror)
        write(10,'(A75)',advance="no") 'WingName                 ControlPoint(x)          ControlPoint(y)          '
        write(10,'(A75)',advance="no") 'ControlPoint(z)          Chord                    Twist(deg)               '
        write(10,'(A75)',advance="no") 'Sweep(deg)               Dihedral(deg)            Area                     '
        write(10,'(A75)',advance="no") 'Section_Alpha(deg)       Section_CL               Section_CD_parasitic     '
        write(10,'(A75)')              'Section_Cm               CL(Ref)                  Section_alpha_L0(deg)    '
        do iwing=1,t%nwings
            do isec=1,t%wings(iwing)%nSec
                si => t%wings(iwing)%sec(isec)
                write(10,120) t%wings(iwing)%name,si%PC(:),&
                              & si%chord_c,si%twist*180.0/pi,si%sweep*180.0/pi,si%dihedral*180.0/pi,&
                              & si%ds,&
                              & si%alpha*180.0/pi,sec_CL(si),sec_CD(si),sec_Cm(si),sec_CL(si)*si%chord_c/t%long_r,&
                              & sec_alpha_L0(si)*180.0/pi
            end do
!            write(10,*)
        end do
        close(10)
    end if

    write(*,'(A)') 'Distribution results written to: '//trim(filename)

end subroutine sf_distributions

!-----------------------------------------------------------------------------------------------------------
subroutine sf_derivs_stability(t)
    type(plane_t) :: t
    character(100) :: name
    character(1000) :: note
    real,allocatable,dimension(:,:,:) :: dF,dM !Global Forces and Moments
    real,allocatable,dimension(:,:) :: dL,dD,dLoD !Global Lift, Drag, and L/D
    real :: delta,alpha,beta
    real :: CLa,Cma

    allocate(dF(t%nrealwings+1,3,3))
    allocate(dM(t%nrealwings+1,3,3))
    allocate(dL(t%nrealwings+1,3))
    allocate(dD(t%nrealwings+1,3))
    allocate(dLoD(t%nrealwings+1,3))

    delta = 0.5*pi/180.0 !radians
    alpha = t%alpha
    beta = t%beta
    note = "Derivatives with respect to angle in radians."

    write(*,*)
    write(*,*) '---------- Stability Derivatives -----------'

    ! Deriv w/respect to alpha
    name = 'd_alpha'
    t%alpha = alpha + delta
    call plane_run_current(t)
    dF = t%GF; dM = t%GM; dL = t%GL; dD = t%GD; dLoD = t%GLoD;
    t%alpha = alpha - delta;
    call plane_run_current(t)
    dF = 0.5/delta*(dF - t%GF); dM = 0.5/delta*(dM - t%GM);
    dL = 0.5/delta*(dL - t%GL); dD = 0.5/delta*(dD - t%GD); dLoD = 0.5/delta*(dLoD - t%GLoD);
    call sf_append_json_file(t,name,dF,dM,dL,dD,dLoD,note);
    t%alpha = alpha;

    CLa = dL(t%nrealwings+1,3)
    Cma = dM(t%nrealwings+1,3,2)

    name = 'alpha'
    call sf_write_derivs_to_screen(dL(t%nrealwings+1,3),dD(t%nrealwings+1,3),dF(t%nrealwings+1,3,:),dM(t%nrealwings+1,3,:),name)

    ! Deriv w/respect to beta
    name = 'd_beta'
    t%beta = beta + delta
    call plane_run_current(t)
    dF = t%GF; dM = t%GM; dL = t%GL; dD = t%GD; dLoD = t%GLoD;
    t%beta = beta - delta;
    call plane_run_current(t)
    dF = 0.5/delta*(dF - t%GF); dM = 0.5/delta*(dM - t%GM);
    dL = 0.5/delta*(dL - t%GL); dD = 0.5/delta*(dD - t%GD); dLoD = 0.5/delta*(dLoD - t%GLoD);
    call sf_append_json_file(t,name,dF,dM,dL,dD,dLoD,note);
    t%beta = beta;

    name = 'beta'
    call sf_write_derivs_to_screen(dL(t%nrealwings+1,3),dD(t%nrealwings+1,3),dF(t%nrealwings+1,3,:),dM(t%nrealwings+1,3,:),name)

    write(*,*) 'Static Margin = ',-Cma/CLa
    write(*,*) '--------------------------------------------'
    write(*,*)

end subroutine sf_derivs_stability

!-----------------------------------------------------------------------------------------------------------
subroutine sf_derivs_control(t)
    type(plane_t) :: t
    character(100) :: name
    character(1000) :: note
    real,allocatable,dimension(:,:,:) :: dF,dM !Global Forces and Moments
    real,allocatable,dimension(:,:) :: dL,dD,dLoD !Global Lift, Drag, and L/D
    integer :: icontrol
    real :: delta, deflection

    allocate(dF(t%nrealwings+1,3,3))
    allocate(dM(t%nrealwings+1,3,3))
    allocate(dL(t%nrealwings+1,3))
    allocate(dD(t%nrealwings+1,3))
    allocate(dLoD(t%nrealwings+1,3))

    delta = 0.25*pi/180.0 !radians
    note = "Derivatives with respect to angle in radians."

    write(*,*) '------------ Control Derivatives ---------------'
    do icontrol = 1,t%ncontrols
        name = trim(t%controls(icontrol)%name)
        name = 'd_'//trim(name)
        deflection = t%controls(icontrol)%deflection
        t%controls(icontrol)%deflection = deflection + delta;
        call plane_run_current(t)
        dF = t%GF; dM = t%GM; dL = t%GL; dD = t%GD; dLoD = t%GLoD;
        t%controls(icontrol)%deflection = deflection - delta;
        call plane_run_current(t)
        dF = 0.5/delta*(dF - t%GF); dM = 0.5/delta*(dM - t%GM);
        dL = 0.5/delta*(dL - t%GL); dD = 0.5/delta*(dD - t%GD); dLoD = 0.5/delta*(dLoD - t%GLoD);
        call sf_append_json_file(t,name,dF,dM,dL,dD,dLoD,note);
        t%controls(icontrol)%deflection = deflection;
        name = trim(t%controls(icontrol)%name)
    call sf_write_derivs_to_screen(dL(t%nrealwings+1,3),dD(t%nrealwings+1,3),dF(t%nrealwings+1,3,:),dM(t%nrealwings+1,3,:),name)
    end do
    write(*,*) '------------------------------------------------'

end subroutine sf_derivs_control

!-----------------------------------------------------------------------------------------------------------
subroutine sf_derivs_damping(t)
    type(plane_t) :: t
    character(100) :: name
    character(1000) :: note
    real,allocatable,dimension(:,:,:) :: dF,dM !Global Forces and Moments
    real,allocatable,dimension(:,:) :: dL,dD,dLoD !Global Lift, Drag, and L/D
    real :: delta,omega_orig(3),lref

    allocate(dF(t%nrealwings+1,3,3))
    allocate(dM(t%nrealwings+1,3,3))
    allocate(dL(t%nrealwings+1,3))
    allocate(dD(t%nrealwings+1,3))
    allocate(dLoD(t%nrealwings+1,3))

    omega_orig(:) = t%omega(:)

    write(*,*)
    write(*,*) '----------- Damping Derivatives ------------'

    ! Deriv w/respect to roll
    lref = t%lat_r
    delta = 0.005*(2.0/lref) ! scaled by lateral_length/2.0
    note = "p_bar = p*lateral_length/(2*V_infinity)"
    name = 'd_p_bar'
    t%omega(1) = omega_orig(1) + delta
    call plane_run_current(t)
    dF = t%GF; dM = t%GM; dL = t%GL; dD = t%GD; dLoD = t%GLoD;
    t%omega(1) = omega_orig(1) - delta
    call plane_run_current(t)
    dF = 0.5/delta*(dF - t%GF)*2.0/lref; dM = 0.5/delta*(dM - t%GM)*2.0/lref;
    dL = 0.5/delta*(dL - t%GL)*2.0/lref; dD = 0.5/delta*(dD - t%GD)*2.0/lref; dLoD = 0.5/delta*(dLoD - t%GLoD)*2.0/lref;
    call sf_append_json_file(t,name,dF,dM,dL,dD,dLoD,note);
    t%omega(1) = omega_orig(1)

    name = 'pbar'
    call sf_write_derivs_to_screen(dL(t%nrealwings+1,3),dD(t%nrealwings+1,3),dF(t%nrealwings+1,3,:),dM(t%nrealwings+1,3,:),name)

    ! Deriv w/respect to pitch
    lref = t%long_r
    delta = 0.005*(2.0/lref) ! scaled by longitudinal_length/2.0
    note = "q_bar = q*longitudinal_length/(2*V_infinity)"
    name = 'd_q_bar'
    t%omega(2) = omega_orig(2) + delta
    call plane_run_current(t)
    dF = t%GF; dM = t%GM; dL = t%GL; dD = t%GD; dLoD = t%GLoD;
    t%omega(2) = omega_orig(2) - delta
    call plane_run_current(t)
    dF = 0.5/delta*(dF - t%GF)*2.0/lref; dM = 0.5/delta*(dM - t%GM)*2.0/lref;
    dL = 0.5/delta*(dL - t%GL)*2.0/lref; dD = 0.5/delta*(dD - t%GD)*2.0/lref; dLoD = 0.5/delta*(dLoD - t%GLoD)*2.0/lref;
    call sf_append_json_file(t,name,dF,dM,dL,dD,dLoD,note);
    t%omega(2) = omega_orig(2)

    name = 'qbar'
    call sf_write_derivs_to_screen(dL(t%nrealwings+1,3),dD(t%nrealwings+1,3),dF(t%nrealwings+1,3,:),dM(t%nrealwings+1,3,:),name)

    ! Deriv w/respect to yaw
    lref = t%lat_r
    delta = 0.005*(2.0/lref) ! scaled by lateral_length/2.0
    note = "r_bar = r*lateral_length/(2*V_infinity)"
    name = 'd_r_bar'
    t%omega(3) = omega_orig(3) + delta
    call plane_run_current(t)
    dF = t%GF; dM = t%GM; dL = t%GL; dD = t%GD; dLoD = t%GLoD;
    t%omega(3) = omega_orig(3) - delta
    call plane_run_current(t)
    dF = 0.5/delta*(dF - t%GF)*2.0/lref; dM = 0.5/delta*(dM - t%GM)*2.0/lref;
    dL = 0.5/delta*(dL - t%GL)*2.0/lref; dD = 0.5/delta*(dD - t%GD)*2.0/lref; dLoD = 0.5/delta*(dLoD - t%GLoD)*2.0/lref;
    call sf_append_json_file(t,name,dF,dM,dL,dD,dLoD,note);
    t%omega(3) = omega_orig(3)

    name = 'rbar'
    call sf_write_derivs_to_screen(dL(t%nrealwings+1,3),dD(t%nrealwings+1,3),dF(t%nrealwings+1,3,:),dM(t%nrealwings+1,3,:),name)

    write(*,*) '--------------------------------------------'
    write(*,*)

end subroutine sf_derivs_damping

!-----------------------------------------------------------------------------------------------------------
subroutine sf_write_derivs_to_screen(dL,dD,dF,dM,name)
    character(100) :: name
    real:: dL, dD, dF(3), dM(3)

    write(*,*) '    CL,',trim(name),' = ',dL
    write(*,*) '    CD,',trim(name),' = ',dD
    write(*,*) '    CY,',trim(name),' = ',dF(2)
    write(*,*) '    Cl,',trim(name),' = ',dM(1)
    write(*,*) '    Cm,',trim(name),' = ',dM(2)
    write(*,*) '    Cn,',trim(name),' = ',dM(3)
    write(*,*)

end subroutine sf_write_derivs_to_screen

!-----------------------------------------------------------------------------------------------------------
subroutine sf_start_json_file(t,json_command)
    type(plane_t) :: t
    type(json_value),intent(in),pointer :: json_command
    character(len=:),allocatable :: cval
    character(100) :: filename

    !Get filename if specified
    filename = trim(adjustl(t%master_filename))//'_derivatives.json'
    call myjson_get(json_command, 'filename', cval, filename);
    filename = trim(adjustl(cval))

    call json_value_create(t%p_json)           ! create the value and associate the pointer
    call to_object(t%p_json,trim(filename))    ! add the file name as the name of the overall structure

end subroutine sf_start_json_file

!-----------------------------------------------------------------------------------------------------------
subroutine sf_append_json_file(t,name,F,M,L,D,LoD,note)
    type(plane_t) :: t
    character(100) :: name
    character(1000) :: note
    real :: F(t%nrealwings+1,3,3),M(t%nrealwings+1,3,3),L(t%nrealwings+1,3),D(t%nrealwings+1,3),LoD(t%nrealwings+1,3)
    type(json_value),pointer    :: p_name, p_object
    integer :: iwing

    !Write JSON File
    call json_value_create(p_name)             !an object
    call to_object(p_name,trim(name))
    call json_value_add(t%p_json, p_name)
    call json_value_add(p_name, 'note',   trim(note))
    iwing = 0
    do iwing =1, t%nrealwings+1
        call json_value_create(p_object)
        if(iwing>t%nrealwings) then
            call to_object(p_object,trim(t%name))
        else
            call to_object(p_object,trim(t%wings(iwing)%name)//'_'//trim(t%wings(iwing)%side))
        end if
        call json_value_add(p_name,p_object)
        call plane_add_json_force(p_object,F(iwing,3,:),M(iwing,3,:),L(iwing,3),D(iwing,3),LoD(iwing,3))
        nullify(p_object)
    end do
    nullify(p_name)

end subroutine sf_append_json_file

!-----------------------------------------------------------------------------------------------------------
subroutine sf_end_json_file(t)
    type(plane_t) :: t
    character(100) :: filename
    integer :: iunit

    !Write File
    filename = trim(adjustl(t%master_filename))//'_derivatives.json'
!    write(*,'(A)') 'writing file '//trim(filename)//'...'
    open(newunit=iunit, file=filename, status='REPLACE')
    call json_print(t%p_json,iunit)
    close(iunit)
    call json_destroy(t%p_json)
!    write(*,'(A)') 'done.'
    write(*,'(A)') 'Derivative results written to: '//trim(filename)

end subroutine sf_end_json_file

!-----------------------------------------------------------------------------------------------------------
subroutine sf_aerocenter(t,json_command)
    type(plane_t) :: t
    type(json_value),intent(in),pointer :: json_command
    character(len=:),allocatable :: cval
    real :: delta,alpha,CG(3)
    real :: CA(3), CN(3), Cm0(3)
    real :: CAa,CAaa,CNa,CNaa,Cm0a,Cm0aa,Cmac
    real :: xac, zac,denom
    type(json_value),pointer    :: p_root
    character(100) :: filename
    integer :: iunit

    delta = 0.5*pi/180.0 !radians
    alpha = t%alpha
    CG(:) = t%CG(:)
    t%CG(:) = 0.0

    write(*,*)
    write(*,*) '---------- Aerodynamic Center  -----------'

    ! Deriv w/respect to alpha
    t%alpha = alpha - delta
    call plane_run_current(t)
    CA(1)  = -t%GF(t%nrealwings+1,3,1)
    CN(1)  = -t%GF(t%nrealwings+1,3,3)
    Cm0(1) =  t%GM(t%nrealwings+1,3,2)

    t%alpha = alpha
    call plane_run_current(t)
    CA(2)  = -t%GF(t%nrealwings+1,3,1)
    CN(2)  = -t%GF(t%nrealwings+1,3,3)
    Cm0(2) =  t%GM(t%nrealwings+1,3,2)

    t%alpha = alpha + delta
    call plane_run_current(t)
    CA(3)  = -t%GF(t%nrealwings+1,3,1)
    CN(3)  = -t%GF(t%nrealwings+1,3,3)
    Cm0(3) =  t%GM(t%nrealwings+1,3,2)

    !Second order derivatives
    CAa = (CA(3) - CA(1))/2.0/delta
    CNa = (CN(3) - CN(1))/2.0/delta
    Cm0a = (Cm0(3) - Cm0(1))/2.0/delta

    CAaa = (CA(1) - 2.0*CA(2) + CA(3))/delta/delta
    CNaa = (CN(1) - 2.0*CN(2) + CN(3))/delta/delta
    Cm0aa = (Cm0(1) - 2.0*Cm0(2) + Cm0(3))/delta/delta

    denom = CNa*CAaa - CAa*CNaa

    xac = (CAa*Cm0aa - Cm0a*CAaa)/denom
    zac = (CNa*Cm0aa - Cm0a*CNaa)/denom
    Cmac = Cm0(2) + xac*CN(2) - zac*CA(2)

    xac = -xac*t%long_r    !negative because in stability coordinates, not aero coordinates
    zac = -zac*t%long_r  !negative because in stability coordinates, not aero coordinates

    write(*,*) '         x_ac = ',xac
    write(*,*) '         z_ac = ',zac
    write(*,*)
    write(*,*) '        Cm_ac = ',Cmac
    write(*,*) '------------------------------------------'
    write(*,*)

    t%alpha = alpha
    t%CG(:) = CG(:)

    !Get filename if specified
    filename = trim(adjustl(t%master_filename))//'_aerocenter.json'
    call myjson_get(json_command, 'filename', cval, filename)
    filename = trim(adjustl(cval))

    call json_value_create(p_root)           ! create the value and associate the pointer
    call to_object(p_root,trim(filename))    ! add the file name as the name of the overall structure

    call json_value_add( p_root, 'x_ac',   xac)
    call json_value_add( p_root, 'z_ac',   zac)
    call json_value_add( p_root, 'Cm_ac', Cmac)

    !Write File
    open(newunit=iunit, file=filename, status='REPLACE')
    call json_print(p_root,iunit)
    close(iunit)
    call json_destroy(p_root)
    write(*,'(A)') 'Aerodynamic Center results written to: '//trim(filename)

end subroutine sf_aerocenter

!-----------------------------------------------------------------------------------------------------------
subroutine sf_stallonset(t,json_command)
    type(plane_t) :: t
    type(json_value),intent(in),pointer :: json_command
    character(len=:),allocatable :: cval
    type(section_t),pointer :: si
    character(1000) :: note
    real :: alpha,maxdiff,start_alpha
    integer :: ialpha,iwing,isec,stalled,swing,ssec
    type(json_value),pointer    :: p_root
    character(100) :: filename
    integer :: iunit

    alpha = t%alpha !store to reset after finished
    note = "Linear Stall Onset."

    call myjson_get(json_command, 'start_alpha', start_alpha, 0.0)
    start_alpha = (start_alpha+1.0)*pi/180.0  ! Add 1 degree, but ialpha loop starts at i=-10

    write(*,*)
    write(*,*) '---------- Evaluating Stall Onset -----------'
    write(*,*)
    write(*,*) 'alpha [deg]        Aircraft CL'

    stalled = 0
    maxdiff = 0.0
    swing = 0
    ssec = 0
    do ialpha = -10, 300, 1
        t%alpha = REAL(ialpha)/10.0*pi/180.0 + start_alpha
        call plane_run_current(t)
        write(*,*) t%alpha*180.0/pi,t%GL(t%nrealwings+1,3)
        do iwing=1,t%nwings
            do isec = 1,t%wings(iwing)%nSec
                si => t%wings(iwing)%sec(isec)
                if(sec_stalled(si).gt.0.0) then
                    write(*,*) 'Wing ',trim(t%wings(iwing)%name),' stalled at spanwise location = ',si%percent_c*100.0,'%'
                    stalled = 1
                    if(sec_stalled(si) > maxdiff) then
                        swing = iwing
                        ssec = isec
                        maxdiff = sec_stalled(si)
                    end if
                end if
            end do
        end do
        if(stalled.eq.1) exit
    end do

    write(*,*)
    si => t%wings(swing)%sec(ssec)
    write(*,*) 'Stall onset on wing = ',trim(t%wings(swing)%name),' at spanwise location = ',si%percent_c*100.0,'%'
    write(*,*)

    !Get filename if specified
    filename = trim(adjustl(t%master_filename))//'_stallonset.json'
    call myjson_get(json_command, 'filename', cval, filename)
    filename = trim(adjustl(cval))

    call json_value_create(p_root)           ! create the value and associate the pointer
    call to_object(p_root,trim(filename))    ! add the file name as the name of the overall structure

    call json_value_add( p_root, 'alpha',   t%alpha*180.0/pi)
    call json_value_add( p_root, 'CL',   t%GL(t%nrealwings+1,3))
    call json_value_add( p_root, 'wing',   trim(t%wings(swing)%name))
    call json_value_add( p_root, 'spanLocation',   si%percent_c)

    !Write File
    open(newunit=iunit, file=filename, status='REPLACE')
    call json_print(p_root,iunit)
    close(iunit)
    call json_destroy(p_root)
    write(*,'(A)') 'Stall onset results written to: '//trim(filename)

    write(*,*) '--------------------------------------------'
    write(*,*)

    t%alpha = alpha
end subroutine sf_stallonset

!-----------------------------------------------------------------------------------------------------------
subroutine sf_pitch_trim(t,json_command)
    type(plane_t) :: t
    type(json_value),intent(in),pointer :: json_command
    character(len=:),allocatable :: cval
    real :: CL_target,Cm_target,CW_target,climb,thrust_x,thrust_z,thrust_a,alpha_temp,de_temp
    real :: delta,relaxation,alpha,de,maxres,residual(2),ans1(2),ans2(2),gradient(2,2),dG(2)
    integer :: trimType
    integer :: icontrol,found,iter,maxiter
    type(json_value),pointer    :: p_root
    character(100) :: controlname, filename
    integer :: iunit

    !Get filename if specified
    filename = trim(adjustl(t%master_filename))//'_pitchtrim.json'
    call myjson_get(json_command, 'filename', cval, filename)
    filename = trim(adjustl(cval))

    write(*,*)
    write(*,*) '---------- Trimming Aircraft in Pitch -----------'

    !Get control surface name
    controlname = 'elevator'
    call myjson_get(json_command, 'control', cval, controlname)
    controlname = trim(adjustl(cval))

    found = 0
    do icontrol=1,t%ncontrols
        if(t%controls(icontrol)%name == trim(controlname)) then
            found = icontrol
            exit
        end if
    end do

    if(found.eq.0) then
        write(*,*)
        write(*,*) '---------------------------------------------------------------'
        write(*,*) 'Error: No control surface found. Exiting trim command.'
        write(*,*) '---------------------------------------------------------------'
        write(*,*)
        open(newunit=iunit, file=filename, status='REPLACE')
        write(iunit,*) '{ "error" : "No elevator found."}'
        close(iunit)
        return
    end if

    write(*,*) 'Using control surface: ',trim(controlname)

    call myjson_get(json_command, 'CL', CL_target, -1.0)  !!??? Is it ever ok for CL < 0? Use NOTDEFD?
!    call json_get(json_command,'CL',CL_target,json_found)

    if(CL_target < 0.0) then !Trim using Cw, and thrust !Not done yet. Fix this!
        call json_clear_exceptions()
        trimType = 2

        call myjson_get(json_command, 'CW', CW_target)
        call myjson_get(json_command, 'climb', climb, 0.0); climb = climb*pi/180.0
        call myjson_get(json_command, 'thrust.x', thrust_x, 0.0);
        call myjson_get(json_command, 'thrust.z', thrust_z, 0.0);
        call myjson_get(json_command, 'thrust.angle', thrust_a, 0.0);

    else !trim using CL, Cm
        trimType = 1
        call myjson_get(json_command, 'Cm', Cm_target, 0.0);
    end if

    !store alpha and de in case no solution is found
    alpha_temp = t%alpha !radians
    de_temp = t%controls(icontrol)%deflection !radians

    icontrol = found
    alpha = t%alpha !radians
    de = t%controls(icontrol)%deflection !radians

    call myjson_get(json_command, 'delta', delta, 0.5);
    call myjson_get(json_command, 'convergence', maxres, 1.0e-10);
    call myjson_get(json_command, 'relaxation', relaxation, 1.0);
    call myjson_get(json_command, 'maxiter', maxiter, 50);

    write(*,*)
    write(*,*) '    derivative step size [deg] = ',delta
    write(*,*) '                   convergence = ',maxres
    delta = 0.5*pi/180.0 !radians

    write(*,*)
    write(*,*) '   alpha[deg]                ',trim(controlname),'[deg]            CL                        &
                  &Cm                        CL_Residual              Cm_Residual'
    call sf_pitch_trim_residual(t,icontrol,alpha,de,CL_target,Cm_target,residual(:))
    write(*,*) alpha*180.0/pi,de*180.0/pi,t%GL(t%nrealwings+1,3), t%GM(t%nrealwings+1,3,2), residual(:)

    iter = 0
    do while ((abs(residual(1)) + abs(residual(2))) > maxres)

        alpha = alpha + delta
        call sf_pitch_trim_residual(t,icontrol,alpha,de,CL_target,Cm_target,ans1(:))
        alpha = alpha - 2.0*delta
        call sf_pitch_trim_residual(t,icontrol,alpha,de,CL_target,Cm_target,ans2(:))
        alpha = alpha + delta

        gradient(:,1) = (ans1(:) - ans2(:))/2.0/delta

        de = de + delta
        call sf_pitch_trim_residual(t,icontrol,alpha,de,CL_target,Cm_target,ans1(:))
        de = de - 2.0*delta
        call sf_pitch_trim_residual(t,icontrol,alpha,de,CL_target,Cm_target,ans2(:))
        de = de + delta

        gradient(:,2) = (ans1(:) - ans2(:))/2.0/delta

!        write(*,*) 'gradient'
!        write(*,*) gradient(1,:)
!        write(*,*) gradient(2,:)

        call math_snyder_ludcmp(gradient,2)
        call math_snyder_lusolv(gradient,-residual,dG,2)

!        write(*,*) 'dG'
!        write(*,*) dG(1)*180.0/pi
!        write(*,*) dG(2)*180.0/pi
!        write(*,*)

        if((dG(1).ne.dG(1)) .or. (dG(2).ne.dG(2))) then
            write(*,*) 'Error: Control surface not set to affect pitch trim.'
            t%alpha  = alpha_temp !radians
            t%controls(icontrol)%deflection = de_temp !radians
            open(newunit=iunit, file=filename, status='REPLACE')
            write(iunit,*) '{ "error" : "No elevator mixing found."}'
            close(iunit)
            return
        end if

        alpha = alpha + relaxation*dG(1)
        de    = de    + relaxation*dG(2)
        call sf_pitch_trim_residual(t,icontrol,alpha,de,CL_target,Cm_target,residual(:))
        write(*,*) alpha*180.0/pi,de*180.0/pi,t%GL(t%nrealwings+1,3), t%GM(t%nrealwings+1,3,2),residual(:)
        iter = iter + 1
        if(iter > maxiter) then
            write(*,*) 'Error. Maximum iterations reached. Trim point not found. Exiting trim.'
            t%alpha  = alpha_temp !radians
            t%controls(icontrol)%deflection = de_temp !radians
            open(newunit=iunit, file=filename, status='REPLACE')
            write(iunit,*) '{ "error" : "We are having difficulty trimming. &
                           &Try making your elevator larger or adjusting your center of gravity."}'
            close(iunit)
            return
        end if
    end do

    write(*,*)
    write(*,*) '------------- Results --------------------'
    write(*,*) '      alpha [deg] = ',alpha*180.0/pi
    write(*,*) '   ',trim(controlname),' [deg] = ',de*180.0/pi
    write(*,*) '------------------------------------------'
    write(*,*)

    if(abs(de*180.0/pi) > 90.0) then
        write(*,*) 'Error. Maximum elevator deflection is extremely large. Exiting trim.'
        t%alpha  = alpha_temp !radians
        t%controls(icontrol)%deflection = de_temp !radians
        open(newunit=iunit, file=filename, status='REPLACE')
        write(iunit,*) '{ "error" : "We are having difficulty trimming. &
                       &Try making your elevator larger or adjusting your center of gravity."}'
        close(iunit)
        return
    end if

    t%alpha = alpha
    t%controls(icontrol)%deflection = de
    call json_value_create(p_root)           ! create the value and associate the pointer
    call to_object(p_root,trim(filename))    ! add the file name as the name of the overall structure

    call json_value_add( p_root, 'alpha',   alpha*180.0/pi)
    call json_value_add( p_root, 'elevator',   de*180.0/pi)

    !Write File
    open(newunit=iunit, file=filename, status='REPLACE')
    call json_print(p_root,iunit)
    close(iunit)
    call json_destroy(p_root)
    write(*,'(A)') 'Trim results written to: '//trim(filename)

end subroutine sf_pitch_trim

!-----------------------------------------------------------------------------------------------------------
subroutine sf_pitch_trim_residual(t,icontrol,alpha,de,CL_target,Cm_target,ans)
    type(plane_t) :: t
    integer :: icontrol
    real :: alpha,de,CL_target,Cm_target,ans(2)

    t%alpha = alpha
    t%controls(icontrol)%deflection = de;
    call plane_run_current(t)
    ans(1) = t%GL(t%nrealwings+1,3) - CL_target
    ans(2) = t%GM(t%nrealwings+1,3,2) - Cm_target
!write(*,*) alpha*180.0/pi,de*180.0/pi,t%GL(t%nrealwings+1,3), t%GM(t%nrealwings+1,3,2),ans(:)

end subroutine sf_pitch_trim_residual

!-----------------------------------------------------------------------------------------------------------
subroutine sf_target_CL(t,json_command)
    type(plane_t) :: t
    type(json_value),intent(in),pointer :: json_command
    character(len=:),allocatable :: cval
    real :: CL_target,alpha_temp
    real :: delta,relaxation,maxres,residual,ans1,ans2,gradient
    integer :: iter,maxiter
    type(json_value),pointer    :: p_root
    character(100) :: filename
    integer :: iunit

    !Get filename if specified
    filename = trim(adjustl(t%master_filename))//'_targetCL.json'
    call myjson_get(json_command, 'filename', cval, filename)
    filename = trim(adjustl(cval))

    write(*,*)
    write(*,*) '---------- Finding alpha to target CL -----------'

    call myjson_get(json_command, 'CL', CL_target);

    !store alpha and de in case no solution is found
    alpha_temp = t%alpha !radians

    call myjson_get(json_command, 'delta', delta, 0.5);
    call myjson_get(json_command, 'convergence', maxres, 1.0e-10);
    call myjson_get(json_command, 'relaxation', relaxation, 1.0);
    call myjson_get(json_command, 'maxiter', maxiter, 50);

    write(*,*)
#ifndef dnad
    write(*,*) '    derivative step size [deg] = ',delta
    write(*,*) '                   convergence = ',maxres
    delta = 0.5*pi/180.0 !radians

    write(*,*)
    write(*,*) '   alpha[deg]                CL                        CL_Residual'
    call sf_target_CL_residual(t,CL_target,residual)
    write(*,*) t%alpha*180.0/pi,t%GL(t%nrealwings+1,3), residual
#else
    write(*,*) '                   convergence = ',maxres

    write(*,*)
    write(*,*) '   alpha[deg]                CL                        d(CL)/d(alpha)'
    call plane_run_current(t)
    write(*,*) t%alpha%x*180.0/pi, t%GL(t%nrealwings+1,3)%x, t%GL(t%nrealwings+1,3)%dx(1)
    residual = t%GL(t%nrealwings+1,3) - CL_target
#endif

    iter = 0
    do while (abs(residual) > maxres)
#ifndef dnad
        t%alpha = t%alpha + delta
        call sf_target_CL_residual(t,CL_target,ans1)
        t%alpha = t%alpha - 2.0*delta
        call sf_target_CL_residual(t,CL_target,ans2)
        t%alpha = t%alpha + delta

        gradient = (ans1 - ans2)/2.0/delta
        t%alpha = t%alpha - relaxation*residual/gradient

        call sf_target_CL_residual(t,CL_target,residual)
        write(*,*) t%alpha*180.0/pi, t%GL(t%nrealwings+1,3), residual
#else
        t%alpha = t%alpha + (CL_target - t%GL(t%nrealwings+1, 3)) / t%GL(t%nrealwings+1, 3)%dx(1)
        t%alpha%dx(1) = t%GL(t%nrealwings+1, 3)%dx(1) / t%GL(t%nrealwings+1, 3)%dx(1)

        call plane_run_current(t)
        residual = t%GL(t%nrealwings+1, 3) - CL_target
        write(*,*) t%alpha%x*180.0/pi, t%GL(t%nrealwings+1,3)%x, t%GL(t%nrealwings+1,3)%dx(1)
#endif

        iter = iter + 1
        if(iter > maxiter) then
            write(*,*) 'Error. Maximum iterations reached. Trim point not found. Exiting trim.'
            t%alpha  = alpha_temp !radians
            open(newunit=iunit, file=filename, status='REPLACE')
            write(iunit,*) '{ "error" : "We are having difficulty hitting that target CL. &
                           &Try checking the lift slope of your main wing."}'
            close(iunit)
            return
        end if
    end do

    write(*,*)
    write(*,*) '------------- Results --------------------'
    write(*,*) '      alpha [deg] = ',t%alpha*180.0/pi
    write(*,*) '               CL = ',t%GL(t%nrealwings+1,3)
    write(*,*) '               CD = ',t%GD(t%nrealwings+1,3)
    write(*,*) '------------------------------------------'
    write(*,*)

    call json_value_create(p_root)           ! create the value and associate the pointer
    call to_object(p_root,trim(filename))    ! add the file name as the name of the overall structure

    call json_value_add( p_root, 'alpha',   t%alpha*180.0/pi)
    call json_value_add( p_root, 'CL',   t%GL(t%nrealwings+1,3))
    call json_value_add( p_root, 'CD',   t%GD(t%nrealwings+1,3))

    !Write File
    open(newunit=iunit, file=filename, status='REPLACE')
    call json_print(p_root,iunit)
    close(iunit)
    call json_destroy(p_root)
    write(*,'(A)') 'Trim results written to: '//trim(filename)

end subroutine sf_target_CL

!-----------------------------------------------------------------------------------------------------------
subroutine sf_target_CL_residual(t,CL_target,ans)
    type(plane_t) :: t
    real :: CL_target,ans

    call plane_run_current(t)
    ans = t%GL(t%nrealwings+1,3) - CL_target
end subroutine sf_target_CL_residual

!-----------------------------------------------------------------------------------------------------------
subroutine sf_target(t)
    type(plane_t) :: t
    character(len=:), allocatable :: target_var, change_var
    real :: result,x0,x1,xnew,f0,f1, target_val, tolerance

    !Read json target info
    call myjson_get(t%json, 'run.target.variable', target_var)
    call myjson_get(t%json, 'run.target.value', target_val)
    call myjson_get(t%json, 'run.target.tolerance', tolerance )
    call myjson_get(t%json, 'run.target.change', change_var)

    write(*,*) '    target variable = ',trim(target_var)
    write(*,*) '       target value = ',target_val
    write(*,*) '          tolerance = ',tolerance
    write(*,*) '    change variable = ',trim(change_var)

    run_type = 'forces'

    x0 = 0.0
    call sf_set_variable(t,change_var,x0)
    call sf_run_single(t,target_var,result)
    f0 = result-target_val

    x1 = x0 + 0.01
    call sf_set_variable(t,change_var,x1)
    call sf_run_single(t,target_var,result)
    f1 = result-target_val

write(*,*)
write(*,*) '  '//trim(change_var)//'     |    '//trim(target_var)
write(*,*) '  --------------------------------------------'
    do while (abs(f1)>tolerance)
        xnew = x1-f1*(x1-x0)/(f1-f0)
        call sf_set_variable(t,change_var,xnew)
        call sf_run_single(t,target_var,result)
        x0 = x1
        f0 = f1
        x1 = xnew
        f1 = result-target_val
write(*,*) x1,result
    end do
write(*,*)

    !Run the final configuration
    t%verbose = 1
    call plane_run_current(t)
end subroutine sf_target

!-----------------------------------------------------------------------------------------------------------
subroutine sf_report(t,json_command)
    type(plane_t) :: t
    type(json_value),intent(in),pointer :: json_command
    type(json_file) :: f_json    !the JSON structure read from the file
    type(json_value), pointer :: c_var,p_root
    character(len=:),allocatable :: cval
    character(100) :: filename,var_name, var_file
    real :: var_value
    integer :: ios,i,nvars,iunit,save_file

    !Get filename if specified
    filename = trim(adjustl(t%master_filename))//'_report.json'
    call myjson_get(json_command, 'filename', cval, filename)
    filename = trim(adjustl(cval))

    call json_value_create(p_root)           ! create the value and associate the pointer
    call to_object(p_root,trim(filename))    ! add the file name as the name of the overall structure

    !Include variables
    nvars = json_value_count(json_command)
    do i=1,nvars
        call json_value_get(json_command,i,c_var)
        if(trim(c_var%name).eq.'filename') cycle
        if(trim(c_var%name).eq.'run') cycle

        call myjson_get(json_command,trim(c_var%name)//'.name', cval); var_name = trim(cval)
        call myjson_get(json_command,trim(c_var%name)//'.file', cval); var_file = trim(cval)
        call myjson_get(json_command,trim(c_var%name)//'.save', save_file, 0);

        call f_json%load_file(filename = var_file);       call json_check()
        call myjson_get(f_json, trim(var_name), var_value)
!        call f_json%get(trim(var_name),   var_value);     call json_check()
        call json_value_add(p_root, trim(c_var%name), var_value)
        write(*,*) '   ',trim(c_var%name),' set from ',trim(var_name),' in file ',trim(var_file),' = ',var_value
        if(save_file .eq. 1) then
            open(unit = 10, File = trim(c_var%name)//'.txt', action = 'write', iostat = ios)
            write(10,'(2ES25.16)') var_value
            close(10)
            write(*,*) '   ',trim(c_var%name),' saved to ',trim(c_var%name)//'.txt'
        end if

    end do
    write(*,*)

    !Write File
    write(*,*) 'Saving Report to ',trim(filename)
    open(newunit=iunit, file=filename, status='REPLACE')
    call json_print(p_root,iunit)
    close(iunit)
    call json_destroy(p_root)
end subroutine sf_report

!-----------------------------------------------------------------------------------------------------------
subroutine sf_fitness_old(t)
    type(plane_t) :: t
    type(json_file) :: f_json    !the JSON structure read from the file:
    character(len=:),allocatable :: cval
    character(100) :: fitness_var, fn
    real :: fitness_val
    integer :: ios

    !Evaluate Fitness Function if present
    call json_clear_exceptions()
    call t%json%get('run.fitness.variable',    cval);   call json_check();    fitness_var = trim(cval)
    call t%json%get('run.fitness.file',        cval);   call json_check();    fn = trim(cval)

    call f_json%load_file(filename = fn);                  call json_check()
    call myjson_get(f_json, trim(fitness_var), fitness_val)

    open(unit = 10, File = 'fitness.txt', action = 'write', iostat = ios)
    write(10,'(2ES25.16)') fitness_val
    close(10)
    write(*,*) 'Fitness = ',fitness_val
end subroutine sf_fitness_old

!-----------------------------------------------------------------------------------------------------------
subroutine sf_optimize(t)
    type(plane_t) :: t
    type(json_file) :: f_json    !the JSON structure read from the file:
    character(len=:),allocatable :: cval
    character(100) :: fn
    real :: fitness_val, Cl, CD, Cn
    integer :: ios, opt_type

    !Evaluate Fitness Function if present
    call json_clear_exceptions()
    call t%json%get('run.optimize.type',    opt_type);   call json_check();
    call t%json%get('run.fitness.file',        cval);   call json_check();    fn = trim(cval)

    call f_json%load_file(filename = fn);                  call json_check()
    call myjson_get(f_json, 'total.MyAirplane.CD', CD)
    call myjson_get(f_json, 'total.MyAirplane.Cl', Cl)
    call myjson_get(f_json, 'total.MyAirplane.Cn', Cn)

    if(opt_type.eq.1) then
        fitness_val = 1.0*sqrt((Cl-0.052950639579748)**2)+CD
    end if


    open(unit = 10, File = 'fitness.txt', action = 'write', iostat = ios)
    write(10,'(2ES25.16)') fitness_val
    close(10)
    write(*,*) 'Fitness = ',fitness_val
end subroutine sf_optimize

!-----------------------------------------------------------------------------------------------------------
subroutine sf_set_variable(t,set_var,value)
    type(plane_t) :: t
    character(100) :: set_var
    real :: value

    if(set_var.eq.'condition.alpha') t%alpha = value*pi/180.0
    if(set_var.eq.'condition.beta') t%beta = value*pi/180.0

    if(set_var.eq.'cg_x') t%CG(1) = value
    if(set_var.eq.'cg_y') t%CG(2) = value
    if(set_var.eq.'cg_z') t%CG(3) = value
end subroutine sf_set_variable

!-----------------------------------------------------------------------------------------------------------
subroutine sf_run_single(t,read_var,result) !this can only handle run force calls. no derivatives yet.
    type(plane_t) :: t
    type(json_file) :: f_json    !the JSON structure read from the file:
    character(100) :: read_var,fn
    character(120) :: command
    real :: result

    fn = trim(adjustl(t%master_filename))//'_sub.json'
    call plane_write_json_file(t,fn)
    command = 'MachUp.out '//fn
    command = trim(adjustl(command))//' > trash.txt'
    call system(command)
    fn = trim(adjustl(t%master_filename))//'_sub_forces.json'
    call f_json%load_file(filename = fn);          call json_check()
    call myjson_get(f_json, trim(read_var), result)

end subroutine sf_run_single

end module special_functions_m
!---notes---
!you have to connect to a wing with a lower id than yourself. So for loads calcs, you can sum
!by starting with the last wing and working backwards, adding the forces and moments
!to the parent wing
