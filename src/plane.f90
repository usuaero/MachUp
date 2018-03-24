module plane_m
#ifdef dnad
    use dnadmod
#define real type(dual)
#endif
    use myjson_m
    use wing_m
    implicit none

    type section_pointer_t
        type(section_t),pointer :: myp
    end type section_pointer_t

    type control_surface_t
        character(20) :: name
        integer :: is_symmetric
        real :: deflection
    end type control_surface_t

    type plane_t
        character(20) :: name
        character(100) :: master_filename

        type(json_file) :: json    !the JSON structure read from the file:
        type(json_value),pointer    :: p_json

        integer :: nSize !number of horseshoe vortices
        integer :: nwings !total number of lifting surfaces
        integer :: nrealwings !number of wings not including reflections in ground plane
        integer :: groundplane
        integer :: ncontrols

        real :: Sr
        real :: long_r,lat_r
        real :: CG(3)
        real :: alpha
        real :: beta
        real :: omega(3)
        real :: Uinf(3)
        real :: hag !height above ground

        integer :: verbose !=1 for write info, =0 for don't write info

        real,allocatable,dimension(:,:,:) :: vij
        real,allocatable,dimension(:,:) :: Amat,Ainv
        real,allocatable,dimension(:) :: Bvec
        real,allocatable,dimension(:) :: Gammas

        type(section_pointer_t),allocatable,dimension(:) :: sec
        type(wing_t),allocatable,dimension(:) :: wings
        type(control_surface_t),allocatable,dimension(:) :: controls
        real,allocatable,dimension(:,:,:) :: GF,GM !Global Forces and Moments
        real,allocatable,dimension(:,:) :: GL,GD,GLoD !Global Lift, Drag, and L/D
        type(dataset_t),allocatable,dimension(:) :: external_forces
        real,allocatable,dimension(:) :: external_force_mult

    end type plane_t

        integer :: allocated

        character(100) :: run_type
        character(20) :: solver
        !Solver Settings
        real :: jacobian_converged
        real :: jacobian_omega
        integer :: nonlinear_maxiter
        integer :: solution_exists

        !Airfoil Database
        character(200) :: DB_Airfoil
        type(airfoil_t),pointer,dimension(:) :: airfoils

contains

!-----------------------------------------------------------------------------------------------------------
subroutine plane_allocate(t)
    type(plane_t) :: t
    if(allocated.eq.1) call plane_deallocate(t)
!    write(*,*) 'Size of computational matrix: ',t%nSize

    allocate(t%sec(t%nSize))
    allocate(t%vij(t%nSize,t%nSize,3))
    allocate(t%Amat(t%nSize,t%nSize))
    allocate(t%Ainv(t%nSize,t%nSize))
    allocate(t%Bvec(t%nSize))
    allocate(t%Gammas(t%nSize))

    allocate(t%GF(t%nrealwings+1,3,3))
    allocate(t%GM(t%nrealwings+1,3,3))
    allocate(t%GL(t%nrealwings+1,3))
    allocate(t%GD(t%nrealwings+1,3))
    allocate(t%GLoD(t%nrealwings+1,3))

    allocated = 1
    t%vij = 0.0
    t%Gammas = 0.0

end subroutine plane_allocate

!-----------------------------------------------------------------------------------------------------------
subroutine plane_deallocate(t)
    type(plane_t) :: t
    deallocate(t%sec)
    deallocate(t%vij)
    deallocate(t%Amat)
    deallocate(t%Ainv)
    deallocate(t%Bvec)
    deallocate(t%Gammas)

    deallocate(t%GF)
    deallocate(t%GM)
    deallocate(t%GL)
    deallocate(t%GD)
    deallocate(t%GLoD)

    allocated = 0
end subroutine plane_deallocate

!-----------------------------------------------------------------------------------------------------------
subroutine plane_set_defaults(t)
    type(plane_t) :: t

    allocated = 0
    t%groundplane = 0
    t%hag = 0.0

    solver = 'linear'
    jacobian_converged = 1.0e-11
    jacobian_omega = 1.0
    nonlinear_maxiter = 100

end subroutine plane_set_defaults

!-----------------------------------------------------------------------------------------------------------
subroutine plane_load_json(t)
    type(plane_t) :: t
    type(json_value),pointer :: j_this, j_wing, j_cont, j_afprop, j_afread, json_command
    character(len=:),allocatable :: cval
    character(5) :: side
    integer :: loc,i,iwing,nreadwings,iforce,icontrol,iairfoil,nairfoils,iaf
    real :: sweep,dihedral,mount,washout

    call json_initialize()

!    write(*,*) 'reading input file: ',t%master_filename
    call t%json%load_file(filename = t%master_filename); call json_check()
    loc = index(t%master_filename,'.json')
    t%master_filename = t%master_filename(1:loc-1) !deletes the .json file extension

!    call t%json%print_file()

    write(*,*) 'Reading Input File...'
    !Read Data
    call t%json%get('plane.name',                        cval);     call json_check(); t%name = trim(cval)
    call myjson_get(t%json,'plane.CGx', t%CG(1))
    call myjson_get(t%json,'plane.CGy', t%CG(2))
    call myjson_get(t%json,'plane.CGz', t%CG(3))

    call myjson_get(t%json,'reference.area', t%Sr)
    call myjson_get(t%json,'reference.longitudinal_length', t%long_r)
    call myjson_get(t%json,'reference.lateral_length', t%lat_r)

    call myjson_get(t%json,'condition.alpha', t%alpha); t%alpha = t%alpha*pi/180.0
#ifdef dnad
    ! Always make alpha a design variable so it can be used in the target_CL function
    t%alpha%dx(1) = 1.0
#endif
    call myjson_get(t%json,'condition.beta', t%beta, 0.0); t%beta = t%beta*pi/180.0

    call myjson_get(t%json,'condition.omega.roll', t%omega(1), 0.0)
    call myjson_get(t%json,'condition.omega.pitch', t%omega(2), 0.0)
    call myjson_get(t%json,'condition.omega.yaw', t%omega(3), 0.0)

    call myjson_get(t%json,'condition.ground', t%hag, 0.0)
    if(t%hag.gt.0.0) t%groundplane = 1

    call t%json%get('solver.type',                       cval);     call json_check(); solver = trim(cval)
    call myjson_get(t%json,'solver.convergence', jacobian_converged, 1.0e-6)
    call myjson_get(t%json,'solver.relaxation', jacobian_omega, 0.9)
    call myjson_get(t%json,'solver.maxiter', nonlinear_maxiter, 100)

    call t%json%get('airfoil_DB',                        cval);     call json_check(); DB_Airfoil = trim(cval)
    write(*,*) 'Airfoil database located at: ',trim(DB_Airfoil)

    ! Read Controls
    write(*,*)
    call t%json%get('controls', j_this)
    if(json_failed()) then
        write(*,*) 'No controls specified.'
        call json_clear_exceptions()
        t%ncontrols = 0
    else
        t%ncontrols = json_value_count(j_this)
        write(*,*) 'Number of controls : ',t%ncontrols
        allocate(t%controls(t%ncontrols))

        do icontrol=1,t%ncontrols
            call json_value_get(j_this,icontrol,j_cont)
            t%controls(icontrol)%name = trim(j_cont%name)
            call t%json%get('controls.'//trim(j_cont%name)//'.is_symmetric',  t%controls(icontrol)%is_symmetric);
            call json_check();
            call myjson_get(t%json, 'controls.'//trim(j_cont%name)//'.deflection', t%controls(icontrol)%deflection)
!            call t%json%get('controls.'//trim(j_cont%name)//'.deflection',    t%controls(icontrol)%deflection);
!            call json_check()
            t%controls(icontrol)%deflection = pi/180.0*t%controls(icontrol)%deflection
            write(*,*) '  : ',trim(t%controls(icontrol)%name), ' ',t%controls(icontrol)%deflection*180.0/pi, ' (deg)'
        end do
        write(*,*)
    end if

    ! Read Wings
    call t%json%get('wings', j_this)
    nreadwings = json_value_count(j_this)
    t%nrealwings = nreadwings
    nairfoils = 0
    do i=1,nreadwings
        call json_value_get(j_this,i,j_wing)
        call t%json%get('wings.'//trim(j_wing%name)//'.side', cval); call json_check(); side = trim(cval)
        if(side == 'both') t%nrealwings = t%nrealwings + 1

        call t%json%get('wings.'//trim(j_wing%name)//'.airfoils', j_afread)
        nairfoils = nairfoils + json_value_count(j_afread)
    end do

    if(t%groundplane == 1) then
        t%nwings = 2*t%nrealwings
    else
        t%nwings = t%nrealwings
    end if

    allocate(t%wings(t%nwings))
    allocate(airfoils(nairfoils))
    iforce = 0
    allocate(t%external_forces(iforce))
    allocate(t%external_force_mult(iforce))

    iwing = 0
    iairfoil = 0
    do i=1,nreadwings
        iwing = iwing + 1
        call json_value_get(j_this,i,j_wing)
        t%wings(iwing)%name = trim(j_wing%name)



        call t%json%get('wings.'//trim(j_wing%name)//'.ID',                 t%wings(iwing)%ID); call json_check()
        call t%json%get('wings.'//trim(j_wing%name)//'.side',                            cval); call json_check();
                                                                                t%wings(iwing)%orig_side = trim(cval)
        call t%json%get('wings.'//trim(j_wing%name)//'.connect.ID',  t%wings(iwing)%connectid); call json_check()
        call t%json%get('wings.'//trim(j_wing%name)//'.connect.location',                cval); call json_check();
                                                                                t%wings(iwing)%connectend = trim(cval)
        call myjson_get(t%json, 'wings.'//trim(j_wing%name)//'.connect.dx', t%wings(iwing)%doffset(1))
!        call t%json%get('wings.'//trim(j_wing%name)//'.connect.dx', t%wings(iwing)%doffset(1)); call json_check()
        call myjson_get(t%json, 'wings.'//trim(j_wing%name)//'.connect.dy', t%wings(iwing)%doffset(2))
!        call t%json%get('wings.'//trim(j_wing%name)//'.connect.dy', t%wings(iwing)%doffset(2)); call json_check()
        call myjson_get(t%json, 'wings.'//trim(j_wing%name)//'.connect.dz', t%wings(iwing)%doffset(3))
!        call t%json%get('wings.'//trim(j_wing%name)//'.connect.dz', t%wings(iwing)%doffset(3)); call json_check()
        call myjson_get(t%json, 'wings.'//trim(j_wing%name)//'.connect.yoffset', t%wings(iwing)%dy)
!        call t%json%get('wings.'//trim(j_wing%name)//'.connect.yoffset',    t%wings(iwing)%dy); call json_check()
        call myjson_get(t%json, 'wings.'//trim(j_wing%name)//'.span', t%wings(iwing)%span)
!        call t%json%get('wings.'//trim(j_wing%name)//'.span',             t%wings(iwing)%span); call json_check()
        call myjson_get(t%json, 'wings.'//trim(j_wing%name)//'.sweep', sweep)
!        call t%json%get('wings.'//trim(j_wing%name)//'.sweep',                          sweep); call json_check()
        call myjson_get(t%json, 'wings.'//trim(j_wing%name)//'.dihedral', dihedral)
!        call t%json%get('wings.'//trim(j_wing%name)//'.dihedral',                    dihedral); call json_check()
        call myjson_get(t%json, 'wings.'//trim(j_wing%name)//'.mounting_angle', mount)
!        call t%json%get('wings.'//trim(j_wing%name)//'.mounting_angle',                 mount); call json_check()
        call myjson_get(t%json, 'wings.'//trim(j_wing%name)//'.washout', washout)
!        call t%json%get('wings.'//trim(j_wing%name)//'.washout',                      washout); call json_check()
        call myjson_get(t%json, 'wings.'//trim(j_wing%name)//'.root_chord', t%wings(iwing)%chord_1)
!        call t%json%get('wings.'//trim(j_wing%name)//'.root_chord',    t%wings(iwing)%chord_1); call json_check()
        call myjson_get(t%json, 'wings.'//trim(j_wing%name)//'.tip_chord', t%wings(iwing)%chord_2)
!        call t%json%get('wings.'//trim(j_wing%name)//'.tip_chord',     t%wings(iwing)%chord_2); call json_check()

        call t%json%get('wings.'//trim(j_wing%name)//'.sweep_definition', t%wings(iwing)%sweep_definition);
        if(json_failed()) t%wings(iwing)%sweep_definition=1
        call json_clear_exceptions()

        ! Read Airfoils
        call t%json%get('wings.'//trim(j_wing%name)//'.airfoils', j_afread)
        nairfoils = json_value_count(j_afread) !just for this wing
        t%wings(iwing)%nairfoils = nairfoils
        call wing_allocate_airfoils(t%wings(iwing))

        do iaf = 1, nairfoils
            iairfoil = iairfoil + 1
            call json_value_get(j_afread,iaf,json_command)
            airfoils(iairfoil)%name = trim(json_command%name)
            call t%json%get('wings.'//trim(j_wing%name)//'.airfoils.'//trim(json_command%name)//'.properties',  j_afprop);
            if(json_failed()) then !Read from airfoil database
                call json_clear_exceptions()
                call plane_load_airfoil(t,iairfoil,trim(json_command%name)//'.',0)
            else !Read from local file
                call plane_load_airfoil(t,iairfoil,'wings.'//trim(j_wing%name)//'.airfoils.'//trim(json_command%name)//'.',1)
            end if
            t%wings(iwing)%airfoils(iaf)%p => airfoils(iairfoil)
        end do


        call t%json%get('wings.'//trim(j_wing%name)//'.grid',             t%wings(iwing)%nSec); call json_check()

        !Grid clustering parameters
        call myjson_get(t%json, 'wings.'//trim(j_wing%name)//'.root_clustering', t%wings(iwing)%root_clustering, 1)
        call myjson_get(t%json, 'wings.'//trim(j_wing%name)//'.tip_clustering', t%wings(iwing)%tip_clustering, 1)

        !control surface defs
        call myjson_get(t%json, 'wings.'//trim(j_wing%name)//'.control.span_root', t%wings(iwing)%control_span_root, -1.0)
!        call t%json%get('wings.'//trim(j_wing%name)//'.control.span_root', t%wings(iwing)%control_span_root);

        if(t%wings(iwing)%control_span_root < 0.0) then
            call json_clear_exceptions()
            t%wings(iwing)%has_control_surface = 0
        else
            t%wings(iwing)%has_control_surface = 1
!            call t%json%get('wings.'//trim(j_wing%name)//'.control.span_root', t%wings(iwing)%control_span_root);
!            call json_check()
            call myjson_get(t%json, 'wings.'//trim(j_wing%name)//'.control.span_tip', t%wings(iwing)%control_span_tip)
!            call t%json%get('wings.'//trim(j_wing%name)//'.control.span_tip',  t%wings(iwing)%control_span_tip );
!            call json_check()
            call myjson_get(t%json, 'wings.'//trim(j_wing%name)//'.control.chord_root', t%wings(iwing)%control_chord_root)
!            call t%json%get('wings.'//trim(j_wing%name)//'.control.chord_root', t%wings(iwing)%control_chord_root);
!            call json_check()
            call myjson_get(t%json, 'wings.'//trim(j_wing%name)//'.control.chord_tip', t%wings(iwing)%control_chord_tip)
!            call t%json%get('wings.'//trim(j_wing%name)//'.control.chord_tip',  t%wings(iwing)%control_chord_tip );
!            call json_check()
            call t%json%get('wings.'//trim(j_wing%name)//'.control.is_sealed',  t%wings(iwing)%control_is_sealed);
            call json_check()
        end if

        t%wings(iwing)%sweep = sweep*pi/180.0
        t%wings(iwing)%dihedral = dihedral*pi/180.0
        t%wings(iwing)%mount = mount*pi/180.0
        t%wings(iwing)%washout = washout*pi/180.0

        !Distributions defined by files
        t%wings(iwing)%f_sweep         = 'none'
        t%wings(iwing)%f_dihedral      = 'none'
        t%wings(iwing)%f_washout       = 'none'
        t%wings(iwing)%f_chord         = 'none'
        t%wings(iwing)%f_af_ratio      = 'none'
        t%wings(iwing)%f_EIGJ          = 'none'
        t%wings(iwing)%f_elastic_twist = 'none'

        !Read Distribution Files
        call t%json%get('wings.'//trim(j_wing%name)//'.chord_file', cval);
        if(.NOT.json_failed()) t%wings(iwing)%f_chord = cval
        call json_clear_exceptions()
        call t%json%get('wings.'//trim(j_wing%name)//'.sweep_file', cval);
        if(.NOT.json_failed()) t%wings(iwing)%f_sweep = cval
        call json_clear_exceptions()
        call t%json%get('wings.'//trim(j_wing%name)//'.dihedral_file', cval);
        if(.NOT.json_failed()) t%wings(iwing)%f_dihedral=cval
        call json_clear_exceptions()
        call t%json%get('wings.'//trim(j_wing%name)//'.washout_file', cval);
        if(.NOT.json_failed()) t%wings(iwing)%f_washout = cval
        call json_clear_exceptions()
        call t%json%get('wings.'//trim(j_wing%name)//'.af_ratio_file', cval);
        if(.NOT.json_failed()) t%wings(iwing)%f_af_ratio = cval
        call json_clear_exceptions()

        !make sure wing is totally set up before continuing. Here it will mirror it to the other side
        if(t%wings(iwing)%orig_side .eq. 'both') then
            t%wings(iwing)%side = 'right'
            t%wings(iwing+1) = t%wings(iwing)
            t%wings(iwing+1)%side = 'left'
            iwing = iwing + 1
        else
            t%wings(iwing)%side = t%wings(iwing)%orig_side
        end if
    end do

    if(t%groundplane .eq. 1) then
        do iwing=1,t%nrealwings
            t%wings(iwing+t%nrealwings) = t%wings(iwing) !make a copy
        end do
    end if
    write(*,*) 'Successfully read input file.'
    write(*,*)
end subroutine plane_load_json

!-----------------------------------------------------------------------------------------------------------
subroutine plane_load_airfoil(t,i,prefix,local)
    type(plane_t) :: t
    type(json_file) :: f_json    !the JSON structure read from the file:
    character(len=*) :: prefix
    integer :: local
    integer :: i
    character(100) :: fn,datafilename
    character(len=:),allocatable :: cval

    if(local.eq.1) then
        fn = trim(t%master_filename)//'.json'
    else
        fn = trim(adjustl(DB_Airfoil))//'/'//trim(adjustl(airfoils(i)%name))//'.json'
    end if

    write(*,*) 'Reading airfoil properties from file: ',trim(fn)
    call f_json%load_file(filename = trim(fn)); call json_check()
!    call f_json%print_file()

    call f_json%get(trim(prefix)//'properties.type', cval);  call json_check(); airfoils(i)%properties_type = trim(cval)

    select case (airfoils(i)%properties_type)
        case ('linear')
            call myjson_get(f_json, trim(prefix)//'properties.alpha_L0', airfoils(i)%aL0);
            call myjson_get(f_json, trim(prefix)//'properties.CL_alpha', airfoils(i)%CLa);
            call myjson_get(f_json, trim(prefix)//'properties.Cm_L0', airfoils(i)%CmL0);
            call myjson_get(f_json, trim(prefix)//'properties.Cm_alpha', airfoils(i)%Cma);
            call myjson_get(f_json, trim(prefix)//'properties.CD0', airfoils(i)%CD0);
            call myjson_get(f_json, trim(prefix)//'properties.CD0_L', airfoils(i)%CD0L);
            call myjson_get(f_json, trim(prefix)//'properties.CD0_L2', airfoils(i)%CD0L2);
            call myjson_get(f_json, trim(prefix)//'properties.CL_max', airfoils(i)%CLmax, -1.0);
            airfoils(i)%has_data_file = 0
        case ('datafile')
            call f_json%get(trim(prefix)//'properties.filename', cval); call json_check(); datafilename = trim(cval)
            datafilename = trim(adjustl(DB_Airfoil))//'/'//trim(adjustl(datafilename))
            call af_create_from_data_file(airfoils(i),datafilename)
            airfoils(i)%has_data_file = 1
        case default
            write(*,*) 'Invalid airfoil properties type! Aborting program.'
            stop
    end select

    write(*,*) 'Loaded airfoil: ',airfoils(i)%name
end subroutine plane_load_airfoil

!-----------------------------------------------------------------------------------------------------------
subroutine plane_write_json_file(t,filename)
    type(plane_t) :: t
    type(json_value),pointer    :: p_root, p_plane, p_reference, p_condition, p_wings, p_solver, p_run, p_wingi, p_connect

    character(100) :: filename
    integer :: ios,iwing,iairfoil,iunit

    iairfoil = 0
    iwing = 0
    ios = 0

    !root
    call json_value_create(p_root)           ! create the value and associate the pointer
    call to_object(p_root,trim(filename))    ! add the file name as the name of the overall structure

    !plane structure:
    call json_value_create(p_plane)             !an object
    call to_object(p_plane,'plane')
    call json_value_add(p_root, p_plane)
    call json_value_add(        p_plane, 'name',            trim(t%name))
    call json_value_add(        p_plane, 'CGx',            t%CG(1))
    call json_value_add(        p_plane, 'CGy',            t%CG(2))
    call json_value_add(        p_plane, 'CGz',            t%CG(3))
    nullify(p_plane)

    !reference structure:
    call json_value_create(p_reference)             !an object
    call to_object(p_reference,'reference')
    call json_value_add(p_root, p_reference)
    call json_value_add(        p_reference, 'area',                    t%Sr)
    call json_value_add(        p_reference, 'longitudinal_length',     t%long_r)
    call json_value_add(        p_reference, 'lateral_length',          t%lat_r)
    nullify(p_reference)

    !condition structure:
    call json_value_create(p_condition)             !an object
    call to_object(p_condition,'condition')
    call json_value_add(p_root, p_condition)
    call json_value_add(        p_condition, 'alpha',           t%alpha*180.0/pi)
    call json_value_add(        p_condition, 'beta',            t%beta*180.0/pi)
    call json_value_add(        p_condition, 'ground',          t%hag)
    nullify(p_condition)

    !solver structure:
    call json_value_create(p_solver)             !an object
    call to_object(p_solver,'solver')
    call json_value_add(p_root, p_solver)
    call json_value_add(        p_solver, 'type',            trim(solver))
    call json_value_add(        p_solver, 'convergence',     jacobian_converged)
    call json_value_add(        p_solver, 'relaxation',      jacobian_omega)
    nullify(p_solver)

    !airfoil database structure:
    call json_value_add(        p_root, 'airfoil_DB',            trim(DB_Airfoil))

    !run structure:
    call json_value_create(p_run)             !an object
    call to_object(p_run,'run')
    call json_value_add(p_root, p_run)
    call json_value_add(        p_run, trim(run_type), '')
!fix this
!    if(trim(run_type).eq."target") then
!        call json_value_add(        p_run, 'output',     jacobian_converged)
!        call json_value_add(        p_run, 'value',      jacobian_omega)
!    end if
    nullify(p_run)

    !wings structure:
    call json_value_create(p_wings)             !an object
    call to_object(p_wings,'wings')
    call json_value_add(p_root, p_wings)

    iwing = 0
    do while (iwing < t%nrealwings)
        iwing = iwing + 1
        call json_value_create(p_wingi)
        call to_object(p_wingi,trim(t%wings(iwing)%name))
        call json_value_add(p_wings,p_wingi)
        call json_value_add(        p_wingi, 'ID',            t%wings(iwing)%ID)
        call json_value_add(        p_wingi, 'side',     trim(t%wings(iwing)%orig_side))

        call json_value_create(p_connect)
        call to_object(p_connect,'connect')
        call json_value_add(p_wingi,p_connect)
        call json_value_add(p_connect, 'ID', t%wings(iwing)%connectid)
        call json_value_add(p_connect, 'location', trim(t%wings(iwing)%connectend))
        call json_value_add(p_connect, 'dx', t%wings(iwing)%doffset(1))
        call json_value_add(p_connect, 'dy', t%wings(iwing)%doffset(2))
        call json_value_add(p_connect, 'dz', t%wings(iwing)%doffset(3))
        nullify(p_connect)

        call json_value_add(        p_wingi, 'span',             t%wings(iwing)%span)
        call json_value_add(        p_wingi, 'sweep',            t%wings(iwing)%sweep*180.0/pi)
        call json_value_add(        p_wingi, 'dihedral',         t%wings(iwing)%dihedral*180.0/pi)
        call json_value_add(        p_wingi, 'mounting_angle',   t%wings(iwing)%mount*180.0/pi)
        call json_value_add(        p_wingi, 'washout',          t%wings(iwing)%washout*180.0/pi)
        call json_value_add(        p_wingi, 'root_chord',       t%wings(iwing)%chord_1)
        call json_value_add(        p_wingi, 'tip_chord',        t%wings(iwing)%chord_2)
!        call json_value_add(        p_wingi, 'root_airfoil',     trim(t%wings(iwing)%af1_text))
!        call json_value_add(        p_wingi, 'tip_airfoil',      trim(t%wings(iwing)%af2_text)) !Fix this
        call json_value_add(        p_wingi, 'grid',             t%wings(iwing)%nSec)

        !Write Distribution Files
        if(t%wings(iwing)%f_chord.ne.'none')    call json_value_add( p_wingi, 'chord_file',     trim(t%wings(iwing)%f_chord))
        if(t%wings(iwing)%f_sweep.ne.'none')    call json_value_add( p_wingi, 'sweep_file',     trim(t%wings(iwing)%f_sweep))
        if(t%wings(iwing)%f_dihedral.ne.'none') call json_value_add( p_wingi, 'dihedral_file',  trim(t%wings(iwing)%f_dihedral))
        if(t%wings(iwing)%f_washout.ne.'none')  call json_value_add( p_wingi, 'washout_file',   trim(t%wings(iwing)%f_washout))
        if(t%wings(iwing)%f_af_ratio.ne.'none')  call json_value_add( p_wingi, 'af_ratio_file',   trim(t%wings(iwing)%f_af_ratio))
        nullify(p_wingi)
        if(t%wings(iwing)%orig_side .eq. 'both') iwing = iwing + 1
    end do
    nullify(p_wings)

    !Write File
    open(newunit=iunit, file=filename, status='REPLACE')
    call json_print(p_root,iunit)
    close(iunit)
    call json_destroy(p_root)

end subroutine plane_write_json_file

!-----------------------------------------------------------------------------------------------------------
subroutine plane_init_setup(t)
    type(plane_t) :: t
    integer :: i,iwing,isec,idc,jwing
    real :: start(3)
    type(section_t), pointer :: seci

!    write(*,*) 'Setting up case'
    call plane_set_control_deflections(t,0)

    do iwing=1,t%nwings
        t%wings(iwing)%connect_actual = 0
        idc = t%wings(iwing)%connectid
        if(idc > 0) then
            jwing = 1
            do while (idc .ne. t%wings(jwing)%ID)
                jwing = jwing + 1
            end do
            if((t%wings(iwing)%side == 'left') .and. (t%wings(jwing)%orig_side == 'both')) then
                jwing = jwing + 1
            end if
            t%wings(iwing)%connect_actual = jwing
            if(t%wings(iwing)%connectend.eq.'root') then
                start = t%wings(jwing)%root
            else
                start = t%wings(jwing)%tip
            end if
        else
            start = 0.0
        end if
        call wing_setup(t%wings(iwing),start)
    end do

    !flip ground plane wings
    if(t%groundplane .eq. 1) then
        do iwing=t%nrealwings+1,t%nwings
            call wing_flip_groundplane(t%wings(iwing),t%alpha,t%CG,t%hag) !fix this. Cannot do imbedded derivatives
        end do
    end if

    t%nSize = 0
    do iwing=1,t%nwings
        t%nSize = t%nSize + t%wings(iwing)%nSec
    end do
    call plane_allocate(t)
    i = 1
    do iwing=1,t%nwings
        do isec = 1,t%wings(iwing)%nSec
            t%sec(i)%myp => t%wings(iwing)%sec(isec)
            seci => t%sec(i)%myp
            i = i + 1
        end do
    end do

!    call plane_write_attributes(t)
    solution_exists = 0

end subroutine plane_init_setup

!-----------------------------------------------------------------------------------------------------------
subroutine plane_set_control_deflections(t,include_sections)
    type(plane_t) :: t
    type(json_value),pointer :: j_this, j_mix
    type(section_t),pointer :: si
    integer :: include_sections ! 0 = no. 1 = yes. Must be 0 until wing is allocated
    integer :: imix,iwing,isec,nmix,icontrol
    real :: mydeflection, ratio, adddeflection

    !this subroutine uses the current control deflections stored in t%controls, and mixes them with each
    !wing surface by reading the json input file "mix" parameter for each wing.

    if(t%verbose.eq.1) write(*,*) '--------- Setting Control Deflections ---------'
    do iwing=1,t%nwings
        if(t%verbose.eq.1) write(*,*) 'Wing: ',trim(t%wings(iwing)%name), ' (',trim(t%wings(iwing)%side), ' side)'
        mydeflection = 0.0
        if(t%wings(iwing)%has_control_surface.eq.1) then
            call t%json%get('wings.'//trim(t%wings(iwing)%name)//'.control.mix',  j_this); call json_check()
            nmix = json_value_count(j_this)
            if(t%verbose.eq.1) write(*,*) '    Mixing with ',nmix, ' controls.'
            do imix=1,nmix
                call json_value_get(j_this,imix,j_mix)
                do icontrol = 1, t%ncontrols
                    if(trim(j_mix%name).eq.trim(t%controls(icontrol)%name)) then
                        call myjson_get(t%json, 'wings.'//trim(t%wings(iwing)%name)//'.control.mix.'//trim(j_mix%name), ratio)
!                        call t%json%get('wings.'//trim(t%wings(iwing)%name)//'.control.mix.'//trim(j_mix%name), ratio);
!                        call json_check()
                        if(t%controls(icontrol)%is_symmetric.eq.1) then
                            if(trim(t%wings(iwing)%side).eq.'right') then
                                adddeflection = t%controls(icontrol)%deflection * ratio
                            else
                                adddeflection = t%controls(icontrol)%deflection * ratio
                            end if
                        else
                            if(trim(t%wings(iwing)%side).eq.'right') then
                                adddeflection = t%controls(icontrol)%deflection * ratio
                            else
                                adddeflection = - t%controls(icontrol)%deflection * ratio
                            end if
                        end if
                        if(t%verbose.eq.1) write(*,*) '        ',trim(j_mix%name), ' : ratio = ',ratio, ', deflection = ',&
                                    &adddeflection*180.0/pi, ' (deg)'
                        mydeflection = mydeflection + adddeflection
                    end if
                end do
            end do
        end if

        t%wings(iwing)%control_deflection = mydeflection
        if(t%verbose.eq.1) write(*,*) '    Total Control Deflection = ', mydeflection*180.0/pi, ' (deg)'

        ! Set wing and appropriate sections to this deflection
        if(include_sections.eq.1) then
            do isec = 1,t%wings(iwing)%nSec
                si => t%wings(iwing)%sec(isec)
                if((si%percent_c > t%wings(iwing)%control_span_root) .and. &
                  &(si%percent_c < t%wings(iwing)%control_span_tip)) then
                    si%control_deflection = mydeflection
                else
                    si%control_deflection = 0.0
                end if
            end do
        end if

    end do
    if(t%verbose.eq.1) write(*,*) '-------------------------------------'

end subroutine plane_set_control_deflections

!-----------------------------------------------------------------------------------------------------------
subroutine plane_run_current(t)
    type(plane_t) :: t
    integer :: i
    real :: sumGamma

    if(t%verbose.eq.1) then
        write(*,*) '----------- Case Description ----------'
        write(*,*) '        alpha (deg) = ',t%alpha*180.0/pi
        write(*,*) '        beta (deg) = ',t%beta*180.0/pi
        do i=1,t%ncontrols
            write(*,*) '        ',trim(t%controls(i)%name), ' (deg) = ',t%controls(i)%deflection*180.0/pi
        end do
        write(*,*) '---------------------------------------'
    end if

    sumGamma = sum(t%Gammas(:))
    if(sumGamma.ne.sumGamma) then
        t%Gammas = 0.0
        solution_exists = 0
    end if
    call plane_case_setup(t)
    call plane_solve(t)
    call plane_global_forces(t)
end subroutine plane_run_current

!-----------------------------------------------------------------------------------------------------------
subroutine plane_case_setup(t)
    type(plane_t) :: t
    integer :: i,j
    call plane_set_Uinf(t)
    call plane_set_control_deflections(t,1)
    do i=1,t%nSize
        do j=1,t%nSize
            call plane_influence(t,i,j,t%sec(j)%myp%PC(:),t%vij(i,j,:))
        end do
    end do
end subroutine plane_case_setup

!-----------------------------------------------------------------------------------------------------------
subroutine plane_set_Uinf(t)
    type(plane_t) :: t
    real :: denom
    denom = sqrt(1.0 - sin(t%alpha)**2 * sin(t%beta)**2)
    t%Uinf(1) = -(cos(t%alpha) * cos(t%beta))/denom
    t%Uinf(2) = -(cos(t%alpha) * sin(t%beta))/denom
    t%Uinf(3) = -(sin(t%alpha) * cos(t%beta))/denom
!old code - assumes all angles are less than 90 deg
!    denom = sqrt(1.0 + tan(t%alpha)**2 + tan(t%beta)**2)
!    t%Uinf(1) = -1.0/denom
!    t%Uinf(2) = -tan(t%beta)/denom
!    t%Uinf(3) = -tan(t%alpha)/denom
end subroutine plane_set_Uinf

!-----------------------------------------------------------------------------------------------------------
subroutine plane_influence(t,i,j,P,ans) !Influence of horseshoe i on point P (pass in j=0 if external point)
    type(plane_t) :: t
    integer :: i,j
    real :: P(3),ans(3),r1(3),r2(3),mr1,mr2,vec1(3),vec2(3),vec(3),denom

    r1 = P(:) - t%sec(i)%myp%P1(:)
    r2 = P(:) - t%sec(i)%myp%P2(:)
    mr1 = sqrt(dot_product(r1,r1))
    mr2 = sqrt(dot_product(r2,r2))
    call math_cross_product(t%Uinf(:),r1(:),vec1(:)) !trailing vortices aligned with freestream
    call math_cross_product(t%Uinf(:),r2(:),vec2(:))
    call math_cross_product(r1(:),r2(:),vec(:))
    ans(:) = vec2(:)/(mr2*(mr2-dot_product(t%Uinf,r2))) - vec1(:)/(mr1*(mr1-dot_product(t%Uinf,r1)))
    denom = (mr1*mr2*(mr1*mr2+dot_product(r1,r2)))
    if(i.ne.j) ans(:) = ans(:) + (mr1+mr2)*vec(:)/denom
    ans(:) = t%sec(i)%myp%chord_c/(4.0*pi)*ans(:)
end subroutine plane_influence

!-----------------------------------------------------------------------------------------------------------
subroutine plane_solve(t)
    type(plane_t) :: t
    integer :: i
    REAL :: time1,time2
    call cpu_time(time1)

    if((solver.eq.'linear') .or. (solution_exists.eq.0)) call plane_solve_linear(t)
    if(solver.eq.'nonlinear') call plane_solve_Jacobian(t)

    !store gammas in section info for distributed loads
    do i=1,t%nSize
        t%sec(i)%myp%Gamma = t%Gammas(i)
!write(*,*) i,t%Gammas(i)
    end do

    solution_exists = 1

    call cpu_time(time2)
    if(t%verbose.eq.1) write(*,*) 'CPU time to solve (sec): ',time2-time1
end subroutine plane_solve

!-----------------------------------------------------------------------------------------------------------
subroutine plane_solve_linear(t)
    type(plane_t) :: t
    type(section_t),pointer :: si
    integer :: i,j
    real :: vec(3),CLa

    t%Gammas = 0.0
    if(t%verbose.eq.1) write(*,*) 'Running the linear solver.'
    call plane_update_alphas(t)
    do i=1,t%nSize
        si => t%sec(i)%myp
        CLa = sec_CLa(si)
        do j=1,t%nSize
            t%Amat(i,j) = -CLa*dot_product(t%vij(j,i,:),si%un(:))
        end do
        call math_cross_product(t%Uinf(:),si%zeta(:),vec(:))
        t%Amat(i,i) = t%Amat(i,i) + 2.0*math_mag(3,vec)
!write(*,*) i,t%Amat(i,i)
!write(*,*) si%PC(:)
        t%Bvec(i) = sec_CL(si) !this is slightly different than in the paper, but is the way Phillips actually does it
        !And it appears to be correct to me. The way in the paper uses linear and small-angle approximations in order to
        !simplify the expression for CL of the section. But this way just looks it up directly so can account for
        !non-linearities and more accurate flap deflection.
!write(*,*) 'Bvec ',i,t%Bvec(i), 'alpha ',si%alpha
    end do
!    call math_matinv(t%nSize,t%Amat,t%Ainv) !This is very slow!
!    t%Gammas = matmul(t%Ainv,t%Bvec)
    call math_snyder_ludcmp(t%Amat,t%nSize)
    call math_snyder_lusolv(t%Amat,t%Bvec,t%Gammas,t%nSize)
end subroutine plane_solve_linear

!-----------------------------------------------------------------------------------------------------------
subroutine plane_update_alphas(t)
    type(plane_t) :: t
    type(section_t),pointer :: si
    integer :: i,j
    real :: vec(3),Rcg(3)
    do i=1,t%nSize
        si => t%sec(i)%myp
        Rcg(:) = si%PC(:) - t%CG
        call math_cross_product(Rcg(:),t%omega(:),vec(:))
        vec(:) = vec(:) + t%Uinf(:)
        si%viom = math_mag(3,vec)**2    ! This is used to scale the local forces by the local dynamic pressure.
                                        ! Does not include downwash.
                                        ! This is how Phillips does it, although I would like to check it someday. !Fix This?
        do j=1,t%nSize
            vec(:) = vec(:) + t%vij(j,i,:)*t%Gammas(j)
        end do
        si%v = vec
        si%alpha = atan2(dot_product(vec(:),si%un(:)),dot_product(vec(:),si%ua(:)))
    end do
end subroutine plane_update_alphas

!-----------------------------------------------------------------------------------------------------------
subroutine plane_solve_jacobian(t)
    type(plane_t) :: t
    type(section_t),pointer :: si
    integer :: i,j,iter
    real :: w(3),vn,va,vi(3),vec(3),error,CLa,dGamma(t%nSize)
    110 FORMAT (1X, I10, 100ES25.16)

    if(t%verbose.eq.1) write(*,*) 'Running the Jacobian solver.'
    if(t%verbose.eq.1) write(*,*) '---------------- Solver Settings ----------------'
    if(t%verbose.eq.1) write(*,*) 'Convergence Criteria = ',jacobian_converged
    if(t%verbose.eq.1) write(*,*) '   Relaxation Factor = ',jacobian_omega
    if(t%verbose.eq.1) write(*,*) '-------------------------------------------------'
    if(t%verbose.eq.1) write(*,*)
    if(t%verbose.eq.1) write(*,*) ' iteration   residual'
    iter = 0
    error = 100.0
    do while (error > jacobian_converged)
        iter = iter + 1
        call plane_update_alphas(t)
        do i=1,t%nSize
            si => t%sec(i)%myp
            vi = si%v
            call math_cross_product(vi(:),si%zeta(:),w(:))
            vn = dot_product(vi(:),si%un(:))
            va = dot_product(vi(:),si%ua(:))
            CLa = sec_CLa(si)
            do j=1,t%nSize
                call math_cross_product(t%vij(j,i,:),si%zeta(:),vec(:))
                t%Amat(i,j) = 2.0*dot_product(w,vec)/math_mag(3,w)*t%Gammas(i) &
                          & - CLa*(va*dot_product(t%vij(j,i,:),si%un(:)) &
                          & - vn*dot_product(t%vij(j,i,:),si%ua(:)))/(va**2+vn**2)
            end do
            t%Amat(i,i) = t%Amat(i,i) + 2.0*math_mag(3,w)
            t%Bvec(i) = 2.0*math_mag(3,w)*t%Gammas(i) - sec_CL(si)
        end do
        call math_snyder_ludcmp(t%Amat,t%nSize)
        call math_snyder_lusolv(t%Amat,-t%Bvec,dGamma,t%nSize)
!        call math_matinv(t%nSize,t%Amat,t%Ainv)
!        dGamma = matmul(t%Ainv,-t%Bvec)
        error = sqrt(dot_product(t%Bvec,t%Bvec))
        t%Gammas(:) = t%Gammas(:) + jacobian_omega*dGamma(:)
        if(t%verbose.eq.1) write(*,110) iter,error
        if(iter > nonlinear_maxiter) then
            exit
        end if
    end do
    if(t%verbose.eq.1) write(*,*) 'Solution successfully converged.'
end subroutine plane_solve_jacobian

!-----------------------------------------------------------------------------------------------------------
subroutine plane_global_forces(t)
    type(plane_t) :: t
    type(section_t),pointer :: si
    type(json_value),pointer    :: p_root, p_forcetype, p_object
    character(100) :: filename
    integer :: iwing,isec,i,j,numwings,iunit
    real :: vi(3),vec(3),vm(3),ui(3),vmp(3),Rcg(3)
!    real :: test_CL, test_CD, test_CN, test_CA, test_fvec(3), test_mvec(3)
    120 Format(A15, A5, 100ES25.13)

    numwings = t%nrealwings !only sum over real wings, not reflected wings


    t%GF = 0.0
    t%GM = 0.0
    t%GL = 0.0
    t%GD = 0.0
    t%GLoD = 0.0

    call plane_update_alphas(t)
    i = 1 !needed to reference gammas
    do iwing=1,numwings
        do isec=1,t%wings(iwing)%nSec
            si => t%wings(iwing)%sec(isec)
            vi = si%v
            ui(:) = vi(:)/math_mag(3,vi)
            Rcg(:) = si%PC(:) - t%CG(:)

            call math_cross_product(vi(:),si%zeta(:),vec(:))
            call math_cross_product(Rcg(:),vec(:),vm(:))
            call math_cross_product(Rcg(:),ui(:),vmp(:))

            !Calculate force and moment from vortex theory of lift
            t%GF(iwing,1,:) = t%GF(iwing,1,:) + 2.0*t%Gammas(i)*vec(:)*si%ds/t%Sr*si%viom
            t%GF(iwing,2,:) = t%GF(iwing,2,:) + sec_CD(si)*si%ds/t%Sr*ui(:)*si%viom    !Parasitic
            t%GM(iwing,1,:) = t%GM(iwing,1,:) + (2.0*t%Gammas(i)*vm(:) - sec_Cm(si)*si%chord_c*si%us(:))*si%ds/t%Sr*si%viom
            t%GM(iwing,2,:) = t%GM(iwing,2,:) + sec_CD(si)*si%ds/t%Sr*vmp(:)*si%viom   !Parasitic

            ! Store local force and moment for later.
            si%F(:) = (2.0*t%Gammas(i)*vec(:) + sec_CD(si)*ui(:))*si%ds*si%viom
            si%M(:) = - sec_Cm(si)*si%chord_c*si%us(:)*si%ds*si%viom

            !Calculate force and moment from section theory of lift
!            test_CL = sec_CL(si)*si%ds/t%Sr*si%viom
!            test_CD = sec_CD(si)*si%ds/t%Sr*si%viom
!            test_CN = test_CL*cos(si%alpha) + test_CD*sin(si%alpha)
!            test_CA = test_CD*cos(si%alpha) - test_CL*sin(si%alpha)
!            test_fvec(:) = test_CN*si%un(:) + test_CA*si%ua(:)
!            t%GF(iwing,1,:) = t%GF(iwing,1,:) + test_fvec(:)
!            t%GF(iwing,2,:) = t%GF(iwing,2,:) + sec_CD(si)*si%ds/t%Sr*ui(:)*si%viom    !Parasitic
!            call math_cross_product(Rcg(:),test_fvec(:),test_mvec(:))
!            t%GM(iwing,1,:) = t%GM(iwing,1,:) + test_mvec - sec_Cm(si)*si%chord_c*si%us(:)*si%ds/t%Sr*si%viom
!            t%GM(iwing,2,:) = t%GM(iwing,2,:) + sec_CD(si)*si%ds/t%Sr*vmp(:)*si%viom   !Parasitic


            i = i + 1
        end do
        t%GF(iwing,3,:) = t%GF(iwing,1,:) + t%GF(iwing,2,:) !total Forces for that wing
        t%GM(iwing,3,:) = t%GM(iwing,1,:) + t%GM(iwing,2,:) !total Moment for that wing

        t%GF(numwings+1,:,:) = t%GF(numwings+1,:,:) + t%GF(iwing,:,:) !totals
        t%GM(numwings+1,:,:) = t%GM(numwings+1,:,:) + t%GM(iwing,:,:) !totals
    end do

    !calculate actual wing area
!    area = 0.0
!    do i=1,t%nrealwings
!        area = area + t%wings(i)%area
!    end do
!    F = F*t%Sr/area
!    M = M*t%Sr/area


    t%GM(:,:,1) = t%GM(:,:,1)/t%lat_r
    t%GM(:,:,2) = t%GM(:,:,2)/t%long_r
    t%GM(:,:,3) = t%GM(:,:,3)/t%lat_r

    t%GL(:,:) = -t%GF(:,:,3)*cos(t%alpha) + t%GF(:,:,1)*sin(t%alpha) !fix this? should I rotate lift and drag by beta as well?
    t%GD(:,:) = -t%GF(:,:,1)*cos(t%alpha) - t%GF(:,:,3)*sin(t%alpha)
    t%GLoD = t%GL(:,:)/t%GD(:,:)

    !Remove NANs from LoD calcs
    do i=1,numwings+1
        do j=1,3
            if(t%GLoD(i,j).ne.t%GLoD(i,j)) t%GLoD(i,j) = 0.0
        end do
    end do

    if(t%verbose.eq.1) then
        !write forces to the screen
        write(*,*)
        write(*,*) 'Contributions to Forces and Moments at CG ------------------------'
        write(*,*) 'Wing Name     side       CX                       CY                       CZ                       &
                   &Cl                       Cm                       Cn                       &
                   &CL                       CD                       L/D'
        write(*,*)
        write(*,*) 'Inviscid'
        do iwing=1,numwings
            write(*,120) t%wings(iwing)%name,t%wings(iwing)%side,t%GF(iwing,1,:),t%GM(iwing,1,:),t%GL(iwing,1),t%GD(iwing,1),&
                    &t%GLoD(iwing,1)
        end do
        write(*,120) t%name,'',t%GF(numwings+1,1,:),t%GM(numwings+1,1,:),t%GL(numwings+1,1),t%GD(numwings+1,1),t%GLoD(numwings+1,1)

        write(*,*)
        write(*,*) 'Viscous'
        do iwing=1,numwings
            write(*,120) t%wings(iwing)%name,t%wings(iwing)%side,t%GF(iwing,2,:),t%GM(iwing,2,:),t%GL(iwing,2),t%GD(iwing,2),&
                    &t%GLoD(iwing,2)
        end do
        write(*,120) t%name,'',t%GF(numwings+1,2,:),t%GM(numwings+1,2,:),t%GL(numwings+1,2),t%GD(numwings+1,2),t%GLoD(numwings+1,2)

        write(*,*)
        write(*,*) 'Total'
        do iwing=1,numwings
            write(*,120) t%wings(iwing)%name,t%wings(iwing)%side,t%GF(iwing,3,:),t%GM(iwing,3,:),t%GL(iwing,3),t%GD(iwing,3),&
                    &t%GLoD(iwing,3)
        end do
        write(*,120) t%name,'',t%GF(numwings+1,3,:),t%GM(numwings+1,3,:),t%GL(numwings+1,3),t%GD(numwings+1,3),t%GLoD(numwings+1,3)
        write(*,*) '------------------------------------------------------------------'


        !Write JSON File
        !root
        filename = trim(adjustl(t%master_filename))//'_forces.json'
        call json_value_create(p_root)           ! create the value and associate the pointer
        call to_object(p_root,trim(filename))    ! add the file name as the name of the overall structure

        do i=1,3
            call json_value_create(p_forcetype)             !an object
            select case (i)
                case (1)
                    call to_object(p_forcetype,'inviscid')
                case (2)
                    call to_object(p_forcetype,'viscous')
                case (3)
                    call to_object(p_forcetype,'total')
            end select
            call json_value_add(p_root, p_forcetype)
            iwing = 0
            do iwing =1, t%nrealwings+1
                call json_value_create(p_object)
                if(iwing>t%nrealwings) then
                    call to_object(p_object,trim(t%name))
                else
                    call to_object(p_object,trim(t%wings(iwing)%name)//'_'//trim(t%wings(iwing)%side))
                end if
                call json_value_add(p_forcetype,p_object)
                call plane_add_json_force(p_object,t%GF(iwing,i,:),t%GM(iwing,i,:),t%GL(iwing,i),t%GD(iwing,i),t%GLoD(iwing,i))
                nullify(p_object)
            end do
            nullify(p_forcetype)
        end do

        !Write File
    !    write(*,'(A)') 'writing file '//trim(filename)//'...'
        open(newunit=iunit, file=filename, status='REPLACE')
        call json_print(p_root,iunit)
        close(iunit)
        call json_destroy(p_root)
    !    write(*,'(A)') 'done.'
        write(*,'(A)') 'Force and moment results written to: '//trim(filename)
    end if
end subroutine plane_global_forces

!-----------------------------------------------------------------------------------------------------------
subroutine plane_add_json_force(this,F,M,L,D,LoD)
    type(json_value),pointer    :: this
    real :: F(3), M(3), L,D,LoD

    call json_value_add(        this, 'CX',            F(1))
    call json_value_add(        this, 'CY',            F(2))
    call json_value_add(        this, 'CZ',            F(3))

    call json_value_add(        this, 'Cl',            M(1))
    call json_value_add(        this, 'Cm',            M(2))
    call json_value_add(        this, 'Cn',            M(3))

    call json_value_add(        this, 'CL',            L)
    call json_value_add(        this, 'CD',            D)
    call json_value_add(        this, 'L/D',           LoD)

end subroutine plane_add_json_force
!-----------------------------------------------------------------------------------------------------------
subroutine plane_distributed_loads(t)
    type(plane_t) :: t
    type(section_t),pointer :: si
    character(100) :: filename
    integer :: iwing,isec,i,numwings,ierror,jwing,iforce,ipoint,min_isec
    real :: vi(3),vec(3),F_tip(3),F_root(3),M_tip(3),M_root(3),prevP(3),prevF(3),prevM(3),myF(3),myM(3)
    real :: min_dist,force_dir(3),force_pos(3),force_mag,dist

    numwings = t%nrealwings !only sum over real wings, not reflected wings

    !clear all Loads variables
    do isec=1,t%nSize
        t%sec(isec)%myp%Load_F(:) = 0.0
        t%sec(isec)%myp%Load_M(:) = 0.0
    end do
    !tie in external forces to closest structure !this should be same as in subroutine for viewing
    min_isec = 0
    do iforce = 1,size(t%external_forces)
        do ipoint=1,t%external_forces(iforce)%datasize
            force_pos(:) = t%external_forces(iforce)%RawData(ipoint,1:3)
            force_dir(:) = t%external_forces(iforce)%RawData(ipoint,4:6)
            force_mag    = t%external_forces(iforce)%RawData(ipoint,7)
            min_dist = 1.0e16
            do isec=1,t%nSize
                si => t%sec(isec)%myp
                dist = math_length(3,force_pos(:),si%PC(:))
                if(dist < min_dist) then
                    min_dist = dist
                    min_isec = isec
                end if
            end do
            si => t%sec(min_isec)%myp
            si%Load_F(:) = si%Load_F(:) + force_dir(:)*force_mag
            call math_cross_product(force_pos(:)-si%PC(:),force_dir(:)*force_mag,vec(:))
            si%Load_M(:) = si%Load_M(:) + vec(:)
        end do
    end do

    call plane_update_alphas(t)
    i = 1 !needed to reference gammas
    do iwing=numwings,1,-1
        F_tip = 0.0
        F_root = 0.0
        M_tip = 0.0
        M_root = 0.0
        !find connecting forces and moments
        do jwing = 1,numwings
            if(t%wings(jwing)%connect_actual .eq. iwing) then
                if(t%wings(jwing)%connectend .eq. 'tip') then
                    F_tip = F_tip + t%wings(jwing)%root_F
                    call math_cross_product(t%wings(jwing)%root(:)-t%wings(iwing)%tip(:),t%wings(jwing)%root_F(:),vec(:))
                    M_tip = M_tip + t%wings(jwing)%root_M + vec
                else
                    F_root = F_root + t%wings(jwing)%root_F
                    call math_cross_product(t%wings(jwing)%root(:)-t%wings(iwing)%root(:),t%wings(jwing)%root_F(:),vec(:))
                    M_root = M_root + t%wings(jwing)%root_M + vec
                end if
            end if
        end do
        t%wings(iwing)%tip_F = F_tip
        t%wings(iwing)%tip_M = M_tip
        prevP = t%wings(iwing)%tip
        prevF = F_tip
        prevM = M_tip
        do isec=t%wings(iwing)%nSec,1,-1
            si => t%wings(iwing)%sec(isec)
            vi = si%v
            call math_cross_product(vi(:),si%zeta(:),vec(:))
            myF = 2.0*si%Gamma*vec(:)*si%ds !ignores parasitic force Fix this
            myM = - sec_Cm(si)*si%chord_c*si%us(:)*si%ds !ignores parasitic moment Fix this!
            call math_cross_product(prevP(:)-si%PC(:),prevF(:),vec(:))
            si%Load_F = si%Load_F + myF + prevF
            si%Load_M = si%Load_M + myM + prevM + vec
            prevP = si%PC
            prevF = si%Load_F
            prevM = si%Load_M
        end do
        call math_cross_product(prevP(:)-t%wings(iwing)%root(:),prevF(:),vec(:))
        t%wings(iwing)%root_F = F_root + prevF
        t%wings(iwing)%root_M = M_root + prevM + vec
    end do

    !write forces to a file
    150 Format(9ES26.16)
    filename = trim(adjustl(t%master_filename))//'_distributed.txt'
    open(unit = 10, File = filename, action = 'write', iostat = ierror)
    write(10,*) 'x y z Fx Fy Fz Mx My Mz'
    do iwing=1,numwings
        write(10,150) t%wings(iwing)%root,t%wings(iwing)%root_F,t%wings(iwing)%root_M
        do isec=1,t%wings(iwing)%nSec
            si => t%wings(iwing)%sec(isec)
            write(10,150) si%PC(:),si%Load_F(:),si%Load_M(:)
        end do
        write(10,150) t%wings(iwing)%tip,t%wings(iwing)%tip_F,t%wings(iwing)%tip_M
    end do
    close(10)

end subroutine plane_distributed_loads

!-----------------------------------------------------------------------------------------------------------
subroutine plane_loads_old(t)
    type(plane_t) :: t
    character(100) :: filename
    integer :: iwing,isec,ierror,iforce
    type(section_t),pointer :: si
    type(dataset_t) :: dist
    real :: rawdata(30,2)
    real :: A(3,3),Atemp(3,3),X(3),B(3),span,span1,span2,span3,F0(3),F1(3),P(3),percent
    real :: zero
    120 Format(A15, 100ES25.13)

    zero = 0.0

    filename = trim(adjustl(t%master_filename))//'_loads.txt'
    open(unit = 10, File = filename, action = 'write', iostat = ierror)
    write(10,*) ' WingName           ControlPoint(x)          ControlPoint(y)          ControlPoint(z)          &
           &SpanwiseCoordinate       Chord                    &
           &F(x)/q/span              F(y)/q/span              F(z)/q/span'

    do iwing=1,1
        do isec=1,t%wings(iwing)%nSec
            si => t%wings(iwing)%sec(isec)
            span = (si%percent_2-si%percent_1)*t%wings(iwing)%span
            write(*,*) isec,si%percent_1*t%wings(iwing)%span,si%F(3)/span
            write(*,*) isec,si%percent_2*t%wings(iwing)%span,si%F(3)/span
            rawdata(isec,1) = si%percent_c*t%wings(iwing)%span
            rawdata(isec,2) = si%F(3)/span
        end do
    end do

    call ds_create_from_data(dist,30,2,rawdata)

    call ds_cubic_setup(dist, 1, 2, zero, 2, zero)
    call ds_print_data(dist)

    do iwing=1,1
        !Solve for end point by fitting parabola to first three points
        si => t%wings(iwing)%sec(1)
        span1 = abs(si%percent_2-si%percent_1)*t%wings(iwing)%span
        A(1,1) = (si%percent_c)**2
        A(1,2) = (si%percent_c)
        A(1,3) = 1.0

        si => t%wings(iwing)%sec(2)
        span2 = abs(si%percent_2-si%percent_1)*t%wings(iwing)%span
        A(2,1) = (si%percent_c)**2
        A(2,2) = (si%percent_c)
        A(2,3) = 1.0

        si => t%wings(iwing)%sec(3)
        span3 = abs(si%percent_2-si%percent_1)*t%wings(iwing)%span
        A(3,1) = (si%percent_c)**2
        A(3,2) = (si%percent_c)
        A(3,3) = 1.0

        Atemp(:,:) = A(:,:)
        do iforce=1,3
            B(1) = t%wings(iwing)%sec(1)%F(iforce)/span1
            B(2) = t%wings(iwing)%sec(2)%F(iforce)/span2
            B(3) = t%wings(iwing)%sec(3)%F(iforce)/span3
            A(:,:) = Atemp(:,:)
            call math_snyder_ludcmp(A,3)
            call math_snyder_lusolv(A,B,X,3)
            F0(iforce) = X(3)
        end do

!        F0(3) = -1.78
        !Write end point to file
        si => t%wings(iwing)%sec(1)
        span = (si%percent_2-si%percent_1)*t%wings(iwing)%span
        P(:) = si%P1(:)
        percent = si%percent_1
        if(t%wings(iwing)%side.eq.'left') then
            P(:) = si%P2(:)
            percent = si%percent_2
        end if
        write(10,120) t%wings(iwing)%name, P(:), &
                      & percent*t%wings(iwing)%span, &
                      & si%chord_1, &
                      & F0(:)

call ds_cubic_interpolate(dist, zero, 0, F1(1:2))
write(*,*) 0.0,F1(2)

        !Solve for all other cell vertex points
        do isec=1,t%wings(iwing)%nSec
            si => t%wings(iwing)%sec(isec)
            span = (si%percent_2-si%percent_1)*t%wings(iwing)%span
            F1(:) = 2.0*si%F(:)/span - F0(:)

            P(:) = si%P2(:)
            percent = si%percent_2
            if(t%wings(iwing)%side.eq.'left') then
                P(:) = si%P1(:)
                percent = si%percent_1
            end if

            call ds_cubic_interpolate(dist,percent*t%wings(iwing)%span,0,F1(1:2))


write(*,*) percent*t%wings(iwing)%span,F1(2)
!write(*,*) isec,si%percent_c*t%wings(iwing)%span,si%F(3)/span,si%percent_1*t%wings(iwing)%span,&
!& span,F0(3),F1(3),si%F(3)

            write(10,120) t%wings(iwing)%name, P(:), &
                          & percent*t%wings(iwing)%span, &
                          & si%chord_2, &
                          & F1(:)
            F0(:) = F1(:)
        end do
        write(10,*)
    end do
    close(10)
    write(*,'(A)') 'Load results written to: '//trim(filename)

end subroutine plane_loads_old
!-----------------------------------------------------------------------------------------------------------
subroutine plane_loads_orig(t)
    type(plane_t) :: t
    character(100) :: filename
    integer :: iwing,isec,ierror
    type(section_t),pointer :: si
    real :: span
    120 Format(A15, 100ES25.13)

    filename = trim(adjustl(t%master_filename))//'_loads.txt'
    open(unit = 10, File = filename, action = 'write', iostat = ierror)
    write(10,*) ' WingName           ControlPoint(x)          ControlPoint(y)          ControlPoint(z)          &
           &SpanwiseCoordinate       Chord                    Area                     &
           &F(x)/q                   F(y)/q                   F(z)/q                   &
           &F(x)/q/span              F(y)/q/span              F(z)/q/span'
    do iwing=1,t%nwings
        do isec=1,t%wings(iwing)%nSec
            si => t%wings(iwing)%sec(isec)
            span = abs(si%percent_2-si%percent_1)*t%wings(iwing)%span
            write(10,120) t%wings(iwing)%name, 0.5*(si%P1(:)+si%P2(:)), &
                          & 0.5*(si%percent_2+si%percent_1)*t%wings(iwing)%span, &
                          & 0.5*(si%chord_1+si%chord_2), si%ds, &
                          & si%F(:), si%F(:)/span
        end do
        write(10,*)
    end do
    close(10)
    write(*,'(A)') 'Load results written to: '//trim(filename)

end subroutine plane_loads_orig

!-----------------------------------------------------------------------------------------------------------
subroutine plane_write_attributes(t)
    type(plane_t) :: t
    integer :: iwing
    write(*,*) '-------------------------------------------------------------'
    write(*,*) '                 Airplane name : ',t%name
    write(*,*) '               Number of wings : ',t%nrealwings
    write(*,*) '                Reference area : [',t%Sr,']'
    write(*,*) ' Longitudinal reference length : [',t%long_r,']'
    write(*,*) '      Lateral reference length : [',t%lat_r,']'
    write(*,*)
    do iwing = 1,t%nrealwings
        call wing_write_attributes(t%wings(iwing))
    end do

end subroutine plane_write_attributes

end module plane_m
