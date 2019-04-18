module wing_m
#ifdef dnad
    use dnadmod
#define real type(dual)
#endif
    use section_m
    implicit none

    type wing_t
        integer :: ID
        character(20) :: name

        character(5) :: side !left/right/both
        character(5) :: orig_side !original spec from input file
        character(4) :: connectend !root/tip
        integer :: connectid,connect_actual
        integer :: has_control_surface
        integer :: control_is_sealed
        integer :: is_linear
        integer :: sweep_definition

        real :: area
        real :: doffset(3) !offset from connecting wing
        real :: dy !y-offset from connecting wing
        real :: root(3)
        real :: tip(3)
        real :: span
        real :: sweep
        real :: dihedral
        real :: chord_1
        real :: chord_2
        real :: mount
        real :: washout
        real :: control_span_root
        real :: control_span_tip
        real :: control_chord_root
        real :: control_chord_tip
        real :: control_deflection

!        character(100) :: af1_text,af2_text !original text from file
        type(airfoil_p),allocatable,dimension(:) :: airfoils

        integer :: nSec
        integer :: nairfoils
        integer :: root_clustering  !Root clustering parameter (1 = clustered, 0 = linear)
        integer :: tip_clustering  !Tip clustering parameter (1 = clustered, 0 = linear)

        real :: root_F(3), root_M(3) !force and moment load at root of the wing
        real :: tip_F(3), tip_M(3) !force and moment load at the tip of the wing
        real :: root_theta, tip_theta !twisted shape at root and tip


        type(section_t),pointer,dimension(:) :: sec

        !possible files
        character(100) :: f_sweep, f_dihedral, f_washout, f_chord, f_af_ratio
        character(100) :: f_EIGJ, f_elastic_twist
    end type wing_t

contains

!-----------------------------------------------------------------------------------------------------------
subroutine wing_allocate(t)
    type(wing_t) :: t
    allocate(t%sec(t%nSec))
end subroutine wing_allocate

!-----------------------------------------------------------------------------------------------------------
subroutine wing_allocate_airfoils(t)
    type(wing_t) :: t
    allocate(t%airfoils(t%nairfoils))
end subroutine wing_allocate_airfoils

!-----------------------------------------------------------------------------------------------------------
subroutine wing_deallocate(t)
    type(wing_t) :: t
    deallocate(t%sec)
end subroutine wing_deallocate

!-----------------------------------------------------------------------------------------------------------
subroutine wing_setup(t,start)
    type(wing_t) :: t
    type(section_t),pointer :: si
    integer :: i, isec
    REAL :: dtheta
    real :: start(3),qvec(3),nvec(3),avec(3),fvec(3),percent_1,percent_2,percent_c,chord_1,chord_2,RA,span
    real :: my_sweep,my_dihedral,my_twist,temp
    real :: sweep1,sweep2,dihedral1,dihedral2,twist1,twist2
    real :: cfc,thetaf,efi,etah,percent_c_af,percent_1_af,percent_2_af
    type(dataset_t) :: data_sweep,data_dihedral,data_washout,data_chord,data_af_ratio,data_elastic_twist
    logical :: linear
    REAL :: theta_start, theta_end, fac, off

    t%is_linear = 1
    if(t%nairfoils > 2) t%is_linear = 0

    !allocate arrays from files
    if(t%f_chord .ne. 'none') then
        t%is_linear = 0; call ds_create_from_file(data_chord,t%f_chord,2)
    end if
    if(t%f_sweep .ne. 'none') then
        t%is_linear = 0; call ds_create_from_file(data_sweep,t%f_sweep,2)
    end if
    if(t%f_dihedral .ne. 'none') then
        t%is_linear = 0; call ds_create_from_file(data_dihedral,t%f_dihedral,2)
    end if
    if(t%f_washout .ne. 'none') then
        t%is_linear = 0; call ds_create_from_file(data_washout,t%f_washout,2)
    end if
    if(t%f_af_ratio .ne. 'none') then
        t%is_linear = 0; call ds_create_from_file(data_af_ratio,t%f_af_ratio,2)
        do i = 1, data_af_ratio%datasize
            if(data_af_ratio%RawData(i, 2) < 0.0 .or. &
             & data_af_ratio%RawData(i, 2) > t%nairfoils - 1) then
                write(*,*) 'Error: Airfoil ratio file "', trim(t%f_af_ratio), &
                    & '" specifies airfoil ratios that are outside the', &
                    & ' allowable limits of 0 to ', t%nairfoils - 1
                STOP
            end if
        end do
    end if
    if(t%f_elastic_twist .ne. 'none') then
        t%is_linear = 0; call ds_create_from_file(data_elastic_twist,t%f_elastic_twist,2)
    end if

    t%root = start + t%doffset
    t%root(2) = t%root(2) + t%dy
    if(t%side.eq.'left') then
        t%root(2) = t%root(2) - 2.0*t%dy
    end if
    start = t%root

    call wing_allocate(t)

    linear = .false.
    theta_start = 0.0
    theta_end = pi
    off = 1.0
    fac = 0.5
    if(t%root_clustering .ne. 0 .and. t%tip_clustering .ne. 0) then
        theta_start = 0.0
        theta_end = pi
        off = 1.0
        fac = 0.5
    else if(t%root_clustering .ne. 0) then
        theta_start = 0.0
        theta_end = pi / 2.0
        off = 1.0
        fac = 1.0
    else if(t%tip_clustering .ne. 0) then
        theta_start = pi / 2.0
        theta_end = pi
        off = 0.0
        fac = 1.0
    else
        linear = .true.
    end if

    dtheta = (theta_end - theta_start)/REAL(t%nSec)

    t%area = 0.0
    span = 0.0
    do isec=1,t%nSec
        si => t%sec(isec)
        if(linear) then
            percent_1 = REAL(isec - 1) / REAL(t%nSec)
            percent_2 = REAL(isec) / REAL(t%nSec)
            percent_c = (REAL(isec) - 0.5) / REAL(t%nSec)
        else
            percent_1 = fac * (off - cos(theta_start + dtheta * REAL(isec-1)))
            percent_2 = fac * (off - cos(theta_start + dtheta * REAL(isec)))
            percent_c = fac * (off - cos(theta_start + dtheta * (REAL(isec) - 0.5)))
        end if
        if(t%side.eq.'left') then !must handle differently for left wing
            temp = percent_1
            percent_1 = percent_2
            percent_2 = temp
        end if

        !data from wing file
        if(t%chord_2 >= 0.0) then
            chord_1 = t%chord_1 + percent_1*(t%chord_2 - t%chord_1)
            chord_2 = t%chord_1 + percent_2*(t%chord_2 - t%chord_1)
        else !elliptic wing
            RA = 8.0*t%span/pi/t%chord_1
            chord_1 = 8.0*t%span/pi/RA*sqrt(1.0-percent_1**2)
            chord_2 = 8.0*t%span/pi/RA*sqrt(1.0-percent_2**2)
        end if
        my_sweep = t%sweep
        my_dihedral = t%dihedral
        my_twist = t%mount - percent_c*t%washout

        ! % airfoil 0=airfoil1, 1=airfoil2, 2=airfoil3
        percent_c_af = percent_c*REAL(t%nairfoils-1)
        percent_1_af = percent_1*REAL(t%nairfoils-1)
        percent_2_af = percent_2*REAL(t%nairfoils-1)

        !For geometry purposes
        sweep1 = my_sweep
        sweep2 = my_sweep
        dihedral1 = my_dihedral
        dihedral2 = my_dihedral
        twist1 = t%mount - percent_1*t%washout
        twist2 = t%mount - percent_2*t%washout

        ! data from input files
        if(t%f_chord .ne. 'none') then
            chord_1 = ds_linear_interpolate_col(data_chord,percent_1,1,2)
            chord_2 = ds_linear_interpolate_col(data_chord,percent_2,1,2)
        end if
        if(t%f_sweep .ne. 'none') then
            my_sweep = ds_linear_interpolate_col(data_sweep,percent_c,1,2)*pi/180.0
            sweep1   = ds_linear_interpolate_col(data_sweep,percent_1,1,2)*pi/180.0
            sweep2   = ds_linear_interpolate_col(data_sweep,percent_2,1,2)*pi/180.0
        end if
        if(t%f_dihedral .ne. 'none') then
            my_dihedral = ds_linear_interpolate_col(data_dihedral,percent_c,1,2)*pi/180.0
            dihedral1   = ds_linear_interpolate_col(data_dihedral,percent_1,1,2)*pi/180.0
            dihedral2   = ds_linear_interpolate_col(data_dihedral,percent_2,1,2)*pi/180.0
        end if
        if(t%f_washout .ne. 'none') then
            my_twist   = t%mount - ds_linear_interpolate_col(data_washout,percent_c,1,2)*pi/180.0
            twist1     = t%mount - ds_linear_interpolate_col(data_washout,percent_1,1,2)*pi/180.0
            twist2     = t%mount - ds_linear_interpolate_col(data_washout,percent_2,1,2)*pi/180.0
        end if
        if(t%f_af_ratio .ne. 'none') then
            percent_c_af = ds_linear_interpolate_col(data_af_ratio,percent_c,1,2)
            percent_1_af = ds_linear_interpolate_col(data_af_ratio,percent_1,1,2)
            percent_2_af = ds_linear_interpolate_col(data_af_ratio,percent_2,1,2)
        end if
        if(t%f_elastic_twist .ne.  'none') then
            my_twist   = my_twist + ds_linear_interpolate_col(data_elastic_twist,percent_c,1,2)*pi/180.0
        end if

        !for geometry purposes
        si%twist = my_twist
        si%twist1 = twist1
        si%twist2 = twist2
        si%dihedral = my_dihedral
        si%dihedral1 = dihedral1
        si%dihedral2 = dihedral2
        si%sweep = my_sweep

        !operate!
        qvec(1) = 0.0; qvec(2) = 1.0; qvec(3) = 0.0
        nvec(1) = 0.0; nvec(2) = 0.0; nvec(3) =-1.0
        avec(1) =-1.0; avec(2) = 0.0; avec(3) = 0.0
        fvec(1) =-1.0; fvec(2) = 0.0; fvec(3) = 0.0 !we can set up fvec here because it doesn't change for derivatives
        call math_rot_z(qvec,my_sweep)
        call math_rot_x(qvec,-my_dihedral)

        call math_rot_y(nvec,my_twist)
        call math_rot_x(nvec,-my_dihedral)
        call math_rot_y(avec,my_twist)
        call math_rot_x(avec,-my_dihedral)
        call math_rot_y(fvec,my_twist+t%control_deflection)
        call math_rot_x(fvec,-my_dihedral)
        qvec(:) = qvec(:)/cos(my_sweep)
        if(t%side.eq.'left') then !must handle differently for left wing
            qvec(2) = -qvec(2)
            nvec(2) = -nvec(2)
            avec(2) = -avec(2)
            fvec(2) = -fvec(2)
        end if
        si%percent_1 = percent_1
        si%percent_2 = percent_2
        si%percent_c = percent_c
        si%chord_1 = chord_1
        si%chord_2 = chord_2
        si%chord_c = 2.0/3.0*(chord_1**2 + chord_1*chord_2 + chord_2**2)/(chord_1 + chord_2)
        si%ds = abs(0.5*(chord_1+chord_2)*(percent_2-percent_1)*t%span)
        si%un(:) = nvec(:)
        si%ua(:) = avec(:)
        si%uf(:) = fvec(:)
        call math_cross_product(avec(:),nvec(:),si%us)

        if(t%side.eq.'left') then !must handle differently for left wing
            si%P2(:) = start(:)
            si%P1(:) = start(:) + qvec(:)*(percent_1-percent_2)*t%span
        else
            si%P1(:) = start(:)
            si%P2(:) = start(:) + qvec(:)*(percent_2-percent_1)*t%span
        end if
        si%dl(:) = si%P2(:) - si%P1(:)
        si%PC(:) = si%P1(:) + si%dl(:)*(percent_c - percent_1)/(percent_2 - percent_1)

        if(t%sweep_definition .eq. 0) then !remove sweep from control points
            si%temp_P1(:) = si%P1(:)
            si%temp_P2(:) = si%P2(:)

            si%P1(1) = si%PC(1)
            si%P2(1) = si%PC(1)
            si%dl(:) = si%P2(:) - si%P1(:)
        end if

        si%zeta(:) = si%chord_c*si%dl(:)/si%ds
        si%Rroot(:) = si%PC(:) - t%root(:)

        !setup airfoils for section
        call wing_get_airfoil_weighting(t,percent_c_af,si%af_weight_c,si%afc_a,si%afc_b)
!        write(*,*) si%percent_c,percent_c_af,si%af_weight_c,trim(si%afc_a%name),trim(si%afc_b%name)
        call wing_get_airfoil_weighting(t,percent_1_af,si%af_weight_1,si%af1_a,si%af1_b)
        call wing_get_airfoil_weighting(t,percent_2_af,si%af_weight_2,si%af2_a,si%af2_b)

        !Add flaps
        si%cf_c = 0.0
        if(t%has_control_surface .eq. 1) then
            si%control_deflection = t%control_deflection !This is set here, but can be updated later for derivs.
            if((percent_c > t%control_span_root) .and. (percent_c < t%control_span_tip)) then
                cfc = t%control_chord_root + (percent_c - t%control_span_root)/(t%control_span_tip - t%control_span_root)*&
                      & (t%control_chord_tip - t%control_chord_root) !Fix this! Should force straight hinge line
                thetaf = acos(2.0*cfc-1.0)
                efi = 1.0 - (thetaf - sin(thetaf))/pi
                etah = 3.9598*atan((cfc + 0.006527)*89.2574 + 4.898015) - 5.18786
                if(t%control_is_sealed .eq. 0) then
                    etah = 0.8*etah
                end if
                si%cf_c = cfc
                si%ef = etah*efi  !adds in etad later because it is a function of deflection
                si%Cmdelta = 0.25*(sin(2.0*thetaf) - 2.0*sin(thetaf))
            else
                si%control_deflection = 0.0
                si%cf_c = 0.0
                si%ef = 0.0
                si%Cmdelta = 0.0
            end if
        end if
!************* this may need fixing.
        !adjust start point
        if(t%side.eq.'left') then
            start(:) = si%P1(:)
        else
            start(:) = si%P2(:)
        end if
        span = span + sqrt( (si%P2(2)-si%P1(2))**2 + (si%P2(3)-si%P1(3))**2)
        t%area = t%area + si%ds
    end do
    t%tip(:) = start(:)
    write(*,*) 'Wing properties for wing named: ',t%name
    write(*,*) '   Span  : ',span
    write(*,*) '   Area  : ',t%area
    write(*,*) 'Using Sweep Definition: ',t%sweep_definition
    write(*,*)
end subroutine wing_setup


!-----------------------------------------------------------------------------------------------------------
subroutine wing_get_airfoil_weighting(t,percent,weight,af1,af2)
    type(wing_t) :: t
    real, intent(in) :: percent ! percent within total airfoil vector
    real, intent(out) :: weight ! weighting of bounding airfoils
    type(airfoil_t),pointer,intent(out) :: af1,af2 !pointers to bounding airfoils
    integer :: lo,hi
    real :: pct

    pct = min(max(percent, 0.0), REAL(t%nairfoils) - 1.0)
    lo = floor(pct) + 1
    hi = min(lo + 1, t%nairfoils)
    weight = pct - REAL(lo - 1)

    if(hi.eq.lo) weight = 0.0

    af1 => t%airfoils(lo)%p
    af2 => t%airfoils(hi)%p

end subroutine wing_get_airfoil_weighting
!-----------------------------------------------------------------------------------------------------------
subroutine wing_flip_groundplane(t,alpha,CG,hag)
    type(wing_t) :: t
    type(section_t),pointer :: si
    real :: alpha,CG(3),hag
    integer :: isec
    real :: A,B,C,D,P1(3),P2(3),P3(3),tempv(3),tempr,zero
    zero = 0.0
    P1(1) = CG(1) - hag*sin(alpha) !offset from CG
    P1(2) = CG(2)
    P1(3) = CG(3) + hag*cos(alpha)
    P2(:) = P1(:)
    P2(2) = P2(2) + 1.0
    P3(1) = P1(1) - cos(alpha)
    P3(2) = P1(2)
    P3(3) = P1(3) - sin(alpha)

    !equation for the plane
    A = P1(2) * (P2(3) - P3(3)) + P2(2) * (P3(3) - P1(3)) + P3(2) * (P1(3) - P2(3))
    B = P1(3) * (P2(1) - P3(1)) + P2(3) * (P3(1) - P1(1)) + P3(3) * (P1(1) - P2(1))
    C = P1(1) * (P2(2) - P3(2)) + P2(1) * (P3(2) - P1(2)) + P3(1) * (P1(2) - P2(2))
    D = -(P1(1)*(P2(2)*P3(3) - P3(2)*P2(3)) + P2(1)*(P3(2)*P1(3) - P1(2)*P3(3)) + P3(1)*(P1(2)*P2(3) - P2(2)*P1(3)))


    do isec=1,t%nSec
        si => t%sec(isec)

        tempr = si%chord_1
        si%chord_1 = si%chord_2
        si%chord_2 = tempr

        call math_reflect_point(A,B,C,zero,si%un,tempv)
        si%un = tempv
        call math_reflect_point(A,B,C,zero,si%ua,tempv)
        si%ua = tempv
        call math_reflect_point(A,B,C,zero,si%uf,tempv)
        si%uf = tempv

        call math_cross_product(si%ua,si%un,si%us)

        tempv = si%P1
        si%P1 = si%P2
        si%P2 = tempv
        call math_reflect_point(A,B,C,D,si%P1,tempv)
        si%P1 = tempv
        call math_reflect_point(A,B,C,D,si%P2,tempv)
        si%P2 = tempv
        si%dl(:) = si%P2(:) - si%P1(:)
        call math_reflect_point(A,B,C,D,si%PC,tempv)
        si%PC = tempv
        si%zeta(:) = si%chord_c*si%dl(:)/si%ds
        call math_reflect_point(A,B,C,D,t%root,tempv)
        si%Rroot(:) = si%PC(:) - tempv(:)
!        si%af1 => t%af2
!        si%af2 => t%af1
!        temp_percent = si%percent_1
!        si%percent
!        si
    end do

end subroutine wing_flip_groundplane

!-----------------------------------------------------------------------------------------------------------
subroutine wing_write_attributes(t)
    type(wing_t) :: t

    write(*,*) '-------------------------------------------------------------'
    write(*,*) '  Wing Name : ',trim(t%name)
    write(*,*) '       side : ',t%side
    write(*,*) '          x : [',t%root(1),']'
    write(*,*) '          y : [',t%root(2),']'
    write(*,*) '          z : [',t%root(3),']'
    write(*,*) '       span : [',t%span,']'
end subroutine wing_write_attributes


end module wing_m
