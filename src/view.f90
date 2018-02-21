module view_m
#ifdef dnad
    use dnadmod
#define real type(dual)
#endif
    use plane_m
    use wing_m
    use section_m
    use airfoil_m
    implicit none

contains

!-----------------------------------------------------------------------------------------------------------
subroutine view_plotmtv(t)
    type(plane_t) :: t
    type(section_t),pointer :: si
    character(100) :: filename
    character(140) :: command
    integer :: i,iwing,isec,gnum,iforce,ipoint,min_isec
    integer :: ierror = 0
    real :: P0(3),P1(3),P2(3),gpsize,delta,force_pos(3),force_dir(3),force_mag,dist,min_dist

    filename = trim(adjustl(t%master_filename))//'_PlotData.txt'
    open(unit = 10, File = filename, action = 'write', iostat = ierror)

    write(10,*) ' $ DATA=CURVE3D'
    write(10,*) ' % toplabel   = "Planform"'
    write(10,*) ' % xlabel     = "-x"'
    write(10,*) ' % ylabel     = "-y"'
    write(10,*) ' % zlabel     = "-z"'
    write(10,*) ' % grid       = False'
!    write(10,*) ' % axislabel  = False'
    write(10,*) ' % equalscale = True'
    write(10,*) ' % axisscale  = False'
!    write(10,*) ' % axisguides = False'
    write(10,*) ' % eyepos.x   = 0.50'
    write(10,*) ' % eyepos.y   = 0.75'
    write(10,*) ' % eyepos.z   = 0.25'
!    write(10,*) ' % xmin   = -0.60'
!    write(10,*) ' % xmax   =  0.60'
!    write(10,*) ' % ymin   = -0.70'
!    write(10,*) ' % ymax   =  0.2'
    write(10,*) ' % dlinecolor = 6'
    write(10,*) ' '

    do iwing=1,t%nrealwings !real wings
        do isec=1,t%wings(iwing)%nSec
            si => t%wings(iwing)%sec(isec)
            write(10,*) ' % linecolor = 5' !Quarter Chord
            write(10,*) -si%P1(:)
            write(10,*) -si%P2(:)
            write(10,*) ' '

            write(10,*) ' % linecolor = 5' !leading edge
            write(10,*) -si%P1(:)
            write(10,*) -(si%P1(:) - 0.25*si%chord_1*si%ua(:))
            write(10,*) -(si%P2(:) - 0.25*si%chord_2*si%ua(:))
            write(10,*) -si%P2(:)
            write(10,*) ' '

            write(10,*) ' % linecolor = 5' !trailing edge
            write(10,*) -si%P1(:)
            write(10,*) -(si%P1(:) + (0.75-si%cf_c)*si%chord_1*si%ua(:))
            write(10,*) -(si%P2(:) + (0.75-si%cf_c)*si%chord_2*si%ua(:))
            write(10,*) -si%P2(:)
            write(10,*) ' '

            if(si%cf_c > 0.0) then
                write(10,*) ' % linecolor = 7' !flaps
                P1 = -(si%P1(:) + (0.75-si%cf_c)*si%chord_1*si%ua(:))
                P2 = -(si%P2(:) + (0.75-si%cf_c)*si%chord_2*si%ua(:))
                write(10,*) P1
                write(10,*) P1 - si%cf_c*si%chord_1*si%uf(:)
                write(10,*) P2 - si%cf_c*si%chord_2*si%uf(:)
                write(10,*) P2
                write(10,*) ' '
            end if
        end do
    end do

    do iwing=t%nrealwings+1,t%nwings !reflected wings
        do isec=1,t%wings(iwing)%nSec
            si => t%wings(iwing)%sec(isec)
            write(10,*) ' % linecolor = 3' !Quarter Chord
            write(10,*) -si%P1(:)
            write(10,*) -si%P2(:)
            write(10,*) ' '

            write(10,*) ' % linecolor = 3' !leading edge
            write(10,*) -si%P1(:)
            write(10,*) -(si%P1(:) - 0.25*si%chord_1*si%ua(:))
            write(10,*) -(si%P2(:) - 0.25*si%chord_2*si%ua(:))
            write(10,*) -si%P2(:)
            write(10,*) ' '

            write(10,*) ' % linecolor = 3' !trailing edge
            write(10,*) -si%P1(:)
            write(10,*) -(si%P1(:) + (0.75-si%cf_c)*si%chord_1*si%ua(:))
            write(10,*) -(si%P2(:) + (0.75-si%cf_c)*si%chord_2*si%ua(:))
            write(10,*) -si%P2(:)
            write(10,*) ' '

            if(si%cf_c > 0.0) then
                write(10,*) ' % linecolor = 9' !flaps
                P1 = -(si%P1(:) + (0.75-si%cf_c)*si%chord_1*si%ua(:))
                P2 = -(si%P2(:) + (0.75-si%cf_c)*si%chord_2*si%ua(:))
                write(10,*) P1
                write(10,*) P1 - si%cf_c*si%chord_1*si%uf(:)
                write(10,*) P2 - si%cf_c*si%chord_2*si%uf(:)
                write(10,*) P2
                write(10,*) ' '
            end if
        end do
    end do

    if(t%Gammas(1) .ne. 0.0) then
        do i=1,t%nSize ! %normal
            si => t%sec(i)%myp
            write(10,*) ' % linecolor = 6'
            write(10,*) -si%PC(:)
!            write(10,*) -(si%PC(:) + 10.0*sec_CL(si)*si%un(:))
            write(10,*) -(si%PC(:) + 10.0*t%Gammas(i)*si%un(:))
            write(10,*) ' '
        end do
    end if

    if(t%groundplane .eq. 1) then !draw a ground plane
        gpsize = 10.0 !length of ground plane
        gnum = 5 !number of squares on plane
        delta = gpsize/REAL(gnum)

        P0(1) = t%CG(1) - t%hag*sin(t%alpha) !offset from CG !this code copied from wing.f90
        P0(2) = t%CG(2)
        P0(3) = t%CG(3) + t%hag*cos(t%alpha)

        do i=-gnum,gnum
            write(10,*) ' % linecolor = 8'
            P1(1) = P0(1) - gpsize*cos(t%alpha); P1(2) = P0(2) + REAL(i)*delta; P1(3) = P0(3) - gpsize*sin(t%alpha)
            P2(1) = P0(1) + gpsize*cos(t%alpha); P2(2) = P0(2) + REAL(i)*delta; P2(3) = P0(3) + gpsize*sin(t%alpha)
            write(10,*) -P1(:)
            write(10,*) -P2(:)
            write(10,*) ' '
        end do
        do i=-gnum,gnum
            write(10,*) ' % linecolor = 8'
            P1(1) = P0(1) - REAL(i)*delta*cos(t%alpha); P1(2) = P0(2) + gpsize; P1(3) = P0(3) - REAL(i)*delta*sin(t%alpha)
            P2(1) = P0(1) - REAL(i)*delta*cos(t%alpha); P2(2) = P0(2) - gpsize; P2(3) = P0(3) - REAL(i)*delta*sin(t%alpha)
            write(10,*) -P1(:)
            write(10,*) -P2(:)
            write(10,*) ' '
        end do
    end if

    !show external forces tied to closest structure
    do iforce = 1,size(t%external_forces)
        do ipoint=1,t%external_forces(iforce)%datasize
            force_pos(:) = t%external_forces(iforce)%RawData(ipoint,1:3)
            force_dir(:) = t%external_forces(iforce)%RawData(ipoint,4:6)
            force_mag    = t%external_forces(iforce)%RawData(ipoint,7)
            min_dist = 1.0e16
            min_isec = 0
            do isec=1,t%nSize
                si => t%sec(isec)%myp
                dist = math_length(3,force_pos(:),si%PC(:))
                if(dist < min_dist) then
                    min_dist = dist
                    min_isec = isec
                end if
            end do
            si => t%sec(min_isec)%myp
            write(10,*) ' % linecolor = 9'
            write(10,*) -si%PC(:)
            write(10,*) -force_pos(:)
            write(10,*)
            write(10,*) ' % linecolor = 10'
            write(10,*) -force_pos(:)
            write(10,*) -(force_pos(:) + force_dir(:))
            write(10,*)
        end do
    end do

    close(10)
    command = '/sw/bin/Plotmtv -fg white -bg black '//filename
    call system(command)
    command = 'Plotmtv -fg white -bg black '//filename
    call system(command)
end subroutine view_plotmtv

!-----------------------------------------------------------------------------------------------------------

subroutine view_vtk(t)
    type(plane_t) :: t
    type(section_t),pointer :: si
    character(100) :: filename,str1
    integer :: iwing,isec
    integer :: ierror = 0
    120 Format(3ES25.13)
    130 Format(3I5)

    filename = trim(adjustl(t%master_filename))//'_view.vtk'
    open(unit = 10, File = filename, action = 'write', iostat = ierror)

    write(10,'(A)') '# vtk DataFile Version 3.0'
    write(10,'(A)') 'vtk output'
    write(10,'(A)') 'ASCII'
    write(10,'(A)') 'DATASET POLYDATA'
    write(10,'(A)')
    write(10,'(A,I5,A)') 'POINTS ',t%nSize*2,' float'

    do iwing=1,t%nrealwings !real wings
        do isec=1,t%wings(iwing)%nSec
            si => t%wings(iwing)%sec(isec)
            write(str1,120) si%PC(:); write(10,'(A)') trim(adjustl(str1))
            write(str1,120) si%PC(:) + (0.75-si%cf_c)*si%chord_c*si%ua(:); write(10,'(A)') trim(adjustl(str1))
        end do
    end do

    write(10,*)
    write(10,'(A,I5,I5)') 'LINES ',t%nSize,t%nSize*3
    do iwing=0,t%nSize-1
        write(str1,130) 2,2*iwing,2*iwing+1; write(10,'(A)') trim(adjustl(str1))
    end do

    close(10)
end subroutine view_vtk

!-----------------------------------------------------------------------------------------------------------

subroutine view_stl(t)
    type(plane_t) :: t
    type(section_t),pointer :: si
    real,allocatable,dimension(:,:) :: af_points1,af_points2
    character(100) :: filename
    integer :: i,iwing,isec,af_datasize
    integer :: ierror = 0
    real :: percent_1, percent_2, percent_c, temp
    real :: chord_1, chord_2, RA, dtheta

    do i=1,size(airfoils)
        call af_create_geom_from_file(airfoils(i),DB_Airfoil)
    end do

    filename = trim(adjustl(t%master_filename))//'_view.stl'
    open(unit = 10, File = filename, action = 'write', iostat = ierror)
!    filename = trim(adjustl(t%master_filename))//'_view.web'
!    open(unit = 20, File = filename, action = 'write', iostat = ierror)

    write(10,'(A)') 'solid geom'

    do iwing=1,t%nrealwings !real wings

        af_datasize = t%wings(iwing)%airfoils(1)%p%geom%datasize
        do i=1,t%wings(iwing)%nairfoils
            if(t%wings(iwing)%airfoils(i)%p%geom%datasize .ne. af_datasize) then
                write(*,*) 'All airfoils for wing ',t%wings(iwing)%name,' must have same number of nodes.'
                stop
            end if
        end do

        allocate(af_points1(af_datasize,3))
        allocate(af_points2(af_datasize,3))

        if(t%wings(iwing)%is_linear.eq.1) then
            do isec=1,t%wings(iwing)%nSec
                si => t%wings(iwing)%sec(isec)

                dtheta = pi / REAL(t%wings(iwing)%nSec)

                percent_1 = 0.5*(1.0-cos(dtheta*REAL(isec-1)))
                percent_2 = 0.5*(1.0-cos(dtheta*REAL(isec)))
                percent_c = 0.5*(1.0-cos(dtheta*(REAL(isec)-0.5)))
                if(t%wings(iwing)%side.eq.'left') then
                    temp = percent_1
                    percent_1 = percent_2
                    percent_2 = temp
                end if

                if(t%wings(iwing)%chord_2 >= 0.0) then
                    chord_1 = t%wings(iwing)%chord_1 + percent_1*(t%wings(iwing)%chord_2 - t%wings(iwing)%chord_1)
                    chord_2 = t%wings(iwing)%chord_1 + percent_2*(t%wings(iwing)%chord_2 - t%wings(iwing)%chord_1)
                else
                    RA = 8.0 * t%wings(iwing)%span / pi / t%wings(iwing)%chord_1
                    chord_1 = 8.0*t%wings(iwing)%span / pi / RA * sqrt(1.0 - percent_1**2)
                    chord_2 = 8.0*t%wings(iwing)%span / pi / RA * sqrt(1.0 - percent_2**2)
                end if

                call view_create_local_airfoil(si%af1_a,si%af1_b,t%wings(iwing)%side,si%af_weight_1,chord_1,&
                                            & si%twist1,si%dihedral1,af_datasize,si%P1(:),af_points1)
                call view_create_local_airfoil(si%af2_a,si%af2_b,t%wings(iwing)%side,si%af_weight_2,chord_2,&
                                            & si%twist2,si%dihedral2,af_datasize,si%P2(:),af_points2)
                call view_create_stl_shell(af_datasize,af_points1,af_points2)
            end do
        end if

        if(t%wings(iwing)%side == 'right') then
            si => t%wings(iwing)%sec(1)
            call view_create_local_airfoil(si%af1_a,si%af1_b,t%wings(iwing)%side,si%af_weight_1,si%chord_1,&
            & si%twist1,si%dihedral1,af_datasize,si%P1(:),af_points1)
            si => t%wings(iwing)%sec(t%wings(iwing)%nSec)
            call view_create_local_airfoil(si%af2_a,si%af2_b,t%wings(iwing)%side,si%af_weight_2,si%chord_2,&
            & si%twist2,si%dihedral2,af_datasize,si%P2(:),af_points2)
        else
            si => t%wings(iwing)%sec(1)
            call view_create_local_airfoil(si%af2_a,si%af2_b,t%wings(iwing)%side,si%af_weight_2,si%chord_2,&
                                        & si%twist2,si%dihedral2,af_datasize,si%P2(:),af_points2)
            si => t%wings(iwing)%sec(t%wings(iwing)%nSec)
            call view_create_local_airfoil(si%af1_a,si%af1_b,t%wings(iwing)%side,si%af_weight_1,si%chord_1,&
                                        & si%twist1,si%dihedral1,af_datasize,si%P1(:),af_points1)
        end if

        if(t%wings(iwing)%is_linear.ne.1) then
            call view_create_stl_shell(af_datasize,af_points1,af_points2)
        end if
!        call view_create_stl_rib(af_datasize,af_points1)
!        call view_create_stl_rib(af_datasize,af_points2)
        deallocate(af_points1)
        deallocate(af_points2)
    end do

    write(10,*) 'endsolid geom'
    close(10)
!    close(20)
end subroutine view_stl

!-----------------------------------------------------------------------------------------------------------
subroutine view_create_stl_shell(datasize,af_points1,af_points2)
    integer :: datasize,i
    real :: af_points1(datasize,3),af_points2(datasize,3)
    real :: P1(3),P2(3),P3(3)

    !Outer surface
    do i=1,datasize-1
        P1(:) = af_points1(i,:)
        P2(:) = af_points2(i,:)
        P3(:) = af_points1(i+1,:)

        call view_add_stl_triangle(P1,P3,P2)

        P1 = P2
        P2(:) = af_points2(i+1,:)
        call view_add_stl_triangle(P1,P3,P2)
    end do
    !Close Surface
    i = datasize
    P1(:) = af_points1(i,:)
    P2(:) = af_points1(1,:)
    P3(:) = af_points2(i,:)

    call view_add_stl_triangle(P1,P3,P2)

    P1 = P2
    P2(:) = af_points2(1,:)
    call view_add_stl_triangle(P1,P3,P2)
end subroutine view_create_stl_shell

!-----------------------------------------------------------------------------------------------------------
subroutine view_create_stl_rib(datasize,af_points)
    integer :: datasize,i
    real :: af_points(datasize,3)
    real :: P1(3),P2(3),P3(3)

    do i=1,(datasize-2)/2
        P1(:) = af_points(i,:)
        P2(:) = af_points(datasize-i+1,:)
        P3(:) = af_points(i+1,:)
        call view_add_stl_triangle(P1,P3,P2)

        P1 = P2
        P2(:) = af_points(datasize-i,:)
        call view_add_stl_triangle(P1,P3,P2)
    end do
end subroutine view_create_stl_rib

!-----------------------------------------------------------------------------------------------------------
subroutine view_add_stl_triangle(P1,P2,P3)
    real :: P1(3),P2(3),P3(3),norm(3)
    110 Format(A15, 3ES25.13)

    call math_plane_normal(P1,P2,P3,norm)

#ifndef dnad
    write(10,110) 'facet normal ',norm(:)
    write(10,*) 'outer loop'
    write(10,110) 'vertex ',P1(:)
    write(10,110) 'vertex ',P2(:)
    write(10,110) 'vertex ',P3(:)
    write(10,*) 'endloop'
    write(10,*) 'endfacet'
#else
    write(10,110) 'facet normal ',norm(:)%x
    write(10,*) 'outer loop'
    write(10,110) 'vertex ',P1(:)%x
    write(10,110) 'vertex ',P2(:)%x
    write(10,110) 'vertex ',P3(:)%x
    write(10,*) 'endloop'
    write(10,*) 'endfacet'
#endif

!    write(20,120) P1(:),P2(:),P3(:)
end subroutine view_add_stl_triangle

!-----------------------------------------------------------------------------------------------------------
subroutine view_create_local_airfoil(af1,af2,side,percent,chord,twist,dihedral,datasize,point,output)
    type(airfoil_t),pointer :: af1
    type(airfoil_t),pointer :: af2
    character(5) :: side
    real :: percent,chord,twist,dihedral
    integer :: datasize,i
    real :: point(3),output(datasize,3)

    output(:,1:2) = af1%geom%RawData(:,:) + percent*(af2%geom%RawData(:,:) - af1%geom%RawData(:,:))
    output(:,1) = output(:,1) - 0.25
    output = chord*output
    output(:,1) = -output(:,1) !flip x
    output(:,3) = -output(:,2) !assign y to z
    output(:,2) = 0.0 !set y to zero


    if(side.eq.'left') then
!        twist = -twist
        dihedral = -dihedral
    end if

    do i=1,datasize
        call math_rot_y(output(i,:),twist)
        call math_rot_x(output(i,:),-dihedral)
        output(i,:) = output(i,:) + point(:)
    end do

end subroutine view_create_local_airfoil


!-----------------------------------------------------------------------------------------------------------
subroutine view_panair(t)
    type(plane_t) :: t
    real,allocatable,dimension(:,:) :: af_points
    character(100) :: filename
    integer :: i,iwing,af_datasize
    integer :: ierror = 0

    do i=1,size(airfoils)
        call af_create_geom_from_file(airfoils(i),DB_Airfoil)
    end do

    do iwing=1,t%nrealwings !real wings
        af_datasize = t%wings(iwing)%airfoils(1)%p%geom%datasize
        do i=1, t%wings(iwing)%nairfoils
            if(t%wings(iwing)%airfoils(i)%p%geom%datasize .ne. af_datasize) then
                write(*,*) 'All airfoils for wing ',t%wings(iwing)%name,' must have same number of nodes.'
                stop
            end if
        end do

        ! Open the file
        if(t%wings(iwing)%orig_side .eq. "both") then
            write(filename, *) trim(adjustl(t%master_filename))//'_'//trim(adjustl(t%wings(iwing)%name))//"-"//trim(adjustl(t%wings(iwing)%side))//'.panair'
        else
            write(filename, *) trim(adjustl(t%master_filename))//'_'//trim(adjustl(t%wings(iwing)%name))//'.panair'
        end if
        open(unit = 10, File = filename, action = 'write', iostat = ierror)
        allocate(af_points(af_datasize,3))
        call view_write_panair_header(iwing)
        call view_write_panair_upper(t%wings(iwing), af_datasize, af_points)
        call view_write_panair_lower(t%wings(iwing), af_datasize, af_points)
        call view_write_panair_wake()
        write(10, "(A)") "$END"
        deallocate(af_points)
        close(10)
    end do

end subroutine view_panair


subroutine view_write_panair_header(iwing)
    integer, intent(in) :: iwing

    ! Write header info
    write(10, "(A, I0)") "$POINTS for wing ", iwing
    write(10, "(A)") "=kn"  ! Number of networks in $POINTS block
    write(10, "(I0)") 2
    write(10, "(A)") "=kt"  ! Boundary condition (1 = solid surface)
    write(10, "(I0)") 1

end subroutine view_write_panair_header


subroutine view_write_panair_upper(wi, af_datasize, af_points)
    type(wing_t), intent(in) :: wi
    integer, intent(in) :: af_datasize
    real, allocatable, dimension(:,:), intent(inout) :: af_points

    call view_write_panair_network(wi, af_datasize, af_points, "upper", af_datasize, af_datasize / 2 + 1)
end subroutine view_write_panair_upper


subroutine view_write_panair_lower(wi, af_datasize, af_points)
    type(wing_t), intent(in) :: wi
    integer, intent(in) :: af_datasize
    real, allocatable, dimension(:,:), intent(inout) :: af_points

    call view_write_panair_network(wi, af_datasize, af_points, "lower", af_datasize / 2 + 1, 1)
end subroutine view_write_panair_lower


subroutine view_write_panair_network(wi, af_datasize, af_points, network, istart, iend)
    type(wing_t), intent(in) :: wi
    integer, intent(in) :: af_datasize
    real, allocatable, dimension(:,:), intent(inout) :: af_points
    character(len=*), intent(in) :: network
    integer, intent(in) :: istart, iend

    integer :: isec, isec_start, isec_end, isec_inc, isec_side1, isec_side2
    type(section_t), pointer :: si

    write(10, "(A, T11, A)") "=nm", "nn"
    write(10, "(I0, T11, I0, T71, A)") istart - iend + 1, wi%nSec + 1, network

    if(wi%side .eq. "left") then
        isec_start = wi%nSec
        isec_end = 1
        isec_inc = -1
        isec_side1 = 1
        isec_side2 = 2
    else
        isec_start = 1
        isec_end = wi%nSec
        isec_inc = 1
        isec_side1 = 1
        isec_side2 = 2
    end if

    do isec = isec_start, isec_end, isec_inc
        si => wi%sec(isec)
        call view_create_local_airfoil_panair(wi, si, isec_side1, af_datasize, af_points)
        call view_write_panair_points(istart, iend, af_points)
    end do

    call view_create_local_airfoil_panair(wi, si, isec_side2, af_datasize, af_points)
    call view_write_panair_points(istart, iend, af_points)

end subroutine view_write_panair_network


subroutine view_create_local_airfoil_panair(wi, si, sec_side, af_datasize, af_points)
    type(wing_t), intent(in) :: wi
    type(section_t), pointer, intent(in) :: si
    integer, intent(in) :: sec_side, af_datasize

    real, allocatable, dimension(:,:), intent(inout) :: af_points

    real :: percent, chord, RA
    integer :: i

    if(sec_side .eq. 1) then
        percent = si%percent_1
    else if(sec_side .eq. 2) then
        percent = si%percent_2
    else
        percent = si%percent_c
    end if

    if(wi%chord_2 >= 0.0) then
        chord = wi%chord_1 + percent*(wi%chord_2 - wi%chord_1)
    else
        RA = 8.0 * wi%span / pi / wi%chord_1
        chord = 8.0 * wi%span / pi / RA * sqrt(1.0 - percent**2)
    end if

    if(sec_side .eq. 1) then
        call view_create_local_airfoil(si%af1_a, si%af1_b, wi%side, si%af_weight_1, chord, &
                & si%twist1, si%dihedral1, af_datasize, si%P1, af_points)
    else if(sec_side .eq. 2) then
        call view_create_local_airfoil(si%af2_a, si%af2_b, wi%side, si%af_weight_2, chord, &
                & si%twist2, si%dihedral2, af_datasize, si%P2, af_points)
    else
        call view_create_local_airfoil(si%afc_a, si%afc_b, wi%side, si%af_weight_c, chord, &
                & si%twist, si%dihedral, af_datasize, si%PC, af_points)
    end if

    do i = 1, af_datasize
        af_points(i, 1) = -af_points(i, 1)
        af_points(i, 3) = -af_points(i, 3)
    end do

end subroutine view_create_local_airfoil_panair


subroutine view_write_panair_points(istart, iend, af_points)
    integer, intent(in) :: istart, iend
    real, dimension(:, :), intent(in) :: af_points

    integer :: ipt

    do ipt = istart, iend + 1, -2
        write(10, "(6F10.6)") af_points(ipt, 1:3), af_points(ipt - 1, 1:3)
    end do
    if (ipt == iend) then
        write(10, "(3F10.6)") af_points(ipt, 1:3)
    end if
end subroutine view_write_panair_points


subroutine view_write_panair_wake()
    write(10, "(A)") "$TRAILING matchw=0"
    write(10, "(A)") "=kn"
    write(10, "(I0)") 1
    write(10, "(A, T11, A)") "=kt", "matchw"
    write(10, "(I0, T11, I0)") 18, 0
    write(10, "(A, T11, A, T21, A, T31, A)") "=inat", "insd", "xwake", "twake"
    write(10, "(A, T11, I0, T21, I0, T31, I0, T71, A)") "upper", 1, 10, 0, "wake"
end subroutine


end module view_m
