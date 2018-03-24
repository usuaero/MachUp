module dataset_m
#ifdef dnad
    use dnadmod
#define real type(dual)
#endif

    use math_m
    use myjson_m
    implicit none
    type dataset_t
        character(100) :: filename
        integer :: dim !dimensions of the dataset (3 for xyz points)
        integer :: datasize !length of dataset

        real,allocatable,dimension(:,:) :: RawData,Deriv
        real,allocatable,dimension(:) :: Xi ! Xi = array the length of datasize that contains the xi values for each point in the raw data set
        real,allocatable,dimension(:) :: xdata !independent variable array used for spline calcs
    end type

    interface myjson_get
        module procedure :: myjson_file_get_dataset
    end interface

    interface json_value_add
        module procedure :: myjson_value_add_dataset
    end interface


contains

!-----------------------------------------------------------------------------------------------------------
subroutine ds_allocate(t)
    type(dataset_t) :: t

    if(allocated(t%RawData)) call ds_deallocate(t)

    allocate(t%RawData(t%datasize,t%dim))
    allocate(t%Deriv(t%datasize,t%dim))
    allocate(t%Xi(t%datasize))
    allocate(t%xdata(t%datasize))
end subroutine ds_allocate

subroutine ds_deallocate(t)
    type(dataset_t) :: t

    deallocate(t%RawData)
    deallocate(t%Deriv)
    deallocate(t%Xi)
    deallocate(t%xdata)
end subroutine ds_deallocate

!-----------------------------------------------------------------------------------------------------------
subroutine ds_create_from_file(t,filename,dim)
    type(dataset_t) :: t
    character(100) :: filename
    integer :: dim

    if(index(filename, '.json') .ne. 0) then
        !Read as JSON file
        call ds_read_json(t,filename)
        call ds_calc_Xi(t)
    else
        !Read as text file
        call ds_read_file(t,filename,dim)
        call ds_calc_Xi(t)
    end if
end subroutine ds_create_from_file

!-----------------------------------------------------------------------------------------------------------
subroutine ds_create_from_data(t,datasize,dim,rawdata)
    type(dataset_t) :: t
    integer :: datasize,dim
    real :: rawdata(datasize,dim)

    t%datasize = datasize
    t%dim = dim

    call ds_allocate(t)

    t%RawData = rawdata
    call ds_calc_Xi(t)

end subroutine ds_create_from_data

!-----------------------------------------------------------------------------------------------------------
subroutine ds_read_file(t,filename,dim)
    type(dataset_t) :: t
    character(100) :: filename
    character(1000) :: buffer
    integer :: dim,ios,irow

    write(*,*) 'Reading file: ',trim(filename)
    t%dim = dim
    t%datasize = 0
    t%filename = trim(filename)
    open(unit = 100, File = trim(filename), action = "read", iostat = ios)

    ! Find file length
    read(100,*) !header line
    do while (ios == 0)
        read(100,'(a)',iostat=ios) buffer
        if (ios == 0) then
            if(buffer .ne. '') t%datasize = t%datasize + 1
        end if
    end do
    close(100)

    write(*,*) 'data size = ',t%datasize
    call ds_allocate(t)

    open(unit = 100, File = filename, action = "read", iostat = ios)
    read(100,*) !header line
    do irow=1,t%datasize
#ifndef dnad
        read(100,*,iostat=ios) t%RawData(irow,:)
#else
        read(100,*,iostat=ios) t%RawData(irow,:)%x
        t%RawData(irow,:) = t%RawData(irow,:)%x
#endif
    end do
    close(100)

end subroutine ds_read_file

!-----------------------------------------------------------------------------------------------------------
subroutine ds_read_json(t,filename)
    type(dataset_t), intent(out) :: t
    character(len=*), intent(in) :: filename

    type(json_file) :: json    !the JSON structure read from the file

    write(*,*) 'Reading JSON file: ',trim(filename)
    call json%load_file(filename = filename); call json_check()
    call myjson_get(json, '', t)

    !Clean up the JSON file
    call json%destroy()

    write(*,*) 'data size = ',t%datasize
end subroutine ds_read_json

!-----------------------------------------------------------------------------------------------------------
subroutine myjson_file_get_dataset(json, name, value)
    type(json_file) :: json
    character(len=*), intent(in) :: name
    type(dataset_t), intent(out) :: value

    type(json_value), pointer :: json_table, json_row
    integer :: irow, icol, ncol
    character(len=10) :: row_name, col_name
    character(len=10) :: fmt_string

    ! Get the table object from the JSON file
    call json%get(name, json_table)
    if(json_failed()) then
        write(*,*) 'Error: Unable to read required table: ',name
        STOP
    end if

    ! Get the number of rows in the data table
    value%datasize = json_value_count(json_table)

    ! Get the number of columns (dimensions) in the data table
    call json%get('r1', json_row)
    value%dim = json_value_count(json_row)

    ! Allocate dataset memory
    call ds_allocate(value)

    ! Read the data table and assign to dataset cells
    do irow = 1, value%datasize
        write(fmt_string, '(A5,I1,A1)') '(A1,I', int(log10(REAL(irow))) + 1, ')'
        write(row_name, fmt_string) "r", irow
        call json_value_get(json_table, row_name, json_row)
        if(json_failed()) then
            write(*,*) 'Error: Unable to read row ', irow, 'of table ', name
            STOP
        end if

        !Make sure every row has the same number of columns
        ncol = json_value_count(json_row)
        if(ncol .ne. value%dim) then
            write(*,*) 'Error: Row ', irow, ' has ', ncol, ' columns (expected ', value%dim, ')'
            STOP
        end if
        do icol = 1, value%dim
            write(fmt_string, '(A5,I1,A1)') '(A1,I', int(log10(REAL(icol))) + 1, ')'
            write(col_name, fmt_string) "c", icol
            call myjson_get(json_row, col_name, value%RawData(irow, icol))
        end do
    end do
end subroutine myjson_file_get_dataset

!-----------------------------------------------------------------------------------------------------------
subroutine myjson_value_add_dataset(me, name, val)

    implicit none

    type(json_value), pointer   :: me
    character(len=*),intent(in) :: name
    type(dataset_t),intent(in)  :: val

    type(json_value),pointer :: var, row
    integer :: irow, icol
    character(len=2) :: row_name, col_name
    character(len=7) :: fmt_string

    !create the variable as an array:
    irow = len(name)
    if(len(name) .eq. 0) then
        var => me
    else
        call json_value_create(var)
        call to_object(var, name)
    end if

    !populate the table:
    do irow=1, val%datasize
        !format the row name
        write(fmt_string, '(A5,I1,A1)') '(A1,I', int(log10(REAL(irow))) + 1, ')'
        write(row_name, fmt_string) 'r', irow

        !create the row
        call json_value_create(row)
        call to_object(row, row_name)
        do icol=1, val%dim
            !format the column name
            write(fmt_string, '(A5,I1,A1)') '(A1,I', int(log10(REAL(irow))) + 1, ')'
            write(col_name, fmt_string) 'c', icol

            !add the column name/value to the row
            call json_value_add(row, col_name, val%RawData(irow, icol))
        end do

        !add the row to the dataset
        call json_value_add(var, row)
    end do

    !add the dataset
    if(len(name) .ne. 0) then
        call json_value_add(me, var)
        nullify(var)
    end if
end subroutine myjson_value_add_dataset

!-----------------------------------------------------------------------------------------------------------
subroutine ds_calc_Xi(t)
    type(dataset_t) :: t
    integer :: i

    !Calculate Xi based on distance between points
    t%Xi(1) = 0
    do i = 2, t%datasize, 1
        t%Xi(i) = t%Xi(i-1) + math_length(t%dim,t%RawData(i,:),t%RawData(i-1,:))
    end do

    !default xdata to Xi
    call ds_set_xcol(t,0)
end subroutine ds_calc_Xi

!-----------------------------------------------------------------------------------------------------------
subroutine ds_set_xcol(t,xcol)
    type(dataset_t) :: t
    integer :: xcol !column to use as independent variable (0 will use Xi value which is distance to each point in space)

    if(xcol > 0) then
        t%xdata(:) = t%RawData(:,xcol)
    else
        t%xdata(:) = t%Xi(:)
    end if

end subroutine ds_set_xcol
!-----------------------------------------------------------------------------------------------------------
subroutine ds_cubic_setup(t,xcol,bc1,bcval1,bc2,bcval2)
    type(dataset_t) :: t
    integer :: xcol !column to use as independent variable (0 will use Xi value which is distance to each point in space)
    integer :: bc1,bc2 !boundary condition flags at each end (1 = 1st deriv, 2 = 2nd deriv)
    real :: bcval1,bcval2 !boundary condition values at each end

    !Declare Local Variables
    integer :: i,j,n
    real :: dt1,dt2,Amat(t%datasize,t%datasize), Bvec(t%datasize), x(t%datasize), y(t%datasize,t%dim)

    call ds_set_xcol(t,xcol)

    n = t%datasize
    x = t%xdata
    y = t%RawData

    do j=1,t%dim
        Amat = 0.0; Bvec = 0.0 !Amat is actually the same for each j, but I'm too lazy to change it right now
        select case (bc1)
            case (1) !first derivative specified
                Amat(1,1) = 1.0; Bvec(1) = bcval1
            case (2) !second derivative specified
                Amat(1,1) = 2.0; Amat(1,2) = 1.0; Bvec(1) = 3.0*(y(2,j) - y(1,j))/(x(2) - x(1)) - 0.5*(x(2) - x(1))*bcval1
            case (3) !first derivative equal at endpoints
                Amat(1,1) = 1.0; Amat(1,n) = -1.0; Bvec(1) = 0.0
            case (4) !second derivative equal at endpoints
                dt1 = x(2) - x(1)
                dt2 = x(n) - x(n-1)
                Amat(1,1) = 2.0*dt2; Amat(1,2) = dt2; Amat(1,n-1) = dt1; Amat(1,n) = 2.0*dt1;
                Bvec(1) = 3.0*dt2/dt1*(y(2,j) - y(1,j)) + 3.0*dt1/dt2*(y(n,j) - y(n-1,j))
        end select
        do i=2,n-1
            dt1 = (x(i+1) - x(i))/(x(i+1) - x(i-1))
            dt2 = (x(i) - x(i-1))/(x(i+1) - x(i-1))
            Amat(i,i-1) = dt1; Amat(i,i) = 2.0; Amat(i,i+1) = dt2
            Bvec(i) = 3.0*dt1*(y(i,j) - y(i-1,j))/(x(i) - x(i-1)) + 3.0*dt2*(y(i+1,j) - y(i,j))/(x(i+1) - x(i))
        end do

        select case (bc2)
            case (1) !first derivative specified
                Amat(n,n) = 1.0; Bvec(n) = bcval2
            case (2) !second derivative specified
                Amat(n,n-1) = 1.0; Amat(n,n) = 2.0; Bvec(n) = 3.0*(y(n,j) - y(n-1,j))/(x(n) - x(n-1)) + 0.5*(x(n) - x(n-1))*bcval2
            case (3) !first derivative equal at endpoints
                Amat(n,1) = 1.0; Amat(n,n) = -1.0; Bvec(n) = 0.0
            case (4) !second derivative equal at endpoints
                dt1 = x(2) - x(1)
                dt2 = x(n) - x(n-1)
                Amat(n,1) = 2.0*dt2; Amat(n,2) = dt2; Amat(n,n-1) = dt1; Amat(n,n) = 2.0*dt1;
                Bvec(n) = 3.0*dt2/dt1*(y(2,j) - y(1,j)) + 3.0*dt1/dt2*(y(n,j) - y(n-1,j))
        end select

        call math_snyder_ludcmp(Amat,n)
        call math_snyder_lusolv(Amat,Bvec,t%Deriv(:,j),n)
    end do

end subroutine ds_cubic_setup

!-----------------------------------------------------------------------------------------------------------
subroutine ds_weighted_interpolate(t,value,weight,ans)
    type(dataset_t) :: t
    real :: value !value of independent variable
    real :: weight ! = 0 for linear, = 1 for cubic spline
    real :: ans(t%dim) !returns the answer through this array the size of t%dim
    real, allocatable,dimension(:) :: ans_linear, ans_cubic

    allocate(ans_linear(t%dim))
    allocate(ans_cubic(t%dim))

    call ds_linear_interpolate(t,value,ans_linear)
    call ds_cubic_interpolate(t,value,0,ans_cubic)

    ans(:) = ans_linear(:) + weight*(ans_cubic(:) - ans_linear(:))

    deallocate(ans_linear)
    deallocate(ans_cubic)

end subroutine ds_weighted_interpolate

!-----------------------------------------------------------------------------------------------------------
subroutine ds_linear_interpolate(t,value,ans)
    type(dataset_t) :: t
    real :: value !value of independent variable
    real :: ans(t%dim) !returns the answer through this array the size of t%dim

    integer :: n,j,i,ival
    real :: zeta, dt
    real :: x(t%datasize), y(t%datasize,t%dim)

    n = t%datasize
    x = t%xdata
    y = t%RawData

    !find correct interval
    ival = 0
    if(value < x(2)) then
        ival = 1
    elseif(value > x(n-1)) then
        ival = n-1
    else
        do i=1,n-1
            if(value > x(i)) then
                ival = i
            end if
        end do
    end if
    i = ival

    dt = x(i+1) - x(i)
    zeta = (value - x(i))/dt

    do j=1,t%dim
        ans(j) = y(i,j) + zeta*(y(i+1,j)-y(i,j))
    end do

end subroutine ds_linear_interpolate

!-----------------------------------------------------------------------------------------------------------
subroutine ds_cubic_interpolate(t,value,flag,ans)
    type(dataset_t) :: t
    real :: value !value of independent variable
    integer :: flag !0 = return spline value, 1 = return 1st deriv, 2 = return 2nd deriv
    real :: ans(t%dim) !returns the answer through this array the size of t%dim

    integer :: n,j,i,ival
    real :: zeta, dt
    real, allocatable, dimension(:) :: x
    real, allocatable, dimension(:,:) :: y

    allocate(x(t%datasize))
    allocate(y(t%datasize,t%dim))

    n = t%datasize
    x = t%xdata
    y = t%RawData

    !find correct interval
    ival = 0
    if(value < x(2)) then
        ival = 1
    elseif(value > x(n-1)) then
        ival = n-1
    else
        do i=1,n-1
            if(value > x(i)) then
                ival = i
            end if
        end do
    end if
    i = ival

    dt = x(i+1) - x(i)
    zeta = (value - x(i))/dt

    select case (flag)
        case (0) !return spline value
            do j=1,t%dim
                ans(j) = y(i,j)*(1.0-3.0*zeta**2 + 2.0*zeta**3)   + y(i+1,j)*(3.0*zeta**2 - 2.0*zeta**3) + &
                       & t%Deriv(i,j)*dt*(zeta - 2.0*zeta**2 + zeta**3) + t%Deriv(i+1,j)*dt*(zeta**3 - zeta**2)
            end do
        case (1) !return 1st Deriv
            do j=1,t%dim
                ans(j) = y(i,j)/dt*(6.0*zeta**2 - 6.0*zeta)  + y(i+1,j)/dt*(6.0*zeta - 6.0*zeta**2) + &
                       & t%Deriv(i,j)*(1.0 - 4.0*zeta + 3.0*zeta**2) + t%Deriv(i+1,j)*(3.0*zeta**2 - 2.0*zeta)
            end do
        case (2) !return 2nd Deriv
            do j=1,t%dim
                ans(j) = y(i,j)/dt**2*(12.0*zeta - 6.0)  + y(i+1,j)/dt**2*(6.0 - 12.0*zeta) + &
                       & t%Deriv(i,j)/dt*(6.0*zeta - 4.0) + t%Deriv(i+1,j)/dt*(6.0*zeta - 2.0)
            end do
    end select

    deallocate(x)
    deallocate(y)

end subroutine ds_cubic_interpolate

!-----------------------------------------------------------------------------------------------------------
subroutine ds_print_data(t)
    type(dataset_t) :: t
    integer :: i

    write(*,*) '          i             Xi(i)               xdata(i)              RawData(i) ...'
    do i=1,t%datasize
        write(*,*) i,t%Xi(i), t%xdata(i), t%RawData(i,:)
    end do

end subroutine ds_print_data

!-----------------------------------------------------------------------------------------------------------
real function ds_linear_interpolate_col(t,value,value_col,ans_col)
    type(dataset_t) :: t
    real :: value, percent
    integer :: value_col,ans_col,i

    real,allocatable,dimension(:) :: col1,col2

    allocate(col1(t%datasize))
    allocate(col2(t%datasize))

    if(value_col .eq. 0) then
        col1 = t%Xi/t%Xi(t%datasize)
    else
        col1 = t%RawData(:,value_col)
    end if

    if(ans_col .eq. 0) then
        col2 = t%Xi/t%Xi(t%datasize)
    else
        col2 = t%RawData(:,ans_col)
    end if

    if(value <= col1(1)) then
        ds_linear_interpolate_col = col2(1)
    else if(value >= col1(t%datasize)) then
        ds_linear_interpolate_col = col2(t%datasize)
    else
        i = 1
        do while(value > col1(i))
            i = i + 1
        end do
        percent = (value-col1(i-1))/(col1(i) - col1(i-1))
        ds_linear_interpolate_col = col2(i-1) + percent*(col2(i) - col2(i-1))
    end if

    deallocate(col1)
    deallocate(col2)

end function ds_linear_interpolate_col


end module dataset_m
