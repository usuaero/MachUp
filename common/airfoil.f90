module airfoil_m
#ifdef dnad
    use dnadmod
#define real type(dual)
#endif

    use dataset_m
    implicit none

    type airfoil_t
        character(100) :: name
        character(20) :: properties_type
        integer :: has_data_file
        integer :: has_geom_file

        real :: aL0
        real :: CLa
        real :: CmL0
        real :: Cma
        real :: CD0,CD0L,CD0L2
        real :: CLmax

        type(dataset_t) :: data2D
        type(dataset_t) :: geom
    end type airfoil_t

    type airfoil_p
        type(airfoil_t),pointer :: p
    end type airfoil_p

contains

!-----------------------------------------------------------------------------------------------------------
subroutine af_create_from_data_file(t,filename)
    type(airfoil_t) :: t
    character(100) :: filename
    t%has_data_file = 1
    call ds_create_from_file(t%data2D,filename,4)
    t%CLmax = math_max(t%data2D%datasize, t%data2D%RawData(:,2))
end subroutine af_create_from_data_file

!-----------------------------------------------------------------------------------------------------------
subroutine af_create_geom_from_file(t,DB_Airfoil)
    type(airfoil_t) :: t
    character(200) :: DB_Airfoil
    character(100) :: filename

    if(.not. allocated(t%geom%RawData)) then
        filename = trim(adjustl(DB_Airfoil))//'/'//trim(adjustl(t%name))//'_profile.txt'
        call ds_create_from_file(t%geom,filename,2)
    end if

end subroutine af_create_geom_from_file

!-----------------------------------------------------------------------------------------------------------
real function af_CLa(t,alpha)
    type(airfoil_t),pointer :: t
    real :: alpha,CL1,CL2
    if(t%has_data_file==1) then
        CL1 = ds_linear_interpolate_col(t%data2D,180.0/pi*(alpha-0.005),1,2)
        CL2 = ds_linear_interpolate_col(t%data2D,180.0/pi*(alpha+0.005),1,2)
        af_CLa = (CL2-CL1)/0.01
    else
        af_CLa = t%CLa
    end if
!    af_CLa = 6.9034609677143 !4412 linear
!    af_CLa = 6.9152325411828 - 2.0*0.3793150634847*alpha !4412 nonlinear
!    write(*,*) alpha*180.0/pi,af_CLa
end function af_CLa

!-----------------------------------------------------------------------------------------------------------
real function af_CL(t,alpha)
    type(airfoil_t),pointer :: t
    real :: alpha
    if(t%has_data_file==1) then
        af_CL = ds_linear_interpolate_col(t%data2D,180.0/pi*alpha,1,2)
    else
        af_CL = t%CLa*(alpha-t%aL0)
    end if
!    af_CL = 0.5173557093757 + 6.9034609677143*alpha !4412 linear
!    af_CL = 0.5183846299204 +  6.9152325411828*alpha - 0.3793150634847*(alpha**2) !4412 nonlinear
!    write(*,*) alpha*180.0/pi,af_CL
end function af_CL

!-----------------------------------------------------------------------------------------------------------
real function af_CD(t,alpha)
    type(airfoil_t),pointer :: t
    real :: alpha,CL
    if(t%has_data_file==1) then
        af_CD = ds_linear_interpolate_col(t%data2D,180.0/pi*alpha,1,3)
    else
        CL = af_CL(t,alpha)
        af_CD = t%CD0 + t%CD0L*CL + t%CD0L2*CL**2
    end if
end function af_CD

!-----------------------------------------------------------------------------------------------------------
real function af_Cm(t,alpha)
    type(airfoil_t),pointer :: t
    real :: alpha
    if(t%has_data_file==1) then
        af_Cm = ds_linear_interpolate_col(t%data2D,180.0/pi*alpha,1,4)
    else
        af_Cm = (t%CmL0 + t%Cma*(alpha-t%aL0))
    end if
end function af_Cm

!-----------------------------------------------------------------------------------------------------------
end module airfoil_m

