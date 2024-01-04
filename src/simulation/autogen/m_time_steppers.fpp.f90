# 1 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
!>
!! @file m_time_steppers.f90
!! @brief Contains module m_time_steppers

# 1 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/common/include/macros.fpp" 1
# 11 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/common/include/macros.fpp"

# 21 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/common/include/macros.fpp"

# 27 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/common/include/macros.fpp"

#define t_vec3   real(kind(0d0)), dimension(1:3)
#define t_mat4x4 real(kind(0d0)), dimension(1:4,1:4)
# 6 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp" 2

!> @brief The following module features a variety of time-stepping schemes.
!!              Currently, it includes the following Runge-Kutta (RK) algorithms:
!!                   1) 1st Order TVD RK
!!                   2) 2nd Order TVD RK
!!                   3) 3rd Order TVD RK
!!              where TVD designates a total-variation-diminishing time-stepper.
module m_time_steppers

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_rhs                  !< Right-hand-side (RHS) evaluation procedures

    use m_data_output          !< Run-time info & solution data output procedures

    use m_bubbles              !< Bubble dynamics routines

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_fftw

    use m_nvtx
    ! ==========================================================================

    implicit none

    type(vector_field), allocatable, dimension(:) :: q_cons_ts !<
    !! Cell-average conservative variables at each time-stage (TS)

    type(scalar_field), allocatable, dimension(:) :: q_prim_vf !<
    !! Cell-average primitive variables at the current time-stage

    type(scalar_field), allocatable, dimension(:) :: rhs_vf !<
    !! Cell-average RHS variables at the current time-stage

    type(vector_field), allocatable, dimension(:) :: q_prim_ts !<
    !! Cell-average primitive variables at consecutive TIMESTEPS

    real(kind(0d0)), allocatable, dimension(:, :, :, :, :) :: rhs_pb

    real(kind(0d0)), allocatable, dimension(:, :, :, :, :) :: rhs_mv

    integer, private :: num_ts !<
    !! Number of time stages in the time-stepping scheme

!$acc declare create(q_cons_ts,q_prim_vf,rhs_vf,q_prim_ts, rhs_mv, rhs_pb)

contains

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_time_steppers_module() ! -----------------------

        type(int_bounds_info) :: ix_t, iy_t, iz_t !<
            !! Indical bounds in the x-, y- and z-directions

        integer :: i, j !< Generic loop iterators

        ! Setting number of time-stages for selected time-stepping scheme
        if (time_stepper == 1) then
            num_ts = 1
        elseif (any(time_stepper == (/2, 3/))) then
            num_ts = 2
        end if

        ! Setting the indical bounds in the x-, y- and z-directions
        ix_t%beg = -buff_size; ix_t%end = m + buff_size

        if (n > 0) then
            iy_t%beg = -buff_size; iy_t%end = n + buff_size

            if (p > 0) then
                iz_t%beg = -buff_size; iz_t%end = p + buff_size
            else
                iz_t%beg = 0; iz_t%end = 0
            end if
        else
            iy_t%beg = 0; iy_t%end = 0
            iz_t%beg = 0; iz_t%end = 0
        end if

        ! Allocating the cell-average conservative variables
#ifdef MFC_DEBUG
# 92 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 92 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 92 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 92 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:92: ', '@:ALLOCATE(q_cons_ts(1:num_ts))'
# 92 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 92 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 92 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 92 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(q_cons_ts(1:num_ts))
# 92 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(q_cons_ts(1:num_ts))

        do i = 1, num_ts
#ifdef MFC_DEBUG
# 95 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 95 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 95 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 95 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:95: ', '@:ALLOCATE(q_cons_ts(i)%vf(1:sys_size))'
# 95 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 95 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 95 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 95 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(q_cons_ts(i)%vf(1:sys_size))
# 95 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(q_cons_ts(i)%vf(1:sys_size))
        end do

        do i = 1, num_ts
            do j = 1, sys_size
#ifdef MFC_DEBUG
# 100 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 100 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 100 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 100 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:100: ', '@:ALLOCATE(q_cons_ts(i)%vf(j)%sf(ix_t%beg:ix_t%end, iy_t%beg:iy_t%end, iz_t%beg:iz_t%end))'
# 100 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 100 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 100 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 100 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(q_cons_ts(i)%vf(j)%sf(ix_t%beg:ix_t%end,                     iy_t%beg:iy_t%end,                     iz_t%beg:iz_t%end))
# 100 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(q_cons_ts(i)%vf(j)%sf(ix_t%beg:ix_t%end,                     iy_t%beg:iy_t%end,                     iz_t%beg:iz_t%end))
# 103 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
            end do
        end do

        ! Allocating the cell-average primitive ts variables
        if (probe_wrt) then
#ifdef MFC_DEBUG
# 108 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 108 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 108 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 108 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:108: ', '@:ALLOCATE(q_prim_ts(0:3))'
# 108 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 108 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 108 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 108 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(q_prim_ts(0:3))
# 108 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(q_prim_ts(0:3))

            do i = 0, 3
#ifdef MFC_DEBUG
# 111 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 111 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 111 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 111 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:111: ', '@:ALLOCATE(q_prim_ts(i)%vf(1:sys_size))'
# 111 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 111 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 111 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 111 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(q_prim_ts(i)%vf(1:sys_size))
# 111 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(q_prim_ts(i)%vf(1:sys_size))
            end do

            do i = 0, 3
                do j = 1, sys_size
#ifdef MFC_DEBUG
# 116 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 116 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 116 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 116 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:116: ', '@:ALLOCATE(q_prim_ts(i)%vf(j)%sf(ix_t%beg:ix_t%end, iy_t%beg:iy_t%end, iz_t%beg:iz_t%end))'
# 116 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 116 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 116 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 116 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(q_prim_ts(i)%vf(j)%sf(ix_t%beg:ix_t%end,                         iy_t%beg:iy_t%end,                         iz_t%beg:iz_t%end))
# 116 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(q_prim_ts(i)%vf(j)%sf(ix_t%beg:ix_t%end,                         iy_t%beg:iy_t%end,                         iz_t%beg:iz_t%end))
# 119 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
                end do
            end do
        end if

        ! Allocating the cell-average primitive variables
#ifdef MFC_DEBUG
# 124 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 124 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 124 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 124 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:124: ', '@:ALLOCATE(q_prim_vf(1:sys_size))'
# 124 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 124 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 124 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 124 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(q_prim_vf(1:sys_size))
# 124 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(q_prim_vf(1:sys_size))

        do i = 1, adv_idx%end
#ifdef MFC_DEBUG
# 127 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 127 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 127 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 127 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:127: ', '@:ALLOCATE(q_prim_vf(i)%sf(ix_t%beg:ix_t%end, iy_t%beg:iy_t%end, iz_t%beg:iz_t%end))'
# 127 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 127 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 127 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 127 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(q_prim_vf(i)%sf(ix_t%beg:ix_t%end,                 iy_t%beg:iy_t%end,                 iz_t%beg:iz_t%end))
# 127 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(q_prim_vf(i)%sf(ix_t%beg:ix_t%end,                 iy_t%beg:iy_t%end,                 iz_t%beg:iz_t%end))
# 130 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        end do

        if (bubbles) then
            do i = bub_idx%beg, bub_idx%end
#ifdef MFC_DEBUG
# 134 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 134 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 134 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 134 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:134: ', '@:ALLOCATE(q_prim_vf(i)%sf(ix_t%beg:ix_t%end, iy_t%beg:iy_t%end, iz_t%beg:iz_t%end))'
# 134 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 134 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 134 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 134 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(q_prim_vf(i)%sf(ix_t%beg:ix_t%end,                     iy_t%beg:iy_t%end,                     iz_t%beg:iz_t%end))
# 134 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(q_prim_vf(i)%sf(ix_t%beg:ix_t%end,                     iy_t%beg:iy_t%end,                     iz_t%beg:iz_t%end))
# 137 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
            end do
        end if

#ifdef MFC_DEBUG
# 140 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 140 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 140 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 140 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:140: ', '@:ALLOCATE(pb_ts(1:2))'
# 140 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 140 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 140 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 140 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(pb_ts(1:2))
# 140 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(pb_ts(1:2))
        !Initialize bubble variables pb and mv at all quadrature nodes for all R0 bins
        if (qbmm .and. (.not. polytropic)) then
#ifdef MFC_DEBUG
# 143 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 143 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 143 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 143 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:143: ', '@:ALLOCATE(pb_ts(1)%sf(ix_t%beg:ix_t%end, iy_t%beg:iy_t%end, iz_t%beg:iz_t%end, 1:nnode, 1:nb))'
# 143 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 143 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 143 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 143 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(pb_ts(1)%sf(ix_t%beg:ix_t%end,                 iy_t%beg:iy_t%end,                 iz_t%beg:iz_t%end, 1:nnode, 1:nb))
# 143 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(pb_ts(1)%sf(ix_t%beg:ix_t%end,                 iy_t%beg:iy_t%end,                 iz_t%beg:iz_t%end, 1:nnode, 1:nb))
# 146 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#ifdef MFC_DEBUG
# 146 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 146 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 146 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 146 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:146: ', '@:ALLOCATE(pb_ts(2)%sf(ix_t%beg:ix_t%end, iy_t%beg:iy_t%end, iz_t%beg:iz_t%end, 1:nnode, 1:nb))'
# 146 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 146 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 146 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 146 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(pb_ts(2)%sf(ix_t%beg:ix_t%end,                 iy_t%beg:iy_t%end,                 iz_t%beg:iz_t%end, 1:nnode, 1:nb))
# 146 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(pb_ts(2)%sf(ix_t%beg:ix_t%end,                 iy_t%beg:iy_t%end,                 iz_t%beg:iz_t%end, 1:nnode, 1:nb))
# 149 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#ifdef MFC_DEBUG
# 149 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 149 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 149 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 149 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:149: ', '@:ALLOCATE(rhs_pb(ix_t%beg:ix_t%end, iy_t%beg:iy_t%end, iz_t%beg:iz_t%end, 1:nnode, 1:nb))'
# 149 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 149 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 149 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 149 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(rhs_pb(ix_t%beg:ix_t%end,                 iy_t%beg:iy_t%end,                 iz_t%beg:iz_t%end, 1:nnode, 1:nb))
# 149 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(rhs_pb(ix_t%beg:ix_t%end,                 iy_t%beg:iy_t%end,                 iz_t%beg:iz_t%end, 1:nnode, 1:nb))
# 152 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        else if (qbmm .and. polytropic) then
#ifdef MFC_DEBUG
# 153 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 153 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 153 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 153 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:153: ', '@:ALLOCATE(pb_ts(1)%sf(ix_t%beg:ix_t%beg + 1, iy_t%beg:iy_t%beg + 1, iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))'
# 153 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 153 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 153 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 153 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(pb_ts(1)%sf(ix_t%beg:ix_t%beg + 1,                 iy_t%beg:iy_t%beg + 1,                 iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
# 153 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(pb_ts(1)%sf(ix_t%beg:ix_t%beg + 1,                 iy_t%beg:iy_t%beg + 1,                 iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
# 156 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#ifdef MFC_DEBUG
# 156 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 156 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 156 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 156 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:156: ', '@:ALLOCATE(pb_ts(2)%sf(ix_t%beg:ix_t%beg + 1, iy_t%beg:iy_t%beg + 1, iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))'
# 156 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 156 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 156 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 156 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(pb_ts(2)%sf(ix_t%beg:ix_t%beg + 1,                 iy_t%beg:iy_t%beg + 1,                 iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
# 156 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(pb_ts(2)%sf(ix_t%beg:ix_t%beg + 1,                 iy_t%beg:iy_t%beg + 1,                 iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
# 159 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#ifdef MFC_DEBUG
# 159 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 159 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 159 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 159 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:159: ', '@:ALLOCATE(rhs_pb(ix_t%beg:ix_t%beg + 1, iy_t%beg:iy_t%beg + 1, iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))'
# 159 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 159 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 159 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 159 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(rhs_pb(ix_t%beg:ix_t%beg + 1,                 iy_t%beg:iy_t%beg + 1,                 iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
# 159 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(rhs_pb(ix_t%beg:ix_t%beg + 1,                 iy_t%beg:iy_t%beg + 1,                 iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
# 162 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        end if

#ifdef MFC_DEBUG
# 164 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 164 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 164 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 164 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:164: ', '@:ALLOCATE(mv_ts(1:2))'
# 164 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 164 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 164 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 164 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(mv_ts(1:2))
# 164 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(mv_ts(1:2))

        if (qbmm .and. (.not. polytropic)) then
#ifdef MFC_DEBUG
# 167 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 167 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 167 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 167 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:167: ', '@:ALLOCATE(mv_ts(1)%sf(ix_t%beg:ix_t%end, iy_t%beg:iy_t%end, iz_t%beg:iz_t%end, 1:nnode, 1:nb))'
# 167 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 167 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 167 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 167 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(mv_ts(1)%sf(ix_t%beg:ix_t%end,                 iy_t%beg:iy_t%end,                 iz_t%beg:iz_t%end, 1:nnode, 1:nb))
# 167 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(mv_ts(1)%sf(ix_t%beg:ix_t%end,                 iy_t%beg:iy_t%end,                 iz_t%beg:iz_t%end, 1:nnode, 1:nb))
# 170 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#ifdef MFC_DEBUG
# 170 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 170 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 170 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 170 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:170: ', '@:ALLOCATE(mv_ts(2)%sf(ix_t%beg:ix_t%end, iy_t%beg:iy_t%end, iz_t%beg:iz_t%end, 1:nnode, 1:nb))'
# 170 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 170 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 170 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 170 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(mv_ts(2)%sf(ix_t%beg:ix_t%end,                 iy_t%beg:iy_t%end,                 iz_t%beg:iz_t%end, 1:nnode, 1:nb))
# 170 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(mv_ts(2)%sf(ix_t%beg:ix_t%end,                 iy_t%beg:iy_t%end,                 iz_t%beg:iz_t%end, 1:nnode, 1:nb))
# 173 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#ifdef MFC_DEBUG
# 173 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 173 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 173 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 173 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:173: ', '@:ALLOCATE(rhs_mv(ix_t%beg:ix_t%end, iy_t%beg:iy_t%end, iz_t%beg:iz_t%end, 1:nnode, 1:nb))'
# 173 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 173 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 173 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 173 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(rhs_mv(ix_t%beg:ix_t%end,                 iy_t%beg:iy_t%end,                 iz_t%beg:iz_t%end, 1:nnode, 1:nb))
# 173 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(rhs_mv(ix_t%beg:ix_t%end,                 iy_t%beg:iy_t%end,                 iz_t%beg:iz_t%end, 1:nnode, 1:nb))
# 176 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        else if (qbmm .and. polytropic) then
#ifdef MFC_DEBUG
# 177 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 177 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 177 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 177 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:177: ', '@:ALLOCATE(mv_ts(1)%sf(ix_t%beg:ix_t%beg + 1, iy_t%beg:iy_t%beg + 1, iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))'
# 177 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 177 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 177 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 177 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(mv_ts(1)%sf(ix_t%beg:ix_t%beg + 1,                 iy_t%beg:iy_t%beg + 1,                 iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
# 177 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(mv_ts(1)%sf(ix_t%beg:ix_t%beg + 1,                 iy_t%beg:iy_t%beg + 1,                 iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
# 180 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#ifdef MFC_DEBUG
# 180 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 180 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 180 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 180 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:180: ', '@:ALLOCATE(mv_ts(2)%sf(ix_t%beg:ix_t%beg + 1, iy_t%beg:iy_t%beg + 1, iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))'
# 180 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 180 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 180 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 180 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(mv_ts(2)%sf(ix_t%beg:ix_t%beg + 1,                 iy_t%beg:iy_t%beg + 1,                 iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
# 180 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(mv_ts(2)%sf(ix_t%beg:ix_t%beg + 1,                 iy_t%beg:iy_t%beg + 1,                 iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
# 183 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#ifdef MFC_DEBUG
# 183 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 183 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 183 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 183 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:183: ', '@:ALLOCATE(rhs_mv(ix_t%beg:ix_t%beg + 1, iy_t%beg:iy_t%beg + 1, iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))'
# 183 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 183 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 183 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 183 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(rhs_mv(ix_t%beg:ix_t%beg + 1,                 iy_t%beg:iy_t%beg + 1,                 iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
# 183 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(rhs_mv(ix_t%beg:ix_t%beg + 1,                 iy_t%beg:iy_t%beg + 1,                 iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
# 186 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        end if

        if (hypoelasticity) then

            do i = stress_idx%beg, stress_idx%end
#ifdef MFC_DEBUG
# 191 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 191 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 191 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 191 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:191: ', '@:ALLOCATE(q_prim_vf(i)%sf(ix_t%beg:ix_t%end, iy_t%beg:iy_t%end, iz_t%beg:iz_t%end))'
# 191 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 191 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 191 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 191 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(q_prim_vf(i)%sf(ix_t%beg:ix_t%end,                     iy_t%beg:iy_t%end,                     iz_t%beg:iz_t%end))
# 191 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(q_prim_vf(i)%sf(ix_t%beg:ix_t%end,                     iy_t%beg:iy_t%end,                     iz_t%beg:iz_t%end))
# 194 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
            end do
        end if

        if (model_eqns == 3) then
            do i = internalEnergies_idx%beg, internalEnergies_idx%end
#ifdef MFC_DEBUG
# 199 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 199 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 199 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 199 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:199: ', '@:ALLOCATE(q_prim_vf(i)%sf(ix_t%beg:ix_t%end, iy_t%beg:iy_t%end, iz_t%beg:iz_t%end))'
# 199 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 199 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 199 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 199 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(q_prim_vf(i)%sf(ix_t%beg:ix_t%end,                     iy_t%beg:iy_t%end,                     iz_t%beg:iz_t%end))
# 199 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(q_prim_vf(i)%sf(ix_t%beg:ix_t%end,                     iy_t%beg:iy_t%end,                     iz_t%beg:iz_t%end))
# 202 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
            end do
        end if

        ! Allocating the cell-average RHS variables
#ifdef MFC_DEBUG
# 206 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 206 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 206 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 206 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:206: ', '@:ALLOCATE(rhs_vf(1:sys_size))'
# 206 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 206 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 206 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 206 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(rhs_vf(1:sys_size))
# 206 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(rhs_vf(1:sys_size))

        do i = 1, sys_size
#ifdef MFC_DEBUG
# 209 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 209 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 209 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 209 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:209: ', '@:ALLOCATE(rhs_vf(i)%sf(0:m, 0:n, 0:p))'
# 209 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 209 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 209 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 209 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    allocate(rhs_vf(i)%sf(0:m, 0:n, 0:p))
# 209 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc enter data create(rhs_vf(i)%sf(0:m, 0:n, 0:p))
        end do

        ! Opening and writing the header of the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_open_run_time_information_file()
        end if

    end subroutine s_initialize_time_steppers_module ! ---------------------

    !> 1st order TVD RK time-stepping algorithm
        !! @param t_step Current time step
    subroutine s_1st_order_tvd_rk(t_step, time_avg) ! --------------------------------

        integer, intent(IN) :: t_step
        real(kind(0d0)), intent(INOUT) :: time_avg

        integer :: i, j, k, l, q!< Generic loop iterator
        real(kind(0d0)) :: start, finish

        ! Stage 1 of 1 =====================================================

        call cpu_time(start)

        call nvtxStartRange("Time_Step")

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step)

#ifdef DEBUG
        print *, 'got rhs'
#endif

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

#ifdef DEBUG
        print *, 'wrote runtime info'
#endif

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

!$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                            q_cons_ts(1)%vf(i)%sf(j, k, l) &
                            + dt*rhs_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do
        !Evolve pb and mv for non-polytropic qbmm
        if (qbmm .and. (.not. polytropic)) then
!$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(1)%sf(j, k, l, q, i) = &
                                    pb_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_pb(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
!$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(1)%sf(j, k, l, q, i) = &
                                    mv_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_mv(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        call nvtxEndRange

        call cpu_time(finish)

        if (t_step >= 4) then
            time_avg = (abs(finish - start) + (t_step - 4)*time_avg)/(t_step - 3)
        else
            time_avg = 0d0
        end if

        ! ==================================================================

    end subroutine s_1st_order_tvd_rk ! ------------------------------------

    !> 2nd order TVD RK time-stepping algorithm
        !! @param t_step Current time-step
    subroutine s_2nd_order_tvd_rk(t_step, time_avg) ! --------------------------------

        integer, intent(IN) :: t_step
        real(kind(0d0)), intent(INOUT) :: time_avg

        integer :: i, j, k, l, q!< Generic loop iterator
        real(kind(0d0)) :: start, finish

        ! Stage 1 of 2 =====================================================

        call cpu_time(start)

        call nvtxStartRange("Time_Step")

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step)

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

!$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(2)%vf(i)%sf(j, k, l) = &
                            q_cons_ts(1)%vf(i)%sf(j, k, l) &
                            + dt*rhs_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do
        !Evolve pb and mv for non-polytropic qbmm
        if (qbmm .and. (.not. polytropic)) then
!$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(2)%sf(j, k, l, q, i) = &
                                    pb_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_pb(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
!$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(2)%sf(j, k, l, q, i) = &
                                    mv_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_mv(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)

        ! ==================================================================

        ! Stage 2 of 2 =====================================================

        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, pb_ts(2)%sf, rhs_pb, mv_ts(2)%sf, rhs_mv, t_step)

!$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                            (q_cons_ts(1)%vf(i)%sf(j, k, l) &
                             + q_cons_ts(2)%vf(i)%sf(j, k, l) &
                             + dt*rhs_vf(i)%sf(j, k, l))/2d0
                    end do
                end do
            end do
        end do

        if (qbmm .and. (.not. polytropic)) then
!$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(1)%sf(j, k, l, q, i) = &
                                    (pb_ts(1)%sf(j, k, l, q, i) &
                                     + pb_ts(2)%sf(j, k, l, q, i) &
                                     + dt*rhs_pb(j, k, l, q, i))/2d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
!$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(1)%sf(j, k, l, q, i) = &
                                    (mv_ts(1)%sf(j, k, l, q, i) &
                                     + mv_ts(2)%sf(j, k, l, q, i) &
                                     + dt*rhs_mv(j, k, l, q, i))/2d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        call nvtxEndRange

        call cpu_time(finish)

        if (t_step >= 4) then
            time_avg = (abs(finish - start) + (t_step - 4)*time_avg)/(t_step - 3)
        else
            time_avg = 0d0
        end if

        ! ==================================================================

    end subroutine s_2nd_order_tvd_rk ! ------------------------------------

    !> 3rd order TVD RK time-stepping algorithm
        !! @param t_step Current time-step
    subroutine s_3rd_order_tvd_rk(t_step, time_avg) ! --------------------------------

        integer, intent(IN) :: t_step
        real(kind(0d0)), intent(INOUT) :: time_avg

        integer :: i, j, k, l, q
        real(kind(0d0)) :: ts_error, denom, error_fraction, time_step_factor !< Generic loop iterator
        real(kind(0d0)) :: start, finish

        ! Stage 1 of 3 =====================================================

        call cpu_time(start)

        call nvtxStartRange("Time_Step")

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step)

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

!$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(2)%vf(i)%sf(j, k, l) = &
                            q_cons_ts(1)%vf(i)%sf(j, k, l) &
                            + dt*rhs_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do
        !Evolve pb and mv for non-polytropic qbmm
        if (qbmm .and. (.not. polytropic)) then
!$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(2)%sf(j, k, l, q, i) = &
                                    pb_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_pb(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
!$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(2)%sf(j, k, l, q, i) = &
                                    mv_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_mv(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)
        if ( model_eqns == 3 .and. (.not.relax ) ) then
            call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
        end if

        ! ==================================================================

        ! Stage 2 of 3 =====================================================

        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, pb_ts(2)%sf, rhs_pb, mv_ts(2)%sf, rhs_mv, t_step)

!$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(2)%vf(i)%sf(j, k, l) = &
                            (3d0*q_cons_ts(1)%vf(i)%sf(j, k, l) &
                             + q_cons_ts(2)%vf(i)%sf(j, k, l) &
                             + dt*rhs_vf(i)%sf(j, k, l))/4d0
                    end do
                end do
            end do
        end do

        if (qbmm .and. (.not. polytropic)) then
!$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(2)%sf(j, k, l, q, i) = &
                                    (3d0*pb_ts(1)%sf(j, k, l, q, i) &
                                     + pb_ts(2)%sf(j, k, l, q, i) &
                                     + dt*rhs_pb(j, k, l, q, i))/4d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
!$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(2)%sf(j, k, l, q, i) = &
                                    (3d0*mv_ts(1)%sf(j, k, l, q, i) &
                                     + mv_ts(2)%sf(j, k, l, q, i) &
                                     + dt*rhs_mv(j, k, l, q, i))/4d0
                            end do
                        end do
                    end do
                end do
            end do
        end if
        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)
        if ( model_eqns == 3 .and. (.not.relax ) ) then
            call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
        end if
        
        ! ==================================================================

        ! Stage 3 of 3 =====================================================
        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, pb_ts(2)%sf, rhs_pb, mv_ts(2)%sf, rhs_mv, t_step)

!$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                            (q_cons_ts(1)%vf(i)%sf(j, k, l) &
                             + 2d0*q_cons_ts(2)%vf(i)%sf(j, k, l) &
                             + 2d0*dt*rhs_vf(i)%sf(j, k, l))/3d0
                    end do
                end do
            end do
        end do

        if (qbmm .and. (.not. polytropic)) then
!$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(1)%sf(j, k, l, q, i) = &
                                    (pb_ts(1)%sf(j, k, l, q, i) &
                                     + 2d0*pb_ts(2)%sf(j, k, l, q, i) &
                                     + 2d0*dt*rhs_pb(j, k, l, q, i))/3d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
!$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(1)%sf(j, k, l, q, i) = &
                                    (mv_ts(1)%sf(j, k, l, q, i) &
                                     + 2d0*mv_ts(2)%sf(j, k, l, q, i) &
                                     + 2d0*dt*rhs_mv(j, k, l, q, i))/3d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)      
        if ( model_eqns == 3 .and. (.not.relax ) ) then
            call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)
        end if

        call nvtxEndRange

        call cpu_time(finish)

        time = time + (finish - start)

        if (t_step >= 4) then
            time_avg = (abs(finish - start) + (t_step - 4)*time_avg)/(t_step - 3)
        else
            time_avg = 0d0
        end if

        ! ==================================================================

    end subroutine s_3rd_order_tvd_rk ! ------------------------------------

    !> This subroutine saves the temporary q_prim_vf vector
        !!      into the q_prim_ts vector that is then used in p_main
        !! @param t_step current time-step
    subroutine s_time_step_cycling(t_step) ! ----------------------------

        integer, intent(IN) :: t_step

        integer :: i !< Generic loop iterator

        do i = 1, sys_size
!$acc update host(q_prim_vf(i)%sf)
        end do

        if (t_step == t_step_start) then
            do i = 1, sys_size
                q_prim_ts(3)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 1) then
            do i = 1, sys_size
                q_prim_ts(2)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 2) then
            do i = 1, sys_size
                q_prim_ts(1)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 3) then
            do i = 1, sys_size
                q_prim_ts(0)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        else ! All other timesteps
            do i = 1, sys_size
                q_prim_ts(3)%vf(i)%sf(:, :, :) = q_prim_ts(2)%vf(i)%sf(:, :, :)
                q_prim_ts(2)%vf(i)%sf(:, :, :) = q_prim_ts(1)%vf(i)%sf(:, :, :)
                q_prim_ts(1)%vf(i)%sf(:, :, :) = q_prim_ts(0)%vf(i)%sf(:, :, :)
                q_prim_ts(0)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        end if

    end subroutine s_time_step_cycling ! -----------------------------------

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_time_steppers_module() ! -------------------------

        integer :: i, j !< Generic loop iterators

        ! Deallocating the cell-average conservative variables
        do i = 1, num_ts

            do j = 1, sys_size
#ifdef MFC_DEBUG
# 733 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 733 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 733 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 733 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:733: ', '@:DEALLOCATE(q_cons_ts(i)%vf(j)%sf)'
# 733 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 733 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 733 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 733 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    deallocate(q_cons_ts(i)%vf(j)%sf)
# 733 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc exit data delete(q_cons_ts(i)%vf(j)%sf)
            end do

#ifdef MFC_DEBUG
# 736 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 736 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 736 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 736 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:736: ', '@:DEALLOCATE(q_cons_ts(i)%vf)'
# 736 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 736 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 736 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 736 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    deallocate(q_cons_ts(i)%vf)
# 736 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc exit data delete(q_cons_ts(i)%vf)

        end do

#ifdef MFC_DEBUG
# 740 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 740 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 740 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 740 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:740: ', '@:DEALLOCATE(q_cons_ts)'
# 740 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 740 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 740 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 740 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    deallocate(q_cons_ts)
# 740 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc exit data delete(q_cons_ts)

        ! Deallocating the cell-average primitive ts variables
        if (probe_wrt) then
            do i = 0, 3
                do j = 1, sys_size
#ifdef MFC_DEBUG
# 746 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 746 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 746 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 746 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:746: ', '@:DEALLOCATE(q_prim_ts(i)%vf(j)%sf)'
# 746 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 746 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 746 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 746 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    deallocate(q_prim_ts(i)%vf(j)%sf)
# 746 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc exit data delete(q_prim_ts(i)%vf(j)%sf)
                end do
#ifdef MFC_DEBUG
# 748 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 748 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 748 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 748 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:748: ', '@:DEALLOCATE(q_prim_ts(i)%vf)'
# 748 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 748 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 748 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 748 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    deallocate(q_prim_ts(i)%vf)
# 748 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc exit data delete(q_prim_ts(i)%vf)
            end do
#ifdef MFC_DEBUG
# 750 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 750 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 750 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 750 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:750: ', '@:DEALLOCATE(q_prim_ts)'
# 750 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 750 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 750 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 750 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    deallocate(q_prim_ts)
# 750 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc exit data delete(q_prim_ts)
        end if

        ! Deallocating the cell-average primitive variables
        do i = 1, adv_idx%end
#ifdef MFC_DEBUG
# 755 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 755 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 755 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 755 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:755: ', '@:DEALLOCATE(q_prim_vf(i)%sf)'
# 755 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 755 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 755 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 755 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    deallocate(q_prim_vf(i)%sf)
# 755 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc exit data delete(q_prim_vf(i)%sf)
        end do

        if (hypoelasticity) then
            do i = stress_idx%beg, stress_idx%end
#ifdef MFC_DEBUG
# 760 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 760 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 760 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 760 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:760: ', '@:DEALLOCATE(q_prim_vf(i)%sf)'
# 760 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 760 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 760 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 760 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    deallocate(q_prim_vf(i)%sf)
# 760 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc exit data delete(q_prim_vf(i)%sf)
            end do
        end if

        if (bubbles) then
            do i = bub_idx%beg, bub_idx%end
#ifdef MFC_DEBUG
# 766 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 766 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 766 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 766 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:766: ', '@:DEALLOCATE(q_prim_vf(i)%sf)'
# 766 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 766 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 766 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 766 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    deallocate(q_prim_vf(i)%sf)
# 766 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc exit data delete(q_prim_vf(i)%sf)
            end do
        end if

        if (model_eqns == 3) then
            do i = internalEnergies_idx%beg, internalEnergies_idx%end
#ifdef MFC_DEBUG
# 772 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 772 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 772 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 772 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:772: ', '@:DEALLOCATE(q_prim_vf(i)%sf)'
# 772 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 772 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 772 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 772 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    deallocate(q_prim_vf(i)%sf)
# 772 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc exit data delete(q_prim_vf(i)%sf)
            end do
        end if

#ifdef MFC_DEBUG
# 776 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 776 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 776 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 776 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:776: ', '@:DEALLOCATE(q_prim_vf)'
# 776 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 776 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 776 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 776 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    deallocate(q_prim_vf)
# 776 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc exit data delete(q_prim_vf)

        ! Deallocating the cell-average RHS variables
        do i = 1, sys_size
#ifdef MFC_DEBUG
# 780 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 780 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 780 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 780 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:780: ', '@:DEALLOCATE(rhs_vf(i)%sf)'
# 780 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 780 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 780 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 780 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    deallocate(rhs_vf(i)%sf)
# 780 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc exit data delete(rhs_vf(i)%sf)
        end do

#ifdef MFC_DEBUG
# 783 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    block
# 783 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        use iso_fortran_env, only: output_unit
# 783 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"

# 783 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        print *, 'm_time_steppers.fpp:783: ', '@:DEALLOCATE(rhs_vf)'
# 783 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
        call flush(output_unit)
# 783 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    end block
# 783 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
#endif
# 783 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    deallocate(rhs_vf)
# 783 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/simulation/m_time_steppers.fpp"
    !$acc exit data delete(rhs_vf)

        ! Writing the footer of and closing the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_close_run_time_information_file()
        end if

    end subroutine s_finalize_time_steppers_module ! -----------------------

end module m_time_steppers
