# 1 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
# 10 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"

# 23 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"

# 32 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"

# 38 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"

# 44 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"

program syscheck

#ifdef MFC_MPI
# 47 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    use mpi
# 47 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
#ifdef MFC_OpenACC
# 48 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    use openacc
# 48 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif

    implicit none

    integer :: ierr, rank = 0, nRanks = 1

#ifdef MFC_OpenACC
# 54 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    integer(acc_device_kind) :: devtype
# 54 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
#ifdef MFC_OpenACC
# 55 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    integer :: i,num_devices
# 55 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
#ifdef MFC_OpenACC
# 56 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    real(kind(0d0)),allocatable,dimension(:) :: arr
# 56 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
#ifdef MFC_OpenACC
# 57 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    integer,parameter :: N = 100
# 57 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif

#ifdef MFC_MPI
# 59 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 59 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 59 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 59 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[TEST] MPI: call mpi_init(ierr)"
# 59 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 59 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 59 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 59 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    call mpi_init(ierr)
# 59 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (ierr /= MPI_SUCCESS) then
# 59 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
        print*, " -> Error: ", ierr
# 59 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
        stop ierr
# 59 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 59 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#else
# 59 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 59 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 59 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 59 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[SKIP] MPI: call mpi_init(ierr)"
# 59 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 59 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 59 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 59 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
#ifdef MFC_MPI
# 60 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 60 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 60 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 60 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[TEST] MPI: call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)"
# 60 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 60 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 60 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 60 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
# 60 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (ierr /= MPI_SUCCESS) then
# 60 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
        print*, " -> Error: ", ierr
# 60 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
        stop ierr
# 60 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 60 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#else
# 60 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 60 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 60 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 60 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[SKIP] MPI: call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)"
# 60 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 60 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 60 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 60 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
#ifdef MFC_MPI
# 61 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 61 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 61 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 61 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[TEST] MPI: call mpi_barrier(MPI_COMM_WORLD, ierr)"
# 61 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 61 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 61 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 61 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    call mpi_barrier(MPI_COMM_WORLD, ierr)
# 61 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (ierr /= MPI_SUCCESS) then
# 61 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
        print*, " -> Error: ", ierr
# 61 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
        stop ierr
# 61 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 61 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#else
# 61 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 61 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 61 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 61 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[SKIP] MPI: call mpi_barrier(MPI_COMM_WORLD, ierr)"
# 61 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 61 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 61 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 61 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
#ifdef MFC_MPI
# 62 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 62 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 62 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 62 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[TEST] MPI: call assert(rank >= 0)"
# 62 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 62 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 62 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 62 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    call assert(rank >= 0)
# 62 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (ierr /= MPI_SUCCESS) then
# 62 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
        print*, " -> Error: ", ierr
# 62 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
        stop ierr
# 62 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 62 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#else
# 62 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 62 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 62 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 62 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[SKIP] MPI: call assert(rank >= 0)"
# 62 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 62 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 62 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 62 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
#ifdef MFC_MPI
# 63 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 63 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 63 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 63 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[TEST] MPI: call mpi_comm_size(MPI_COMM_WORLD, nRanks, ierr)"
# 63 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 63 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 63 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 63 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    call mpi_comm_size(MPI_COMM_WORLD, nRanks, ierr)
# 63 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (ierr /= MPI_SUCCESS) then
# 63 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
        print*, " -> Error: ", ierr
# 63 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
        stop ierr
# 63 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 63 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#else
# 63 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 63 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 63 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 63 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[SKIP] MPI: call mpi_comm_size(MPI_COMM_WORLD, nRanks, ierr)"
# 63 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 63 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 63 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 63 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
#ifdef MFC_MPI
# 64 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 64 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 64 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 64 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[TEST] MPI: call assert(nRanks > 0 .and. rank < nRanks)"
# 64 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 64 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 64 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 64 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    call assert(nRanks > 0 .and. rank < nRanks)
# 64 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (ierr /= MPI_SUCCESS) then
# 64 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
        print*, " -> Error: ", ierr
# 64 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
        stop ierr
# 64 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 64 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#else
# 64 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 64 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 64 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 64 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[SKIP] MPI: call assert(nRanks > 0 .and. rank < nRanks)"
# 64 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 64 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 64 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 64 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif

#ifdef MFC_OpenACC
# 66 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 66 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 66 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 66 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[TEST] ACC: devtype = acc_get_device_type()"
# 66 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 66 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 66 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 66 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    devtype = acc_get_device_type()
# 66 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#else
# 66 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 66 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 66 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 66 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[SKIP] ACC: devtype = acc_get_device_type()"
# 66 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 66 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 66 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 66 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
#ifdef MFC_OpenACC
# 67 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 67 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 67 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 67 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[TEST] ACC: num_devices = acc_get_num_devices(devtype)"
# 67 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 67 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 67 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 67 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    num_devices = acc_get_num_devices(devtype)
# 67 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#else
# 67 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 67 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 67 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 67 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[SKIP] ACC: num_devices = acc_get_num_devices(devtype)"
# 67 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 67 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 67 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 67 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
#ifdef MFC_OpenACC
# 68 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 68 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 68 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 68 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[TEST] ACC: call assert(num_devices > 0)"
# 68 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 68 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 68 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 68 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    call assert(num_devices > 0)
# 68 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#else
# 68 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 68 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 68 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 68 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[SKIP] ACC: call assert(num_devices > 0)"
# 68 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 68 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 68 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 68 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
#ifdef MFC_OpenACC
# 69 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 69 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 69 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 69 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[TEST] ACC: call acc_set_device_num(mod(rank, nRanks), devtype)"
# 69 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 69 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 69 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 69 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    call acc_set_device_num(mod(rank, nRanks), devtype)
# 69 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#else
# 69 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 69 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 69 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 69 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[SKIP] ACC: call acc_set_device_num(mod(rank, nRanks), devtype)"
# 69 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 69 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 69 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 69 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
#ifdef MFC_OpenACC
# 70 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 70 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 70 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 70 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[TEST] ACC: allocate(arr(1:N))"
# 70 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 70 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 70 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 70 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    allocate(arr(1:N))
# 70 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#else
# 70 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 70 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 70 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 70 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[SKIP] ACC: allocate(arr(1:N))"
# 70 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 70 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 70 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 70 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
#ifdef MFC_OpenACC
# 71 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 71 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 71 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 71 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[TEST] ACC: !$acc enter data create(arr(1:N))"
# 71 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 71 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 71 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 71 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    !$acc enter data create(arr(1:N))
# 71 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#else
# 71 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 71 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 71 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 71 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[SKIP] ACC: !$acc enter data create(arr(1:N))"
# 71 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 71 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 71 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 71 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
#ifdef MFC_OpenACC
# 72 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 72 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 72 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 72 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[TEST] ACC: !$acc parallel loop"
# 72 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 72 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 72 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 72 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    !$acc parallel loop
# 72 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#else
# 72 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 72 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 72 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 72 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[SKIP] ACC: !$acc parallel loop"
# 72 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 72 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 72 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 72 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
#ifdef MFC_OpenACC
# 73 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    do i = 1,N
# 73 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
#ifdef MFC_OpenACC
# 74 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    arr(i) = i
# 74 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
#ifdef MFC_OpenACC
# 75 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end do
# 75 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
#ifdef MFC_OpenACC
# 76 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 76 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 76 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 76 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[TEST] ACC: !$acc update host(arr(1:N))"
# 76 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 76 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 76 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 76 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    !$acc update host(arr(1:N))
# 76 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#else
# 76 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 76 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 76 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 76 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[SKIP] ACC: !$acc update host(arr(1:N))"
# 76 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 76 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 76 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 76 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
#ifdef MFC_OpenACC
# 77 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 77 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 77 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 77 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[TEST] ACC: !$acc exit data delete(arr)"
# 77 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 77 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 77 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 77 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    !$acc exit data delete(arr)
# 77 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#else
# 77 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 77 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 77 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 77 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[SKIP] ACC: !$acc exit data delete(arr)"
# 77 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 77 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 77 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 77 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif

#ifdef MFC_MPI
# 79 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 79 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 79 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 79 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[TEST] MPI: call mpi_barrier(MPI_COMM_WORLD, ierr)"
# 79 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 79 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 79 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 79 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    call mpi_barrier(MPI_COMM_WORLD, ierr)
# 79 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (ierr /= MPI_SUCCESS) then
# 79 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
        print*, " -> Error: ", ierr
# 79 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
        stop ierr
# 79 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 79 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#else
# 79 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 79 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 79 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 79 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[SKIP] MPI: call mpi_barrier(MPI_COMM_WORLD, ierr)"
# 79 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 79 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 79 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 79 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
#ifdef MFC_MPI
# 80 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 80 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 80 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 80 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[TEST] MPI: call mpi_finalize(ierr)"
# 80 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 80 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 80 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 80 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    call mpi_finalize(ierr)
# 80 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (ierr /= MPI_SUCCESS) then
# 80 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
        print*, " -> Error: ", ierr
# 80 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
        stop ierr
# 80 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 80 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#else
# 80 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 80 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 80 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 80 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "[SKIP] MPI: call mpi_finalize(ierr)"
# 80 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 80 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 80 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 80 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif

#ifdef MFC_MPI
# 82 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 82 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 82 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, ""
# 82 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 82 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 82 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
#ifdef MFC_MPI
# 83 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    if (rank == 0) then
# 83 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif
# 83 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    print *, "Syscheck: PASSED."
# 83 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#ifdef MFC_MPI
# 83 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
    end if
# 83 "/ocean/projects/phy230019p/jrchreim/MFC-JRChreim/src/syscheck/syscheck.fpp"
#endif

end program syscheck

subroutine assert(condition)
    
    use iso_fortran_env, only: output_unit, error_unit
    
    logical, intent(in) :: condition
    
    if (.not. condition) then
        call flush(int(output_unit))
        call flush(int(error_unit))
        stop 1
    end if

end subroutine assert
