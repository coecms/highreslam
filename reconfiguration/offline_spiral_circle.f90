!> 
!! Copyright 2019 Scott Wales
!!
!! \author  Scott Wales <scott.wales@unimelb.edu.au>
!!
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!! You may obtain a copy of the License at
!!
!!     http://www.apache.org/licenses/LICENSE-2.0
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.

module offline_spiral_circle
    use iso_c_binding
    public

    type, public :: spiral_circle_data
        logical(kind=c_bool), allocatable :: lsm(:), unres(:)
        integer(kind=c_int64_t), allocatable :: input_idx(:)
        real(kind=c_double), allocatable :: lat(:), lon(:)

        logical(kind=c_bool) :: is_land_field, constrained
        logical(kind=c_bool) :: cyclic_domain

        real(kind=c_double) :: constrained_max_dist, dist_step
        real(kind=c_double) :: planet_radius
    end type

    interface read_var
        module procedure :: read_var_r
        module procedure :: read_var_i
        module procedure :: read_var_l
    end interface

contains
    subroutine error(message)
        use ifcore
        use mpi_f08
        character(len=*), intent(in) :: message
        call tracebackqq(message, user_exit_code=0)
        call mpi_abort(MPI_COMM_WORLD, 10)
    end subroutine

    subroutine nfe(err)
        use netcdf4_f03
        integer, intent(in) :: err

        if (err /= NF_NOERR) then
            call error(nf_strerror(err))
        end if
    end subroutine

    subroutine read_var_i(ncid, dimname, varname, var)
        use netcdf4_f03
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: dimname, varname
        integer(kind=8), allocatable, intent(out) :: var(:)

        integer :: dimid, len, varid
        
        call nfe(nf_inq_dimid(ncid, dimname, dimid))
        call nfe(nf_inq_dimlen(ncid, dimid, len))
        allocate(var(len))
        call nfe(nf_inq_varid(ncid, varname, varid))
        call nfe(nf_get_var(ncid, varid, var))
    end subroutine
    subroutine read_var_r(ncid, dimname, varname, var)
        use netcdf4_f03
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: dimname, varname
        real(kind=8), allocatable, intent(out) :: var(:)

        integer :: dimid, len, varid
        
        call nfe(nf_inq_dimid(ncid, dimname, dimid))
        call nfe(nf_inq_dimlen(ncid, dimid, len))
        allocate(var(len))
        call nfe(nf_inq_varid(ncid, varname, varid))
        call nfe(nf_get_var(ncid, varid, var))
    end subroutine
    subroutine read_var_l(ncid, dimname, varname, var)
        use netcdf4_f03
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: dimname, varname
        logical(kind=c_bool), allocatable, intent(out) :: var(:)

        integer :: dimid, len, varid
        
        call nfe(nf_inq_dimid(ncid, dimname, dimid))
        call nfe(nf_inq_dimlen(ncid, dimid, len))
        allocate(var(len))
        call nfe(nf_inq_varid(ncid, varname, varid))
        call nfe(nf_get_var(ncid, varid, var))
    end subroutine

    subroutine read_netcdf(path, sc)
        use netcdf4_f03
        character(len=*), intent(in) :: path
        type(spiral_circle_data), intent(out) :: sc

        integer :: ncid

        call nfe(nf_open(path, NF_NOWRITE, ncid))

        call read_var(ncid, 'phi', 'lat', sc%lat)
        call read_var(ncid, 'lambda', 'lon', sc%lon)
        call read_var(ncid, 'flat', 'lsm', sc%lsm)
        call read_var(ncid, 'flat', 'unres_mask', sc%unres)
        call read_var(ncid, 'unres', 'input', sc%input_idx)

        call nfe(nf_get_att(ncid, NF_GLOBAL, 'is_land_field', sc%is_land_field))
        call nfe(nf_get_att(ncid, NF_GLOBAL, 'constrained', sc%constrained))
        call nfe(nf_get_att(ncid, NF_GLOBAL, 'cyclic_domain', sc%cyclic_domain))

        call nfe(nf_get_att(ncid, NF_GLOBAL, 'constrained_max_dist', sc%constrained_max_dist))
        call nfe(nf_get_att(ncid, NF_GLOBAL, 'dist_step', sc%dist_step))
        call nfe(nf_get_att(ncid, NF_GLOBAL, 'planet_radius', sc%planet_radius))

        sc%constrained_max_dist = 1000

        call nfe(nf_close(ncid))
        
    end subroutine

    subroutine spiral_circle(sc, output_idx)
        use mpi_f08
        use f_shum_spiral_search_mod
        type(spiral_circle_data), intent(in) :: sc
        integer(kind=8), allocatable, intent(out) :: output_idx(:)
        integer(kind=8), allocatable :: buffer(:)
        integer(kind=8) :: err, nlat, nlon

        character(len=1024) :: cmessage

        integer :: npes, rank, offset, i
        integer(kind=8) :: mycount
        integer, allocatable :: recvcounts(:), displs(:)

        call mpi_comm_size(mpi_comm_world, npes)
        call mpi_comm_rank(mpi_comm_world, rank)

        if (rank == 0) then
            write(*,*) sc%lsm(1:5)
            write(*,*) sc%unres(1:5)
            write(*,*) sc%input_idx(1:5)
            write(*,*) sc%lat(1:5)
            write(*,*) sc%lon(1:5)
            write(*,*) count(sc%lsm), size(sc%lsm)
            write(*,*) count(sc%unres), size(sc%unres)
            write(*,*) sc%is_land_field, sc%constrained, sc%cyclic_domain
        end if

        allocate(recvcounts(npes), displs(npes))
        recvcounts = size(sc%input_idx) / npes 
        displs(1) = 1
        do i=2,npes
            displs(i) = displs(i-1) + recvcounts(i-1)
        end do
        recvcounts(npes) = size(sc%input_idx) - displs(npes)

        if (size(sc%lat) * size(sc%lon) /= size(sc%lsm)) then
            call error("Bad sizes")
        end if

        offset = displs(rank+1)
        mycount = recvcounts(rank+1)

        allocate(buffer(mycount))

        nlat = size(sc%lat)
        nlon = size(sc%lon)

        err = f_shum_spiral_search_algorithm( &
            sc%lsm(:), sc%input_idx(offset:offset+mycount), mycount, &
            nlat, nlon, sc%lat, sc%lon, &
            sc%is_land_field, sc%constrained, sc%constrained_max_dist, &
            sc%dist_step, sc%cyclic_domain, sc%unres, &
            buffer, sc%planet_radius, cmessage &
        )
        write(*,*) rank, 'done'

        if (err /= 0) then
            write(*,*) rank, trim(cmessage)
        end if

        allocate(output_idx(size(sc%input_idx)))
        call MPI_Gatherv(buffer, recvcounts(rank+1), MPI_INTEGER8, output_idx, recvcounts, &
            displs, MPI_INTEGER8, 0, MPI_COMM_WORLD)

        deallocate(buffer)
    end subroutine

    subroutine write_netcdf(path, output_idx)
        use mpi_f08
        use netcdf4_f03
        character(len=*), intent(in) :: path
        integer(kind=8), intent(in) :: output_idx(:)
        integer :: ncid, varid
        integer :: rank

        call mpi_comm_rank(mpi_comm_world, rank)

        if (rank == 0) then
            call nfe(nf_open(path, NF_WRITE, ncid))
            call nfe(nf_inq_varid(ncid, 'output', varid))
            call nfe(nf_put_var(ncid, varid, output_idx))
            call nfe(nf_close(ncid))
        end if
    end subroutine

end module

