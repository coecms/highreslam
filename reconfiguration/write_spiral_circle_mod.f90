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

module external_spiral_circle_mod
contains

subroutine nferr(err)
USE Ereport_Mod, ONLY: Ereport
use netcdf
integer(kind=4), intent(in) :: err
integer :: ex
if (err/=nf90_NOERR) then
    ex = 10
    call ereport('external_spiral_circle', ex, nf90_strerror(err))
end if
end subroutine

subroutine write_int_mask(filename, int_lsm_in, int_lsm_out, lsm_in, lsm_out, local_lsm_out)
use netcdf
use um_parcore, only: mype
use rcf_field_type_mod, only: field_type
use mpl
implicit none
character(len=*), intent(in) :: filename
integer(kind=8), intent(in) :: int_lsm_in(:), int_lsm_out(:)
logical(kind=8), intent(in) :: local_lsm_out(:)
type(field_type), intent(in) :: lsm_in, lsm_out

integer(kind=8) :: err
integer(kind=4) :: points, dim_x, dim_y, var_in, var_out, var_local, ncid
integer(kind=4), allocatable :: data_local(:)
character(len=4) :: mpi_ext

write (mpi_ext, '(I0.4)') mype

call nferr(nf90_create(trim(filename) // mpi_ext, ior( nf90_CLOBBER, nf90_NETCDF4), ncid))

points = lsm_in%glob_row_len
call nferr(nf90_def_dim(ncid, 'x_in', points, dim_x))
points = lsm_in%glob_rows
call nferr(nf90_def_dim(ncid, 'y_in', points, dim_y))
call nferr(nf90_def_var(ncid, 'mask_in', nf90_byte, [dim_x, dim_y], var_in))

points = lsm_out%glob_row_len
call nferr(nf90_def_dim(ncid, 'x_out', points, dim_x))
points = lsm_out%glob_rows
call nferr(nf90_def_dim(ncid, 'y_out', points, dim_y))
call nferr(nf90_def_var(ncid, 'mask_out', nf90_byte, [dim_x, dim_y], var_out))

points = lsm_out%row_len
call nferr(nf90_def_dim(ncid, 'x_local', points, dim_x))
points = lsm_out%rows
call nferr(nf90_def_dim(ncid, 'y_local', points, dim_y))
call nferr(nf90_def_var(ncid, 'mask_local', nf90_byte, [dim_x, dim_y], var_local))

call nferr(nf90_enddef(ncid))

call nferr(nf90_put_var(ncid, var_in, reshape(int_lsm_in,[lsm_in%glob_row_len, lsm_in%glob_rows])))
call nferr(nf90_put_var(ncid, var_out, reshape(int_lsm_out,[lsm_out%glob_row_len, lsm_out%glob_rows])))

allocate(data_local(lsm_out%row_len*lsm_out%rows))
data_local = 0
where(local_lsm_out) data_local = 1
call nferr(nf90_put_var(ncid, var_local, reshape(data_local,[lsm_out%row_len, lsm_out%rows])))
call nferr(nf90_close(ncid))

call mpl_barrier(mpl_comm_world, err)
end subroutine

subroutine write_spiral_circle_s                                                   &
           (filename, lsm, index_unres, no_point_unres,                                  &
            points_phi, points_lambda, lats, lons,                             &
            is_land_field, constrained,                   &
             cyclic_domain, unres_mask,                              &
            indices)

use iso_c_binding, only: c_bool
use netcdf
use um_parcore, only: mype
use planet_constants_mod, only: planet_radius
use mpl
IMPLICIT NONE
character(len=*), intent(in) :: filename

                                          ! Number of rows in grid
INTEGER(KIND=8), INTENT(IN)    :: points_phi
                                          ! Number of columns in grid
INTEGER(KIND=8), INTENT(IN)    :: points_lambda
                                          ! Land sea mask
LOGICAL(KIND=8),  INTENT(IN)    :: lsm(points_lambda*points_phi)
                                          ! Number of unresolved points
INTEGER(KIND=8), INTENT(IN)    :: no_point_unres
                                          ! Index to unresolved pts
INTEGER(KIND=8), INTENT(IN)    :: index_unres(no_point_unres)
                                          ! Latitudes
REAL(KIND=8),   INTENT(IN)    :: lats(points_phi)
                                          ! Longitudes
REAL(KIND=8),   INTENT(IN)    :: lons(points_lambda)
                                          ! False for sea, True for land field
LOGICAL(KIND=8),  INTENT(IN)    :: is_land_field
                                    ! True to apply constraint distance
LOGICAL(KIND=8),  INTENT(IN)    :: constrained
                                    ! Contraint distance (m)
LOGICAL(KIND=8),  INTENT(IN)    :: cyclic_domain
                                          ! False for a point that is resolved,
                                          ! True for an unresolved point
LOGICAL(KIND=8),  INTENT(IN)    :: unres_mask(points_lambda*points_phi)
                                          ! Indices to resolved pts
INTEGER(KIND=8), INTENT(OUT)   :: indices(no_point_unres)
                                          ! Radius of planet (in m)

integer(kind=8) :: err
integer(kind=4) :: ncid, dim_phi, dim_lambda, dim_unres, dim_flat
integer(kind=4) :: var_lsm, var_lat, var_lon, var_input, var_output, var_unres
integer(kind=4) :: deflate, chunk_1d(1), chunk_2d(2)
logical(kind=4) :: shuffle
integer(kind=4) :: points_phi4, points_lambda4, no_point_unres4, points_flat

integer(kind=4), allocatable :: mask(:)
integer(kind=4) :: switch
character(len=4) :: mpi_ext

write (mpi_ext, '(I0.4)') mype

points_phi4 = points_phi
points_lambda4 = points_lambda
no_point_unres4 = no_point_unres
points_flat = points_phi * points_lambda

deflate = 4
chunk_1d = 100000
chunk_2d = [1000, 1000]
shuffle = .true.

call nferr(nf90_create(trim(filename) // mpi_ext, ior( nf90_CLOBBER, nf90_NETCDF4), ncid))
call nferr(nf90_def_dim(ncid, 'phi', points_phi4, dim_phi))
call nferr(nf90_def_dim(ncid, 'lambda', points_lambda4, dim_lambda))
call nferr(nf90_def_dim(ncid, 'unres', no_point_unres4, dim_unres))
call nferr(nf90_def_dim(ncid, 'flat', points_flat, dim_flat))

call nferr(nf90_def_var(ncid, 'lsm', nf90_BYTE, dim_flat, var_lsm &
    ))!, contiguous=nf90_CHUNKED))!, shuffle=shuffle, chunksizes=chunk_2d, deflate_level=deflate)
call nferr(nf90_def_var(ncid, 'unres_mask', nf90_BYTE, dim_flat, var_unres &
    ))!contiguous=nf90_CHUNKED, shuffle=shuffle, chunksizes=chunk_2d, deflate_level=deflate)
call nferr(nf90_def_var(ncid, 'lat', nf90_DOUBLE, dim_phi, var_lat &
    ))!contiguous=nf90_CHUNKED, shuffle=shuffle, chunksizes=chunk_1d, deflate_level=deflate)
call nferr(nf90_def_var(ncid, 'lon', nf90_DOUBLE, dim_lambda, var_lon &
    ))!contiguous=nf90_CHUNKED, shuffle=shuffle, chunksizes=chunk_1d, deflate_level=deflate)
call nferr(nf90_def_var(ncid, 'input', nf90_INT64, dim_unres, var_input &
    ))!contiguous=nf90_CHUNKED, shuffle=shuffle, chunksizes=chunk_1d, deflate_level=deflate)
call nferr(nf90_def_var(ncid, 'output', nf90_INT64, dim_unres, var_output &
    ))!contiguous=nf90_CHUNKED, shuffle=shuffle, chunksizes=chunk_1d, deflate_level=deflate)

call nferr(nf90_enddef(ncid))

allocate(mask(size(lsm)))
mask = 0
where(lsm) mask = 1
call nferr(nf90_put_var(ncid, var_lsm, mask))
mask = 0
where(unres_mask) mask = 1
call nferr(nf90_put_var(ncid, var_unres, mask))
deallocate(mask)

call nferr(nf90_put_var(ncid, var_lat, lats))
call nferr(nf90_put_var(ncid, var_lon, lons))
call nferr(nf90_put_var(ncid, var_input, index_unres))

switch = 0
if (is_land_field) switch = 1
call nferr(nf90_put_att(ncid, nf90_GLOBAL, 'is_land_field', switch))
switch = 0
if (constrained) switch = 1
call nferr(nf90_put_att(ncid, nf90_GLOBAL, 'constrained', switch))
switch = 0
if (cyclic_domain) switch = 1
call nferr(nf90_put_att(ncid, nf90_GLOBAL, 'cyclic_domain', switch))

call nferr(nf90_put_att(ncid, nf90_GLOBAL, 'constrained_max_dist', 200000.0))
call nferr(nf90_put_att(ncid, nf90_GLOBAL, 'dist_step', 3.0))
call nferr(nf90_put_att(ncid, nf90_GLOBAL, 'planet_radius', planet_radius))

call nferr(nf90_close(ncid))

call mpl_barrier(mpl_comm_world, err)

end subroutine

subroutine read_spiral_circle_s &
           (filename, lsm, index_unres, no_point_unres,                                  &
            points_phi, points_lambda, lats, lons,                             &
            is_land_field, constrained,                   &
             cyclic_domain, unres_mask,                              &
            indices)

use netcdf
use um_parcore, only: mype
use mpl
use iso_c_binding, only: c_bool
USE Ereport_Mod, ONLY: Ereport
IMPLICIT NONE
character(len=*), intent(in) :: filename

                                          ! Number of rows in grid
INTEGER(KIND=8), INTENT(IN)    :: points_phi
                                          ! Number of columns in grid
INTEGER(KIND=8), INTENT(IN)    :: points_lambda
                                          ! Land sea mask
LOGICAL(KIND=8),  INTENT(IN)    :: lsm(points_lambda*points_phi)
                                          ! Number of unresolved points
INTEGER(KIND=8), INTENT(IN)    :: no_point_unres
                                          ! Index to unresolved pts
INTEGER(KIND=8), INTENT(IN)    :: index_unres(no_point_unres)
                                          ! Latitudes
REAL(KIND=8),   INTENT(IN)    :: lats(points_phi)
                                          ! Longitudes
REAL(KIND=8),   INTENT(IN)    :: lons(points_lambda)
                                          ! False for sea, True for land field
LOGICAL(KIND=8),  INTENT(IN)    :: is_land_field
                                          ! True to apply constraint distance
LOGICAL(KIND=8),  INTENT(IN)    :: constrained
                                          ! Contraint distance (m)
LOGICAL(KIND=8),  INTENT(IN)    :: cyclic_domain
                                          ! False for a point that is resolved,
                                          ! True for an unresolved point
LOGICAL(KIND=8),  INTENT(IN)    :: unres_mask(points_lambda*points_phi)
                                          ! Indices to resolved pts
INTEGER(KIND=8), INTENT(OUT)   :: indices(no_point_unres)
                                          ! Radius of planet (in m)

integer(kind=8) :: err
integer(kind=4) :: ncid
integer(kind=4) :: var_output, dimid, unres_len
character(len=4) :: mpi_ext
character(len=1024) :: cmessage

write (mpi_ext, '(I0.4)') mype

call nferr(nf90_open(trim(filename) // mpi_ext, nf90_NOWRITE, ncid))

call nferr(nf90_inq_dimid(ncid, 'unres', dimid))
call nferr(nf90_inquire_dimension(ncid, dimid, len=unres_len))

if (unres_len /= no_point_unres) then
    err = 1
    write (cmessage, '(I,I,A)') unres_len, no_point_unres, "Dimension size mismatch"
    call ereport('read_spiral_circle', err, cmessage)
end if 

call nferr(nf90_inq_varid(ncid, 'output', var_output))
call nferr(nf90_get_var(ncid, var_output, indices))

call nferr(nf90_close(ncid))

end subroutine

end module
