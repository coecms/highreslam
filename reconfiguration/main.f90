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

program main
    use mpi_f08
    use offline_spiral_circle
    
    character(len=:), allocatable :: path
    integer :: path_len
    type(spiral_circle_data) :: sc
    integer(kind=8), allocatable :: output_idx(:)

    call mpi_init

    call get_command_argument(1, length=path_len)
    allocate(character(len=path_len) :: path)
    call get_command_argument(1, value=path)

    call read_netcdf(path, sc)
    call spiral_circle(sc, output_idx)
    call write_netcdf(path, output_idx)

    call mpi_finalize

end program
