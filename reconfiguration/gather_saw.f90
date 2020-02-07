!> 
!! Copyright 2020 Scott Wales
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

module gather_saw_mod
    private
    public gather_field_saw
contains

    subroutine gather_field_saw(l_field, g_field, &
        l_nx, l_ny, g_nx, g_ny, &
        field_type, halo_type, root, group)

        use mpi_f08

        integer, parameter :: max_ranks = 1024

        logical(kind=8), intent(in) :: l_field(:) ! Local data
        logical(kind=8), intent(out) :: g_field(:) ! Global data

        integer(kind=8), intent(in) :: l_nx, l_ny ! Local size
        integer(kind=8), intent(in) :: g_nx, g_ny ! Local size

        integer(kind=8), intent(in) :: root ! Gather process
        integer(kind=8), intent(in) :: field_type, halo_type, group ! Ignored

        integer(kind=4) comm_rank
        integer :: handle

        ! Get MPI info
        call MPI_comm_rank(MPI_COMM_WORLD, comm_rank)

        if (comm_rank == dest) then
            open(newunit=handle, file='/g/data/w35/saw562/HighResLAM/mask.bin', FORM='BINARY')
            read(handle) g_field
            close(handle)
        end if

    end subroutine

end module
