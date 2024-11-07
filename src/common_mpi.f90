! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 The CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_common_mpi
  implicit none
  public
  integer :: myid,ierr
  integer :: halo(3)
  integer :: ipencil_axis
end module mod_common_mpi
