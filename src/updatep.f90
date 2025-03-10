! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 The CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_updatep
  use mod_types
  implicit none
  private
  public updatep,extrapl_p
  contains
  subroutine updatep(pp,p)
    !
    ! updates the final pressure
    !
    implicit none
    real(rp), intent(in   ), dimension(0:,0:,0:) :: pp
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    !
    !$acc kernels default(present) async(1)
    p(:,:,:) = p(:,:,:) + pp(:,:,:)
    !$acc end kernels
  end subroutine updatep
  !
  subroutine extrapl_p(dt,dto,pp,po,p)
    !
    ! updates the final pressure
    !
    implicit none
    real(rp), intent(in   ) :: dt,dto
    real(rp), intent(in   ), dimension(0:,0:,0:) :: p
    real(rp), intent(inout), dimension(0:,0:,0:) :: po
    real(rp), intent(out  ), dimension(0:,0:,0:) :: pp
    real(rp) :: factor
    !
    factor = dt/dto
    !$acc kernels default(present) async(1)
    pp(:,:,:) = (1.+factor)*p(:,:,:) - factor*po(:,:,:)
    po(:,:,:) = p(:,:,:)
    !$acc end kernels
  end subroutine extrapl_p
end module mod_updatep
