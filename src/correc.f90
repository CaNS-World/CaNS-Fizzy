! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_correc
  use mod_types
  implicit none
  private
  public correc
  contains
  subroutine correc(n,dli,dzci,rho0,rho12,dt,p,psi,u,v,w)
    !
    ! corrects the velocity so that it is divergence free
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3 ) :: dli
    real(rp), intent(in), dimension(0:) :: dzci
    real(rp), intent(in)                :: rho0,rho12(2)
    real(rp), intent(in) :: dt
    real(rp), intent(in   ), dimension(0:,0:,0:) :: p,psi
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    real(rp) :: factori,factorj
    integer :: i,j,k
    real(rp) :: rhop,rho,drho
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    rhop = 1.
    !
    !factor = rkcoeffab(rkiter)*dt
    !
    factori = dt*dli(1)
    factorj = dt*dli(2)
    !factork = dt*dzci!dli(3)
    !$acc parallel loop collapse(3) default(present) firstprivate(rhop) async(1)
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)
#if !defined(_CONSTANT_COEFFS_POISSON)
          rhop = rho + drho*0.5*(psi(i+1,j,k)+psi(i,j,k))
#endif
          u(i,j,k) = u(i,j,k) - factori*(p(i+1,j,k)-p( i,j,k))/rhop
        end do
      end do
    end do
    !$acc parallel loop collapse(3) default(present) firstprivate(rhop) async(1)
    do k=0,n(3)+1
      do j=0,n(2)
        do i=0,n(1)+1
#if !defined(_CONSTANT_COEFFS_POISSON)
          rhop = rho + drho*0.5*(psi(i,j+1,k)+psi(i,j,k))
#endif
          v(i,j,k) = v(i,j,k) - factorj*(p(i,j+1,k)-p(i,j,k))/rhop
        end do
      end do
    end do
    !$acc parallel loop collapse(3) default(present) firstprivate(rhop) async(1)
    do k=0,n(3)
      do j=0,n(2)+1
        do i=0,n(1)+1
#if !defined(_CONSTANT_COEFFS_POISSON)
          rhop = rho + drho*0.5*(psi(i,j,k+1)+psi(i,j,k))
#endif
          w(i,j,k) = w(i,j,k) - dt*dzci(k)*(p(i,j,k+1)-p(i,j,k))/rhop
        end do
      end do
    end do
  end subroutine correc
end module mod_correc
