! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_correc
  use mod_types
  implicit none
  private
  public correc
  contains
  subroutine correc(n,nh_v,nh_p,nh_psi,nh_g,dli,dzci,wzc,rho12,dt,p,psi,u,v,w)
    !
    ! corrects the velocity so that it is divergence free;
    ! with the variable-coefficient Poisson solver, the pressure gradient is divided by the density at the
    ! cell faces (from `rho12` and `psi`, with the volume weights `wzc` along z), while with the
    ! constant-coefficient one the reference density is part of `dt`, and `nh_psi`, `wzc`, `rho12` and `psi`
    ! are not used
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    integer , intent(in)               :: nh_v,nh_p,nh_psi,nh_g
    real(rp), intent(in), dimension(3 ) :: dli
    real(rp), intent(in), dimension(1-nh_g:) :: dzci,wzc
    real(rp), intent(in), dimension(2) :: rho12
    real(rp), intent(in) :: dt
    real(rp), intent(in   ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
#if !defined(_CONSTANT_COEFFS_POISSON)
    real(rp) :: rho,drho
#endif
    real(rp) :: factori,factorj
    !real(rp), dimension(0:n(3)+1) :: factork
    integer :: i,j,k,m
    !
    !factor = rkcoeffab(rkiter)*dt
    !
    ! the velocity is also corrected in the halo layers covered by `p`, since the
    ! predicted velocity halos are updated before; only the upper halo layers along
    ! the normal direction are left for the halo exchange that follows
    !
    m = min(nh_v,nh_p)
#if !defined(_CONSTANT_COEFFS_POISSON)
    m = min(m,nh_psi)
    rho = rho12(2); drho = rho12(1)-rho12(2)
#endif
    !
    factori = dt*dli(1)
    factorj = dt*dli(2)
    !factork = dt*dzci!dli(3)
    !$acc parallel     loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3)
    do k=1-m,n(3)+m
      do j=1-m,n(2)+m
        do i=1-m,n(1)
#if defined(_CONSTANT_COEFFS_POISSON)
          u(i,j,k) = u(i,j,k) - factori*(   p(i+1,j,k)-p(i,j,k))
#else
          u(i,j,k) = u(i,j,k) - factori*(   p(i+1,j,k)-p(i,j,k))/(rho + drho*0.5*(psi(i+1,j,k)+psi(i,j,k)))
#endif
        end do
      end do
    end do
    !$acc parallel     loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3)
    do k=1-m,n(3)+m
      do j=1-m,n(2)
        do i=1-m,n(1)+m
#if defined(_CONSTANT_COEFFS_POISSON)
          v(i,j,k) = v(i,j,k) - factorj*(   p(i,j+1,k)-p(i,j,k))
#else
          v(i,j,k) = v(i,j,k) - factorj*(   p(i,j+1,k)-p(i,j,k))/(rho + drho*0.5*(psi(i,j+1,k)+psi(i,j,k)))
#endif
        end do
      end do
    end do
    !$acc parallel     loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3)
    do k=1-m,n(3)
      do j=1-m,n(2)+m
        do i=1-m,n(1)+m
#if defined(_CONSTANT_COEFFS_POISSON)
          w(i,j,k) = w(i,j,k) - dt*dzci(k)*(p(i,j,k+1)-p(i,j,k))
#else
          w(i,j,k) = w(i,j,k) - dt*dzci(k)*(p(i,j,k+1)-p(i,j,k))/(rho + drho*(wzc(k)*psi(i,j,k)+(1.-wzc(k))*psi(i,j,k+1)))
#endif
        end do
      end do
    end do
  end subroutine correc
end module mod_correc
