! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_forcing
  use mod_param, only: pi
  use mod_types
  implicit none
  private
  public lscale_forcing
  contains
  subroutine lscale_forcing(ftype,lo,hi,nh_v,nh_g,alpha,dt,l,dl,zc,u,v,w)
    !
    ! adds a steady large-scale acceleration to the velocity field of a triply-periodic domain, e.g., to
    ! sustain homogeneous isotropic turbulence; with k = 2*pi*k0/l along each direction:
    !
    !   'tgv', Taylor-Green              : f0*( sin(k*x)*cos(k*y)*cos(k*z),
    !                                          -cos(k*x)*sin(k*y)*cos(k*z),
    !                                           0                          )
    !   'abc', Arnold-Beltrami-Childress : f0*( a*sin(k*z) + c*cos(k*y),
    !                                           b*sin(k*x) + a*cos(k*z),
    !                                           c*sin(k*y) + b*cos(k*x)    )
    !
    ! with f0 = alpha*(2*pi/l(3))**2, and alpha the input `forced_hit_amp`
    !
    implicit none
    real(rp), parameter :: a = 1._rp, b = 1._rp, c = 1._rp, k0 = 2._rp
    character(len=*), intent(in   ) :: ftype
    integer , intent(in   ), dimension(3) :: lo,hi
    integer , intent(in   ) :: nh_v,nh_g
    real(rp), intent(in   ) :: alpha,dt
    real(rp), intent(in   ), dimension(3) :: l,dl
    real(rp), intent(in   ), dimension(lo(3)-nh_g:) :: zc
    real(rp), intent(inout), dimension(lo(1)-nh_v:,lo(2)-nh_v:,lo(3)-nh_v:) :: u,v,w
    !
    real(rp) :: f0,factor
    real(rp) :: xxc,yyc,zzc,xxf,yyf
    integer  :: i,j,k
    !
    f0 = alpha*(2.*pi/l(3))**2
    factor = dt
    !
    select case(ftype)
    case('abc')
      !$acc parallel     loop collapse(3) default(present) private(xxc,yyc,zzc) async(1)
      !$omp target teams loop collapse(3)                  private(xxc,yyc,zzc)
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            xxc = 2.*pi*k0*(i-0.5)*dl(1)/l(1)
            yyc = 2.*pi*k0*(j-0.5)*dl(2)/l(2)
            zzc = 2.*pi*k0*zc(k)/l(3)
            u(i,j,k) = u(i,j,k) + factor*( f0*(a*sin(zzc) + c*cos(yyc)) )
            v(i,j,k) = v(i,j,k) + factor*( f0*(b*sin(xxc) + a*cos(zzc)) )
            w(i,j,k) = w(i,j,k) + factor*( f0*(c*sin(yyc) + b*cos(xxc)) )
          end do
        end do
      end do
    case default ! 'tgv'
      !$acc parallel     loop collapse(3) default(present) private(xxc,yyc,zzc,xxf,yyf) async(1)
      !$omp target teams loop collapse(3)                  private(xxc,yyc,zzc,xxf,yyf)
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            xxc = 2.*pi*k0*(i-0.5)*dl(1)/l(1)
            xxf = 2.*pi*k0*(i    )*dl(1)/l(1)
            yyc = 2.*pi*k0*(j-0.5)*dl(2)/l(2)
            yyf = 2.*pi*k0*(j    )*dl(2)/l(2)
            zzc = 2.*pi*k0*zc(k)/l(3)
            u(i,j,k) = u(i,j,k) + factor*( f0*sin(xxf)*cos(yyc)*cos(zzc) )
            v(i,j,k) = v(i,j,k) - factor*( f0*cos(xxc)*sin(yyf)*cos(zzc) )
          end do
        end do
      end do
    end select
  end subroutine lscale_forcing
end module mod_forcing
