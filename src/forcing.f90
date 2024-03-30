! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_forcing
  use mod_types
  implicit none
  private
  public lscale_forcing
  contains
  subroutine lscale_forcing(ftype,lo,hi,alpha,dt,l,dl,zc,zf,u,v,w)
    use mod_param, only: pi
    implicit none
    integer, parameter :: FTYPE_ABC          = 1, &
                          FTYPE_TAYLOR_GREEN = 2
    real(rp), parameter :: a = 1._rp, b = 1._rp, c = 1._rp, k0 = 1
    integer , intent(in) :: ftype
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in) :: alpha,dt
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: zc,zf
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w
    !
    real(rp) :: f0,factor
    real(rp) :: xxc,yyc,zzc,xxf,yyf,zzf
    integer  :: i,j,k
    !
    f0 = alpha*(2.*pi/l(3))**2
    factor = dt
    !
    select case(ftype)
    case(FTYPE_ABC)
      !$acc parallel loop collapse(3) default(present), private(xxc,yyc,zzc) async(1)
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            xxc = (i-0.5)*dl(1)/l(1)
            yyc = (j-0.5)*dl(2)/l(2)
            zzc = zc(k)/l(3)
            u(i,j,k) = u(i,j,k) + factor*( f0*(a*sin(k0*zzc) + c*cos(k0*yyc)) )
            v(i,j,k) = v(i,j,k) + factor*( f0*(b*sin(k0*xxc) + a*cos(k0*zzc)) )
            w(i,j,k) = w(i,j,k) + factor*( f0*(c*sin(k0*yyc) + b*cos(k0*xxc)) )
          enddo
        enddo
      enddo
    case default ! FTYPE_TAYLOR_GREEN
      !$acc parallel loop collapse(3) default(present), private(xxc,yyc,zzc,xxf,yyf,zzf) async(1)
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            xxc = (i-0.5)*dl(1)/l(1)
            xxf = xxc + 0.5*dl(1)/l(1)
            yyc = (j-0.5)*dl(2)/l(2)
            yyf = yyc + 0.5*dl(2)/l(2)
            zzc = zc(k)/l(3)
            zzf = zf(k)/l(3)
            u(i,j,k) = u(i,j,k) + factor*( f0*sin(k0*xxf)*cos(k0*yyc)*cos(k0*zzc) )
            v(i,j,k) = v(i,j,k) - factor*( f0*cos(k0*xxc)*sin(k0*yyf)*cos(k0*zzc) )
          enddo
        enddo
      enddo
    end select
  end subroutine lscale_forcing
end module mod_forcing
