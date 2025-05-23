! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 The CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_chkdt
  use mpi
  use mod_common_mpi, only:ierr
  use mod_types
  implicit none
  private
  public chkdt
  contains
  subroutine chkdt(n,dl,dzci,dzfi,is_solve_ns,is_track_interface,mu12,rho12,sigma,gacc,u,v,w,dtmax,gam,seps,ka12,cp12)
    !
    ! computes maximum allowed time step
    ! see Kang et al. Journal of Scientific Computing 15, 323–360
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    logical,  intent(in) :: is_solve_ns,is_track_interface
    real(rp), intent(in) :: mu12(2),rho12(2),sigma,gacc(3),gam,seps
    real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out) :: dtmax
    real(rp), intent(in ), optional :: ka12(2),cp12(2)
    real(rp) :: dxi,dyi,dzi
    real(rp) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
    real(rp) :: dtix,dtiy,dtiz,dtiv,dtig,dtik,dtipsi,dti
    integer :: i,j,k
    real(rp), save :: dlmin
    logical , save :: is_first = .true.
    !
    dxi = 1./dl(1)
    dyi = 1./dl(2)
    dzi = 1./dl(3)
    if(is_first) then ! calculate dlmin only once
      is_first = .false.
      dlmin     = minval(dl(1:2))
      dlmin     = min(dlmin,minval(1./dzfi))
      call MPI_ALLREDUCE(MPI_IN_PLACE,dlmin,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
    end if
    !
    dti = 0.
    !$acc data copy(dti) async(1)
    !$acc parallel loop collapse(3) default(present) private(ux,uy,uz,vx,vy,vz,wx,wy,wz,dtix,dtiy,dtiz) reduction(max:dti) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          ux = abs(u(i,j,k))
          vx = 0.25*abs( v(i,j,k)+v(i,j-1,k)+v(i+1,j,k)+v(i+1,j-1,k) )
          wx = 0.25*abs( w(i,j,k)+w(i,j,k-1)+w(i+1,j,k)+w(i+1,j,k-1) )
          dtix = ux*dxi+vx*dyi+wx*dzfi(k)
          uy = 0.25*abs( u(i,j,k)+u(i,j+1,k)+u(i-1,j+1,k)+u(i-1,j,k) )
          vy = abs(v(i,j,k))
          wy = 0.25*abs( w(i,j,k)+w(i,j+1,k)+w(i,j+1,k-1)+w(i,j,k-1) )
          dtiy = uy*dxi+vy*dyi+wy*dzfi(k)
          uz = 0.25*abs( u(i,j,k)+u(i-1,j,k)+u(i-1,j,k+1)+u(i,j,k+1) )
          vz = 0.25*abs( v(i,j,k)+v(i,j-1,k)+v(i,j-1,k+1)+v(i,j,k+1) )
          wz = abs(w(i,j,k))
          dtiz = uz*dxi+vz*dyi+wz*dzci(k)
          dti = max(dti,dtix,dtiy,dtiz)
        end do
      end do
    end do
    !$acc end data
    !$acc wait(1)
    call MPI_ALLREDUCE(MPI_IN_PLACE,dti,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
#if !defined(_INTERFACE_CAPTURING_VOF)
    dtipsi = gam*seps/dlmin**2
    if(dtipsi == 0.) dtipsi = 1.
#else
    dtipsi = dti
#endif
    if(.not.is_solve_ns) then
      dtmax = min(dti**(-1),dtipsi**(-1))
    else
      dtiv = maxval(mu12(:)/rho12(:))*2.*(3./dlmin**2)
      dtik = sqrt(sigma/(minval(rho12(:)))/dlmin**3)
      dtig = maxval(abs(gacc))/dlmin
      dti = 0.5*(dti+dtiv+sqrt((dti+dtiv)**2+4.*(dtig**2+dtik**2)))
      if(dti    == 0.) dti    = 1.
      dtmax = dti**(-1)
      if(is_track_interface) dtmax = min(dtmax,dtipsi**(-1))
    end if
#if defined(_SCALAR)
    if(present(ka12) .and. present(cp12)) &
      dti = maxval(ka12(:)/(rho12(:)*cp12(:)))/dlmin**2
    if(dti == 0.) dti = 1.
    dtmax = min(dtmax,dti**(-1))
#endif
  end subroutine chkdt
end module mod_chkdt
