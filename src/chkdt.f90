! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_chkdt
  use mpi
  use mod_common_mpi, only:ierr
  use mod_param     , only:nscal,is_impdiff,is_impdiff_1d,is_ifc_acdi
  use mod_types
  implicit none
  private
  public chkdt
  contains
  subroutine chkdt(n,nh,nh_g,dl,dzci,dzfi,is_solve_ns,is_track_interface,mu12,rho12,sigma,gacc,u,v,w,dtmax,gam,seps,ka12,cp12)
    !
    ! computes maximum allowed time step
    ! see Kang et al. Journal of Scientific Computing 15, 323–360
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    integer , intent(in)               :: nh,nh_g
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(1-nh_g:) :: dzci,dzfi
    logical , intent(in) :: is_solve_ns,is_track_interface
    real(rp), intent(in) :: mu12(2),rho12(2),sigma,gacc(3),gam,seps
    real(rp), intent(in), dimension(1-nh:,1-nh:,1-nh:) :: u,v,w
    real(rp), intent(out) :: dtmax
    real(rp), intent(in ), optional :: ka12(:,:),cp12(:,:) ! (2,nscal)
    real(rp) :: dxi,dyi,dzi
    real(rp) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
    real(rp) :: dtix,dtiy,dtiz,dtiv,dtig,dtik,dtipsi,dti,dl2i
    integer :: i,j,k,is
    real(rp), save :: dlmin,dzmin
    logical , save :: is_first = .true.
    !
    dxi = 1./dl(1)
    dyi = 1./dl(2)
    dzi = 1./dl(3)
    if(is_first) then ! calculate dlmin only once
      is_first = .false.
      dzmin = minval(1./dzfi)
      call MPI_ALLREDUCE(MPI_IN_PLACE,dzmin,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
      dlmin = min(minval(dl(1:2)),dzmin)
    end if
    !
    !
    ! largest eigenvalue of the discrete Laplacian, times the inverse of the stability interval of the
    ! AB2 scheme along the negative real axis (one): the diffusion terms limit the time step to
    ! 1/(diffusivity*dl2i); without the direction integrated implicitly, if any
    !
    dl2i = 4.*(dxi**2+dyi**2+dzmin**(-2))
    if(is_impdiff_1d) dl2i = 4.*(dxi**2+dyi**2)
    !
    dti = 0.
    !$acc data copy(dti) async(1)
    !$acc parallel loop collapse(3) default(present) private(ux,uy,uz,vx,vy,vz,wx,wy,wz,dtix,dtiy,dtiz) reduction(max:dti) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(ux,uy,uz,vx,vy,vz,wx,wy,wz,dtix,dtiy,dtiz) REDUCTION(max:dti)
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
    if(is_ifc_acdi) then
      dtipsi = gam*seps/dlmin**2
      if(dtipsi == 0.) dtipsi = 1.
    else
      dtipsi = dti
    end if
    if(.not.is_solve_ns) then
      dtmax = min(dti**(-1),dtipsi**(-1))
    else
      dtiv = maxval(mu12(:)/rho12(:))*dl2i
      if(is_impdiff.and..not.is_impdiff_1d) dtiv = 0.
      dtik = sqrt(sigma/(minval(rho12(:)))/dlmin**3)
      dtig = maxval(abs(gacc))/dlmin
      dti = 0.5*(dti+dtiv+sqrt((dti+dtiv)**2+4.*(dtig**2+dtik**2)))
      if(dti    == 0.) dti    = 1.
      dtmax = dti**(-1)
      if(is_track_interface.and.dtipsi > 0.) dtmax = min(dtmax,dtipsi**(-1))
    end if
    if(nscal > 0) then
      if(present(ka12) .and. present(cp12)) then
        dti = 0.
        do is=1,nscal
          dti = max(dti,maxval(ka12(:,is)/(rho12(:)*cp12(:,is)))*dl2i)
        end do
      end if
      if(is_impdiff.and..not.is_impdiff_1d) dti = 0.
      if(dti == 0.) dti = 1.
      dtmax = min(dtmax,dti**(-1))
    end if
  end subroutine chkdt
end module mod_chkdt
