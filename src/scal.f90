! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 The CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_scal
  use mod_types
  implicit none
  private
  public scal_ad,cmpt_scalflux
  contains
  subroutine scal_ad(nx,ny,nz,dxi,dyi,dzci,dzfi,ka12,rhocp12,psiflx_x,psiflx_y,psiflx_z,psi,u,v,w,s,dsdt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), intent(in), dimension(2) :: ka12,rhocp12
    real(rp), intent(in), dimension(0:,0:,0:), optional :: psiflx_x,psiflx_y,psiflx_z
    real(rp), intent(in), dimension(0:,0:,0:) :: psi,u,v,w,s
    real(rp), dimension(:,:,:), intent(out) :: dsdt
    integer :: i,j,k
    real(rp) :: usip,usim,vsjp,vsjm,wskp,wskm
    real(rp) :: dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
    real(rp) :: kaxp,kaxm,kayp,kaym,kazp,kazm,rhocp_c
    real(rp) :: ka,dka,rhocp,drhocp
    !
    ka    = ka12(2)   ; dka    = ka12(1)-ka12(2)
    rhocp = rhocp12(2); drhocp = rhocp12(1)-rhocp12(2)
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
#if defined(_CONSISTENT_ADVECTION)
          usip = 0.5*(rhocp*u(i  ,j,k)+drhocp*psiflx_x(i  ,j,k))*(s(i+1,j,k)+s(i,j,k))
          usim = 0.5*(rhocp*u(i-1,j,k)+drhocp*psiflx_x(i-1,j,k))*(s(i-1,j,k)+s(i,j,k))
          vsjm = 0.5*(rhocp*v(i,j  ,k)+drhocp*psiflx_y(i,j  ,k))*(s(i,j-1,k)+s(i,j,k))
          vsjp = 0.5*(rhocp*v(i,j-1,k)+drhocp*psiflx_y(i,j-1,k))*(s(i,j+1,k)+s(i,j,k))
          wskp = 0.5*(rhocp*w(i,j,k  )+drhocp*psiflx_z(i,j,k  ))*(s(i,j,k+1)+s(i,j,k))
          wskm = 0.5*(rhocp*w(i,j,k-1)+drhocp*psiflx_z(i,j,k-1))*(s(i,j,k-1)+s(i,j,k))
#else
          usip = 0.5*u(i  ,j,k)*(s(i+1,j,k)+s(i,j,k))
          usim = 0.5*u(i-1,j,k)*(s(i-1,j,k)+s(i,j,k))
          vsjm = 0.5*v(i,j-1,k)*(s(i,j-1,k)+s(i,j,k))
          vsjp = 0.5*v(i,j  ,k)*(s(i,j+1,k)+s(i,j,k))
          wskp = 0.5*w(i,j,k  )*(s(i,j,k+1)+s(i,j,k))
          wskm = 0.5*w(i,j,k-1)*(s(i,j,k-1)+s(i,j,k))
#endif
          !
          kaxp = ka+dka*0.5*(psi(i+1,j,k)+psi(i  ,j,k))
          kaxm = ka+dka*0.5*(psi(i  ,j,k)+psi(i-1,j,k))
          kayp = ka+dka*0.5*(psi(i,j+1,k)+psi(i,j  ,k))
          kaym = ka+dka*0.5*(psi(i,j  ,k)+psi(i,j-1,k))
          kazp = ka+dka*0.5*(psi(i,j,k+1)+psi(i,j,k  ))
          kazm = ka+dka*0.5*(psi(i,j,k  )+psi(i,j,k-1))
          !
          dsdxp = (s(i+1,j,k)-s(i  ,j,k))*dxi
          dsdxm = (s(i  ,j,k)-s(i-1,j,k))*dxi
          dsdyp = (s(i,j+1,k)-s(i,j  ,k))*dyi
          dsdym = (s(i,j  ,k)-s(i,j-1,k))*dyi
          dsdzp = (s(i,j,k+1)-s(i,j,k  ))*dzci(k  )
          dsdzm = (s(i,j,k  )-s(i,j,k-1))*dzci(k-1)
          !
#if defined(_CONSISTENT_ADVECTION)
          dsdt(i,j,k) = dxi*(     -usip + usim ) + &
                        dyi*(     -vsjp + vsjm ) + &
                        dzfi(k)*( -wskp + wskm ) + &
#else
          rhocp_c = rhocp+drhocp*psi(i,j,k)
          dsdt(i,j,k) = ( dxi*(     -usip + usim ) + &
                          dyi*(     -vsjp + vsjm ) + &
                          dzfi(k)*( -wskp + wskm ) )*rhocp_c + &
#endif
                        ( (kaxp*dsdxp-kaxm*dsdxm)*dxi + &
                          (kayp*dsdyp-kaym*dsdym)*dyi + &
                          (kazp*dsdzp-kazm*dsdzm)*dzfi(k) )
        end do
      end do
    end do
  end subroutine scal_ad
  !
  subroutine cmpt_scalflux(n,is_bound,l,dli,dzci,dzfi,ka12,psi,s,flux)
    use mpi
    implicit none
    integer , intent(in ), dimension(3) :: n
    logical , intent(in ), dimension(0:1,3) :: is_bound
    real(rp), intent(in ), dimension(3)     :: l,dli
    real(rp), intent(in ), dimension(0:)    :: dzci,dzfi
    real(rp), intent(in ), dimension(2) :: ka12
    real(rp), intent(in ), dimension(0:,0:,0:) :: psi,s
    real(rp), intent(out), dimension(3) :: flux
    real(rp) :: dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
    real(rp) :: kaxp,kaxm,kayp,kaym,kazp,kazm
    real(rp) :: flux_x,flux_y,flux_z
    integer :: i,j,k,nx,ny,nz
    real(rp) :: dxi,dyi,lx,ly,lz
    real(rp) :: ka,dka
    integer :: ierr
    !
    ka = ka12(2); dka = ka12(1)-ka12(2)
    !
    nx = n(1); ny = n(2); nz = n(3)
    dxi = dli(1); dyi = dli(2)
    lx = l(1); ly = l(2); lz = l(3)
    flux_x = 0._rp
    !$acc data copy(flux_x) async(1)
    if(is_bound(0,1)) then
      i = 0
      !$acc parallel loop collapse(2) default(present) private(dsdxp,kaxp) reduction(+:flux_x) async(1)
      do k=1,nz
        do j=1,ny
          dsdxp = (s(i+1,j,k)-s(i  ,j,k))*dxi
          kaxp = ka+dka*0.5*(psi(i+1,j,k)+psi(i,j,k))
          flux_x = flux_x + kaxp*dsdxp/(dyi*dzfi(k)*ly*lz)
        end do
      end do
    end if
    if(is_bound(1,1)) then
      i = nx+1
      !$acc parallel loop collapse(2) default(present) private(dsdxm,kaxm) reduction(+:flux_x) async(1)
      do k=1,nz
        do j=1,ny
          dsdxm = (s(i  ,j,k)-s(i-1,j,k))*dxi
          kaxm = ka+dka*0.5*(psi(i,j,k)+psi(i-1,j,k))
          flux_x = flux_x - kaxm*dsdxm/(dyi*dzfi(k)*ly*lz)
        end do
      end do
    end if
    !$acc end data
    flux_y = 0._rp
    !$acc data copy(flux_y) async(1)
    if(is_bound(0,2)) then
      j = 0
      !$acc parallel loop collapse(2) default(present) private(dsdyp,kayp) reduction(+:flux_y) async(1)
      do k=1,nz
        do i=1,nx
          dsdyp = (s(i,j+1,k)-s(i,j  ,k))*dyi
          kayp = ka+dka*0.5*(psi(i,j+1,k)+psi(i,j,k))
          flux_y = flux_y + kayp*dsdyp/(dxi*dzfi(k)*lx*lz)
        end do
      end do
    end if
    if(is_bound(1,2)) then
      j = ny+1
      !$acc parallel loop collapse(2) default(present) private(dsdym,kaym) reduction(+:flux_y) async(1)
      do k=1,nz
        do i=1,nx
          dsdym = (s(i,j  ,k)-s(i,j-1,k))*dyi
          kaym = ka+dka*0.5*(psi(i,j,k)+psi(i,j-1,k))
          flux_y = flux_y - kaym*dsdym/(dxi*dzfi(k)*lx*lz)
        end do
      end do
    end if
    !$acc end data
    flux_z = 0._rp
    !$acc data copy(flux_z) async(1)
    if(is_bound(0,3)) then
      k = 0
      !$acc parallel loop collapse(2) default(present) private(dsdzp,kazp) reduction(+:flux_z) async(1)
      do j=1,ny
        do i=1,nx
          dsdzp = (s(i,j,k+1)-s(i,j,k  ))*dzci(k  )
          kazp = ka+dka*0.5*(psi(i,j,k+1)+psi(i,j,k))
          flux_z = flux_z + kazp*dsdzp/(dxi*dyi*lx*ly)
        end do
      end do
    end if
    if(is_bound(1,3)) then
      k = nz+1
      !$acc parallel loop collapse(2) default(present) private(dsdzm,kazm) reduction(+:flux_z) async(1)
      do j=1,ny
        do i=1,nx
          dsdzm = (s(i,j,k  )-s(i,j,k-1))*dzci(k-1)
          kazm = ka+dka*0.5*(psi(i,j,k)+psi(i,j,k-1))
          flux_z = flux_z - kazm*dsdzm/(dxi*dyi*lx*ly)
        end do
      end do
    end if
    !$acc end data
    !$acc wait(1)
    flux(:) = [flux_x,flux_y,flux_z]
    call MPI_ALLREDUCE(MPI_IN_PLACE,flux(3),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  end subroutine cmpt_scalflux
end module mod_scal
