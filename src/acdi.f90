! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_acdi
  use mpi
  use mod_common_mpi, only:ierr
  use mod_types
  use mod_param , only: eps
  implicit none
  private
  public acdi_set_epsilon,acdi_set_gamma,acdi_transport_pf
  contains
  subroutine acdi_set_epsilon(dl,dzfi,seps_factor,seps)
    !
    ! computes the acdi epsilon parameter
    !
    implicit none
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: dzfi
    real(rp), intent(in)  :: seps_factor
    real(rp), intent(out) :: seps
    real(rp), save :: dlmax
    !
    dlmax  = max(maxval(dl(1:2)),maxval(1./dzfi))
    call MPI_ALLREDUCE(MPI_IN_PLACE,dlmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    seps = dlmax*seps_factor
  end subroutine acdi_set_epsilon
  !
  subroutine acdi_set_gamma(n,gam_factor,u,v,w,gam)
    !
    ! computes the acdi gamma parameter
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(in) :: gam_factor
    real(rp), intent(out) :: gam
    real(rp) :: uc,vc,wc,vel,velmax
    integer :: i,j,k
    !
    velmax = 0.
    !$acc data copy(vel) async(1)
    !$acc parallel loop collapse(3) default(present) private(uc,vc,wc,vel) reduction(max:velmax) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          uc  = 0.5*(u(i,j,k) + u(i-1,j,k))
          vc  = 0.5*(v(i,j,k) + v(i,j-1,k))
          wc  = 0.5*(w(i,j,k) + w(i,j,k-1))
          vel = uc**2 + vc**2 + wc**2
          velmax = max(velmax,vel)
        end do
      end do
    end do
    !$acc end data
    !$acc wait(1)
    call MPI_ALLREDUCE(MPI_IN_PLACE,velmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    gam = sqrt(velmax)*gam_factor
  end subroutine acdi_set_gamma
  !
  subroutine acdi_transport_pf(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,gam,seps,u,v,w,psi,dpsidt)
    !
    !
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,dzi,gam,seps
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w
    real(rp), dimension(0:,0:,0:), intent(inout) :: psi
    real(rp), dimension(:,:,:), intent(out) :: dpsidt
    integer :: i,j,k
    real(rp) :: sepsi
    real(rp) :: wghtpp,wghtpm,wghtmp,wghtmm
    real(rp) :: adv,diff,sharp
    real(rp) :: upsiip,upsiim,vpsijp,vpsijm,wpsikp,wpsikm
    real(rp) :: dpsidxp,dpsidxm,dpsidyp,dpsidym,dpsidzp,dpsidzm
    real(rp) :: phimcm,phicmm,phiccm,phicpm,phipcm,phimmc,phimcc,phimpc,phicmc,phiccc, &
                phicpc,phipmc,phipcc,phippc,phimcp,phicmp,phiccp,phicpp,phipcp
    real(rp) :: dphidx_01,dphidy_01,dphidz_01,dphidx_11,dphidy_11,dphidz_11, &
                dphidx_02,dphidy_02,dphidz_02,dphidx_12,dphidy_12,dphidz_12, &
                dphidx_03,dphidy_03,dphidz_03,dphidx_13,dphidy_13,dphidz_13
    real(rp) :: rn_01,rn_11,rn_02,rn_12,rn_03,rn_13
    !
    sepsi = 1./seps
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ! weights along the stretched grid direction
          !
          wghtpp = dzci(k)/dzfi(k+1)
          wghtpm = dzci(k)/dzfi(k)
          wghtmp = dzci(k-1)/dzfi(k)
          wghtmm = dzci(k-1)/dzfi(k-1)
          !
          ! advection term
          !
          upsiim  = 0.5*(psi(i-1,j,k)       +psi(i,j,k)       )*u(i-1,j,k)
          upsiip  = 0.5*(psi(i+1,j,k)       +psi(i,j,k)       )*u(i  ,j,k)
          vpsijm  = 0.5*(psi(i,j-1,k)       +psi(i,j,k)       )*v(i,j-1,k)
          vpsijp  = 0.5*(psi(i,j+1,k)       +psi(i,j,k)       )*v(i,j  ,k)
          wpsikm  = 0.5*(psi(i,j,k-1)*wghtmp+psi(i,j,k)*wghtmm)*w(i,j,k-1)
          wpsikp  = 0.5*(psi(i,j,k+1)*wghtpm+psi(i,j,k)*wghtpp)*w(i,j,k  )
          adv     = dxi*(upsiip-upsiim) + dyi*(vpsijp-vpsijm) + dzfi(k)*(wpsikp-wpsikm)
          !
          ! diffusion term
          !
          dpsidxp = (psi(i+1,j,k)-psi(i  ,j,k))*dxi
          dpsidxm = (psi(i  ,j,k)-psi(i-1,j,k))*dxi
          dpsidyp = (psi(i,j+1,k)-psi(i,j  ,k))*dyi
          dpsidym = (psi(i,j  ,k)-psi(i,j-1,k))*dyi
          dpsidzp = (psi(i,j,k+1)-psi(i,j,k  ))*dzci(k  )
          dpsidzm = (psi(i,j,k  )-psi(i,j,k-1))*dzci(k-1)
          diff    = gam*seps*( (dpsidxp-dpsidxm)*dxi + (dpsidyp-dpsidym)*dyi + (dpsidzp-dpsidzm)*dzfi(k) )
          !
          ! sharpening term
          !
          !phimmm = seps*log((psi(i-1,j-1,k-1)+eps)/(1.-psi(i-1,j-1,k-1)+eps))
          phimcm = seps*log((psi(i-1,j  ,k-1)+eps)/(1.-psi(i-1,j  ,k-1)+eps))
          !phimpm = seps*log((psi(i-1,j+1,k-1)+eps)/(1.-psi(i-1,j+1,k-1)+eps))
          phicmm = seps*log((psi(i  ,j-1,k-1)+eps)/(1.-psi(i  ,j-1,k-1)+eps))
          phiccm = seps*log((psi(i  ,j  ,k-1)+eps)/(1.-psi(i  ,j  ,k-1)+eps))
          phicpm = seps*log((psi(i  ,j+1,k-1)+eps)/(1.-psi(i  ,j+1,k-1)+eps))
          !phipmm = seps*log((psi(i+1,j-1,k-1)+eps)/(1.-psi(i+1,j-1,k-1)+eps))
          phipcm = seps*log((psi(i+1,j  ,k-1)+eps)/(1.-psi(i+1,j  ,k-1)+eps))
          !phippm = seps*log((psi(i+1,j+1,k-1)+eps)/(1.-psi(i+1,j+1,k-1)+eps))
          phimmc = seps*log((psi(i-1,j-1,k  )+eps)/(1.-psi(i-1,j-1,k  )+eps))
          phimcc = seps*log((psi(i-1,j  ,k  )+eps)/(1.-psi(i-1,j  ,k  )+eps))
          phimpc = seps*log((psi(i-1,j+1,k  )+eps)/(1.-psi(i-1,j+1,k  )+eps))
          phicmc = seps*log((psi(i  ,j-1,k  )+eps)/(1.-psi(i  ,j-1,k  )+eps))
          phiccc = seps*log((psi(i  ,j  ,k  )+eps)/(1.-psi(i  ,j  ,k  )+eps))
          phicpc = seps*log((psi(i  ,j+1,k  )+eps)/(1.-psi(i  ,j+1,k  )+eps))
          phipmc = seps*log((psi(i+1,j-1,k  )+eps)/(1.-psi(i+1,j-1,k  )+eps))
          phipcc = seps*log((psi(i+1,j  ,k  )+eps)/(1.-psi(i+1,j  ,k  )+eps))
          phippc = seps*log((psi(i+1,j+1,k  )+eps)/(1.-psi(i+1,j+1,k  )+eps))
          !phimmp = seps*log((psi(i-1,j-1,k+1)+eps)/(1.-psi(i-1,j-1,k+1)+eps))
          phimcp = seps*log((psi(i-1,j  ,k+1)+eps)/(1.-psi(i-1,j  ,k+1)+eps))
          !phimpp = seps*log((psi(i-1,j+1,k+1)+eps)/(1.-psi(i-1,j+1,k+1)+eps))
          phicmp = seps*log((psi(i  ,j-1,k+1)+eps)/(1.-psi(i  ,j-1,k+1)+eps))
          phiccp = seps*log((psi(i  ,j  ,k+1)+eps)/(1.-psi(i  ,j  ,k+1)+eps))
          phicpp = seps*log((psi(i  ,j+1,k+1)+eps)/(1.-psi(i  ,j+1,k+1)+eps))
          !phipmp = seps*log((psi(i+1,j-1,k+1)+eps)/(1.-psi(i+1,j-1,k+1)+eps))
          phipcp = seps*log((psi(i+1,j  ,k+1)+eps)/(1.-psi(i+1,j  ,k+1)+eps))
          !phippp = seps*log((psi(i+1,j+1,k+1)+eps)/(1.-psi(i+1,j+1,k+1)+eps))
          !
          dphidx_01 =      (phiccc       -phimcc       )*dxi
          dphidy_01 = 0.25*(phicpc+phimpc-phicmc-phimmc)*dyi
          dphidz_01 = 0.50*(phiccp+phimcp-phiccm-phimcm)*(dzci(k-1)*dzci(k))/(dzci(k-1)+dzci(k))
          dphidx_11 =      (phipcc       -phiccc       )*dxi
          dphidy_11 = 0.25*(phippc+phicpc-phipmc-phicmc)*dyi
          dphidz_11 = 0.50*(phipcp+phiccp-phipcm-phiccm)*(dzci(k-1)*dzci(k))/(dzci(k-1)+dzci(k))
          dphidx_02 = 0.25*(phipcc+phipmc-phimcc-phimmc)*dxi
          dphidy_02 =      (phiccc       -phicmc       )*dyi
          dphidz_02 = 0.50*(phiccp+phicmp-phiccm-phicmm)*(dzci(k-1)*dzci(k))/(dzci(k-1)+dzci(k))
          dphidx_12 = 0.25*(phippc+phipcc-phimpc-phimcc)*dxi
          dphidy_12 =      (phicpc       -phiccc       )*dyi
          dphidz_12 = 0.50*(phicpp+phiccp-phicpm-phiccm)*(dzci(k-1)*dzci(k))/(dzci(k-1)+dzci(k))
          dphidx_03 = 0.25*((phipcm-phimcm)*wghtmp+(phipcc-phimcc)*wghtmm)*dxi
          dphidy_03 = 0.25*((phicpm-phicmm)*wghtmp+(phicpc-phicmc)*wghtmm)*dyi
          dphidz_03 =      ( phiccc               -phiccm                )*dzci(k-1)
          dphidx_13 = 0.25*((phipcp-phimcp)*wghtpm+(phipcc-phimcc)*wghtpp)*dxi
          dphidy_13 = 0.25*((phicpp-phicmp)*wghtpm+(phicpc-phicmc)*wghtpp)*dyi
          dphidz_13 =      ( phiccp               -phiccc                )*dzci(k)
          !
          rn_01 = dphidx_01/(sqrt(dphidx_01**2+dphidy_01**2+dphidz_01**2)+eps)
          rn_11 = dphidx_11/(sqrt(dphidx_11**2+dphidy_11**2+dphidz_11**2)+eps)
          rn_02 = dphidy_02/(sqrt(dphidx_02**2+dphidy_02**2+dphidz_02**2)+eps)
          rn_12 = dphidy_12/(sqrt(dphidx_12**2+dphidy_12**2+dphidz_12**2)+eps)
          rn_03 = dphidz_03/(sqrt(dphidx_03**2+dphidy_03**2+dphidz_03**2)+eps)
          rn_13 = dphidz_13/(sqrt(dphidx_13**2+dphidy_13**2+dphidz_13**2)+eps)
          !
          sharp = 0.25*gam*( (1.-(tanh(0.25*(phipcc       +phiccc       )*sepsi))**2)*rn_11 - &
                             (1.-(tanh(0.25*(phiccc       +phimcc       )*sepsi))**2)*rn_01 )*dxi + &
                  0.25*gam*( (1.-(tanh(0.25*(phicpc       +phiccc       )*sepsi))**2)*rn_12 - &
                             (1.-(tanh(0.25*(phiccc       +phicmc       )*sepsi))**2)*rn_02 )*dyi + &
                  0.25*gam*( (1.-(tanh(0.25*(phiccp*wghtpm+phiccc*wghtpp)*sepsi))**2)*rn_13 - &
                             (1.-(tanh(0.25*(phiccc*wghtmm+phiccm*wghtmp)*sepsi))**2)*rn_03 )*dzfi(k)
          !
          ! full transport
          !
          dpsidt(i,j,k) = -adv+diff-sharp
        end do
      end do
    end do
  end subroutine acdi_transport_pf
  !
  pure elemental real(rp) function acdi_cmpt_phi(psi,seps) result(phi)
    use mod_param, only:eps
    !$acc routine seq
    !
    ! computes the signed-distance field phi from the smoothed volume fraction field psi
    !
    implicit none
    !
    real(rp), intent(in) :: psi,seps
    phi = seps*log((psi+eps)/(1.-psi+eps))
  end function acdi_cmpt_phi
end module mod_acdi
