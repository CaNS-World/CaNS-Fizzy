! -
!
! SPDX-FileCopy_11Text: Copy_11 (c) 2017-2022 Pedro Costa and the CaNS contributors. All _11s reserved.
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
    !$acc data copy(vel) async(1)
    !$acc parallel loop collapse(3) default(present) private(uc,vc,wc) reduction(max:vel) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(uc,vc,wc) REDUCTION(max:vel)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          uc  = 0.5*(u(i,j,k) + u(i-1,j,k))
          vc  = 0.5*(v(i,j,k) + v(i,j-1,k))
          wc  = 0.5*(w(i,j,k) + w(i,j,k-1))
          vel = sqrt(uc**2 + vc**2 + wc**2)
        end do
      end do
    end do
    !$acc end data
    !$acc wait(1)
    call MPI_ALLREDUCE(MPI_IN_PLACE,vel,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    velmax = vel
    gam = velmax*gam_factor
  end subroutine acdi_set_gamma
  !
  subroutine acdi_transport_pf(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,gam,seps,u,v,w,phi,dphidt)
    !
    !
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,dzi,gam,seps
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w
    real(rp), dimension(0:,0:,0:), intent(inout) :: phi
    real(rp), dimension(:,:,:), intent(out) :: dphidt
    integer :: i,j,k
    real(rp) :: sepsi
    real(rp) :: wghtpp,wghtpm,wghtmp,wghtmm
    real(rp) :: ADV,DIFF,SHARP
    real(rp) :: uphiip,uphiim,vphijp,vphijm,wphikp,wphikm
    real(rp) :: dphidxp,dphidxm,dphidyp,dphidym,dphidzp,dphidzm
    real(rp) :: psimcm,psicmm,psiccm,psicpm,psipcm,psimmc,psimcc,psimpc,psicmc,psiccc, &
                psicpc,psipmc,psipcc,psippc,psimcp,psicmp,psiccp,psicpp,psipcp
    real(rp) :: dpsidx_01,dpsidy_01,dpsidz_01,dpsidx_11,dpsidy_11,dpsidz_11, &
                dpsidx_02,dpsidy_02,dpsidz_02,dpsidx_12,dpsidy_12,dpsidz_12, &
                dpsidx_03,dpsidy_03,dpsidz_03,dpsidx_13,dpsidy_13,dpsidz_13
    real(rp) :: rn_01,rn_11,rn_02,rn_12,rn_03,rn_13

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
          uphiim  = 0.5*(phi(i-1,j,k)       +phi(i,j,k)       )*u(i-1,j,k)
          uphiip  = 0.5*(phi(i+1,j,k)       +phi(i,j,k)       )*u(i  ,j,k)
          vphijm  = 0.5*(phi(i,j-1,k)       +phi(i,j,k)       )*v(i,j-1,k)
          vphijp  = 0.5*(phi(i,j+1,k)       +phi(i,j,k)       )*v(i,j  ,k)
          wphikm  = 0.5*(phi(i,j,k-1)*wghtmp+phi(i,j,k)*wghtmm)*w(i,j,k-1)
          wphikp  = 0.5*(phi(i,j,k+1)*wghtpm+phi(i,j,k)*wghtpp)*w(i,j,k  )
          adv     = dxi*(uphiip-uphiim) + dyi*(vphijp-vphijm) + dzfi(k)*(wphikp-wphikm)
          !
          ! diffusion term
          !
          dphidxp = (phi(i+1,j,k)-phi(i  ,j,k))*dxi
          dphidxm = (phi(i  ,j,k)-phi(i-1,j,k))*dxi
          dphidyp = (phi(i,j+1,k)-phi(i,j  ,k))*dyi
          dphidym = (phi(i,j  ,k)-phi(i,j-1,k))*dyi
          dphidzp = (phi(i,j,k+1)-phi(i,j,k  ))*dzci(k  )
          dphidzm = (phi(i,j,k  )-phi(i,j,k-1))*dzci(k-1)
          diff    = gam*seps*( (dphidxp-dphidxm)*dxi + (dphidyp-dphidym)*dyi + (dphidzp-dphidzm)*dzfi(k) )
          !
          ! sharpening term
          !
          !psimmm = seps*log((phi(i-1,j-1,k-1)+eps)/(1.-phi(i-1,j-1,k-1)+eps))
          psimcm = seps*log((phi(i-1,j  ,k-1)+eps)/(1.-phi(i-1,j  ,k-1)+eps))
          !psimpm = seps*log((phi(i-1,j+1,k-1)+eps)/(1.-phi(i-1,j+1,k-1)+eps))
          psicmm = seps*log((phi(i  ,j-1,k-1)+eps)/(1.-phi(i  ,j-1,k-1)+eps))
          psiccm = seps*log((phi(i  ,j  ,k-1)+eps)/(1.-phi(i  ,j  ,k-1)+eps))
          psicpm = seps*log((phi(i  ,j+1,k-1)+eps)/(1.-phi(i  ,j+1,k-1)+eps))
          !psipmm = seps*log((phi(i+1,j-1,k-1)+eps)/(1.-phi(i+1,j-1,k-1)+eps))
          psipcm = seps*log((phi(i+1,j  ,k-1)+eps)/(1.-phi(i+1,j  ,k-1)+eps))
          !psippm = seps*log((phi(i+1,j+1,k-1)+eps)/(1.-phi(i+1,j+1,k-1)+eps))
          psimmc = seps*log((phi(i-1,j-1,k  )+eps)/(1.-phi(i-1,j-1,k  )+eps))
          psimcc = seps*log((phi(i-1,j  ,k  )+eps)/(1.-phi(i-1,j  ,k  )+eps))
          psimpc = seps*log((phi(i-1,j+1,k  )+eps)/(1.-phi(i-1,j+1,k  )+eps))
          psicmc = seps*log((phi(i  ,j-1,k  )+eps)/(1.-phi(i  ,j-1,k  )+eps))
          psiccc = seps*log((phi(i  ,j  ,k  )+eps)/(1.-phi(i  ,j  ,k  )+eps))
          psicpc = seps*log((phi(i  ,j+1,k  )+eps)/(1.-phi(i  ,j+1,k  )+eps))
          psipmc = seps*log((phi(i+1,j-1,k  )+eps)/(1.-phi(i+1,j-1,k  )+eps))
          psipcc = seps*log((phi(i+1,j  ,k  )+eps)/(1.-phi(i+1,j  ,k  )+eps))
          psippc = seps*log((phi(i+1,j+1,k  )+eps)/(1.-phi(i+1,j+1,k  )+eps))
          !psimmp = seps*log((phi(i-1,j-1,k+1)+eps)/(1.-phi(i-1,j-1,k+1)+eps))
          psimcp = seps*log((phi(i-1,j  ,k+1)+eps)/(1.-phi(i-1,j  ,k+1)+eps))
          !psimpp = seps*log((phi(i-1,j+1,k+1)+eps)/(1.-phi(i-1,j+1,k+1)+eps))
          psicmp = seps*log((phi(i  ,j-1,k+1)+eps)/(1.-phi(i  ,j-1,k+1)+eps))
          psiccp = seps*log((phi(i  ,j  ,k+1)+eps)/(1.-phi(i  ,j  ,k+1)+eps))
          psicpp = seps*log((phi(i  ,j+1,k+1)+eps)/(1.-phi(i  ,j+1,k+1)+eps))
          !psipmp = seps*log((phi(i+1,j-1,k+1)+eps)/(1.-phi(i+1,j-1,k+1)+eps))
          psipcp = seps*log((phi(i+1,j  ,k+1)+eps)/(1.-phi(i+1,j  ,k+1)+eps))
          !psippp = seps*log((phi(i+1,j+1,k+1)+eps)/(1.-phi(i+1,j+1,k+1)+eps))
          !
          dpsidx_01 =      (psiccc       -psimcc       )*dxi
          dpsidy_01 = 0.25*(psicpc+psimpc-psicmc-psimmc)*dyi
          dpsidz_01 = 0.50*(psiccp+psimcp-psiccm-psimcm)*(dzci(k-1)*dzci(k))/(dzci(k-1)+dzci(k))
          dpsidx_11 =      (psipcc       -psiccc       )*dxi
          dpsidy_11 = 0.25*(psippc+psicpc-psipmc-psicmc)*dyi
          dpsidz_11 = 0.50*(psipcp+psiccp-psipcm-psiccm)*(dzci(k-1)*dzci(k))/(dzci(k-1)+dzci(k))
          dpsidx_02 = 0.25*(psippc+psipcc-psimpc-psimcc)*dxi
          dpsidy_02 =      (psicpc       -psiccc       )*dyi
          dpsidz_02 = 0.50*(psicpp+psiccp-psicpm-psiccm)*(dzci(k-1)*dzci(k))/(dzci(k-1)+dzci(k))
          dpsidx_12 = 0.25*(psipcc+psipmc-psimcc-psimmc)*dxi
          dpsidy_12 =      (psiccc       -psicmc       )*dyi
          dpsidz_12 = 0.50*(psiccp+psicmp-psiccm-psicmm)*(dzci(k-1)*dzci(k))/(dzci(k-1)+dzci(k))
          dpsidx_03 = 0.25*((psipcm-psimcm)*wghtmp+(psipcc-psimcc)*wghtmm)*dxi
          dpsidy_03 = 0.25*((psicpm-psicmm)*wghtmp+(psicpc-psicmc)*wghtmm)*dyi
          dpsidz_03 =      ( psiccc               -psiccm                )*dzci(k-1)
          dpsidx_13 = 0.25*((psipcp-psimcp)*wghtpm+(psipcc-psimcc)*wghtpp)*dxi
          dpsidy_13 = 0.25*((psicpp-psicmp)*wghtpm+(psicpc-psicmc)*wghtpp)*dyi
          dpsidz_13 =      ( psiccp               -psiccc                )*dzci(k)
          !
          rn_01 = dpsidx_01/(sqrt(dpsidx_01**2+dpsidy_01**2+dpsidz_01**2)+eps)
          rn_11 = dpsidx_11/(sqrt(dpsidx_11**2+dpsidy_11**2+dpsidz_11**2)+eps)
          rn_02 = dpsidy_02/(sqrt(dpsidx_02**2+dpsidy_02**2+dpsidz_02**2)+eps)
          rn_12 = dpsidy_12/(sqrt(dpsidx_12**2+dpsidy_12**2+dpsidz_12**2)+eps)
          rn_03 = dpsidz_03/(sqrt(dpsidx_03**2+dpsidy_03**2+dpsidz_03**2)+eps)
          rn_13 = dpsidz_13/(sqrt(dpsidx_13**2+dpsidy_13**2+dpsidz_13**2)+eps)
          !
          sharp = 0.25*gam*( (1.-(tanh(0.25*(psipcc       +psiccc       )*sepsi))**2)*rn_11 - &
                             (1.-(tanh(0.25*(psiccc       +psimcc       )*sepsi))**2)*rn_01 )*dxi + &
                  0.25*gam*( (1.-(tanh(0.25*(psicpc       +psiccc       )*sepsi))**2)*rn_02 - &
                             (1.-(tanh(0.25*(psiccc       +psicmc       )*sepsi))**2)*rn_12 )*dyi + &
                  0.25*gam*( (1.-(tanh(0.25*(psiccp*wghtpm+psiccc*wghtpp)*sepsi))**2)*rn_13 - &
                             (1.-(tanh(0.25*(psiccc*wghtmm+psiccm*wghtmp)*sepsi))**2)*rn_03 )*dzfi(k)
          dphidt(i,j,k) = -adv+diff-sharp
        end do
      end do
    end do
  end subroutine acdi_transport_pf
  pure elemental real(rp) function acdi_cmpt_psi(phi,seps) result(res)
    use mod_param, only:eps
    !$acc routine seq
    !
    ! smooth impulse Dirac delta function using trigonometric functions
    !
    implicit none
    !
    real(rp), intent(in) :: phi,seps
    res = seps*log((phi+eps)/(1.-phi+eps))
  end function acdi_cmpt_psi
end module mod_acdi
