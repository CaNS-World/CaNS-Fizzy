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
  public acdi_set_epsilon,acdi_set_gamma,acdi_transport_pf,acdi_cmpt_norm_curv,acdi_cmpt_rglr,acdi_cmpt_phi
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
  subroutine acdi_transport_pf(n,dli,dzci,dzfi,gam,seps,u,v,w,normx,normy,normz,psi,dpsidt)
    !
    !
    !
    implicit none
    integer , intent(in ), dimension(3) :: n
    real(rp), intent(in ), dimension(3) :: dli
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), intent(in)                :: gam,seps
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w
    real(rp), intent(in), dimension(0:,0:,0:) :: normx,normy,normz
    real(rp), dimension(0:,0:,0:), intent(inout) :: psi
    real(rp), dimension(:,:,:), intent(out) :: dpsidt
    integer :: i,j,k
    real(rp) :: dxi,dyi
    real(rp) :: sepsi
    real(rp) :: wghtpp,wghtpm,wghtmp,wghtmm
    real(rp) :: adv,diff,sharp
    real(rp) :: upsiip,upsiim,vpsijp,vpsijm,wpsikp,wpsikm
    real(rp) :: dpsidxp,dpsidxm,dpsidyp,dpsidym,dpsidzp,dpsidzm
    real(rp) :: phiccm,phimcc,phicmc,phiccc,phicpc,phipcc,phiccp
    real(rp) :: rn_01,rn_11,rn_02,rn_12,rn_03,rn_13
    !
    dxi = dli(1)
    dyi = dli(2)
    sepsi = 1./seps
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      !
      ! weights along the stretched grid direction
      !
      wghtpp = dzci(k)/dzfi(k+1)
      wghtpm = dzci(k)/dzfi(k)
      wghtmp = dzci(k-1)/dzfi(k)
      wghtmm = dzci(k-1)/dzfi(k-1)
      !
      do j=1,n(2)
        do i=1,n(1)
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
          !phimcm = seps*log((psi(i-1,j  ,k-1)+eps)/(1.-psi(i-1,j  ,k-1)+eps))
          !phimpm = seps*log((psi(i-1,j+1,k-1)+eps)/(1.-psi(i-1,j+1,k-1)+eps))
          !phicmm = seps*log((psi(i  ,j-1,k-1)+eps)/(1.-psi(i  ,j-1,k-1)+eps))
          phiccm = seps*log((psi(i  ,j  ,k-1)+eps)/(1.-psi(i  ,j  ,k-1)+eps))
          !phicpm = seps*log((psi(i  ,j+1,k-1)+eps)/(1.-psi(i  ,j+1,k-1)+eps))
          !phipmm = seps*log((psi(i+1,j-1,k-1)+eps)/(1.-psi(i+1,j-1,k-1)+eps))
          !phipcm = seps*log((psi(i+1,j  ,k-1)+eps)/(1.-psi(i+1,j  ,k-1)+eps))
          !phippm = seps*log((psi(i+1,j+1,k-1)+eps)/(1.-psi(i+1,j+1,k-1)+eps))
          !phimmc = seps*log((psi(i-1,j-1,k  )+eps)/(1.-psi(i-1,j-1,k  )+eps))
          phimcc = seps*log((psi(i-1,j  ,k  )+eps)/(1.-psi(i-1,j  ,k  )+eps))
          !phimpc = seps*log((psi(i-1,j+1,k  )+eps)/(1.-psi(i-1,j+1,k  )+eps))
          phicmc = seps*log((psi(i  ,j-1,k  )+eps)/(1.-psi(i  ,j-1,k  )+eps))
          phiccc = seps*log((psi(i  ,j  ,k  )+eps)/(1.-psi(i  ,j  ,k  )+eps))
          phicpc = seps*log((psi(i  ,j+1,k  )+eps)/(1.-psi(i  ,j+1,k  )+eps))
          !phipmc = seps*log((psi(i+1,j-1,k  )+eps)/(1.-psi(i+1,j-1,k  )+eps))
          phipcc = seps*log((psi(i+1,j  ,k  )+eps)/(1.-psi(i+1,j  ,k  )+eps))
          !phippc = seps*log((psi(i+1,j+1,k  )+eps)/(1.-psi(i+1,j+1,k  )+eps))
          !phimmp = seps*log((psi(i-1,j-1,k+1)+eps)/(1.-psi(i-1,j-1,k+1)+eps))
          !phimcp = seps*log((psi(i-1,j  ,k+1)+eps)/(1.-psi(i-1,j  ,k+1)+eps))
          !phimpp = seps*log((psi(i-1,j+1,k+1)+eps)/(1.-psi(i-1,j+1,k+1)+eps))
          !phicmp = seps*log((psi(i  ,j-1,k+1)+eps)/(1.-psi(i  ,j-1,k+1)+eps))
          phiccp = seps*log((psi(i  ,j  ,k+1)+eps)/(1.-psi(i  ,j  ,k+1)+eps))
          !phicpp = seps*log((psi(i  ,j+1,k+1)+eps)/(1.-psi(i  ,j+1,k+1)+eps))
          !phipmp = seps*log((psi(i+1,j-1,k+1)+eps)/(1.-psi(i+1,j-1,k+1)+eps))
          !phipcp = seps*log((psi(i+1,j  ,k+1)+eps)/(1.-psi(i+1,j  ,k+1)+eps))
          !phippp = seps*log((psi(i+1,j+1,k+1)+eps)/(1.-psi(i+1,j+1,k+1)+eps))
          !
          rn_01 = 0.5*(normx(i-1,j,k)+normx(i,j,k))/(sqrt((0.5*(normx(i-1,j,k)+normx(i,j,k)))**2+ &
                 (0.5*(normy(i-1,j,k)+normy(i,j,k)))**2+(0.5*(normz(i-1,j,k)+normz(i,j,k)))**2)+eps)
          rn_11 = 0.5*(normx(i,j,k)+normx(i+1,j,k))/(sqrt((0.5*(normx(i,j,k)+normx(i+1,j,k)))**2+ &
                 (0.5*(normy(i,j,k)+normy(i+1,j,k)))**2+(0.5*(normz(i,j,k)+normz(i+1,j,k)))**2)+eps)
          rn_02 = 0.5*(normy(i,j-1,k)+normy(i,j,k))/(sqrt((0.5*(normx(i,j-1,k)+normx(i,j,k)))**2+ &
                 (0.5*(normy(i,j-1,k)+normy(i,j,k)))**2+(0.5*(normz(i,j-1,k)+normz(i,j,k)))**2)+eps)
          rn_12 = 0.5*(normy(i,j,k)+normy(i,j+1,k))/(sqrt((0.5*(normx(i,j,k)+normx(i,j+1,k)))**2+ &
                 (0.5*(normy(i,j,k)+normy(i,j+1,k)))**2+(0.5*(normz(i,j,k)+normz(i,j+1,k)))**2)+eps)
          rn_03 = 0.5*(normz(i,j,k-1)*wghtmp+normz(i,j,k)*wghtmm)/(sqrt((0.5*(normx(i,j,k-1)*wghtmp+normx(i,j,k)*wghtmm))**2+ &
                 (0.5*(normy(i,j,k-1)*wghtmp+normy(i,j,k)*wghtmm))**2+(0.5*(normz(i,j,k-1)*wghtmp+normz(i,j,k)*wghtmm))**2)+eps)
          rn_13 = 0.5*(normz(i,j,k)*wghtpp+normz(i,j,k+1)*wghtpm)/(sqrt((0.5*(normx(i,j,k)*wghtpp+normx(i,j,k+1)*wghtpm))**2+ &
                 (0.5*(normy(i,j,k)*wghtpp+normy(i,j,k+1)*wghtpm))**2+(0.5*(normz(i,j,k)*wghtpp+normz(i,j,k+1)*wghtpm))**2)+eps)
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
  subroutine acdi_cmpt_norm_curv(n,dli,dzci,dzfi,seps,psi,kappa,normx,normy,normz)
    !
    ! computes the normals and curvature based on a phase field
    ! using finite-differences based on Youngs method
    !
    implicit none
    real(rp), parameter :: eps = epsilon(1._rp)
    integer , intent(in ), dimension(3) :: n
    real(rp), intent(in ), dimension(3) :: dli
    real(rp), intent(in ), dimension(0:)          :: dzci,dzfi
    real(rp), intent(in)                          :: seps
    real(rp), intent(in ), dimension(0:,0:,0:)    :: psi
    real(rp), intent(out), dimension(0:,0:,0:)    :: kappa
    real(rp), intent(out), dimension(0:,0:,0:)    :: normx,normy,normz
    real(rp) :: phimmm,phimcm,phimpm,phicmm,phiccm,phicpm,phipmm,phipcm,phippm, &
                phimmc,phimcc,phimpc,phicmc,phiccc,phicpc,phipmc,phipcc,phippc, &
                phimmp,phimcp,phimpp,phicmp,phiccp,phicpp,phipmp,phipcp,phippp
    real(rp) :: wghtpp,wghtpm,wghtmp,wghtmm
    real(rp) :: norm
    logical , save :: is_first = .true.
    real(rp), allocatable, dimension(:), save :: mx,my,mz
    integer  :: i,j,k,q
    !
    if(is_first) then
      is_first = .false.
      allocate(mx(8),my(8),mz(8))
      !$acc enter data create(mx,my,mz)
    end if
    !
    !$acc parallel loop collapse(3) default(present) private(mx,my,mz,norm) async(1)
    do k=1,n(3)
      !
      ! weights along the stretched grid direction
      !
      wghtpp = dzci(k)/dzfi(k+1)
      wghtpm = dzci(k)/dzfi(k)
      wghtmp = dzci(k-1)/dzfi(k)
      wghtmm = dzci(k-1)/dzfi(k-1)
      !
      do j=1,n(2)
        do i=1,n(1)
          !
          phimmm = seps*log((psi(i-1,j-1,k-1)+eps)/(1.-psi(i-1,j-1,k-1)+eps))
          phimcm = seps*log((psi(i-1,j  ,k-1)+eps)/(1.-psi(i-1,j  ,k-1)+eps))
          phimpm = seps*log((psi(i-1,j+1,k-1)+eps)/(1.-psi(i-1,j+1,k-1)+eps))
          phicmm = seps*log((psi(i  ,j-1,k-1)+eps)/(1.-psi(i  ,j-1,k-1)+eps))
          phiccm = seps*log((psi(i  ,j  ,k-1)+eps)/(1.-psi(i  ,j  ,k-1)+eps))
          phicpm = seps*log((psi(i  ,j+1,k-1)+eps)/(1.-psi(i  ,j+1,k-1)+eps))
          phipmm = seps*log((psi(i+1,j-1,k-1)+eps)/(1.-psi(i+1,j-1,k-1)+eps))
          phipcm = seps*log((psi(i+1,j  ,k-1)+eps)/(1.-psi(i+1,j  ,k-1)+eps))
          phippm = seps*log((psi(i+1,j+1,k-1)+eps)/(1.-psi(i+1,j+1,k-1)+eps))
          phimmc = seps*log((psi(i-1,j-1,k  )+eps)/(1.-psi(i-1,j-1,k  )+eps))
          phimcc = seps*log((psi(i-1,j  ,k  )+eps)/(1.-psi(i-1,j  ,k  )+eps))
          phimpc = seps*log((psi(i-1,j+1,k  )+eps)/(1.-psi(i-1,j+1,k  )+eps))
          phicmc = seps*log((psi(i  ,j-1,k  )+eps)/(1.-psi(i  ,j-1,k  )+eps))
          phiccc = seps*log((psi(i  ,j  ,k  )+eps)/(1.-psi(i  ,j  ,k  )+eps))
          phicpc = seps*log((psi(i  ,j+1,k  )+eps)/(1.-psi(i  ,j+1,k  )+eps))
          phipmc = seps*log((psi(i+1,j-1,k  )+eps)/(1.-psi(i+1,j-1,k  )+eps))
          phipcc = seps*log((psi(i+1,j  ,k  )+eps)/(1.-psi(i+1,j  ,k  )+eps))
          phippc = seps*log((psi(i+1,j+1,k  )+eps)/(1.-psi(i+1,j+1,k  )+eps))
          phimmp = seps*log((psi(i-1,j-1,k+1)+eps)/(1.-psi(i-1,j-1,k+1)+eps))
          phimcp = seps*log((psi(i-1,j  ,k+1)+eps)/(1.-psi(i-1,j  ,k+1)+eps))
          phimpp = seps*log((psi(i-1,j+1,k+1)+eps)/(1.-psi(i-1,j+1,k+1)+eps))
          phicmp = seps*log((psi(i  ,j-1,k+1)+eps)/(1.-psi(i  ,j-1,k+1)+eps))
          phiccp = seps*log((psi(i  ,j  ,k+1)+eps)/(1.-psi(i  ,j  ,k+1)+eps))
          phicpp = seps*log((psi(i  ,j+1,k+1)+eps)/(1.-psi(i  ,j+1,k+1)+eps))
          phipmp = seps*log((psi(i+1,j-1,k+1)+eps)/(1.-psi(i+1,j-1,k+1)+eps))
          phipcp = seps*log((psi(i+1,j  ,k+1)+eps)/(1.-psi(i+1,j  ,k+1)+eps))
          phippp = seps*log((psi(i+1,j+1,k+1)+eps)/(1.-psi(i+1,j+1,k+1)+eps))
          !
          mx(1) = 0.25*((phipcc+phippc-phiccc-phicpc)*wghtpp+(phipcp+phippp-phiccp-phicpp)*wghtpm)*dli(1)
          mx(2) = 0.25*((phipcc+phipmc-phiccc-phicmc)*wghtpp+(phipcp+phipmp-phiccp-phicmp)*wghtpm)*dli(1)
          mx(3) = 0.25*((phipcc+phippc-phiccc-phicpc)*wghtmm+(phipcm+phippm-phiccm-phicpm)*wghtmp)*dli(1)
          mx(4) = 0.25*((phipcc+phipmc-phiccc-phicmc)*wghtmm+(phipcm+phipmm-phiccm-phicmm)*wghtmp)*dli(1)
          mx(5) = 0.25*((phiccc+phicpc-phimcc-phimpc)*wghtpp+(phiccp+phicpp-phimcp-phimpp)*wghtpm)*dli(1)
          mx(6) = 0.25*((phiccc+phicmc-phimcc-phimmc)*wghtpp+(phiccp+phicmp-phimcp-phimmp)*wghtpm)*dli(1)
          mx(7) = 0.25*((phiccc+phicpc-phimcc-phimpc)*wghtmm+(phiccm+phicpm-phimcm-phimpm)*wghtmp)*dli(1)
          mx(8) = 0.25*((phiccc+phicmc-phimcc-phimmc)*wghtmm+(phiccm+phicmm-phimcm-phimmm)*wghtmp)*dli(1)
          !
          my(1) = 0.25*((phicpc+phippc-phiccc-phipcc)*wghtpp+(phicpp+phippp-phiccp-phipcp)*wghtpm)*dli(2)
          my(2) = 0.25*((phiccc+phipcc-phicmc-phipmc)*wghtpp+(phiccp+phipcp-phicmp-phipmp)*wghtpm)*dli(2)
          my(3) = 0.25*((phicpc+phippc-phiccc-phipcc)*wghtmm+(phicpm+phippm-phiccm-phipcm)*wghtmp)*dli(2)
          my(4) = 0.25*((phiccc+phipcc-phicmc-phipmc)*wghtmm+(phiccm+phipcm-phicmm-phipmm)*wghtmp)*dli(2)
          my(5) = 0.25*((phicpc+phimpc-phiccc-phimcc)*wghtpp+(phicpp+phimpp-phiccp-phimcp)*wghtpm)*dli(2)
          my(6) = 0.25*((phiccc+phimcc-phicmc-phimmc)*wghtpp+(phiccp+phimcp-phicmp-phimmp)*wghtpm)*dli(2)
          my(7) = 0.25*((phicpc+phimpc-phiccc-phimcc)*wghtmm+(phicpm+phimpm-phiccm-phimcm)*wghtmp)*dli(2)
          my(8) = 0.25*((phiccc+phimcc-phicmc-phimmc)*wghtmm+(phiccm+phimcm-phicmm-phimmm)*wghtmp)*dli(2)
          !
          mz(1) = 0.25*((phiccp+phipcp+phicpp+phippp)       -(phiccc+phipcc+phicpc+phippc))*dzci(k  )
          mz(2) = 0.25*((phiccp+phipcp+phicmp+phipmp)       -(phiccc+phipcc+phicmc+phipmc))*dzci(k  )
          mz(3) = 0.25*((phiccc+phipcc+phicpc+phippc)       -(phiccm+phipcm+phicpm+phippm))*dzci(k-1)
          mz(4) = 0.25*((phiccc+phipcc+phicmc+phipmc)       -(phiccm+phipcm+phicmm+phipmm))*dzci(k-1)
          mz(5) = 0.25*((phiccp+phimcp+phicpp+phimpp)       -(phiccc+phimcc+phicpc+phimpc))*dzci(k  )
          mz(6) = 0.25*((phiccp+phimcp+phicmp+phimmp)       -(phiccc+phimcc+phicmc+phimmc))*dzci(k  )
          mz(7) = 0.25*((phiccc+phimcc+phicpc+phimpc)       -(phiccm+phimcm+phicpm+phimpm))*dzci(k-1)
          mz(8) = 0.25*((phiccc+phimcc+phicmc+phimmc)       -(phiccm+phimcm+phicmm+phimmm))*dzci(k-1)
          !
          !$acc loop seq
          do q=1,8
            norm = sqrt(mx(q)**2+my(q)**2+mz(q)**2)+eps
            mx(q) = mx(q)/norm
            my(q) = my(q)/norm
            mz(q) = mz(q)/norm
          end do
          !
          ! compute the normal vector
          !
          normx(i,j,k) = .125*(mx(1)+mx(2)+mx(3)+mx(4)+mx(5)+mx(6)+mx(7)+mx(8))
          normy(i,j,k) = .125*(my(1)+my(2)+my(3)+my(4)+my(5)+my(6)+my(7)+my(8))
          normz(i,j,k) = .125*(mz(1)+mz(2)+mz(3)+mz(4)+mz(5)+mz(6)+mz(7)+mz(8))
          norm = sqrt(normx(i,j,k)**2+normy(i,j,k)**2+normz(i,j,k)**2)+eps
          normx(i,j,k) = normx(i,j,k)/norm
          normy(i,j,k) = normy(i,j,k)/norm
          normz(i,j,k) = normz(i,j,k)/norm
          !
          ! compute the curvature
          !
          kappa(i,j,k) = - ( 0.25*((mx(1)+mx(2)+mx(3)+mx(4))-(mx(5)+mx(6)+mx(7)+mx(8)))*dli(1) + &
                             0.25*((my(1)+my(3)+my(5)+my(7))-(my(2)+my(4)+my(6)+my(8)))*dli(2) + &
                             0.25*((mz(1)+mz(2)+mz(5)+mz(6))-(mz(3)+mz(4)+mz(7)+mz(8)))*dzfi(k) )
        end do
      end do
    end do
  end subroutine acdi_cmpt_norm_curv
  !
  subroutine acdi_cmpt_rglr(n,dli,dzci,dzfi,gam,seps,normx,normy,normz,psi,rglrx,rglry,rglrz)
    !
    ! compute the components of the interface regularization vector (staggered like the velocity)
    !
    implicit none
    integer , intent(in), dimension(3)  :: n
    real(rp), intent(in), dimension(3)  :: dli
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), intent(in)                :: gam,seps
    real(rp), intent(in ), dimension(0:,0:,0:) :: normx,normy,normz
    real(rp), intent(in ), dimension(0:,0:,0:) :: psi
    real(rp), intent(out), dimension(0:,0:,0:) :: rglrx,rglry,rglrz
    integer :: i,j,k
    real(rp) :: dxi,dyi
    real(rp) :: sepsi
    real(rp) :: wghtpp,wghtpm,wghtmp,wghtmm
    real(rp) :: diffx,diffy,diffz,sharpx,sharpy,sharpz
    real(rp) :: phiccc,phicpc,phipcc,phiccp
    real(rp) :: rn_11,rn_12,rn_13
    !
    dxi = dli(1)
    dyi = dli(2)
    sepsi = 1./seps
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      !
      ! weights along the stretched grid direction
      !
      wghtpp = dzci(k)/dzfi(k+1)
      wghtpm = dzci(k)/dzfi(k)
      wghtmp = dzci(k-1)/dzfi(k)
      wghtmm = dzci(k-1)/dzfi(k-1)
      !
      do j=1,n(2)
        do i=1,n(1)
          !
          ! diffusion term
          !
          diffx = gam*seps*(psi(i+1,j,k)-psi(i,j,k))*dxi
          diffy = gam*seps*(psi(i,j+1,k)-psi(i,j,k))*dyi
          diffz = gam*seps*(psi(i,j,k+1)-psi(i,j,k))*dzci(k)
          !
          ! sharpening term
          !
          phiccc = seps*log((psi(i  ,j  ,k  )+eps)/(1.-psi(i  ,j  ,k  )+eps))
          phicpc = seps*log((psi(i  ,j+1,k  )+eps)/(1.-psi(i  ,j+1,k  )+eps))
          phipcc = seps*log((psi(i+1,j  ,k  )+eps)/(1.-psi(i+1,j  ,k  )+eps))
          phiccp = seps*log((psi(i  ,j  ,k+1)+eps)/(1.-psi(i  ,j  ,k+1)+eps))
          !
          rn_11 = 0.5*(normx(i,j,k)+normx(i+1,j,k))/(sqrt((0.5*(normx(i,j,k)+normx(i+1,j,k)))**2+ &
                 (0.5*(normy(i,j,k)+normy(i+1,j,k)))**2+(0.5*(normz(i,j,k)+normz(i+1,j,k)))**2)+eps)
          rn_12 = 0.5*(normy(i,j,k)+normy(i,j+1,k))/(sqrt((0.5*(normx(i,j,k)+normx(i,j+1,k)))**2+ &
                 (0.5*(normy(i,j,k)+normy(i,j+1,k)))**2+(0.5*(normz(i,j,k)+normz(i,j+1,k)))**2)+eps)
          rn_13 = 0.5*(normz(i,j,k)*wghtpp+normz(i,j,k+1)*wghtpm)/(sqrt((0.5*(normx(i,j,k)*wghtpp+normx(i,j,k+1)*wghtpm))**2+ &
                 (0.5*(normy(i,j,k)*wghtpp+normy(i,j,k+1)*wghtpm))**2+(0.5*(normz(i,j,k)*wghtpp+normz(i,j,k+1)*wghtpm))**2)+eps)
          !  
          sharpx = 0.25*gam*((1.-(tanh(0.25*(phipcc       +phiccc       )*sepsi))**2)*rn_11)
          sharpy = 0.25*gam*((1.-(tanh(0.25*(phicpc       +phiccc       )*sepsi))**2)*rn_12)
          sharpz = 0.25*gam*((1.-(tanh(0.25*(phiccp*wghtpm+phiccc*wghtpp)*sepsi))**2)*rn_13)
          !
          rglrx(i,j,k) = diffx-sharpx
          rglry(i,j,k) = diffy-sharpy
          rglrz(i,j,k) = diffz-sharpz
        end do
      end do
    end do
  end subroutine acdi_cmpt_rglr
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
