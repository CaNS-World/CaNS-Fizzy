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
  public acdi_set_epsilon,acdi_set_gamma,acdi_transport_pf,acdi_cmpt_norm_curv
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
    !$acc data copy(velmax) async(1)
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
  subroutine acdi_transport_pf(n,dli,dzci,dzfi,gam,seps,u,v,w,normx,normy,normz,psi,dpsidt,rglrx,rglry,rglrz)
    !
    !
    !
    implicit none
    integer , intent(in   ), dimension(3)        :: n
    real(rp), intent(in   ), dimension(3)        :: dli
    real(rp), intent(in   ), dimension(0:)       :: dzci,dzfi
    real(rp), intent(in   )                      :: gam,seps
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(in   ), dimension(0:,0:,0:) :: normx,normy,normz
    real(rp), intent(inout), dimension(0:,0:,0:) :: psi
    real(rp), intent(out  ), dimension(: ,: ,: ) :: dpsidt
    real(rp), intent(out  ), dimension(0:,0:,0:), optional :: rglrx,rglry,rglrz
    integer :: i,j,k
    real(rp) :: dxi,dyi
    real(rp) :: adv,diff,sharp,rglr
    real(rp) :: upsiip,upsiim,vpsijp,vpsijm,wpsikp,wpsikm
    real(rp) :: dpsidxp,dpsidxm,dpsidyp,dpsidym,dpsidzp,dpsidzm
    real(rp) :: sharpxp,sharpxm,sharpyp,sharpym,sharpzp,sharpzm
    real(rp) :: rglrxp,rglrxm,rglryp,rglrym,rglrzp,rglrzm
    real(rp) :: dzfi_c,dzci_c,dzci_m
    real(rp) :: psi_ccm,psi_mcc,psi_cmc,psi_ccc,psi_cpc,psi_pcc,psi_ccp
    real(rp) :: u_mcc,u_ccc,v_cmc,v_ccc,w_ccm,w_ccc
    real(rp) :: normx_xm,normx_xp,normx_ym,normx_yp,normx_zm,normx_zp, &
                normy_xm,normy_xp,normy_ym,normy_yp,normy_zm,normy_zp, &
                normz_xm,normz_xp,normz_ym,normz_yp,normz_zm,normz_zp
    real(rp) :: phi_ccm,phi_mcc,phi_cmc,phi_ccc,phi_cpc,phi_pcc,phi_ccp
    real(rp) :: rn_01,rn_11,rn_02,rn_12,rn_03,rn_13
    !
    dxi = dli(1)
    dyi = dli(2)
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          ! touch each array at a time to maximize cache efficiency
          !
          dzfi_c  = dzfi(k)
          dzci_c  = dzci(k)
          dzci_m  = dzci(k-1)
          !
          psi_ccm = psi(i  ,j  ,k-1)
          psi_cmc = psi(i  ,j-1,k  )
          psi_mcc = psi(i-1,j  ,k  )
          psi_ccc = psi(i  ,j  ,k  )
          psi_pcc = psi(i+1,j  ,k  )
          psi_cpc = psi(i  ,j+1,k  )
          psi_ccp = psi(i  ,j  ,k+1)
          !
          u_mcc = u(i-1,j  ,k  )
          u_ccc = u(i  ,j  ,k  )
          v_cmc = v(i  ,j-1,k  )
          v_ccc = v(i  ,j  ,k  )
          w_ccm = w(i  ,j  ,k-1)
          w_ccc = w(i  ,j  ,k  )
          !
          normx_xm = 0.5*(normx(i-1,j,k)+normx(i,j,k))
          normx_xp = 0.5*(normx(i+1,j,k)+normx(i,j,k))
          normx_ym = 0.5*(normx(i,j-1,k)+normx(i,j,k))
          normx_yp = 0.5*(normx(i,j+1,k)+normx(i,j,k))
          normx_zm = 0.5*(normx(i,j,k-1)+normx(i,j,k))
          normx_zp = 0.5*(normx(i,j,k+1)+normx(i,j,k))
          normy_xm = 0.5*(normy(i-1,j,k)+normy(i,j,k))
          normy_xp = 0.5*(normy(i+1,j,k)+normy(i,j,k))
          normy_ym = 0.5*(normy(i,j-1,k)+normy(i,j,k))
          normy_yp = 0.5*(normy(i,j+1,k)+normy(i,j,k))
          normy_zm = 0.5*(normy(i,j,k-1)+normy(i,j,k))
          normy_zp = 0.5*(normy(i,j,k+1)+normy(i,j,k))
          normz_xm = 0.5*(normz(i-1,j,k)+normz(i,j,k))
          normz_xp = 0.5*(normz(i+1,j,k)+normz(i,j,k))
          normz_ym = 0.5*(normz(i,j-1,k)+normz(i,j,k))
          normz_yp = 0.5*(normz(i,j+1,k)+normz(i,j,k))
          normz_zm = 0.5*(normz(i,j,k-1)+normz(i,j,k))
          normz_zp = 0.5*(normz(i,j,k+1)+normz(i,j,k))
          !
          ! advection term
          !
          upsiim  = 0.5*(psi_mcc+psi_ccc)*u_mcc
          upsiip  = 0.5*(psi_pcc+psi_ccc)*u_ccc
          vpsijm  = 0.5*(psi_cmc+psi_ccc)*v_cmc
          vpsijp  = 0.5*(psi_cpc+psi_ccc)*v_ccc
          wpsikm  = 0.5*(psi_ccm+psi_ccc)*w_ccm
          wpsikp  = 0.5*(psi_ccp+psi_ccc)*w_ccc
          adv     = dxi*(upsiip-upsiim) + dyi*(vpsijp-vpsijm) + dzfi_c*(wpsikp-wpsikm)
          !
          ! diffusion term
          !
          dpsidxp = gam*seps*(psi_pcc-psi_ccc)*dxi
          dpsidxm = gam*seps*(psi_ccc-psi_mcc)*dxi
          dpsidyp = gam*seps*(psi_cpc-psi_ccc)*dyi
          dpsidym = gam*seps*(psi_ccc-psi_cmc)*dyi
          dpsidzp = gam*seps*(psi_ccp-psi_ccc)*dzci_c
          dpsidzm = gam*seps*(psi_ccc-psi_ccm)*dzci_m
          !
          ! sharpening term
          !
          phi_ccm = seps*log((psi_ccm+eps)/(1.-psi_ccm+eps))
          phi_cmc = seps*log((psi_cmc+eps)/(1.-psi_cmc+eps))
          phi_mcc = seps*log((psi_mcc+eps)/(1.-psi_mcc+eps))
          phi_ccc = seps*log((psi_ccc+eps)/(1.-psi_ccc+eps))
          phi_pcc = seps*log((psi_pcc+eps)/(1.-psi_pcc+eps))
          phi_cpc = seps*log((psi_cpc+eps)/(1.-psi_cpc+eps))
          phi_ccp = seps*log((psi_ccp+eps)/(1.-psi_ccp+eps))
          !
          rn_01 = normx_xm/sqrt(normx_xm**2+normy_xm**2+normz_xm**2+eps)
          rn_11 = normx_xp/sqrt(normx_xp**2+normy_xp**2+normz_xp**2+eps)
          rn_02 = normy_ym/sqrt(normx_ym**2+normy_ym**2+normz_ym**2+eps)
          rn_12 = normy_yp/sqrt(normx_yp**2+normy_yp**2+normz_yp**2+eps)
          rn_03 = normz_zm/sqrt(normx_zm**2+normy_zm**2+normz_zm**2+eps)
          rn_13 = normz_zp/sqrt(normx_zp**2+normy_zp**2+normz_zp**2+eps)
          !
          sharpxp = 0.25*gam*((1.-(tanh(0.25*(phi_pcc+phi_ccc)/seps))**2)*rn_11)
          sharpxm = 0.25*gam*((1.-(tanh(0.25*(phi_ccc+phi_mcc)/seps))**2)*rn_01)
          sharpyp = 0.25*gam*((1.-(tanh(0.25*(phi_cpc+phi_ccc)/seps))**2)*rn_12)
          sharpym = 0.25*gam*((1.-(tanh(0.25*(phi_ccc+phi_cmc)/seps))**2)*rn_02)
          sharpzp = 0.25*gam*((1.-(tanh(0.25*(phi_ccp+phi_ccc)/seps))**2)*rn_13)
          sharpzm = 0.25*gam*((1.-(tanh(0.25*(phi_ccc+phi_ccm)/seps))**2)*rn_03)
          !
          ! transport
          !
          rglrxp = dpsidxp - sharpxp
          rglrxm = dpsidxm - sharpxm
          rglryp = dpsidyp - sharpyp
          rglrym = dpsidym - sharpym
          rglrzp = dpsidzp - sharpzp
          rglrzm = dpsidzm - sharpzm
          !
          rglr = (rglrxp-rglrxm)*dxi + (rglryp-rglrym)*dyi + (rglrzp-rglrzm)*dzfi_c
          !
          dpsidt(i,j,k) = -adv + rglr
#if defined(_CONSERVATIVE_MOMENTUM)
          rglrx(i,j,k) = rglrxp
          rglry(i,j,k) = rglryp
          rglrz(i,j,k) = rglrzp
#endif
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
    !$acc parallel loop collapse(3) default(present) private(mx,my,mz) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
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
          mx(1) = 0.25*((phipcc+phippc+phipcp+phippp)-(phiccc+phicpc+phiccp+phicpp))*dli(1)
          mx(2) = 0.25*((phipcc+phipmc+phipcp+phipmp)-(phiccc+phicmc+phiccp+phicmp))*dli(1)
          mx(3) = 0.25*((phipcc+phippc+phipcm+phippm)-(phiccc+phicpc+phiccm+phicpm))*dli(1)
          mx(4) = 0.25*((phipcc+phipmc+phipcm+phipmm)-(phiccc+phicmc+phiccm+phicmm))*dli(1)
          mx(5) = 0.25*((phiccc+phicpc+phiccp+phicpp)-(phimcc+phimpc+phimcp+phimpp))*dli(1)
          mx(6) = 0.25*((phiccc+phicmc+phiccp+phicmp)-(phimcc+phimmc+phimcp+phimmp))*dli(1)
          mx(7) = 0.25*((phiccc+phicpc+phiccm+phicpm)-(phimcc+phimpc+phimcm+phimpm))*dli(1)
          mx(8) = 0.25*((phiccc+phicmc+phiccm+phicmm)-(phimcc+phimmc+phimcm+phimmm))*dli(1)
          !
          my(1) = 0.25*((phicpc+phippc+phicpp+phippp)-(phiccc+phipcc+phiccp+phipcp))*dli(2)
          my(2) = 0.25*((phiccc+phipcc+phiccp+phipcp)-(phicmc+phipmc+phicmp+phipmp))*dli(2)
          my(3) = 0.25*((phicpc+phippc+phicpm+phippm)-(phiccc+phipcc+phiccm+phipcm))*dli(2)
          my(4) = 0.25*((phiccc+phipcc+phiccm+phipcm)-(phicmc+phipmc+phicmm+phipmm))*dli(2)
          my(5) = 0.25*((phicpc+phimpc+phicpp+phimpp)-(phiccc+phimcc+phiccp+phimcp))*dli(2)
          my(6) = 0.25*((phiccc+phimcc+phiccp+phimcp)-(phicmc+phimmc+phicmp+phimmp))*dli(2)
          my(7) = 0.25*((phicpc+phimpc+phicpm+phimpm)-(phiccc+phimcc+phiccm+phimcm))*dli(2)
          my(8) = 0.25*((phiccc+phimcc+phiccm+phimcm)-(phicmc+phimmc+phicmm+phimmm))*dli(2)
          !
          mz(1) = 0.25*((phiccp+phipcp+phicpp+phippp)-(phiccc+phipcc+phicpc+phippc))*dzci(k  )
          mz(2) = 0.25*((phiccp+phipcp+phicmp+phipmp)-(phiccc+phipcc+phicmc+phipmc))*dzci(k  )
          mz(3) = 0.25*((phiccc+phipcc+phicpc+phippc)-(phiccm+phipcm+phicpm+phippm))*dzci(k-1)
          mz(4) = 0.25*((phiccc+phipcc+phicmc+phipmc)-(phiccm+phipcm+phicmm+phipmm))*dzci(k-1)
          mz(5) = 0.25*((phiccp+phimcp+phicpp+phimpp)-(phiccc+phimcc+phicpc+phimpc))*dzci(k  )
          mz(6) = 0.25*((phiccp+phimcp+phicmp+phimmp)-(phiccc+phimcc+phicmc+phimmc))*dzci(k  )
          mz(7) = 0.25*((phiccc+phimcc+phicpc+phimpc)-(phiccm+phimcm+phicpm+phimpm))*dzci(k-1)
          mz(8) = 0.25*((phiccc+phimcc+phicmc+phimmc)-(phiccm+phimcm+phicmm+phimmm))*dzci(k-1)
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
