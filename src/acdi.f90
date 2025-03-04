! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 The CaNS contributors. All rights reserved.
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
  public acdi_set_epsilon,acdi_set_gamma,acdi_transport_pf,acdi_cmpt_phi
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
  subroutine acdi_transport_pf(n,dli,dzci,dzfi,gam,seps,u,v,w,normx,normy,normz,phi,psi,dpsidt,flux_x,flux_y,flux_z)
    !
    ! compute right-hand side of the phase field transport equation
    !
    implicit none
    integer , intent(in   ), dimension(3)        :: n
    real(rp), intent(in   ), dimension(3)        :: dli
    real(rp), intent(in   ), dimension(0:)       :: dzci,dzfi
    real(rp), intent(in   )                      :: gam,seps
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(in   ), dimension(0:,0:,0:) :: normx,normy,normz
    real(rp), intent(in   ), dimension(0:,0:,0:) :: phi
    real(rp), intent(inout), dimension(0:,0:,0:) :: psi
    real(rp), intent(out  ), dimension(: ,: ,: ) :: dpsidt
    real(rp), intent(out  ), dimension(0:,0:,0:), optional :: flux_x,flux_y,flux_z
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
          phi_ccm = phi(i  ,j  ,k-1)
          phi_cmc = phi(i  ,j-1,k  )
          phi_mcc = phi(i-1,j  ,k  )
          phi_ccc = phi(i  ,j  ,k  )
          phi_pcc = phi(i+1,j  ,k  )
          phi_cpc = phi(i  ,j+1,k  )
          phi_ccp = phi(i  ,j  ,k+1)
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
          dpsidxm = gam*seps*(psi_ccc-psi_mcc)*dxi
          dpsidxp = gam*seps*(psi_pcc-psi_ccc)*dxi
          dpsidym = gam*seps*(psi_ccc-psi_cmc)*dyi
          dpsidyp = gam*seps*(psi_cpc-psi_ccc)*dyi
          dpsidzm = gam*seps*(psi_ccc-psi_ccm)*dzci_m
          dpsidzp = gam*seps*(psi_ccp-psi_ccc)*dzci_c
          !
          ! sharpening term
          !
          rn_01 = normx_xm/(sqrt(normx_xm**2+normy_xm**2+normz_xm**2)+eps)
          rn_11 = normx_xp/(sqrt(normx_xp**2+normy_xp**2+normz_xp**2)+eps)
          rn_02 = normy_ym/(sqrt(normx_ym**2+normy_ym**2+normz_ym**2)+eps)
          rn_12 = normy_yp/(sqrt(normx_yp**2+normy_yp**2+normz_yp**2)+eps)
          rn_03 = normz_zm/(sqrt(normx_zm**2+normy_zm**2+normz_zm**2)+eps)
          rn_13 = normz_zp/(sqrt(normx_zp**2+normy_zp**2+normz_zp**2)+eps)
          !
          sharpxm = 0.25*gam*((1.-(tanh(0.25*(phi_ccc+phi_mcc)/seps))**2)*rn_01)
          sharpxp = 0.25*gam*((1.-(tanh(0.25*(phi_pcc+phi_ccc)/seps))**2)*rn_11)
          sharpym = 0.25*gam*((1.-(tanh(0.25*(phi_ccc+phi_cmc)/seps))**2)*rn_02)
          sharpyp = 0.25*gam*((1.-(tanh(0.25*(phi_cpc+phi_ccc)/seps))**2)*rn_12)
          sharpzm = 0.25*gam*((1.-(tanh(0.25*(phi_ccc+phi_ccm)/seps))**2)*rn_03)
          sharpzp = 0.25*gam*((1.-(tanh(0.25*(phi_ccp+phi_ccc)/seps))**2)*rn_13)
          !
          ! transport
          !
          rglrxm = dpsidxm - sharpxm
          rglrxp = dpsidxp - sharpxp
          rglrym = dpsidym - sharpym
          rglryp = dpsidyp - sharpyp
          rglrzm = dpsidzm - sharpzm
          rglrzp = dpsidzp - sharpzp
          !
          rglr = (rglrxp-rglrxm)*dxi + (rglryp-rglrym)*dyi + (rglrzp-rglrzm)*dzfi_c
          !
          dpsidt(i,j,k) = -adv + rglr
#if defined(_CONSISTENT_ADVECTION)
          flux_x(i,j,k) = upsiip - rglrxp
          flux_y(i,j,k) = vpsijp - rglryp
          flux_z(i,j,k) = wpsikp - rglrzp
#endif
        end do
      end do
    end do
  end subroutine acdi_transport_pf
  !
  subroutine acdi_cmpt_phi(n,seps,psi,phi)
    !
    ! computes the signed-distance field phi from the smoothed volume fraction field psi
    !
    implicit none
    integer , intent(in ), dimension(3)           :: n
    real(rp), intent(in )                         :: seps
    real(rp), intent(in ), dimension(0:,0:,0:)    :: psi
    real(rp), intent(out), dimension(0:,0:,0:)    :: phi
    real(rp) :: psi_aux
    integer  :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          psi_aux = psi(i,j,k)
          phi(i,j,k) = seps*log((psi_aux+eps)/(1.-psi_aux+eps))
        end do
      end do
    end do
  end subroutine acdi_cmpt_phi
  !
  pure elemental real(rp) function acdi_phi(psi,seps) result(phi)
    use mod_param, only:eps
    !$acc routine seq
    !
    ! computes the signed-distance field phi from the smoothed volume fraction field psi
    !
    implicit none
    real(rp), intent(in) :: psi,seps
    phi = seps*log((psi+eps)/(1.-psi+eps))
  end function acdi_phi
end module mod_acdi
