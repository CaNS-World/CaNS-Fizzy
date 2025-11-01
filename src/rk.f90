! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 The CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
#define _FAST_MOM_KERNELS
module mod_rk
  use mod_mom  , only: mom_xyz_ad,mom_xyz_oth
  use mod_utils, only: swap
  use mod_types
  implicit none
  private
  public rk,rk_scal,rk_2fl
  contains
  subroutine rk(rkpar,n,dli,dzci,dzfi,dt,dt_r, &
                bforce,gacc,sigma,rho_av,rho12,mu12,beta12,rho0,psi,kappa,p,pn,po,s, &
                psio,psiflx_x,psiflx_y,psiflx_z,u,v,w)
    !
    ! RK3 scheme for time integration of the momentum equations
    !
    implicit none
    real(rp), intent(in), dimension(2) :: rkpar
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: dli
    real(rp), intent(in) :: dt,dt_r
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi
    real(rp), intent(in   ), dimension(3)  :: bforce,gacc
    real(rp), intent(in   )                :: sigma,rho_av(3)
    real(rp), intent(in   ), dimension(2)  :: rho12,mu12,beta12
    real(rp), intent(in   )                :: rho0
    real(rp), intent(in   ), dimension(0:,0:,0:)           :: psi,kappa
    real(rp), intent(in   ), dimension(0:,0:,0:)           :: p
    real(rp), intent(in   ), dimension(0:,0:,0:), optional :: pn,po,s
    real(rp), intent(in   ), dimension(0:,0:,0:), optional :: psio,psiflx_x,psiflx_y,psiflx_z
    real(rp), intent(inout), dimension(0:,0:,0:)           :: u,v,w
    real(rp), target       , allocatable, dimension(:,:,:), save :: dudtrk_t ,dvdtrk_t ,dwdtrk_t , &
                                                                    dudtrko_t,dvdtrko_t,dwdtrko_t
    real(rp), pointer      , contiguous , dimension(:,:,:), save :: dudtrk   ,dvdtrk   ,dwdtrk   , &
                                                                    dudtrko  ,dvdtrko  ,dwdtrko
    real(rp) :: rho,drho,rhox_n,rhox_p,rhoy_n,rhoy_p,rhoz_n,rhoz_p
    logical, save :: is_first = .true.
    real(rp) :: factor1,factor2,factor12
    integer :: i,j,k
    !
    factor1  = rkpar(1)*dt
    factor2  = rkpar(2)*dt
    factor12 = factor1+factor2
    !
    ! initialization
    !
    if(is_first) then ! leverage save attribute to allocate these arrays on the device only once
      allocate(dudtrk_t( n(1),n(2),n(3)),dvdtrk_t( n(1),n(2),n(3)),dwdtrk_t( n(1),n(2),n(3)))
      allocate(dudtrko_t(n(1),n(2),n(3)),dvdtrko_t(n(1),n(2),n(3)),dwdtrko_t(n(1),n(2),n(3)))
      !$acc enter data create(dudtrk_t ,dvdtrk_t ,dwdtrk_t ) async(1)
      !$acc enter data create(dudtrko_t,dvdtrko_t,dwdtrko_t) async(1)
      dudtrk  => dudtrk_t
      dvdtrk  => dvdtrk_t
      dwdtrk  => dwdtrk_t
      dudtrko => dudtrko_t
      dvdtrko => dvdtrko_t
      dwdtrko => dwdtrko_t
    end if
    !
    ! advection, diffusion and regularization terms
    !
    call mom_xyz_ad(n,dli,dzci,dzfi,rho12,mu12,psiflx_x,psiflx_y,psiflx_z,u,v,w,psi,dudtrk,dvdtrk,dwdtrk)
    !
    if(is_first) then
      !$acc kernels default(present) async(1)
      dudtrko(:,:,:) = dudtrk(:,:,:)
      dvdtrko(:,:,:) = dvdtrk(:,:,:)
      dwdtrko(:,:,:) = dwdtrk(:,:,:)
      !$acc end kernels
      is_first = .false.
    end if
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          rhox_n = rho + drho*0.5*(psio(i,j,k)+psio(i+1,j,k))
          rhoy_n = rho + drho*0.5*(psio(i,j,k)+psio(i,j+1,k))
          rhoz_n = rho + drho*0.5*(psio(i,j,k)+psio(i,j,k+1))
          rhox_p = rho + drho*0.5*(psi( i,j,k)+psi( i+1,j,k))
          rhoy_p = rho + drho*0.5*(psi( i,j,k)+psi( i,j+1,k))
          rhoz_p = rho + drho*0.5*(psi( i,j,k)+psi( i,j,k+1))
          u(i,j,k) = u(i,j,k)*rhox_n/rhox_p + (factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k))/rhox_p
          v(i,j,k) = v(i,j,k)*rhoy_n/rhoy_p + (factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k))/rhoy_p
          w(i,j,k) = w(i,j,k)*rhoz_n/rhoz_p + (factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k))/rhoz_p
        end do
      end do
    end do
    !
    ! swap d?dtrk <-> d?dtrko
    !
    call swap(dudtrk,dudtrko)
    call swap(dvdtrk,dvdtrko)
    call swap(dwdtrk,dwdtrko)
    !
    ! pressure, surface tension, and buoyancy terms
    !
    call mom_xyz_oth(n,dli,dzci,dzfi,dt_r,rho12,beta12,bforce,gacc,sigma,rho0,rho_av, &
                     p,psi,kappa,s,pn,po,dudtrk,dvdtrk,dwdtrk)
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          rhox_p = rho + drho*0.5*(psi(i,j,k)+psi(i+1,j,k))
          rhoy_p = rho + drho*0.5*(psi(i,j,k)+psi(i,j+1,k))
          rhoz_p = rho + drho*0.5*(psi(i,j,k)+psi(i,j,k+1))
          u(i,j,k) = u(i,j,k) + factor12*dudtrk(i,j,k)/rhox_p
          v(i,j,k) = v(i,j,k) + factor12*dvdtrk(i,j,k)/rhoy_p
          w(i,j,k) = w(i,j,k) + factor12*dwdtrk(i,j,k)/rhoz_p
        end do
      end do
    end do
  end subroutine rk
  !
  subroutine rk_scal(rkpar,n,dli,dzci,dzfi,dt,ssource,rho12,ka12,cp12,psi,u,v,w,psio,psiflx_x,psiflx_y,psiflx_z,s)
    use mod_scal, only: scal_ad
    !
    ! RK3 scheme for time integration of the scalar field
    !
    implicit none
    real(rp), intent(in   ), dimension(2) :: rkpar
    integer , intent(in   ), dimension(3) :: n
    real(rp), intent(in   ), dimension(3) :: dli
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi
    real(rp), intent(in   )                :: dt
    real(rp), intent(in   )                :: ssource
    real(rp), intent(in   ), dimension(2)  :: rho12,ka12,cp12
    real(rp), intent(in   ), dimension(0:,0:,0:) :: psi
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(in   ), dimension(0:,0:,0:) :: psio,psiflx_x,psiflx_y,psiflx_z
    real(rp), intent(inout), dimension(0:,0:,0:) :: s
    real(rp), target     , allocatable, dimension(:,:,:), save :: dsdtrk_t,dsdtrko_t
    real(rp), pointer    , contiguous , dimension(:,:,:), save :: dsdtrk  ,dsdtrko
    real(rp), dimension(2), save :: rhocp12
    real(rp) :: rhocp,drhocp,rhocp_n,rhocp_p
    logical, save :: is_first = .true.
    real(rp) :: factor1,factor2,factor12
    integer :: i,j,k
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    if(is_first) then ! leverage save attribute to allocate these arrays on the device only once
      allocate(dsdtrk_t(n(1),n(2),n(3)),dsdtrko_t(n(1),n(2),n(3)))
      !$acc enter data create(dsdtrk_t,dsdtrko_t) async(1)
      dsdtrk  => dsdtrk_t
      dsdtrko => dsdtrko_t
      rhocp12(:) = rho12(:)*cp12(:)
      !$acc wait(1)
    end if
    call scal_ad(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,ka12,rhocp12,psiflx_x,psiflx_y,psiflx_z,psi,u,v,w,s,dsdtrk)
    if(is_first) then
      !$acc kernels default(present) async(1)
      dsdtrko(:,:,:) = dsdtrk(:,:,:)
      !$acc end kernels
      is_first = .false.
    end if
    rhocp = rhocp12(2); drhocp = rhocp12(1)-rhocp12(2)
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          rhocp_p = rhocp + drhocp*psi( i,j,k)
          rhocp_n = rhocp + drhocp*psio(i,j,k)
          s(i,j,k) = s(i,j,k)*rhocp_n/rhocp_p + (factor1*dsdtrk(i,j,k) + factor2*dsdtrko(i,j,k) + factor12*ssource)/rhocp_p
        end do
      end do
    end do
    !
    ! swap d?dtrk <-> d?dtrko
    !
    call swap(dsdtrk,dsdtrko)
  end subroutine rk_scal
  !
  subroutine rk_2fl(rkpar,n,dli,dzci,dzfi,dt,gam,seps,beta,u,v,w,normx,normy,normz,phi,psi,psiflx_x,psiflx_y,psiflx_z)
#if !defined(_INTERFACE_CAPTURING_VOF)
    use mod_acdi        , only: acdi_transport_pf
#else
    use mod_vof_thinc_qq, only: vof_thinc_transport_psi
#endif
    use mod_two_fluid   , only: clip_field
    !
    ! RK3 scheme for time integration of the phase field
    !
    implicit none
    real(rp), intent(in   ), dimension(2) :: rkpar
    integer , intent(in   ), dimension(3) :: n
    real(rp), intent(in   ), dimension(3) :: dli
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi
    real(rp), intent(in   ) :: gam,seps,beta,dt
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(in   ), dimension(0:,0:,0:) :: normx,normy,normz
    real(rp), intent(in   ), dimension(0:,0:,0:), optional :: phi
    real(rp), intent(inout), dimension(0:,0:,0:) :: psi
    real(rp), intent(out  ), dimension(0:,0:,0:) :: psiflx_x,psiflx_y,psiflx_z
    real(rp), target     , allocatable, dimension(:,:,:), save :: dpsidtrk_t,dpsidtrko_t
    real(rp), pointer    , contiguous , dimension(:,:,:), save :: dpsidtrk  ,dpsidtrko
    logical, save :: is_first = .true.
    real(rp) :: factor1,factor2,factor12
    integer :: i,j,k
    !
    factor1  = rkpar(1)*dt
    factor2  = rkpar(2)*dt
    factor12 = factor1 + factor2 ! unused but left here in case a source term is needed in the future
    if(is_first) then ! leverage save attribute to allocate these arrays on the device only once
      allocate(dpsidtrk_t(n(1),n(2),n(3)),dpsidtrko_t(n(1),n(2),n(3)))
      !$acc enter data create(dpsidtrk_t,dpsidtrko_t) async(1)
      dpsidtrk  => dpsidtrk_t
      dpsidtrko => dpsidtrko_t
    end if
#if !defined(_INTERFACE_CAPTURING_VOF)
    call acdi_transport_pf(n,dli,dzci,dzfi,gam,seps,u,v,w,normx,normy,normz,phi,psi,dpsidtrk,psiflx_x,psiflx_y,psiflx_z)
#else
    call vof_thinc_transport_psi(n,dli,dzfi,beta,u,v,w,normx,normy,normz,psi,dpsidtrk,psiflx_x,psiflx_y,psiflx_z)
#endif
    if(is_first) then
      !$acc kernels default(present) async(1)
      dpsidtrko(:,:,:) = dpsidtrk(:,:,:)
      !$acc end kernels
      is_first = .false.
    end if
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          psi(i,j,k) = psi(i,j,k) + factor1*dpsidtrk(i,j,k) + factor2*dpsidtrko(i,j,k)
        end do
      end do
    end do
    call clip_field([1,1,1],[0._rp,1._rp],psi)
    !
    ! swap d?dtrk <-> d?dtrko
    !
    call swap(dpsidtrk,dpsidtrko)
    !
  end subroutine rk_2fl
end module mod_rk
