! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
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
  subroutine rk(rkpar,n,dli,dzci,dzfi,dt, &
                bforce,gacc,sigma,rho_av,rho12,mu12,beta12,rho0,psi,kappa,s,p,pp, &
                acdi_rgx,acdi_rgy,acdi_rgz,u,v,w)
    !
    ! Adams-Bashforth scheme for time integration of the momentum equations
    !
    implicit none
    real(rp), intent(in), dimension(2) :: rkpar
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: dli
    real(rp), intent(in) :: dt
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi
    real(rp), intent(in   ), dimension(3)  :: bforce,gacc
    real(rp), intent(in   )                :: sigma,rho_av
    real(rp), intent(in   ), dimension(2)  :: rho12,mu12,beta12
    real(rp), intent(in   )                :: rho0
    real(rp), intent(in   ), dimension(0:,0:,0:) :: psi,kappa,s,p,pp
    real(rp), intent(in   ), dimension(0:,0:,0:) :: acdi_rgx,acdi_rgy,acdi_rgz
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    real(rp), target       , allocatable, dimension(:,:,:), save :: dudtrk_t ,dvdtrk_t ,dwdtrk_t , &
                                                                    dudtrko_t,dvdtrko_t,dwdtrko_t
    real(rp), pointer      , contiguous , dimension(:,:,:), save :: dudtrk   ,dvdtrk   ,dwdtrk   , &
                                                                    dudtrko  ,dvdtrko  ,dwdtrko
    logical, save :: is_first = .true.
    real(rp) :: factor1,factor2,factor12,dt_r
    integer :: i,j,k
    !
    factor1  = rkpar(1)*dt
    factor2  = rkpar(2)*dt
    factor12 = factor1+factor2
    dt_r     = -2.*rkpar(2) ! N.B.: For AB2 only, rkpar(2) = -(dt/dto)/2 = dt_r/2
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
    call mom_xyz_ad(n,dli,dzci,dzfi,rho12,mu12,acdi_rgx,acdi_rgy,acdi_rgz,u,v,w,psi,dudtrk,dvdtrk,dwdtrk)
    !
    if(is_first) then ! use Euler forward
      !$acc kernels default(present) async(1)
      dudtrko(:,:,:) = dudtrk(:,:,:)
      dvdtrko(:,:,:) = dvdtrk(:,:,:)
      dwdtrko(:,:,:) = dwdtrk(:,:,:)
      !$acc end kernels
      is_first = .false.
    end if
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          u(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k)
          v(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k)
          w(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k)
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
                     p,pp,psi,kappa,s,dudtrk,dvdtrk,dwdtrk)
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          u(i,j,k) = u(i,j,k) + factor12*dudtrk(i,j,k)
          v(i,j,k) = v(i,j,k) + factor12*dvdtrk(i,j,k)
          w(i,j,k) = w(i,j,k) + factor12*dwdtrk(i,j,k)
        end do
      end do
    end do
  end subroutine rk
  !
  subroutine rk_scal(rkpar,n,dli,dzci,dzfi,dt, &
                     ssource,rho12,ka12,cp12,psi,u,v,w,s)
    use mod_scal, only: scal_ad
    !
    ! Adams-Bashfroth scheme for time integration of the scalar field
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
    real(rp), intent(inout), dimension(0:,0:,0:) :: s
    real(rp), target     , allocatable, dimension(:,:,:), save :: dsdtrk_t,dsdtrko_t
    real(rp), pointer    , contiguous , dimension(:,:,:), save :: dsdtrk  ,dsdtrko
    logical, save :: is_first = .true.
    real(rp) :: factor1,factor2,factor12
    integer :: i,j,k
    real(rp), dimension(2) :: rhocp12
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    rhocp12 = rho12(:)*cp12(:)
    if(is_first) then ! leverage save attribute to allocate these arrays on the device only once
      allocate(dsdtrk_t(n(1),n(2),n(3)),dsdtrko_t(n(1),n(2),n(3)))
      !$acc enter data create(dsdtrk_t,dsdtrko_t) async(1)
      dsdtrk  => dsdtrk_t
      dsdtrko => dsdtrko_t
    end if
    call scal_ad(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,ssource,ka12,rhocp12,psi,u,v,w,s,dsdtrk)
    if(is_first) then ! use Euler forward
      !$acc kernels default(present) async(1)
      dsdtrko(:,:,:) = dsdtrk(:,:,:)
      !$acc end kernels
      is_first = .false.
    end if
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          s(i,j,k) = s(i,j,k) + factor1*dsdtrk(i,j,k) + factor2*dsdtrko(i,j,k) + factor12*ssource
        end do
      end do
    end do
    !
    ! swap d?dtrk <-> d?dtrko
    !
    call swap(dsdtrk,dsdtrko)
  end subroutine rk_scal
  !
  subroutine rk_2fl(rkpar,n,dli,dzci,dzfi,dt,gam,seps,u,v,w,normx,normy,normz,psi,rglrx,rglry,rglrz)
    use mod_acdi     , only: acdi_transport_pf
    use mod_two_fluid, only: clip_field
    !
    ! Adams-Bashforth scheme for time integration of the phase field
    !
    implicit none
    logical , parameter :: is_cmpt_wallflux = .false.
    real(rp), intent(in   ), dimension(2) :: rkpar
    integer , intent(in   ), dimension(3) :: n
    real(rp), intent(in   ), dimension(3) :: dli
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi
    real(rp), intent(in   ) :: gam,seps,dt
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(in   ), dimension(0:,0:,0:) :: normx,normy,normz
    real(rp), intent(inout), dimension(0:,0:,0:) :: psi
    real(rp), intent(out  ), dimension(0:,0:,0:) :: rglrx,rglry,rglrz
    real(rp), target     , allocatable, dimension(:,:,:), save :: dpsidtrk_t,dpsidtrko_t
    real(rp), pointer    , contiguous , dimension(:,:,:), save :: dpsidtrk  ,dpsidtrko
    logical, save :: is_first = .true.
    real(rp) :: factor1,factor2,factor12
    integer :: i,j,k
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    if(is_first) then ! leverage save attribute to allocate these arrays on the device only once
      allocate(dpsidtrk_t(n(1),n(2),n(3)),dpsidtrko_t(n(1),n(2),n(3)))
      !$acc enter data create(dpsidtrk_t,dpsidtrko_t) async(1)
      dpsidtrk  => dpsidtrk_t
      dpsidtrko => dpsidtrko_t
    end if
    call acdi_transport_pf(n,dli,dzci,dzfi,gam,seps,u,v,w,normx,normy,normz,psi,dpsidtrk,rglrx,rglry,rglrz)
    if(is_first) then ! use Euler forward
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
