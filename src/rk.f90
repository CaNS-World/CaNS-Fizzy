! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
#define _FAST_MOM_KERNELS
module mod_rk
  use mod_mom  , only: mom_xyz_all
  use mod_utils, only: swap
  use mod_types
  implicit none
  private
  public rk
  contains
  subroutine rk(rkpar,n,dli,dzci,dzfi,dt, &
                bforce,gacc,sigma,rho_av,rho12,mu12,beta12,rho0,psi,normx,normy,normz,kappa,s, &
                p,pp,u,v,w)
    !
    ! low-storage 3rd-order Runge-Kutta scheme
    ! for time integration of the momentum equations.
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
    real(rp), intent(in   ), dimension(0:,0:,0:) :: psi,normx,normy,normz,kappa,s,p,pp
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    real(rp), target     , allocatable, dimension(:,:,:), save :: dudtrk_t ,dvdtrk_t ,dwdtrk_t , &
                                                                  dudtrko_t,dvdtrko_t,dwdtrko_t
    real(rp), pointer    , contiguous , dimension(:,:,:), save :: dudtrk   ,dvdtrk   ,dwdtrk   , &
                                                                  dudtrko  ,dvdtrko  ,dwdtrko
    logical, save :: is_first = .true.
    real(rp) :: factor1,factor2,factor12
    integer :: i,j,k
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    !
    ! initialization
    !
    if(is_first) then ! leverage save attribute to allocate these arrays on the device only once
      is_first = .false.
      allocate(dudtrk_t( n(1),n(2),n(3)),dvdtrk_t( n(1),n(2),n(3)),dwdtrk_t( n(1),n(2),n(3)))
      allocate(dudtrko_t(n(1),n(2),n(3)),dvdtrko_t(n(1),n(2),n(3)),dwdtrko_t(n(1),n(2),n(3)))
      !$acc enter data create(dudtrk_t ,dvdtrk_t ,dwdtrk_t ) async(1)
      !$acc enter data create(dudtrko_t,dvdtrko_t,dwdtrko_t) async(1)
      !$acc kernels default(present) async(1) ! not really necessary
      dudtrko_t(:,:,:) = dudtrk_t(:,:,:)
      dvdtrko_t(:,:,:) = dvdtrk_t(:,:,:)
      dwdtrko_t(:,:,:) = dwdtrk_t(:,:,:)
      !$acc end kernels
      dudtrk  => dudtrk_t
      dvdtrk  => dvdtrk_t
      dwdtrk  => dwdtrk_t
      dudtrko => dudtrko_t
      dvdtrko => dvdtrko_t
      dwdtrko => dwdtrko_t
    end if
    !
    call mom_xyz_all(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,rho12,mu12,beta12,bforce,gacc,sigma,rho0,rho_av, &
                     u,v,w,p,pp,psi,kappa,s,dudtrk,dvdtrk,dwdtrk)
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
  end subroutine rk
  !
  subroutine rk_scal(rkpar,n,dli,dzci,dzfi,dt, &
                     ssource,rho12,ka12,cp12,psi,u,v,w,s)
  use mod_scal, only: scal_ad
    !
    ! low-storage 3rd-order Runge-Kutta scheme
    ! for time integration of the scalar field.
    !
    ! n.b.: since we leverage the `save` attribute for dsdtrk*, this subroutine only supports
    !       transport of a single scalar; extension to n arbtritrary scalars could be done,
    !       e.g., using a loop through an array of scalar derived types, with an ASSOCIATE
    !       statement to keep the same piece of code to for each scalar
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
    if(is_first) then ! leverage save attribute to allocate these arrays on the device only once
      is_first = .false.
      allocate(dsdtrk_t(n(1),n(2),n(3)),dsdtrko_t(n(1),n(2),n(3)))
      !$acc enter data create(dsdtrk_t,dsdtrko_t) async(1)
      !$acc kernels default(present) async(1) ! not really necessary
      dsdtrko_t(:,:,:) = 0._rp
      !$acc end kernels
      dsdtrk  => dsdtrk_t
      dsdtrko => dsdtrko_t
    end if
    call scal_ad(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,ssource,ka12,rhocp12,psi,u,v,w,s,dsdtrk)
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
end module mod_rk
