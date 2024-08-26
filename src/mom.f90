! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_mom
  use mpi
  use mod_param, only: nh
  use mod_types
  implicit none
  private
  public mom_xyz_ad,mom_xyz_oth
  contains
  !
  subroutine momx_a(nx,ny,nz,dxi,dyi,dzfi,u,v,w,dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(1-nh:) :: dzfi
    real(rp), dimension(1-nh:,1-nh:,1-nh:), intent(in   ) :: u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    integer :: i,j,k
    real(rp) :: uuip,uuim,vujp,vujm,wukp,wukm
    !
    !$acc parallel loop collapse(3) default(present) private(uuip,uuim,vujp,vujm,wukp,wukm) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          uuip  = 0.25*( u(i+1,j  ,k  )+u(i,j  ,k  ) )*( u(i,j,k)+u(i+1,j,k) )
          uuim  = 0.25*( u(i-1,j  ,k  )+u(i,j  ,k  ) )*( u(i,j,k)+u(i-1,j,k) )
          vujp  = 0.25*( v(i+1,j  ,k  )+v(i,j  ,k  ) )*( u(i,j,k)+u(i,j+1,k) )
          vujm  = 0.25*( v(i+1,j-1,k  )+v(i,j-1,k  ) )*( u(i,j,k)+u(i,j-1,k) )
          wukp  = 0.25*( w(i+1,j  ,k  )+w(i,j  ,k  ) )*( u(i,j,k)+u(i,j,k+1) )
          wukm  = 0.25*( w(i+1,j  ,k-1)+w(i,j  ,k-1) )*( u(i,j,k)+u(i,j,k-1) )
          !
          dudt(i,j,k) = dudt(i,j,k) + &
                        dxi*(     -uuip + uuim ) + &
                        dyi*(     -vujp + vujm ) + &
                        dzfi(k)*( -wukp + wukm )
        end do
      end do
    end do
  end subroutine momx_a
  !
  subroutine momy_a(nx,ny,nz,dxi,dyi,dzfi,u,v,w,dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(1-nh:) :: dzfi
    real(rp), dimension(1-nh:,1-nh:,1-nh:), intent(in   ) :: u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    integer :: i,j,k
    real(rp) :: uvip,uvim,vvjp,vvjm,wvkp,wvkm
    !
    !$acc parallel loop collapse(3) default(present) private(uvip,uvim,vvjp,vvjm,wvkp,wvkm) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          uvip  = 0.25*( u(i  ,j,k  )+u(i  ,j+1,k  ) )*( v(i,j,k)+v(i+1,j,k) )
          uvim  = 0.25*( u(i-1,j,k  )+u(i-1,j+1,k  ) )*( v(i,j,k)+v(i-1,j,k) )
          vvjp  = 0.25*( v(i  ,j,k  )+v(i  ,j+1,k  ) )*( v(i,j,k)+v(i,j+1,k) )
          vvjm  = 0.25*( v(i  ,j,k  )+v(i  ,j-1,k  ) )*( v(i,j,k)+v(i,j-1,k) )
          wvkp  = 0.25*( w(i  ,j,k  )+w(i  ,j+1,k  ) )*( v(i,j,k)+v(i,j,k+1) )
          wvkm  = 0.25*( w(i  ,j,k-1)+w(i  ,j+1,k-1) )*( v(i,j,k)+v(i,j,k-1) )
          !
          dvdt(i,j,k) = dvdt(i,j,k) + &
                        dxi*(     -uvip + uvim ) + &
                        dyi*(     -vvjp + vvjm ) + &
                        dzfi(k)*( -wvkp + wvkm )
        end do
      end do
    end do
  end subroutine momy_a
  !
  subroutine momz_a(nx,ny,nz,dxi,dyi,dzci,u,v,w,dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(1-nh:) :: dzci
    real(rp), dimension(1-nh:,1-nh:,1-nh:), intent(in   ) :: u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    integer :: i,j,k
    real(rp) :: uwip,uwim,vwjp,vwjm,wwkp,wwkm
    !
    !$acc parallel loop collapse(3) default(present) private(uwip,uwim,vwjp,vwjm,wwkp,wwkm) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          uwip  = 0.25*( u(i  ,j  ,k)+u(i  ,j  ,k+1) )*( w(i,j,k)+w(i+1,j,k) )
          uwim  = 0.25*( u(i-1,j  ,k)+u(i-1,j  ,k+1) )*( w(i,j,k)+w(i-1,j,k) )
          vwjp  = 0.25*( v(i  ,j  ,k)+v(i  ,j  ,k+1) )*( w(i,j,k)+w(i,j+1,k) )
          vwjm  = 0.25*( v(i  ,j-1,k)+v(i  ,j-1,k+1) )*( w(i,j,k)+w(i,j-1,k) )
          wwkp  = 0.25*( w(i  ,j  ,k)+w(i  ,j  ,k+1) )*( w(i,j,k)+w(i,j,k+1) )
          wwkm  = 0.25*( w(i  ,j  ,k)+w(i  ,j  ,k-1) )*( w(i,j,k)+w(i,j,k-1) )
          !
          dwdt(i,j,k) = dwdt(i,j,k) + &
                        dxi*(     -uwip + uwim ) + &
                        dyi*(     -vwjp + vwjm ) + &
                        dzci(k)*( -wwkp + wwkm )
        end do
      end do
    end do
  end subroutine momz_a
  !
  subroutine momx_d(nx,ny,nz,dxi,dyi,dzci,dzfi,rho12,mu12,psi,u,v,w,dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(1-nh:) :: dzci,dzfi
    real(rp), intent(in), dimension(2) :: rho12,mu12
    real(rp), dimension(1-nh:,1-nh:,1-nh:), intent(in   ) :: psi,u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm, &
                dvdxp,dvdxm,dwdxp,dwdxm
    real(rp) :: muxp,muxm,muyp,muym,muzp,muzm,rhop
    integer :: i,j,k
    real(rp) :: rho,drho,mu,dmu
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    mu  = mu12(2) ; dmu  = mu12(1)-mu12(2)
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(dudxp,dudxm,dudyp,dudym,dudzp,dudzm,dvdxp,dvdxm,dwdxp,dwdxm) &
    !$acc private(muxp,muxm,muyp,muym,muzp,muzm,rhop) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dudxp = (u(i+1,j  ,k  )-u(i  ,j  ,k  ))*dxi
          dudxm = (u(i  ,j  ,k  )-u(i-1,j  ,k  ))*dxi
          dvdxp = (v(i+1,j  ,k  )-v(i  ,j  ,k  ))*dxi
          dvdxm = (v(i+1,j-1,k  )-v(i  ,j-1,k  ))*dxi
          dudyp = (u(i  ,j+1,k  )-u(i  ,j  ,k  ))*dyi
          dudym = (u(i  ,j  ,k  )-u(i  ,j-1,k  ))*dyi
          dudzp = (u(i  ,j  ,k+1)-u(i  ,j  ,k  ))*dzci(k  )
          dudzm = (u(i  ,j  ,k  )-u(i  ,j  ,k-1))*dzci(k-1)
          dwdxp = (w(i+1,j  ,k  )-w(i  ,j  ,k  ))*dxi
          dwdxm = (w(i+1,j  ,k-1)-w(i  ,j  ,k-1))*dxi
          !
          muxp = mu + dmu*psi(i+1,j,k)
          muxm = mu + dmu*psi(i  ,j,k)
          muyp = mu + dmu*0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i+1,j+1,k)+psi(i+1,j,k))
          muym = mu + dmu*0.25*(psi(i,j,k)+psi(i,j-1,k)+psi(i+1,j-1,k)+psi(i+1,j,k))
          muzp = mu + dmu*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i+1,j,k+1)+psi(i+1,j,k))
          muzm = mu + dmu*0.25*(psi(i,j,k)+psi(i,j,k-1)+psi(i+1,j,k-1)+psi(i+1,j,k))
          rhop = rho + drho*0.5*(psi(i+1,j,k)+psi(i,j,k))
          !
          dudt(i,j,k) = dudt(i,j,k) + &
                        dxi*(    (dudxp+dudxp)*muxp-(dudxm+dudxm)*muxm)/rhop + &
                        dyi*(    (dudyp+dvdxp)*muyp-(dudym+dvdxm)*muym)/rhop + &
                        dzfi(k)*((dudzp+dwdxp)*muzp-(dudzm+dwdxm)*muzm)/rhop
        end do
      end do
    end do
  end subroutine momx_d
  !
  subroutine momy_d(nx,ny,nz,dxi,dyi,dzci,dzfi,rho12,mu12,psi,u,v,w,dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(1-nh:) :: dzci,dzfi
    real(rp), intent(in), dimension(2) :: rho12,mu12
    real(rp), dimension(1-nh:,1-nh:,1-nh:), intent(in   ) :: psi,u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    real(rp) :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm, &
                dudyp,dudym,dwdyp,dwdym
    real(rp) :: muxp,muxm,muyp,muym,muzp,muzm,rhop
    integer :: i,j,k
    real(rp) :: rho,drho,mu,dmu
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    mu  = mu12(2) ; dmu  = mu12(1)-mu12(2)
    !
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm,dudyp,dudym,dwdyp,dwdym) &
    !$acc private(muxp,muxm,muyp,muym,muzp,muzm,rhop) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dvdxp = (v(i+1,j  ,k  )-v(i  ,j  ,k  ))*dxi
          dvdxm = (v(i  ,j  ,k  )-v(i-1,j  ,k  ))*dxi
          dudyp = (u(i  ,j+1,k  )-u(i  ,j  ,k  ))*dyi
          dudym = (u(i-1,j+1,k  )-u(i-1,j  ,k  ))*dyi
          dvdyp = (v(i  ,j+1,k  )-v(i  ,j  ,k  ))*dyi
          dvdym = (v(i  ,j  ,k  )-v(i  ,j-1,k  ))*dyi
          dvdzp = (v(i  ,j  ,k+1)-v(i  ,j  ,k  ))*dzci(k  )
          dvdzm = (v(i  ,j  ,k  )-v(i  ,j  ,k-1))*dzci(k-1)
          dwdyp = (w(i  ,j+1,k  )-w(i  ,j  ,k  ))*dyi
          dwdym = (w(i  ,j+1,k-1)-w(i  ,j  ,k-1))*dyi
          !
          muxp = mu + dmu*0.25*(psi(i,j,k)+psi(i+1,j,k)+psi(i+1,j+1,k)+psi(i,j+1,k))
          muxm = mu + dmu*0.25*(psi(i,j,k)+psi(i-1,j,k)+psi(i-1,j+1,k)+psi(i,j+1,k))
          muyp = mu + dmu*psi(i,j+1,k)
          muym = mu + dmu*psi(i,j  ,k)
          muzp = mu + dmu*0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i,j+1,k+1)+psi(i,j,k+1))
          muzm = mu + dmu*0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i,j+1,k-1)+psi(i,j,k-1))
          rhop = rho + drho*0.5*(psi(i,j+1,k)+psi(i,j,k))
          !
          dvdt(i,j,k) = dvdt(i,j,k) + &
                        dxi*(    (dvdxp+dudyp)*muxp-(dvdxm+dudym)*muxm)/rhop + &
                        dyi*(    (dvdyp+dvdyp)*muyp-(dvdym+dvdym)*muym)/rhop + &
                        dzfi(k)*((dvdzp+dwdyp)*muzp-(dvdzm+dwdym)*muzm)/rhop
        end do
      end do
    end do
  end subroutine momy_d
  !
  subroutine momz_d(nx,ny,nz,dxi,dyi,dzci,dzfi,rho12,mu12,psi,u,v,w,dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(1-nh:) :: dzci,dzfi
    real(rp), intent(in), dimension(2) :: rho12,mu12
    real(rp), dimension(1-nh:,1-nh:,1-nh:), intent(in   ) :: psi,u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    integer :: i,j,k
    real(rp) :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm, &
                dudzp,dudzm,dvdzp,dvdzm
    real(rp) :: muxp,muxm,muyp,muym,muzp,muzm,rhop
    real(rp) :: rho,drho,mu,dmu
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    mu  = mu12(2) ; dmu  = mu12(1)-mu12(2)
    !
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm,dudzp,dudzm,dvdzp,dvdzm) &
    !$acc private(muxp,muxm,muyp,muym,muzp,muzm,rhop) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dwdxp = (w(i+1,j  ,k  )-w(i  ,j  ,k))*dxi
          dwdxm = (w(i  ,j  ,k  )-w(i-1,j  ,k))*dxi
          dudzp = (u(i  ,j  ,k+1)-u(i  ,j  ,k))*dzci(k  )
          dudzm = (u(i-1,j  ,k+1)-u(i-1,j  ,k))*dzci(k  )
          dwdyp = (w(i  ,j+1,k  )-w(i  ,j  ,k))*dyi
          dwdym = (w(i  ,j  ,k  )-w(i  ,j-1,k))*dyi
          dvdzp = (v(i  ,j  ,k+1)-v(i  ,j  ,k))*dzci(k  )
          dvdzm = (v(i  ,j-1,k+1)-v(i  ,j-1,k))*dzci(k  )
          dwdzp = (w(i  ,j  ,k+1)-w(i  ,j,k  ))*dzfi(k+1)
          dwdzm = (w(i  ,j  ,k  )-w(i  ,j,k-1))*dzfi(k  )
          !
          muxp = mu + dmu*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i+1,j,k+1)+psi(i+1,j,k) )
          muxm = mu + dmu*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i-1,j,k+1)+psi(i-1,j,k) )
          muyp = mu + dmu*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i,j+1,k+1)+psi(i,j+1,k) )
          muym = mu + dmu*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i,j-1,k+1)+psi(i,j-1,k) )
          muzp = mu + dmu*psi(i,j,k+1)
          muzm = mu + dmu*psi(i,j,k  )
          rhop = rho + drho*0.5*(psi(i,j,k+1)+psi(i,j,k))
          !
          dwdt(i,j,k) = dwdt(i,j,k) + &
                        dxi*(    (dwdxp+dudzp)*muxp-(dwdxm+dudzm)*muxm)/rhop + &
                        dyi*(    (dwdyp+dvdzp)*muyp-(dwdym+dvdzm)*muym)/rhop + &
                        dzci(k)*((dwdzp+dwdzp)*muzp-(dwdzm+dwdzm)*muzm)/rhop
        end do
      end do
    end do
  end subroutine momz_d
  !
  subroutine momx_p(nx,ny,nz,dxi,dt_r,bforce,gacc,rho0,rho_av,rho12,psi,p,pp,dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dt_r
    real(rp), intent(in) :: bforce,gacc,rho0,rho_av,rho12(2)
    real(rp), dimension(1-nh:,1-nh:,1-nh:), intent(in   ) :: psi,p,pp
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    real(rp) :: rhop,dpdl
    integer :: i,j,k
    real(rp) :: rho,drho
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    !
    !$acc parallel loop collapse(3) default(present) private(rhop,dpdl) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          rhop = rho + drho*0.5*(psi(i+1,j,k)+psi(i,j,k))
          dpdl = (p(i+1,j,k)-p(i,j,k))*dxi
          !
          dudt(i,j,k) = dudt(i,j,k) + bforce/rhop + gacc*(1.-rho_av/rhop) &
#if defined(_CONSTANT_COEFFS_POISSON)
                        -dpdl/rho0 - (1./rhop-1./rho0)*( ((1.+dt_r)*p(i+1,j,k)-dt_r*(pp(i+1,j,k))) - &
                                                         ((1.+dt_r)*p(i  ,j,k)-dt_r*(pp(i  ,j,k))) &
                                                        )*dxi
#else
                        -dpdl/rhop
#endif
        end do
      end do
    end do
  end subroutine momx_p
  !
  subroutine momy_p(nx,ny,nz,dyi,dt_r,bforce,gacc,rho0,rho_av,rho12,psi,p,pp,dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dyi,dt_r
    real(rp), intent(in) :: bforce,gacc,rho0,rho_av,rho12(2)
    real(rp), dimension(1-nh:,1-nh:,1-nh:), intent(in   ) :: psi,p,pp
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    integer :: i,j,k
    real(rp) :: rhop,dpdl
    real(rp) :: rho,drho
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    !
    !$acc parallel loop collapse(3) default(present) private(rhop,dpdl) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          rhop = rho + drho*0.5*(psi(i,j+1,k)+psi(i,j,k))
          dpdl = (p(i,j+1,k)-p(i,j,k))*dyi
          !
          dvdt(i,j,k) = dvdt(i,j,k) + bforce/rhop + gacc*(1.-rho_av/rhop) &
#if defined(_CONSTANT_COEFFS_POISSON)
                        -dpdl/rho0 - (1./rhop-1./rho0)*( ((1.+dt_r)*p(i,j+1,k)-dt_r*(pp(i,j+1,k))) - &
                                                         ((1.+dt_r)*p(i,j  ,k)-dt_r*(pp(i,j  ,k))) &
                                                       )*dyi
#else
                        -dpdl/rhop
#endif
        end do
      end do
    end do
  end subroutine momy_p
  !
  subroutine momz_p(nx,ny,nz,dzci,dt_r,bforce,gacc,rho0,rho_av,rho12,psi,p,pp,dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in), dimension(1-nh:) :: dzci
    real(rp), intent(in) :: dt_r
    real(rp), intent(in) :: bforce,gacc,rho0,rho_av,rho12(2)
    real(rp), dimension(1-nh:,1-nh:,1-nh:), intent(in   ) :: psi,p,pp
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    real(rp) :: rhop,dpdl
    integer :: i,j,k
    real(rp) :: rho,drho
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    !
    !$acc parallel loop collapse(3) default(present) private(rhop,dpdl) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          rhop = rho + drho*0.5*(psi(i,j,k+1)+psi(i,j,k))
          dpdl = (p(i,j,k+1)-p(i,j,k))*dzci(k)
          !
          dwdt(i,j,k) = dwdt(i,j,k) + bforce/rhop + gacc*(1.-rho_av/rhop) &
#if defined(_CONSTANT_COEFFS_POISSON)
                        -dpdl/rho0 - (1./rhop-1./rho0)*( ((1.+dt_r)*p(i,j,k+1)-dt_r*(pp(i,j,k+1))) - &
                                                         ((1.+dt_r)*p(i,j,k  )-dt_r*(pp(i,j,k  ))) &
                                                       )*dzci(k)
#else
                        -dpdl/rhop
#endif
        end do
      end do
    end do
  end subroutine momz_p
  !
  subroutine momx_sigma(nx,ny,nz,dxi,sigma,rho12,kappa,psi,dudt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi
    real(rp), intent(in) :: sigma,rho12(2)
    real(rp), dimension(1-nh:,1-nh:,1-nh:), intent(in   ) :: psi,kappa
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    real(rp) :: rhop,kappap
    integer :: i,j,k
    real(rp) :: rho,drho
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    !
    !$acc parallel loop collapse(3) default(present) private(rhop,kappap) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          rhop   = rho + drho*0.5*(psi(i+1,j,k)+psi(i,j,k))
          kappap = 0.5*(kappa(i+1,j,k)+kappa(i,j,k))
          !
          dudt(i,j,k) = dudt(i,j,k) + &
                        sigma*kappap*(psi(i+1,j,k)-psi(i,j,k))*dxi/rhop
        end do
      end do
    end do
  end subroutine momx_sigma
  !
  subroutine momy_sigma(nx,ny,nz,dyi,sigma,rho12,psi,kappa,dvdt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dyi
    real(rp), intent(in) :: sigma,rho12(2)
    real(rp), dimension(1-nh:,1-nh:,1-nh:), intent(in   ) :: psi,kappa
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    real(rp) :: rhop,kappap
    integer :: i,j,k
    real(rp) :: rho,drho
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    !
    !$acc parallel loop collapse(3) default(present) private(rhop,kappap) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          rhop   = rho + drho*0.5*(psi(i,j+1,k)+psi(i,j,k))
          kappap = 0.5*(kappa(i,j+1,k)+kappa(i,j,k))
          !
          dvdt(i,j,k) = dvdt(i,j,k) + &
                        sigma*kappap*(psi(i,j+1,k)-psi(i,j,k))*dyi/rhop
        end do
      end do
    end do
  end subroutine momy_sigma
  !
  subroutine momz_sigma(nx,ny,nz,dzci,sigma,rho12,psi,kappa,dwdt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in), dimension(1-nh:) :: dzci
    real(rp), intent(in) :: sigma,rho12(2)
    real(rp), dimension(1-nh:,1-nh:,1-nh:), intent(in   ) :: psi,kappa
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    real(rp) :: rhop,kappap
    integer :: i,j,k
    real(rp) :: rho,drho
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    !
    !$acc parallel loop collapse(3) default(present) private(rhop,kappap) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          rhop   = rho + drho*0.5*(psi(i,j,k+1)+psi(i,j,k))
          kappap = 0.5*(kappa(i,j,k+1)+kappa(i,j,k))
          !
          dwdt(i,j,k) = dwdt(i,j,k) + &
                        sigma*kappap*(psi(i,j,k+1)-psi(i,j,k))*dzci(k)/rhop
        end do
      end do
    end do
  end subroutine momz_sigma
  !
  subroutine momx_buoy(nx,ny,nz,gacc,rho12,beta12,psi,s,dudt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: gacc,rho12(2),beta12(2)
    real(rp), dimension(1-nh:,1-nh:,1-nh:), intent(in   ) :: psi,s
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    real(rp) :: psip,rhop,factorp
    integer :: i,j,k
    real(rp) :: rho,drho,rhobeta,drhobeta
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    rhobeta = rho12(2)*beta12(2); drhobeta = rho12(1)*beta12(1)-rho12(2)*beta12(2)
    !
    !$acc parallel loop collapse(3) default(present) private(psip,rhop,factorp) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          psip = 0.5*(psi(i+1,j,k)+psi(i,j,k))
          rhop = rho + drho*psip
          factorp = rhobeta + drhobeta*psip
          !
          dudt(i,j,k) = dudt(i,j,k) - &
                        gacc*factorp*0.5*(s(i+1,j,k)+s(i,j,k))/rhop
        end do
      end do
    end do
  end subroutine momx_buoy
  !
  subroutine momy_buoy(nx,ny,nz,gacc,rho12,beta12,psi,s,dvdt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: gacc,rho12(2),beta12(2)
    real(rp), dimension(1-nh:,1-nh:,1-nh:), intent(in   ) :: psi,s
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    real(rp) :: psip,rhop,factorp
    integer :: i,j,k
    real(rp) :: rho,drho,rhobeta,drhobeta
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    rhobeta = rho12(2)*beta12(2); drhobeta = rho12(1)*beta12(1)-rho12(2)*beta12(2)
    !
    !$acc parallel loop collapse(3) default(present) private(factorp) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          psip = 0.5*(psi(i,j+1,k)+psi(i,j,k))
          rhop = rho + drho*psip
          factorp = rhobeta + drhobeta*psip
          !
          dvdt(i,j,k) = dvdt(i,j,k) - &
                        gacc*factorp*0.5*(s(i,j+1,k)+s(i,j,k))/rhop
        end do
      end do
    end do
  end subroutine momy_buoy
  !
  subroutine momz_buoy(nx,ny,nz,gacc,rho12,beta12,psi,s,dwdt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: gacc,rho12(2),beta12(2)
    real(rp), dimension(1-nh:,1-nh:,1-nh:), intent(in   ) :: psi,s
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    real(rp) :: psip,rhop,factorp
    integer :: i,j,k
    real(rp) :: rho,drho,rhobeta,drhobeta
    !
    rho = rho12(2); drho = rho12(2)
    rhobeta = rho12(2)*beta12(2); drhobeta = rho12(1)*beta12(1)-rho12(2)*beta12(2)
    !
    !$acc parallel loop collapse(3) default(present) private(rhop,factorp) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          psip = 0.5*(psi(i,j,k+1)+psi(i,j,k))
          rhop = rho + drho*psip
          factorp = rhobeta + drhobeta*psip
          !
          dwdt(i,j,k) = dwdt(i,j,k) - &
                        gacc*factorp*0.5*(s(i,j,k+1)+s(i,j,k))/rhop
        end do
      end do
    end do
  end subroutine momz_buoy
  !
  subroutine mom_xyz_ad(n,dli,dzci,dzfi,rho12,mu12,acdi_rglrx,acdi_rglry,acdi_rglrz, &
                        u,v,w,psi,psio,dudt,dvdt,dwdt)
    implicit none
    integer , intent(in ), dimension(3) :: n
    real(rp), intent(in ), dimension(3) :: dli
    real(rp), intent(in ), dimension(1-nh:) :: dzci,dzfi
    real(rp), intent(in ), dimension(2) :: rho12,mu12
    real(rp), intent(in ), dimension(1-nh:,1-nh:,1-nh:) :: u,v,w
    real(rp), intent(in ), dimension(1-nh:,1-nh:,1-nh:) :: psi
    real(rp), intent(in ), dimension(1-nh:,1-nh:,1-nh:), optional :: psio,acdi_rglrx,acdi_rglry,acdi_rglrz
    real(rp), intent(out), dimension( :, :, :) :: dudt,dvdt,dwdt
    integer :: i,j,k
    real(rp) :: rho,drho,mu,dmu,rhobeta,drhobeta
    real(rp) :: dxi,dyi
    real(rp) :: c_ccm,c_pcm,c_cpm,c_cmc,c_pmc,c_mcc,c_ccc,c_pcc,c_mpc,c_cpc,c_cmp,c_mcp,c_ccp,c_cpp,c_ppc,c_pcp, &
                d_ccm,d_pcm,d_cpm,d_cmc,d_pmc,d_mcc,d_ccc,d_pcc,d_mpc,d_cpc,d_cmp,d_mcp,d_ccp,d_cpp,d_ppc,d_pcp, &
                u_ccl,u_ccm,u_pcm,u_cpm,u_clc,u_cmc,u_pmc,u_lcc,u_mcc,u_ccc,u_pcc,u_qcc,u_mpc,u_cpc,u_cqc,u_cmp,u_mcp,u_ccp,u_ccq, &
                v_ccl,v_ccm,v_pcm,v_cpm,v_clc,v_cmc,v_pmc,v_lcc,v_mcc,v_ccc,v_pcc,v_qcc,v_mpc,v_cpc,v_cqc,v_cmp,v_mcp,v_ccp,v_ccq, &
                w_ccl,w_ccm,w_pcm,w_cpm,w_clc,w_cmc,w_pmc,w_lcc,w_mcc,w_ccc,w_pcc,w_qcc,w_mpc,w_cpc,w_cqc,w_cmp,w_mcp,w_ccp,w_ccq, &
                rglrx_mcc,rglrx_ccc,rglrx_pcc,rglrx_mpc,rglrx_cpc,rglrx_mcp,rglrx_ccp, &
                rglry_cmc,rglry_ccc,rglry_cpc,rglry_pmc,rglry_pcc,rglry_cmp,rglry_ccp, &
                rglrz_ccm,rglrz_ccc,rglrz_ccp,rglrz_pcm,rglrz_pcc,rglrz_cpm,rglrz_cpc, &
                dzci_c,dzci_m,dzfi_c,dzfi_p, &
                psixp,psiyp,psizp
    real(rp) :: uuip,uuim,vujp,vujm,wukp,wukm,uvip,uvim,vvjp,vvjm,wvkp,wvkm,uwip,uwim,vwjp,vwjm,wwkp,wwkm
    real(rp) :: ududx,vdudy,wdudz,udvdx,vdvdy,wdvdz,udwdx,vdwdy,wdwdz
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm,dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm,dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    real(rp) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
    real(rp) :: omxp,omxm,omyp,omym,omzp,omzm,domzdy,domydz,domxdz,domzdx,domydx,domxdy
    real(rp) :: dpsidx,dpsidy,dpsidz,dpsidxp,dpsidxm,dpsidyp,dpsidym,dpsidzp,dpsidzm
    real(rp) :: mux,muy,muz,muxp,muxm,muyp,muym,muzp,muzm,dmudx,dmudy,dmudz
    real(rp) :: rxup,rxum,ryup,ryum,rzup,rzum,rxvp,rxvm,ryvp,ryvm,rzvp,rzvm,rxwp,rxwm,rywp,rywm,rzwp,rzwm
    real(rp) :: rhox,rhoy,rhoz,rhoxp,rhoxm,rhoyp,rhoym,rhozp,rhozm
    real(rp) :: dudt_aux,dvdt_aux,dwdt_aux
    !
    real(rp), parameter, dimension(2,2) :: c = 1.d0/2.d0*reshape((/-1.d0,3.d0,&
                                                                  1.d0, 1.d0/),shape(c))
    real(rp), parameter, dimension(2)   :: sigma = 1.d0/3.d0*(/1.d0,2.d0/)
    real(rp), parameter :: eps = 1.e-40
    integer :: a
    real(rp), dimension(-1:1) :: f
    real(rp), dimension(2) :: beta,we,dudlh,dvdlh,dwdlh
    real(rp) :: tauP
    real(rp) :: lmbd
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    mu  = mu12(2);  dmu  = mu12(1)-mu12(2)
    dxi = dli(1)
    dyi = dli(2)
    !
    lmbd = 1.d0/6.d0
    !
    !
    ! making an exception for this kernel -- private variables not explicitly mentioned for the sake of conciseness
    !                                        all scalars should be firstprivate/private
    !
    !$acc parallel loop collapse(3) default(present) private(f,beta,we,dudlh,dvdlh,dwdlh) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          u_ccl = u(i  ,j  ,k-2)
          u_ccm = u(i  ,j  ,k-1)
          u_pcm = u(i+1,j  ,k-1)
          u_cpm = u(i  ,j+1,k-1)
          u_clc = u(i  ,j-2,k  )
          u_cmc = u(i  ,j-1,k  )
          u_pmc = u(i+1,j-1,k  )
          u_lcc = u(i-2,j  ,k  )
          u_mcc = u(i-1,j  ,k  )
          u_ccc = u(i  ,j  ,k  )
          u_pcc = u(i+1,j  ,k  )
          u_qcc = u(i+2,j  ,k  )
          u_mpc = u(i-1,j+1,k  )
          u_cpc = u(i  ,j+1,k  )
          u_cqc = u(i  ,j+2,k  )
          u_cmp = u(i  ,j-1,k+1)
          u_mcp = u(i-1,j  ,k+1)
          u_ccp = u(i  ,j  ,k+1)
          u_ccq = u(i  ,j  ,k+2)
          !
          v_ccl = v(i  ,j  ,k-2)
          v_ccm = v(i  ,j  ,k-1)
          v_pcm = v(i+1,j  ,k-1)
          v_cpm = v(i  ,j+1,k-1)
          v_clc = v(i  ,j-2,k  )
          v_cmc = v(i  ,j-1,k  )
          v_pmc = v(i+1,j-1,k  )
          v_lcc = v(i-2,j  ,k  )
          v_mcc = v(i-1,j  ,k  )
          v_ccc = v(i  ,j  ,k  )
          v_pcc = v(i+1,j  ,k  )
          v_qcc = v(i+2,j  ,k  )
          v_mpc = v(i-1,j+1,k  )
          v_cpc = v(i  ,j+1,k  )
          v_cqc = v(i  ,j+2,k  )
          v_cmp = v(i  ,j-1,k+1)
          v_mcp = v(i-1,j  ,k+1)
          v_ccp = v(i  ,j  ,k+1)
          v_ccq = v(i  ,j  ,k+2)
          !
          w_ccl = w(i  ,j  ,k-2)
          w_ccm = w(i  ,j  ,k-1)
          w_pcm = w(i+1,j  ,k-1)
          w_cpm = w(i  ,j+1,k-1)
          w_clc = w(i  ,j-2,k  )
          w_cmc = w(i  ,j-1,k  )
          w_pmc = w(i+1,j-1,k  )
          w_lcc = w(i-2,j  ,k  )
          w_mcc = w(i-1,j  ,k  )
          w_ccc = w(i  ,j  ,k  )
          w_pcc = w(i+1,j  ,k  )
          w_qcc = w(i+2,j  ,k  )
          w_mpc = w(i-1,j+1,k  )
          w_cpc = w(i  ,j+1,k  )
          w_cqc = w(i  ,j+2,k  )
          w_cmp = w(i  ,j-1,k+1)
          w_mcp = w(i-1,j  ,k+1)
          w_ccp = w(i  ,j  ,k+1)
          w_ccq = w(i  ,j  ,k+2)
          !
          c_ccm = psi(i  ,j  ,k-1)
          c_pcm = psi(i+1,j  ,k-1)
          c_cpm = psi(i  ,j+1,k-1)
          c_cmc = psi(i  ,j-1,k  )
          c_pmc = psi(i+1,j-1,k  )
          c_mcc = psi(i-1,j  ,k  )
          c_ccc = psi(i  ,j  ,k  )
          c_pcc = psi(i+1,j  ,k  )
          c_mpc = psi(i-1,j+1,k  )
          c_cpc = psi(i  ,j+1,k  )
          c_cmp = psi(i  ,j-1,k+1)
          c_mcp = psi(i-1,j  ,k+1)
          c_ccp = psi(i  ,j  ,k+1)
          c_cpp = psi(i  ,j+1,k+1)
          c_ppc = psi(i+1,j+1,k  )
          c_pcp = psi(i+1,j  ,k+1)
          !
#if defined(_CONSERVATIVE_MOMENTUM)
          d_ccm = psio(i  ,j  ,k-1)
          d_pcm = psio(i+1,j  ,k-1)
          d_cpm = psio(i  ,j+1,k-1)
          d_cmc = psio(i  ,j-1,k  )
          d_pmc = psio(i+1,j-1,k  )
          d_mcc = psio(i-1,j  ,k  )
          d_ccc = psio(i  ,j  ,k  )
          d_pcc = psio(i+1,j  ,k  )
          d_mpc = psio(i-1,j+1,k  )
          d_cpc = psio(i  ,j+1,k  )
          d_cmp = psio(i  ,j-1,k+1)
          d_mcp = psio(i-1,j  ,k+1)
          d_ccp = psio(i  ,j  ,k+1)
          d_cpp = psio(i  ,j+1,k+1)
          d_ppc = psio(i+1,j+1,k  )
          d_pcp = psio(i+1,j  ,k+1)
          !
          rglrx_mcc = acdi_rglrx(i-1,j  ,k  )
          rglrx_ccc = acdi_rglrx(i  ,j  ,k  )
          rglrx_pcc = acdi_rglrx(i+1,j  ,k  )
          rglrx_mpc = acdi_rglrx(i-1,j+1,k  )
          rglrx_cpc = acdi_rglrx(i  ,j+1,k  )
          rglrx_mcp = acdi_rglrx(i-1,j  ,k+1)
          rglrx_ccp = acdi_rglrx(i  ,j  ,k+1)
          !
          rglry_cmc = acdi_rglry(i  ,j-1,k  )
          rglry_ccc = acdi_rglry(i  ,j  ,k  )
          rglry_cpc = acdi_rglry(i  ,j+1,k  )
          rglry_pmc = acdi_rglry(i+1,j-1,k  )
          rglry_pcc = acdi_rglry(i+1,j  ,k  )
          rglry_cmp = acdi_rglry(i  ,j-1,k+1)
          rglry_ccp = acdi_rglry(i  ,j  ,k+1)
          !
          rglrz_ccm = acdi_rglrz(i  ,j  ,k-1)
          rglrz_ccc = acdi_rglrz(i  ,j  ,k  )
          rglrz_ccp = acdi_rglrz(i  ,j  ,k+1)
          rglrz_pcm = acdi_rglrz(i+1,j  ,k-1)
          rglrz_pcc = acdi_rglrz(i+1,j  ,k  )
          rglrz_cpm = acdi_rglrz(i  ,j+1,k-1)
          rglrz_cpc = acdi_rglrz(i  ,j+1,k  )
#endif
          !
          dzci_c = dzci(k  )
          dzci_m = dzci(k-1)
          dzfi_c = dzfi(k  )
          dzfi_p = dzfi(k+1)
          !
          psixp = 0.5*(c_pcc+c_ccc)
          psiyp = 0.5*(c_cpc+c_ccc)
          psizp = 0.5*(c_ccp+c_ccc)
          !
          ! advection
          !
#if defined(_CONSERVATIVE_MOMENTUM)
          !
          ! n.b.: interpolations can be recycled from viscosity
          !
          rhoxp = rho + drho*d_pcc
          rhoxm = rho + drho*d_ccc
          rhoyp = rho + drho*0.25*(d_ccc+d_cpc+d_ppc+d_pcc)
          rhoym = rho + drho*0.25*(d_ccc+d_cmc+d_pmc+d_pcc)
          rhozp = rho + drho*0.25*(d_ccc+d_pcc+d_ccp+d_pcp)
          rhozm = rho + drho*0.25*(d_ccc+d_pcc+d_ccm+d_pcm)
          rhox  = rho + drho*psixp
          uuip  = 0.25*(u_pcc+u_ccc)*(u_ccc+u_pcc)*rhoxp
          uuim  = 0.25*(u_mcc+u_ccc)*(u_ccc+u_mcc)*rhoxm
          vujp  = 0.25*(v_pcc+v_ccc)*(u_ccc+u_cpc)*rhoyp
          vujm  = 0.25*(v_pmc+v_cmc)*(u_ccc+u_cmc)*rhoym
          wukp  = 0.25*(w_pcc+w_ccc)*(u_ccc+u_ccp)*rhozp
          wukm  = 0.25*(w_pcm+w_ccm)*(u_ccc+u_ccm)*rhozm
          dudt_aux = (dxi*( -uuip + uuim ) + dyi*( -vujp + vujm ) + dzfi_c*( -wukp + wukm ))/rhox
          !
          rhoxp = rho + drho*0.25*(d_ccc+d_pcc+d_ppc+d_cpc)
          rhoxm = rho + drho*0.25*(d_ccc+d_mcc+d_mpc+d_cpc)
          rhoyp = rho + drho*d_cpc
          rhoym = rho + drho*d_ccc
          rhozp = rho + drho*0.25*(d_ccc+d_cpc+d_cpp+d_ccp)
          rhozm = rho + drho*0.25*(d_ccc+d_cpc+d_cpm+d_ccm)
          rhoy  = rho + drho*psiyp
          uvip  = 0.25*(u_ccc+u_cpc)*(v_ccc+v_pcc)*rhoxp
          uvim  = 0.25*(u_mcc+u_mpc)*(v_ccc+v_mcc)*rhoxm
          vvjp  = 0.25*(v_ccc+v_cpc)*(v_ccc+v_cpc)*rhoyp
          vvjm  = 0.25*(v_ccc+v_cmc)*(v_ccc+v_cmc)*rhoym
          wvkp  = 0.25*(w_ccc+w_cpc)*(v_ccc+v_ccp)*rhozp
          wvkm  = 0.25*(w_ccm+w_cpm)*(v_ccc+v_ccm)*rhozm
          dvdt_aux = (dxi*( -uvip + uvim ) + dyi*( -vvjp + vvjm ) + dzfi_c*( -wvkp + wvkm ))/rhoy
          !
          rhoxp = rho + drho*0.25*(d_ccc+d_pcc+d_ccp+d_pcp)
          rhoxm = rho + drho*0.25*(d_ccc+d_mcc+d_ccp+d_mcp)
          rhoyp = rho + drho*0.25*(d_ccc+d_cpc+d_ccp+d_cpp)
          rhoym = rho + drho*0.25*(d_ccc+d_cmc+d_ccp+d_cmp)
          rhozp = rho + drho*d_ccp
          rhozm = rho + drho*d_ccc
          rhoz  = rho + drho*psizp
          uwip  = 0.25*(u_ccc+u_ccp)*(w_ccc+w_pcc)*rhoxp
          uwim  = 0.25*(u_mcc+u_mcp)*(w_ccc+w_mcc)*rhoxm
          vwjp  = 0.25*(v_ccc+v_ccp)*(w_ccc+w_cpc)*rhoyp
          vwjm  = 0.25*(v_cmc+v_cmp)*(w_ccc+w_cmc)*rhoym
          wwkp  = 0.25*(w_ccc+w_ccp)*(w_ccc+w_ccp)*rhozp
          wwkm  = 0.25*(w_ccc+w_ccm)*(w_ccc+w_ccm)*rhozm
          dwdt_aux = (dxi*( -uwip + uwim ) + dyi*( -vwjp + vwjm ) + dzci_c*( -wwkp + wwkm ))/rhoz
#else
          !!!QUICK
          !!if ((u_mcc+u_ccc+u_ccc+u_pcc) >= 0.) then
          !!  ududx = 0.25*(u_mcc+u_ccc+u_ccc+u_pcc)*(2.*u_pcc + 3.*u_ccc - 6.*u_mcc + 1.*u_lcc)*dxi/6.
          !!else
          !!  ududx = 0.25*(u_mcc+u_ccc+u_ccc+u_pcc)*(1.*u_qcc - 6.*u_pcc + 3.*u_ccc + 2.*u_mcc)*dxi/6.
          !!end if
          !if ((u_ccc) >= 0.) then
          !  !ududx = u_ccc*(2.*u_pcc + 3.*u_ccc - 6.*u_mcc + 1.*u_lcc)*dxi/6.
          !  ududx = u_ccc*(3.*u_pcc + 3.*u_ccc - 7.*u_mcc + 1.*u_lcc)*dxi/8.
          !  !ududx = u_ccc*(0.5*(u_ccc+u_pcc)-0.5*(u_mcc+u_ccc))*dxi
          !else
          !  !ududx = u_ccc*(1.*u_qcc - 6.*u_pcc + 3.*u_ccc + 2.*u_mcc)*dxi/6.
          !  ududx = u_ccc*(1.*u_qcc - 7.*u_pcc + 3.*u_ccc + 3.*u_mcc)*dxi/8.
          !  !ududx = u_ccc*(0.5*(u_ccc+u_pcc)-0.5*(u_mcc+u_ccc))*dxi
          !end if
          !if ((v_cmc+v_pmc+v_ccc+v_pcc) >= 0.) then
          !  vdudy = 0.25*(v_cmc+v_pmc+v_ccc+v_pcc)*(2.*u_cpc + 3.*u_ccc - 6.*u_cmc + 1.*u_clc)*dyi/6.
          !else
          !  vdudy = 0.25*(v_cmc+v_pmc+v_ccc+v_pcc)*(1.*u_cqc - 6.*u_cpc + 3.*u_ccc + 2.*u_cmc)*dyi/6.
          !end if
          !if ((w_ccm+w_pcm+w_ccc+w_pcc) >= 0.) then
          !  !wdudz = 0.25*(w_ccm+w_pcm+w_ccc+w_pcc)*(2.*u_ccp + 3.*u_ccc - 6.*u_ccm + 1.*u_ccl)*dzfi_c/6.
          !  wdudz = 0.25*(w_ccm+w_pcm+w_ccc+w_pcc)*(3.*u_ccp + 3.*u_ccc - 7.*u_ccm + 1.*u_ccl)*dzfi_c/8.
          !  !wdudz = 0.25*(w_ccm+w_pcm+w_ccc+w_pcc)*(0.5*(u_ccc+u_ccp)-0.5*(u_ccm+u_ccc))*dzfi_c
          !else
          !  !wdudz = 0.25*(w_ccm+w_pcm+w_ccc+w_pcc)*(1.*u_ccq - 6.*u_ccp + 3.*u_ccc + 2.*u_ccm)*dzfi_c/6.
          !  wdudz = 0.25*(w_ccm+w_pcm+w_ccc+w_pcc)*(1.*u_ccq - 7.*u_ccp + 3.*u_ccc + 3.*u_ccm)*dzfi_c/8.
          !  !wdudz = 0.25*(w_ccm+w_pcm+w_ccc+w_pcc)*(0.5*(u_ccc+u_ccp)-0.5*(u_ccm+u_ccc))*dzfi_c
          !end if
          !dudt_aux = -ududx -vdudy -wdudz
          !!
          !if ((u_mcc+u_ccc+u_mpc+u_cpc) >= 0.) then
          !  udvdx = 0.25*(u_mcc+u_ccc+u_mpc+u_cpc)*(2.*v_pcc + 3.*v_ccc - 6.*v_mcc + 1.*v_lcc)*dxi/6.
          !else
          !  udvdx = 0.25*(u_mcc+u_ccc+u_mpc+u_cpc)*(1.*v_qcc - 6.*v_pcc + 3.*v_ccc + 2.*v_mcc)*dxi/6.
          !end if
          !if ((v_ccc) >= 0.) then
          !  vdvdy = v_ccc*(2.*v_cpc + 3.*v_ccc - 6.*v_cmc + 1.*v_clc)*dyi/6.
          !else
          !  vdvdy = v_ccc*(1.*v_cqc - 6.*v_cpc + 3.*v_ccc + 2.*v_cmc)*dyi/6.
          !end if
          !!if ((v_cmc+v_ccc+v_ccc+v_cpc) >= 0.) then
          !!  vdvdy = 0.25*(v_cmc+v_ccc+v_ccc+v_cpc)*(2.*v_cpc + 3.*v_ccc - 6.*v_cmc + 1.*v_clc)*dyi/6.
          !!else
          !!  vdvdy = 0.25*(v_cmc+v_ccc+v_ccc+v_cpc)*(1.*v_cqc - 6.*v_cpc + 3.*v_ccc + 2.*v_cmc)*dyi/6.
          !!end if
          !if ((w_ccm+w_cpm+w_ccc+w_cpc) >= 0.) then
          !  wdvdz = 0.25*(w_ccm+w_cpm+w_ccc+w_cpc)*(2.*v_ccp + 3.*v_ccc - 6.*v_ccm + 1.*v_ccl)*dzfi_c/6.
          !else
          !  wdvdz = 0.25*(w_ccm+w_cpm+w_ccc+w_cpc)*(1.*v_ccq - 6.*v_ccp + 3.*v_ccc + 2.*v_ccm)*dzfi_c/6.
          !end if
          !dvdt_aux = -udvdx -vdvdy -wdvdz
          !!
          !if ((u_mcc+u_mcp+u_ccc+u_ccp) >= 0.) then
          !  !udwdx = 0.25*(u_mcc+u_mcp+u_ccc+u_ccp)*(2.*w_pcc + 3.*w_ccc - 6.*w_mcc + 1.*w_lcc)*dxi/6.
          !  udwdx = 0.25*(u_mcc+u_mcp+u_ccc+u_ccp)*(3.*w_pcc + 3.*w_ccc - 7.*w_mcc + 1.*w_lcc)*dxi/8.
          !  !udwdx = 0.25*(u_mcc+u_mcp+u_ccc+u_ccp)*(0.5*(w_ccc+w_pcc)-0.5*(w_mcc+w_ccc))*dxi
          !else
          !  !udwdx = 0.25*(u_mcc+u_mcp+u_ccc+u_ccp)*(1.*w_qcc - 6.*w_pcc + 3.*w_ccc + 2.*w_mcc)*dxi/6.
          !  udwdx = 0.25*(u_mcc+u_mcp+u_ccc+u_ccp)*(1.*w_qcc - 7.*w_pcc + 3.*w_ccc + 3.*w_mcc)*dxi/8.
          !  !udwdx = 0.25*(u_mcc+u_mcp+u_ccc+u_ccp)*(0.5*(w_ccc+w_pcc)-0.5*(w_mcc+w_ccc))*dxi
          !end if
          !if ((v_cmc+v_cmp+v_ccc+v_ccp) >= 0.) then
          !  vdwdy = 0.25*(v_cmc+v_cmp+v_ccc+v_ccp)*(2.*w_cpc + 3.*w_ccc - 6.*w_cmc + 1.*w_clc)*dyi/6.
          !else
          !  vdwdy = 0.25*(v_cmc+v_cmp+v_ccc+v_ccp)*(1.*w_cqc - 6.*w_cpc + 3.*w_ccc + 2.*w_cmc)*dyi/6.
          !end if
          !if ((w_ccc) >= 0.) then
          !  !wdwdz = w_ccc*(2.*w_ccp + 3.*w_ccc - 6.*w_ccm + 1.*w_ccl)*dzfi_c/6.
          !  wdwdz = w_ccc*(3.*w_ccp + 3.*w_ccc - 7.*w_ccm + 1.*w_ccl)*dzfi_c/8.
          !  !wdwdz = w_ccc*(0.5*(w_ccc+w_ccp)-0.5*(w_ccm+w_ccc))*dzfi_c
          !else
          !  !wdwdz = w_ccc*(1.*w_ccq - 6.*w_ccp + 3.*w_ccc + 2.*w_ccm)*dzfi_c/6.
          !  wdwdz = w_ccc*(1.*w_ccq - 7.*w_ccp + 3.*w_ccc + 3.*w_ccm)*dzfi_c/8.
          !  !wdwdz = w_ccc*(0.5*(w_ccc+w_ccp)-0.5*(w_ccm+w_ccc))*dzfi_c
          !end if
          !!if ((w_ccm+w_ccc+w_ccc+w_ccp) >= 0.) then
          !!  wdwdz = 0.25*(w_ccm+w_ccc+w_ccc+w_ccp)*(2.*w_ccp + 3.*w_ccc - 6.*w_ccm + 1.*w_ccl)*dzfi_c/6.
          !!else
          !!  wdwdz = 0.25*(w_ccm+w_ccc+w_ccc+w_ccp)*(1.*w_ccq - 6.*w_ccp + 3.*w_ccc + 2.*w_ccm)*dzfi_c/6.
          !!end if
          !dwdt_aux = -udwdx -vdwdy -wdwdz
!!!!!!!!!!!!!!!!!
          !!weno3
          !
          a = nint(sign(1.d0,u_ccc))
          f(-1) = a*(u(i-1*a,j,k) - u(i-2*a,j,k))*dli(1)
          f( 0) = a*(u(i+0*a,j,k) - u(i-1*a,j,k))*dli(1)
          f( 1) = a*(u(i+1*a,j,k) - u(i+0*a,j,k))*dli(1)
          beta(1) = (f(-1) - f(0))**2
          beta(2) = (f( 0) - f(1))**2
          tauP  = abs(0.5d0*(beta(1)+beta(2))-0.25d0*(f(-1)-f(1))**2)
          we(:) = sigma(:)*(1.d0+tauP/(beta(:)+eps)+1.d0/dli(1)**(lmbd)*((beta(:)+eps)/(tauP+eps)))
          we(:) = we(:)/sum(we(:))
          dudlh(1) = sum(c(:,1)*f(-1:0))
          dudlh(2) = sum(c(:,2)*f( 0:1))
          ududx    = u_ccc*sum(we(:)*dudlh(:))
          !
          a = nint(sign(1.d0,v_cmc+v_pmc+v_ccc+v_pcc))
          f(-1) = a*(u(i,j-1*a,k) - u(i,j-2*a,k))*dli(2)
          f( 0) = a*(u(i,j+0*a,k) - u(i,j-1*a,k))*dli(2)
          f( 1) = a*(u(i,j+1*a,k) - u(i,j+0*a,k))*dli(2)
          beta(1) = (f(-1) - f(0))**2
          beta(2) = (f( 0) - f(1))**2
          tauP  = abs(0.5d0*(beta(1)+beta(2))-0.25d0*(f(-1)-f(1))**2)
          we(:) = sigma(:)*(1.d0+tauP/(beta(:)+eps)+1.d0/dli(2)**(lmbd)*((beta(:)+eps)/(tauP+eps)))
          we(:) = we(:)/sum(we(:))
          dudlh(1) = sum(c(:,1)*f(-1:0))
          dudlh(2) = sum(c(:,2)*f( 0:1))
          vdudy    = 0.25*(v_cmc+v_pmc+v_ccc+v_pcc)*sum(we(:)*dudlh(:))
          !
          a = nint(sign(1.d0,w_ccm+w_pcm+w_ccc+w_pcc))
          if(a > 0.) then
            f(-1) = a*(u(i,j,k-1*a) - u(i,j,k-2*a))*dzci(k-2)
            f( 0) = a*(u(i,j,k+0*a) - u(i,j,k-1*a))*dzci(k-1)
            f( 1) = a*(u(i,j,k+1*a) - u(i,j,k+0*a))*dzci(k  )
          else
            f(-1) = a*(u(i,j,k-1*a) - u(i,j,k-2*a))*dzci(k+1)
            f( 0) = a*(u(i,j,k+0*a) - u(i,j,k-1*a))*dzci(k  )
            f( 1) = a*(u(i,j,k+1*a) - u(i,j,k+0*a))*dzci(k-1)
          endif
          beta(1) = (f(-1) - f(0))**2
          beta(2) = (f( 0) - f(1))**2
          tauP  = abs(0.5d0*(beta(1)+beta(2))-0.25d0*(f(-1)-f(1))**2)
          we(:) = sigma(:)*(1.d0+tauP/(beta(:)+eps)+1.d0/dzci(k)**(lmbd)*((beta(:)+eps)/(tauP+eps)))
          we(:) = we(:)/sum(we(:))
          dudlh(1) = sum(c(:,1)*f(-1:0))
          dudlh(2) = sum(c(:,2)*f( 0:1))
          wdudz    = 0.25*(w_ccm+w_pcm+w_ccc+w_pcc)*sum(we(:)*dudlh(:))
          !
          dudt_aux = -ududx -vdudy -wdudz
          !
          a = nint(sign(1.d0,u_mcc+u_ccc+u_mpc+u_cpc))
          f(-1) = a*(v(i-1*a,j,k) - v(i-2*a,j,k))*dli(1)
          f( 0) = a*(v(i+0*a,j,k) - v(i-1*a,j,k))*dli(1)
          f( 1) = a*(v(i+1*a,j,k) - v(i+0*a,j,k))*dli(1)
          beta(1) = (f(-1) - f(0))**2
          beta(2) = (f( 0) - f(1))**2
          tauP  = abs(0.5d0*(beta(1)+beta(2))-0.25d0*(f(-1)-f(1))**2)
          we(:) = sigma(:)*(1.d0+tauP/(beta(:)+eps)+1.d0/dli(1)**(lmbd)*((beta(:)+eps)/(tauP+eps)))
          we(:) = we(:)/sum(we(:))
          dvdlh(1) = sum(c(:,1)*f(-1:0))
          dvdlh(2) = sum(c(:,2)*f( 0:1))
          udvdx    = 0.25*(u_mcc+u_ccc+u_mpc+u_cpc)*sum(we(:)*dvdlh(:))
          !
          a = nint(sign(1.d0,v_ccc))
          f(-1) = a*(v(i,j-1*a,k) - v(i,j-2*a,k))*dli(2)
          f( 0) = a*(v(i,j+0*a,k) - v(i,j-1*a,k))*dli(2)
          f( 1) = a*(v(i,j+1*a,k) - v(i,j+0*a,k))*dli(2)
          beta(1) = (f(-1) - f(0))**2
          beta(2) = (f( 0) - f(1))**2
          tauP  = abs(0.5d0*(beta(1)+beta(2))-0.25d0*(f(-1)-f(1))**2)
          we(:) = sigma(:)*(1.d0+tauP/(beta(:)+eps)+1.d0/dli(2)**(lmbd)*((beta(:)+eps)/(tauP+eps)))
          we(:) = we(:)/sum(we(:))
          dvdlh(1) = sum(c(:,1)*f(-1:0))
          dvdlh(2) = sum(c(:,2)*f( 0:1))
          vdvdy    = v_ccc*sum(we(:)*dvdlh(:))
          !
          a = nint(sign(1.d0,w_ccm+w_cpm+w_ccc+w_cpc))
          if(a > 0.) then
            f(-1) = a*(v(i,j,k-1*a) - v(i,j,k-2*a))*dzci(k-2)
            f( 0) = a*(v(i,j,k+0*a) - v(i,j,k-1*a))*dzci(k-1)
            f( 1) = a*(v(i,j,k+1*a) - v(i,j,k+0*a))*dzci(k  )
          else
            f(-1) = a*(v(i,j,k-1*a) - v(i,j,k-2*a))*dzci(k+1)
            f( 0) = a*(v(i,j,k+0*a) - v(i,j,k-1*a))*dzci(k  )
            f( 1) = a*(v(i,j,k+1*a) - v(i,j,k+0*a))*dzci(k-1)
          endif
          beta(1) = (f(-1) - f(0))**2
          beta(2) = (f( 0) - f(1))**2
          tauP  = abs(0.5d0*(beta(1)+beta(2))-0.25d0*(f(-1)-f(1))**2)
          we(:) = sigma(:)*(1.d0+tauP/(beta(:)+eps)+1.d0/dzci(k)**(lmbd)*((beta(:)+eps)/(tauP+eps)))
          we(:) = we(:)/sum(we(:))
          dvdlh(1) = sum(c(:,1)*f(-1:0))
          dvdlh(2) = sum(c(:,2)*f( 0:1))
          wdvdz    = 0.25*(w_ccm+w_cpm+w_ccc+w_cpc)*sum(we(:)*dvdlh(:))
          !
          dvdt_aux = -udvdx -vdvdy -wdvdz
          !
          a = nint(sign(1.d0,u_mcc+u_mcp+u_ccc+u_ccp))
          f(-1) = a*(w(i-1*a,j,k) - w(i-2*a,j,k))*dli(1)
          f( 0) = a*(w(i+0*a,j,k) - w(i-1*a,j,k))*dli(1)
          f( 1) = a*(w(i+1*a,j,k) - w(i+0*a,j,k))*dli(1)
          beta(1) = (f(-1) - f(0))**2
          beta(2) = (f( 0) - f(1))**2
          tauP  = abs(0.5d0*(beta(1)+beta(2))-0.25d0*(f(-1)-f(1))**2)
          we(:) = sigma(:)*(1.d0+tauP/(beta(:)+eps)+1.d0/dli(1)**(lmbd)*((beta(:)+eps)/(tauP+eps)))
          we(:) = we(:)/sum(we(:))
          dwdlh(1) = sum(c(:,1)*f(-1:0))
          dwdlh(2) = sum(c(:,2)*f( 0:1))
          udwdx    = 0.25*(u_mcc+u_mcp+u_ccc+u_ccp)*sum(we(:)*dwdlh(:))
          !
          a = nint(sign(1.d0,v_cmc+v_cmp+v_ccc+v_ccp))
          f(-1) = a*(w(i,j-1*a,k) - w(i,j-2*a,k))*dli(2)
          f( 0) = a*(w(i,j+0*a,k) - w(i,j-1*a,k))*dli(2)
          f( 1) = a*(w(i,j+1*a,k) - w(i,j+0*a,k))*dli(2)
          beta(1) = (f(-1) - f(0))**2
          beta(2) = (f( 0) - f(1))**2
          tauP  = abs(0.5d0*(beta(1)+beta(2))-0.25d0*(f(-1)-f(1))**2)
          we(:) = sigma(:)*(1.d0+tauP/(beta(:)+eps)+1.d0/dli(2)**(lmbd)*((beta(:)+eps)/(tauP+eps)))
          we(:) = we(:)/sum(we(:))
          dwdlh(1) = sum(c(:,1)*f(-1:0))
          dwdlh(2) = sum(c(:,2)*f( 0:1))
          vdwdy    = 0.25*(v_cmc+v_cmp+v_ccc+v_ccp)*sum(we(:)*dwdlh(:))
          !
          a = nint(sign(1.d0,w_ccc))
          if(a > 0.) then
            f(-1) = a*(w(i,j,k-1*a) - w(i,j,k-2*a))*dzci(k-2)
            f( 0) = a*(w(i,j,k+0*a) - w(i,j,k-1*a))*dzci(k-1)
            f( 1) = a*(w(i,j,k+1*a) - w(i,j,k+0*a))*dzci(k  )
          else
            f(-1) = a*(w(i,j,k-1*a) - w(i,j,k-2*a))*dzci(k+1)
            f( 0) = a*(w(i,j,k+0*a) - w(i,j,k-1*a))*dzci(k  )
            f( 1) = a*(w(i,j,k+1*a) - w(i,j,k+0*a))*dzci(k-1)
          endif
          beta(1) = (f(-1) - f(0))**2
          beta(2) = (f( 0) - f(1))**2
          tauP  = abs(0.5d0*(beta(1)+beta(2))-0.25d0*(f(-1)-f(1))**2)
          we(:) = sigma(:)*(1.d0+tauP/(beta(:)+eps)+1.d0/dzci(k)**(lmbd)*((beta(:)+eps)/(tauP+eps)))
          we(:) = we(:)/sum(we(:))
          dwdlh(1) = sum(c(:,1)*f(-1:0))
          dwdlh(2) = sum(c(:,2)*f( 0:1))
          wdwdz    = w_ccc*sum(we(:)*dwdlh(:))
          !
          dwdt_aux = -udwdx -vdwdy -wdwdz
#endif
          !
          ! diffusion
          !
#if defined(_SPLIT_VISCOUS_DIFFUSION)
          dpsidyp = 0.5*(c_ppc+c_cpc-c_pcc-c_ccc)*dyi
          dpsidym = 0.5*(c_pcc+c_ccc-c_pmc-c_cmc)*dyi
          dpsidzp = 0.5*(c_pcp+c_ccp-c_pcc-c_ccc)*dzci_c
          dpsidzm = 0.5*(c_pcc+c_ccc-c_pcm-c_ccm)*dzci_m
          dpsidx  = (c_pcc-c_ccc)*dxi
          dpsidy  = 0.5*(dpsidyp+dpsidym)
          dpsidz  = 0.5*(dpsidzp+dpsidzm)
          dmudx = dmu*dpsidx
          dmudy = dmu*dpsidy
          dmudz = dmu*dpsidz
          mux  = mu + dmu*psixp
          rhoxp = rho + drho*psixp
          dudxp = (u_pcc-u_ccc)*dxi
          dudxm = (u_ccc-u_mcc)*dxi
          dudyp = (u_cpc-u_ccc)*dyi
          dudym = (u_ccc-u_cmc)*dyi
          dudzp = (u_ccp-u_ccc)*dzci_c
          dudzm = (u_ccc-u_ccm)*dzci_m
          dvdxp = (v_pcc-v_ccc)*dxi
          dvdxm = (v_pmc-v_cmc)*dxi
          dwdxp = (w_pcc-w_ccc)*dxi
          dwdxm = (w_pcm-w_ccm)*dxi
          omyp = (dudzp-dwdxp)
          omym = (dudzm-dwdxm)
          omzp = (dvdxp-dudyp)
          omzm = (dvdxm-dudym)
          domydz = (omyp-omym)*dzfi_c
          domzdy = (omzp-omzm)*dyi
          dudx = 0.5*(dudxp+dudxm)
          dudy = 0.5*(dudyp+dudym)
          dudz = 0.5*(dudzp+dudzm)
          dvdx = 0.5*(dvdxp+dvdxm)
          dwdx = 0.5*(dwdxp+dwdxm)
          dudt_aux = dudt_aux + (dmudx*(dudx+dudx) + dmudy*(dudy+dvdx) + dmudz*(dudz+dwdx) - &
                                mux*(domzdy-domydz))/rhoxp ! vorticity curl
          !dudt_aux = dudt_aux + (dmudx*(dudx+dudx) + dmudy*(dudy+dvdx) + dmudz*(dudz+dwdx) + &
          !                      mux*((dudxp-dudxm)*dxi+(dudyp-dudym)*dyi+(dudzp-dudzm)*dzfi_c))/rhoxp ! laplacian
          !
          dpsidxp = 0.5*(c_ppc+c_pcc-c_cpc-c_ccc)*dxi
          dpsidxm = 0.5*(c_cpc+c_ccc-c_mpc-c_mcc)*dxi
          dpsidzp = 0.5*(c_cpp+c_ccp-c_cpc-c_ccc)*dzci_c
          dpsidzm = 0.5*(c_cpc+c_ccc-c_cpm-c_ccm)*dzci_m
          dpsidx  = 0.5*(dpsidxp+dpsidxm)
          dpsidy  = (c_cpc-c_ccc)*dyi
          dpsidz  = 0.5*(dpsidzp+dpsidzm)
          dmudx = dmu*dpsidx
          dmudy = dmu*dpsidy
          dmudz = dmu*dpsidz
          muy  = mu + dmu*psiyp
          rhoyp = rho + drho*psiyp
          dvdxp = (v_pcc-v_ccc)*dxi
          dvdxm = (v_ccc-v_mcc)*dxi
          dvdyp = (v_cpc-v_ccc)*dyi
          dvdym = (v_ccc-v_cmc)*dyi
          dvdzp = (v_ccp-v_ccc)*dzci_c
          dvdzm = (v_ccc-v_ccm)*dzci_m
          dudyp = (u_cpc-u_ccc)*dyi
          dudym = (u_mpc-u_mcc)*dyi
          dwdyp = (w_cpc-w_ccc)*dyi
          dwdym = (w_cpm-w_ccm)*dyi
          omxp = (dwdyp-dvdzp)
          omxm = (dwdym-dvdzm)
          omzp = (dvdxp-dudyp)
          omzm = (dvdxm-dudym)
          domxdz = (omxp-omxm)*dzfi_c
          domzdx = (omzp-omzm)*dxi
          dvdx = 0.5*(dvdxp+dvdxm)
          dvdy = 0.5*(dvdyp+dvdym)
          dvdz = 0.5*(dvdzp+dvdzm)
          dudy = 0.5*(dudyp+dudym)
          dwdy = 0.5*(dwdyp+dwdym)
          dvdt_aux = dvdt_aux + (dmudx*(dudy+dvdx) + dmudy*(dvdy+dvdy) + dmudz*(dvdz+dwdy) - &
                                muy*(domxdz-domzdx))/rhoyp ! vorticity curl
          !dvdt_aux = dvdt_aux + (dmudx*(dudy+dvdx) + dmudy*(dvdy+dvdy) + dmudz*(dvdz+dwdy) + &
          !                      muy*((dvdxp-dvdxm)*dxi+(dvdyp-dvdym)*dyi+(dvdzp-dvdzm)*dzfi_c))/rhoyp ! laplacian
          !
          dpsidxp = 0.5*(c_pcp+c_pcc-c_ccp-c_ccc)*dxi
          dpsidxm = 0.5*(c_ccp+c_ccc-c_mcp-c_mcc)*dxi
          dpsidyp = 0.5*(c_cpp+c_cpc-c_ccp-c_ccc)*dyi
          dpsidym = 0.5*(c_ccp+c_ccc-c_cmp-c_cmc)*dyi
          dpsidx  = 0.5*(dpsidxp+dpsidxm)
          dpsidy  = 0.5*(dpsidyp+dpsidym)
          dpsidz  = (c_ccp-c_ccc)*dzci_c
          dmudx = dmu*dpsidx
          dmudy = dmu*dpsidy
          dmudz = dmu*dpsidz
          muz  = mu + dmu*psizp
          rhozp = rho + drho*psizp
          dwdxp = (w_pcc-w_ccc)*dxi
          dwdxm = (w_ccc-w_mcc)*dxi
          dwdyp = (w_cpc-w_ccc)*dyi
          dwdym = (w_ccc-w_cmc)*dyi
          dwdzp = (w_ccp-w_ccc)*dzfi_p
          dwdzm = (w_ccc-w_ccm)*dzfi_c
          dudzp = (u_ccp-u_ccc)*dzci_c
          dudzm = (u_mcp-u_mcc)*dzci_c
          dvdzp = (v_ccp-v_ccc)*dzci_c
          dvdzm = (v_cmp-v_cmc)*dzci_c
          omxp = (dwdyp-dvdzp)
          omxm = (dwdym-dvdzm)
          omyp = (dudzp-dwdxp)
          omym = (dudzm-dwdxm)
          domxdy = (omxp-omxm)*dyi
          domydx = (omyp-omym)*dxi
          dwdx = 0.5*(dwdxp+dwdxm)
          dwdy = 0.5*(dwdyp+dwdym)
          dwdz = 0.5*(dwdzp+dwdzm)
          dudz = 0.5*(dudzp+dudzm)
          dvdz = 0.5*(dvdzp+dvdzm)
          dwdt_aux = dwdt_aux + (dmudx*(dudz+dwdx) + dmudy*(dvdz+dwdy) + dmudz*(dwdz+dwdz) - &
                                muz*(domydx-domxdy))/rhozp ! vorticity curl
          !dwdt_aux = dwdt_aux + (dmudx*(dudz+dwdx) + dmudy*(dvdz+dwdy) + dmudz*(dwdz+dwdz) + &
          !                      muz*((dwdxp-dwdxm)*dxi+(dwdyp-dwdym)*dyi+(dwdzp-dwdzm)*dzci_c))/rhozp ! laplacian
#else
          dudxp = (u_pcc-u_ccc)*dxi
          dudxm = (u_ccc-u_mcc)*dxi
          dvdxp = (v_pcc-v_ccc)*dxi
          dvdxm = (v_pmc-v_cmc)*dxi
          dudyp = (u_cpc-u_ccc)*dyi
          dudym = (u_ccc-u_cmc)*dyi
          dudzp = (u_ccp-u_ccc)*dzci_c
          dudzm = (u_ccc-u_ccm)*dzci_m
          dwdxp = (w_pcc-w_ccc)*dxi
          dwdxm = (w_pcm-w_ccm)*dxi
          muxp = mu + dmu*c_pcc
          muxm = mu + dmu*c_ccc
          muyp = mu + dmu*0.25*(c_ccc+c_cpc+c_ppc+c_pcc)
          muym = mu + dmu*0.25*(c_ccc+c_cmc+c_pmc+c_pcc)
          muzp = mu + dmu*0.25*(c_ccc+c_pcc+c_ccp+c_pcp)
          muzm = mu + dmu*0.25*(c_ccc+c_pcc+c_ccm+c_pcm)
          rhoxp = rho + drho*psixp
          dudt_aux = dudt_aux + dxi*(   (dudxp+dudxp)*muxp-(dudxm+dudxm)*muxm)/rhoxp + &
                                dyi*(   (dudyp+dvdxp)*muyp-(dudym+dvdxm)*muym)/rhoxp + &
                                dzfi_c*((dudzp+dwdxp)*muzp-(dudzm+dwdxm)*muzm)/rhoxp
          !
          dvdxp = (v_pcc-v_ccc)*dxi
          dvdxm = (v_ccc-v_mcc)*dxi
          dudyp = (u_cpc-u_ccc)*dyi
          dudym = (u_mpc-u_mcc)*dyi
          dvdyp = (v_cpc-v_ccc)*dyi
          dvdym = (v_ccc-v_cmc)*dyi
          dvdzp = (v_ccp-v_ccc)*dzci_c
          dvdzm = (v_ccc-v_ccm)*dzci_m
          dwdyp = (w_cpc-w_ccc)*dyi
          dwdym = (w_cpm-w_ccm)*dyi
          muxp = mu + dmu*0.25*(c_ccc+c_pcc+c_ppc+c_cpc)
          muxm = mu + dmu*0.25*(c_ccc+c_mcc+c_mpc+c_cpc)
          muyp = mu + dmu*c_cpc
          muym = mu + dmu*c_ccc
          muzp = mu + dmu*0.25*(c_ccc+c_cpc+c_cpp+c_ccp)
          muzm = mu + dmu*0.25*(c_ccc+c_cpc+c_cpm+c_ccm)
          rhoyp = rho + drho*psiyp
          dvdt_aux = dvdt_aux + dxi*(   (dvdxp+dudyp)*muxp-(dvdxm+dudym)*muxm)/rhoyp + &
                                dyi*(   (dvdyp+dvdyp)*muyp-(dvdym+dvdym)*muym)/rhoyp + &
                                dzfi_c*((dvdzp+dwdyp)*muzp-(dvdzm+dwdym)*muzm)/rhoyp
          !
          dwdxp = (w_pcc-w_ccc)*dxi
          dwdxm = (w_ccc-w_mcc)*dxi
          dudzp = (u_ccp-u_ccc)*dzci_c
          dudzm = (u_mcp-u_mcc)*dzci_c
          dwdyp = (w_cpc-w_ccc)*dyi
          dwdym = (w_ccc-w_cmc)*dyi
          dvdzp = (v_ccp-v_ccc)*dzci_c
          dvdzm = (v_cmp-v_cmc)*dzci_c
          dwdzp = (w_ccp-w_ccc)*dzfi_p
          dwdzm = (w_ccc-w_ccm)*dzfi_c
          muxp = mu + dmu*0.25*(c_ccc+c_pcc+c_ccp+c_pcp)
          muxm = mu + dmu*0.25*(c_ccc+c_mcc+c_ccp+c_mcp)
          muyp = mu + dmu*0.25*(c_ccc+c_cpc+c_ccp+c_cpp)
          muym = mu + dmu*0.25*(c_ccc+c_cmc+c_ccp+c_cmp)
          muzp = mu + dmu*c_ccp
          muzm = mu + dmu*c_ccc
          rhozp = rho + drho*psizp
          dwdt_aux = dwdt_aux + dxi*(   (dwdxp+dudzp)*muxp-(dwdxm+dudzm)*muxm)/rhozp + &
                                dyi*(   (dwdyp+dvdzp)*muyp-(dwdym+dvdzm)*muym)/rhozp + &
                                dzci_c*((dwdzp+dwdzp)*muzp-(dwdzm+dwdzm)*muzm)/rhozp
#endif
#if defined(_CONSERVATIVE_MOMENTUM)
          !
          ! acdi interface regularization term
          ! (WIP: to be validated)
          !
          rxup = 0.25*drho*(rglrx_ccc+rglrx_pcc)*(u_ccc+u_pcc)
          rxum = 0.25*drho*(rglrx_mcc+rglrx_ccc)*(u_mcc+u_ccc)
          ryup = 0.25*drho*(rglry_ccc+rglry_pcc)*(u_ccc+u_cpc)
          ryum = 0.25*drho*(rglry_cmc+rglry_pmc)*(u_cmc+u_ccc)
          rzup = 0.25*drho*(rglrz_ccc+rglrz_pcc)*(u_ccc+u_ccp)
          rzum = 0.25*drho*(rglrz_ccm+rglrz_pcm)*(u_ccm+u_ccc)
          dudt_aux = dudt_aux + dxi*(   rxup-rxum)/rhoxp + &
                                dyi*(   ryup-ryum)/rhoxp + &
                                dzfi_c*(rzup-rzum)/rhoxp
          !
          rxvp = 0.25*drho*(rglrx_ccc+rglrx_cpc)*(v_ccc+v_pcc)
          rxvm = 0.25*drho*(rglrx_mcc+rglrx_mpc)*(v_mcc+v_ccc)
          ryvp = 0.25*drho*(rglry_ccc+rglry_cpc)*(v_ccc+v_cpc)
          ryvm = 0.25*drho*(rglry_cmc+rglry_ccc)*(v_cmc+v_ccc)
          rzvp = 0.25*drho*(rglrz_ccc+rglrz_cpc)*(v_ccc+v_ccp)
          rzvm = 0.25*drho*(rglrz_ccm+rglrz_cpm)*(v_ccm+v_ccc)
          dvdt_aux = dvdt_aux + dxi*(   rxvp-rxvm)/rhoyp + &
                                dyi*(   ryvp-ryvm)/rhoyp + &
                                dzfi_c*(rzvp-rzvm)/rhoyp
          !
          rxwp = 0.25*drho*(rglrx_ccc+rglrx_ccp)*(w_ccc+w_pcc)
          rxwm = 0.25*drho*(rglrx_mcc+rglrx_mcp)*(w_mcc+w_ccc)
          rywp = 0.25*drho*(rglry_ccc+rglry_ccp)*(w_ccc+w_cpc)
          rywm = 0.25*drho*(rglry_cmc+rglry_cmp)*(w_cmc+w_ccc)
          rzwp = 0.25*drho*(rglrz_ccc+rglrz_ccp)*(w_ccc+w_ccp)
          rzwm = 0.25*drho*(rglrz_ccm+rglrz_ccc)*(w_ccm+w_ccc)
          dwdt_aux = dwdt_aux + dxi*(   rxwp-rxwm)/rhozp + &
                                dyi*(   rywp-rywm)/rhozp + &
                                dzci_c*(rzwp-rzwm)/rhozp
#endif
          !
          dudt(i,j,k) = dudt_aux
          dvdt(i,j,k) = dvdt_aux
          dwdt(i,j,k) = dwdt_aux
        end do
      end do
    end do
  end subroutine mom_xyz_ad
  !
  subroutine mom_xyz_oth(n,dli,dzci,dzfi,rkpar,dt_r,rho12,beta12,bforce,gacc,sigma,rho0,rho_av, &
                         p,pn,po,psi,kappa,s,psio,kappao,dudt,dvdt,dwdt)
    implicit none
    integer , intent(in   ), dimension(3) :: n
    real(rp), intent(in   ), dimension(3) :: dli
    real(rp), intent(in   ), dimension(1-nh:) :: dzci,dzfi
    real(rp), intent(in   ), dimension(2) :: rkpar
    real(rp), intent(in   ) :: dt_r
    real(rp), intent(in   ), dimension(2) :: rho12,beta12
    real(rp), intent(in   ), dimension(3) :: bforce,gacc
    real(rp), intent(in   ) :: sigma
    real(rp), intent(in   ) :: rho0,rho_av
    real(rp), intent(in   ), dimension(1-nh:,1-nh:,1-nh:)    :: p,pn,po,psi,kappa,s
    real(rp), intent(in   ), dimension(1-nh:,1-nh:,1-nh:,1:), optional :: psio,kappao
    real(rp), intent(inout), dimension( :, :, :)    :: dudt,dvdt,dwdt
    integer :: i,j,k
    real(rp) :: rho,drho,rhobeta,drhobeta
    real(rp) :: dxi,dyi,sum_rkpar
    real(rp) :: c_ccm,c_pcm,c_cpm,c_cmc,c_pmc,c_mcc,c_ccc,c_pcc,c_mpc,c_cpc,c_cmp,c_mcp,c_ccp,c_cpp,c_ppc,c_pcp, &
                p_ccc,p_pcc,p_cpc,p_ccp, &
                pn_ccc,pn_pcc,pn_cpc,pn_ccp, &
                po_ccc,po_pcc,po_cpc,po_ccp, &
                q_ccc,q_pcc,q_cpc,q_ccp, &
                s_ccc,s_pcc,s_cpc,s_ccp, &
                k_ccc,k_pcc,k_cpc,k_ccp, &
                bforcex,bforcey,bforcez,gaccx,gaccy,gaccz, &
                dzci_c,dzci_m,dzfi_c,dzfi_p, &
                psixp,psiyp,psizp
    real(rp) :: rhoxp,rhoyp,rhozp, &
                kappaxp,kappayp,kappazp, &
                factorxp,factoryp,factorzp, &
                dpdx  ,dpdy  ,dpdz  , &
                dpdx_e,dpdy_e,dpdz_e, &
                surfx  ,surfy  ,surfz
#if defined(_SURFACE_TENSION_SPLITTING)
    real(rp) :: surfx_1,surfy_1,surfz_1, &
                surfx_2,surfy_2,surfz_2, &
                surfx_e,surfy_e,surfz_e
#endif
    real(rp) :: dudt_aux,dvdt_aux,dwdt_aux
    !
    rho     = rho12(2); drho = rho12(1)-rho12(2)
    rhobeta = rho12(2)*beta12(2); drhobeta = rho12(1)*beta12(1)- rho12(2)*beta12(2)
    dxi = dli(1)
    dyi = dli(2)
    sum_rkpar = sum(rkpar)
    !
    ! making an exception for this kernel -- private variables not explicitly mentioned for the sake of conciseness
    !                                        all scalars should be firstprivate/private
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          dzci_c = dzci(k  )
          !
          p_ccc = p(i  ,j  ,k  )
          p_pcc = p(i+1,j  ,k  )
          p_cpc = p(i  ,j+1,k  )
          p_ccp = p(i  ,j  ,k+1)
          pn_ccc = pn(i  ,j  ,k  )
          pn_pcc = pn(i+1,j  ,k  )
          pn_cpc = pn(i  ,j+1,k  )
          pn_ccp = pn(i  ,j  ,k+1)
          po_ccc = po(i  ,j  ,k  )
          po_pcc = po(i+1,j  ,k  )
          po_cpc = po(i  ,j+1,k  )
          po_ccp = po(i  ,j  ,k+1)
          !
#if defined(_CONSTANT_COEFFS_POISSON)
          q_ccc = p_ccc + dt_r*sum_rkpar*(pn_ccc-po_ccc)
          q_pcc = p_pcc + dt_r*sum_rkpar*(pn_pcc-po_pcc)
          q_cpc = p_cpc + dt_r*sum_rkpar*(pn_cpc-po_cpc)
          q_ccp = p_ccp + dt_r*sum_rkpar*(pn_ccp-po_ccp)
#endif
          !
#if defined(_SCALAR) && defined(_BOUSSINESQ_BUOYANCY)
          s_ccc = s(i  ,j  ,k  )
          s_pcc = s(i+1,j  ,k  )
          s_cpc = s(i  ,j+1,k  )
          s_ccp = s(i  ,j  ,k+1)
#endif
          !
#if defined(_CONSTANT_COEFFS_POISSON)
#if defined(_SURFACE_TENSION_SPLITTING)
          k_ccc = kappao(i  ,j  ,k  ,2)
          k_pcc = kappao(i+1,j  ,k  ,2)
          k_cpc = kappao(i  ,j+1,k  ,2)
          k_ccp = kappao(i  ,j  ,k+1,2)
          c_ccc = psio(i  ,j  ,k  ,2)
          c_pcc = psio(i+1,j  ,k  ,2)
          c_cpc = psio(i  ,j+1,k  ,2)
          c_ccp = psio(i  ,j  ,k+1,2)
          kappaxp = 0.5*(k_pcc+k_ccc)
          kappayp = 0.5*(k_cpc+k_ccc)
          kappazp = 0.5*(k_ccp+k_ccc)
          surfx_2 = sigma*kappaxp*(c_pcc-c_ccc)*dxi
          surfy_2 = sigma*kappayp*(c_cpc-c_ccc)*dyi
          surfz_2 = sigma*kappazp*(c_ccp-c_ccc)*dzci_c
          k_ccc = kappao(i  ,j  ,k  ,1)
          k_pcc = kappao(i+1,j  ,k  ,1)
          k_cpc = kappao(i  ,j+1,k  ,1)
          k_ccp = kappao(i  ,j  ,k+1,1)
          c_ccc = psio(i  ,j  ,k  ,1)
          c_pcc = psio(i+1,j  ,k  ,1)
          c_cpc = psio(i  ,j+1,k  ,1)
          c_ccp = psio(i  ,j  ,k+1,1)
          kappaxp = 0.5*(k_pcc+k_ccc)
          kappayp = 0.5*(k_cpc+k_ccc)
          kappazp = 0.5*(k_ccp+k_ccc)
          surfx_1 = sigma*kappaxp*(c_pcc-c_ccc)*dxi
          surfy_1 = sigma*kappayp*(c_cpc-c_ccc)*dyi
          surfz_1 = sigma*kappazp*(c_ccp-c_ccc)*dzci_c
          surfx_e = ((1.+dt_r)*surfx_1-dt_r*surfx_2)
          surfy_e = ((1.+dt_r)*surfy_1-dt_r*surfy_2)
          surfz_e = ((1.+dt_r)*surfz_1-dt_r*surfz_2)
#endif
#endif
          k_ccc = kappa(i  ,j  ,k  )
          k_pcc = kappa(i+1,j  ,k  )
          k_cpc = kappa(i  ,j+1,k  )
          k_ccp = kappa(i  ,j  ,k+1)
          !
          c_ccc = psi(i  ,j  ,k  )
          c_pcc = psi(i+1,j  ,k  )
          c_cpc = psi(i  ,j+1,k  )
          c_ccp = psi(i  ,j  ,k+1)
          !
          bforcex = bforce(1)
          bforcey = bforce(2)
          bforcez = bforce(3)
          gaccx = gacc(1)
          gaccy = gacc(2)
          gaccz = gacc(3)
          !
          psixp = 0.5*(c_pcc+c_ccc)
          psiyp = 0.5*(c_cpc+c_ccc)
          psizp = 0.5*(c_ccp+c_ccc)
          !
          rhoxp = rho + drho*psixp
          rhoyp = rho + drho*psiyp
          rhozp = rho + drho*psizp
          !
          ! pressure gradient and surface tension
          !
          dpdx = (p_pcc-p_ccc)*dxi
          dpdy = (p_cpc-p_ccc)*dyi
          dpdz = (p_ccp-p_ccc)*dzci_c
          dudt_aux = bforcex/rhoxp + gaccx*(1.-rho_av/rhoxp)
          dvdt_aux = bforcey/rhoyp + gaccy*(1.-rho_av/rhoyp)
          dwdt_aux = bforcez/rhozp + gaccz*(1.-rho_av/rhozp)
          !
          kappaxp = 0.5*(k_pcc+k_ccc)
          kappayp = 0.5*(k_cpc+k_ccc)
          kappazp = 0.5*(k_ccp+k_ccc)
          surfx   = sigma*kappaxp*(c_pcc-c_ccc)*dxi
          surfy   = sigma*kappayp*(c_cpc-c_ccc)*dyi
          surfz   = sigma*kappazp*(c_ccp-c_ccc)*dzci_c
#if defined(_CONSTANT_COEFFS_POISSON)
          dpdx_e = (q_pcc-q_ccc)*dxi
          dpdy_e = (q_cpc-q_ccc)*dyi
          dpdz_e = (q_ccp-q_ccc)*dzci_c
#if defined(_SURFACE_TENSION_SPLITTING)
          dudt_aux = dudt_aux + (-dpdx + surfx)/rho0 + (1./rhoxp-1./rho0)*(-dpdx_e + surfx_e)
          dvdt_aux = dvdt_aux + (-dpdy + surfy)/rho0 + (1./rhoyp-1./rho0)*(-dpdy_e + surfy_e)
          dwdt_aux = dwdt_aux + (-dpdz + surfz)/rho0 + (1./rhozp-1./rho0)*(-dpdz_e + surfz_e)
#else
          dudt_aux = dudt_aux - dpdx/rho0 + surfx/rhoxp - (1./rhoxp-1./rho0)*dpdx_e
          dvdt_aux = dvdt_aux - dpdy/rho0 + surfy/rhoyp - (1./rhoyp-1./rho0)*dpdy_e
          dwdt_aux = dwdt_aux - dpdz/rho0 + surfz/rhozp - (1./rhozp-1./rho0)*dpdz_e
#endif
#else
          dudt_aux = dudt_aux + (-dpdx + surfx)/rhoxp
          dvdt_aux = dvdt_aux + (-dpdy + surfy)/rhoyp
          dwdt_aux = dwdt_aux + (-dpdz + surfz)/rhozp
#endif
          !
          ! buoyancy
          !
#if defined(_SCALAR) && defined(_BOUSSINESQ_BUOYANCY)
          factorxp = rhobeta + drhobeta*psixp
          factoryp = rhobeta + drhobeta*psiyp
          factorzp = rhobeta + drhobeta*psizp
          dudt_aux = dudt_aux - gaccx*factorxp*0.5*(s_pcc+s_ccc)/rhoxp
          dvdt_aux = dvdt_aux - gaccy*factoryp*0.5*(s_cpc+s_ccc)/rhoyp
          dwdt_aux = dwdt_aux - gaccz*factorzp*0.5*(s_ccp+s_ccc)/rhozp
#endif
          dudt(i,j,k) = dudt_aux
          dvdt(i,j,k) = dvdt_aux
          dwdt(i,j,k) = dwdt_aux
        end do
      end do
    end do
  end subroutine mom_xyz_oth
  !
  subroutine cmpt_wallshear(n,is_bound,l,dli,dzci,dzfi,mu12,psi,u,v,w,taux,tauy,tauz)
    !
    ! n.b.: `is_bound` should exclude periodic BCs
    !
    implicit none
    integer , intent(in ), dimension(3) :: n
    logical , intent(in ), dimension(0:1,3) :: is_bound
    real(rp), intent(in ), dimension(3)     :: l,dli
    real(rp), intent(in ), dimension(1-nh:)    :: dzci,dzfi
    real(rp), intent(in )                   :: mu12(2)
    real(rp), intent(in ), dimension(1-nh:,1-nh:,1-nh:) :: psi,u,v,w
    real(rp), intent(out), dimension(3) :: taux,tauy,tauz
    real(rp) :: dudyp,dudym,dudzp,dudzm, &
                dvdxp,dvdxm,dvdzp,dvdzm, &
                dwdxp,dwdxm,dwdyp,dwdym, &
                muxp,muxm,muyp,muym,muzp,muzm
    !
    real(rp) :: tau21,tau31,tau12,tau32,tau13,tau23
    integer :: i,j,k,nx,ny,nz
    real(rp) :: dxi,dyi,lx,ly,lz,mu,dmu
    real(rp) :: tau(3,3)
    integer :: ierr
    !
    mu  = mu12(2); dmu  = mu12(1)-mu12(2)
    nx = n(1); ny = n(2); nz = n(3)
    dxi = dli(1); dyi = dli(2)
    lx = l(1); ly = l(2); lz = l(3)
    tau21 = 0._rp
    !$acc data copy(tau21) async(1)
    if(is_bound(0,2)) then
      j = 0
      !$acc parallel loop collapse(2) default(present) private(dvdxp,dudyp,muyp) &
      !$acc reduction(+:tau21) async(1)
      do k=1,nz
        do i=1,nx
          dvdxp = (v(i+1,j  ,k  )-v(i  ,j  ,k  ))*dxi
          dudyp = (u(i  ,j+1,k  )-u(i  ,j  ,k  ))*dyi
          muyp = mu + dmu*0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i+1,j+1,k)+psi(i+1,j,k))
          !
          tau21 = tau21 + (dudyp+dvdxp)*muyp/(dxi*dzfi(k)*lx*lz)
        end do
      end do
    end if
    if(is_bound(1,2)) then
      j = ny+1
      !$acc parallel loop collapse(2) default(present) private(dvdxm,dudym,muym) &
      !$acc reduction(+:tau21) async(1)
      do k=1,nz
        do i=1,nx
          dvdxm = (v(i+1,j-1,k  )-v(i  ,j-1,k  ))*dxi
          dudym = (u(i  ,j  ,k  )-u(i  ,j-1,k  ))*dyi
          muym = mu + dmu*0.25*(psi(i,j,k)+psi(i,j-1,k)+psi(i+1,j-1,k)+psi(i+1,j,k))
          !
          tau21 = tau21 - (dudym+dvdxm)*muym/(dxi*dzfi(k)*lx*lz)
        end do
      end do
    end if
    !$acc end data
    tau31 = 0._rp
    !$acc data copy(tau31) async(1)
    if(is_bound(0,3)) then
      k = 0
      !$acc parallel loop collapse(2) default(present) private(dudzp,dwdxp,muzp) &
      !$acc reduction(+:tau31) async(1)
      do j=1,ny
        do i=1,nx
          dudzp = (u(i  ,j  ,k+1)-u(i  ,j  ,k  ))*dzci(k  )
          dwdxp = (w(i+1,j  ,k  )-w(i  ,j  ,k  ))*dxi
          muzp = mu + dmu*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i+1,j,k+1)+psi(i+1,j,k))
          !
          tau31 = tau31 + (dudzp+dwdxp)*muzp/(dxi*dyi*lx*ly)
        end do
      end do
    end if
    if(is_bound(1,3)) then
      k = nz+1
      !$acc parallel loop collapse(2) default(present) private(dudzm,dwdxm,muzm) &
      !$acc reduction(+:tau31) async(1)
      do j=1,ny
        do i=1,nx
          dudzm = (u(i  ,j  ,k  )-u(i  ,j  ,k-1))*dzci(k-1)
          dwdxm = (w(i+1,j  ,k-1)-w(i  ,j  ,k-1))*dxi
          muzm = mu + dmu*0.25*(psi(i,j,k)+psi(i,j,k-1)+psi(i+1,j,k-1)+psi(i+1,j,k))
          !
          tau31 = tau31 - (dudzm+dwdxm)*muzm/(dxi*dyi*lx*ly)
        end do
      end do
    end if
    !$acc end data
    !
    tau12 = 0._rp
    !$acc data copy(tau12) async(1)
    if(is_bound(0,1)) then
      i = 0
      !$acc parallel loop collapse(2) default(present) private(dvdxp,dudyp,muxp) &
      !$acc reduction(+:tau12) async(1)
      do k=1,nz
        do j=1,ny
          dvdxp = (v(i+1,j  ,k  )-v(i  ,j  ,k  ))*dxi
          dudyp = (u(i  ,j+1,k  )-u(i  ,j  ,k  ))*dyi
          muxp = mu + dmu*0.25*(psi(i,j,k)+psi(i+1,j,k)+psi(i+1,j+1,k)+psi(i,j+1,k))
          !
          tau12 = tau12 + (dvdxp+dudyp)*muxp/(dyi*dzfi(k)*ly*lz)
        end do
      end do
    end if
    if(is_bound(1,1)) then
      i = nx+1
      !$acc parallel loop collapse(2) default(present) private(dvdxm,dudym,muxm) &
      !$acc reduction(+:tau12) async(1)
      do k=1,nz
        do j=1,ny
          dvdxm = (v(i  ,j  ,k  )-v(i-1,j  ,k  ))*dxi
          dudym = (u(i-1,j+1,k  )-u(i-1,j  ,k  ))*dyi
          muxm = mu + dmu*0.25*(psi(i,j,k)+psi(i-1,j,k)+psi(i-1,j+1,k)+psi(i,j+1,k))
          !
          tau12 = tau12 - (dvdxm+dudym)*muxm/(dyi*dzfi(k)*ly*lz)
        end do
      end do
    end if
    !$acc end data
    tau32 = 0._rp
    !$acc data copy(tau32) async(1)
    if(is_bound(0,3)) then
      k = 0
      !$acc parallel loop collapse(2) default(present) private(dvdzp,dwdyp,muzp) &
      !$acc reduction(+:tau32) async(1)
      do j=1,ny
        do i=1,nx
          dvdzp = (v(i  ,j  ,k+1)-v(i  ,j  ,k  ))*dzci(k  )
          dwdyp = (w(i  ,j+1,k  )-w(i  ,j  ,k  ))*dyi
          muzp = mu + dmu*0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i,j+1,k+1)+psi(i,j,k+1))
          !
          tau32 = tau32 + (dvdzp+dwdyp)*muzp/(dxi*dyi*lx*ly)
        end do
      end do
    end if
    if(is_bound(1,3)) then
      k = nz+1
      !$acc parallel loop collapse(2) default(present) private(dvdzm,dwdym,muzm) &
      !$acc reduction(+:tau32) async(1)
      do j=1,ny
        do i=1,nx
          dvdzm = (v(i  ,j  ,k  )-v(i  ,j  ,k-1))*dzci(k-1)
          dwdym = (w(i  ,j+1,k-1)-w(i  ,j  ,k-1))*dyi
          muzm = mu + dmu*0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i,j+1,k-1)+psi(i,j,k-1))
          !
          tau32 = tau32 - (dvdzm+dwdym)*muzm/(dxi*dyi*lx*ly)
        end do
      end do
    end if
    !$acc end data
    !
    tau13 = 0._rp
    !$acc data copy(tau13) async(1)
    if(is_bound(0,1)) then
      i = 0
      !$acc parallel loop collapse(2) default(present) private(dwdxp,dudzp,muxp) &
      !$acc reduction(+:tau13) async(1)
      do k=1,nz
        do j=1,ny
          dwdxp = (w(i+1,j  ,k  )-w(i  ,j  ,k))*dxi
          dudzp = (u(i  ,j  ,k+1)-u(i  ,j  ,k))*dzci(k  )
          muxp = mu + dmu*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i+1,j ,k+1)+psi(i+1,j ,k) )
          !
          tau13 = tau13 + (dwdxp+dudzp)*muxp/(dyi*dzfi(k)*ly*lz)
        end do
      end do
    end if
    if(is_bound(1,1)) then
      i = nx+1
      !$acc parallel loop collapse(2) default(present) private(dwdxm,dudzm,muxm) &
      !$acc reduction(+:tau13) async(1)
      do k=1,nz
        do j=1,ny
          dwdxm = (w(i  ,j  ,k  )-w(i-1,j  ,k))*dxi
          dudzm = (u(i-1,j  ,k+1)-u(i-1,j  ,k))*dzci(k  )
          muxm = mu + dmu*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i-1,j ,k+1)+psi(i-1,j ,k) )
          !
          tau13 = tau13 - (dwdxm+dudzm)*muxm/(dyi*dzfi(k)*ly*lz)
        end do
      end do
    end if
    !$acc end data
    tau23 = 0._rp
    !$acc data copy(tau23) async(1)
    if(is_bound(0,2)) then
      j = 0
      !$acc parallel loop collapse(2) default(present) private(dwdyp,dvdzp,muyp) &
      !$acc reduction(+:tau23) async(1)
      do k=1,nz
        do i=1,nx
          dwdyp = (w(i  ,j+1,k  )-w(i  ,j  ,k))*dyi
          dvdzp = (v(i  ,j  ,k+1)-v(i  ,j  ,k))*dzci(k  )
          muyp = mu + dmu*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i ,j+1,k+1)+psi(i ,j+1,k) )
          !
          tau23 = tau23 + (dwdyp+dvdzp)*muyp/(dxi*dzfi(k)*lx*lz)
        end do
      end do
    end if
    if(is_bound(1,2)) then
      j = ny+1
      !$acc parallel loop collapse(2) default(present) private(dwdym,dvdzm,muym) &
      !$acc reduction(+:tau23) async(1)
      do k=1,nz
        do i=1,nx
          dwdym = (w(i  ,j  ,k  )-w(i  ,j-1,k))*dyi
          dvdzm = (v(i  ,j-1,k+1)-v(i  ,j-1,k))*dzci(k  )
          muym = mu + dmu*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i ,j-1,k+1)+psi(i ,j-1,k) )
          !
          tau23 = tau23 - (dwdym+dvdzm)*muym/(dxi*dzfi(k)*lx*lz)
        end do
      end do
    end if
    !$acc end data
    !$acc wait(1)
    tau(:,:) = 0._rp
    tau(2,1) = tau21
    tau(3,1) = tau31
    tau(1,2) = tau12
    tau(3,2) = tau32
    tau(1,3) = tau13
    tau(2,3) = tau23
    call MPI_ALLREDUCE(MPI_IN_PLACE,tau(1,1),9,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    taux(:) = tau(:,1)
    tauy(:) = tau(:,2)
    tauz(:) = tau(:,3)
  end subroutine cmpt_wallshear
end module mod_mom
