! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 The CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_mom
  use mpi
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
    real(rp), intent(in), dimension(0:) :: dzfi
    real(rp), dimension(0:,0:,0:), intent(in   ) :: u,v,w
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
    real(rp), intent(in), dimension(0:) :: dzfi
    real(rp), dimension(0:,0:,0:), intent(in   ) :: u,v,w
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
    real(rp), intent(in), dimension(0:) :: dzci
    real(rp), dimension(0:,0:,0:), intent(in   ) :: u,v,w
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
  subroutine momx_d(nx,ny,nz,dxi,dyi,dzci,dzfi,mu12,psi,u,v,w,dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), intent(in), dimension(2) :: mu12
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm, &
                dvdxp,dvdxm,dwdxp,dwdxm
    real(rp) :: muxp,muxm,muyp,muym,muzp,muzm
    integer :: i,j,k
    real(rp) :: mu,dmu
    !
    mu  = mu12(2) ; dmu  = mu12(1)-mu12(2)
    !
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(dudxp,dudxm,dudyp,dudym,dudzp,dudzm,dvdxp,dvdxm,dwdxp,dwdxm) &
    !$acc private(muxp,muxm,muyp,muym,muzp,muzm) async(1)
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
          !
          dudt(i,j,k) = dudt(i,j,k) + &
                        dxi*(    (dudxp+dudxp)*muxp-(dudxm+dudxm)*muxm) + &
                        dyi*(    (dudyp+dvdxp)*muyp-(dudym+dvdxm)*muym) + &
                        dzfi(k)*((dudzp+dwdxp)*muzp-(dudzm+dwdxm)*muzm)
        end do
      end do
    end do
  end subroutine momx_d
  !
  subroutine momy_d(nx,ny,nz,dxi,dyi,dzci,dzfi,mu12,psi,u,v,w,dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), intent(in), dimension(2) :: mu12
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    real(rp) :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm, &
                dudyp,dudym,dwdyp,dwdym
    real(rp) :: muxp,muxm,muyp,muym,muzp,muzm
    integer :: i,j,k
    real(rp) :: mu,dmu
    !
    mu  = mu12(2) ; dmu  = mu12(1)-mu12(2)
    !
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm,dudyp,dudym,dwdyp,dwdym) &
    !$acc private(muxp,muxm,muyp,muym,muzp,muzm) async(1)
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
          !
          dvdt(i,j,k) = dvdt(i,j,k) + &
                        dxi*(    (dvdxp+dudyp)*muxp-(dvdxm+dudym)*muxm) + &
                        dyi*(    (dvdyp+dvdyp)*muyp-(dvdym+dvdym)*muym) + &
                        dzfi(k)*((dvdzp+dwdyp)*muzp-(dvdzm+dwdym)*muzm)
        end do
      end do
    end do
  end subroutine momy_d
  !
  subroutine momz_d(nx,ny,nz,dxi,dyi,dzci,dzfi,mu12,psi,u,v,w,dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), intent(in), dimension(2) :: mu12
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    integer :: i,j,k
    real(rp) :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm, &
                dudzp,dudzm,dvdzp,dvdzm
    real(rp) :: muxp,muxm,muyp,muym,muzp,muzm
    real(rp) :: mu,dmu
    !
    mu  = mu12(2) ; dmu  = mu12(1)-mu12(2)
    !
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm,dudzp,dudzm,dvdzp,dvdzm) &
    !$acc private(muxp,muxm,muyp,muym,muzp,muzm) async(1)
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
          !
          dwdt(i,j,k) = dwdt(i,j,k) + &
                        dxi*(    (dwdxp+dudzp)*muxp-(dwdxm+dudzm)*muxm) + &
                        dyi*(    (dwdyp+dvdzp)*muyp-(dwdym+dvdzm)*muym) + &
                        dzci(k)*((dwdzp+dwdzp)*muzp-(dwdzm+dwdzm)*muzm)
        end do
      end do
    end do
  end subroutine momz_d
  !
  subroutine momx_p(nx,ny,nz,dxi,dt_r,bforce,gacc,rho0,rhox_av,rho12,psi,p,pp,dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dt_r
    real(rp), intent(in) :: bforce,gacc,rho0,rhox_av,rho12(2)
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,p,pp
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    real(rp) :: rhop,dpdl
    integer :: i,j,k
    real(rp) :: rho,drho
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    !
    !$acc parallel loop collapse(3) default(present) private(dpdl,rhop) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          rhop = rho + drho*0.5*(psi(i+1,j,k)+psi(i,j,k))
          dpdl = (p(i+1,j,k)-p(i,j,k))*dxi
          !
          dudt(i,j,k) = dudt(i,j,k) + bforce + gacc*(rhop-rhox_av) &
#if defined(_CONSTANT_COEFFS_POISSON)
                        -dpdl*rhop/rho0 - (1.-rhop/rho0)*( ((1.+dt_r)*p(i+1,j,k)-dt_r*(pp(i+1,j,k))) - &
                                                           ((1.+dt_r)*p(i  ,j,k)-dt_r*(pp(i  ,j,k))) &
                                                         )*dxi
#else
                        -dpdl
#endif
        end do
      end do
    end do
  end subroutine momx_p
  !
  subroutine momy_p(nx,ny,nz,dyi,dt_r,bforce,gacc,rho0,rhoy_av,rho12,psi,p,pp,dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dyi,dt_r
    real(rp), intent(in) :: bforce,gacc,rho0,rhoy_av,rho12(2)
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,p,pp
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
          dvdt(i,j,k) = dvdt(i,j,k) + bforce + gacc*(rhop-rhoy_av) &
#if defined(_CONSTANT_COEFFS_POISSON)
                        -dpdl*rhop/rho0 - (1.-1.*rhop/rho0)*( ((1.+dt_r)*p(i,j+1,k)-dt_r*(pp(i,j+1,k))) - &
                                                              ((1.+dt_r)*p(i,j  ,k)-dt_r*(pp(i,j  ,k))) &
                                                            )*dyi
#else
                        -dpdl
#endif
        end do
      end do
    end do
  end subroutine momy_p
  !
  subroutine momz_p(nx,ny,nz,dzci,dt_r,bforce,gacc,rho0,rhoz_av,rho12,psi,p,pp,dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in), dimension(0:) :: dzci
    real(rp), intent(in) :: dt_r
    real(rp), intent(in) :: bforce,gacc,rho0,rhoz_av,rho12(2)
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,p,pp
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
          dwdt(i,j,k) = dwdt(i,j,k) + bforce + gacc*(rhop-rhoz_av) &
#if defined(_CONSTANT_COEFFS_POISSON)
                        -dpdl*rhop/rho0 - (1.-1.*rhop/rho0)*( ((1.+dt_r)*p(i,j,k+1)-dt_r*(pp(i,j,k+1))) - &
                                                              ((1.+dt_r)*p(i,j,k  )-dt_r*(pp(i,j,k  ))) &
                                                            )*dzci(k)
#else
                        -dpdl
#endif
        end do
      end do
    end do
  end subroutine momz_p
  !
  subroutine momx_sigma(nx,ny,nz,dxi,sigma,kappa,psi,dudt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi
    real(rp), intent(in) :: sigma
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,kappa
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    real(rp) :: kappap
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) private(kappap) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          kappap = 0.5*(kappa(i+1,j,k)+kappa(i,j,k))
          !
          dudt(i,j,k) = dudt(i,j,k) + &
                        sigma*kappap*(psi(i+1,j,k)-psi(i,j,k))*dxi
        end do
      end do
    end do
  end subroutine momx_sigma
  !
  subroutine momy_sigma(nx,ny,nz,dyi,sigma,psi,kappa,dvdt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dyi
    real(rp), intent(in) :: sigma
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,kappa
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    real(rp) :: kappap
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) private(kappap) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          kappap = 0.5*(kappa(i,j+1,k)+kappa(i,j,k))
          !
          dvdt(i,j,k) = dvdt(i,j,k) + &
                        sigma*kappap*(psi(i,j+1,k)-psi(i,j,k))*dyi
        end do
      end do
    end do
  end subroutine momy_sigma
  !
  subroutine momz_sigma(nx,ny,nz,dzci,sigma,psi,kappa,dwdt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in), dimension(0:) :: dzci
    real(rp), intent(in) :: sigma
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,kappa
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    real(rp) :: kappap
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) private(kappap) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          kappap = 0.5*(kappa(i,j,k+1)+kappa(i,j,k))
          !
          dwdt(i,j,k) = dwdt(i,j,k) + &
                        sigma*kappap*(psi(i,j,k+1)-psi(i,j,k))*dzci(k)
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
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,s
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    real(rp) :: psip,factorp
    integer :: i,j,k
    real(rp) :: rhobeta,drhobeta
    !
    rhobeta = rho12(2)*beta12(2); drhobeta = rho12(1)*beta12(1)-rho12(2)*beta12(2)
    !
    !$acc parallel loop collapse(3) default(present) private(psip,factorp) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          psip = 0.5*(psi(i+1,j,k)+psi(i,j,k))
          factorp = rhobeta + drhobeta*psip
          !
          dudt(i,j,k) = dudt(i,j,k) - &
                        gacc*factorp*0.5*(s(i+1,j,k)+s(i,j,k))
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
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,s
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    real(rp) :: psip,factorp
    integer :: i,j,k
    real(rp) :: rhobeta,drhobeta
    !
    rhobeta = rho12(2)*beta12(2); drhobeta = rho12(1)*beta12(1)-rho12(2)*beta12(2)
    !
    !$acc parallel loop collapse(3) default(present) private(factorp) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          psip = 0.5*(psi(i,j+1,k)+psi(i,j,k))
          factorp = rhobeta + drhobeta*psip
          !
          dvdt(i,j,k) = dvdt(i,j,k) - &
                        gacc*factorp*0.5*(s(i,j+1,k)+s(i,j,k))
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
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,s
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    real(rp) :: psip,factorp
    integer :: i,j,k
    real(rp) :: rhobeta,drhobeta
    !
    rhobeta = rho12(2)*beta12(2); drhobeta = rho12(1)*beta12(1)-rho12(2)*beta12(2)
    !
    !$acc parallel loop collapse(3) default(present) private(factorp) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          psip = 0.5*(psi(i,j,k+1)+psi(i,j,k))
          factorp = rhobeta + drhobeta*psip
          !
          dwdt(i,j,k) = dwdt(i,j,k) - &
                        gacc*factorp*0.5*(s(i,j,k+1)+s(i,j,k))
        end do
      end do
    end do
  end subroutine momz_buoy
  !
  subroutine mom_xyz_ad(n,dli,dzci,dzfi,rho12,mu12,psiflx_x,psiflx_y,psiflx_z,u,v,w,psi,dudt,dvdt,dwdt)
    implicit none
    integer , intent(in ), dimension(3) :: n
    real(rp), intent(in ), dimension(3) :: dli
    real(rp), intent(in ), dimension(0:) :: dzci,dzfi
    real(rp), intent(in ), dimension(2) :: rho12,mu12
    real(rp), intent(in ), dimension(0:,0:,0:) :: u,v,w,psi
    real(rp), intent(in ), dimension(0:,0:,0:), optional :: psiflx_x,psiflx_y,psiflx_z
    real(rp), intent(out), dimension( :, :, :) :: dudt,dvdt,dwdt
    integer :: i,j,k
    real(rp) :: rho,drho,mu,dmu,rhobeta,drhobeta
    real(rp) :: dxi,dyi
    real(rp) :: c_ccm,c_pcm,c_cpm,c_cmc,c_pmc,c_mcc,c_ccc,c_pcc,c_mpc,c_cpc,c_cmp,c_mcp,c_ccp,c_cpp,c_ppc,c_pcp, &
                d_ccm,d_pcm,d_cpm,d_cmc,d_pmc,d_mcc,d_ccc,d_pcc,d_mpc,d_cpc,d_cmp,d_mcp,d_ccp,d_cpp,d_ppc,d_pcp, &
                u_ccm,u_pcm,u_cpm,u_cmc,u_pmc,u_mcc,u_ccc,u_pcc,u_mpc,u_cpc,u_cmp,u_mcp,u_ccp, &
                v_ccm,v_pcm,v_cpm,v_cmc,v_pmc,v_mcc,v_ccc,v_pcc,v_mpc,v_cpc,v_cmp,v_mcp,v_ccp, &
                w_ccm,w_pcm,w_cpm,w_cmc,w_pmc,w_mcc,w_ccc,w_pcc,w_mpc,w_cpc,w_cmp,w_mcp,w_ccp, &
                flx_x_mcc,flx_x_ccc,flx_x_pcc,flx_x_mpc,flx_x_cpc,flx_x_mcp,flx_x_ccp, &
                flx_y_cmc,flx_y_ccc,flx_y_cpc,flx_y_pmc,flx_y_pcc,flx_y_cmp,flx_y_ccp, &
                flx_z_ccm,flx_z_ccc,flx_z_ccp,flx_z_pcm,flx_z_pcc,flx_z_cpm,flx_z_cpc, &
                dzci_c,dzci_m,dzfi_c,dzfi_p, &
                psixp,psiyp,psizp
    real(rp) :: uuip,uuim,vujp,vujm,wukp,wukm,uvip,uvim,vvjp,vvjm,wvkp,wvkm,uwip,uwim,vwjp,vwjm,wwkp,wwkm
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm,dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm,dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    real(rp) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
    real(rp) :: omxp,omxm,omyp,omym,omzp,omzm,domzdy,domydz,domxdz,domzdx,domydx,domxdy
    real(rp) :: dpsidx,dpsidy,dpsidz,dpsidxp,dpsidxm,dpsidyp,dpsidym,dpsidzp,dpsidzm
    real(rp) :: mux,muy,muz,muxp,muxm,muyp,muym,muzp,muzm,dmudx,dmudy,dmudz
    real(rp) :: rhoxp,rhoyp,rhozp
    real(rp) :: dudt_aux,dvdt_aux,dwdt_aux
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    mu  = mu12(2);  dmu  = mu12(1)-mu12(2)
    dxi = dli(1)
    dyi = dli(2)
    !
    !
    ! making an exception for this kernel -- private variables not explicitly mentioned for the sake of conciseness
    !                                        all scalars should be firstprivate/private
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          u_ccm = u(i  ,j  ,k-1)
          u_pcm = u(i+1,j  ,k-1)
          u_cpm = u(i  ,j+1,k-1)
          u_cmc = u(i  ,j-1,k  )
          u_pmc = u(i+1,j-1,k  )
          u_mcc = u(i-1,j  ,k  )
          u_ccc = u(i  ,j  ,k  )
          u_pcc = u(i+1,j  ,k  )
          u_mpc = u(i-1,j+1,k  )
          u_cpc = u(i  ,j+1,k  )
          u_cmp = u(i  ,j-1,k+1)
          u_mcp = u(i-1,j  ,k+1)
          u_ccp = u(i  ,j  ,k+1)
          !
          v_ccm = v(i  ,j  ,k-1)
          v_pcm = v(i+1,j  ,k-1)
          v_cpm = v(i  ,j+1,k-1)
          v_cmc = v(i  ,j-1,k  )
          v_pmc = v(i+1,j-1,k  )
          v_mcc = v(i-1,j  ,k  )
          v_ccc = v(i  ,j  ,k  )
          v_pcc = v(i+1,j  ,k  )
          v_mpc = v(i-1,j+1,k  )
          v_cpc = v(i  ,j+1,k  )
          v_cmp = v(i  ,j-1,k+1)
          v_mcp = v(i-1,j  ,k+1)
          v_ccp = v(i  ,j  ,k+1)
          !
          w_ccm = w(i  ,j  ,k-1)
          w_pcm = w(i+1,j  ,k-1)
          w_cpm = w(i  ,j+1,k-1)
          w_cmc = w(i  ,j-1,k  )
          w_pmc = w(i+1,j-1,k  )
          w_mcc = w(i-1,j  ,k  )
          w_ccc = w(i  ,j  ,k  )
          w_pcc = w(i+1,j  ,k  )
          w_mpc = w(i-1,j+1,k  )
          w_cpc = w(i  ,j+1,k  )
          w_cmp = w(i  ,j-1,k+1)
          w_mcp = w(i-1,j  ,k+1)
          w_ccp = w(i  ,j  ,k+1)
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
          flx_x_mcc = rho*u_mcc + drho*psiflx_x(i-1,j  ,k  )
          flx_x_ccc = rho*u_ccc + drho*psiflx_x(i  ,j  ,k  )
          flx_x_pcc = rho*u_pcc + drho*psiflx_x(i+1,j  ,k  )
          flx_x_mpc = rho*u_mpc + drho*psiflx_x(i-1,j+1,k  )
          flx_x_cpc = rho*u_cpc + drho*psiflx_x(i  ,j+1,k  )
          flx_x_mcp = rho*u_mcp + drho*psiflx_x(i-1,j  ,k+1)
          flx_x_ccp = rho*u_ccp + drho*psiflx_x(i  ,j  ,k+1)
          !
          flx_y_cmc = rho*v_cmc + drho*psiflx_y(i  ,j-1,k  )
          flx_y_ccc = rho*v_ccc + drho*psiflx_y(i  ,j  ,k  )
          flx_y_cpc = rho*v_cpc + drho*psiflx_y(i  ,j+1,k  )
          flx_y_pmc = rho*v_pmc + drho*psiflx_y(i+1,j-1,k  )
          flx_y_pcc = rho*v_pcc + drho*psiflx_y(i+1,j  ,k  )
          flx_y_cmp = rho*v_cmp + drho*psiflx_y(i  ,j-1,k+1)
          flx_y_ccp = rho*v_ccp + drho*psiflx_y(i  ,j  ,k+1)
          !
          flx_z_ccm = rho*w_ccm + drho*psiflx_z(i  ,j  ,k-1)
          flx_z_ccc = rho*w_ccc + drho*psiflx_z(i  ,j  ,k  )
          flx_z_ccp = rho*w_ccp + drho*psiflx_z(i  ,j  ,k+1)
          flx_z_pcm = rho*w_pcm + drho*psiflx_z(i+1,j  ,k-1)
          flx_z_pcc = rho*w_pcc + drho*psiflx_z(i+1,j  ,k  )
          flx_z_cpm = rho*w_cpm + drho*psiflx_z(i  ,j+1,k-1)
          flx_z_cpc = rho*w_cpc + drho*psiflx_z(i  ,j+1,k  )
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
          !
          ! mass--momentum-consistent flx
          !
          uuip = 0.25*(flx_x_pcc+flx_x_ccc)*(u_ccc+u_pcc)
          uuim = 0.25*(flx_x_mcc+flx_x_ccc)*(u_ccc+u_mcc)
          vujp = 0.25*(flx_y_pcc+flx_y_ccc)*(u_ccc+u_cpc)
          vujm = 0.25*(flx_y_pmc+flx_y_cmc)*(u_ccc+u_cmc)
          wukp = 0.25*(flx_z_pcc+flx_z_ccc)*(u_ccc+u_ccp)
          wukm = 0.25*(flx_z_pcm+flx_z_ccm)*(u_ccc+u_ccm)
          dudt_aux = dxi*( -uuip + uuim ) + dyi*( -vujp + vujm ) + dzfi_c*( -wukp + wukm )
          !
          uvip = 0.25*(flx_x_ccc+flx_x_cpc)*(v_ccc+v_pcc)
          uvim = 0.25*(flx_x_mcc+flx_x_mpc)*(v_ccc+v_mcc)
          vvjp = 0.25*(flx_y_ccc+flx_y_cpc)*(v_ccc+v_cpc)
          vvjm = 0.25*(flx_y_ccc+flx_y_cmc)*(v_ccc+v_cmc)
          wvkp = 0.25*(flx_z_ccc+flx_z_cpc)*(v_ccc+v_ccp)
          wvkm = 0.25*(flx_z_ccm+flx_z_cpm)*(v_ccc+v_ccm)
          dvdt_aux = dxi*( -uvip + uvim ) + dyi*( -vvjp + vvjm ) + dzfi_c*( -wvkp + wvkm )
          !
          uwip = 0.25*(flx_x_ccc+flx_x_ccp)*(w_ccc+w_pcc)
          uwim = 0.25*(flx_x_mcc+flx_x_mcp)*(w_ccc+w_mcc)
          vwjp = 0.25*(flx_y_ccc+flx_y_ccp)*(w_ccc+w_cpc)
          vwjm = 0.25*(flx_y_cmc+flx_y_cmp)*(w_ccc+w_cmc)
          wwkp = 0.25*(flx_z_ccc+flx_z_ccp)*(w_ccc+w_ccp)
          wwkm = 0.25*(flx_z_ccc+flx_z_ccm)*(w_ccc+w_ccm)
          dwdt_aux = dxi*( -uwip + uwim ) + dyi*( -vwjp + vwjm ) + dzci_c*( -wwkp + wwkm )
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
                                mux*(domzdy-domydz)) ! vorticity curl
          !dudt_aux = dudt_aux + (dmudx*(dudx+dudx) + dmudy*(dudy+dvdx) + dmudz*(dudz+dwdx) + &
          !                      mux*((dudxp-dudxm)*dxi+(dudyp-dudym)*dyi+(dudzp-dudzm)*dzfi_c)) ! laplacian
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
                                muy*(domxdz-domzdx)) ! vorticity curl
          !dvdt_aux = dvdt_aux + (dmudx*(dudy+dvdx) + dmudy*(dvdy+dvdy) + dmudz*(dvdz+dwdy) + &
          !                      muy*((dvdxp-dvdxm)*dxi+(dvdyp-dvdym)*dyi+(dvdzp-dvdzm)*dzfi_c)) ! laplacian
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
                                muz*(domydx-domxdy)) ! vorticity curl
          !dwdt_aux = dwdt_aux + (dmudx*(dudz+dwdx) + dmudy*(dvdz+dwdy) + dmudz*(dwdz+dwdz) + &
          !                      muz*((dwdxp-dwdxm)*dxi+(dwdyp-dwdym)*dyi+(dwdzp-dwdzm)*dzci_c)) ! laplacian
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
          dudt_aux = dudt_aux + dxi*(   (dudxp+dudxp)*muxp-(dudxm+dudxm)*muxm) + &
                                dyi*(   (dudyp+dvdxp)*muyp-(dudym+dvdxm)*muym) + &
                                dzfi_c*((dudzp+dwdxp)*muzp-(dudzm+dwdxm)*muzm)
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
          dvdt_aux = dvdt_aux + dxi*(   (dvdxp+dudyp)*muxp-(dvdxm+dudym)*muxm) + &
                                dyi*(   (dvdyp+dvdyp)*muyp-(dvdym+dvdym)*muym) + &
                                dzfi_c*((dvdzp+dwdyp)*muzp-(dvdzm+dwdym)*muzm)
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
          dwdt_aux = dwdt_aux + dxi*(   (dwdxp+dudzp)*muxp-(dwdxm+dudzm)*muxm) + &
                                dyi*(   (dwdyp+dvdzp)*muyp-(dwdym+dvdzm)*muym) + &
                                dzci_c*((dwdzp+dwdzp)*muzp-(dwdzm+dwdzm)*muzm)
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
  subroutine mom_xyz_oth(n,dli,dzci,dzfi,dt_r,rho12,beta12,bforce,gacc,sigma,rho0,rho_av, &
                         p,psi,kappa,s,pn,po,dudt,dvdt,dwdt)
    implicit none
    integer , intent(in   ), dimension(3) :: n
    real(rp), intent(in   ), dimension(3) :: dli
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi
    real(rp), intent(in   ) :: dt_r
    real(rp), intent(in   ), dimension(2) :: rho12,beta12
    real(rp), intent(in   ), dimension(3) :: bforce,gacc
    real(rp), intent(in   ) :: sigma
    real(rp), intent(in   ) :: rho0,rho_av(3)
    real(rp), intent(in   ), dimension(0:,0:,0:)           :: p,psi,kappa
    real(rp), intent(in   ), dimension(0:,0:,0:), optional :: s,pn,po
    real(rp), intent(inout), dimension( :, :, :)           :: dudt,dvdt,dwdt
    integer :: i,j,k
    real(rp) :: rho,drho,rhobeta,drhobeta
    real(rp) :: dxi,dyi
    real(rp) :: c_ccm,c_pcm,c_cpm,c_cmc,c_pmc,c_mcc,c_ccc,c_pcc,c_mpc,c_cpc,c_cmp,c_mcp,c_ccp,c_cpp,c_ppc,c_pcp, &
                p_ccc,p_pcc,p_cpc,p_ccp, &
                q_ccc,q_pcc,q_cpc,q_ccp, &
                s_ccc,s_pcc,s_cpc,s_ccp, &
                k_ccc,k_pcc,k_cpc,k_ccp, &
                bforcex,bforcey,bforcez,gaccx,gaccy,gaccz,rhox_av,rhoy_av,rhoz_av, &
                dzci_c,dzci_m,dzfi_c,dzfi_p, &
                psixp,psiyp,psizp
    real(rp) :: rhoxp,rhoyp,rhozp, &
                skappaxp,skappayp,skappazp, &
                factorxp,factoryp,factorzp, &
                dpdx  ,dpdy  ,dpdz  , &
                dpdx_e,dpdy_e,dpdz_e, &
                surfx  ,surfy  ,surfz
    real(rp) :: dudt_aux,dvdt_aux,dwdt_aux
    real(rp) :: surf_factor
    !
    rho     = rho12(2); drho = rho12(1)-rho12(2)
    rhobeta = rho12(2)*beta12(2); drhobeta = rho12(1)*beta12(1)- rho12(2)*beta12(2)
    dxi = dli(1)
    dyi = dli(2)
    surf_factor = sigma*2./(rho12(1)+rho12(2))
    rhox_av = rho_av(1); rhoy_av = rho_av(2); rhoz_av = rho_av(3)
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
          !
#if defined(_CONSTANT_COEFFS_POISSON)
          q_ccc = (1.+dt_r)*pn(i  ,j  ,k  )-dt_r*po(i  ,j  ,k  )
          q_pcc = (1.+dt_r)*pn(i+1,j  ,k  )-dt_r*po(i+1,j  ,k  )
          q_cpc = (1.+dt_r)*pn(i  ,j+1,k  )-dt_r*po(i  ,j+1,k  )
          q_ccp = (1.+dt_r)*pn(i  ,j  ,k+1)-dt_r*po(i  ,j  ,k+1)
#endif
          !
#if defined(_SCALAR) && defined(_BOUSSINESQ_BUOYANCY)
          s_ccc = s(i  ,j  ,k  )
          s_pcc = s(i+1,j  ,k  )
          s_cpc = s(i  ,j+1,k  )
          s_ccp = s(i  ,j  ,k+1)
#endif
          !
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
          dudt_aux = bforcex + (rhoxp-rhox_av)*gaccx
          dvdt_aux = bforcey + (rhoyp-rhoy_av)*gaccy
          dwdt_aux = bforcez + (rhozp-rhoz_av)*gaccz
          !
          skappaxp = 0.5*(k_pcc+k_ccc)*surf_factor
          skappayp = 0.5*(k_cpc+k_ccc)*surf_factor
          skappazp = 0.5*(k_ccp+k_ccc)*surf_factor
          surfx = skappaxp*rhoxp*(c_pcc-c_ccc)*dxi
          surfy = skappayp*rhoyp*(c_cpc-c_ccc)*dyi
          surfz = skappazp*rhozp*(c_ccp-c_ccc)*dzci_c
          dudt_aux = dudt_aux + surfx
          dvdt_aux = dvdt_aux + surfy
          dwdt_aux = dwdt_aux + surfz
#if defined(_CONSTANT_COEFFS_POISSON)
          dpdx_e = (q_pcc-q_ccc)*dxi
          dpdy_e = (q_cpc-q_ccc)*dyi
          dpdz_e = (q_ccp-q_ccc)*dzci_c
          dudt_aux = dudt_aux - (dpdx_e + rhoxp/rho0*(dpdx-dpdx_e))
          dvdt_aux = dvdt_aux - (dpdy_e + rhoyp/rho0*(dpdy-dpdy_e))
          dwdt_aux = dwdt_aux - (dpdz_e + rhozp/rho0*(dpdz-dpdz_e))
#else
          dudt_aux = dudt_aux - dpdx
          dvdt_aux = dvdt_aux - dpdy
          dwdt_aux = dwdt_aux - dpdz
#endif
          !
          ! buoyancy
          !
#if defined(_SCALAR) && defined(_BOUSSINESQ_BUOYANCY)
          factorxp = rhobeta + drhobeta*psixp
          factoryp = rhobeta + drhobeta*psiyp
          factorzp = rhobeta + drhobeta*psizp
          dudt_aux = dudt_aux - gaccx*factorxp*0.5*(s_pcc+s_ccc)
          dvdt_aux = dvdt_aux - gaccy*factoryp*0.5*(s_cpc+s_ccc)
          dwdt_aux = dwdt_aux - gaccz*factorzp*0.5*(s_ccp+s_ccc)
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
    real(rp), intent(in ), dimension(0:)    :: dzci,dzfi
    real(rp), intent(in ), dimension(2)     :: mu12
    real(rp), intent(in ), dimension(0:,0:,0:) :: psi,u,v,w
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
