! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_mom
  use mpi
  use mod_types
  implicit none
  private
  public momx_a
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
    real(rp) :: uuip,uuim,uvjp,uvjm,uwkp,uwkm
    !
    !$acc parallel loop collapse(3) default(present) private(uuip,uuim,uvjp,uvjm,uwkp,uwkm) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(uuip,uuim,uvjp,uvjm,uwkp,uwkm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,u,v,w,dudt,dzfi)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          uuip  = 0.25*( u(i+1,j,k)+u(i,j,k) )*( u(i+1,j  ,k  )+u(i,j  ,k  ) )
          uuim  = 0.25*( u(i-1,j,k)+u(i,j,k) )*( u(i-1,j  ,k  )+u(i,j  ,k  ) )
          uvjp  = 0.25*( u(i,j+1,k)+u(i,j,k) )*( v(i+1,j  ,k  )+v(i,j  ,k  ) )
          uvjm  = 0.25*( u(i,j-1,k)+u(i,j,k) )*( v(i+1,j-1,k  )+v(i,j-1,k  ) )
          uwkp  = 0.25*( u(i,j,k+1)+u(i,j,k) )*( w(i+1,j  ,k  )+w(i,j  ,k  ) )
          uwkm  = 0.25*( u(i,j,k-1)+u(i,j,k) )*( w(i+1,j  ,k-1)+w(i,j  ,k-1) )
          !
          dudt(i,j,k) = dudt(i,j,k) + &
                        dxi*(     -uuip + uuim ) + &
                        dyi*(     -uvjp + uvjm ) + &
                        dzfi(k)*( -uwkp + uwkm )
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
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(uvip,uvim,vvjp,vvjm,wvkp,wvkm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzfi,u,v,w,dvdt)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          uvip  = 0.25*( u(i  ,j,k  )+u(i  ,j+1,k  ) )*( v(i,j,k)+v(i+1,j,k) )
          uvim  = 0.25*( u(i-1,j,k  )+u(i-1,j+1,k  ) )*( v(i,j,k)+v(i-1,j,k) )
          vvjp  = 0.25*( v(i  ,j,k  )+v(i  ,j+1,k  ) )*( v(i,j,k)+v(i,j+1,k) )
          vvjm  = 0.25*( v(i  ,j,k  )+v(i  ,j-1,k  ) )*( v(i,j,k)+v(i,j-1,k) )
          wvkp  = 0.25*( w(i  ,j,k  )+w(i  ,j+1,k  ) )*( v(i,j,k+1)+v(i,j,k) )
          wvkm  = 0.25*( w(i  ,j,k-1)+w(i  ,j+1,k-1) )*( v(i,j,k-1)+v(i,j,k) )
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
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(uwip,uwim,vwjp,vwjm,wwkp,wwkm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzci,u,v,w,dwdt)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          uwip  = 0.25*( w(i,j,k)+w(i+1,j,k) )*( u(i  ,j  ,k)+u(i  ,j  ,k+1) )
          uwim  = 0.25*( w(i,j,k)+w(i-1,j,k) )*( u(i-1,j  ,k)+u(i-1,j  ,k+1) )
          vwjp  = 0.25*( w(i,j,k)+w(i,j+1,k) )*( v(i  ,j  ,k)+v(i  ,j  ,k+1) )
          vwjm  = 0.25*( w(i,j,k)+w(i,j-1,k) )*( v(i  ,j-1,k)+v(i  ,j-1,k+1) )
          wwkp  = 0.25*( w(i,j,k)+w(i,j,k+1) )*( w(i  ,j  ,k)+w(i  ,j  ,k+1) )
          wwkm  = 0.25*( w(i,j,k)+w(i,j,k-1) )*( w(i  ,j  ,k)+w(i  ,j  ,k-1) )
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
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), intent(in), dimension(2) :: rho12,mu12
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm, &
                dvdxp,dvdxm,dwdxp,dwdxm
    real(rp) :: muxp,muxm,muyp,muym,muzp,muzm,rhop
    integer :: i,j,k
    real(rp) :: rho1,rho2,mu1,mu2
    !
    rho1 = rho12(1); rho2 = rho12(2)
    mu1  = mu12(1) ; mu2  = mu12(2)
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(dudxp,dudxm,dudyp,dudym,dudzp,dudzm,dvdxp,dvdxm,dwdxp,dwdxm) &
    !$acc private(muxp,muxm,muyp,muym,muzp,muzm,rhop) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(dudxp,dudxm,dudyp,dudym,dudzp,dudzm,dvdxp,dvdxm,dwdxp,dwdxm) &
    !$OMP PRIVATE(muxp,muxm,muyp,muym,muzp,muzm,rhop) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzci,dzfi,rho1,rho2,mu1,mu2,psi,u,v,w,dudt)
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
          muxp = mu2 + (mu1-mu2)*psi(i+1,j,k)
          muxm = mu2 + (mu1-mu2)*psi(i  ,j,k)
          muyp = mu2 + (mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i+1,j+1,k)+psi(i+1,j,k))
          muym = mu2 + (mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j-1,k)+psi(i+1,j-1,k)+psi(i+1,j,k))
          muzp = mu2 + (mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i+1,j,k+1)+psi(i+1,j,k))
          muzm = mu2 + (mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j,k-1)+psi(i+1,j,k-1)+psi(i+1,j,k))
          rhop = rho2 + (rho1-rho2)*0.5*(psi(i+1,j,k)+psi(i,j,k))
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
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), intent(in), dimension(2) :: rho12,mu12
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    real(rp) :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm, &
                dudyp,dudym,dwdyp,dwdym
    real(rp) :: muxp,muxm,muyp,muym,muzp,muzm,rhop
    integer :: i,j,k
    real(rp) :: rho1,rho2,mu1,mu2
    !
    rho1 = rho12(1); rho2 = rho12(2)
    mu1  = mu12(1) ; mu2  = mu12(2)
    !
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm,dudyp,dudym,dwdyp,dwdym) &
    !$acc private(muxp,muxm,muyp,muym,muzp,muzm,rhop) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm,dudyp,dudym,dwdyp,dwdym) &
    !$OMP PRIVATE(muxp,muxm,muyp,muym,muzp,muzm,rhop) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzci,dzfi,rho1,rho2,mu1,mu2,psi,u,v,w,dvdt)
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
          muxp = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i+1,j,k)+psi(i+1,j+1,k)+psi(i,j+1,k))
          muxm = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i-1,j,k)+psi(i-1,j+1,k)+psi(i,j+1,k))
          muyp = mu2+(mu1-mu2)*psi(i,j+1,k)
          muym = mu2+(mu1-mu2)*psi(i,j  ,k)
          muzp = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i,j+1,k+1)+psi(i,j,k+1))
          muzm = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i,j+1,k-1)+psi(i,j,k-1))
          rhop = rho2+(rho1-rho2)*0.5*(psi(i,j+1,k)+psi(i,j,k))
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
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), intent(in), dimension(2) :: rho12,mu12
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    integer :: i,j,k
    real(rp) :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm, &
                dudzp,dudzm,dvdzp,dvdzm
    real(rp) :: muxp,muxm,muyp,muym,muzp,muzm,rhop
    real(rp) :: rho1,rho2,mu1,mu2
    !
    rho1 = rho12(1); rho2 = rho12(2)
    mu1  = mu12(1) ; mu2  = mu12(2)
    !
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm,dudzp,dudzm,dvdzp,dvdzm) &
    !$acc private(muxp,muxm,muyp,muym,muzp,muzm,rhop) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm,dudzp,dudzm,dvdzp,dvdzm) &
    !$OMP PRIVATE(muxp,muxm,muyp,muym,muzp,muzm,rhop) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzci,dzfi,rho1,rho2,mu1,mu2,psi,u,v,w,dwdt)
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
          muxp = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i+1,j ,k+1)+psi(i+1,j ,k) )
          muxm = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i-1,j ,k+1)+psi(i-1,j ,k) )
          muyp = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i ,j+1,k+1)+psi(i ,j+1,k) )
          muym = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i ,j-1,k+1)+psi(i ,j-1,k) )
          muzp = mu2+(mu1-mu2)*psi(i,j,k+1)
          muzm = mu2+(mu1-mu2)*psi(i,j,k  )
          rhop = rho2+(rho1-rho2)*0.5*(psi(i,j,k+1)+psi(i,j,k))
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
  subroutine momx_p(nx,ny,nz,dxi,bforce,gacc,rho0,rho_av,rho12,psi,p,pp,dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi
    real(rp), intent(in) :: bforce,gacc,rho0,rho_av,rho12(2)
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,p,pp
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    real(rp) :: rhop,dpdl
    integer :: i,j,k
    real(rp) :: rho1,rho2
    !
    rho1 = rho12(1); rho2 = rho12(2)
    !
    !$acc parallel loop collapse(3) default(present) private(rhop,dpdl) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(nx,ny,nz,dxi,bforce,gacc,rho0,rho_av,rho1,rho2,psi,p,dudt) &
    !$OMP PRIVATE(rhop,dpdl)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          rhop = rho2+(rho1-rho2)*0.5*(psi(i+1,j,k)+psi(i,j,k))
          dpdl = - (p(i+1,j,k)-p(i,j,k))*dxi + bforce
          dudt(i,j,k) = dudt(i,j,k) + &
#if defined(_CONSTANT_COEFFS_POISSON)
                        dpdl/rho0 + (1./rhop-1./rho0)*(pp(i+1,j,k)-pp(i,j,k))*dxi + gacc
#else                                 
                        dpdl/rhop + gacc*(1.-rho_av/rhop)
#endif
        end do
      end do
    end do
  end subroutine momx_p
  !
  subroutine momy_p(nx,ny,nz,dyi,bforce,gacc,rho0,rho_av,rho12,psi,p,pp,dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dyi
    real(rp), intent(in) :: bforce,gacc,rho0,rho_av,rho12(2)
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,p,pp
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    integer :: i,j,k
    real(rp) :: rhop,dpdl
    real(rp) :: rho1,rho2
    !
    rho1 = rho12(1); rho2 = rho12(2)
    !
    !$acc parallel loop collapse(3) default(present) private(rhop,dpdl) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(nx,ny,nz,dyi,bforce,gacc,rho0,rho_av,rho1,rho2,psi,p,dvdt) &
    !$OMP PRIVATE(rhop,dpdl)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          rhop = rho2+(rho1-rho2)*0.5*(psi(i,j+1,k)+psi(i,j,k))
          dpdl = - (p(i,j+1,k)-p(i,j,k))*dyi + bforce
          dvdt(i,j,k) = dvdt(i,j,k) + &
#if defined(_CONSTANT_COEFFS_POISSON)
                        dpdl/rho0 + (1./rhop-1./rho0)*(pp(i,j+1,k)-pp(i,j,k))*dyi + gacc
#else                                 
                        dpdl/rhop + gacc*(1.-rho_av/rhop)
#endif
        end do
      end do
    end do
  end subroutine momy_p
  !
  subroutine momz_p(nx,ny,nz,dzci,bforce,gacc,rho0,rho_av,rho12,psi,p,pp,dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in), dimension(0:) :: dzci
    real(rp), intent(in) :: bforce,gacc,rho0,rho_av,rho12(2)
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,p,pp
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    real(rp) :: rhop,dpdl
    integer :: i,j,k
    real(rp) :: rho1,rho2
    !
    rho1 = rho12(1); rho2 = rho12(2)
    !
    !$acc parallel loop collapse(3) default(present) private(rhop,dpdl) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(nx,ny,nz,dzci,bforce,gacc,rho0,rho_av,rho1,rho2,psi,p,dwdt) &
    !$OMP PRIVATE(rhop,dpdl)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          rhop = rho2+(rho1-rho2)*0.5*(psi(i,j,k+1)+psi(i,j,k))
          dpdl = - (p(i,j,k+1)-p(i,j,k))*dzci(k) + bforce
          dwdt(i,j,k) = dwdt(i,j,k) + &
#if defined(_CONSTANT_COEFFS_POISSON)
                        dpdl/rho0 + (1./rhop-1./rho0)*(pp(i,j,k+1)-pp(i,j,k))*dzci(k) + gacc
#else
                        dpdl/rhop + gacc*(1.-rho_av/rhop)
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
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,kappa
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    real(rp) :: rhop,kappap
    integer :: i,j,k
    real(rp) :: rho1,rho2
    !
    rho1 = rho12(1); rho2 = rho12(2)
    !
    !$acc parallel loop collapse(3) default(present) private(rhop,kappap) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(nx,ny,nz,dxi,sigma,rho1,rho2,psi,kappa,dudt) &
    !$OMP PRIVATE(rhop,kappap)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          rhop   = rho2+(rho1-rho2)*0.5*(psi(i+1,j,k)+psi(i,j,k))
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
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,kappa
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    real(rp) :: rhop,kappap
    integer :: i,j,k
    real(rp) :: rho1,rho2
    !
    rho1 = rho12(1); rho2 = rho12(2)
    !
    !$acc parallel loop collapse(3) default(present) private(rhop,kappap) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(nx,ny,nz,dyi,sigma,rho1,rho2,psi,kappa,dvdt) &
    !$OMP PRIVATE(rhop,kappap)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          rhop   = rho2+(rho1-rho2)*0.5*(psi(i,j+1,k)+psi(i,j,k))
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
    real(rp), intent(in), dimension(0:) :: dzci
    real(rp), intent(in) :: sigma,rho12(2)
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,kappa
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    real(rp) :: rhop,kappap
    integer :: i,j,k
    real(rp) :: rho1,rho2
    !
    rho1 = rho12(1); rho2 = rho12(2)
    !
    !$acc parallel loop collapse(3) default(present) private(rhop,kappap) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(nx,ny,nz,dzci,sigma,rho1,rho2,psi,kappa,dwdt) &
    !$OMP PRIVATE(rhop,kappap)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          rhop   = rho2+(rho1-rho2)*0.5*(psi(i,j,k+1)+psi(i,j,k))
          kappap = 0.5*(kappa(i,j,k+1)+kappa(i,j,k))
          !
          dwdt(i,j,k) = dwdt(i,j,k) + &
                        sigma*kappap*(psi(i,j,k+1)-psi(i,j,k))*dzci(k)/rhop
        end do
      end do
    end do
  end subroutine momz_sigma
  !
  subroutine momx_buoy(nx,ny,nz,gacc,rhobeta12,psi,s,dudt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: gacc,rhobeta12(2)
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,s
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    real(rp) :: factorp
    integer :: i,j,k
    real(rp) :: rhobeta1,rhobeta2
    !
    rhobeta1 = rhobeta12(1); rhobeta2 = rhobeta12(2)
    !
    !$acc parallel loop collapse(3) default(present) private(factorp) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(nx,ny,nz,gacc,rhobeta1,rhobeta2,psi,s,dudt) &
    !$OMP PRIVATE(factorp)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          factorp = rhobeta2-(rhobeta1-rhobeta2)*0.5*(psi(i+1,j,k)+psi(i,j,k))
          !
          dudt(i,j,k) = dudt(i,j,k) + &
                        gacc*factorp*0.5*(s(i+1,j,k)+s(i,j,k))
        end do
      end do
    end do
  end subroutine momx_buoy
  !
  subroutine momy_buoy(nx,ny,nz,gacc,rhobeta12,psi,s,dvdt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: gacc,rhobeta12(2)
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,s
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    real(rp) :: factorp
    integer :: i,j,k
    real(rp) :: rhobeta1,rhobeta2
    !
    rhobeta1 = rhobeta12(1); rhobeta2 = rhobeta12(2)
    !
    !$acc parallel loop collapse(3) default(present) private(factorp) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(nx,ny,nz,gacc,rhobeta1,rhobeta2,psi,s,dvdt) &
    !$OMP PRIVATE(factorp)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          factorp = rhobeta2-(rhobeta1-rhobeta2)*0.5*(psi(i,j+1,k)+psi(i,j,k))
          !
          dvdt(i,j,k) = dvdt(i,j,k) + &
                        gacc*factorp*0.5*(s(i,j+1,k)+s(i,j,k))
        end do
      end do
    end do
  end subroutine momy_buoy
  !
  subroutine momz_buoy(nx,ny,nz,gacc,rhobeta12,psi,s,dwdt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: gacc,rhobeta12(2)
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi,s
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    real(rp) :: factorp
    integer :: i,j,k
    real(rp) :: rhobeta1,rhobeta2
    !
    rhobeta1 = rhobeta12(1); rhobeta2 = rhobeta12(2)
    !
    !$acc parallel loop collapse(3) default(present) private(factorp) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(nx,ny,nz,gacc,rhobeta1,rhobeta2,psi,s,dwdt) &
    !$OMP PRIVATE(factorp)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          factorp = rhobeta2-(rhobeta1-rhobeta2)*0.5*(psi(i,j,k+1)+psi(i,j,k))
          !
          dwdt(i,j,k) = dwdt(i,j,k) + &
                        gacc*factorp*0.5*(s(i,j,k+1)+s(i,j,k))
        end do
      end do
    end do
  end subroutine momz_buoy
  !
#if 1
  subroutine cmpt_wallshear(n,is_bound,l,dli,dzci,dzfi,mu12,psi,u,v,w,taux,tauy,tauz)
    !
    ! n.b.: is_bound should exclude periodic BCs
    !
    implicit none
    integer , intent(in ), dimension(3) :: n
    logical , intent(in ), dimension(0:1,3) :: is_bound
    real(rp), intent(in ), dimension(3)     :: l,dli
    real(rp), intent(in ), dimension(0:)    :: dzci,dzfi
    real(rp), intent(in )                   :: mu12(2)
    real(rp), intent(in ), dimension(0:,0:,0:) :: psi,u,v,w
    real(rp), intent(out), dimension(3) :: taux,tauy,tauz
    real(rp) :: dudyp,dudym,dudzp,dudzm, &
                dvdxp,dvdxm,dvdzp,dvdzm, &
                dwdxp,dwdxm,dwdyp,dwdym, &
                muxp,muxm,muyp,muym,muzp,muzm
    !
    real(rp) :: tau21,tau31,tau12,tau32,tau13,tau23
    integer :: i,j,k,nx,ny,nz
    real(rp) :: dxi,dyi,lx,ly,lz,mu1,mu2
    real(rp) :: tau(3,3)
    integer :: ierr
    !
    mu1  = mu12(1) ; mu2  = mu12(2)
    nx = n(1); ny = n(2); nz = n(3)
    dxi = dli(1); dyi = dli(2)
    lx = l(1); ly = l(2); lz = l(3)
    tau21 = 0._rp
    !$acc data copy(tau21) async(1)
    if(is_bound(0,2)) then
      j = 1
      !$acc parallel loop collapse(2) default(present) private(dvdxp,dvdxm,dudyp,dudym,muyp,muym) &
      !$acc reduction(+:tau21) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(dvdxp,dvdxm,dudyp,dudym,muyp,muym) REDUCTION(+:tau21)
      do k=1,nz
        do i=1,nx
          dvdxp = (v(i+1,j  ,k  )-v(i  ,j  ,k  ))*dxi
          dvdxm = (v(i+1,j-1,k  )-v(i  ,j-1,k  ))*dxi
          dudyp = (u(i  ,j+1,k  )-u(i  ,j  ,k  ))*dyi
          dudym = (u(i  ,j  ,k  )-u(i  ,j-1,k  ))*dyi
          !
          muyp = mu2 + (mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i+1,j+1,k)+psi(i+1,j,k))
          muym = mu2 + (mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j-1,k)+psi(i+1,j-1,k)+psi(i+1,j,k))
          !
          tau21 = tau21 + ((dudyp+dvdxp)*muyp-(dudym+dvdxm)*muym)/(dxi*dzfi(k)*lx*lz)
        end do
      end do
    end if
    if(is_bound(1,2)) then
      j = ny+1
      !$acc parallel loop collapse(2) default(present) private(dvdxp,dvdxm,dudyp,dudym,muyp,muym) &
      !$acc reduction(+:tau21) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(dvdxp,dvdxm,dudyp,dudym,muyp,muym) REDUCTION(+:tau21)
      do k=1,nz
        do i=1,nx
          dvdxp = (v(i+1,j  ,k  )-v(i  ,j  ,k  ))*dxi
          dvdxm = (v(i+1,j-1,k  )-v(i  ,j-1,k  ))*dxi
          dudyp = (u(i  ,j+1,k  )-u(i  ,j  ,k  ))*dyi
          dudym = (u(i  ,j  ,k  )-u(i  ,j-1,k  ))*dyi
          !
          muyp = mu2 + (mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i+1,j+1,k)+psi(i+1,j,k))
          muym = mu2 + (mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j-1,k)+psi(i+1,j-1,k)+psi(i+1,j,k))
          !
          tau21 = tau21 - ((dudyp+dvdxp)*muyp-(dudym+dvdxm)*muym)/(dxi*dzfi(k)*lx*lz)
        end do
      end do
    end if
    !$acc end data
    tau31 = 0._rp
    !$acc data copy(tau31) async(1)
    if(is_bound(0,3)) then
      k = 1
      !$acc parallel loop collapse(2) default(present) private(dudzp,dudzm,dwdxp,dwdxm,muzp,muzm) &
      !$acc reduction(+:tau31) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(dudzp,dudzm,dwdxp,dwdxm,muzp,muzm) REDUCTION(+:tau31)
      do j=1,ny
        do i=1,nx
          dudzp = (u(i  ,j  ,k+1)-u(i  ,j  ,k  ))*dzci(k  )
          dudzm = (u(i  ,j  ,k  )-u(i  ,j  ,k-1))*dzci(k-1)
          dwdxp = (w(i+1,j  ,k  )-w(i  ,j  ,k  ))*dxi
          dwdxm = (w(i+1,j  ,k-1)-w(i  ,j  ,k-1))*dxi
          !
          muzp = mu2 + (mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i+1,j,k+1)+psi(i+1,j,k))
          muzm = mu2 + (mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j,k-1)+psi(i+1,j,k-1)+psi(i+1,j,k))
          !
          tau31 = tau31 + ((dudzp+dwdxp)*muzp-(dudzm+dwdxm)*muzm)/(dxi*dyi*lx*ly)
        end do
      end do
    end if
    if(is_bound(1,3)) then
      k = nz+1
      !$acc parallel loop collapse(2) default(present) private(dudzp,dudzm,dwdxp,dwdxm,muzp,muzm) &
      !$acc reduction(+:tau31) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(dudzp,dudzm,dwdxp,dwdxm,muzp,muzm) REDUCTION(+:tau31)
      do j=1,ny
        do i=1,nx
          dudzp = (u(i  ,j  ,k+1)-u(i  ,j  ,k  ))*dzci(k  )
          dudzm = (u(i  ,j  ,k  )-u(i  ,j  ,k-1))*dzci(k-1)
          dwdxp = (w(i+1,j  ,k  )-w(i  ,j  ,k  ))*dxi
          dwdxm = (w(i+1,j  ,k-1)-w(i  ,j  ,k-1))*dxi
          !
          muzp = mu2 + (mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i+1,j,k+1)+psi(i+1,j,k))
          muzm = mu2 + (mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j,k-1)+psi(i+1,j,k-1)+psi(i+1,j,k))
          !
          tau31 = tau31 - ((dudzp+dwdxp)*muzp-(dudzm+dwdxm)*muzm)/(dxi*dyi*lx*ly)
        end do
      end do
    end if
    !$acc end data
    !
    tau12 = 0._rp
    !$acc data copy(tau12) async(1)
    if(is_bound(0,1)) then
      i = 1
      !$acc parallel loop collapse(2) default(present) private(dvdxp,dvdxm,dudyp,dudym,muxp,muxm) &
      !$acc reduction(+:tau12) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(dvdxp,dvdxm,dudyp,dudym,muxp,muxm) REDUCTION(+:tau12)
      do k=1,nz
        do j=1,ny
          dvdxp = (v(i+1,j  ,k  )-v(i  ,j  ,k  ))*dxi
          dvdxm = (v(i  ,j  ,k  )-v(i-1,j  ,k  ))*dxi
          dudyp = (u(i  ,j+1,k  )-u(i  ,j  ,k  ))*dyi
          dudym = (u(i-1,j+1,k  )-u(i-1,j  ,k  ))*dyi
          !
          muxp = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i+1,j,k)+psi(i+1,j+1,k)+psi(i,j+1,k))
          muxm = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i-1,j,k)+psi(i-1,j+1,k)+psi(i,j+1,k))
          !
          tau12 = tau12 + ((dvdxp+dudyp)*muxp-(dvdxm+dudym)*muxm)/(dyi*dzfi(k)*ly*lz)
        end do
      end do
    end if
    if(is_bound(1,1)) then
      i = nx+1
      !$acc parallel loop collapse(2) default(present) private(dvdxp,dvdxm,dudyp,dudym,muxp,muxm) &
      !$acc reduction(+:tau12) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(dvdxp,dvdxm,dudyp,dudym,muxp,muxm) REDUCTION(+:tau12)
      do k=1,nz
        do j=1,ny
          dvdxp = (v(i+1,j  ,k  )-v(i  ,j  ,k  ))*dxi
          dvdxm = (v(i  ,j  ,k  )-v(i-1,j  ,k  ))*dxi
          dudyp = (u(i  ,j+1,k  )-u(i  ,j  ,k  ))*dyi
          dudym = (u(i-1,j+1,k  )-u(i-1,j  ,k  ))*dyi
          !
          muxp = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i+1,j,k)+psi(i+1,j+1,k)+psi(i,j+1,k))
          muxm = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i-1,j,k)+psi(i-1,j+1,k)+psi(i,j+1,k))
          !
          tau12 = tau12 - ((dvdxp+dudyp)*muxp-(dvdxm+dudym)*muxm)/(dyi*dzfi(k)*ly*lz)
        end do
      end do
    end if
    !$acc end data
    tau32 = 0._rp
    !$acc data copy(tau32) async(1)
    if(is_bound(0,3)) then
      k = 1
      !$acc parallel loop collapse(2) default(present) private(dvdzp,dvdzm,dwdyp,dwdym,muzp,muzm) &
      !$acc reduction(+:tau32) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(dvdzp,dvdzm,dwdyp,dwdym,muzp,muzm) REDUCTION(+:tau32)
      do j=1,ny
        do i=1,nx
          dvdzp = (v(i  ,j  ,k+1)-v(i  ,j  ,k  ))*dzci(k  )
          dvdzm = (v(i  ,j  ,k  )-v(i  ,j  ,k-1))*dzci(k-1)
          dwdyp = (w(i  ,j+1,k  )-w(i  ,j  ,k  ))*dyi
          dwdym = (w(i  ,j+1,k-1)-w(i  ,j  ,k-1))*dyi
          !
          muzp = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i,j+1,k+1)+psi(i,j,k+1))
          muzm = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i,j+1,k-1)+psi(i,j,k-1))
          !
          tau32 = tau32 + ((dvdzp+dwdyp)*muzp-(dvdzm+dwdym)*muzm)/(dxi*dyi*lx*ly)
        end do
      end do
    end if
    if(is_bound(1,3)) then
      k = nz+1
      !$acc parallel loop collapse(2) default(present) private(dvdzp,dvdzm,dwdyp,dwdym,muzp,muzm) &
      !$acc reduction(+:tau32) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(dvdzp,dvdzm,dwdyp,dwdym,muzp,muzm) REDUCTION(+:tau32)
      do j=1,ny
        do i=1,nx
          dvdzp = (v(i  ,j  ,k+1)-v(i  ,j  ,k  ))*dzci(k  )
          dvdzm = (v(i  ,j  ,k  )-v(i  ,j  ,k-1))*dzci(k-1)
          dwdyp = (w(i  ,j+1,k  )-w(i  ,j  ,k  ))*dyi
          dwdym = (w(i  ,j+1,k-1)-w(i  ,j  ,k-1))*dyi
          !
          muzp = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i,j+1,k+1)+psi(i,j,k+1))
          muzm = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i,j+1,k-1)+psi(i,j,k-1))
          !
          tau32 = tau32 - ((dvdzp+dwdyp)*muzp-(dvdzm+dwdym)*muzm)/(dxi*dyi*lx*ly)
        end do
      end do
    end if
    !$acc end data
    !
    tau13 = 0._rp
    !$acc data copy(tau13) async(1)
    if(is_bound(0,1)) then
      i = 1
      !$acc parallel loop collapse(2) default(present) private(dwdxp,dwdxm,dudzp,dudzm,muxp,muxm) &
      !$acc reduction(+:tau13) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(dwdxp,dwdxm,dudzp,dudzm,muxp,muxm) REDUCTION(+:tau13)
      do k=1,nz
        do j=1,ny
          dwdxp = (w(i+1,j  ,k  )-w(i  ,j  ,k))*dxi
          dwdxm = (w(i  ,j  ,k  )-w(i-1,j  ,k))*dxi
          dudzp = (u(i  ,j  ,k+1)-u(i  ,j  ,k))*dzci(k  )
          dudzm = (u(i-1,j  ,k+1)-u(i-1,j  ,k))*dzci(k  )
          !
          muxp = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i+1,j ,k+1)+psi(i+1,j ,k) )
          muxm = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i-1,j ,k+1)+psi(i-1,j ,k) )
          !
          tau13 = tau13 + ((dwdxp+dudzp)*muxp-(dwdxm+dudzm)*muxm)/(dyi*dzfi(k)*ly*lz)
        end do
      end do
    end if
    if(is_bound(1,1)) then
      i = nx+1
      !$acc parallel loop collapse(2) default(present) private(dwdxp,dwdxm,dudzp,dudzm,muxp,muxm) &
      !$acc reduction(+:tau13) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(dwdxp,dwdxm,dudzp,dudzm,muxp,muxm) REDUCTION(+:tau13)
      do k=1,nz
        do j=1,ny
          dwdxp = (w(i+1,j  ,k  )-w(i  ,j  ,k))*dxi
          dwdxm = (w(i  ,j  ,k  )-w(i-1,j  ,k))*dxi
          dudzp = (u(i  ,j  ,k+1)-u(i  ,j  ,k))*dzci(k  )
          dudzm = (u(i-1,j  ,k+1)-u(i-1,j  ,k))*dzci(k  )
          !
          muxp = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i+1,j ,k+1)+psi(i+1,j ,k) )
          muxm = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i-1,j ,k+1)+psi(i-1,j ,k) )
          !
          tau13 = tau13 - ((dwdxp+dudzp)*muxp-(dwdxm+dudzm)*muxm)/(dyi*dzfi(k)*ly*lz)
        end do
      end do
    end if
    !$acc end data
    tau23 = 0._rp
    !$acc data copy(tau23) async(1)
    if(is_bound(0,2)) then
      j = 1
      !$acc parallel loop collapse(2) default(present) private(dwdyp,dwdym,dvdzp,dvdzm,muyp,muym) &
      !$acc reduction(+:tau23) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(dwdyp,dwdym,dvdzp,dvdzm,muyp,muym) REDUCTION(+:tau32)
      do k=1,nz
        do i=1,nx
          dwdyp = (w(i  ,j+1,k  )-w(i  ,j  ,k))*dyi
          dwdym = (w(i  ,j  ,k  )-w(i  ,j-1,k))*dyi
          dvdzp = (v(i  ,j  ,k+1)-v(i  ,j  ,k))*dzci(k  )
          dvdzm = (v(i  ,j-1,k+1)-v(i  ,j-1,k))*dzci(k  )
          !
          muyp = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i ,j+1,k+1)+psi(i ,j+1,k) )
          muym = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i ,j-1,k+1)+psi(i ,j-1,k) )
          !
          tau23 = tau23 + ((dwdyp+dvdzp)*muyp-(dwdym+dvdzm)*muym)/(dxi*dzfi(k)*lx*lz)
        end do
      end do
    end if
    if(is_bound(1,2)) then
      j = ny+1
      !$acc parallel loop collapse(2) default(present) private(dwdyp,dwdym,dvdzp,dvdzm,muyp,muym) &
      !$acc reduction(+:tau23) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(dwdyp,dwdym,dvdzp,dvdzm,muyp,muym) REDUCTION(+:tau32)
      do k=1,nz
        do i=1,nx
          dwdyp = (w(i  ,j+1,k  )-w(i  ,j  ,k))*dyi
          dwdym = (w(i  ,j  ,k  )-w(i  ,j-1,k))*dyi
          dvdzp = (v(i  ,j  ,k+1)-v(i  ,j  ,k))*dzci(k  )
          dvdzm = (v(i  ,j-1,k+1)-v(i  ,j-1,k))*dzci(k  )
          !
          muyp = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i ,j+1,k+1)+psi(i ,j+1,k) )
          muym = mu2+(mu1-mu2)*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i ,j-1,k+1)+psi(i ,j-1,k) )
          !
          tau23 = tau23 - ((dwdyp+dvdzp)*muyp-(dwdym+dvdzm)*muym)/(dxi*dzfi(k)*lx*lz)
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
#endif
end module mod_mom
