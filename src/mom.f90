! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_mom
  use mpi
  use mod_param, only: is_boussinesq_buoyancy,is_harmonic_interpolation,is_ifc_split,is_visc_limiter, &
                       is_brackbill_normalization
  use mod_types
  implicit none
  private
  public mom_xyz_ad,mom_xyz_a_cons,mom_xyz_oth,mom_xyz_marangoni,cmpt_wallshear,cmpt_bodyforce,cmpt_momentum,cmpt_capillary_force
  !
  type, public :: marker_fields
    real(rp), pointer, contiguous, dimension(:,:,:) :: psi => null(),kappa => null() ! phase indicator and curvature of a
                                                                                     ! marker (see `rk`)
  end type marker_fields
  contains
  !
  subroutine momx_a(nx,ny,nz,nh_v,nh_g,dxi,dyi,dzfi,u,v,w,dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    integer , intent(in) :: nh_v,nh_g
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(1-nh_g:) :: dzfi
    real(rp), dimension(1-nh_v:,1-nh_v:,1-nh_v:), intent(in   ) :: u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    integer :: i,j,k
    real(rp) :: uuip,uuim,vujp,vujm,wukp,wukm
    !
    !$acc parallel     loop collapse(3) default(present) private(uuip,uuim,vujp,vujm,wukp,wukm) async(1)
    !$omp target teams loop collapse(3)                  private(uuip,uuim,vujp,vujm,wukp,wukm)
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
  subroutine momy_a(nx,ny,nz,nh_v,nh_g,dxi,dyi,dzfi,u,v,w,dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    integer , intent(in) :: nh_v,nh_g
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(1-nh_g:) :: dzfi
    real(rp), dimension(1-nh_v:,1-nh_v:,1-nh_v:), intent(in   ) :: u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    integer :: i,j,k
    real(rp) :: uvip,uvim,vvjp,vvjm,wvkp,wvkm
    !
    !$acc parallel     loop collapse(3) default(present) private(uvip,uvim,vvjp,vvjm,wvkp,wvkm) async(1)
    !$omp target teams loop collapse(3)                  private(uvip,uvim,vvjp,vvjm,wvkp,wvkm)
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
  subroutine momz_a(nx,ny,nz,nh_v,nh_g,dxi,dyi,dzci,u,v,w,dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    integer , intent(in) :: nh_v,nh_g
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(1-nh_g:) :: dzci
    real(rp), dimension(1-nh_v:,1-nh_v:,1-nh_v:), intent(in   ) :: u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    integer :: i,j,k
    real(rp) :: uwip,uwim,vwjp,vwjm,wwkp,wwkm
    !
    !$acc parallel     loop collapse(3) default(present) private(uwip,uwim,vwjp,vwjm,wwkp,wwkm) async(1)
    !$omp target teams loop collapse(3)                  private(uwip,uwim,vwjp,vwjm,wwkp,wwkm)
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
  subroutine momx_d(nx,ny,nz,nh_v,nh_psi,nh_g,dxi,dyi,dzci,dzfi,mu12,psi,u,v,w,dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    integer , intent(in) :: nh_v,nh_psi,nh_g
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(1-nh_g:) :: dzci,dzfi
    real(rp), intent(in), dimension(2) :: mu12
    real(rp), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:), intent(in   ) :: psi
    real(rp), dimension(1-nh_v:,1-nh_v:,1-nh_v:), intent(in   ) :: u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm, &
                dvdxp,dvdxm,dwdxp,dwdxm
    real(rp) :: muxp,muxm,muyp,muym,muzp,muzm
    integer :: i,j,k
    real(rp) :: mu,dmu
    logical  :: is_harm
    !
    is_harm = is_harmonic_interpolation
    if(is_harm) then
      mu  = 1._rp/mu12(2); dmu = 1._rp/mu12(1)-1._rp/mu12(2)
    else
      mu  = mu12(2)      ; dmu = mu12(1)-mu12(2)
    end if
    !
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(dudxp,dudxm,dudyp,dudym,dudzp,dudzm,dvdxp,dvdxm,dwdxp,dwdxm) &
    !$acc private(muxp,muxm,muyp,muym,muzp,muzm) async(1)
    !$omp target teams loop collapse(3) &
    !$omp private(dudxp,dudxm,dudyp,dudym,dudzp,dudzm,dvdxp,dvdxm,dwdxp,dwdxm) &
    !$omp private(muxp,muxm,muyp,muym,muzp,muzm)
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
          muxp = mu+dmu*psi(i+1,j,k)
          muxm = mu+dmu*psi(i  ,j,k)
          muyp = mu+dmu*0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i+1,j+1,k)+psi(i+1,j,k))
          muym = mu+dmu*0.25*(psi(i,j,k)+psi(i,j-1,k)+psi(i+1,j-1,k)+psi(i+1,j,k))
          muzp = mu+dmu*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i+1,j,k+1)+psi(i+1,j,k))
          muzm = mu+dmu*0.25*(psi(i,j,k)+psi(i,j,k-1)+psi(i+1,j,k-1)+psi(i+1,j,k))
          if(is_harm) then
            muxp = 1._rp/muxp
            muxm = 1._rp/muxm
            muyp = 1._rp/muyp
            muym = 1._rp/muym
            muzp = 1._rp/muzp
            muzm = 1._rp/muzm
          end if
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
  subroutine momy_d(nx,ny,nz,nh_v,nh_psi,nh_g,dxi,dyi,dzci,dzfi,mu12,psi,u,v,w,dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    integer , intent(in) :: nh_v,nh_psi,nh_g
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(1-nh_g:) :: dzci,dzfi
    real(rp), intent(in), dimension(2) :: mu12
    real(rp), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:), intent(in   ) :: psi
    real(rp), dimension(1-nh_v:,1-nh_v:,1-nh_v:), intent(in   ) :: u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    real(rp) :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm, &
                dudyp,dudym,dwdyp,dwdym
    real(rp) :: muxp,muxm,muyp,muym,muzp,muzm
    integer :: i,j,k
    real(rp) :: mu,dmu
    logical  :: is_harm
    !
    is_harm = is_harmonic_interpolation
    if(is_harm) then
      mu  = 1._rp/mu12(2); dmu = 1._rp/mu12(1)-1._rp/mu12(2)
    else
      mu  = mu12(2)      ; dmu = mu12(1)-mu12(2)
    end if
    !
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm,dudyp,dudym,dwdyp,dwdym) &
    !$acc private(muxp,muxm,muyp,muym,muzp,muzm) async(1)
    !$omp target teams loop collapse(3) &
    !$omp private(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm,dudyp,dudym,dwdyp,dwdym) &
    !$omp private(muxp,muxm,muyp,muym,muzp,muzm)
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
          muxp = mu+dmu*0.25*(psi(i,j,k)+psi(i+1,j,k)+psi(i+1,j+1,k)+psi(i,j+1,k))
          muxm = mu+dmu*0.25*(psi(i,j,k)+psi(i-1,j,k)+psi(i-1,j+1,k)+psi(i,j+1,k))
          muyp = mu+dmu*psi(i,j+1,k)
          muym = mu+dmu*psi(i,j  ,k)
          muzp = mu+dmu*0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i,j+1,k+1)+psi(i,j,k+1))
          muzm = mu+dmu*0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i,j+1,k-1)+psi(i,j,k-1))
          if(is_harm) then
            muxp = 1._rp/muxp
            muxm = 1._rp/muxm
            muyp = 1._rp/muyp
            muym = 1._rp/muym
            muzp = 1._rp/muzp
            muzm = 1._rp/muzm
          end if
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
  subroutine momz_d(nx,ny,nz,nh_v,nh_psi,nh_g,dxi,dyi,dzci,dzfi,mu12,psi,u,v,w,dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    integer , intent(in) :: nh_v,nh_psi,nh_g
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(1-nh_g:) :: dzci,dzfi
    real(rp), intent(in), dimension(2) :: mu12
    real(rp), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:), intent(in   ) :: psi
    real(rp), dimension(1-nh_v:,1-nh_v:,1-nh_v:), intent(in   ) :: u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    integer :: i,j,k
    real(rp) :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm, &
                dudzp,dudzm,dvdzp,dvdzm
    real(rp) :: muxp,muxm,muyp,muym,muzp,muzm
    real(rp) :: mu,dmu
    logical  :: is_harm
    !
    is_harm = is_harmonic_interpolation
    if(is_harm) then
      mu  = 1._rp/mu12(2); dmu = 1._rp/mu12(1)-1._rp/mu12(2)
    else
      mu  = mu12(2)      ; dmu = mu12(1)-mu12(2)
    end if
    !
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm,dudzp,dudzm,dvdzp,dvdzm) &
    !$acc private(muxp,muxm,muyp,muym,muzp,muzm) async(1)
    !$omp target teams loop collapse(3) &
    !$omp private(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm,dudzp,dudzm,dvdzp,dvdzm) &
    !$omp private(muxp,muxm,muyp,muym,muzp,muzm)
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
          muxp = mu+dmu*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i+1,j,k+1)+psi(i+1,j,k) )
          muxm = mu+dmu*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i-1,j,k+1)+psi(i-1,j,k) )
          muyp = mu+dmu*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i,j+1,k+1)+psi(i,j+1,k) )
          muym = mu+dmu*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i,j-1,k+1)+psi(i,j-1,k) )
          muzp = mu+dmu*psi(i,j,k+1)
          muzm = mu+dmu*psi(i,j,k  )
          if(is_harm) then
            muxp = 1._rp/muxp
            muxm = 1._rp/muxm
            muyp = 1._rp/muyp
            muym = 1._rp/muym
            muzp = 1._rp/muzp
            muzm = 1._rp/muzm
          end if
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
  subroutine momx_p(nx,ny,nz,nh_p,nh_psi,dxi,dt_r,bforce,gacc,rho0,rhox_av,rho12,psi,p,pp,dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    integer , intent(in) :: nh_p,nh_psi
    real(rp), intent(in) :: dxi,dt_r
    real(rp), intent(in) :: bforce,gacc,rho0,rhox_av,rho12(2)
    real(rp), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:), intent(in   ) :: psi
    real(rp), dimension(1-nh_p:,1-nh_p:,1-nh_p:), intent(in   ) :: p,pp
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    real(rp) :: rhop,dpdl
    integer :: i,j,k
    real(rp) :: rho,drho
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    !
    !$acc parallel     loop collapse(3) default(present) private(dpdl,rhop) async(1)
    !$omp target teams loop collapse(3)                  private(dpdl,rhop)
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
  subroutine momy_p(nx,ny,nz,nh_p,nh_psi,dyi,dt_r,bforce,gacc,rho0,rhoy_av,rho12,psi,p,pp,dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    integer , intent(in) :: nh_p,nh_psi
    real(rp), intent(in) :: dyi,dt_r
    real(rp), intent(in) :: bforce,gacc,rho0,rhoy_av,rho12(2)
    real(rp), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:), intent(in   ) :: psi
    real(rp), dimension(1-nh_p:,1-nh_p:,1-nh_p:), intent(in   ) :: p,pp
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    integer :: i,j,k
    real(rp) :: rhop,dpdl
    real(rp) :: rho,drho
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    !
    !$acc parallel     loop collapse(3) default(present) private(rhop,dpdl) async(1)
    !$omp target teams loop collapse(3)                  private(rhop,dpdl)
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
  subroutine momz_p(nx,ny,nz,nh_p,nh_psi,nh_g,dzci,dt_r,bforce,gacc,rho0,rhoz_av,rho12,psi,p,pp,dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    integer , intent(in) :: nh_p,nh_psi,nh_g
    real(rp), intent(in), dimension(1-nh_g:) :: dzci
    real(rp), intent(in) :: dt_r
    real(rp), intent(in) :: bforce,gacc,rho0,rhoz_av,rho12(2)
    real(rp), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:), intent(in   ) :: psi
    real(rp), dimension(1-nh_p:,1-nh_p:,1-nh_p:), intent(in   ) :: p,pp
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    real(rp) :: rhop,dpdl
    integer :: i,j,k
    real(rp) :: rho,drho
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    !
    !$acc parallel     loop collapse(3) default(present) private(rhop,dpdl) async(1)
    !$omp target teams loop collapse(3)                  private(rhop,dpdl)
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
  subroutine momx_sigma(nx,ny,nz,nh_psi,nh_ifc,dxi,sigma,kappa,psi,dudt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    integer , intent(in) :: nh_psi,nh_ifc
    real(rp), intent(in) :: dxi
    real(rp), intent(in) :: sigma
    real(rp), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:), intent(in   ) :: psi
    real(rp), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:), intent(in   ) :: kappa
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    real(rp) :: kappap
    integer :: i,j,k
    !
    !$acc parallel     loop collapse(3) default(present) private(kappap) async(1)
    !$omp target teams loop collapse(3)                  private(kappap)
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
  subroutine momy_sigma(nx,ny,nz,nh_psi,nh_ifc,dyi,sigma,psi,kappa,dvdt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    integer , intent(in) :: nh_psi,nh_ifc
    real(rp), intent(in) :: dyi
    real(rp), intent(in) :: sigma
    real(rp), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:), intent(in   ) :: psi
    real(rp), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:), intent(in   ) :: kappa
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    real(rp) :: kappap
    integer :: i,j,k
    !
    !$acc parallel     loop collapse(3) default(present) private(kappap) async(1)
    !$omp target teams loop collapse(3)                  private(kappap)
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
  subroutine momz_sigma(nx,ny,nz,nh_psi,nh_ifc,nh_g,dzci,sigma,psi,kappa,dwdt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    integer , intent(in) :: nh_psi,nh_ifc,nh_g
    real(rp), intent(in), dimension(1-nh_g:) :: dzci
    real(rp), intent(in) :: sigma
    real(rp), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:), intent(in   ) :: psi
    real(rp), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:), intent(in   ) :: kappa
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    real(rp) :: kappap
    integer :: i,j,k
    !
    !$acc parallel     loop collapse(3) default(present) private(kappap) async(1)
    !$omp target teams loop collapse(3)                  private(kappap)
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
  subroutine momx_buoy(nx,ny,nz,nh_s,nh_psi,gacc,rho12,beta12,psi,s,dudt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    integer , intent(in) :: nh_s,nh_psi
    real(rp), intent(in) :: gacc,rho12(2),beta12(2)
    real(rp), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:), intent(in   ) :: psi
    real(rp), dimension(1-nh_s:,1-nh_s:,1-nh_s:), intent(in   ) :: s
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    real(rp) :: psip,factorp
    integer :: i,j,k
    real(rp) :: rhobeta,drhobeta
    !
    rhobeta = rho12(2)*beta12(2); drhobeta = rho12(1)*beta12(1)-rho12(2)*beta12(2)
    !
    !$acc parallel     loop collapse(3) default(present) private(psip,factorp) async(1)
    !$omp target teams loop collapse(3)                  private(psip,factorp)
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
  subroutine momy_buoy(nx,ny,nz,nh_s,nh_psi,gacc,rho12,beta12,psi,s,dvdt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    integer , intent(in) :: nh_s,nh_psi
    real(rp), intent(in) :: gacc,rho12(2),beta12(2)
    real(rp), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:), intent(in   ) :: psi
    real(rp), dimension(1-nh_s:,1-nh_s:,1-nh_s:), intent(in   ) :: s
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    real(rp) :: psip,factorp
    integer :: i,j,k
    real(rp) :: rhobeta,drhobeta
    !
    rhobeta = rho12(2)*beta12(2); drhobeta = rho12(1)*beta12(1)-rho12(2)*beta12(2)
    !
    !$acc parallel     loop collapse(3) default(present) private(factorp) async(1)
    !$omp target teams loop collapse(3)                  private(factorp)
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
  subroutine momz_buoy(nx,ny,nz,nh_s,nh_psi,gacc,rho12,beta12,psi,s,dwdt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    integer , intent(in) :: nh_s,nh_psi
    real(rp), intent(in) :: gacc,rho12(2),beta12(2)
    real(rp), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:), intent(in   ) :: psi
    real(rp), dimension(1-nh_s:,1-nh_s:,1-nh_s:), intent(in   ) :: s
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    real(rp) :: psip,factorp
    integer :: i,j,k
    real(rp) :: rhobeta,drhobeta
    !
    rhobeta = rho12(2)*beta12(2); drhobeta = rho12(1)*beta12(1)-rho12(2)*beta12(2)
    !
    !$acc parallel     loop collapse(3) default(present) private(factorp) async(1)
    !$omp target teams loop collapse(3)                  private(factorp)
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
  subroutine mom_xyz_ad(n,nh_v,nh_psi,nh_ifc,nh_g,dli,dzci,dzfi,wzc,rho12,mu12,psiflx_x,psiflx_y,psiflx_z,u,v,w,psi,dudt,dvdt,dwdt)
    implicit none
    integer , intent(in ), dimension(3) :: n
    integer , intent(in) :: nh_v,nh_psi,nh_ifc,nh_g
    real(rp), intent(in ), dimension(3) :: dli
    real(rp), intent(in ), dimension(1-nh_g:) :: dzci,dzfi,wzc
    real(rp), intent(in ), dimension(2) :: rho12,mu12
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(in ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in ), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:), optional :: psiflx_x,psiflx_y,psiflx_z
    real(rp), intent(out), dimension( :, :, :) :: dudt,dvdt,dwdt
    integer :: i,j,k
    real(rp) :: rho,drho,mu,dmu,rhobeta,drhobeta,nu1,nu2
    real(rp) :: dxi,dyi
    real(rp) :: c_ccm,c_pcm,c_cpm,c_cmc,c_pmc,c_mcc,c_ccc,c_pcc,c_mpc,c_cpc,c_cmp,c_mcp,c_ccp,c_cpp,c_ppc,c_pcp, &
                d_ccm,d_pcm,d_cpm,d_cmc,d_pmc,d_mcc,d_ccc,d_pcc,d_mpc,d_cpc,d_cmp,d_mcp,d_ccp,d_cpp,d_ppc,d_pcp, &
                u_ccm,u_pcm,u_cpm,u_cmc,u_pmc,u_mcc,u_ccc,u_pcc,u_mpc,u_cpc,u_cmp,u_mcp,u_ccp, &
                v_ccm,v_pcm,v_cpm,v_cmc,v_pmc,v_mcc,v_ccc,v_pcc,v_mpc,v_cpc,v_cmp,v_mcp,v_ccp, &
                w_ccm,w_pcm,w_cpm,w_cmc,w_pmc,w_mcc,w_ccc,w_pcc,w_mpc,w_cpc,w_cmp,w_mcp,w_ccp, &
                flx_x_mcc,flx_x_ccc,flx_x_pcc,flx_x_mpc,flx_x_cpc,flx_x_mcp,flx_x_ccp, &
                flx_y_cmc,flx_y_ccc,flx_y_cpc,flx_y_pmc,flx_y_pcc,flx_y_cmp,flx_y_ccp, &
                flx_z_ccm,flx_z_ccc,flx_z_ccp,flx_z_pcm,flx_z_pcc,flx_z_cpm,flx_z_cpc, &
                dzci_c,dzci_m,dzfi_c,dzfi_p,wgt_c,wgt_p,wgt_m
    real(rp) :: uuip,uuim,vujp,vujm,wukp,wukm,uvip,uvim,vvjp,vvjm,wvkp,wvkm,uwip,uwim,vwjp,vwjm,wwkp,wwkm
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm,dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm,dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    real(rp) :: muxp,muxm,muyp,muym,muzp,muzm
    real(rp) :: rhoxp,rhoyp,rhozp
    real(rp) :: dudt_aux,dvdt_aux,dwdt_aux
    logical  :: is_adv,is_harm,is_lim,is_lim_c
    real(rp) :: c_min,c_max
    real(rp), parameter :: c_eps = max(1.e-12_rp,1.e3_rp*epsilon(1._rp)) ! full and empty cells for the viscosity
                                                                         ! limiter (as for the reconstruction)
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    is_harm = is_harmonic_interpolation
    if(is_harm) then
      mu  = 1._rp/mu12(2); dmu = 1._rp/mu12(1)-1._rp/mu12(2)
    else
      mu  = mu12(2)      ; dmu = mu12(1)-mu12(2)
    end if
    nu1 = mu12(1)/rho12(1)
    nu2 = mu12(2)/rho12(2)
    is_adv = .not.is_ifc_split ! (otherwise, advection is computed in `mom_xyz_a_cons` or along with the phase indicator)
    is_lim = is_visc_limiter
    dxi = dli(1)
    dyi = dli(2)
    !
    !
    ! making an exception for this kernel -- private variables not explicitly mentioned for the sake of conciseness
    !                                        all scalars should be firstprivate/private
    !
    !$acc parallel     loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3)
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
          ! the viscosity limiter only matters where the two phases meet: it is skipped where all the cells of
          ! the stencil are full, or all empty, to the tolerance `c_eps`
          !
          is_lim_c = is_lim
          if(is_lim) then
            c_min = min(c_ccm,c_pcm,c_cpm,c_cmc,c_pmc,c_mcc,c_ccc,c_pcc,c_mpc,c_cpc,c_cmp,c_mcp,c_ccp,c_cpp,c_ppc,c_pcp)
            c_max = max(c_ccm,c_pcm,c_cpm,c_cmc,c_pmc,c_mcc,c_ccc,c_pcc,c_mpc,c_cpc,c_cmp,c_mcp,c_ccp,c_cpp,c_ppc,c_pcp)
            is_lim_c = c_min < 1._rp-c_eps .and. c_max > c_eps
          end if
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
          wgt_c  = wzc(k)      ! fractions of the control volume of w within the cells k and k+1
          wgt_p  = 1.-wzc(k)
          wgt_m  = wzc(k-1)
          !
          ! advection
          !
          if(is_adv) then
            !
            ! mass--momentum-consistent flx (the mass fluxes of the staggered control volumes are averages of
            ! those of the cell faces weighted by the overlapped volumes, as their densities are)
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
            uwip = 0.5*(wgt_c*flx_x_ccc+wgt_p*flx_x_ccp)*(w_ccc+w_pcc)
            uwim = 0.5*(wgt_c*flx_x_mcc+wgt_p*flx_x_mcp)*(w_ccc+w_mcc)
            vwjp = 0.5*(wgt_c*flx_y_ccc+wgt_p*flx_y_ccp)*(w_ccc+w_cpc)
            vwjm = 0.5*(wgt_c*flx_y_cmc+wgt_p*flx_y_cmp)*(w_ccc+w_cmc)
            wwkp = 0.25*(flx_z_ccc+flx_z_ccp)*(w_ccc+w_ccp)
            wwkm = 0.25*(flx_z_ccc+flx_z_ccm)*(w_ccc+w_ccm)
            dwdt_aux = dxi*( -uwip + uwim ) + dyi*( -vwjp + vwjm ) + dzci_c*( -wwkp + wwkm )
          else
            !
            ! (computed elsewhere)
            !
            dudt_aux = 0.
            dvdt_aux = 0.
            dwdt_aux = 0.
          end if
          !
          ! diffusion
          !
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
          muxp = mu+dmu*c_pcc
          muxm = mu+dmu*c_ccc
          muyp = mu+dmu*0.25*(c_ccc+c_cpc+c_ppc+c_pcc)
          muym = mu+dmu*0.25*(c_ccc+c_cmc+c_pmc+c_pcc)
          muzp = mu+dmu*0.25*(c_ccc+c_pcc+c_ccp+c_pcp)
          muzm = mu+dmu*0.25*(c_ccc+c_pcc+c_ccm+c_pcm)
          if(is_harm) then
            muxp = 1._rp/muxp
            muxm = 1._rp/muxm
            muyp = 1._rp/muyp
            muym = 1._rp/muym
            muzp = 1._rp/muzp
            muzm = 1._rp/muzm
          end if
          if(is_lim_c) then ! kinematic viscosity bounded at the interface, see `mu_lim_edge`
            muxp = mu_lim_cell(muxp,c_pcc,rho,drho,nu1,nu2)
            muxm = mu_lim_cell(muxm,c_ccc,rho,drho,nu1,nu2)
            muyp = mu_lim_edge(muyp,c_ccc,c_pcc,c_cpc,c_ppc,0.5_rp,rho,drho,nu1,nu2)
            muym = mu_lim_edge(muym,c_cmc,c_pmc,c_ccc,c_pcc,0.5_rp,rho,drho,nu1,nu2)
            muzp = mu_lim_edge(muzp,c_ccc,c_pcc,c_ccp,c_pcp,wgt_c ,rho,drho,nu1,nu2)
            muzm = mu_lim_edge(muzm,c_ccm,c_pcm,c_ccc,c_pcc,wgt_m ,rho,drho,nu1,nu2)
          end if
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
          muxp = mu+dmu*0.25*(c_ccc+c_pcc+c_ppc+c_cpc)
          muxm = mu+dmu*0.25*(c_ccc+c_mcc+c_mpc+c_cpc)
          muyp = mu+dmu*c_cpc
          muym = mu+dmu*c_ccc
          muzp = mu+dmu*0.25*(c_ccc+c_cpc+c_cpp+c_ccp)
          muzm = mu+dmu*0.25*(c_ccc+c_cpc+c_cpm+c_ccm)
          if(is_harm) then
            muxp = 1._rp/muxp
            muxm = 1._rp/muxm
            muyp = 1._rp/muyp
            muym = 1._rp/muym
            muzp = 1._rp/muzp
            muzm = 1._rp/muzm
          end if
          if(is_lim_c) then
            muxp = mu_lim_edge(muxp,c_ccc,c_pcc,c_cpc,c_ppc,0.5_rp,rho,drho,nu1,nu2)
            muxm = mu_lim_edge(muxm,c_mcc,c_ccc,c_mpc,c_cpc,0.5_rp,rho,drho,nu1,nu2)
            muyp = mu_lim_cell(muyp,c_cpc,rho,drho,nu1,nu2)
            muym = mu_lim_cell(muym,c_ccc,rho,drho,nu1,nu2)
            muzp = mu_lim_edge(muzp,c_ccc,c_cpc,c_ccp,c_cpp,wgt_c ,rho,drho,nu1,nu2)
            muzm = mu_lim_edge(muzm,c_ccm,c_cpm,c_ccc,c_cpc,wgt_m ,rho,drho,nu1,nu2)
          end if
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
          muxp = mu+dmu*0.25*(c_ccc+c_pcc+c_ccp+c_pcp)
          muxm = mu+dmu*0.25*(c_ccc+c_mcc+c_ccp+c_mcp)
          muyp = mu+dmu*0.25*(c_ccc+c_cpc+c_ccp+c_cpp)
          muym = mu+dmu*0.25*(c_ccc+c_cmc+c_ccp+c_cmp)
          muzp = mu+dmu*c_ccp
          muzm = mu+dmu*c_ccc
          if(is_harm) then
            muxp = 1._rp/muxp
            muxm = 1._rp/muxm
            muyp = 1._rp/muyp
            muym = 1._rp/muym
            muzp = 1._rp/muzp
            muzm = 1._rp/muzm
          end if
          if(is_lim_c) then
            muxp = mu_lim_edge(muxp,c_ccc,c_pcc,c_ccp,c_pcp,wgt_c ,rho,drho,nu1,nu2)
            muxm = mu_lim_edge(muxm,c_mcc,c_ccc,c_mcp,c_ccp,wgt_c ,rho,drho,nu1,nu2)
            muyp = mu_lim_edge(muyp,c_ccc,c_cpc,c_ccp,c_cpp,wgt_c ,rho,drho,nu1,nu2)
            muym = mu_lim_edge(muym,c_cmc,c_ccc,c_cmp,c_ccp,wgt_c ,rho,drho,nu1,nu2)
            muzp = mu_lim_cell(muzp,c_ccp,rho,drho,nu1,nu2)
            muzm = mu_lim_cell(muzm,c_ccc,rho,drho,nu1,nu2)
          end if
          dwdt_aux = dwdt_aux + dxi*(   (dwdxp+dudzp)*muxp-(dwdxm+dudzm)*muxm) + &
                                dyi*(   (dwdyp+dvdzp)*muyp-(dwdym+dvdzm)*muym) + &
                                dzci_c*((dwdzp+dwdzp)*muzp-(dwdzm+dwdzm)*muzm)
          !
          dudt(i,j,k) = dudt_aux
          dvdt(i,j,k) = dvdt_aux
          dwdt(i,j,k) = dwdt_aux
        end do
      end do
    end do
  end subroutine mom_xyz_ad
  !
  pure real(rp) function mu_lim_edge(mu_e,ca,cb,cc,cd,wz,rho,drho,nu1,nu2) result(mu_l)
    !
    ! limits the viscosity on a cell edge such that its kinematic viscosity, relative to the lightest of the
    ! four staggered control volumes that share the edge, does not exceed that of the phase that dominates
    ! the four cells around the edge (T.-Y. Huang et al., arXiv:2606.02467, 2026): with a sharp interface,
    ! the viscosity interpolated to an edge can otherwise be that of the heavy phase next to a control volume
    ! with the density of the light one, with a kinematic viscosity far above those of both phases;
    ! ca, cb, cc, cd are the volume fractions of the cells (0,0), (1,0), (0,1), (1,1) around the edge, and
    ! wz the weight of the cells (.,0) in the control volumes along the second direction (see `wzc`)
    !
    implicit none
    real(rp), intent(in) :: mu_e,ca,cb,cc,cd,wz,rho,drho,nu1,nu2
    real(rp) :: rho_min
    !$acc routine seq
    !$omp declare target
    rho_min = min(rho+drho*0.5_rp*(ca+cb),rho+drho*0.5_rp*(cc+cd), &
                  rho+drho*(wz*ca+(1._rp-wz)*cc),rho+drho*(wz*cb+(1._rp-wz)*cd))
    mu_l = min(mu_e,rho_min*merge(nu1,nu2,ca+cb+cc+cd >= 2._rp))
  end function mu_lim_edge
  !
  pure real(rp) function mu_lim_cell(mu_c,c,rho,drho,nu1,nu2) result(mu_l)
    !
    ! same as `mu_lim_edge`, for the viscosity at a cell centre (relative to the density of the cell)
    !
    implicit none
    real(rp), intent(in) :: mu_c,c,rho,drho,nu1,nu2
    !$acc routine seq
    !$omp declare target
    mu_l = min(mu_c,(rho+drho*c)*merge(nu1,nu2,c >= 0.5_rp))
  end function mu_lim_cell
  !
  subroutine mom_xyz_a_cons(n,nh_v,nh_psi,nh_ifc,nh_g,dli,dzci,dzfi,wzc,rho12,psi,psio,psiflx_x,psiflx_y,psiflx_z,ua,va,wa, &
                            factor1,factor2,factor12,dudt,dvdt,dwdt,dudto,dvdto,dwdto,u,v,w)
    !
    ! advection terms in conservative form, consistent with the transport of the geometric volume-of-fluid
    ! method, integrated in one step (over the time `factor12`), together with the AB2 update of the velocity
    ! with the other terms (`dudt` and `dudto`, with weights `factor1` and `factor2`), in a single kernel:
    !   u <- (rho_old*u + factor12*advection + factor1*dudt + factor2*dudto)/rho_new,   dudto <- dudt
    !
    ! the cell-face mass fluxes are rho_2*ua + (rho_1-rho_2)*psiflx, with psiflx the (time-averaged) phase-1
    ! volume fluxes of the volume-of-fluid step and ua its divergence-free advecting velocity, which is also
    ! the transported one; the mass fluxes of the staggered control volumes are averages of those of the cell
    ! faces weighted by the overlapped volumes (`wzc` along z), as the staggered densities are averages of
    ! those of the cells: the discrete continuity equation of the staggered control volumes is then satisfied
    ! exactly, and the central interpolation of the transported velocity makes the advection operator
    ! skew-symmetric (kinetic energy conserving)
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   ) :: nh_v,nh_psi,nh_ifc,nh_g
    real(rp), intent(in   ), dimension(3) :: dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzci,dzfi,wzc
    real(rp), intent(in   ), dimension(2) :: rho12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi,psio
    real(rp), intent(in   ), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: psiflx_x,psiflx_y,psiflx_z
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: ua,va,wa
    real(rp), intent(in   ) :: factor1,factor2,factor12
    real(rp), intent(in   ), dimension(:,:,:) :: dudt,dvdt,dwdt
    real(rp), intent(inout), dimension(:,:,:) :: dudto,dvdto,dwdto
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    integer :: i,j,k
    real(rp) :: rho,drho,dxi,dyi,dzci_c,dzfi_c,wgt_c,wgt_p
    real(rp) :: u_ccm,u_cmc,u_mcc,u_ccc,u_pcc,u_cpc,u_ccp, &
                v_ccm,v_cmc,v_mcc,v_ccc,v_pcc,v_cpc,v_ccp, &
                w_ccm,w_cmc,w_mcc,w_ccc,w_pcc,w_cpc,w_ccp, &
                flx_x_mcc,flx_x_ccc,flx_x_pcc,flx_x_mpc,flx_x_cpc,flx_x_mcp,flx_x_ccp, &
                flx_y_cmc,flx_y_ccc,flx_y_cpc,flx_y_pmc,flx_y_pcc,flx_y_cmp,flx_y_ccp, &
                flx_z_ccm,flx_z_ccc,flx_z_ccp,flx_z_pcm,flx_z_pcc,flx_z_cpm,flx_z_cpc
    real(rp) :: uuip,uuim,vujp,vujm,wukp,wukm,uvip,uvim,vvjp,vvjm,wvkp,wvkm,uwip,uwim,vwjp,vwjm,wwkp,wwkm
    real(rp) :: rhoxp,rhoyp,rhozp,rhoxo,rhoyo,rhozo,dudt_a,dvdt_a,dwdt_a
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    dxi = dli(1)
    dyi = dli(2)
    !
    ! making an exception for this kernel -- private variables not explicitly mentioned for the sake of conciseness
    !                                        all scalars should be firstprivate/private
    !
    !$acc parallel     loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          u_ccm = ua(i  ,j  ,k-1)
          u_cmc = ua(i  ,j-1,k  )
          u_mcc = ua(i-1,j  ,k  )
          u_ccc = ua(i  ,j  ,k  )
          u_pcc = ua(i+1,j  ,k  )
          u_cpc = ua(i  ,j+1,k  )
          u_ccp = ua(i  ,j  ,k+1)
          !
          v_ccm = va(i  ,j  ,k-1)
          v_cmc = va(i  ,j-1,k  )
          v_mcc = va(i-1,j  ,k  )
          v_ccc = va(i  ,j  ,k  )
          v_pcc = va(i+1,j  ,k  )
          v_cpc = va(i  ,j+1,k  )
          v_ccp = va(i  ,j  ,k+1)
          !
          w_ccm = wa(i  ,j  ,k-1)
          w_cmc = wa(i  ,j-1,k  )
          w_mcc = wa(i-1,j  ,k  )
          w_ccc = wa(i  ,j  ,k  )
          w_pcc = wa(i+1,j  ,k  )
          w_cpc = wa(i  ,j+1,k  )
          w_ccp = wa(i  ,j  ,k+1)
          !
          flx_x_mcc = rho*u_mcc           + drho*psiflx_x(i-1,j  ,k  )
          flx_x_ccc = rho*u_ccc           + drho*psiflx_x(i  ,j  ,k  )
          flx_x_pcc = rho*u_pcc           + drho*psiflx_x(i+1,j  ,k  )
          flx_x_mpc = rho*ua(i-1,j+1,k  ) + drho*psiflx_x(i-1,j+1,k  )
          flx_x_cpc = rho*u_cpc           + drho*psiflx_x(i  ,j+1,k  )
          flx_x_mcp = rho*ua(i-1,j  ,k+1) + drho*psiflx_x(i-1,j  ,k+1)
          flx_x_ccp = rho*u_ccp           + drho*psiflx_x(i  ,j  ,k+1)
          !
          flx_y_cmc = rho*v_cmc           + drho*psiflx_y(i  ,j-1,k  )
          flx_y_ccc = rho*v_ccc           + drho*psiflx_y(i  ,j  ,k  )
          flx_y_cpc = rho*v_cpc           + drho*psiflx_y(i  ,j+1,k  )
          flx_y_pmc = rho*va(i+1,j-1,k  ) + drho*psiflx_y(i+1,j-1,k  )
          flx_y_pcc = rho*v_pcc           + drho*psiflx_y(i+1,j  ,k  )
          flx_y_cmp = rho*va(i  ,j-1,k+1) + drho*psiflx_y(i  ,j-1,k+1)
          flx_y_ccp = rho*v_ccp           + drho*psiflx_y(i  ,j  ,k+1)
          !
          flx_z_ccm = rho*w_ccm           + drho*psiflx_z(i  ,j  ,k-1)
          flx_z_ccc = rho*w_ccc           + drho*psiflx_z(i  ,j  ,k  )
          flx_z_ccp = rho*w_ccp           + drho*psiflx_z(i  ,j  ,k+1)
          flx_z_pcm = rho*wa(i+1,j  ,k-1) + drho*psiflx_z(i+1,j  ,k-1)
          flx_z_pcc = rho*w_pcc           + drho*psiflx_z(i+1,j  ,k  )
          flx_z_cpm = rho*wa(i  ,j+1,k-1) + drho*psiflx_z(i  ,j+1,k-1)
          flx_z_cpc = rho*w_cpc           + drho*psiflx_z(i  ,j+1,k  )
          !
          dzci_c = dzci(k)
          dzfi_c = dzfi(k)
          wgt_c  = wzc(k)      ! fractions of the control volume of w within the cells k and k+1
          wgt_p  = 1.-wzc(k)
          !
          uuip = 0.25*(flx_x_pcc+flx_x_ccc)*(u_ccc+u_pcc)
          uuim = 0.25*(flx_x_mcc+flx_x_ccc)*(u_ccc+u_mcc)
          vujp = 0.25*(flx_y_pcc+flx_y_ccc)*(u_ccc+u_cpc)
          vujm = 0.25*(flx_y_pmc+flx_y_cmc)*(u_ccc+u_cmc)
          wukp = 0.25*(flx_z_pcc+flx_z_ccc)*(u_ccc+u_ccp)
          wukm = 0.25*(flx_z_pcm+flx_z_ccm)*(u_ccc+u_ccm)
          !
          uvip = 0.25*(flx_x_ccc+flx_x_cpc)*(v_ccc+v_pcc)
          uvim = 0.25*(flx_x_mcc+flx_x_mpc)*(v_ccc+v_mcc)
          vvjp = 0.25*(flx_y_ccc+flx_y_cpc)*(v_ccc+v_cpc)
          vvjm = 0.25*(flx_y_ccc+flx_y_cmc)*(v_ccc+v_cmc)
          wvkp = 0.25*(flx_z_ccc+flx_z_cpc)*(v_ccc+v_ccp)
          wvkm = 0.25*(flx_z_ccm+flx_z_cpm)*(v_ccc+v_ccm)
          !
          uwip = 0.5*(wgt_c*flx_x_ccc+wgt_p*flx_x_ccp)*(w_ccc+w_pcc)
          uwim = 0.5*(wgt_c*flx_x_mcc+wgt_p*flx_x_mcp)*(w_ccc+w_mcc)
          vwjp = 0.5*(wgt_c*flx_y_ccc+wgt_p*flx_y_ccp)*(w_ccc+w_cpc)
          vwjm = 0.5*(wgt_c*flx_y_cmc+wgt_p*flx_y_cmp)*(w_ccc+w_cmc)
          wwkp = 0.25*(flx_z_ccc+flx_z_ccp)*(w_ccc+w_ccp)
          wwkm = 0.25*(flx_z_ccc+flx_z_ccm)*(w_ccc+w_ccm)
          !
          dudt_a = dxi*( -uuip + uuim ) + dyi*( -vujp + vujm ) + dzfi_c*( -wukp + wukm )
          dvdt_a = dxi*( -uvip + uvim ) + dyi*( -vvjp + vvjm ) + dzfi_c*( -wvkp + wvkm )
          dwdt_a = dxi*( -uwip + uwim ) + dyi*( -vwjp + vwjm ) + dzci_c*( -wwkp + wwkm )
          !
          ! update, with the densities of the staggered control volumes before and after the volume-of-fluid step
          !
          rhoxo = rho + drho*0.5*(psio(i+1,j,k)+psio(i,j,k))
          rhoyo = rho + drho*0.5*(psio(i,j+1,k)+psio(i,j,k))
          rhozo = rho + drho*(wgt_c*psio(i,j,k)+wgt_p*psio(i,j,k+1))
          rhoxp = rho + drho*0.5*(psi( i+1,j,k)+psi( i,j,k))
          rhoyp = rho + drho*0.5*(psi( i,j+1,k)+psi( i,j,k))
          rhozp = rho + drho*(wgt_c*psi( i,j,k)+wgt_p*psi( i,j,k+1))
          u(i,j,k) = (u(i,j,k)*rhoxo + factor12*dudt_a + factor1*dudt(i,j,k) + factor2*dudto(i,j,k))/rhoxp
          v(i,j,k) = (v(i,j,k)*rhoyo + factor12*dvdt_a + factor1*dvdt(i,j,k) + factor2*dvdto(i,j,k))/rhoyp
          w(i,j,k) = (w(i,j,k)*rhozo + factor12*dwdt_a + factor1*dwdt(i,j,k) + factor2*dwdto(i,j,k))/rhozp
          dudto(i,j,k) = dudt(i,j,k)
          dvdto(i,j,k) = dvdt(i,j,k)
          dwdto(i,j,k) = dwdt(i,j,k)
        end do
      end do
    end do
  end subroutine mom_xyz_a_cons
  !
  subroutine mom_xyz_oth(n,nh_v,nh_p,nh_s,nh_psi,nh_ifc,nh_g,dli,dzci,dzfi,wzc,dt_r,rho12,beta12,bforce,gacc,sigma,rho0,rho_av, &
                         p,psi,kappa,s,pn,po, &
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
                         surfx_n,surfy_n,surfz_n,surfx_o,surfy_o,surfz_o, &
#endif
                         factor,u,v,w,fsx,fsy,fsz)
    !
    ! pressure, surface tension and body forces: integrated over the time `factor` and added to the velocity,
    ! divided by the face density; if present, fsx, fsy and fsz are the surface tension force at the faces,
    ! computed beforehand (several markers, see `cmpt_capillary_force`), which replaces the one computed here
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in) :: nh_v,nh_p,nh_s,nh_psi,nh_ifc,nh_g
    real(rp), intent(in   ), dimension(3) :: dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzci,dzfi,wzc
    real(rp), intent(in   ) :: dt_r
    real(rp), intent(in   ), dimension(2) :: rho12,beta12
    real(rp), intent(in   ), dimension(3) :: bforce,gacc
    real(rp), intent(in   ) :: sigma
    real(rp), intent(in   ) :: rho0,rho_av(3)
    real(rp), intent(in   ), dimension(1-nh_p:,1-nh_p:,1-nh_p:)           :: p
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:)           :: psi
    real(rp), intent(in   ), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:)           :: kappa
    real(rp), intent(in   ), dimension(1-nh_s:,1-nh_s:,1-nh_s:), optional :: s
    real(rp), intent(in   ), dimension(1-nh_p:,1-nh_p:,1-nh_p:), optional :: pn,po
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
    real(rp), intent(inout), dimension(:,:,:) :: surfx_n,surfy_n,surfz_n, &
                                                 surfx_o,surfy_o,surfz_o
#endif
    real(rp), intent(in   ) :: factor
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(in   ), dimension(:,:,:), optional :: fsx,fsy,fsz
    integer :: i,j,k
    real(rp) :: rho,drho,rhobeta,drhobeta,rho0i
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
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
    real(rp) :: surfx_e,surfy_e,surfz_e
#endif
    real(rp) :: dudt_aux,dvdt_aux,dwdt_aux
    real(rp) :: surf_factor
    logical  :: is_buoyant,is_brack,is_fs
    !
    is_buoyant = is_boussinesq_buoyancy.and.present(s)
    is_fs = present(fsx)
    is_brack = is_brackbill_normalization
    rho     = rho12(2); drho = rho12(1)-rho12(2)
    rho0i   = 1._rp/rho0
    rhobeta = rho12(2)*beta12(2); drhobeta = rho12(1)*beta12(1)- rho12(2)*beta12(2)
    dxi = dli(1)
    dyi = dli(2)
    surf_factor = sigma
    if(is_brack) surf_factor = sigma*2._rp/(rho12(1)+rho12(2)) ! density-based scaling of Brackbill et al.
    rhox_av = rho_av(1); rhoy_av = rho_av(2); rhoz_av = rho_av(3)
    !
    ! making an exception for this kernel -- private variables not explicitly mentioned for the sake of conciseness
    !                                        all scalars should be firstprivate/private
    !
    !$acc parallel     loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3)
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
          if(is_buoyant) then
            s_ccc = s(i  ,j  ,k  )
            s_pcc = s(i+1,j  ,k  )
            s_cpc = s(i  ,j+1,k  )
            s_ccp = s(i  ,j  ,k+1)
          end if
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
          psizp = wzc(k)*c_ccc+(1.-wzc(k))*c_ccp ! volume average
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
          if(is_brack) then
            skappaxp = skappaxp*rhoxp
            skappayp = skappayp*rhoyp
            skappazp = skappazp*rhozp
          end if
          surfx = skappaxp*(c_pcc-c_ccc)*dxi
          surfy = skappayp*(c_cpc-c_ccc)*dyi
          surfz = skappazp*(c_ccp-c_ccc)*dzci_c
          if(is_fs) then
            surfx = fsx(i,j,k)
            surfy = fsy(i,j,k)
            surfz = fsz(i,j,k)
          end if
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
          !
          ! split and extrapolate the balanced residual grad(p)-f_sigma, instead of grad(p) alone, which jumps
          ! across the interface (see P. Cifani, J. Comput. Phys. 398 (2019) 108904)
          !
          surfx_e = (1._rp+dt_r)*surfx_n(i,j,k)-dt_r*surfx_o(i,j,k)
          surfy_e = (1._rp+dt_r)*surfy_n(i,j,k)-dt_r*surfy_o(i,j,k)
          surfz_e = (1._rp+dt_r)*surfz_n(i,j,k)-dt_r*surfz_o(i,j,k)
          dpdx = dpdx-surfx
          dpdy = dpdy-surfy
          dpdz = dpdz-surfz
#else
          dudt_aux = dudt_aux + surfx
          dvdt_aux = dvdt_aux + surfy
          dwdt_aux = dwdt_aux + surfz
#endif
#if defined(_CONSTANT_COEFFS_POISSON)
          dpdx_e = (q_pcc-q_ccc)*dxi
          dpdy_e = (q_cpc-q_ccc)*dyi
          dpdz_e = (q_ccp-q_ccc)*dzci_c
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
          dpdx_e = dpdx_e-surfx_e
          dpdy_e = dpdy_e-surfy_e
          dpdz_e = dpdz_e-surfz_e
#endif
          dudt_aux = dudt_aux - (dpdx_e + rhoxp*rho0i*(dpdx-dpdx_e))
          dvdt_aux = dvdt_aux - (dpdy_e + rhoyp*rho0i*(dpdy-dpdy_e))
          dwdt_aux = dwdt_aux - (dpdz_e + rhozp*rho0i*(dpdz-dpdz_e))
#else
          dudt_aux = dudt_aux - dpdx
          dvdt_aux = dvdt_aux - dpdy
          dwdt_aux = dwdt_aux - dpdz
#endif
          !
          ! buoyancy
          !
          if(is_buoyant) then
            factorxp = rhobeta + drhobeta*psixp
            factoryp = rhobeta + drhobeta*psiyp
            factorzp = rhobeta + drhobeta*psizp
            dudt_aux = dudt_aux - gaccx*factorxp*0.5*(s_pcc+s_ccc)
            dvdt_aux = dvdt_aux - gaccy*factoryp*0.5*(s_cpc+s_ccc)
            dwdt_aux = dwdt_aux - gaccz*factorzp*0.5*(s_ccp+s_ccc)
          end if
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
          surfx_o(i,j,k) = surfx_n(i,j,k)
          surfy_o(i,j,k) = surfy_n(i,j,k)
          surfz_o(i,j,k) = surfz_n(i,j,k)
          surfx_n(i,j,k) = surfx
          surfy_n(i,j,k) = surfy
          surfz_n(i,j,k) = surfz
#endif
          u(i,j,k) = u(i,j,k) + factor*dudt_aux/rhoxp
          v(i,j,k) = v(i,j,k) + factor*dvdt_aux/rhoyp
          w(i,j,k) = w(i,j,k) + factor*dwdt_aux/rhozp
        end do
      end do
    end do
  end subroutine mom_xyz_oth
  !
  subroutine cmpt_bodyforce(n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,is_cmpt,dli,dzci,wzc,grid_vol_ratio_c,grid_vol_ratio_f,dt_r, &
                            rho12,beta12,gacc,sigma,rho0,rho_av,psi,psio,kappa,s,u,v,w, &
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
                            surfx_n,surfy_n,surfz_n,surfx_o,surfy_o,surfz_o, &
#endif
                            fgrav,fsurf,fbuoy,momentum,fsx,fsy,fsz)
    !
    ! computes signed volume averages of the non-pressure body-force
    ! densities exactly as they are applied in mom_xyz_oth.
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    integer , intent(in) :: nh_v,nh_s,nh_psi,nh_ifc,nh_g
    logical , intent(in), dimension(3) :: is_cmpt
    real(rp), intent(in), dimension(3) :: dli
    real(rp), intent(in), dimension(1-nh_g:) :: dzci,wzc,grid_vol_ratio_c,grid_vol_ratio_f
    real(rp), intent(in) :: dt_r
    real(rp), intent(in), dimension(2) :: rho12,beta12
    real(rp), intent(in), dimension(3) :: gacc,rho_av
    real(rp), intent(in) :: sigma,rho0
    real(rp), intent(in), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi,psio
    real(rp), intent(in), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: kappa
    real(rp), intent(in), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(in), dimension(1-nh_s:,1-nh_s:,1-nh_s:), optional :: s
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
    real(rp), intent(in), dimension(:,:,:) :: surfx_n,surfy_n,surfz_n, &
                                              surfx_o,surfy_o,surfz_o
#endif
    real(rp), intent(out), dimension(3) :: fgrav,fsurf,fbuoy,momentum
    real(rp), intent(in ), dimension(:,:,:), optional :: fsx,fsy,fsz ! (see `mom_xyz_oth`)
    real(rp) :: rho,drho,rhobeta,drhobeta,surf_factor,rho0i
    real(rp) :: dxi,dyi,dzci_c,weight_f,weight_c
    real(rp) :: c_ccc,c_pcc,c_cpc,c_ccp,k_ccc,k_pcc,k_cpc,k_ccp
    real(rp) :: co_ccc,co_pcc,co_cpc,co_ccp
    real(rp) :: s_ccc,s_pcc,s_cpc,s_ccp
    real(rp) :: psixp,psiyp,psizp,rhoxp,rhoyp,rhozp
    real(rp) :: rhox_old,rhoy_old,rhoz_old
    real(rp) :: skappaxp,skappayp,skappazp
    real(rp) :: surfx,surfy,surfz,surfx_f,surfy_f,surfz_f
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
    real(rp) :: surfx_e,surfy_e,surfz_e
#endif
    real(rp) :: factorxp,factoryp,factorzp
    real(rp) :: gaccx,gaccy,gaccz,rhox_av,rhoy_av,rhoz_av
    real(rp) :: gravx,gravy,gravz,surfsumx,surfsumy,surfsumz,buoyx,buoyy,buoyz
    real(rp) :: momx,momy,momz,force(3,4)
    logical :: is_cmptx,is_cmpty,is_cmptz
    logical :: is_buoyant,is_brack,is_fs
    integer :: i,j,k,ierr
    !
    is_buoyant = is_boussinesq_buoyancy.and.present(s)
    is_brack = is_brackbill_normalization
    is_fs = present(fsx)
    rho = rho12(2); drho = rho12(1)-rho12(2)
    rho0i = 1._rp/rho0
    rhobeta = rho12(2)*beta12(2)
    drhobeta = rho12(1)*beta12(1)-rho12(2)*beta12(2)
    dxi = dli(1); dyi = dli(2)
    gaccx = gacc(1); gaccy = gacc(2); gaccz = gacc(3)
    rhox_av = rho_av(1); rhoy_av = rho_av(2); rhoz_av = rho_av(3)
    is_cmptx = is_cmpt(1); is_cmpty = is_cmpt(2); is_cmptz = is_cmpt(3)
    surf_factor = sigma
    if(is_brack) surf_factor = sigma*2._rp/(rho12(1)+rho12(2)) ! density-based scaling of Brackbill et al.
    gravx = 0._rp; gravy = 0._rp; gravz = 0._rp
    surfsumx = 0._rp; surfsumy = 0._rp; surfsumz = 0._rp
    buoyx = 0._rp; buoyy = 0._rp; buoyz = 0._rp
    momx = 0._rp; momy = 0._rp; momz = 0._rp
    !$acc        data copy(gravx,gravy,gravz,surfsumx,surfsumy,surfsumz,buoyx,buoyy,buoyz,momx,momy,momz) async(1)
    !$omp target data map(tofrom:gravx,gravy,gravz,surfsumx,surfsumy,surfsumz,buoyx,buoyy,buoyz,momx,momy,momz)
    !$acc parallel loop collapse(3) default(present) firstprivate(is_cmptx,is_cmpty,is_cmptz) &
    !$acc reduction(+:gravx,gravy,gravz,surfsumx,surfsumy,surfsumz,buoyx,buoyy,buoyz,momx,momy,momz) async(1)
    !$omp target teams loop collapse(3) firstprivate(is_cmptx,is_cmpty,is_cmptz) &
    !$omp reduction(+:gravx,gravy,gravz,surfsumx,surfsumy,surfsumz,buoyx,buoyy,buoyz,momx,momy,momz)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          dzci_c = dzci(k)
          weight_f = grid_vol_ratio_f(k)
          weight_c = grid_vol_ratio_c(k)
          !
          c_ccc = psi(i  ,j  ,k  )
          c_pcc = psi(i+1,j  ,k  )
          c_cpc = psi(i  ,j+1,k  )
          c_ccp = psi(i  ,j  ,k+1)
          co_ccc = psio(i  ,j  ,k  )
          co_pcc = psio(i+1,j  ,k  )
          co_cpc = psio(i  ,j+1,k  )
          co_ccp = psio(i  ,j  ,k+1)
          k_ccc = kappa(i  ,j  ,k  )
          k_pcc = kappa(i+1,j  ,k  )
          k_cpc = kappa(i  ,j+1,k  )
          k_ccp = kappa(i  ,j  ,k+1)
          !
          psixp = 0.5_rp*(c_pcc+c_ccc)
          psiyp = 0.5_rp*(c_cpc+c_ccc)
          psizp = wzc(k)*c_ccc+(1.-wzc(k))*c_ccp
          rhoxp = rho+drho*psixp
          rhoyp = rho+drho*psiyp
          rhozp = rho+drho*psizp
          rhox_old = rho+drho*0.5_rp*(co_pcc+co_ccc)
          rhoy_old = rho+drho*0.5_rp*(co_cpc+co_ccc)
          rhoz_old = rho+drho*(wzc(k)*co_ccc+(1.-wzc(k))*co_ccp)
          !
          skappaxp = 0.5_rp*(k_pcc+k_ccc)*surf_factor
          skappayp = 0.5_rp*(k_cpc+k_ccc)*surf_factor
          skappazp = 0.5_rp*(k_ccp+k_ccc)*surf_factor
          if(is_brack) then
            skappaxp = skappaxp*rhoxp
            skappayp = skappayp*rhoyp
            skappazp = skappazp*rhozp
          end if
          surfx = skappaxp*(c_pcc-c_ccc)*dxi
          surfy = skappayp*(c_cpc-c_ccc)*dyi
          surfz = skappazp*(c_ccp-c_ccc)*dzci_c
          if(is_fs) then
            surfx = fsx(i,j,k)
            surfy = fsy(i,j,k)
            surfz = fsz(i,j,k)
          end if
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
          surfx_e = (1._rp+dt_r)*surfx_n(i,j,k)-dt_r*surfx_o(i,j,k)
          surfy_e = (1._rp+dt_r)*surfy_n(i,j,k)-dt_r*surfy_o(i,j,k)
          surfz_e = (1._rp+dt_r)*surfz_n(i,j,k)-dt_r*surfz_o(i,j,k)
#if defined(_CONSTANT_COEFFS_POISSON)
          surfx_f = surfx_e+rhoxp*rho0i*(surfx-surfx_e)
          surfy_f = surfy_e+rhoyp*rho0i*(surfy-surfy_e)
          surfz_f = surfz_e+rhozp*rho0i*(surfz-surfz_e)
#else
          surfx_f = surfx
          surfy_f = surfy
          surfz_f = surfz
#endif
#else
          surfx_f = surfx
          surfy_f = surfy
          surfz_f = surfz
#endif
          if(is_buoyant) then
            s_ccc = s(i  ,j  ,k  )
            s_pcc = s(i+1,j  ,k  )
            s_cpc = s(i  ,j+1,k  )
            s_ccp = s(i  ,j  ,k+1)
            factorxp = rhobeta+drhobeta*psixp
            factoryp = rhobeta+drhobeta*psiyp
            factorzp = rhobeta+drhobeta*psizp
          end if
          !
          if(is_cmptx) then
            gravx = gravx+(rhoxp-rhox_av)*gaccx*weight_f
            surfsumx = surfsumx+surfx_f*weight_f
            momx = momx+rhox_old*u(i,j,k)*weight_f
            if(is_buoyant) buoyx = buoyx-gaccx*factorxp*0.5_rp*(s_pcc+s_ccc)*weight_f
          end if
          if(is_cmpty) then
            gravy = gravy+(rhoyp-rhoy_av)*gaccy*weight_f
            surfsumy = surfsumy+surfy_f*weight_f
            momy = momy+rhoy_old*v(i,j,k)*weight_f
            if(is_buoyant) buoyy = buoyy-gaccy*factoryp*0.5_rp*(s_cpc+s_ccc)*weight_f
          end if
          if(is_cmptz) then
            gravz = gravz+(rhozp-rhoz_av)*gaccz*weight_c
            surfsumz = surfsumz+surfz_f*weight_c
            momz = momz+rhoz_old*w(i,j,k)*weight_c
            if(is_buoyant) buoyz = buoyz-gaccz*factorzp*0.5_rp*(s_ccp+s_ccc)*weight_c
          end if
        end do
      end do
    end do
    !$omp end target data
    !$acc end        data
    !$acc wait(1)
    force(:,1) = [gravx,gravy,gravz]
    force(:,2) = [surfsumx,surfsumy,surfsumz]
    force(:,3) = [buoyx,buoyy,buoyz]
    force(:,4) = [momx,momy,momz]
    call MPI_ALLREDUCE(MPI_IN_PLACE,force(1,1),12,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    fgrav(:) = force(:,1)
    fsurf(:) = force(:,2)
    fbuoy(:) = force(:,3)
    momentum(:) = force(:,4)
    where(.not.is_cmpt(:))
      fgrav(:) = 0._rp
      fsurf(:) = 0._rp
      fbuoy(:) = 0._rp
      momentum(:) = 0._rp
    end where
  end subroutine cmpt_bodyforce
  !
  subroutine cmpt_momentum(n,nh_v,nh_psi,nh_g,is_cmpt,grid_vol_ratio_c,grid_vol_ratio_f,wzc,rho12,psi,u,v,w,momentum)
    !
    ! computes the volume-averaged mixture momentum on the staggered grid
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    integer , intent(in) :: nh_v,nh_psi,nh_g
    logical , intent(in), dimension(3) :: is_cmpt
    real(rp), intent(in), dimension(1-nh_g:) :: grid_vol_ratio_c,grid_vol_ratio_f,wzc
    real(rp), intent(in), dimension(2) :: rho12
    real(rp), intent(in), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(out), dimension(3) :: momentum
    real(rp) :: rho,drho,rhoxp,rhoyp,rhozp,weight_f,weight_c
    real(rp) :: momx,momy,momz
    logical :: is_cmptx,is_cmpty,is_cmptz
    integer :: i,j,k,ierr
    !
    rho = rho12(2)
    drho = rho12(1)-rho12(2)
    is_cmptx = is_cmpt(1); is_cmpty = is_cmpt(2); is_cmptz = is_cmpt(3)
    momx = 0._rp; momy = 0._rp; momz = 0._rp
    !$acc        data copy(momx,momy,momz) async(1)
    !$omp target data map(tofrom:momx,momy,momz)
    !$acc parallel loop collapse(3) default(present) firstprivate(is_cmptx,is_cmpty,is_cmptz) &
    !$acc reduction(+:momx,momy,momz) async(1)
    !$omp target teams loop collapse(3) firstprivate(is_cmptx,is_cmpty,is_cmptz) &
    !$omp reduction(+:momx,momy,momz)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          weight_f = grid_vol_ratio_f(k)
          weight_c = grid_vol_ratio_c(k)
          rhoxp = rho+drho*0.5_rp*(psi(i+1,j,k)+psi(i,j,k))
          rhoyp = rho+drho*0.5_rp*(psi(i,j+1,k)+psi(i,j,k))
          rhozp = rho+drho*(wzc(k)*psi(i,j,k)+(1.-wzc(k))*psi(i,j,k+1))
          if(is_cmptx) momx = momx+rhoxp*u(i,j,k)*weight_f
          if(is_cmpty) momy = momy+rhoyp*v(i,j,k)*weight_f
          if(is_cmptz) momz = momz+rhozp*w(i,j,k)*weight_c
        end do
      end do
    end do
    !$omp end target data
    !$acc end        data
    !$acc wait(1)
    momentum(:) = [momx,momy,momz]
    call MPI_ALLREDUCE(MPI_IN_PLACE,momentum,3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    where(.not.is_cmpt(:)) momentum(:) = 0._rp
  end subroutine cmpt_momentum
  !
  subroutine cmpt_wallshear(n,nh_v,nh_psi,nh_g,is_cmpt,is_bound,l,dli,dzci,dzfi,mu12,psi,u,v,w,taux,tauy,tauz)
    !
    ! n.b.: `is_bound` should exclude periodic BCs
    !
    implicit none
    integer , intent(in ), dimension(3) :: n
    integer , intent(in) :: nh_v,nh_psi,nh_g
    logical , intent(in ), dimension(3) :: is_cmpt
    logical , intent(in ), dimension(0:1,3) :: is_bound
    real(rp), intent(in ), dimension(3)     :: l,dli
    real(rp), intent(in ), dimension(1-nh_g:)    :: dzci,dzfi
    real(rp), intent(in ), dimension(2)     :: mu12
    real(rp), intent(in ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(out), dimension(3) :: taux,tauy,tauz
    real(rp) :: dudyp,dudym,dudzp,dudzm, &
                dvdxp,dvdxm,dvdzp,dvdzm, &
                dwdxp,dwdxm,dwdyp,dwdym, &
                muxp,muxm,muyp,muym,muzp,muzm
    !
    real(rp) :: tau21,tau31,tau12,tau32,tau13,tau23
    integer :: i,j,k,nx,ny,nz
    real(rp) :: dxi,dyi,lx,ly,lz,mu,dmu
    logical  :: is_harm
    real(rp) :: tau(3,3)
    integer :: ierr
    !
    is_harm = is_harmonic_interpolation
    if(is_harm) then
      mu  = 1._rp/mu12(2); dmu = 1._rp/mu12(1)-1._rp/mu12(2)
    else
      mu  = mu12(2)      ; dmu = mu12(1)-mu12(2)
    end if
    nx = n(1); ny = n(2); nz = n(3)
    dxi = dli(1); dyi = dli(2)
    lx = l(1); ly = l(2); lz = l(3)
    tau21 = 0._rp
    !$acc        data copy(tau21) async(1)
    !$omp target data map(tofrom:tau21)
    if(is_cmpt(1).and.is_bound(0,2)) then
      j = 0
      !$acc parallel loop collapse(2) default(present) private(dvdxp,dudyp,muyp) &
      !$acc reduction(+:tau21) async(1)
      !$omp target teams loop collapse(2) private(dvdxp,dudyp,muyp) &
      !$omp reduction(+:tau21)
      do k=1,nz
        do i=1,nx
          dvdxp = (v(i+1,j  ,k  )-v(i  ,j  ,k  ))*dxi
          dudyp = (u(i  ,j+1,k  )-u(i  ,j  ,k  ))*dyi
          muyp = mu+dmu*0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i+1,j+1,k)+psi(i+1,j,k))
          if(is_harm) muyp = 1._rp/muyp
          !
          tau21 = tau21 + (dudyp+dvdxp)*muyp/(dxi*dzfi(k)*lx*lz)
        end do
      end do
    end if
    if(is_cmpt(1).and.is_bound(1,2)) then
      j = ny+1
      !$acc parallel loop collapse(2) default(present) private(dvdxm,dudym,muym) &
      !$acc reduction(+:tau21) async(1)
      !$omp target teams loop collapse(2) private(dvdxm,dudym,muym) &
      !$omp reduction(+:tau21)
      do k=1,nz
        do i=1,nx
          dvdxm = (v(i+1,j-1,k  )-v(i  ,j-1,k  ))*dxi
          dudym = (u(i  ,j  ,k  )-u(i  ,j-1,k  ))*dyi
          muym = mu+dmu*0.25*(psi(i,j,k)+psi(i,j-1,k)+psi(i+1,j-1,k)+psi(i+1,j,k))
          if(is_harm) muym = 1._rp/muym
          !
          tau21 = tau21 - (dudym+dvdxm)*muym/(dxi*dzfi(k)*lx*lz)
        end do
      end do
    end if
    !$omp end target data
    !$acc end        data
    tau31 = 0._rp
    !$acc        data copy(tau31) async(1)
    !$omp target data map(tofrom:tau31)
    if(is_cmpt(1).and.is_bound(0,3)) then
      k = 0
      !$acc parallel loop collapse(2) default(present) private(dudzp,dwdxp,muzp) &
      !$acc reduction(+:tau31) async(1)
      !$omp target teams loop collapse(2) private(dudzp,dwdxp,muzp) &
      !$omp reduction(+:tau31)
      do j=1,ny
        do i=1,nx
          dudzp = (u(i  ,j  ,k+1)-u(i  ,j  ,k  ))*dzci(k  )
          dwdxp = (w(i+1,j  ,k  )-w(i  ,j  ,k  ))*dxi
          muzp = mu+dmu*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i+1,j,k+1)+psi(i+1,j,k))
          if(is_harm) muzp = 1._rp/muzp
          !
          tau31 = tau31 + (dudzp+dwdxp)*muzp/(dxi*dyi*lx*ly)
        end do
      end do
    end if
    if(is_cmpt(1).and.is_bound(1,3)) then
      k = nz+1
      !$acc parallel loop collapse(2) default(present) private(dudzm,dwdxm,muzm) &
      !$acc reduction(+:tau31) async(1)
      !$omp target teams loop collapse(2) private(dudzm,dwdxm,muzm) &
      !$omp reduction(+:tau31)
      do j=1,ny
        do i=1,nx
          dudzm = (u(i  ,j  ,k  )-u(i  ,j  ,k-1))*dzci(k-1)
          dwdxm = (w(i+1,j  ,k-1)-w(i  ,j  ,k-1))*dxi
          muzm = mu+dmu*0.25*(psi(i,j,k)+psi(i,j,k-1)+psi(i+1,j,k-1)+psi(i+1,j,k))
          if(is_harm) muzm = 1._rp/muzm
          !
          tau31 = tau31 - (dudzm+dwdxm)*muzm/(dxi*dyi*lx*ly)
        end do
      end do
    end if
    !$omp end target data
    !$acc end        data
    !
    tau12 = 0._rp
    !$acc        data copy(tau12) async(1)
    !$omp target data map(tofrom:tau12)
    if(is_cmpt(2).and.is_bound(0,1)) then
      i = 0
      !$acc parallel loop collapse(2) default(present) private(dvdxp,dudyp,muxp) &
      !$acc reduction(+:tau12) async(1)
      !$omp target teams loop collapse(2) private(dvdxp,dudyp,muxp) &
      !$omp reduction(+:tau12)
      do k=1,nz
        do j=1,ny
          dvdxp = (v(i+1,j  ,k  )-v(i  ,j  ,k  ))*dxi
          dudyp = (u(i  ,j+1,k  )-u(i  ,j  ,k  ))*dyi
          muxp = mu+dmu*0.25*(psi(i,j,k)+psi(i+1,j,k)+psi(i+1,j+1,k)+psi(i,j+1,k))
          if(is_harm) muxp = 1._rp/muxp
          !
          tau12 = tau12 + (dvdxp+dudyp)*muxp/(dyi*dzfi(k)*ly*lz)
        end do
      end do
    end if
    if(is_cmpt(2).and.is_bound(1,1)) then
      i = nx+1
      !$acc parallel loop collapse(2) default(present) private(dvdxm,dudym,muxm) &
      !$acc reduction(+:tau12) async(1)
      !$omp target teams loop collapse(2) private(dvdxm,dudym,muxm) &
      !$omp reduction(+:tau12)
      do k=1,nz
        do j=1,ny
          dvdxm = (v(i  ,j  ,k  )-v(i-1,j  ,k  ))*dxi
          dudym = (u(i-1,j+1,k  )-u(i-1,j  ,k  ))*dyi
          muxm = mu+dmu*0.25*(psi(i,j,k)+psi(i-1,j,k)+psi(i-1,j+1,k)+psi(i,j+1,k))
          if(is_harm) muxm = 1._rp/muxm
          !
          tau12 = tau12 - (dvdxm+dudym)*muxm/(dyi*dzfi(k)*ly*lz)
        end do
      end do
    end if
    !$omp end target data
    !$acc end        data
    tau32 = 0._rp
    !$acc        data copy(tau32) async(1)
    !$omp target data map(tofrom:tau32)
    if(is_cmpt(2).and.is_bound(0,3)) then
      k = 0
      !$acc parallel loop collapse(2) default(present) private(dvdzp,dwdyp,muzp) &
      !$acc reduction(+:tau32) async(1)
      !$omp target teams loop collapse(2) private(dvdzp,dwdyp,muzp) &
      !$omp reduction(+:tau32)
      do j=1,ny
        do i=1,nx
          dvdzp = (v(i  ,j  ,k+1)-v(i  ,j  ,k  ))*dzci(k  )
          dwdyp = (w(i  ,j+1,k  )-w(i  ,j  ,k  ))*dyi
          muzp = mu+dmu*0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i,j+1,k+1)+psi(i,j,k+1))
          if(is_harm) muzp = 1._rp/muzp
          !
          tau32 = tau32 + (dvdzp+dwdyp)*muzp/(dxi*dyi*lx*ly)
        end do
      end do
    end if
    if(is_cmpt(2).and.is_bound(1,3)) then
      k = nz+1
      !$acc parallel loop collapse(2) default(present) private(dvdzm,dwdym,muzm) &
      !$acc reduction(+:tau32) async(1)
      !$omp target teams loop collapse(2) private(dvdzm,dwdym,muzm) &
      !$omp reduction(+:tau32)
      do j=1,ny
        do i=1,nx
          dvdzm = (v(i  ,j  ,k  )-v(i  ,j  ,k-1))*dzci(k-1)
          dwdym = (w(i  ,j+1,k-1)-w(i  ,j  ,k-1))*dyi
          muzm = mu+dmu*0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i,j+1,k-1)+psi(i,j,k-1))
          if(is_harm) muzm = 1._rp/muzm
          !
          tau32 = tau32 - (dvdzm+dwdym)*muzm/(dxi*dyi*lx*ly)
        end do
      end do
    end if
    !$omp end target data
    !$acc end        data
    !
    tau13 = 0._rp
    !$acc        data copy(tau13) async(1)
    !$omp target data map(tofrom:tau13)
    if(is_cmpt(3).and.is_bound(0,1)) then
      i = 0
      !$acc parallel loop collapse(2) default(present) private(dwdxp,dudzp,muxp) &
      !$acc reduction(+:tau13) async(1)
      !$omp target teams loop collapse(2) private(dwdxp,dudzp,muxp) &
      !$omp reduction(+:tau13)
      do k=1,nz
        do j=1,ny
          dwdxp = (w(i+1,j  ,k  )-w(i  ,j  ,k))*dxi
          dudzp = (u(i  ,j  ,k+1)-u(i  ,j  ,k))*dzci(k  )
          muxp = mu+dmu*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i+1,j,k+1)+psi(i+1,j,k))
          if(is_harm) muxp = 1._rp/muxp
          !
          tau13 = tau13 + (dwdxp+dudzp)*muxp/(dyi*dzfi(k)*ly*lz)
        end do
      end do
    end if
    if(is_cmpt(3).and.is_bound(1,1)) then
      i = nx+1
      !$acc parallel loop collapse(2) default(present) private(dwdxm,dudzm,muxm) &
      !$acc reduction(+:tau13) async(1)
      !$omp target teams loop collapse(2) private(dwdxm,dudzm,muxm) &
      !$omp reduction(+:tau13)
      do k=1,nz
        do j=1,ny
          dwdxm = (w(i  ,j  ,k  )-w(i-1,j  ,k))*dxi
          dudzm = (u(i-1,j  ,k+1)-u(i-1,j  ,k))*dzci(k  )
          muxm = mu+dmu*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i-1,j,k+1)+psi(i-1,j,k))
          if(is_harm) muxm = 1._rp/muxm
          !
          tau13 = tau13 - (dwdxm+dudzm)*muxm/(dyi*dzfi(k)*ly*lz)
        end do
      end do
    end if
    !$omp end target data
    !$acc end        data
    tau23 = 0._rp
    !$acc        data copy(tau23) async(1)
    !$omp target data map(tofrom:tau23)
    if(is_cmpt(3).and.is_bound(0,2)) then
      j = 0
      !$acc parallel loop collapse(2) default(present) private(dwdyp,dvdzp,muyp) &
      !$acc reduction(+:tau23) async(1)
      !$omp target teams loop collapse(2) private(dwdyp,dvdzp,muyp) &
      !$omp reduction(+:tau23)
      do k=1,nz
        do i=1,nx
          dwdyp = (w(i  ,j+1,k  )-w(i  ,j  ,k))*dyi
          dvdzp = (v(i  ,j  ,k+1)-v(i  ,j  ,k))*dzci(k  )
          muyp = mu+dmu*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i,j+1,k+1)+psi(i,j+1,k))
          if(is_harm) muyp = 1._rp/muyp
          !
          tau23 = tau23 + (dwdyp+dvdzp)*muyp/(dxi*dzfi(k)*lx*lz)
        end do
      end do
    end if
    if(is_cmpt(3).and.is_bound(1,2)) then
      j = ny+1
      !$acc parallel loop collapse(2) default(present) private(dwdym,dvdzm,muym) &
      !$acc reduction(+:tau23) async(1)
      !$omp target teams loop collapse(2) private(dwdym,dvdzm,muym) &
      !$omp reduction(+:tau23)
      do k=1,nz
        do i=1,nx
          dwdym = (w(i  ,j  ,k  )-w(i  ,j-1,k))*dyi
          dvdzm = (v(i  ,j-1,k+1)-v(i  ,j-1,k))*dzci(k  )
          muym = mu+dmu*0.25*(psi(i,j,k)+psi(i,j,k+1)+psi(i,j-1,k+1)+psi(i,j-1,k))
          if(is_harm) muym = 1._rp/muym
          !
          tau23 = tau23 - (dwdym+dvdzm)*muym/(dxi*dzfi(k)*lx*lz)
        end do
      end do
    end if
    !$omp end target data
    !$acc end        data
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
  !
  subroutine cmpt_capillary_force(n,nh_psi,nh_ifc,nh_g,dli,dzci,wzc,sigma,rho12,psi,psi_m,kappa_m,is_init,fsx,fsy,fsz)
    !
    ! surface tension force at the faces of one marker, sigma*kappa_m*grad(psi_m), with the density weighting
    ! of `mom_xyz_oth` (from the total phase indicator psi) if `is_brackbill_normalization`, exactly as it is
    ! computed there for a single marker; the force is stored in fsx, fsy, fsz if `is_init`, and added to them
    ! otherwise: with several markers, the forces of all of them are summed, such that interfaces of different
    ! markers that share a cell keep their own curvature
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   ) :: nh_psi,nh_ifc,nh_g
    real(rp), intent(in   ), dimension(3) :: dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzci,wzc
    real(rp), intent(in   ) :: sigma
    real(rp), intent(in   ), dimension(2) :: rho12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi,psi_m
    real(rp), intent(in   ), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: kappa_m
    logical , intent(in   ) :: is_init
    real(rp), intent(inout), dimension(:,:,:) :: fsx,fsy,fsz
    real(rp) :: rho,drho,dxi,dyi,dzci_c,surf_factor,skappaxp,skappayp,skappazp,surfx,surfy,surfz
    logical  :: is_brack
    integer  :: i,j,k
    !
    is_brack = is_brackbill_normalization
    rho = rho12(2); drho = rho12(1)-rho12(2)
    dxi = dli(1)
    dyi = dli(2)
    surf_factor = sigma
    if(is_brack) surf_factor = sigma*2._rp/(rho12(1)+rho12(2))
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(dzci_c,skappaxp,skappayp,skappazp,surfx,surfy,surfz) async(1)
    !$omp target teams loop collapse(3) &
    !$omp private(dzci_c,skappaxp,skappayp,skappazp,surfx,surfy,surfz)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          dzci_c = dzci(k)
          skappaxp = 0.5*(kappa_m(i+1,j,k)+kappa_m(i,j,k))*surf_factor
          skappayp = 0.5*(kappa_m(i,j+1,k)+kappa_m(i,j,k))*surf_factor
          skappazp = 0.5*(kappa_m(i,j,k+1)+kappa_m(i,j,k))*surf_factor
          if(is_brack) then
            skappaxp = skappaxp*(rho + drho*0.5*(psi(i+1,j,k)+psi(i,j,k)))
            skappayp = skappayp*(rho + drho*0.5*(psi(i,j+1,k)+psi(i,j,k)))
            skappazp = skappazp*(rho + drho*(wzc(k)*psi(i,j,k)+(1.-wzc(k))*psi(i,j,k+1)))
          end if
          surfx = skappaxp*(psi_m(i+1,j,k)-psi_m(i,j,k))*dxi
          surfy = skappayp*(psi_m(i,j+1,k)-psi_m(i,j,k))*dyi
          surfz = skappazp*(psi_m(i,j,k+1)-psi_m(i,j,k))*dzci_c
          if(is_init) then
            fsx(i,j,k) = surfx
            fsy(i,j,k) = surfy
            fsz(i,j,k) = surfz
          else
            fsx(i,j,k) = fsx(i,j,k) + surfx
            fsy(i,j,k) = fsy(i,j,k) + surfy
            fsz(i,j,k) = fsz(i,j,k) + surfz
          end if
        end do
      end do
    end do
  end subroutine cmpt_capillary_force
  !
  subroutine mom_xyz_marangoni(n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,dli,dzci,wzc,rho12,mu12,dsigma,s_ref,is_brack,psi_p,psi,kappa,s, &
                               factor,u,v,w)
    !
    ! capillary force of the variation of the surface tension with the first scalar,
    ! sigma(s) = sigma + dsigma*(s-s_ref), integrated over the time `factor` and added to the velocity, divided
    ! by the face density; it complements the force of the constant part, sigma, of `mom_xyz_oth`:
    !  - normal force of the variation, (sigma(s)-sigma)*kappa*grad(psi), with the density weighting of the
    !    constant part if `is_brack`;
    !  - tangential (Marangoni) force, (grad(sigma) - n (n.grad(sigma)))*|grad(psi)|, with the gradient of psi
    !    at the faces: the difference across the face along its normal, and the average of the centred
    !    differences of the two cells along the other directions; this force is weighted with the local
    !    viscosity over the mean one: it is balanced by viscous stresses, and an even split between the two
    !    fluids makes the less viscous one slip at the interface (the analogue, for the tangential force,
    !    of the density weighting of the normal force)
    ! the fluid properties follow the phase indicator psi_p, and the interface is that of psi, with curvature
    ! kappa: the same one, or, with several markers, the total phase indicator and that of one marker
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   ) :: nh_v,nh_s,nh_psi,nh_ifc,nh_g
    real(rp), intent(in   ), dimension(3) :: dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzci,wzc
    real(rp), intent(in   ), dimension(2) :: rho12,mu12
    real(rp), intent(in   ) :: dsigma,s_ref
    logical , intent(in   ) :: is_brack
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi_p,psi
    real(rp), intent(in   ), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: kappa
    real(rp), intent(in   ), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s
    real(rp), intent(in   ) :: factor
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp) :: rho,drho,mu,dmu,bfac,mfac,dxi,dyi,hz_c
    real(rp) :: fx,fy,fz,ga,gb,gc,sa,sb,sc,gm2,psif,c_min,c_max
    real(rp), parameter :: c_eps = max(1.e-12_rp,1.e3_rp*epsilon(1._rp)) ! (as for the reconstruction and the limiter)
    integer  :: i,j,k
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    mu  = mu12(2) ; dmu  = mu12(1) -mu12(2)
    bfac = 2._rp/(rho12(1)+rho12(2))
    mfac = 2._rp/(mu12(1) +mu12(2) )
    dxi = dli(1); dyi = dli(2)
    !$acc parallel     loop collapse(3) default(present) private(hz_c,fx,fy,fz,ga,gb,gc,sa,sb,sc,gm2,psif,c_min,c_max) async(1)
    !$omp target teams loop collapse(3)                  private(hz_c,fx,fy,fz,ga,gb,gc,sa,sb,sc,gm2,psif,c_min,c_max)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          ! the force vanishes where the phase indicator is uniform over the stencil (to the tolerance of full
          ! and empty cells): these cells, most of the domain, are skipped
          !
          c_min = min(psi(i  ,j  ,k  ),psi(i+1,j  ,k  ),psi(i-1,j  ,k  ),psi(i  ,j+1,k  ),psi(i  ,j-1,k  ), &
                      psi(i+1,j+1,k  ),psi(i+1,j-1,k  ),psi(i-1,j+1,k  ), &
                      psi(i  ,j  ,k+1),psi(i+1,j  ,k+1),psi(i-1,j  ,k+1),psi(i  ,j+1,k+1),psi(i  ,j-1,k+1), &
                      psi(i  ,j  ,k-1),psi(i+1,j  ,k-1),psi(i  ,j+1,k-1))
          c_max = max(psi(i  ,j  ,k  ),psi(i+1,j  ,k  ),psi(i-1,j  ,k  ),psi(i  ,j+1,k  ),psi(i  ,j-1,k  ), &
                      psi(i+1,j+1,k  ),psi(i+1,j-1,k  ),psi(i-1,j+1,k  ), &
                      psi(i  ,j  ,k+1),psi(i+1,j  ,k+1),psi(i-1,j  ,k+1),psi(i  ,j+1,k+1),psi(i  ,j-1,k+1), &
                      psi(i  ,j  ,k-1),psi(i+1,j  ,k-1),psi(i  ,j+1,k-1))
          if(c_max-c_min <= c_eps) cycle
          hz_c = dzci(k)*dzci(k-1)/(dzci(k)+dzci(k-1)) ! (factor of the centred difference along z)
          !
          ! x face
          !
          psif = 0.5*(psi_p(i+1,j,k)+psi_p(i,j,k))
          ga = (psi(i+1,j,k)-psi(i,j,k))*dxi
          gb = 0.25*((psi(i,j+1,k)-psi(i,j-1,k))+(psi(i+1,j+1,k)-psi(i+1,j-1,k)))*dyi
          gc = 0.5*((psi(i,j,k+1)-psi(i,j,k-1))+(psi(i+1,j,k+1)-psi(i+1,j,k-1)))*hz_c
          sa = (s(i+1,j,k)-s(i,j,k))*dxi
          sb = 0.25*((s(i,j+1,k)-s(i,j-1,k))+(s(i+1,j+1,k)-s(i+1,j-1,k)))*dyi
          sc = 0.5*((s(i,j,k+1)-s(i,j,k-1))+(s(i+1,j,k+1)-s(i+1,j,k-1)))*hz_c
          gm2 = ga**2+gb**2+gc**2
          fx = dsigma*(0.5*(s(i+1,j,k)+s(i,j,k))-s_ref)*0.5*(kappa(i+1,j,k)+kappa(i,j,k))*ga
          if(is_brack) fx = fx*(rho+drho*psif)*bfac
          if(gm2 > 0._rp) fx = fx + dsigma*(sa*gm2-ga*(ga*sa+gb*sb+gc*sc))/sqrt(gm2)*(mu+dmu*psif)*mfac
          u(i,j,k) = u(i,j,k) + factor*fx/(rho+drho*psif)
          !
          ! y face
          !
          psif = 0.5*(psi_p(i,j+1,k)+psi_p(i,j,k))
          ga = (psi(i,j+1,k)-psi(i,j,k))*dyi
          gb = 0.25*((psi(i+1,j,k)-psi(i-1,j,k))+(psi(i+1,j+1,k)-psi(i-1,j+1,k)))*dxi
          gc = 0.5*((psi(i,j,k+1)-psi(i,j,k-1))+(psi(i,j+1,k+1)-psi(i,j+1,k-1)))*hz_c
          sa = (s(i,j+1,k)-s(i,j,k))*dyi
          sb = 0.25*((s(i+1,j,k)-s(i-1,j,k))+(s(i+1,j+1,k)-s(i-1,j+1,k)))*dxi
          sc = 0.5*((s(i,j,k+1)-s(i,j,k-1))+(s(i,j+1,k+1)-s(i,j+1,k-1)))*hz_c
          gm2 = ga**2+gb**2+gc**2
          fy = dsigma*(0.5*(s(i,j+1,k)+s(i,j,k))-s_ref)*0.5*(kappa(i,j+1,k)+kappa(i,j,k))*ga
          if(is_brack) fy = fy*(rho+drho*psif)*bfac
          if(gm2 > 0._rp) fy = fy + dsigma*(sa*gm2-ga*(ga*sa+gb*sb+gc*sc))/sqrt(gm2)*(mu+dmu*psif)*mfac
          v(i,j,k) = v(i,j,k) + factor*fy/(rho+drho*psif)
          !
          ! z face
          !
          psif = (wzc(k)*psi_p(i,j,k)+(1.-wzc(k))*psi_p(i,j,k+1))
          ga = (psi(i,j,k+1)-psi(i,j,k))*dzci(k)
          gb = 0.25*((psi(i+1,j,k)-psi(i-1,j,k))+(psi(i+1,j,k+1)-psi(i-1,j,k+1)))*dxi
          gc = 0.25*((psi(i,j+1,k)-psi(i,j-1,k))+(psi(i,j+1,k+1)-psi(i,j-1,k+1)))*dyi
          sa = (s(i,j,k+1)-s(i,j,k))*dzci(k)
          sb = 0.25*((s(i+1,j,k)-s(i-1,j,k))+(s(i+1,j,k+1)-s(i-1,j,k+1)))*dxi
          sc = 0.25*((s(i,j+1,k)-s(i,j-1,k))+(s(i,j+1,k+1)-s(i,j-1,k+1)))*dyi
          gm2 = ga**2+gb**2+gc**2
          fz = dsigma*(0.5*(s(i,j,k+1)+s(i,j,k))-s_ref)*0.5*(kappa(i,j,k+1)+kappa(i,j,k))*ga
          if(is_brack) fz = fz*(rho+drho*psif)*bfac
          if(gm2 > 0._rp) fz = fz + dsigma*(sa*gm2-ga*(ga*sa+gb*sb+gc*sc))/sqrt(gm2)*(mu+dmu*psif)*mfac
          w(i,j,k) = w(i,j,k) + factor*fz/(rho+drho*psif)
        end do
      end do
    end do
  end subroutine mom_xyz_marangoni
end module mod_mom
