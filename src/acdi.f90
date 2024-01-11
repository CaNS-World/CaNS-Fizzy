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
  use mod_param , only: eps,seps_factor,gam_factor
  implicit none
  private
  public set_seps,set_gam,pf
  contains
  subroutine set_seps(dl,dzfi,seps)
    !
    ! computes seps (only once because the grid is not changing)
    !
    implicit none
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: dzfi
    real(rp), intent(out) :: seps
    real(rp), save :: dlmax
    !
    dlmax  = max(maxval(dl(1:2)),maxval(1./dzfi))
    call MPI_ALLREDUCE(MPI_IN_PLACE,dlmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    seps = dlmax*seps_factor
  end subroutine set_seps
  !
  subroutine set_gam(n,u,v,w,gam)
    !
    ! computes gam (every icheck iterations)
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out) :: gam
    real(rp) :: uc,vc,wc,vel,velmax
    integer :: i,j,k
    !
    !$acc data copy(vel) async(1)
    !$acc parallel loop collapse(3) default(present) private(uc,vc,wc) reduction(max:vel) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(uc,vc,wc) REDUCTION(max:vel)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          uc  = 0.5*(u(i,j,k) + u(i-1,j,k))
          vc  = 0.5*(v(i,j,k) + v(i,j-1,k))
          wc  = 0.5*(w(i,j,k) + w(i,j,k-1))
          vel = sqrt(uc**2 + vc**2 + wc**2)
        end do
      end do
    end do
    !$acc end data
    !$acc wait(1)
    call MPI_ALLREDUCE(MPI_IN_PLACE,vel,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    velmax = vel
    gam = velmax*gam_factor
  end subroutine set_gam
  !
  subroutine pf(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,gam,seps,u,v,w,psi,dpsidt)
    !
    !
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,dzi,gam,seps
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w
    real(rp), dimension(0:,0:,0:), intent(inout) :: psi
    real(rp), dimension(:,:,:), intent(out) :: dpsidt
    integer :: i,j,k
    real(rp) :: sepsi 
    real(rp) :: wghtpp,wghtpm,wghtmp,wghtmm
    real(rp) :: ADV,DIFF,SHARP
    real(rp) :: upsiip,upsiim,vpsijp,vpsijm,wpsikp,wpsikm
    real(rp) :: dpsidxp,dpsidxm,dpsidyp,dpsidym,dpsidzp,dpsidzm
    real(rp) :: phimcm,phicmm,phiccm,phicpm,phipcm,phimmc,phimcc,phimpc,phicmc,phiccc, &
                phicpc,phipmc,phipcc,phippc,phimcp,phicmp,phiccp,phicpp,phipcp
    real(rp) :: dphidxleft,dphidyleft,dphidzleft,dphidxright,dphidyright,dphidzright, &
                dphidxback,dphidyback,dphidzback,dphidxfront,dphidyfront,dphidzfront, &
                dphidxbottom,dphidybottom,dphidzbottom,dphidxtop,dphidytop,dphidztop
    real(rp) :: nleft,nright,nback,nfront,nbottom,ntop
    
    sepsi = 1./seps
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          ! weights on the stretched grid direction
          wghtpp = dzci(k)/dzfi(k+1)
          wghtpm = dzci(k)/dzfi(k)
          wghtmp = dzci(k-1)/dzfi(k)
          wghtmm = dzci(k-1)/dzfi(k-1)
          !wghtpp = 1. 
          !wghtpm = 1. 
          !wghtmp = 1. 
          !wghtmm = 1. 
          !
          ! clip the field
          !psi(i,j,k) = min(max(0.,psi(i,j,k)),1.)
          !
          ! ADVECTION
          !
          upsiim  = 0.5*( psi(i-1,j,k)+psi(i,j,k) )*u(i-1,j,k)
          upsiip  = 0.5*( psi(i+1,j,k)+psi(i,j,k) )*u(i  ,j,k)
          vpsijm  = 0.5*( psi(i,j-1,k)+psi(i,j,k) )*v(i,j-1,k)
          vpsijp  = 0.5*( psi(i,j+1,k)+psi(i,j,k) )*v(i,j  ,k)
          wpsikm  = 0.5*( psi(i,j,k-1)*wghtmp + psi(i,j,k)*wghtmm )*w(i,j,k-1) 
          wpsikp  = 0.5*( psi(i,j,k+1)*wghtpm + psi(i,j,k)*wghtpp )*w(i,j,k  )
          ADV     = dxi*( upsiip - upsiim ) + dyi*( vpsijp - vpsijm ) + dzfi(k)*( wpsikp - wpsikm )
          !
          ! DIFFUSION
          !
          dpsidxp = (psi(i+1,j,k)-psi(i  ,j,k))*dxi
          dpsidxm = (psi(i  ,j,k)-psi(i-1,j,k))*dxi
          dpsidyp = (psi(i,j+1,k)-psi(i,j  ,k))*dyi
          dpsidym = (psi(i,j  ,k)-psi(i,j-1,k))*dyi
          dpsidzp = (psi(i,j,k+1)-psi(i,j,k  ))*dzci(k  )
          dpsidzm = (psi(i,j,k  )-psi(i,j,k-1))*dzci(k-1)
          DIFF    = gam*seps*( (dpsidxp-dpsidxm)*dxi + (dpsidyp-dpsidym)*dyi + (dpsidzp-dpsidzm)*dzfi(k) )
          !
          ! SHARPENING
          ! 
          ! calculate the surface normal at the cell faces
         !phimmm = seps*log((psi(i-1,j-1,k-1)+eps)/(1.-psi(i-1,j-1,k-1)+eps))  
          phimcm = seps*log((psi(i-1,j  ,k-1)+eps)/(1.-psi(i-1,j  ,k-1)+eps))  
         !phimpm = seps*log((psi(i-1,j+1,k-1)+eps)/(1.-psi(i-1,j+1,k-1)+eps))  
          phicmm = seps*log((psi(i  ,j-1,k-1)+eps)/(1.-psi(i  ,j-1,k-1)+eps))  
          phiccm = seps*log((psi(i  ,j  ,k-1)+eps)/(1.-psi(i  ,j  ,k-1)+eps))  
          phicpm = seps*log((psi(i  ,j+1,k-1)+eps)/(1.-psi(i  ,j+1,k-1)+eps))  
         !phipmm = seps*log((psi(i+1,j-1,k-1)+eps)/(1.-psi(i+1,j-1,k-1)+eps))  
          phipcm = seps*log((psi(i+1,j  ,k-1)+eps)/(1.-psi(i+1,j  ,k-1)+eps))  
         !phippm = seps*log((psi(i+1,j+1,k-1)+eps)/(1.-psi(i+1,j+1,k-1)+eps))  
          phimmc = seps*log((psi(i-1,j-1,k  )+eps)/(1.-psi(i-1,j-1,k  )+eps))  
          phimcc = seps*log((psi(i-1,j  ,k  )+eps)/(1.-psi(i-1,j  ,k  )+eps))  
          phimpc = seps*log((psi(i-1,j+1,k  )+eps)/(1.-psi(i-1,j+1,k  )+eps))  
          phicmc = seps*log((psi(i  ,j-1,k  )+eps)/(1.-psi(i  ,j-1,k  )+eps))  
          phiccc = seps*log((psi(i  ,j  ,k  )+eps)/(1.-psi(i  ,j  ,k  )+eps))  
          phicpc = seps*log((psi(i  ,j+1,k  )+eps)/(1.-psi(i  ,j+1,k  )+eps))  
          phipmc = seps*log((psi(i+1,j-1,k  )+eps)/(1.-psi(i+1,j-1,k  )+eps))  
          phipcc = seps*log((psi(i+1,j  ,k  )+eps)/(1.-psi(i+1,j  ,k  )+eps))  
          phippc = seps*log((psi(i+1,j+1,k  )+eps)/(1.-psi(i+1,j+1,k  )+eps))  
         !phimmp = seps*log((psi(i-1,j-1,k+1)+eps)/(1.-psi(i-1,j-1,k+1)+eps))  
          phimcp = seps*log((psi(i-1,j  ,k+1)+eps)/(1.-psi(i-1,j  ,k+1)+eps))  
         !phimpp = seps*log((psi(i-1,j+1,k+1)+eps)/(1.-psi(i-1,j+1,k+1)+eps))  
          phicmp = seps*log((psi(i  ,j-1,k+1)+eps)/(1.-psi(i  ,j-1,k+1)+eps))  
          phiccp = seps*log((psi(i  ,j  ,k+1)+eps)/(1.-psi(i  ,j  ,k+1)+eps))  
          phicpp = seps*log((psi(i  ,j+1,k+1)+eps)/(1.-psi(i  ,j+1,k+1)+eps))  
         !phipmp = seps*log((psi(i+1,j-1,k+1)+eps)/(1.-psi(i+1,j-1,k+1)+eps))  
          phipcp = seps*log((psi(i+1,j  ,k+1)+eps)/(1.-psi(i+1,j  ,k+1)+eps))  
         !phippp = seps*log((psi(i+1,j+1,k+1)+eps)/(1.-psi(i+1,j+1,k+1)+eps))  
          dphidxleft   = 1.00*(phiccc       -phimcc       )*dxi
          dphidyleft   = 0.25*(phicpc+phimpc-phicmc-phimmc)*dyi
          dphidzleft   = 0.50*(phiccp+phimcp-phiccm-phimcm)*(dzci(k-1)*dzci(k))/(dzci(k-1)+dzci(k))
          dphidxright  = 1.00*(phipcc       -phiccc       )*dxi
          dphidyright  = 0.25*(phippc+phicpc-phipmc-phicmc)*dyi
          dphidzright  = 0.50*(phipcp+phiccp-phipcm-phiccm)*(dzci(k-1)*dzci(k))/(dzci(k-1)+dzci(k))
          dphidxback   = 0.25*(phippc+phipcc-phimpc-phimcc)*dxi
          dphidyback   = 1.00*(phicpc       -phiccc       )*dyi
          dphidzback   = 0.50*(phicpp+phiccp-phicpm-phiccm)*(dzci(k-1)*dzci(k))/(dzci(k-1)+dzci(k))
          dphidxfront  = 0.25*(phipcc+phipmc-phimcc-phimmc)*dxi
          dphidyfront  = 1.00*(phiccc       -phicmc       )*dyi
          dphidzfront  = 0.50*(phiccp+phicmp-phiccm-phicmm)*(dzci(k-1)*dzci(k))/(dzci(k-1)+dzci(k))
          dphidxbottom = 0.25*((phipcm-phimcm)*wghtmp+(phipcc-phimcc)*wghtmm)*dxi
          dphidybottom = 0.25*((phicpm-phicmm)*wghtmp+(phicpc-phicmc)*wghtmm)*dyi
          dphidzbottom = 1.00*( phiccc               -phiccm                )*dzci(k-1)
          dphidxtop    = 0.25*((phipcp-phimcp)*wghtpm+(phipcc-phimcc)*wghtpp)*dxi
          dphidytop    = 0.25*((phicpp-phicmp)*wghtpm+(phicpc-phicmc)*wghtpp)*dyi
          dphidztop    = 1.00*( phiccp               -phiccc                )*dzci(k)
          nleft   = dphidxleft  / (sqrt(  dphidxleft**2 +   dphidyleft**2 +   dphidzleft**2)+eps) 
          nright  = dphidxright / (sqrt( dphidxright**2 +  dphidyright**2 +  dphidzright**2)+eps)
          nback   = dphidyback  / (sqrt(  dphidxback**2 +   dphidyback**2 +   dphidzback**2)+eps)
          nfront  = dphidyfront / (sqrt( dphidxfront**2 +  dphidyfront**2 +  dphidzfront**2)+eps)       
          nbottom = dphidzbottom/ (sqrt(dphidxbottom**2 + dphidybottom**2 + dphidzbottom**2)+eps)
          ntop    = dphidztop   / (sqrt(   dphidxtop**2 +    dphidytop**2 +    dphidztop**2)+eps)
          SHARP = 0.25*gam*((1.-(tanh(0.25*(phipcc       +phiccc       )*sepsi))**2)*nright-&
                            (1.-(tanh(0.25*(phiccc       +phimcc       )*sepsi))**2)*nleft  )*dxi + &
                  0.25*gam*((1.-(tanh(0.25*(phicpc       +phiccc       )*sepsi))**2)*nback -&
                            (1.-(tanh(0.25*(phiccc       +phicmc       )*sepsi))**2)*nfront )*dyi + &
                  0.25*gam*((1.-(tanh(0.25*(phiccp*wghtpm+phiccc*wghtpp)*sepsi))**2)*ntop  -&
                            (1.-(tanh(0.25*(phiccc*wghtmm+phiccm*wghtmp)*sepsi))**2)*nbottom)*dzfi(k)
          !
          ! FULL TRANSPORT
          !
          dpsidt(i,j,k) = -ADV+DIFF-SHARP
        end do
      end do
    end do
  end subroutine pf
!
end module mod_acdi
