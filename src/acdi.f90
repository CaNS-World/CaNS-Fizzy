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
  subroutine pf(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,gam,seps,u,v,w,phi,dphidt)
    !
    !
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,dzi,gam,seps
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w
    real(rp), dimension(0:,0:,0:), intent(inout) :: phi
    real(rp), dimension(:,:,:), intent(out) :: dphidt
    integer :: i,j,k
    real(rp) :: sepsi 
    real(rp) :: wghtpp,wghtpm,wghtmp,wghtmm
    real(rp) :: ADV,DIFF,SHARP
    real(rp) :: uphiip,uphiim,vphijp,vphijm,wphikp,wphikm
    real(rp) :: dphidxp,dphidxm,dphidyp,dphidym,dphidzp,dphidzm
    real(rp) :: psimcm,psicmm,psiccm,psicpm,psipcm,psimmc,psimcc,psimpc,psicmc,psiccc, &
                psicpc,psipmc,psipcc,psippc,psimcp,psicmp,psiccp,psicpp,psipcp
    real(rp) :: dpsidxleft,dpsidyleft,dpsidzleft,dpsidxright,dpsidyright,dpsidzright, &
                dpsidxback,dpsidyback,dpsidzback,dpsidxfront,dpsidyfront,dpsidzfront, &
                dpsidxbottom,dpsidybottom,dpsidzbottom,dpsidxtop,dpsidytop,dpsidztop
    real(rp) :: nleft,nright,nback,nfront,nbottom,ntop
    
    sepsi = 1./seps
    !
    !$acc parallel loop collapse(3) default(present) private(wghtpp,wghtpm,wghtmp,wghtmm) &
    !$acc private(uphiip,uphiim,vphijp,vphijm,wphikp,wphikm,ADV) &
    !$acc private(dphidxp,dphidxm,dphidyp,dphidym,dphidzp,dphidzm,DIFF) &
    !$acc private(psimcm,psicmm,psiccm,psicpm,psipcm,psimmc,psimcc,psimpc,psicmc,psiccc) &
    !$acc private(psicpc,psipmc,psipcc,psippc,psimcp,psicmp,psiccp,psicpp,psipcp) &
    !$acc private(dpsidxleft,dpsidyleft,dpsidzleft,dpsidxright,dpsidyright,dpsidzright) &
    !$acc private(dpsidxback,dpsidyback,dpsidzback,dpsidxfront,dpsidyfront,dpsidzfront) &
    !$acc private(dpsidxbottom,dpsidybottom,dpsidzbottom,dpsidxtop,dpsidytop,dpsidztop) &
    !$acc private(nleft,nright,nback,nfront,nbottom,ntop,SHARP) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(wghtpp,wghtpm,wghtmp,wghtmm) &
    !$OMP private(uphiip,uphiim,vphijp,vphijm,wphikp,wphikm,ADV) &
    !$OMP private(dphidxp,dphidxm,dphidyp,dphidym,dphidzp,dphidzm,DIFF) &
    !$OMP private(psimcm,psicmm,psiccm,psicpm,psipcm,psimmc,psimcc,psimpc,psicmc,psiccc) &
    !$OMP private(psicpc,psipmc,psipcc,psippc,psimcp,psicmp,psiccp,psicpp,psipcp) &
    !$OMP private(dpsidxleft,dpsidyleft,dpsidzleft,dpsidxright,dpsidyright,dpsidzright) &
    !$OMP private(dpsidxback,dpsidyback,dpsidzback,dpsidxfront,dpsidyfront,dpsidzfront) &
    !$OMP private(dpsidxbottom,dpsidybottom,dpsidzbottom,dpsidxtop,dpsidytop,dpsidztop) &
    !$OMP private(nleft,nright,nback,nfront,nbottom,ntop,SHARP)
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
          !phi(i,j,k) = min(max(0.,phi(i,j,k)),1.)
          !
          ! ADVECTION
          !
          uphiim  = 0.5*( phi(i-1,j,k)+phi(i,j,k) )*u(i-1,j,k)
          uphiip  = 0.5*( phi(i+1,j,k)+phi(i,j,k) )*u(i  ,j,k)
          vphijm  = 0.5*( phi(i,j-1,k)+phi(i,j,k) )*v(i,j-1,k)
          vphijp  = 0.5*( phi(i,j+1,k)+phi(i,j,k) )*v(i,j  ,k)
          wphikm  = 0.5*( phi(i,j,k-1)*wghtmp + phi(i,j,k)*wghtmm )*w(i,j,k-1) 
          wphikp  = 0.5*( phi(i,j,k+1)*wghtpm + phi(i,j,k)*wghtpp )*w(i,j,k  )
          ADV     = dxi*( uphiip - uphiim ) + dyi*( vphijp - vphijm ) + dzfi(k)*( wphikp - wphikm )
          ! upwind
          !ADV     = dxi*( 1.*(phi(i,j,k)-phi(i-1,j,k)) ) + dyi*( 0.*(phi(i,j,k)-phi(i,j-1,k)) )
          !
          ! DIFFUSION
          !
          dphidxp = (phi(i+1,j,k)-phi(i  ,j,k))*dxi
          dphidxm = (phi(i  ,j,k)-phi(i-1,j,k))*dxi
          dphidyp = (phi(i,j+1,k)-phi(i,j  ,k))*dyi
          dphidym = (phi(i,j  ,k)-phi(i,j-1,k))*dyi
          dphidzp = (phi(i,j,k+1)-phi(i,j,k  ))*dzci(k  )
          dphidzm = (phi(i,j,k  )-phi(i,j,k-1))*dzci(k-1)
          DIFF    = gam*seps*( (dphidxp-dphidxm)*dxi + (dphidyp-dphidym)*dyi + (dphidzp-dphidzm)*dzfi(k) )
          !
          ! SHARPENING
          ! 
          ! calculate the surface normal at the cell faces
          !!!!if (phi(i+1,j+1,k) < 0.) then
          !!!if (ANY (phi < 0.)) then
          !!!  print*, 'indices', i-1,i,i+1,j-1,j,j+1,k-1,k,k+1
          !!!  print*, 'suspicious', (phi(i,j+2,k)), ' ', (phi(i+1,j+2,k)), ' ', (phi(i+2,j+2,k))
          !!!  print*, 'suspicious', (phi(i,j+1,k)), ' ', (phi(i+1,j+1,k)), ' ', (phi(i+2,j+1,k))
          !!!  print*, 'suspicious', (phi(i,j,k  )), ' ', (phi(i+1,j,k  )), ' ', (phi(i+2,j,k  ))
          !!!endif
         !psimmm = seps*log((phi(i-1,j-1,k-1)+eps)/(1.-phi(i-1,j-1,k-1)+eps))  
          psimcm = seps*log((phi(i-1,j  ,k-1)+eps)/(1.-phi(i-1,j  ,k-1)+eps))  
         !psimpm = seps*log((phi(i-1,j+1,k-1)+eps)/(1.-phi(i-1,j+1,k-1)+eps))  
          psicmm = seps*log((phi(i  ,j-1,k-1)+eps)/(1.-phi(i  ,j-1,k-1)+eps))  
          psiccm = seps*log((phi(i  ,j  ,k-1)+eps)/(1.-phi(i  ,j  ,k-1)+eps))  
          psicpm = seps*log((phi(i  ,j+1,k-1)+eps)/(1.-phi(i  ,j+1,k-1)+eps))  
         !psipmm = seps*log((phi(i+1,j-1,k-1)+eps)/(1.-phi(i+1,j-1,k-1)+eps))  
          psipcm = seps*log((phi(i+1,j  ,k-1)+eps)/(1.-phi(i+1,j  ,k-1)+eps))  
         !psippm = seps*log((phi(i+1,j+1,k-1)+eps)/(1.-phi(i+1,j+1,k-1)+eps))  
          psimmc = seps*log((phi(i-1,j-1,k  )+eps)/(1.-phi(i-1,j-1,k  )+eps))  
          psimcc = seps*log((phi(i-1,j  ,k  )+eps)/(1.-phi(i-1,j  ,k  )+eps))  
          psimpc = seps*log((phi(i-1,j+1,k  )+eps)/(1.-phi(i-1,j+1,k  )+eps))  
          psicmc = seps*log((phi(i  ,j-1,k  )+eps)/(1.-phi(i  ,j-1,k  )+eps))  
          psiccc = seps*log((phi(i  ,j  ,k  )+eps)/(1.-phi(i  ,j  ,k  )+eps))  
          psicpc = seps*log((phi(i  ,j+1,k  )+eps)/(1.-phi(i  ,j+1,k  )+eps))  
          psipmc = seps*log((phi(i+1,j-1,k  )+eps)/(1.-phi(i+1,j-1,k  )+eps))  
          psipcc = seps*log((phi(i+1,j  ,k  )+eps)/(1.-phi(i+1,j  ,k  )+eps))  
          psippc = seps*log((phi(i+1,j+1,k  )+eps)/(1.-phi(i+1,j+1,k  )+eps))  
         !psimmp = seps*log((phi(i-1,j-1,k+1)+eps)/(1.-phi(i-1,j-1,k+1)+eps))  
          psimcp = seps*log((phi(i-1,j  ,k+1)+eps)/(1.-phi(i-1,j  ,k+1)+eps))  
         !psimpp = seps*log((phi(i-1,j+1,k+1)+eps)/(1.-phi(i-1,j+1,k+1)+eps))  
          psicmp = seps*log((phi(i  ,j-1,k+1)+eps)/(1.-phi(i  ,j-1,k+1)+eps))  
          psiccp = seps*log((phi(i  ,j  ,k+1)+eps)/(1.-phi(i  ,j  ,k+1)+eps))  
          psicpp = seps*log((phi(i  ,j+1,k+1)+eps)/(1.-phi(i  ,j+1,k+1)+eps))  
         !psipmp = seps*log((phi(i+1,j-1,k+1)+eps)/(1.-phi(i+1,j-1,k+1)+eps))  
          psipcp = seps*log((phi(i+1,j  ,k+1)+eps)/(1.-phi(i+1,j  ,k+1)+eps))  
         !psippp = seps*log((phi(i+1,j+1,k+1)+eps)/(1.-phi(i+1,j+1,k+1)+eps))  
          dpsidxleft   = 1.00*(psiccc       -psimcc       )*dxi
          dpsidyleft   = 0.25*(psicpc+psimpc-psicmc-psimmc)*dyi
          dpsidzleft   = 0.50*(psiccp+psimcp-psiccm-psimcm)*(dzci(k-1)*dzci(k))/(dzci(k-1)+dzci(k))
          dpsidxright  = 1.00*(psipcc       -psiccc       )*dxi
          dpsidyright  = 0.25*(psippc+psicpc-psipmc-psicmc)*dyi
          dpsidzright  = 0.50*(psipcp+psiccp-psipcm-psiccm)*(dzci(k-1)*dzci(k))/(dzci(k-1)+dzci(k))
          dpsidxback   = 0.25*(psippc+psipcc-psimpc-psimcc)*dxi
          dpsidyback   = 1.00*(psicpc       -psiccc       )*dyi
          dpsidzback   = 0.50*(psicpp+psiccp-psicpm-psiccm)*(dzci(k-1)*dzci(k))/(dzci(k-1)+dzci(k))
          dpsidxfront  = 0.25*(psipcc+psipmc-psimcc-psimmc)*dxi
          dpsidyfront  = 1.00*(psiccc       -psicmc       )*dyi
          dpsidzfront  = 0.50*(psiccp+psicmp-psiccm-psicmm)*(dzci(k-1)*dzci(k))/(dzci(k-1)+dzci(k))
          dpsidxbottom = 0.25*((psipcm-psimcm)*wghtmp+(psipcc-psimcc)*wghtmm)*dxi
          dpsidybottom = 0.25*((psicpm-psicmm)*wghtmp+(psicpc-psicmc)*wghtmm)*dyi
          dpsidzbottom = 1.00*( psiccc               -psiccm                )*dzci(k-1)
          dpsidxtop    = 0.25*((psipcp-psimcp)*wghtpm+(psipcc-psimcc)*wghtpp)*dxi
          dpsidytop    = 0.25*((psicpp-psicmp)*wghtpm+(psicpc-psicmc)*wghtpp)*dyi
          dpsidztop    = 1.00*( psiccp               -psiccc                )*dzci(k)
          nleft   = dpsidxleft  / (sqrt(  dpsidxleft**2 +   dpsidyleft**2 +   dpsidzleft**2)+eps) 
          nright  = dpsidxright / (sqrt( dpsidxright**2 +  dpsidyright**2 +  dpsidzright**2)+eps)
          nback   = dpsidyback  / (sqrt(  dpsidxback**2 +   dpsidyback**2 +   dpsidzback**2)+eps)
          nfront  = dpsidyfront / (sqrt( dpsidxfront**2 +  dpsidyfront**2 +  dpsidzfront**2)+eps)       
          nbottom = dpsidzbottom/ (sqrt(dpsidxbottom**2 + dpsidybottom**2 + dpsidzbottom**2)+eps)
          ntop    = dpsidztop   / (sqrt(   dpsidxtop**2 +    dpsidytop**2 +    dpsidztop**2)+eps)
          SHARP = 0.25*gam*((1.-(tanh(0.25*(psipcc       +psiccc       )*sepsi))**2)*nright-&
                            (1.-(tanh(0.25*(psiccc       +psimcc       )*sepsi))**2)*nleft  )*dxi + &
                  0.25*gam*((1.-(tanh(0.25*(psicpc       +psiccc       )*sepsi))**2)*nback -&
                            (1.-(tanh(0.25*(psiccc       +psicmc       )*sepsi))**2)*nfront )*dyi + &
                  0.25*gam*((1.-(tanh(0.25*(psiccp*wghtpm+psiccc*wghtpp)*sepsi))**2)*ntop  -&
                            (1.-(tanh(0.25*(psiccc*wghtmm+psiccm*wghtmp)*sepsi))**2)*nbottom)*dzfi(k)
          dphidt(i,j,k) = -ADV+DIFF-SHARP
        end do
      end do
    end do
  end subroutine pf
!
end module mod_acdi
