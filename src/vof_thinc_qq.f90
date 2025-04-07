! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 The CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_vof_thinc_qq
  !
  ! this module implements the THINC/QQ method from Bin & Xiao. J. of Comput. Phys 349 (2017): 415-440.
  !
  ! Acknowledgements: We are grateful for the discussions with Dr. Naoki Hori, which gave valuable advise
  !                   on the practical implementation of THINC/QQ for multiphase flow simulations,
  !                   with two nice repos describing the method:
  !                     - github.com/NaokiHori/SimpleVOFSolver
  !                     - github.com/NaokiHori/SimpleBubblyFlowSolver
  !
  use mpi
  use mod_types
  use mod_param, only: eps
  implicit none
  private
  public vof_thinc_transport_psi
  !
  ! constants relate to Gaussian quadrature
  !
  integer , parameter :: ngauss = 2 ! two-point Gaussian quadrature
  real(rp), dimension(ngauss) :: gauss_p = [-sqrt(3._rp)/6._rp,sqrt(3._rp)/6._rp], &
                                 gauss_w = [0.5_rp,0.5_rp]
#if 0
!  integer , parameter :: ngauss = 1 ! one-point Gaussian quadrature
!  real(rp), dimension(ngauss) :: gauss_p = [0._rp], &
!                                 gauss_w = [1._rp]
!  integer , parameter :: ngauss = 3 ! three-point Gaussian quadrature
!  real(rp), dimension(ngauss) :: gauss_p = [-sqrt(3._rp/5._rp)/2._rp,0._rp,sqrt(3._rp/5._rp)/2._rp], &
!                                 gauss_w = [5._rp/18._rp,8._rp/18._rp,5._rp/18._rp]
#endif
  !
  contains
  subroutine vof_thinc_transport_psi(n,dli,dzfi,beta,u,v,w,normx,normy,normz,psi,dpsidt,flux_x,flux_y,flux_z)
    !
    ! compute right-hand side of the VoF transport equation
    !
    ! uses mid-point rule (1-point Gaussian quadrature), for which there is an analytical solution
    ! for regularized Heaviside function intercept
    !
    implicit none
    integer , intent(in   ), dimension(3)        :: n
    real(rp), intent(in   ), dimension(3)        :: dli
    real(rp), intent(in   ), dimension(0:)       :: dzfi
    real(rp), intent(in   )                      :: beta
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(in   ), dimension(0:,0:,0:) :: normx,normy,normz
    real(rp), intent(inout), dimension(0:,0:,0:) :: psi
    real(rp), intent(out  ), dimension(: ,: ,: ) :: dpsidt
    real(rp), intent(out  ), dimension(0:,0:,0:) :: flux_x,flux_y,flux_z
    integer  :: i,j,k,ii,jj,kk,q1,q2
    real(rp) :: xx,yy,zz
    real(rp) :: vel_x,vel_y,vel_z,lpsi_x,lpsi_y,lpsi_z, &
                lnormx_x,lnormx_y,lnormx_z, &
                lnormy_x,lnormy_y,lnormy_z, &
                lnormz_x,lnormz_y,lnormz_z, &
                d_x,d_y,d_z,lflux_x,lflux_y,lflux_z
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=0,n(3)
      do j=0,n(2)
        do i=0,n(1)
          vel_x = u(i,j,k)
          vel_y = v(i,j,k)
          vel_z = w(i,j,k)
          !
          if(vel_x > 0.) then
            ii = i
            xx = 0.5
          else
            ii = i+1
            xx = -0.5
          end if
          if(vel_y > 0.) then
            jj = j
            yy = 0.5
          else
            jj = j+1
            yy = -0.5
          end if
          if(vel_z > 0.) then
            kk = k
            zz = 0.5
          else
            kk = k+1
            zz = -0.5
          end if
          !
          lpsi_x = psi(ii,j,k)
          lpsi_y = psi(i,jj,k)
          lpsi_z = psi(i,j,kk)
          lnormx_x = normx(ii,j,k)
          lnormy_y = normy(i,jj,k)
          lnormz_z = normz(i,j,kk)
          !
          d_x = -.5_rp/beta*log(1._rp/(lpsi_x+eps)-1._rp+eps) !d_x = atanh(1._rp/lpsi_x-1._rp)
          d_y = -.5_rp/beta*log(1._rp/(lpsi_y+eps)-1._rp+eps) !d_y = atanh(1._rp/lpsi_y-1._rp)
          d_z = -.5_rp/beta*log(1._rp/(lpsi_z+eps)-1._rp+eps) !d_z = atanh(1._rp/lpsi_z-1._rp)
          !
          lflux_x = 1._rp/(1._rp + exp(-2._rp*beta*(lnormx_x*xx + d_x)))
          lflux_y = 1._rp/(1._rp + exp(-2._rp*beta*(lnormy_y*yy + d_y)))
          lflux_z = 1._rp/(1._rp + exp(-2._rp*beta*(lnormz_z*zz + d_z)))
          flux_x(i,j,k) = vel_x*lflux_x
          flux_y(i,j,k) = vel_y*lflux_y
          flux_z(i,j,k) = vel_z*lflux_z
        end do
      end do
    end do
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          dpsidt(i,j,k) = -(dli(1)*( flux_x(i,j,k)-flux_x(i-1,j,k)) + &
                            dli(2)*( flux_y(i,j,k)-flux_y(i,j-1,k)) + &
                            dzfi(k)*(flux_z(i,j,k)-flux_z(i,j,k-1)))
        end do
      end do
    end do
  end subroutine vof_thinc_transport_psi
  !
  subroutine vof_thinc_transport_psi_gauss(n,dli,dzfi,vofmin,beta,u,v,w,normx,normy,normz,psi,dpsidt,flux_x,flux_y,flux_z)
    !
    ! compute right-hand side of the VoF transport equation
    !
    implicit none
    integer , intent(in   ), dimension(3)        :: n
    real(rp), intent(in   ), dimension(3)        :: dli
    real(rp), intent(in   ), dimension(0:)       :: dzfi
    real(rp), intent(in   )                      :: vofmin,beta
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(in   ), dimension(0:,0:,0:) :: normx,normy,normz
    real(rp), intent(inout), dimension(0:,0:,0:) :: psi
    real(rp), intent(out  ), dimension(: ,: ,: ) :: dpsidt
    real(rp), intent(out  ), dimension(0:,0:,0:) :: flux_x,flux_y,flux_z
    integer  :: i,j,k,ii,jj,kk,q1,q2
    real(rp) :: xx,yy,zz,r1,r2,ww
    real(rp) :: vel_x,vel_y,vel_z,lpsi_x,lpsi_y,lpsi_z, &
                lnormx_x,lnormx_y,lnormx_z, &
                lnormy_x,lnormy_y,lnormy_z, &
                lnormz_x,lnormz_y,lnormz_z, &
                d_x,d_y,d_z,lflux_x,lflux_y,lflux_z
    !
    do k=0,n(3)
      do j=0,n(2)
        do i=0,n(1)
          vel_x = u(i,j,k)
          vel_y = v(i,j,k)
          vel_z = w(i,j,k)
          !
          if(vel_x > 0.) then
            ii = i
            xx = 0.5
          else
            ii = i+1
            xx = -0.5
          end if
          if(vel_y > 0.) then
            jj = j
            yy = 0.5
          else
            jj = j+1
            yy = -0.5
          end if
          if(vel_z > 0.) then
            kk = k
            zz = 0.5
          else
            kk = k+1
            zz = -0.5
          end if
          !
          lpsi_x = psi(ii,j,k)
          lpsi_y = psi(i,jj,k)
          lpsi_z = psi(i,j,kk)
          lnormx_x = normx(ii,j,k)
          lnormx_y = normy(ii,j,k)
          lnormx_z = normz(ii,j,k)
          lnormy_x = normx(i,jj,k)
          lnormy_y = normy(i,jj,k)
          lnormy_z = normz(i,jj,k)
          lnormz_x = normx(i,j,kk)
          lnormz_y = normy(i,j,kk)
          lnormz_z = normz(i,j,kk)
          !
          d_x = cmpt_intercept(lpsi_x,beta,lnormx_x,lnormy_x,lnormz_x)
          d_y = cmpt_intercept(lpsi_y,beta,lnormx_y,lnormy_y,lnormz_y)
          d_z = cmpt_intercept(lpsi_z,beta,lnormx_z,lnormy_z,lnormz_z)
          !
          lflux_x = 0.
          lflux_y = 0.
          lflux_z = 0.
          do q2=1,ngauss
            do q1=1,ngauss
              ww = gauss_w(q1)*gauss_w(q2)
              r1 = gauss_p(q1)
              r2 = gauss_p(q2)
              lflux_x = lflux_x + ww*smooth_step_tanh(beta,lnormx_x,lnormy_x,lnormz_x,xx,r1,r2,d_x)
              lflux_y = lflux_y + ww*smooth_step_tanh(beta,lnormx_y,lnormy_y,lnormz_y,r1,yy,r2,d_y)
              lflux_z = lflux_z + ww*smooth_step_tanh(beta,lnormx_z,lnormy_z,lnormz_z,r1,r2,zz,d_z)
            end do
          end do
          if(lpsi_x < vofmin .or. lpsi_x > (1. - vofmin)) then
            lflux_x = lpsi_x
            !cycle
          end if
          if(lpsi_y < vofmin .or. lpsi_y > (1. - vofmin)) then
            lflux_y = lpsi_y
            !cycle
          end if
          if(lpsi_z < vofmin .or. lpsi_z > (1. - vofmin)) then
            lflux_z = lpsi_z
            !cycle
          end if
          flux_x(i,j,k) = vel_x*lflux_x
          flux_y(i,j,k) = vel_y*lflux_y
          flux_z(i,j,k) = vel_z*lflux_z
        end do
      end do
    end do
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          dpsidt(i,j,k) = -(dli(1)*( flux_x(i,j,k)-flux_x(i-1,j,k)) + &
                            dli(2)*( flux_y(i,j,k)-flux_y(i,j-1,k)) + &
                            dzfi(k)*(flux_z(i,j,k)-flux_z(i,j,k-1)))
        end do
      end do
    end do
  end subroutine vof_thinc_transport_psi_gauss
  !
  pure function cmpt_intercept(psi,beta,normx,normy,normz) result(d)
    !
    ! Newton-Rapshon iteration to compute the d parameter in the THINC/QQ reconstruction function,
    ! intercepting the numerical integration of the multi-dimensional hyperbolic tangent function
    ! with the cell volume fraction
    !
    ! solves for \sum_{g=1}^{NG} \omega_{g} (A_g + D)/(1 + A_g D) = Q,
    ! with A_g = \tanh(\beta P(\xi_{g})); D = \tanh(\beta d); Q = 2(\psi - 1/2)
    ! with some slight changes using the identity (1+tanh(x))/2 == 1/(1+exp(-2x))
    !
    implicit none
    integer , parameter :: itermax = 8
    real(rp), parameter :: resmax = 1.e-12
    real(rp), intent(in) :: psi,beta
    real(rp), intent(in) :: normx,normy,normz
    real(rp) :: d
    !
    real(rp), dimension(ngauss,ngauss,ngauss) :: betap
    integer  :: i,j,k,iter
    real(rp) :: p,w,factor
    real(rp) :: f0,f1
    real(rp) :: res,sol
    !$acc routine seq
    !
    ! pre-compute constant terms tanh(\beta P(\xi_{g}))/2
    !
    do k=1,ngauss
      do j=1,ngauss
        do i=1,ngauss
          betap(i,j,k) = exp(-2._rp*beta*(normx*gauss_p(i) + normy*gauss_p(j) + normz*gauss_p(k)))
          !betap(i,j,k) = tanh(beta*(normx*gauss_p(i) + normy*gauss_p(j) + normz*gauss_p(k)))
        end do
      end do
    end do
    !
    ! initial guess based on the solution with the first-order Guassian quadrature
    ! careful about division by zero - should be avoided in other subroutines
    !
    sol = 1./psi - 1. ! careful about division by zero
    iter = 0
    res = huge(1._rp)
    do while(abs(res) > resmax .or. iter < itermax)
      f0 = -psi
      f1 = 0.
      do k=1,ngauss
        do j=1,ngauss
          do i=1,ngauss
            w = gauss_w(i)*gauss_w(j)*gauss_w(k)
            p = betap(i,j,k)
            factor = 1._rp/(1._rp+p*sol) ! careful about division by zero
            f0 = f0 + w*factor
            f1 = f1 - w*p*factor**2
          end do
        end do
      end do
      sol = sol - f0/f1 ! careful about division by zero
      res = f0
      iter = iter + 1
    end do
    !
    ! compute d from D
    !
    d = -0.5_rp*log(sol)/beta
    !d = atanh(sol)/beta
  end function cmpt_intercept
  !
  pure function smooth_step_tanh(beta,nx,ny,nz,x,y,z,d) result(h)
    !
    ! THINC indicator function
    !
    real(rp), intent(in) :: beta,nx,ny,nz,x,y,z,d
    real(rp) :: h
    !$acc routine seq
    h = 1._rp/(1._rp + exp(-2._rp*beta*(nx*x + ny*y + nz*z + d)))
    !h = 0.5_rp*(1._rp+tanh(beta*(nx*x+ny*y+nz*z+d)))
  end function smooth_step_tanh
end module mod_vof_thinc_qq
