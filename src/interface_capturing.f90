! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_interface_capturing
  !
  ! common interface to the interface-capturing methods, selected at run time (see `mod_param`):
  !   - ACDI diffuse-interface method (`ifc_reconstruction = 'acdi'`, default)
  !   - THINC/QQ volume-of-fluid method (`'thinc'`)
  !   - geometric volume-of-fluid method with conservative direction-split advection (`'plic'`)
  !
  ! the main program calls, in this order:
  !   ifc_init_params   -> method parameters needed by the initial condition
  !   ifc_init_field    -> initial phase indicator
  !   ifc_init          -> initialization from the initial (or loaded) phase indicator
  !   ifc_advance       -> one time step of the phase indicator, returning the (volumetric) phase fluxes
  !                        used for the mass-consistent transport of momentum and scalars, and, for the
  !                        geometric volume-of-fluid method, the advecting velocity of the time step
  !                        (ifc_advecting_velocity returns the latter when the interface is not tracked)
  !   ifc_update_params -> update of the method parameters after each time step
  !   ifc_write_log     -> log file of the method
  !
  use mod_types
  use mod_common_mpi, only: myid
  use mod_param     , only: datadir,ng,nb,is_bound,cbcvel,bcvel,cbcpsi,bcpsi,cbcnor,bcnor,cbccur,bccur, &
                            acdi_gam_factor,acdi_gam_min,vof_thinc_beta, &
                            rho12,cp12,nscal,cbcsca,bcsca, &
                            is_ifc_acdi,is_ifc_plic,is_ifc_split,is_sdf_normals,is_vof_sweep_transport
  use mod_bound     , only: boundp,bounduvw
  use mod_output    , only: out0d
  use mod_two_fluid , only: init2fl,init2fl_sharp,cmpt_norm_curv => cmpt_norm_curv_youngs
  use mod_rk        , only: rk_2fl
  use mod_acdi      , only: acdi_set_epsilon
  use mod_acdi      , only: acdi_set_gamma,acdi_cmpt_phi
  use mod_vof_split , only: vof_split_init,vof_split_advance,vof_split_velocity,vof_split_curvature, &
                            vof_split_diagnostics
  use mod_vof_thinc_qq, only: vof_thinc_cmpt_phi
  implicit none
  private
  public ifc_init_params,ifc_init_field,ifc_init,ifc_advance,ifc_advecting_velocity,ifc_update_params,ifc_write_log
  real(rp) :: vol0 ! initial volume of phase 1 (geometric volume-of-fluid method)
  contains
  subroutine ifc_init_params(dl,nh_g,dzfi,psi_thickness_factor,seps)
    !
    ! method parameters needed by the initial condition: the interface width, used by all methods
    !
    implicit none
    real(rp), intent(in   ), dimension(3)       :: dl
    integer , intent(in   )                     :: nh_g
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzfi
    real(rp), intent(in   ) :: psi_thickness_factor
    real(rp), intent(inout) :: seps
    call acdi_set_epsilon(dl,nh_g,dzfi,psi_thickness_factor,seps)
  end subroutine ifc_init_params
  !
  subroutine ifc_init_field(inipsi,seps,lo,hi,nh_psi,nh_g,l,dl,zc_g,dzf_g,psi)
    !
    ! initial condition for the phase indicator: a smooth profile of width `seps`, or
    ! the exact volume fractions for the geometric volume-of-fluid method
    !
    implicit none
    character(len=*), intent(in) :: inipsi
    real(rp), intent(in ) :: seps
    integer , intent(in ), dimension(3)        :: lo,hi
    integer , intent(in )                      :: nh_psi,nh_g
    real(rp), intent(in ), dimension(3)        :: l,dl
    real(rp), intent(in ), dimension(1-nh_g:)  :: zc_g,dzf_g
    real(rp), intent(out), dimension(lo(1)-nh_psi:,lo(2)-nh_psi:,lo(3)-nh_psi:) :: psi
    if(.not.is_ifc_plic) then
      call init2fl(inipsi,cbcpsi,seps,lo,hi,nh_psi,nh_g,l,dl,zc_g,psi)
    else
      call init2fl_sharp(inipsi,cbcpsi,lo,hi,nh_psi,nh_g,l,dl,zc_g,dzf_g,psi)
    end if
  end subroutine ifc_init_field
  !
  subroutine ifc_init(n,nh_v,nh_psi,nh_ifc,nh_g,dl,dli,dzc,dzf,dzci,dzfi,seps,u,v,w,psi,phi,normx,normy,normz,kappa,gam)
    !
    ! initialization from the initial (or loaded) phase indicator, already present on the device
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v,nh_psi,nh_ifc,nh_g
    real(rp), intent(in   ), dimension(3) :: dl,dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzc,dzf,dzci,dzfi
    real(rp), intent(in   ) :: seps
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)       :: u,v,w
    real(rp), intent(inout), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: normx,normy,normz,kappa
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:), optional :: phi ! (present with `is_sdf_normals`)
    real(rp), intent(inout) :: gam
    real(rp) :: clip,fmin,fmax,cfl(3),kstats(3)
    integer  :: nsub,i,j,k
    !
    call boundp(cbcpsi,n,nh_psi,nh_g,bcpsi,nb,is_bound,dl,dzc,psi)
    !$acc wait
    if(is_ifc_split) then
      call vof_split_init(ng,n,nh_psi,is_vof_sweep_transport,nscal > 0)
      call vof_split_diagnostics(n,nh_psi,nh_g,dl,dzf,psi,vol0,clip,fmin,fmax,cfl,nsub,kstats)
      if(myid == 0) print*, 'Geometric VoF. Initial volume of phase 1: ', vol0
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1-nh_ifc,n(3)+nh_ifc
        do j=1-nh_ifc,n(2)+nh_ifc
          do i=1-nh_ifc,n(1)+nh_ifc
            normx(i,j,k) = 0.
            normy(i,j,k) = 0.
            normz(i,j,k) = 0.
            kappa(i,j,k) = 0.
          end do
        end do
      end do
    end if
    call ifc_geometry(n,nh_psi,nh_ifc,nh_g,dl,dli,dzc,dzf,dzci,dzfi,seps,psi,phi,normx,normy,normz,kappa)
    if(is_ifc_acdi) then
      call acdi_set_gamma(n,nh_v,acdi_gam_factor,u,v,w,gam)
      gam = max(gam,acdi_gam_min)
      if(myid == 0) print*, 'ACDI parameters. Gamma: ', gam, 'Epsilon: ', seps
    end if
  end subroutine ifc_init
  !
  subroutine ifc_advance(tm_coeff,n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,dl,dli,dzc,dzf,dzci,dzfi,wzc,dt,dto,gam,seps,u,v,w, &
                         normx,normy,normz,phi,kappa,psi,psiflx_x,psiflx_y,psiflx_z,is_first,dpsidtrko, &
                         istep,uo,vo,wo,uadv,vadv,wadv,s,sadv)
    !
    ! advances the phase indicator over one time step, returns the phase fluxes (with their boundary
    ! conditions), and updates the interface normals and curvature; with the geometric volume-of-fluid
    ! method and `is_vof_sweep_transport`, the velocity (and the scalar) is transported along with the
    ! phase indicator, and returned in uadv, vadv, wadv (sadv) instead of the advecting velocity
    !
    implicit none
    real(rp), intent(in   ), dimension(2) :: tm_coeff
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v,nh_s,nh_psi,nh_ifc,nh_g
    real(rp), intent(in   ), dimension(3) :: dl,dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzc,dzf,dzci,dzfi,wzc
    real(rp), intent(in   ) :: dt,dto,gam,seps
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)       :: u,v,w
    real(rp), intent(inout), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: normx,normy,normz,kappa
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:), optional :: phi ! (present with `is_sdf_normals`)
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: psiflx_x,psiflx_y,psiflx_z
    logical , intent(in   ) :: is_first
    real(rp), intent(inout), dimension(:,:,:) :: dpsidtrko
    integer , intent(in   ) :: istep
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:), optional :: uo,vo,wo,uadv,vadv,wadv
    real(rp), intent(in   ), dimension(1-nh_s:,1-nh_s:,1-nh_s:), optional :: s
    real(rp), intent(inout), dimension(1-nh_s:,1-nh_s:,1-nh_s:), optional :: sadv
    !
    if(.not.is_ifc_split) then
      call rk_2fl(tm_coeff,n,nh_v,nh_psi,nh_ifc,nh_g,dli,dzci,dzfi,dt,gam,seps,vof_thinc_beta,u,v,w,normx,normy,normz,phi,psi, &
                  psiflx_x,psiflx_y,psiflx_z,is_first,dpsidtrko)
      call bounduvw(cbcvel,n,nh_ifc,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,psiflx_x,psiflx_y,psiflx_z)
      call boundp(cbcpsi,n,nh_psi,nh_g,bcpsi,nb,is_bound,dl,dzc,psi)
    else
      call vof_split_advance(tm_coeff,is_first,.true.,istep,n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,dl,dzc,dzf,dzci,dzfi,wzc,dt,dto, &
                             cbcpsi,bcpsi,nb,is_bound,rho12,rho12*cp12,u,v,w,psi,psiflx_x,psiflx_y,psiflx_z, &
                             uo,vo,wo,uadv,vadv,wadv,s,sadv)
      if(.not.is_vof_sweep_transport) then
        call bounduvw(cbcvel,n,nh_ifc,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,psiflx_x,psiflx_y,psiflx_z)
      else
        call bounduvw(cbcvel,n,nh_v,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,uadv,vadv,wadv)
        if(nscal > 0) call boundp(cbcsca,n,nh_s,nh_g,bcsca,nb,is_bound,dl,dzc,sadv)
      end if
    end if
    call ifc_geometry(n,nh_psi,nh_ifc,nh_g,dl,dli,dzc,dzf,dzci,dzfi,seps,psi,phi,normx,normy,normz,kappa)
  end subroutine ifc_advance
  !
  subroutine ifc_advecting_velocity(tm_coeff,is_first,istep,n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,dl,dzc,dzf,dzci,dzfi,wzc,dt,dto, &
                                    u,v,w,psi,psiflx_x,psiflx_y,psiflx_z,uo,vo,wo,uadv,vadv,wadv,s,sadv)
    !
    ! advecting velocity of a time step in which the interface is not tracked (geometric volume-of-fluid
    ! method only); with `is_vof_sweep_transport`, the velocity (and the scalar) is transported with the
    ! sweeps, without updating the phase indicator, and returned in uadv, vadv, wadv (sadv) instead
    !
    implicit none
    real(rp), intent(in   ), dimension(2) :: tm_coeff
    logical , intent(in   )               :: is_first
    integer , intent(in   )               :: istep
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v,nh_s,nh_psi,nh_ifc,nh_g
    real(rp), intent(in   ), dimension(3) :: dl
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzc,dzf,dzci,dzfi,wzc
    real(rp), intent(in   ) :: dt,dto
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(inout), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: psiflx_x,psiflx_y,psiflx_z
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:), optional :: uo,vo,wo,uadv,vadv,wadv
    real(rp), intent(in   ), dimension(1-nh_s:,1-nh_s:,1-nh_s:), optional :: s
    real(rp), intent(inout), dimension(1-nh_s:,1-nh_s:,1-nh_s:), optional :: sadv
    real(rp), dimension(3) :: cfl
    !
    if(.not.is_ifc_split) return
    if(is_vof_sweep_transport) then
      call vof_split_advance(tm_coeff,is_first,.false.,istep,n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,dl,dzc,dzf,dzci,dzfi,wzc,dt,dto, &
                             cbcpsi,bcpsi,nb,is_bound,rho12,rho12*cp12,u,v,w,psi,psiflx_x,psiflx_y,psiflx_z, &
                             uo,vo,wo,uadv,vadv,wadv,s,sadv)
      call bounduvw(cbcvel,n,nh_v,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,uadv,vadv,wadv)
      if(nscal > 0) call boundp(cbcsca,n,nh_s,nh_g,bcsca,nb,is_bound,dl,dzc,sadv)
    else
      call vof_split_velocity(tm_coeff,is_first,n,nh_v,nh_g,dl,dzfi,dt,u,v,w,uo,vo,wo,uadv,vadv,wadv,cfl)
    end if
  end subroutine ifc_advecting_velocity
  !
  subroutine ifc_geometry(n,nh_psi,nh_ifc,nh_g,dl,dli,dzc,dzf,dzci,dzfi,seps,psi,phi,normx,normy,normz,kappa)
    !
    ! interface normals and curvature, with their boundary conditions
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_psi,nh_ifc,nh_g
    real(rp), intent(in   ), dimension(3) :: dl,dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzc,dzf,dzci,dzfi
    real(rp), intent(in   ) :: seps
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: normx,normy,normz,kappa
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:), optional :: phi ! (present with `is_sdf_normals`)
    character(len=1), dimension(0:1,3) :: cbc_tmp
    !
    if(is_ifc_plic) then
      !
      ! height-function curvature, in stages separated by halo exchanges; the intermediate exchanges mirror
      ! (instead of imposing Dirichlet conditions), such that undefined values are only copied; the interface
      ! normals are not needed with this method, and are left at zero
      !
      cbc_tmp(:,:) = cbccur(:,:)
      where(cbc_tmp(:,:) == 'D') cbc_tmp(:,:) = 'N'
      call vof_split_curvature(1,n,nh_psi,nh_ifc,nh_g,dl,dzf,psi,kappa)
      call boundp(cbc_tmp,n,nh_ifc,nh_g,bccur,nb,is_bound,dl,dzc,kappa)
      call vof_split_curvature(2,n,nh_psi,nh_ifc,nh_g,dl,dzf,psi,kappa)
      call boundp(cbc_tmp,n,nh_ifc,nh_g,bccur,nb,is_bound,dl,dzc,kappa)
      call vof_split_curvature(3,n,nh_psi,nh_ifc,nh_g,dl,dzf,psi,kappa)
      call boundp(cbccur,n,nh_ifc,nh_g,bccur,nb,is_bound,dl,dzc,kappa)
    else
      if(is_sdf_normals) then
        if(is_ifc_acdi) then
          call acdi_cmpt_phi(n,nh_psi,nh_ifc,seps,psi,phi)
        else
          call vof_thinc_cmpt_phi(n,nh_psi,nh_ifc,vof_thinc_beta,psi,phi)
        end if
        call cmpt_norm_curv(n,nh_psi,nh_ifc,nh_g,dli,dzci,dzfi,phi,normx,normy,normz,kappa)
      else
        call cmpt_norm_curv(n,nh_psi,nh_ifc,nh_g,dli,dzci,dzfi,psi,normx,normy,normz,kappa)
      end if
      call boundp(cbccur,n,nh_ifc,nh_g,bccur,nb,is_bound,dl,dzc,kappa)
      call boundp(cbcnor(:,:,1),n,nh_ifc,nh_g,bcnor(:,:,1),nb,is_bound,dl,dzc,normx)
      call boundp(cbcnor(:,:,2),n,nh_ifc,nh_g,bcnor(:,:,2),nb,is_bound,dl,dzc,normy)
      call boundp(cbcnor(:,:,3),n,nh_ifc,nh_g,bcnor(:,:,3),nb,is_bound,dl,dzc,normz)
    end if
  end subroutine ifc_geometry
  !
  subroutine ifc_update_params(n,nh_v,u,v,w,gam)
    !
    ! update of the method parameters after a time step
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(inout) :: gam
    if(is_ifc_acdi) then
      call acdi_set_gamma(n,nh_v,acdi_gam_factor,u,v,w,gam)
      gam = max(gam,acdi_gam_min)
    end if
  end subroutine ifc_update_params
  !
  subroutine ifc_write_log(istep,time,n,nh_psi,nh_g,dl,dzf,gam,seps,psi)
    !
    ! log file of the method
    !
    implicit none
    integer , intent(in) :: istep
    real(rp), intent(in) :: time
    integer , intent(in), dimension(3) :: n
    integer , intent(in)               :: nh_psi,nh_g
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(1-nh_g:) :: dzf
    real(rp), intent(in) :: gam,seps
    real(rp), intent(in), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), dimension(14) :: var
    real(rp) :: vol,clip,fmin,fmax,cfl(3),kstats(3)
    integer  :: nsub
    !
    if(is_ifc_split) then
      !
      ! volume of phase 1 and its relative change, accumulated volume added by clipping, extrema of the
      ! volume fraction before clipping, CFL numbers and number of sub-steps of the last time step, and
      ! numbers of interfacial cells with height-function curvature, neighbour average, and no estimate
      !
      call vof_split_diagnostics(n,nh_psi,nh_g,dl,dzf,psi,vol,clip,fmin,fmax,cfl,nsub,kstats)
      var(1)    = 1.*istep
      var(2)    = time
      var(3)    = vol
      var(4)    = (vol-vol0)/max(vol0,tiny(1._rp))
      var(5)    = clip
      var(6)    = fmin
      var(7)    = fmax
      var(8:10) = cfl(:)
      var(11)   = 1.*nsub
      var(12:14) = kstats(:)
      call out0d(trim(datadir)//'log_vof.out',14,var)
    else if(is_ifc_acdi) then
      var(1) = 1.*istep
      var(2) = time
      var(3) = gam
      var(4) = seps
      call out0d(trim(datadir)//'log_acdi.out',4,var(1:4))
    end if
  end subroutine ifc_write_log
end module mod_interface_capturing
