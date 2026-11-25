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
  !                        used for the mass-consistent transport of momentum and scalars, and, with the
  !                        direction-split advection, the advecting velocity of the time step
  !                        (ifc_advecting_velocity returns the latter when the interface is not tracked)
  !   ifc_update_params -> update of the method parameters after each time step
  !   ifc_write_log     -> log file of the method
  !
  use mpi
  use mod_types
  use mod_common_mpi, only: myid,ierr
  use mod_param     , only: datadir,ng,nb,is_bound,cbcvel,bcvel,cbcvel_halo,bcvel_halo,cbcpsi,bcpsi,cbcnor,bcnor, &
                            cbccur,bccur,contact_angle, &
                            acdi_gam_factor,acdi_gam_min,vof_thinc_beta, &
                            rho12,cp12,nscal,cbcscal,bcscal, &
                            is_ifc_acdi,is_ifc_thinc,is_ifc_plic,is_ifc_split,is_sdf_normals,is_vof_sweep_transport, &
                            is_hf_normals,sigma,inipsi
  use mod_bound     , only: boundp,bounduvw
  use mod_mom       , only: cmpt_capillary_force
  use mod_output    , only: out0d
  use mod_two_fluid , only: init2fl,init2fl_sharp,cmpt_norm_curv => cmpt_norm_curv_youngs
  use mod_rk        , only: rk_2fl
  use mod_acdi      , only: acdi_set_epsilon,acdi_set_gamma,acdi_cmpt_phi
  use mod_scal      , only: scalar
  use mod_vof_split , only: vof_split_init,vof_split_advance,vof_split_velocity,vof_split_curvature,vof_split_contact_angle, &
                            vof_split_diagnostics,vof_split_curvature_reset,vof_marker,markers_sum
  use mod_thinc     , only: thinc_cmpt_phi
  implicit none
  private
  public ifc_init_params,ifc_init_field,ifc_init,ifc_advance,ifc_advecting_velocity,ifc_update_params,ifc_write_log, &
         marker,ifc_markers_init,ifc_markers_advance,ifc_markers_force
  !
  ! several markers: the dispersed phase (drops or bubbles) is distributed over phase indicators that are
  ! advanced separately, each as the phase indicator of a run with its drops only, such that interfaces of
  ! different markers do not merge; the fluid properties follow the total phase indicator, and the surface
  ! tension forces of the markers are summed
  !
  type marker
    real(rp), allocatable, dimension(:,:,:) :: psi,dpsidtrko,kappa,normx,normy,normz,phi
  end type marker
  real(rp), allocatable, dimension(:,:,:) :: wflx_x,wflx_y,wflx_z ! phase fluxes of one marker
  real(rp) :: vol0 ! initial volume of phase 1 (direction-split advection)
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
  subroutine ifc_init_field(inipsi,seps,lo,hi,nh_psi,nh_g,l,dl,zc_g,dzf_g,psi,imarker)
    !
    ! initial condition for the phase indicator: a smooth profile of width `seps`, or
    ! the exact volume fractions for the geometric volume-of-fluid method; with `imarker`, that of the spheres
    ! of that marker alone
    !
    implicit none
    character(len=*), intent(in ) :: inipsi
    real(rp), intent(in ) :: seps
    integer , intent(in ), dimension(3)        :: lo,hi
    integer , intent(in )                      :: nh_psi,nh_g
    real(rp), intent(in ), dimension(3)        :: l,dl
    real(rp), intent(in ), dimension(1-nh_g:)  :: zc_g,dzf_g
    real(rp), intent(out), dimension(lo(1)-nh_psi:,lo(2)-nh_psi:,lo(3)-nh_psi:) :: psi
    integer , intent(in ), optional :: imarker
    if(.not.is_ifc_plic) then
      call init2fl(inipsi,cbcpsi,seps,lo,hi,nh_psi,nh_g,l,dl,zc_g,psi,imarker)
    else
      call init2fl_sharp(inipsi,cbcpsi,lo,hi,nh_psi,nh_g,l,dl,zc_g,dzf_g,psi,imarker)
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
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: kappa
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:), optional :: normx,normy,normz ! (absent with `plic`)
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:), optional :: phi ! (present with `is_sdf_normals`)
    real(rp), intent(inout) :: gam
    real(rp) :: clip,fmin,fmax,cfl(3),kstats(3)
    integer  :: nsub,i,j,k
    !
    call boundp(cbcpsi,n,nh_psi,nh_g,bcpsi,nb,is_bound,dl,dzc,psi)
    !$acc wait
    if(is_ifc_split) then
      call vof_split_init(ng,n,nh_psi,is_ifc_thinc,vof_thinc_beta,is_hf_normals,is_vof_sweep_transport,nscal,contact_angle)
      call vof_split_contact_angle(n,nh_psi,nh_g,dl,dzc,dzf,cbcpsi,bcpsi,nb,is_bound,psi)
      call vof_split_diagnostics(n,nh_psi,nh_g,dl,dzf,psi,vol0,clip,fmin,fmax,cfl,nsub,kstats)
      if(myid == 0) print*, 'Direction-split advection. Initial volume of phase 1: ', vol0
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1-nh_ifc,n(3)+nh_ifc
        do j=1-nh_ifc,n(2)+nh_ifc
          do i=1-nh_ifc,n(1)+nh_ifc
            kappa(i,j,k) = 0.
          end do
        end do
      end do
      if(.not.is_ifc_plic) then
        !$acc parallel loop collapse(3) default(present) async(1)
        do k=1-nh_ifc,n(3)+nh_ifc
          do j=1-nh_ifc,n(2)+nh_ifc
            do i=1-nh_ifc,n(1)+nh_ifc
              normx(i,j,k) = 0.
              normy(i,j,k) = 0.
              normz(i,j,k) = 0.
            end do
          end do
        end do
      end if
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
                         istep,uo,vo,wo,uadv,vadv,wadv,scalars)
    !
    ! advances the phase indicator over one time step, returns the phase fluxes (with their boundary
    ! conditions), and updates the interface normals and curvature; with the direction-split
    ! advection and `is_vof_sweep_transport`, the velocity (and the scalars) is transported along with the
    ! phase indicator, and returned in uadv, vadv, wadv (and in `sadv` of each scalar) instead of the advecting
    ! velocity
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
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: kappa
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:), optional :: normx,normy,normz ! (absent with `plic`)
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:), optional :: phi ! (present with `is_sdf_normals`)
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: psiflx_x,psiflx_y,psiflx_z
    logical , intent(in   ) :: is_first
    real(rp), intent(inout), dimension(:,:,:), optional :: dpsidtrko ! (present with the AB2 transport)
    integer , intent(in   ) :: istep
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:), optional :: uo,vo,wo,uadv,vadv,wadv
    type(scalar), intent(inout), dimension(:), optional :: scalars
    !
    if(.not.is_ifc_split) then
      call rk_2fl(tm_coeff,n,nh_v,nh_psi,nh_ifc,nh_g,dli,dzci,dzfi,dt,gam,seps,vof_thinc_beta,u,v,w,normx,normy,normz,phi,psi, &
                  psiflx_x,psiflx_y,psiflx_z,is_first,dpsidtrko)
      call bounduvw(cbcvel,n,nh_ifc,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,psiflx_x,psiflx_y,psiflx_z)
      call boundp(cbcpsi,n,nh_psi,nh_g,bcpsi,nb,is_bound,dl,dzc,psi)
    else
      call advance_split(tm_coeff,n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,dl,dzc,dzf,dzci,dzfi,wzc,dt,dto,u,v,w,psi, &
                         psiflx_x,psiflx_y,psiflx_z,is_first,istep,uo,vo,wo,uadv,vadv,wadv,scalars)
    end if
    call ifc_geometry(n,nh_psi,nh_ifc,nh_g,dl,dli,dzc,dzf,dzci,dzfi,seps,psi,phi,normx,normy,normz,kappa)
  end subroutine ifc_advance
  !
  subroutine advance_split(tm_coeff,n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,dl,dzc,dzf,dzci,dzfi,wzc,dt,dto,u,v,w,psi, &
                           psiflx_x,psiflx_y,psiflx_z,is_first,istep,uo,vo,wo,uadv,vadv,wadv,scalars,psi_m,nm_off)
    !
    ! direction-split advection of the phase indicator, or of those of several markers (`psi_m`, with the
    ! offset `nm_off` of their total; see `vof_split_advance`), followed by the boundary conditions of the
    ! phase fluxes or, with the transport along with the sweeps, of the transported velocity and scalars
    !
    implicit none
    real(rp), intent(in   ), dimension(2) :: tm_coeff
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v,nh_s,nh_psi,nh_ifc,nh_g
    real(rp), intent(in   ), dimension(3) :: dl
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzc,dzf,dzci,dzfi,wzc
    real(rp), intent(in   ) :: dt,dto
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)       :: u,v,w
    real(rp), intent(inout), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: psiflx_x,psiflx_y,psiflx_z
    logical , intent(in   ) :: is_first
    integer , intent(in   ) :: istep
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: uo,vo,wo,uadv,vadv,wadv
    type(scalar), intent(inout), dimension(:), optional :: scalars
    type(vof_marker), intent(in   ), dimension(:), optional :: psi_m
    integer , intent(in   ), optional :: nm_off
    real(rp), dimension(2,max(nscal,1)) :: rhocp12
    integer :: is
    !
    do is=1,max(nscal,1)
      rhocp12(:,is) = rho12(:)*cp12(:,is)
    end do
    call vof_split_advance(tm_coeff,is_first,.true.,istep,n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,dl,dzc,dzf,dzci,dzfi,wzc,dt,dto, &
                           cbcpsi,bcpsi,nb,is_bound,rho12,rhocp12,u,v,w,psi,psiflx_x,psiflx_y,psiflx_z, &
                           uo,vo,wo,uadv,vadv,wadv,scalars,psi_m,nm_off)
    if(.not.is_vof_sweep_transport) then
      call bounduvw(cbcvel,n,nh_ifc,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,psiflx_x,psiflx_y,psiflx_z)
    else
      call bounduvw(cbcvel_halo,n,nh_v,nh_g,bcvel_halo,nb,is_bound,dl,dzc,dzf,uadv,vadv,wadv)
      do is=1,nscal
        call boundp(cbcscal(:,:,is),n,nh_s,nh_g,bcscal(:,:,is),nb,is_bound,dl,dzc,scalars(is)%sadv)
      end do
    end if
  end subroutine advance_split
  !
  subroutine ifc_markers_init(n,nh_v,nh_psi,nh_ifc,nh_g,dl,dli,dzc,dzf,dzci,dzfi,seps,is_restart,markers,psi)
    !
    ! initialization of the markers from their initial (or loaded) phase indicators, already present on the
    ! device: boundary conditions, normals and curvature of each, and the total phase indicator; in a
    ! restarted run, the total must reproduce the loaded one, which it does not if `inipsi` or `nmarker`
    ! differ from those of the run that saved the checkpoint (see `markers_offset`)
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v,nh_psi,nh_ifc,nh_g
    real(rp), intent(in   ), dimension(3) :: dl,dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzc,dzf,dzci,dzfi
    real(rp), intent(in   ) :: seps
    logical , intent(in   ) :: is_restart
    type(marker), intent(inout), dimension(:) :: markers
    real(rp), intent(inout), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp) :: vol_o,vol
    integer :: m
    !
    if(is_restart) vol_o = psi_sum(n,nh_psi,psi)
    if(is_ifc_plic) call vof_split_curvature_reset
    do m=1,size(markers)
      call boundp(cbcpsi,n,nh_psi,nh_g,bcpsi,nb,is_bound,dl,dzc,markers(m)%psi)
      if(is_ifc_split) call vof_split_contact_angle(n,nh_psi,nh_g,dl,dzc,dzf,cbcpsi,bcpsi,nb,is_bound,markers(m)%psi)
      call ifc_geometry(n,nh_psi,nh_ifc,nh_g,dl,dli,dzc,dzf,dzci,dzfi,seps,markers(m)%psi,markers(m)%phi, &
                        markers(m)%normx,markers(m)%normy,markers(m)%normz,markers(m)%kappa)
      call markers_sum(m == 1,markers(m)%psi,psi)
    end do
    call markers_total(n,nh_v,nh_psi,nh_ifc,size(markers),psi)
    if(is_restart) then
      vol = psi_sum(n,nh_psi,psi)
      if(abs(vol-vol_o) > 1.e-4_rp*product(1._rp*ng(:))) then
        if(myid == 0) print*, 'ERROR: the total phase indicator of the checkpoint is not that of its markers: `inipsi` &
                                      &(`bub*`, or `drp*`/`dis*`) and `nmarker` must be those of the run that saved it.'
        call MPI_FINALIZE(ierr)
        error stop
      end if
    end if
  end subroutine ifc_markers_init
  !
  real(rp) function psi_sum(n,nh_psi,psi)
    !
    ! sum of the phase indicator over the cells of the domain
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    integer , intent(in)               :: nh_psi
    real(rp), intent(in), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp) :: psis
    integer  :: i,j,k
    !
    psis = 0._rp
    !$acc data copy(psis) async(1)
    !$acc parallel loop collapse(3) default(present) reduction(+:psis) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          psis = psis + psi(i,j,k)
        end do
      end do
    end do
    !$acc end data
    !$acc wait(1)
    call MPI_ALLREDUCE(MPI_IN_PLACE,psis,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    psi_sum = psis
  end function psi_sum
  !
  subroutine ifc_markers_advance(tm_coeff,n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,dl,dli,dzc,dzf,dzci,dzfi,wzc,dt,dto,gam,seps,u,v,w, &
                                 is_first,istep,markers,psi,psiflx_x,psiflx_y,psiflx_z,uo,vo,wo,uadv,vadv,wadv,scalars)
    !
    ! advances each marker over one time step, as `ifc_advance` does for a single phase indicator, and returns
    ! the total phase indicator and the total phase fluxes (with their boundary conditions); with the
    ! direction-split advection, all the markers go through the sweeps together (see `vof_split_advance`)
    !
    implicit none
    real(rp), intent(in   ), dimension(2) :: tm_coeff
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v,nh_s,nh_psi,nh_ifc,nh_g
    real(rp), intent(in   ), dimension(3) :: dl,dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzc,dzf,dzci,dzfi,wzc
    real(rp), intent(in   ) :: dt,dto,gam,seps
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    logical , intent(in   ) :: is_first
    integer , intent(in   ) :: istep
    type(marker), intent(inout), target, dimension(:) :: markers
    real(rp), intent(inout), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: psiflx_x,psiflx_y,psiflx_z
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:), optional :: uo,vo,wo,uadv,vadv,wadv
    type(scalar), intent(inout), dimension(:), optional :: scalars
    type(vof_marker), dimension(size(markers)) :: psi_m
    integer :: m,nm
    !
    nm = size(markers)
    if(.not.is_ifc_split) then
      if(.not.allocated(wflx_x)) then
        allocate(wflx_x,wflx_y,wflx_z,mold=psiflx_x)
        !$acc enter data create(wflx_x,wflx_y,wflx_z)
      end if
      do m=1,nm
        call ifc_advance(tm_coeff,n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,dl,dli,dzc,dzf,dzci,dzfi,wzc,dt,dto,gam,seps,u,v,w, &
                         markers(m)%normx,markers(m)%normy,markers(m)%normz,markers(m)%phi,markers(m)%kappa,markers(m)%psi, &
                         wflx_x,wflx_y,wflx_z,is_first,markers(m)%dpsidtrko,istep)
        call markers_sum(m == 1,markers(m)%psi,psi)
        call markers_sum(m == 1,wflx_x,psiflx_x)
        call markers_sum(m == 1,wflx_y,psiflx_y)
        call markers_sum(m == 1,wflx_z,psiflx_z)
      end do
      call markers_total(n,nh_v,nh_psi,nh_ifc,nm,psi,u,v,w,psiflx_x,psiflx_y,psiflx_z)
    else
      do m=1,nm
        psi_m(m)%psi => markers(m)%psi
      end do
      call advance_split(tm_coeff,n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,dl,dzc,dzf,dzci,dzfi,wzc,dt,dto,u,v,w,psi, &
                         psiflx_x,psiflx_y,psiflx_z,is_first,istep,uo,vo,wo,uadv,vadv,wadv,scalars,psi_m,markers_offset(nm))
      do m=1,nm
        call ifc_geometry(n,nh_psi,nh_ifc,nh_g,dl,dli,dzc,dzf,dzci,dzfi,seps,markers(m)%psi,markers(m)%phi, &
                          markers(m)%normx,markers(m)%normy,markers(m)%normz,markers(m)%kappa)
        if(.not.is_vof_sweep_transport) call markers_sum(m == 1,markers(m)%psi,psi)
      end do
      !
      ! (with the transport along with the sweeps, the total phase indicator is formed there, and the phase
      ! fluxes are not used)
      !
      if(.not.is_vof_sweep_transport) then
        call markers_total(n,nh_v,nh_psi,nh_ifc,nm,psi,uadv,vadv,wadv,psiflx_x,psiflx_y,psiflx_z)
      end if
    end if
  end subroutine ifc_markers_advance
  !
  subroutine ifc_markers_force(n,nh_psi,nh_ifc,nh_g,dli,dzci,wzc,markers,psi,fsx,fsy,fsz)
    !
    ! surface tension force at the faces, summed over the markers (see `cmpt_capillary_force`)
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_psi,nh_ifc,nh_g
    real(rp), intent(in   ), dimension(3) :: dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzci,wzc
    type(marker), intent(in   ), dimension(:) :: markers
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(:,:,:) :: fsx,fsy,fsz
    integer :: m
    !
    do m=1,size(markers)
      call cmpt_capillary_force(n,nh_psi,nh_ifc,nh_g,dli,dzci,wzc,sigma,rho12,psi,markers(m)%psi,markers(m)%kappa,m == 1, &
                                fsx,fsy,fsz)
    end do
  end subroutine ifc_markers_force
  !
  integer function markers_offset(nm)
    !
    ! offset of the sum of the phase indicators of nm markers with respect to the total one (see `markers_total`)
    !
    implicit none
    integer, intent(in) :: nm
    !
    markers_offset = 0
    if(inipsi(1:3) == 'bub') markers_offset = nm-1
  end function markers_offset
  !
  subroutine markers_total(n,nh_v,nh_psi,nh_ifc,nm,psi,u,v,w,psiflx_x,psiflx_y,psiflx_z)
    !
    ! total phase indicator (and, if present, total phase fluxes) from the sum over the nm markers: if the
    ! dispersed phase is phase 1 (drops, psi = 1 inside), the sum itself; if it is phase 2 (bubbles, psi = 0
    ! inside), each marker is 1 outside its own bubbles, and the total is the sum minus nm-1 (for the fluxes,
    ! minus nm-1 times the velocity, the advecting one with the direction-split advection); the total phase
    ! indicator, which sets the fluid properties, is clipped to [0,1] where markers overlap, while the markers
    ! themselves are left as they are
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v,nh_psi,nh_ifc,nm
    real(rp), intent(inout), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)      , optional :: u,v,w
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:), optional :: psiflx_x,psiflx_y,psiflx_z
    real(rp) :: offset
    integer  :: i,j,k,nh
    !
    offset = markers_offset(nm)
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1-nh_psi,n(3)+nh_psi
      do j=1-nh_psi,n(2)+nh_psi
        do i=1-nh_psi,n(1)+nh_psi
          psi(i,j,k) = min(max(psi(i,j,k)-offset,0._rp),1._rp)
        end do
      end do
    end do
    if(present(psiflx_x) .and. offset /= 0._rp) then
      nh = min(nh_v,nh_ifc)
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1-nh,n(3)+nh
        do j=1-nh,n(2)+nh
          do i=1-nh,n(1)+nh
            psiflx_x(i,j,k) = psiflx_x(i,j,k) - offset*u(i,j,k)
            psiflx_y(i,j,k) = psiflx_y(i,j,k) - offset*v(i,j,k)
            psiflx_z(i,j,k) = psiflx_z(i,j,k) - offset*w(i,j,k)
          end do
        end do
      end do
    end if
  end subroutine markers_total
  !
  subroutine ifc_advecting_velocity(tm_coeff,is_first,istep,n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,dl,dzc,dzf,dzci,dzfi,wzc,dt,dto, &
                                    u,v,w,psi,psiflx_x,psiflx_y,psiflx_z,uo,vo,wo,uadv,vadv,wadv,scalars)
    !
    ! phase fluxes of a time step in which the interface is not tracked: those of a frozen phase indicator
    ! carried across the faces by the velocity (see `frozen_fluxes`), such that the mass fluxes of the
    ! transport of momentum and scalars are those of the local fluid; with the direction-split advection,
    ! also the advecting velocity, or, with `is_vof_sweep_transport`, the velocity (and the scalars)
    ! transported with the sweeps, without updating the phase indicator, returned in uadv, vadv, wadv (and in
    ! `sadv` of each scalar) instead
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
    type(scalar), intent(inout), dimension(:), optional :: scalars
    real(rp), dimension(2,max(nscal,1)) :: rhocp12
    real(rp), dimension(3) :: cfl
    integer :: is
    !
    if(.not.is_ifc_split) then
      call frozen_fluxes(n,nh_v,nh_psi,nh_ifc,nh_g,wzc,u,v,w,psi,psiflx_x,psiflx_y,psiflx_z)
      call bounduvw(cbcvel,n,nh_ifc,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,psiflx_x,psiflx_y,psiflx_z)
      return
    end if
    do is=1,max(nscal,1)
      rhocp12(:,is) = rho12(:)*cp12(:,is)
    end do
    if(is_vof_sweep_transport) then
      call vof_split_advance(tm_coeff,is_first,.false.,istep,n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,dl,dzc,dzf,dzci,dzfi,wzc,dt,dto, &
                             cbcpsi,bcpsi,nb,is_bound,rho12,rhocp12,u,v,w,psi,psiflx_x,psiflx_y,psiflx_z, &
                             uo,vo,wo,uadv,vadv,wadv,scalars)
      call bounduvw(cbcvel_halo,n,nh_v,nh_g,bcvel_halo,nb,is_bound,dl,dzc,dzf,uadv,vadv,wadv)
      do is=1,nscal
        call boundp(cbcscal(:,:,is),n,nh_s,nh_g,bcscal(:,:,is),nb,is_bound,dl,dzc,scalars(is)%sadv)
      end do
    else
      call vof_split_velocity(tm_coeff,is_first,n,nh_v,nh_g,dl,dzfi,dt,u,v,w,uo,vo,wo,uadv,vadv,wadv,cfl)
      call frozen_fluxes(n,nh_v,nh_psi,nh_ifc,nh_g,wzc,uadv,vadv,wadv,psi,psiflx_x,psiflx_y,psiflx_z)
      call bounduvw(cbcvel,n,nh_ifc,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,psiflx_x,psiflx_y,psiflx_z)
    end if
  end subroutine ifc_advecting_velocity
  !
  subroutine frozen_fluxes(n,nh_v,nh_psi,nh_ifc,nh_g,wzc,u,v,w,psi,psiflx_x,psiflx_y,psiflx_z)
    !
    ! phase fluxes of a phase indicator that is not advanced: its value at each face, interpolated with the
    ! weights of the density of the staggered control volumes, times the face velocity; exact for a uniform
    ! phase indicator (a single fluid)
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v,nh_psi,nh_ifc,nh_g
    real(rp), intent(in   ), dimension(1-nh_g:) :: wzc
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)       :: u,v,w
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: psiflx_x,psiflx_y,psiflx_z
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          psiflx_x(i,j,k) = 0.5*(psi(i,j,k)+psi(i+1,j,k))*u(i,j,k)
          psiflx_y(i,j,k) = 0.5*(psi(i,j,k)+psi(i,j+1,k))*v(i,j,k)
          psiflx_z(i,j,k) = (wzc(k)*psi(i,j,k)+(1.-wzc(k))*psi(i,j,k+1))*w(i,j,k)
        end do
      end do
    end do
  end subroutine frozen_fluxes
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
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: kappa
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:), optional :: normx,normy,normz ! (absent with `plic`)
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
          call thinc_cmpt_phi(n,nh_psi,nh_ifc,vof_thinc_beta,psi,phi)
        end if
        call cmpt_norm_curv(n,nh_ifc,nh_ifc,nh_g,dli,dzci,dzfi,phi,normx,normy,normz,kappa,is_bound,contact_angle) ! (halo
                                                                                                    ! width of phi)
      else
        call cmpt_norm_curv(n,nh_psi,nh_ifc,nh_g,dli,dzci,dzfi,psi,normx,normy,normz,kappa,is_bound,contact_angle)
      end if
      call boundp(cbccur,n,nh_ifc,nh_g,bccur,nb,is_bound,dl,dzc,kappa)
      if(.not.is_ifc_split) then
        !
        ! the halo of the normals is only read by the AB2 transport of the phase indicator
        !
        call boundp(cbcnor(:,:,1),n,nh_ifc,nh_g,bcnor(:,:,1),nb,is_bound,dl,dzc,normx)
        call boundp(cbcnor(:,:,2),n,nh_ifc,nh_g,bcnor(:,:,2),nb,is_bound,dl,dzc,normy)
        call boundp(cbcnor(:,:,3),n,nh_ifc,nh_g,bcnor(:,:,3),nb,is_bound,dl,dzc,normz)
      end if
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
      var(1)    = 1._rp*istep
      var(2)    = time
      var(3)    = vol
      var(4)    = (vol-vol0)/max(vol0,tiny(1._rp))
      var(5)    = clip
      var(6)    = fmin
      var(7)    = fmax
      var(8:10) = cfl(:)
      var(11)   = 1._rp*nsub
      var(12:14) = kstats(:)
      call out0d(trim(datadir)//'log_vof.out',14,var)
    else if(is_ifc_acdi) then
      var(1) = 1._rp*istep
      var(2) = time
      var(3) = gam
      var(4) = seps
      call out0d(trim(datadir)//'log_acdi.out',4,var(1:4))
    end if
  end subroutine ifc_write_log
end module mod_interface_capturing
