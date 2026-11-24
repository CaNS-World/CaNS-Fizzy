! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_vof_split
  !
  ! conservative direction-split advection of the volume fraction, following
  !   G.D. Weymouth & D.K.-P. Yue, J. Comput. Phys. 229 (2010) 2853--2865,
  ! with geometric donor-region fluxes from a PLIC reconstruction (see `mod_plic`), or with those of a THINC
  ! (hyperbolic-tangent) reconstruction along the direction of the sweep (see `mod_thinc`), also exactly
  ! conservative and bounded
  !
  ! each sweep along direction d updates
  !   f <- f - (F_{+} - F_{-})/h_d + c_c*dt*(u_{+} - u_{-})/h_d,
  ! where F is the volume of phase 1 per unit face area crossing a face during the sweep, and
  ! c_c = H(f^n-1/2) is frozen over the (sub-)step; the sum of the dilatation terms over the sweeps
  ! vanishes for a discretely divergence-free velocity, so the update is exactly conservative,
  ! and bounded for a CFL number not larger than 1/2 along each direction; larger time steps are
  ! split into sub-steps
  !
  ! the face fluxes summed over the sweeps (and sub-steps), divided by the time step, are returned
  ! for a mass-consistent transport of momentum and other conserved quantities
  !
  ! alternatively (`is_sweep_transport`), the momentum (and the scalar, with its heat capacity) is transported
  ! here, along with the sweeps of the volume fraction: with the same mass fluxes and dilatation term in each
  ! sweep (T.-Y. Huang, N. Valle, A.K. Lidtke, K. Hendrickson & G.D. Weymouth, arXiv:2606.02467, 2026;
  ! T. Arrufat et al., Comput. Fluids 215, 2021), and a transported velocity averaged over the mass-flux donor
  ! region of a linear reconstruction, limited near the interface (bounded) and unlimited elsewhere (low
  ! dissipation); mass and momentum then stay consistent after each sweep, which keeps the transport robust
  ! at large density ratios
  !
  ! the volume fraction needs two halo layers for the advection (the interface is also reconstructed in the
  ! first one, whose cells are the donors of the fluxes entering the subdomain), and four for the interface
  ! curvature, computed from height functions in nine-cell columns
  !
  use mpi
  use mod_types
  use, intrinsic :: iso_fortran_env, only: int8
  use mod_common_mpi, only: ierr
  use mod_bound     , only: boundp
  use mod_plic      , only: plic_alpha,plic_volume_box,plic_normal_myc,hf_curvature,hf_normal
  use mod_scal      , only: scalar
  use mod_thinc     , only: thinc_intercept,thinc_slab
  implicit none
  private
  public vof_split_init,vof_split_advance,vof_split_velocity,vof_split_curvature,vof_split_diagnostics, &
         vof_split_contact_angle
  !
  real(rp), parameter :: cfl_max = 0.5_rp    ! maximum CFL number of a sweep
  real(rp), parameter :: f_eps   = 1.e-12_rp ! cells with f_eps < f < 1-f_eps are reconstructed
  integer , parameter :: nq_hf   = 4         ! half-length of the height-function columns (nine cells)
  real(rp), parameter :: eps_hf  = 1.e-6_rp  ! tolerance for pure cells at the ends of the height-function columns
  real(rp), parameter :: kappa_undef = huge(1._rp) ! marks cells without a curvature estimate (intermediate stages)
  integer , parameter :: nh_mom     = 2      ! halo width of the fields transported with the sweeps (donor-region slopes)
  integer , parameter :: lim_band   = 3      ! slope limiter in the interface band: van Leer (1), MC (2), superbee (3)
  integer , parameter :: slope_bulk = 2      ! slope elsewhere: Fromm (1), Lax-Wendroff (2)
  !
  real(rp), allocatable, dimension(:,:,:) :: cc       ! dilatation indicator c_c
  real(rp), allocatable, dimension(:,:,:) :: fl       ! face fluxes of a sweep
  real(rp), allocatable, dimension(:,:,:) :: kap      ! work array for the curvature
  integer , allocatable, dimension(:)     :: ilist    ! compacted list of interfacial donor cells (packed indices)
  integer , dimension(1) :: nlist                     ! number of entries of ilist
  logical , dimension(3) :: is_active                 ! false along the directions in which the problem is invariant
  logical  :: is_thinc   = .false.                    ! THINC reconstruction instead of the geometric one
  real(rp) :: thinc_beta = 1._rp                      ! its sharpness
  real(rp), allocatable, dimension(:,:,:) :: mxs,mys,mzs ! normals of the geometric reconstruction from height functions,
                                                         ! in the donor cells of a sweep (allocated if they are used)
  !
  ! contact angles at walls, with the geometric reconstruction (see `vof_split_contact_angle`)
  !
  integer , parameter :: nr_ca = 2                    ! search radius (along the wall) for the plane of a ghost cell
  logical :: is_ca_any = .false.                      ! any wall with a contact angle other than 90 degrees
  logical , dimension(0:1,3) :: is_ca = .false.       ! walls with a contact angle other than 90 degrees
  real(rp), dimension(0:1,3) :: ca_cos = 0._rp,ca_sin = 1._rp ! cosine and sine of the contact angles
  real(rp), allocatable, dimension(:,:,:) :: ca_pl    ! planes of the wall cells (alpha,m; m = 0 without interface)
  !
  ! transport of momentum (and scalar) along with the sweeps
  !
  logical :: is_sweep = .false.                       ! momentum (and scalar) transported along with the sweeps
  integer :: nscal_tr = 0                             ! number of scalars transported too
  logical :: is_scal  = .false.                       ! any
  integer :: ncc = 0                                  ! layers beyond the interior with dilatation indicator and fluxes
  real(rp)     , allocatable, dimension(:,:,:)   :: f_prev      ! volume fraction before a sweep
  integer(int8), allocatable, dimension(:,:,:)   :: band        ! interface band (limited transported value), 0 or 1
  integer(int8), allocatable, dimension(:,:,:)   :: bflag       ! presence of the phases in a cell (bits 0 and 1)
  real(rp)     , allocatable, dimension(:,:,:)   :: ub,vb,wb    ! velocity during the sweeps
  real(rp)     , allocatable, dimension(:,:,:,:) :: sb          ! scalars during the sweeps
  real(rp)     , allocatable, dimension(:,:,:,:) :: qflx        ! fluxes of the transported quantities
  real(rp)     , allocatable, dimension(:,:,:)   :: au,av,aw    ! non-advective rates of change of the previous time step
  real(rp)     , allocatable, dimension(:,:,:,:) :: as          ! the same, for the scalars
  !
  ! diagnostics (local to the task, accumulated on the device; reduced by `vof_split_diagnostics`)
  !
  real(rp), dimension(3) :: clip_stats = [0._rp,0._rp,1._rp] ! accumulated volume added by clipping, and extrema
                                                             ! (minimum, maximum) before clipping in the last step
  real(rp), dimension(3) :: clip_sweep          ! the same, for one sweep
  real(rp), dimension(3) :: cfl_last = 0._rp    ! CFL numbers of the advecting velocity in the last step
  integer  :: nsub_last = 0                     ! number of sub-steps in the last step
  real(rp), dimension(3) :: kappa_stats = 0._rp ! interfacial cells with curvature from height functions,
                                                ! from neighbouring height functions, and without estimate
  contains
  subroutine vof_split_init(ng,n,nh_psi,is_thinc_reconstruction,beta,is_hf_normals,is_sweep_transport,nscal,theta_c)
    !
    ! allocates the work arrays; `is_thinc_reconstruction` selects the THINC reconstruction, of sharpness
    ! `beta`, instead of the geometric one, `is_hf_normals` the normals of the geometric reconstruction from
    ! height functions, and `is_sweep_transport` the transport of momentum (and of the
    ! `nscal` scalars) along with the sweeps; `theta_c` are the contact angles at the walls (in degrees,
    ! measured through phase 1), imposed here with the geometric reconstruction
    !
    implicit none
    integer , intent(in), dimension(3) :: ng,n
    integer , intent(in)               :: nh_psi
    logical , intent(in)               :: is_thinc_reconstruction,is_hf_normals,is_sweep_transport
    integer , intent(in)               :: nscal
    real(rp), intent(in)               :: beta
    real(rp), intent(in), dimension(0:1,3) :: theta_c
    integer :: i,j,k
    !
    is_thinc   = is_thinc_reconstruction
    thinc_beta = beta
    if(nh_psi < nq_hf) error stop 'ERROR: the halo width of the volume fraction is smaller than the height-function columns.'
    !
    ! directions with one or two grid cells are taken as invariant (2D runs): no sweeps,
    ! and no contribution to the interface normal
    !
    is_active(:) = ng(:) > 2
    allocate(cc(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    allocate(fl,kap,mold=cc)
    is_sweep = is_sweep_transport
    nscal_tr = merge(nscal,0,is_sweep)
    is_scal  = nscal_tr > 0
    ncc = merge(1,0,is_sweep)
    allocate(ilist((n(1)+2+ncc)*(n(2)+2+ncc)*(n(3)+2+ncc)))
    nlist(:) = 0
    !$acc enter data create(cc,fl,kap,ilist) copyin(nlist)
    if(is_hf_normals) then
      !
      ! height-function normals: the columns of the donor cells in the halo layers need a wider halo (the
      ! arrays of the normals are otherwise not allocated, and absent arguments of the routines that take them)
      !
      if(nh_psi < nq_hf+1+ncc) error stop 'ERROR: the halo width of the volume fraction is too small for the &
                                          &height-function normals.'
      allocate(mxs(0:n(1)+1+ncc,0:n(2)+1+ncc,0:n(3)+1+ncc))
      allocate(mys,mzs,mold=mxs)
      !$acc enter data create(mxs,mys,mzs)
    end if
    if(is_sweep) then
      allocate(f_prev(1-nh_psi:n(1)+nh_psi,1-nh_psi:n(2)+nh_psi,1-nh_psi:n(3)+nh_psi))
      allocate(bflag( 1-nh_psi:n(1)+nh_psi,1-nh_psi:n(2)+nh_psi,1-nh_psi:n(3)+nh_psi))
      allocate(band,mold=bflag)
      allocate(ub(1-nh_mom:n(1)+nh_mom,1-nh_mom:n(2)+nh_mom,1-nh_mom:n(3)+nh_mom))
      allocate(vb,wb,mold=ub)
      allocate(qflx(0:n(1)+1,0:n(2)+1,0:n(3)+1,merge(0,1,is_scal):3))
      allocate(au,av,aw,mold=cc)
      if(is_scal) then
        allocate(sb(1-nh_mom:n(1)+nh_mom,1-nh_mom:n(2)+nh_mom,1-nh_mom:n(3)+nh_mom,nscal_tr))
        allocate(as(0:n(1)+1,0:n(2)+1,0:n(3)+1,nscal_tr))
      else
        allocate(sb(1,1,1,1),as(1,1,1,1)) ! (not used)
      end if
      !$acc enter data create(f_prev,bflag,band,ub,vb,wb,sb,qflx,au,av,aw,as)
      !
      ! (the outer halo layers are not updated at physical boundaries, and must be finite)
      !
      call set_zero(ub)
      call set_zero(vb)
      call set_zero(wb)
      do i=1,size(sb,4)
        call set_zero(sb(:,:,:,i))
      end do
      call set_zero(f_prev)
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1-nh_psi,n(3)+nh_psi
        do j=1-nh_psi,n(2)+nh_psi
          do i=1-nh_psi,n(1)+nh_psi
            band(i,j,k) = 0_int8
          end do
        end do
      end do
    end if
    !
    ! contact angles (walls with an angle other than 90 degrees within round-off)
    !
    if(.not.is_thinc) then
      ca_cos(:,:) = cos(theta_c(:,:)*acos(-1._rp)/180._rp)
      ca_sin(:,:) = sin(theta_c(:,:)*acos(-1._rp)/180._rp)
      is_ca(:,:)  = abs(ca_cos(:,:)) > 1.e-12_rp
      is_ca_any   = any(is_ca(:,:))
      if(is_ca_any) then
        allocate(ca_pl(0:3,1-nh_psi:maxval(n(:))+nh_psi,1-nh_psi:maxval(n(:))+nh_psi))
        !$acc enter data create(ca_pl)
      end if
    end if
    !$acc enter data copyin(clip_stats,kappa_stats) create(clip_sweep)
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          fl(i,j,k) = 0._rp
        end do
      end do
    end do
  end subroutine vof_split_init
  !
  subroutine vof_split_advance(tm_coeff,is_first,is_track,istep,n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,dl,dzc,dzf,dzci,dzfi,wzc,dt,dto, &
                               cbcpsi,bcpsi,nb,is_bound,rho12,rhocp12,u,v,w,psi,psiflx_x,psiflx_y,psiflx_z, &
                               uo,vo,wo,ua,va,wa,scalars)
    !
    ! advances the volume fraction over one time step, with the face velocity extrapolated to the
    ! middle of the time step from the current and previous velocity fields (both divergence free);
    ! the volume fraction is returned with its halo layers updated, and the advecting velocity in ua, va, wa;
    ! the order of the sweeps is reversed at each (sub-)step, starting from the parity of the time step number
    !
    ! with the transport along with the sweeps (`is_sweep_transport`), the velocity (and each of the
    ! `scalars`, with the heat capacities per unit volume `rhocp12` of the two phases) is transported too, and
    ! ua, va, wa (and `sadv` of each scalar) return the transported velocity (scalar) instead; on entry they
    ! hold that of the previous time step (unless `is_first`), whose time step is dto; the phase fluxes are
    ! then not used; if `is_track` is false, the volume fraction is not updated (transport without interface
    ! tracking)
    !
    implicit none
    real(rp), intent(in   ), dimension(2) :: tm_coeff
    logical , intent(in   )               :: is_first,is_track
    integer , intent(in   )               :: istep
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v,nh_s,nh_psi,nh_ifc,nh_g
    real(rp), intent(in   ), dimension(3) :: dl
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzc,dzf,dzci,dzfi,wzc
    real(rp), intent(in   ) :: dt,dto
    character(len=1), intent(in), dimension(0:1,3) :: cbcpsi
    real(rp), intent(in   ), dimension(0:1,3) :: bcpsi
    integer , intent(in   ), dimension(0:1,3) :: nb
    logical , intent(in   ), dimension(0:1,3) :: is_bound
    real(rp), intent(in   ), dimension(2) :: rho12
    real(rp), intent(in   ), dimension(:,:) :: rhocp12 ! (2,nscal)
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)       :: u,v,w
    real(rp), intent(inout), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: psiflx_x,psiflx_y,psiflx_z
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:)       :: uo,vo,wo,ua,va,wa
    type(scalar), intent(inout), dimension(:), optional :: scalars
    character(len=1), dimension(0:1,3) :: cbc_tr
    logical , dimension(0:1,3) :: is_bound_tr,is_bound_phys
    logical , save :: is_psiflx_set = .false.
    real(rp) :: dts,dtpi,tc
    integer  :: i,j,k,isub,nsub,q,idir,is
    logical  :: is_forward
    !
    if(is_sweep) then
      !
      ! non-advective rates of change of the previous time step (pressure, viscous and body forces), from the
      ! transported velocity (scalar) of the previous step, still in ua, va, wa (sa), used to centre the
      ! transported values in time (see `transport_sweep`); the velocity (scalar) is also copied into the
      ! work arrays of the transport, whose outer halo layers are then updated
      !
      dtpi = 0._rp
      if(.not.is_first) dtpi = 1._rp/dto
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=0,n(3)+1
        do j=0,n(2)+1
          do i=0,n(1)+1
            au(i,j,k) = 0._rp
            av(i,j,k) = 0._rp
            aw(i,j,k) = 0._rp
            if(dtpi > 0._rp) then
              au(i,j,k) = (u(i,j,k)-ua(i,j,k))*dtpi
              av(i,j,k) = (v(i,j,k)-va(i,j,k))*dtpi
              aw(i,j,k) = (w(i,j,k)-wa(i,j,k))*dtpi
            end if
            ub(i,j,k) = u(i,j,k)
            vb(i,j,k) = v(i,j,k)
            wb(i,j,k) = w(i,j,k)
          end do
        end do
      end do
      do is=1,nscal_tr
        call transport_scal_begin(n,nh_s,dtpi,scalars(is)%val,scalars(is)%sadv,as(:,:,:,is),sb(:,:,:,is))
      end do
      !
      ! the halos of the transported fields are only needed across periodic and inter-task boundaries: at
      ! physical boundaries, the transported value uses zero slopes, and the values of the beginning of the
      ! time step are kept
      !
      cbc_tr(:,:) = 'N'
      where(cbcpsi(0,:)//cbcpsi(1,:) == 'PP')
        cbc_tr(0,:) = 'P'
        cbc_tr(1,:) = 'P'
      end where
      is_bound_tr(  :,:) = is_bound(:,:).and.cbc_tr(:,:) == 'P'
      is_bound_phys(:,:) = is_bound(:,:).and.cbc_tr(:,:) /= 'P'
      call transport_halo(cbc_tr,n,nh_g,nb,is_bound_tr,dl,dzc)
      !
      ! the phase fluxes are not used: they are set to zero once
      !
      if(.not.is_psiflx_set) then
        call set_zero(psiflx_x)
        call set_zero(psiflx_y)
        call set_zero(psiflx_z)
        is_psiflx_set = .true.
      end if
      if(.not.is_track) call set_zero(fl)
    end if
    !
    ! advecting velocity, and CFL numbers along each direction
    !
    call vof_split_velocity(tm_coeff,is_first,n,nh_v,nh_g,dl,dzfi,dt,u,v,w,uo,vo,wo,ua,va,wa,cfl_last)
    call MPI_ALLREDUCE(MPI_IN_PLACE,cfl_last,3,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    !
    ! number of sub-steps
    !
    nsub = 1
    if(any(is_active(:))) nsub = max(1,ceiling(maxval(cfl_last(:),mask=is_active(:))/cfl_max))
    nsub_last = nsub
    dts = dt/nsub
    !
    ! the face fluxes are set by the first sweep of each direction (and accumulated over the sub-steps),
    ! except along the directions without sweeps
    !
    if(.not.is_sweep) then
      if(.not.is_active(1)) call set_zero(psiflx_x)
      if(.not.is_active(2)) call set_zero(psiflx_y)
      if(.not.is_active(3)) call set_zero(psiflx_z)
    end if
    !
    if(is_track) then
      !$acc serial default(present) async(1)
      clip_stats(2) =  huge(1._rp)
      clip_stats(3) = -huge(1._rp)
      !$acc end serial
    end if
    do isub=1,nsub
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1-ncc,n(3)+ncc
        do j=1-ncc,n(2)+ncc
          do i=1-ncc,n(1)+ncc
            cc(i,j,k) = merge(1._rp,0._rp,psi(i,j,k) > 0.5_rp)
          end do
        end do
      end do
      is_forward = mod(istep+isub,2) == 0
      tc = (isub-0.5_rp)*dts ! time from the beginning of the step to the middle of the sub-step
      do q=1,3
        idir = merge(q,4-q,is_forward)
        if(.not.is_active(idir)) cycle
        if(is_sweep) call transport_prepare(n,nh_psi,psi)
        if(is_track) then
          select case(idir)
          case(1)
            call sweep(1,n,nh_v,nh_psi,nh_ifc,nh_g,dl,dzf,dzfi,dts,ua,isub == 1,1._rp/dt,psiflx_x,psi, &
                       mxs,mys,mzs)
          case(2)
            call sweep(2,n,nh_v,nh_psi,nh_ifc,nh_g,dl,dzf,dzfi,dts,va,isub == 1,1._rp/dt,psiflx_y,psi, &
                       mxs,mys,mzs)
          case(3)
            call sweep(3,n,nh_v,nh_psi,nh_ifc,nh_g,dl,dzf,dzfi,dts,wa,isub == 1,1._rp/dt,psiflx_z,psi, &
                       mxs,mys,mzs)
          end select
          call boundp(cbcpsi,n,nh_psi,nh_g,bcpsi,nb,is_bound,dl,dzc,psi)
          call vof_split_contact_angle(n,nh_psi,nh_g,dl,dzc,dzf,cbcpsi,bcpsi,nb,is_bound,psi)
        end if
        if(is_sweep) then
          select case(idir)
          case(1)
            call transport_all(1,ua)
          case(2)
            call transport_all(2,va)
          case(3)
            call transport_all(3,wa)
          end select
          call transport_halo(cbc_tr,n,nh_g,nb,is_bound_tr,dl,dzc)
        end if
      end do
    end do
    if(is_sweep) then
      !
      ! transported velocity (and scalar)
      !
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=0,n(3)+1
        do j=0,n(2)+1
          do i=0,n(1)+1
            ua(i,j,k) = ub(i,j,k)
            va(i,j,k) = vb(i,j,k)
            wa(i,j,k) = wb(i,j,k)
          end do
        end do
      end do
      do is=1,nscal_tr
        call transport_scal_end(n,nh_s,sb(:,:,:,is),scalars(is)%sadv)
      end do
    end if
  contains
    subroutine transport_all(idir,vel)
      !
      ! transport along direction idir, with the face velocity vel, of the velocity together with the first
      ! scalar, if any, and then of each of the other scalars
      !
      implicit none
      integer , intent(in) :: idir
      real(rp), intent(in), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: vel
      integer :: is
      if(nscal_tr == 0) then
        call transport_sweep(idir,n,nh_v,nh_s,nh_psi,nh_g,dl,dzf,dzci,dzfi,wzc,dts,tc,rho12,rhocp12(:,1),is_bound_phys, &
                             vel,u,v,w,psi,1,3,sb(:,:,:,1),as(:,:,:,1))
      else
        call transport_sweep(idir,n,nh_v,nh_s,nh_psi,nh_g,dl,dzf,dzci,dzfi,wzc,dts,tc,rho12,rhocp12(:,1),is_bound_phys, &
                             vel,u,v,w,psi,0,3,sb(:,:,:,1),as(:,:,:,1),scalars(1)%val)
        do is=2,nscal_tr
          call transport_sweep(idir,n,nh_v,nh_s,nh_psi,nh_g,dl,dzf,dzci,dzfi,wzc,dts,tc,rho12,rhocp12(:,is),is_bound_phys, &
                               vel,u,v,w,psi,0,0,sb(:,:,:,is),as(:,:,:,is),scalars(is)%val)
        end do
      end if
    end subroutine transport_all
  end subroutine vof_split_advance
  !
  subroutine transport_scal_begin(n,nh_s,dtpi,s,sa,asq,sbq)
    !
    ! non-advective rate of change of a scalar over the previous time step (see `vof_split_advance`), and
    ! its copy into the work array of the transport
    !
    implicit none
    integer , intent(in ), dimension(3) :: n
    integer , intent(in )               :: nh_s
    real(rp), intent(in )               :: dtpi
    real(rp), intent(in ), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s,sa
    real(rp), intent(out), dimension(0:,0:,0:) :: asq
    real(rp), intent(out), dimension(1-nh_mom:,1-nh_mom:,1-nh_mom:) :: sbq
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          asq(i,j,k) = 0._rp
          if(dtpi > 0._rp) asq(i,j,k) = (s(i,j,k)-sa(i,j,k))*dtpi
          sbq(i,j,k) = s(i,j,k)
        end do
      end do
    end do
  end subroutine transport_scal_begin
  !
  subroutine transport_scal_end(n,nh_s,sbq,sa)
    !
    ! transported scalar
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_s
    real(rp), intent(in   ), dimension(1-nh_mom:,1-nh_mom:,1-nh_mom:) :: sbq
    real(rp), intent(inout), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: sa
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          sa(i,j,k) = sbq(i,j,k)
        end do
      end do
    end do
  end subroutine transport_scal_end
  !
  subroutine vof_split_velocity(tm_coeff,is_first,n,nh_v,nh_g,dl,dzfi,dt,u,v,w,uo,vo,wo,ua,va,wa,cfl)
    !
    ! advecting velocity, extrapolated to the middle of the time step from the current and previous
    ! (divergence-free) velocity fields, and its (local) maximum CFL numbers along each direction;
    ! the current velocity is stored in uo, vo, wo for the next time step; if `is_first`, the previous
    ! velocity is not available, and the advecting velocity is the current one
    !
    implicit none
    real(rp), intent(in   ), dimension(2) :: tm_coeff
    logical , intent(in   )               :: is_first
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v,nh_g
    real(rp), intent(in   ), dimension(3) :: dl
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzfi
    real(rp), intent(in   ) :: dt
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: uo,vo,wo,ua,va,wa
    real(rp), intent(out  ), dimension(3) :: cfl
    real(rp) :: c1,c2,dtx,dty,cflx,cfly,cflz,uaux,vaux,waux
    integer  :: i,j,k,n1,n2,n3
    !
    c1 = tm_coeff(1)
    c2 = tm_coeff(2)
    if(is_first) then
      c1 = 1._rp
      c2 = 0._rp
    end if
    dtx = dt/dl(1)
    dty = dt/dl(2)
    n1 = n(1); n2 = n(2); n3 = n(3)
    cflx = 0._rp
    cfly = 0._rp
    cflz = 0._rp
    !$acc parallel loop collapse(3) default(present) private(uaux,vaux,waux) reduction(max:cflx,cfly,cflz) async(1)
    do k=0,n3+1
      do j=0,n2+1
        do i=0,n1+1
          uaux = c1*u(i,j,k)+c2*uo(i,j,k)
          vaux = c1*v(i,j,k)+c2*vo(i,j,k)
          waux = c1*w(i,j,k)+c2*wo(i,j,k)
          ua(i,j,k) = uaux
          va(i,j,k) = vaux
          wa(i,j,k) = waux
          uo(i,j,k) = u(i,j,k)
          vo(i,j,k) = v(i,j,k)
          wo(i,j,k) = w(i,j,k)
          if(i <= n1 .and. j <= n2 .and. k <= n3) then
            cflx = max(cflx,abs(uaux)*dtx)
            cfly = max(cfly,abs(vaux)*dty)
            cflz = max(cflz,abs(waux)*dt*max(dzfi(k),dzfi(k+1)))
          end if
        end do
      end do
    end do
    !$acc wait(1)
    cfl(:) = [cflx,cfly,cflz]
  end subroutine vof_split_velocity
  !
  subroutine sweep(idir,n,nh_v,nh_psi,nh_ifc,nh_g,dl,dzf,dzfi,dt,vel,is_init,scl,psiflx,psi,mx_hf,my_hf,mz_hf)
    !
    ! one directional sweep: geometric donor-region fluxes, and conservative update with the dilatation term;
    ! the face fluxes, times scl, are stored in psiflx if is_init, and accumulated otherwise; the extrema of
    ! the volume fraction before clipping and the clipped volume are accumulated in `clip_stats`;
    ! mx_hf, my_hf, mz_hf are the work arrays of the height-function normals, present if these are used
    !
    implicit none
    integer , intent(in   ) :: idir
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v,nh_psi,nh_ifc,nh_g
    real(rp), intent(in   ), dimension(3) :: dl
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzf,dzfi
    real(rp), intent(in   ) :: dt
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: vel
    logical , intent(in   ) :: is_init
    real(rp), intent(in   ) :: scl
    real(rp), intent(inout), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: psiflx
    real(rp), intent(inout), dimension(0:,0:,0:), optional :: mx_hf,my_hf,mz_hf
    real(rp), dimension(-1:1,-1:1,-1:1) :: c
    real(rp), dimension(-1:1) :: dz3
    logical , dimension(3) :: act
    logical  :: is_own_p,is_own_m,is_act1,is_act2,is_act3,is_flx,is_th,is_hfn
    real(rp) :: dx,dy,dxi,dyi,h,rh,fc,vp,vm,cfl_p,cfl_m,frac_p,frac_m,fnew,fcl,mx,my,mz,beta
    integer  :: i,j,k,di,dj,dk,ii,nd,ndf,ilo,jlo,klo,ihi,jhi,khi,ni,nj,il,lin
    !
    di = merge(1,0,idir == 1)
    dj = merge(1,0,idir == 2)
    dk = merge(1,0,idir == 3)
    nd = n(idir)
    dx = dl(1)
    dy = dl(2)
    dxi = 1._rp/dx
    dyi = 1._rp/dy
    is_act1 = is_active(1)
    is_act2 = is_active(2)
    is_act3 = is_active(3)
    ndf = nd+ncc
    is_flx = .not.is_sweep
    is_th  = is_thinc
    is_hfn = present(mx_hf)
    beta   = thinc_beta
    ilo = 1-di; ihi = n(1)+di+ncc
    jlo = 1-dj; jhi = n(2)+dj+ncc
    klo = 1-dk; khi = n(3)+dk+ncc
    !
    ! fluxes: each face is computed by its donor (upwind) cell; the fluxes of the (few) interfacial donor
    ! cells, which reconstruct the interface, are computed by a second kernel over a compacted list of
    ! these cells, such that the work of the reconstruction is not spread over all the GPU threads; with the
    ! transport along with the sweeps, also the fluxes through the faces 0:n+1 along idir and 1:n+1 along the
    ! other directions are computed
    !
    ni = ihi-ilo+1
    nj = jhi-jlo+1
    nlist(1) = 0
    !$acc update device(nlist) async(1)
    !$acc parallel loop collapse(3) default(present) private(il) async(1)
    do k=klo,khi
      do j=jlo,jhi
        do i=ilo,ihi
          ii = i*di+j*dj+k*dk
          vp = 0._rp
          if(ii <= ndf) vp = vel(i,j,k)
          vm = 0._rp
          if(ii >= 1) vm = vel(i-di,j-dj,k-dk)
          is_own_p = ii <= ndf .and. vp >= 0._rp
          is_own_m = ii >= 1  .and. vm <  0._rp
          if(is_own_p .or. is_own_m) then
            fc = min(max(psi(i,j,k),0._rp),1._rp)
            if(is_own_p) fl(i,j,k) = fc*vp*dt
            if(is_own_m) fl(i-di,j-dj,k-dk) = fc*vm*dt
            if(fc > f_eps .and. fc < 1._rp-f_eps) then
              !$acc atomic capture
              nlist(1) = nlist(1)+1
              il = nlist(1)
              !$acc end atomic
              ilist(il) = (i-ilo)+ni*((j-jlo)+nj*(k-klo))
            end if
          end if
        end do
      end do
    end do
    !
    ! normals of the listed donor cells from height functions, where valid
    !
    if(is_hfn) call hf_list(ilo,jlo,klo,ni,nj,nh_psi,nh_ifc,nh_g,dl,dzf,psi,mx=mx_hf,my=my_hf,mz=mz_hf)
    !
    ! (the number of listed cells is read from the device copy of nlist, without a host synchronization)
    !
    !$acc parallel loop default(present) private(c,dz3,act) async(1)
    do il=1,nlist(1)
      lin = ilist(il)
      i = ilo+mod(lin,ni)
      j = jlo+mod(lin/ni,nj)
      k = klo+lin/(ni*nj)
      ii = i*di+j*dj+k*dk
      h  = merge(dx,merge(dy,dzf(k),idir == 2),idir == 1)
      vp = 0._rp
      if(ii <= ndf) vp = vel(i,j,k)
      vm = 0._rp
      if(ii >= 1) vm = vel(i-di,j-dj,k-dk)
      is_own_p = ii <= ndf .and. vp >= 0._rp
      is_own_m = ii >= 1  .and. vm <  0._rp
      fc = min(max(psi(i,j,k),0._rp),1._rp)
      cfl_p = 0._rp
      cfl_m = 0._rp
      if(is_own_p) cfl_p =  vp*dt/h
      if(is_own_m) cfl_m = -vm*dt/h
      act(:) = [is_act1,is_act2,is_act3]
      if(is_hfn) then
        mx = mx_hf(i,j,k); my = my_hf(i,j,k); mz = mz_hf(i,j,k)
      else
        c(:,:,:) = psi(i-1:i+1,j-1:j+1,k-1:k+1)
        dz3(:) = dzf(k-1:k+1)
        call plic_normal_myc(c,dx,dy,dz3,act,mx,my,mz)
      end if
      call donor_fractions(mx,my,mz,fc,idir,is_th,beta,cfl_p,cfl_m,frac_p,frac_m)
      if(is_own_p) fl(i,j,k) = frac_p*vp*dt
      if(is_own_m) fl(i-di,j-dj,k-dk) = frac_m*vm*dt
    end do
    !
    ! update, with the dilatation term, clipping of round-off errors, and accumulation of the fluxes
    !
    ! (the extrema start from 0 and 1, and only the few values out of these bounds, the ones clipped, update
    ! them and the clipped volume, with atomics instead of reductions over the whole domain)
    !
    !$acc serial default(present) async(1)
    clip_sweep(1) = 0._rp
    clip_sweep(2) = 0._rp
    clip_sweep(3) = 1._rp
    !$acc end serial
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          ii = i*di+j*dj+k*dk
          rh = merge(dxi,merge(dyi,dzfi(k),idir == 2),idir == 1)
          fnew = psi(i,j,k)-(fl(i,j,k)-fl(i-di,j-dj,k-dk))*rh &
                 +cc(i,j,k)*(vel(i,j,k)-vel(i-di,j-dj,k-dk))*(dt*rh)
          fcl  = min(max(fnew,0._rp),1._rp)
          if(fcl /= fnew) then
            !$acc atomic update
            clip_sweep(2) = min(clip_sweep(2),fnew)
            !$acc atomic update
            clip_sweep(3) = max(clip_sweep(3),fnew)
            !$acc atomic update
            clip_sweep(1) = clip_sweep(1)+(fcl-fnew)*dx*dy*dzf(k)
          end if
          psi(i,j,k) = fcl
          if(is_flx) then
            if(is_init) then
              psiflx(i,j,k) = fl(i,j,k)*scl
              if(ii == 1) psiflx(i-di,j-dj,k-dk) = fl(i-di,j-dj,k-dk)*scl
            else
              psiflx(i,j,k) = psiflx(i,j,k)+fl(i,j,k)*scl
              if(ii == 1) psiflx(i-di,j-dj,k-dk) = psiflx(i-di,j-dj,k-dk)+fl(i-di,j-dj,k-dk)*scl
            end if
          end if
        end do
      end do
    end do
    !$acc serial default(present) async(1)
    clip_stats(1) = clip_stats(1)+clip_sweep(1)
    clip_stats(2) = min(clip_stats(2),clip_sweep(2))
    clip_stats(3) = max(clip_stats(3),clip_sweep(3))
    !$acc end serial
  end subroutine sweep
  !
  subroutine transport_halo(cbc,n,nh_g,nb,is_bound,dl,dzc)
    !
    ! halo update of the fields transported with the sweeps, across periodic and inter-task boundaries
    ! only: at physical boundaries, the halo values of the beginning of the time step are kept
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    integer , intent(in), dimension(3) :: n
    integer , intent(in)               :: nh_g
    integer , intent(in), dimension(0:1,3) :: nb
    logical , intent(in), dimension(0:1,3) :: is_bound
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(1-nh_g:) :: dzc
    real(rp), dimension(0:1,3) :: bc
    integer :: is
    !
    bc(:,:) = 0._rp
    call boundp(cbc,n,nh_mom,nh_g,bc,nb,is_bound,dl,dzc,ub)
    call boundp(cbc,n,nh_mom,nh_g,bc,nb,is_bound,dl,dzc,vb)
    call boundp(cbc,n,nh_mom,nh_g,bc,nb,is_bound,dl,dzc,wb)
    do is=1,nscal_tr
      call boundp(cbc,n,nh_mom,nh_g,bc,nb,is_bound,dl,dzc,sb(:,:,:,is))
    end do
  end subroutine transport_halo
  !
  subroutine transport_prepare(n,nh_psi,psi)
    !
    ! volume fraction before a sweep, and interface band: cells whose 3x3x3 neighbourhood contains
    ! interfacial cells, or both full and empty cells, i.e., cells with phase 1 (f > f_eps) and cells with
    ! phase 2 (f < 1-f_eps); the presence of the two phases is first stored as one-byte flags (bits 0 and 1)
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    integer , intent(in)               :: nh_psi
    real(rp), intent(in), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    integer  :: i,j,k,ii,jj,kk
    integer(int8) :: ib
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=-1,n(3)+3
      do j=-1,n(2)+3
        do i=-1,n(1)+3
          bflag(i,j,k) = merge(1_int8,0_int8,psi(i,j,k) > f_eps)+merge(2_int8,0_int8,psi(i,j,k) < 1._rp-f_eps)
          f_prev(i,j,k) = psi(i,j,k)
        end do
      end do
    end do
    !$acc parallel loop collapse(3) default(present) private(ib) async(1)
    do k=0,n(3)+2
      do j=0,n(2)+2
        do i=0,n(1)+2
          ib = 0_int8
          !$acc loop seq
          do kk=-1,1
            !$acc loop seq
            do jj=-1,1
              !$acc loop seq
              do ii=-1,1
                ib = ior(ib,bflag(i+ii,j+jj,k+kk))
              end do
            end do
          end do
          band(i,j,k) = merge(1_int8,0_int8,ib == 3_int8)
        end do
      end do
    end do
  end subroutine transport_prepare
  !
  subroutine transport_sweep(idir,n,nh_v,nh_s,nh_psi,nh_g,dl,dzf,dzci,dzfi,wzc,dt,tc,rho12,rhocp12,is_bound_phys, &
                             vel,u,v,w,psi,iq_lo,iq_hi,sbq,asq,s)
    !
    ! conservative transport along direction idir, after the volume-fraction update of the sweep, of the
    ! quantities per unit mass iq = iq_lo,...,iq_hi: the three velocity components (iq = 1,2,3) on their
    ! staggered control volumes and a cell-centred scalar (iq = 0; its work arrays are sbq and asq, parts of sb
    ! and as, and s its value at the beginning of the time step), with the masses per unit volume of
    ! the two phases rho12 (densities) and rhocp12 (heat capacities), respectively; the transport is
    ! synchronized with the volume-fraction sweep (Huang et al., arXiv:2606.02467):
    !  - the mass fluxes through the control-volume faces are those of the volume fraction, averaged as the
    !    masses of the control volumes (which are averages of the cell masses, see `rk`);
    !  - the dilatation term is that of the volume fraction (averaged likewise) times the value of q at the
    !    beginning of the time step (u, v, w and s), such that it cancels over the sweeps;
    !  - the transported value is the average over the donor region (the mass that crosses the face) of a
    !    linear reconstruction in the upwind control volume, with limited slopes in the interface band
    !    (bounded) and unlimited ones elsewhere (low dissipation), see `donor_value`; the slopes vanish next
    !    to physical boundaries;
    !  - the transported value is centred in time with the non-advective rate of change of q of the previous
    !    time step (pressure, viscous and body forces; au, av, aw, asq), over the time tc from the beginning of
    !    the step to the middle of the sub-step; without it, the transport of q^n followed by the other terms
    !    is a first-order (Lie) splitting, with a significant dissipation of kinetic energy;
    ! the new value is the new momentum divided by the new mass, computed from the updated volume fraction; it
    ! is stored in ub, vb, wb, sbq (the update of a control volume only reads its own value, the fluxes being
    ! computed beforehand); the normal velocity at the last face of a physical boundary (a boundary value)
    ! is kept
    !
    implicit none
    integer , intent(in ) :: idir
    integer , intent(in ), dimension(3) :: n
    integer , intent(in )               :: nh_v,nh_s,nh_psi,nh_g
    real(rp), intent(in ), dimension(3) :: dl
    real(rp), intent(in ), dimension(1-nh_g:) :: dzf,dzci,dzfi,wzc
    real(rp), intent(in ) :: dt,tc
    real(rp), intent(in ), dimension(2) :: rho12,rhocp12
    logical , intent(in ), dimension(0:1,3) :: is_bound_phys
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)       :: vel,u,v,w
    real(rp), intent(in ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    integer , intent(in ) :: iq_lo,iq_hi
    real(rp), intent(inout), dimension(1-nh_mom:,1-nh_mom:,1-nh_mom:) :: sbq
    real(rp), intent(in   ), dimension(0:,0:,0:) :: asq
    real(rp), intent(in ), dimension(1-nh_s:,1-nh_s:,1-nh_s:), optional :: s
    real(rp) :: c2,dc,c2m,dcm,c2s,dcs,dl_d,dli_d,w1,w2,w1d,rext,ext_d,mf,mass,dil1,dil1m,dil1s,dil2,rh1,rh2, &
                rho_o,rho_n,qm,q0,qp,qa,qn
    integer  :: i,j,k,iq,di,dj,dk,qi,qj,qk,off,id,jd,kd,dd,nd,n1,n2,n3
    logical  :: is_stag,is_lo,is_hi,is_zero,is_lim,is_keep,is_keep_x,is_keep_y,is_keep_z
    !
    di = merge(1,0,idir == 1)
    dj = merge(1,0,idir == 2)
    dk = merge(1,0,idir == 3)
    n1 = n(1); n2 = n(2); n3 = n(3)
    nd = n(idir)
    is_keep_x = is_bound_phys(1,1)
    is_keep_y = is_bound_phys(1,2)
    is_keep_z = is_bound_phys(1,3)
    is_lo = is_bound_phys(0,idir)
    is_hi = is_bound_phys(1,idir)
    dl_d = dl(idir)
    dli_d = 1._rp/dl_d
    c2m = rho12(2)
    dcm = rho12(1)-rho12(2)
    c2s = rhocp12(2)
    dcs = rhocp12(1)-rhocp12(2)
    !
    ! fluxes of q (per unit area, over the sweep) through the lower faces of the control volumes 1:n+1 along
    ! idir: mass flux times transported value, from the upwind control volume; each face is computed once
    !
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(iq,qi,qj,qk,is_stag,c2,dc,w1,w2,w1d,ext_d,mf,mass,qm,q0,qp,qa,off,id,jd,kd,dd,is_zero,is_lim) &
    !$acc async(1)
    do k=1,n(3)+dk
      do j=1,n(2)+dj
        do i=1,n(1)+di
          !$acc loop seq
          do iq=iq_lo,iq_hi
            qi = merge(1,0,iq == 1)
            qj = merge(1,0,iq == 2)
            qk = merge(1,0,iq == 3)
            is_stag = iq == idir
            c2 = c2m
            dc = dcm
            if(iq == 0) then
              c2 = c2s
              dc = dcs
            end if
            w1 = 0.5_rp
            if(iq == 0) w1 = 1._rp
            if(iq == 3) w1 = wzc(k)
            w2 = 1._rp-w1
            if(is_stag) then
              mf = 0.5_rp*((c2*vel(i-di,j-dj,k-dk)*dt+dc*fl(i-di,j-dj,k-dk)) + &
                           (c2*vel(i   ,j   ,k   )*dt+dc*fl(i   ,j   ,k   )))
            else
              mf = w1*(c2*vel(i-di   ,j-dj   ,k-dk   )*dt+dc*fl(i-di   ,j-dj   ,k-dk   )) + &
                   w2*(c2*vel(i-di+qi,j-dj+qj,k-dk+qk)*dt+dc*fl(i-di+qi,j-dj+qj,k-dk+qk))
            end if
            off = -merge(1,0,mf >= 0._rp)
            id = i+off*di
            jd = j+off*dj
            kd = k+off*dk
            dd = id*di+jd*dj+kd*dk
            is_zero = (is_lo .and. dd <= 1) .or. (is_hi .and. dd >= nd)
            is_lim  = band(id,jd,kd) > 0_int8 .or. band(id+qi,jd+qj,kd+qk) > 0_int8
            w1d = 0.5_rp
            if(iq == 0) w1d = 1._rp
            if(iq == 3) w1d = wzc(kd)
            ext_d = dl_d
            if(idir == 3) then
              ext_d = dzf(kd)
              if(is_stag) ext_d = 0.5_rp*(dzf(kd)+dzf(kd+1))
            end if
            mass = (c2+dc*(w1d*f_prev(id,jd,kd)+(1._rp-w1d)*f_prev(id+qi,jd+qj,kd+qk)))*ext_d
            select case(iq)
            case(1)
              qm = ub(id-di,jd-dj,kd-dk); q0 = ub(id,jd,kd); qp = ub(id+di,jd+dj,kd+dk); qa = au(id,jd,kd)
            case(2)
              qm = vb(id-di,jd-dj,kd-dk); q0 = vb(id,jd,kd); qp = vb(id+di,jd+dj,kd+dk); qa = av(id,jd,kd)
            case(3)
              qm = wb(id-di,jd-dj,kd-dk); q0 = wb(id,jd,kd); qp = wb(id+di,jd+dj,kd+dk); qa = aw(id,jd,kd)
            case default
              qm = sbq(id-di,jd-dj,kd-dk); q0 = sbq(id,jd,kd); qp = sbq(id+di,jd+dj,kd+dk); qa = asq(id,jd,kd)
            end select
            qflx(i,j,k,iq) = mf*donor_value(qm,q0,qp,mf,mass,tc*qa,is_zero,is_lim)
          end do
        end do
      end do
    end do
    !
    ! update of the control volumes: new momentum (with the dilatation terms of the two cells), and new mass
    !
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(iq,qi,qj,qk,is_stag,is_keep,c2,dc,w1,w2,rext,dil1,dil1m,dil1s,dil2,rh1,rh2,rho_o,rho_n,qn) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          rh1 = dli_d
          if(idir == 3) rh1 = dzfi(k)
          dil1m = (c2m+dcm*cc(i,j,k))*(vel(i,j,k)-vel(i-di,j-dj,k-dk))*(dt*rh1)
          dil1s = 0._rp
          if(iq_lo == 0) dil1s = (c2s+dcs*cc(i,j,k))*(vel(i,j,k)-vel(i-di,j-dj,k-dk))*(dt*rh1)
          !$acc loop seq
          do iq=iq_lo,iq_hi
            qi = merge(1,0,iq == 1)
            qj = merge(1,0,iq == 2)
            qk = merge(1,0,iq == 3)
            is_stag = iq == idir
            c2 = c2m
            dc = dcm
            dil1 = dil1m
            if(iq == 0) then
              c2 = c2s
              dc = dcs
              dil1 = dil1s
            end if
            w1 = 0.5_rp
            if(iq == 0) w1 = 1._rp
            if(iq == 3) w1 = wzc(k)
            w2 = 1._rp-w1
            rext = dli_d
            rh2  = dli_d
            if(idir == 3) then
              rext = dzfi(k)
              if(is_stag) rext = dzci(k)
              rh2 = dzfi(k+qk)
            end if
            dil2 = (c2+dc*cc(i+qi,j+qj,k+qk))*(vel(i+qi,j+qj,k+qk)-vel(i-di+qi,j-dj+qj,k-dk+qk))*(dt*rh2)
            rho_o = c2+dc*(w1*f_prev(i,j,k)+w2*f_prev(i+qi,j+qj,k+qk))
            rho_n = c2+dc*(w1*psi(   i,j,k)+w2*psi(   i+qi,j+qj,k+qk))
            is_keep = (iq == 1 .and. is_keep_x .and. i == n1) .or. &
                      (iq == 2 .and. is_keep_y .and. j == n2) .or. &
                      (iq == 3 .and. is_keep_z .and. k == n3)
            if(.not.is_keep) then
              select case(iq)
              case(1)
                ub(i,j,k) = (rho_o*ub(i,j,k)-(qflx(i+di,j+dj,k+dk,iq)-qflx(i,j,k,iq))*rext+u(i,j,k)*(w1*dil1+w2*dil2))/rho_n
              case(2)
                vb(i,j,k) = (rho_o*vb(i,j,k)-(qflx(i+di,j+dj,k+dk,iq)-qflx(i,j,k,iq))*rext+v(i,j,k)*(w1*dil1+w2*dil2))/rho_n
              case(3)
                wb(i,j,k) = (rho_o*wb(i,j,k)-(qflx(i+di,j+dj,k+dk,iq)-qflx(i,j,k,iq))*rext+w(i,j,k)*(w1*dil1+w2*dil2))/rho_n
              case default
                qn = s(i,j,k)
                sbq(i,j,k) = (rho_o*sbq(i,j,k)-(qflx(i+di,j+dj,k+dk,iq)-qflx(i,j,k,iq))*rext+qn*(w1*dil1+w2*dil2))/rho_n
              end select
            end if
          end do
        end do
      end do
    end do
  end subroutine transport_sweep
  !
  pure real(rp) function donor_value(qm,q0,qp,mf,mass,dq,is_zero,is_lim) result(qtr)
    !
    ! transported value of q through a face: average over the donor region (the part of the upwind control
    ! volume, with value q0 and neighbours qm and qp along the sweep direction, holding the mass mf that
    ! crosses the face; mf is positive if the face is the upper face of the control volume, and mass is
    ! the mass of the control volume) of a linear reconstruction, centred in time with the non-advective
    ! change dq of q over the first half of the (sub-)step; the slope is limited (lim_band) if is_lim,
    ! unlimited (slope_bulk) otherwise, and zero if is_zero
    !
    ! in the interface band (is_lim), the time correction is scaled by 1-sig**2, such that the value that
    ! remains in the control volume, q0-sig*(slope correction)-sig*(1+sig)*dq, stays bounded also when (almost)
    ! all the mass leaves, e.g., when a heavy phase leaves a control volume that is left with a light phase;
    ! the scaling deviates from one by O(sig**2) only, which preserves the second-order accuracy in time
    !
    implicit none
    real(rp), intent(in) :: qm,q0,qp,mf,mass,dq
    logical , intent(in) :: is_zero,is_lim
    real(rp) :: a,b,sl,sig
    !$acc routine seq
    sl = 0._rp
    if(.not.is_zero) then
      a = q0-qm
      b = qp-q0
      if(is_lim) then
        if(a*b > 0._rp) then
          select case(lim_band)
          case(1) ! van Leer
            sl = 2._rp*a*b/(a+b)
          case(2) ! monotonized central
            sl = sign(min(2._rp*abs(a),2._rp*abs(b),0.5_rp*abs(a+b)),a)
          case default ! superbee
            sl = sign(max(min(2._rp*abs(a),abs(b)),min(abs(a),2._rp*abs(b))),a)
          end select
        end if
      else
        select case(slope_bulk)
        case(1) ! Fromm
          sl = 0.5_rp*(a+b)
        case default ! Lax-Wendroff (downwind slope)
          sl = merge(b,a,mf >= 0._rp)
        end select
      end if
    end if
    sig = min(abs(mf)/max(mass,tiny(1._rp)),1._rp)
    qtr = q0+merge(0.5_rp,-0.5_rp,mf >= 0._rp)*sl*(1._rp-sig)
    if(is_lim) then
      qtr = qtr+dq*(1._rp-sig**2)
    else
      qtr = qtr+dq
    end if
  end function donor_value
  !
  subroutine vof_split_contact_angle(n,nh_psi,nh_g,dl,dzc,dzf,cbcpsi,bcpsi,nb,is_bound,psi)
    !
    ! contact angles at walls (approach of S. Afkhami & M. Bussmann, Int. J. Numer. Methods Fluids 57 (2008)
    ! 453--472 and 61 (2009) 827--847, through the volume fractions of the ghost cells): the ghost cells of a
    ! wall with contact angle theta (measured through phase 1) take the volume fractions of the interface
    ! extended into the wall at that angle; in each wall cell with an interface, the plane has the volume
    ! fraction of the cell and the normal n = sin(theta)*t + cos(theta)*e (pointing out of phase 1), with t the
    ! unit wall-tangential part of the MYC normal (computed with the mirrored ghost cells of the Neumann
    ! condition) and e the unit normal to the wall pointing into the domain; each ghost cell takes the plane of
    ! the nearest wall cell with an interface along the wall (within nr_ca cells; otherwise, the mirrored value
    ! is kept); the normals of the cells next to the wall (interface reconstruction) and the height functions
    ! (curvature) then see the contact angle; the ghost cells do not enter the advection (no flux through a
    ! wall), such that the volume is conserved
    !
    ! to be called after the boundary conditions of the volume fraction are imposed; the halo cells of the
    ! ghost layers along the wall are then updated without reimposing the Neumann condition at the walls
    ! with a contact angle (boundary type 'C', which `boundp` leaves as it is); at an edge with another wall,
    ! the ghost cells are mirrored about it
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_psi,nh_g
    real(rp), intent(in   ), dimension(3) :: dl
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzc,dzf
    character(len=1), intent(in), dimension(0:1,3) :: cbcpsi
    real(rp), intent(in   ), dimension(0:1,3) :: bcpsi
    integer , intent(in   ), dimension(0:1,3) :: nb
    logical , intent(in   ), dimension(0:1,3) :: is_bound
    real(rp), intent(inout), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), dimension(-1:1,-1:1,-1:1) :: c3
    real(rp), dimension(-1:1) :: dz3
    real(rp), dimension(3) :: m
    logical , dimension(3) :: act
    character(len=1), dimension(0:1,3) :: cbc_ca
    real(rp) :: dx,dy,dz,cs,sn,fw,tn,d2,d2b,n1,n2,n3,o3,s3
    integer  :: id,is,da,db,iw,se,ia,ib,ip,iq,l,ic,jc,kc,i,j,k,ipb,iqb,kk,na,nb_,gi,gj,gk,pi,pj,pk,ig
    logical  :: is_act1,is_act2,is_act3
    !
    if(.not.is_ca_any) return
    dx = dl(1)
    dy = dl(2)
    is_act1 = is_active(1)
    is_act2 = is_active(2)
    is_act3 = is_active(3)
    do id=1,3
      do is=0,1
        if(.not.(is_ca(is,id).and.is_bound(is,id))) cycle
        da = merge(2,1,id == 1) ! directions along the wall (da < db)
        db = merge(2,3,id == 3)
        na = n(da)
        nb_ = n(db)
        iw = merge(1,n(id),is == 0) ! index of the wall cells
        se = merge(1,-1,is == 0)    ! orientation of the normal to the wall pointing into the domain
        cs = ca_cos(is,id)
        sn = ca_sin(is,id)
        !
        ! planes of the wall cells with an interface
        !
        !$acc parallel loop gang vector collapse(2) default(present) private(c3,dz3,m,act) async(1)
        do ib=2-nh_psi,nb_+nh_psi-1
          do ia=2-nh_psi,na+nh_psi-1
            i = merge(iw,ia,id == 1)
            j = merge(ia,merge(iw,ib,id == 2),id == 1)
            k = merge(iw,ib,id == 3)
            fw = psi(i,j,k)
            ca_pl(0,ia,ib) = 0._rp
            ca_pl(1,ia,ib) = 0._rp
            ca_pl(2,ia,ib) = 0._rp
            ca_pl(3,ia,ib) = 0._rp
            if(fw > f_eps .and. fw < 1._rp-f_eps) then
              act(:) = [is_act1,is_act2,is_act3]
              do kc=-1,1
                do jc=-1,1
                  do ic=-1,1
                    c3(ic,jc,kc) = psi(i+ic,j+jc,k+kc)
                  end do
                end do
                dz3(kc) = dzf(k+kc)
              end do
              call plic_normal_myc(c3,dx,dy,dz3,act,m(1),m(2),m(3))
              dz = dz3(0)
              n1 = merge(0._rp,m(1)/dx,id == 1) ! wall-tangential part of the unit normal
              n2 = merge(0._rp,m(2)/dy,id == 2)
              n3 = merge(0._rp,m(3)/dz,id == 3)
              tn = sqrt(n1**2+n2**2+n3**2)
              if(tn > 0._rp) then
                n1 = merge(cs*se,sn*n1/tn,id == 1)
                n2 = merge(cs*se,sn*n2/tn,id == 2)
                n3 = merge(cs*se,sn*n3/tn,id == 3)
                ca_pl(1,ia,ib) = n1*dx
                ca_pl(2,ia,ib) = n2*dy
                ca_pl(3,ia,ib) = n3*dz
                ca_pl(0,ia,ib) = plic_alpha(n1*dx,n2*dy,n3*dz,fw)
              end if
            end if
          end do
        end do
        !
        ! ghost cells: volume fraction of the plane of the nearest wall cell with an interface
        !
        !$acc parallel loop collapse(3) default(present) async(1)
        do l=1,nh_psi
          do ib=1,nb_
            do ia=1,na
              ipb = 0; iqb = 0
              d2b = huge(1._rp)
              do iq=-nr_ca,nr_ca
                do ip=-nr_ca,nr_ca
                  d2 = real(ip**2+iq**2,rp)
                  if(d2 < d2b .and. abs(ca_pl(1,ia+ip,ib+iq))+abs(ca_pl(2,ia+ip,ib+iq))+abs(ca_pl(3,ia+ip,ib+iq)) > 0._rp) then
                    d2b = d2; ipb = ip; iqb = iq
                  end if
                end do
              end do
              if(d2b < huge(1._rp)) then
                !
                ! ghost cell (gi,gj,gk), and wall cell (pi,pj,pk) with the plane
                !
                ig = iw-se*l
                gi = merge(ig,ia,id == 1)
                gj = merge(ia,merge(ig,ib,id == 2),id == 1)
                gk = merge(ig,ib,id == 3)
                pi = merge(iw,ia+ipb,id == 1)
                pj = merge(ia+ipb,merge(iw,ib+iqb,id == 2),id == 1)
                pk = merge(iw,ib+iqb,id == 3)
                !
                ! box of the ghost cell in the cell-scaled coordinates of the wall cell (non-uniform along z)
                !
                o3 = 0._rp
                do kk=min(gk,pk),max(gk,pk)-1
                  o3 = o3+dzf(kk)
                end do
                o3 = sign(o3,real(gk-pk,rp))/dzf(pk)
                s3 = dzf(gk)/dzf(pk)
                psi(gi,gj,gk) = plic_volume_box(ca_pl(1,ia+ipb,ib+iqb),ca_pl(2,ia+ipb,ib+iqb),ca_pl(3,ia+ipb,ib+iqb), &
                                                ca_pl(0,ia+ipb,ib+iqb),real(gi-pi,rp),real(gj-pj,rp),o3,1._rp,1._rp,s3)
              end if
            end do
          end do
        end do
      end do
    end do
    !
    ! halo cells of the ghost layers along the walls (on all the tasks)
    !
    cbc_ca(:,:) = cbcpsi(:,:)
    where(is_ca(:,:)) cbc_ca(:,:) = 'C'
    call boundp(cbc_ca,n,nh_psi,nh_g,bcpsi,nb,is_bound,dl,dzc,psi)
  end subroutine vof_split_contact_angle
  !
  pure subroutine donor_fractions(mx,my,mz,vol,idir,is_thinc_rec,beta,cfl_p,cfl_m,frac_p,frac_m)
    !
    ! fractions of the donor slabs of an interfacial cell occupied by phase 1, for the slab [1-cfl_p,1]
    ! adjacent to its upper face and the slab [0,cfl_m] adjacent to its lower face along direction idir
    ! (in cell-scaled coordinates), for an interface with normal (mx,my,mz) (cell-scaled, pointing out of
    ! phase 1) and volume fraction vol:
    !  - geometric reconstruction: exact slab volumes cut by the plane;
    !  - THINC reconstruction (`is_thinc_rec`): slab averages of the hyperbolic-tangent profile along idir
    !    with sharpness beta, through the cell centre (the transverse integrals are approximated by the
    !    midpoint rule, as in the THINC/QQ method), whose cell average is vol; the slab fractions lie in
    !    [0,1] and are bounded by the content of the cell, as required by the bounded split advection
    !
    implicit none
    real(rp), intent(in ) :: mx,my,mz,vol
    integer , intent(in ) :: idir
    logical , intent(in ) :: is_thinc_rec
    real(rp), intent(in ) :: beta
    real(rp), intent(in ) :: cfl_p,cfl_m
    real(rp), intent(out) :: frac_p,frac_m
    real(rp), parameter :: a_min = 1.e-12_rp
    real(rp), dimension(3) :: o,s
    real(rp) :: alpha,a,d
    !$acc routine seq
    if(.not.is_thinc_rec) then
      alpha = plic_alpha(mx,my,mz,vol)
      frac_p = 0._rp
      frac_m = 0._rp
      if(cfl_p > 0._rp) then
        o(:) = 0._rp; s(:) = 1._rp
        o(idir) = 1._rp-cfl_p; s(idir) = cfl_p
        frac_p = plic_volume_box(mx,my,mz,alpha,o(1),o(2),o(3),s(1),s(2),s(3))
      end if
      if(cfl_m > 0._rp) then
        o(:) = 0._rp; s(:) = 1._rp
        s(idir) = cfl_m
        frac_m = plic_volume_box(mx,my,mz,alpha,o(1),o(2),o(3),s(1),s(2),s(3))
      end if
    else
      a = merge(mx,merge(my,mz,idir == 2),idir == 1)/max(sqrt(mx**2+my**2+mz**2),tiny(1._rp))
      frac_p = vol
      frac_m = vol
      if(abs(a) > a_min) then
        d = thinc_intercept(abs(a),vol,beta)
        if(a > 0._rp) then
          if(cfl_p > 0._rp) frac_p = thinc_slab(a,d,beta,0.5_rp-cfl_p,0.5_rp)
          if(cfl_m > 0._rp) frac_m = thinc_slab(a,d,beta,-0.5_rp,-0.5_rp+cfl_m)
        else ! mirrored profile, s -> -s
          if(cfl_p > 0._rp) frac_p = thinc_slab(-a,d,beta,-0.5_rp,-0.5_rp+cfl_p)
          if(cfl_m > 0._rp) frac_m = thinc_slab(-a,d,beta,0.5_rp-cfl_m,0.5_rp)
        end if
      end if
    end if
  end subroutine donor_fractions
  !
  subroutine vof_split_curvature(istage,n,nh_psi,nh_ifc,nh_g,dl,dzf,psi,kappa)
    !
    ! interface curvature (sign convention of Fizzy: kappa = -div(n), n pointing into phase 1), in stages
    ! separated by halo exchanges of kappa (without Dirichlet conditions, see `ifc_geometry`):
    !   1. height functions (nine-cell columns) in interfacial cells, along the directions in decreasing order
    !      of the MYC normal components
    !   2. interfacial cells without valid heights: average of the height-function curvatures of the
    !      neighbouring interfacial cells
    !   3. full or empty cells next to interfacial cells: average of the curvatures of the latter (used by the
    !      face interpolation of the surface tension force); zero elsewhere
    !
    implicit none
    integer , intent(in   ) :: istage
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_psi,nh_ifc,nh_g
    real(rp), intent(in   ), dimension(3) :: dl
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzf
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: kappa
    real(rp) :: fc,ksum,kcnt
    integer  :: i,j,k,ii,jj,kk,q,il,lin,n1,n2
    !
    select case(istage)
    case(1)
      !
      ! (the interfacial cells are compacted into a list, processed by a second kernel)
      !
      nlist(1) = 0
      !$acc update device(nlist) async(1)
      !$acc serial default(present) async(1)
      kappa_stats(1) = 0._rp
      !$acc end serial
      !$acc parallel loop collapse(3) default(present) private(il) async(1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            kappa(i,j,k) = 0._rp
            fc = psi(i,j,k)
            if(fc > f_eps .and. fc < 1._rp-f_eps) then
              !$acc atomic capture
              nlist(1) = nlist(1)+1
              il = nlist(1)
              !$acc end atomic
              ilist(il) = (i-1)+n(1)*((j-1)+n(2)*(k-1))
            end if
          end do
        end do
      end do
      call hf_list(1,1,1,n(1),n(2),nh_psi,nh_ifc,nh_g,dl,dzf,psi,kappa)
    case(2)
      !
      ! (only the listed interfacial cells are visited; the few of them without height-function curvature
      ! are counted with atomics, on the device, and take their new value once all are computed)
      !
      n1 = n(1)
      n2 = n(2)
      !$acc serial default(present) async(1)
      kappa_stats(2) = 0._rp
      kappa_stats(3) = 0._rp
      !$acc end serial
      !$acc parallel loop default(present) async(1)
      do il=1,nlist(1)
        lin = ilist(il)
        i = 1+mod(lin,n1)
        j = 1+mod(lin/n1,n2)
        k = 1+lin/(n1*n2)
        if(kappa(i,j,k) == kappa_undef) then
          ksum = 0._rp
          kcnt = 0._rp
          do kk=-1,1
            do jj=-1,1
              do ii=-1,1
                fc = psi(i+ii,j+jj,k+kk)
                if(fc > f_eps .and. fc < 1._rp-f_eps .and. kappa(i+ii,j+jj,k+kk) /= kappa_undef) then
                  ksum = ksum + kappa(i+ii,j+jj,k+kk)
                  kcnt = kcnt + 1._rp
                end if
              end do
            end do
          end do
          if(kcnt > 0._rp) then
            kap(i,j,k) = ksum/kcnt
            !$acc atomic update
            kappa_stats(2) = kappa_stats(2) + 1._rp
          else
            kap(i,j,k) = 0._rp
            !$acc atomic update
            kappa_stats(3) = kappa_stats(3) + 1._rp
          end if
        else
          kap(i,j,k) = kappa(i,j,k)
        end if
      end do
      !$acc parallel loop default(present) async(1)
      do il=1,nlist(1)
        lin = ilist(il)
        i = 1+mod(lin,n1)
        j = 1+mod(lin/n1,n2)
        k = 1+lin/(n1*n2)
        kappa(i,j,k) = kap(i,j,k)
      end do
    case(3)
      !
      ! (in place: only the curvature of interfacial cells is read, and only that of the other cells, zero
      ! since the first stage, is written)
      !
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            fc = psi(i,j,k)
            if(.not.(fc > f_eps .and. fc < 1._rp-f_eps)) then
              ksum = 0._rp
              kcnt = 0._rp
              do q=1,6
                ii = merge(merge(-1,1,q == 1),0,q <= 2)
                jj = merge(merge(-1,1,q == 3),0,q == 3 .or. q == 4)
                kk = merge(merge(-1,1,q == 5),0,q >= 5)
                fc = psi(i+ii,j+jj,k+kk)
                if(fc > f_eps .and. fc < 1._rp-f_eps) then
                  ksum = ksum + kappa(i+ii,j+jj,k+kk)
                  kcnt = kcnt + 1._rp
                end if
              end do
              if(kcnt > 0._rp) kappa(i,j,k) = ksum/kcnt
            end if
          end do
        end do
      end do
    end select
  end subroutine vof_split_curvature
  !
  subroutine hf_list(ilo,jlo,klo,ni,nj,nh_psi,nh_ifc,nh_g,dl,dzf,psi,kappa,mx,my,mz)
    !
    ! height-function geometry of the cells of the compacted list ilist (packed indices, with offsets ilo,
    ! jlo, klo and dimensions ni, nj): MYC normal, and height functions (nine-cell columns) along the
    ! directions in decreasing order of the MYC normal components (cell-scaled), until valid ones are found;
    ! with `kappa`, the curvature is stored in it (kappa_undef where the heights are not valid), and the cells
    ! with valid heights are counted in kappa_stats(1); without it, the normal in cell-scaled coordinates
    ! pointing out of phase 1 (|mx|+|my|+|mz| = 1), for the geometric reconstruction, is stored in mx, my
    ! and mz: that of the height functions if valid, and the MYC one otherwise; the number of listed cells is
    ! read from the device copy of nlist, and the count is kept on the device, without host synchronizations
    !
    ! (a single kernel: calling a device routine per cell, with the volume fraction as argument, is several
    ! times slower on GPUs)
    !
    implicit none
    integer , intent(in   )               :: ilo,jlo,klo,ni,nj
    integer , intent(in   )               :: nh_psi,nh_ifc,nh_g
    real(rp), intent(in   ), dimension(3) :: dl
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzf
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:), optional :: kappa
    real(rp), intent(inout), dimension(0:,0:,0:), optional :: mx,my,mz
    real(rp), dimension(-1:1,-1:1,-1:1) :: c3
    real(rp), dimension(-1:1) :: dz3
    real(rp), dimension(-1:1,-1:1) :: hgt
    real(rp), dimension(3) :: dm,dp,m,h,nrm
    logical , dimension(3) :: act
    logical  :: is_valid,is_act1,is_act2,is_act3,is_curv
    real(rp) :: dx,dy,kap_c,sgn,lo,hi,dqq,vol,fq,f0,f1,na,nb,nq,s
    integer  :: i,j,k,q,il,lin,a,b,r,d,d1,d2,d3,ia,ib,itmp,ic,jc,kc
    !
    is_curv = present(kappa)
    dx = dl(1)
    dy = dl(2)
    is_act1 = is_active(1)
    is_act2 = is_active(2)
    is_act3 = is_active(3)
    !$acc parallel loop default(present) private(m,act,c3,dz3,hgt,dm,dp,h,nrm) async(1)
    do il=1,nlist(1)
      lin = ilist(il)
      i = ilo+mod(lin,ni)
      j = jlo+mod(lin/ni,nj)
      k = klo+lin/(ni*nj)
      act(:) = [is_act1,is_act2,is_act3]
      do kc=-1,1
        do jc=-1,1
          do ic=-1,1
            c3(ic,jc,kc) = psi(i+ic,j+jc,k+kc)
          end do
        end do
        dz3(kc) = dzf(k+kc)
      end do
      call plic_normal_myc(c3,dx,dy,dz3,act,m(1),m(2),m(3))
      kap_c = 0._rp
      is_valid = .false.
      !
      ! directions in decreasing order of the MYC normal components (scalars: a small index array updated in
      ! place and then used in subscripts is miscompiled in device code by nvfortran)
      !
      d1 = 1; d2 = 2; d3 = 3
      if(abs(m(d1)) < abs(m(d2))) then
        itmp = d1; d1 = d2; d2 = itmp
      end if
      if(abs(m(d2)) < abs(m(d3))) then
        itmp = d2; d2 = d3; d3 = itmp
      end if
      if(abs(m(d1)) < abs(m(d2))) then
        itmp = d1; d1 = d2; d2 = itmp
      end if
      dm(:) = [dx,dy,0.5_rp*(dzf(k-1)+dzf(k))]
      dp(:) = [dx,dy,0.5_rp*(dzf(k)+dzf(k+1))]
      do r=1,3
        d = merge(d1,merge(d2,d3,r == 2),r == 1)
        if(.not.act(d)) cycle
        ia = merge(2,1,d == 1)
        ib = merge(2,3,d == 3)
        lo = 0.5_rp*merge(dzf(k),merge(dx,dy,d == 1),d == 3)
        hi = lo
        do q=1,nq_hf
          lo = lo + merge(dzf(k-q),merge(dx,dy,d == 1),d == 3)
          hi = hi + merge(dzf(k+q),merge(dx,dy,d == 1),d == 3)
        end do
        f0 = psi(i+merge(-nq_hf,0,d == 1),j+merge(-nq_hf,0,d == 2),k+merge(-nq_hf,0,d == 3))
        f1 = psi(i+merge( nq_hf,0,d == 1),j+merge( nq_hf,0,d == 2),k+merge( nq_hf,0,d == 3))
        sgn = 0._rp
        if(     f0 >= 1._rp-eps_hf .and. f1 <= eps_hf) then
          sgn =  1._rp
        else if(f0 <= eps_hf .and. f1 >= 1._rp-eps_hf) then
          sgn = -1._rp
        end if
        is_valid = sgn /= 0._rp
        hgt(:,:) = 0._rp
        do b=-1,1
          do a=-1,1
            vol = 0._rp
            do q=-nq_hf,nq_hf
              fq  = psi(i+merge(q,a,d == 1),j+merge(q,merge(a,b,d == 1),d == 2),k+merge(q,b,d == 3))
              dqq = merge(dzf(k+q),merge(dx,dy,d == 1),d == 3)
              if(q == -nq_hf) f0 = fq
              if(q ==  nq_hf) f1 = fq
              vol = vol + fq*dqq
            end do
            if(sgn > 0._rp) then
              is_valid = is_valid .and. f0 >= 1._rp-eps_hf .and. f1 <= eps_hf
            else
              is_valid = is_valid .and. f0 <= eps_hf .and. f1 >= 1._rp-eps_hf
            end if
            hgt(a,b) = merge(vol-lo,hi-vol,sgn > 0._rp)
          end do
        end do
        if(is_valid) then
          if(is_curv) then
            kap_c = hf_curvature(hgt,sgn,dm(ia),dp(ia),dm(ib),dp(ib),act(ia),act(ib))
          else
            call hf_normal(hgt,sgn,dm(ia),dp(ia),dm(ib),dp(ib),act(ia),act(ib),na,nb,nq)
            nrm(ia) = na; nrm(ib) = nb; nrm(d) = nq
          end if
          exit
        end if
      end do
      if(is_curv) then
        if(is_valid) then
          kappa(i,j,k) = kap_c
          !$acc atomic update
          kappa_stats(1) = kappa_stats(1) + 1._rp
        else
          kappa(i,j,k) = kappa_undef
        end if
      else
        if(is_valid) then
          !
          ! (norm stored in a scalar: with m(1:3) on both sides, nvfortran allocates a temporary on the device heap)
          !
          h(:) = [dx,dy,dz3(0)]
          m(:) = nrm(:)*h(:)
          s = abs(m(1))+abs(m(2))+abs(m(3))
          m(:) = m(:)/s
        end if
        mx(i,j,k) = m(1)
        my(i,j,k) = m(2)
        mz(i,j,k) = m(3)
      end if
    end do
  end subroutine hf_list
  !
  subroutine set_zero(p)
    implicit none
    real(rp), intent(inout), dimension(:,:,:) :: p
    integer :: i,j,k
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,size(p,3)
      do j=1,size(p,2)
        do i=1,size(p,1)
          p(i,j,k) = 0._rp
        end do
      end do
    end do
  end subroutine set_zero
  !
  subroutine vof_split_diagnostics(n,nh_psi,nh_g,dl,dzf,psi,vol,clip,fmin,fmax,cfl,nsub,kstats)
    !
    ! total volume of phase 1, accumulated clipped volume, extrema of the volume fraction before
    ! clipping in the last step, CFL numbers and number of sub-steps of the last step, and numbers of
    ! interfacial cells with curvature from height functions, from neighbouring ones, and without estimate
    !
    implicit none
    integer , intent(in ), dimension(3) :: n
    integer , intent(in )               :: nh_psi,nh_g
    real(rp), intent(in ), dimension(3) :: dl
    real(rp), intent(in ), dimension(1-nh_g:) :: dzf
    real(rp), intent(in ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(out) :: vol,clip,fmin,fmax
    real(rp), intent(out), dimension(3) :: cfl
    integer , intent(out) :: nsub
    real(rp), intent(out), dimension(3) :: kstats
    integer :: i,j,k
    !
    vol = 0._rp
    !$acc parallel loop collapse(3) default(present) reduction(+:vol) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          vol = vol+psi(i,j,k)*dl(1)*dl(2)*dzf(k)
        end do
      end do
    end do
    !$acc update self(clip_stats,kappa_stats) async(1)
    !$acc wait(1)
    call MPI_ALLREDUCE(MPI_IN_PLACE,vol,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(clip_stats(1),clip,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(clip_stats(2),fmin,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(clip_stats(3),fmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(kappa_stats,kstats,3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    cfl(:) = cfl_last(:)
    nsub = nsub_last
  end subroutine vof_split_diagnostics
end module mod_vof_split
