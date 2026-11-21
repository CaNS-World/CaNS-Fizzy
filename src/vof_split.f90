! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 The CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_vof_split
  !
  ! conservative direction-split advection of the volume fraction, following
  !   G.D. Weymouth & D.K.-P. Yue, J. Comput. Phys. 229 (2010) 2853--2865,
  ! with geometric donor-region fluxes from a PLIC reconstruction (see `mod_plic`)
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
  ! the volume fraction needs two halo layers: the interface is also reconstructed in the first one,
  ! whose cells are the donors of the fluxes entering the subdomain
  !
  use mpi
  use mod_types
  use mod_common_mpi, only: ierr
  use mod_bound     , only: boundp
  use mod_plic      , only: plic_alpha,plic_volume_box,plic_normal_myc
  implicit none
  private
  public vof_split_init,vof_split_advance,vof_split_diagnostics
  !
  real(rp), parameter :: cfl_max = 0.5_rp    ! maximum CFL number of a sweep
  real(rp), parameter :: f_eps   = 1.e-12_rp ! cells with f_eps < f < 1-f_eps are reconstructed
  !
  real(rp), allocatable, dimension(:,:,:) :: cc       ! dilatation indicator c_c
  real(rp), allocatable, dimension(:,:,:) :: fl       ! face fluxes of a sweep
  real(rp), allocatable, dimension(:,:,:) :: ua,va,wa ! advecting velocity, at the middle of the time step
  real(rp), allocatable, dimension(:,:,:) :: uo,vo,wo ! velocity at the previous time step
  integer , allocatable, dimension(:)     :: ilist    ! compacted list of interfacial donor cells (packed indices)
  integer , dimension(1) :: nlist                     ! number of entries of ilist
  logical , dimension(3) :: is_active                 ! false along the directions in which the problem is invariant
  logical :: is_first   = .true.
  logical :: is_forward = .true.                      ! sweep order, reversed at each (sub-)step
  !
  ! diagnostics (local to the task, accumulated on the device; reduced by `vof_split_diagnostics`)
  !
  real(rp) :: clip_vol = 0._rp                  ! accumulated volume added by clipping
  real(rp) :: fmin_raw = 0._rp,fmax_raw = 1._rp ! extrema before clipping in the last step
  real(rp) :: clip_s,fmin_s,fmax_s              ! the same, for one sweep
  !$acc declare create(clip_vol,fmin_raw,fmax_raw,clip_s,fmin_s,fmax_s)
  real(rp), dimension(3) :: cfl_last = 0._rp    ! CFL numbers of the advecting velocity in the last step
  integer  :: nsub_last = 0                     ! number of sub-steps in the last step
  contains
  subroutine vof_split_init(ng,n,nh_psi)
    !
    ! allocates the work arrays
    !
    implicit none
    integer , intent(in), dimension(3) :: ng,n
    integer , intent(in)               :: nh_psi
    integer :: i,j,k
    !
    if(nh_psi < 2) error stop 'ERROR: the split VoF needs two halo layers of the volume fraction.'
    !
    ! directions with one or two grid cells are taken as invariant (2D runs): no sweeps,
    ! and no contribution to the interface normal
    !
    is_active(:) = ng(:) > 2
    allocate(cc(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    allocate(fl,ua,va,wa,uo,vo,wo,mold=cc)
    allocate(ilist((n(1)+2)*(n(2)+2)*(n(3)+2)))
    nlist(:) = 0
    !$acc enter data create(cc,fl,ua,va,wa,uo,vo,wo,ilist) copyin(nlist)
    !$acc update device(clip_vol,fmin_raw,fmax_raw)
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
  subroutine vof_split_advance(tm_coeff,n,nh_v,nh_psi,nh_g,dl,dzc,dzf,dzfi,dt,cbcpsi,bcpsi,nb,is_bound,u,v,w, &
                               psi,psiflx_x,psiflx_y,psiflx_z)
    !
    ! advances the volume fraction over one time step, with the face velocity extrapolated to the
    ! middle of the time step from the current and previous velocity fields (both divergence free);
    ! the volume fraction is returned with its halo layers updated
    !
    implicit none
    real(rp), intent(in   ), dimension(2) :: tm_coeff
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v,nh_psi,nh_g
    real(rp), intent(in   ), dimension(3) :: dl
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzc,dzf,dzfi
    real(rp), intent(in   ) :: dt
    character(len=1), intent(in), dimension(0:1,3) :: cbcpsi
    real(rp), intent(in   ), dimension(0:1,3) :: bcpsi
    integer , intent(in   ), dimension(0:1,3) :: nb
    logical , intent(in   ), dimension(0:1,3) :: is_bound
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)       :: u,v,w
    real(rp), intent(inout), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psiflx_x,psiflx_y,psiflx_z
    real(rp) :: dts
    integer  :: i,j,k,isub,nsub,q,idir
    !
    ! advecting velocity, and CFL numbers along each direction
    !
    call advecting_velocity(tm_coeff,n,nh_v,nh_g,dl,dzfi,dt,u,v,w,cfl_last)
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
    if(.not.is_active(1)) call set_zero(psiflx_x)
    if(.not.is_active(2)) call set_zero(psiflx_y)
    if(.not.is_active(3)) call set_zero(psiflx_z)
    !
    !$acc serial async(1)
    fmin_raw =  huge(1._rp)
    fmax_raw = -huge(1._rp)
    !$acc end serial
    do isub=1,nsub
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            cc(i,j,k) = merge(1._rp,0._rp,psi(i,j,k) > 0.5_rp)
          end do
        end do
      end do
      do q=1,3
        idir = merge(q,4-q,is_forward)
        if(.not.is_active(idir)) cycle
        select case(idir)
        case(1)
          call sweep(1,n,nh_psi,nh_g,dl,dzf,dzfi,dts,ua,isub == 1,1._rp/dt,psiflx_x,psi)
        case(2)
          call sweep(2,n,nh_psi,nh_g,dl,dzf,dzfi,dts,va,isub == 1,1._rp/dt,psiflx_y,psi)
        case(3)
          call sweep(3,n,nh_psi,nh_g,dl,dzf,dzfi,dts,wa,isub == 1,1._rp/dt,psiflx_z,psi)
        end select
        call boundp(cbcpsi,n,nh_psi,nh_g,bcpsi,nb,is_bound,dl,dzc,psi)
      end do
      is_forward = .not.is_forward
    end do
  end subroutine vof_split_advance
  !
  subroutine advecting_velocity(tm_coeff,n,nh_v,nh_g,dl,dzfi,dt,u,v,w,cfl)
    !
    ! advecting velocity, extrapolated to the middle of the time step from the current and previous
    ! (divergence-free) velocity fields, and its (local) maximum CFL numbers along each direction;
    ! the current velocity is stored for the next time step
    !
    implicit none
    real(rp), intent(in ), dimension(2) :: tm_coeff
    integer , intent(in ), dimension(3) :: n
    integer , intent(in )               :: nh_v,nh_g
    real(rp), intent(in ), dimension(3) :: dl
    real(rp), intent(in ), dimension(1-nh_g:) :: dzfi
    real(rp), intent(in ) :: dt
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(out), dimension(3) :: cfl
    real(rp) :: c1,c2,dtx,dty,cflx,cfly,cflz,uaux,vaux,waux
    integer  :: i,j,k,n1,n2,n3
    !
    c1 = tm_coeff(1)
    c2 = tm_coeff(2)
    if(is_first) then
      c1 = 1._rp
      c2 = 0._rp
      is_first = .false.
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
  end subroutine advecting_velocity
  !
  subroutine sweep(idir,n,nh_psi,nh_g,dl,dzf,dzfi,dt,vel,is_init,scl,psiflx,psi)
    !
    ! one directional sweep: geometric donor-region fluxes, and conservative update with the dilatation term;
    ! the face fluxes, times scl, are stored in psiflx if is_init, and accumulated otherwise; the extrema of
    ! the volume fraction before clipping and the clipped volume are accumulated in fmin_raw, fmax_raw and
    ! clip_vol
    !
    implicit none
    integer , intent(in   ) :: idir
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_psi,nh_g
    real(rp), intent(in   ), dimension(3) :: dl
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzf,dzfi
    real(rp), intent(in   ) :: dt
    real(rp), intent(in   ), dimension(0:,0:,0:) :: vel
    logical , intent(in   ) :: is_init
    real(rp), intent(in   ) :: scl
    real(rp), intent(inout), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psiflx,psi
    real(rp), dimension(-1:1,-1:1,-1:1) :: c
    real(rp), dimension(-1:1) :: dz3
    logical , dimension(3) :: act
    logical  :: is_own_p,is_own_m,is_act1,is_act2,is_act3
    real(rp) :: dx,dy,dxi,dyi,h,rh,fc,vp,vm,cfl_p,cfl_m,frac_p,frac_m,fnew,fcl,mx,my,mz
    integer  :: i,j,k,di,dj,dk,ii,nd,ilo,jlo,klo,ihi,jhi,khi,ni,nj,il,lin
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
    ilo = 1-di; ihi = n(1)+di
    jlo = 1-dj; jhi = n(2)+dj
    klo = 1-dk; khi = n(3)+dk
    !
    ! fluxes: each face is computed by its donor (upwind) cell; the fluxes of the (few) interfacial donor
    ! cells, which reconstruct the interface, are computed by a second kernel over a compacted list of
    ! these cells, such that the work of the reconstruction is not spread over all the GPU threads
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
          if(ii <= nd) vp = vel(i,j,k)
          vm = 0._rp
          if(ii >= 1) vm = vel(i-di,j-dj,k-dk)
          is_own_p = ii <= nd .and. vp >= 0._rp
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
      if(ii <= nd) vp = vel(i,j,k)
      vm = 0._rp
      if(ii >= 1) vm = vel(i-di,j-dj,k-dk)
      is_own_p = ii <= nd .and. vp >= 0._rp
      is_own_m = ii >= 1  .and. vm <  0._rp
      fc = min(max(psi(i,j,k),0._rp),1._rp)
      cfl_p = 0._rp
      cfl_m = 0._rp
      if(is_own_p) cfl_p =  vp*dt/h
      if(is_own_m) cfl_m = -vm*dt/h
      act(:) = [is_act1,is_act2,is_act3]
      c(:,:,:) = psi(i-1:i+1,j-1:j+1,k-1:k+1)
      dz3(:) = dzf(k-1:k+1)
      call plic_normal_myc(c,dx,dy,dz3,act,mx,my,mz)
      call donor_fractions(mx,my,mz,fc,idir,cfl_p,cfl_m,frac_p,frac_m)
      if(is_own_p) fl(i,j,k) = frac_p*vp*dt
      if(is_own_m) fl(i-di,j-dj,k-dk) = frac_m*vm*dt
    end do
    !
    ! update, with the dilatation term, clipping of round-off errors, and accumulation of the fluxes
    !
    ! (the extrema start from 0 and 1, and only the few values out of these bounds, the ones clipped, update
    ! them and the clipped volume, with atomics instead of reductions over the whole domain)
    !
    !$acc serial async(1)
    fmin_s = 0._rp
    fmax_s = 1._rp
    clip_s = 0._rp
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
            fmin_s = min(fmin_s,fnew)
            !$acc atomic update
            fmax_s = max(fmax_s,fnew)
            !$acc atomic update
            clip_s = clip_s+(fcl-fnew)*dx*dy*dzf(k)
          end if
          psi(i,j,k) = fcl
          if(is_init) then
            psiflx(i,j,k) = fl(i,j,k)*scl
            if(ii == 1) psiflx(i-di,j-dj,k-dk) = fl(i-di,j-dj,k-dk)*scl
          else
            psiflx(i,j,k) = psiflx(i,j,k)+fl(i,j,k)*scl
            if(ii == 1) psiflx(i-di,j-dj,k-dk) = psiflx(i-di,j-dj,k-dk)+fl(i-di,j-dj,k-dk)*scl
          end if
        end do
      end do
    end do
    !$acc serial async(1)
    fmin_raw = min(fmin_raw,fmin_s)
    fmax_raw = max(fmax_raw,fmax_s)
    clip_vol = clip_vol+clip_s
    !$acc end serial
  end subroutine sweep
  !
  pure subroutine donor_fractions(mx,my,mz,vol,idir,cfl_p,cfl_m,frac_p,frac_m)
    !
    ! fractions of the donor slabs of an interfacial cell occupied by phase 1, for the slab [1-cfl_p,1]
    ! adjacent to its upper face and the slab [0,cfl_m] adjacent to its lower face along direction idir
    ! (in cell-scaled coordinates), for an interface with normal (mx,my,mz) (cell-scaled, pointing out of
    ! phase 1) and volume fraction vol: exact slab volumes cut by the plane
    !
    implicit none
    real(rp), intent(in ) :: mx,my,mz,vol
    integer , intent(in ) :: idir
    real(rp), intent(in ) :: cfl_p,cfl_m
    real(rp), intent(out) :: frac_p,frac_m
    real(rp), dimension(3) :: o,s
    real(rp) :: alpha
    !$acc routine seq
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
  end subroutine donor_fractions
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
  subroutine vof_split_diagnostics(n,nh_psi,nh_g,dl,dzf,psi,vol,clip,fmin,fmax,cfl,nsub)
    !
    ! total volume of phase 1, accumulated clipped volume, extrema of the volume fraction before
    ! clipping in the last step, CFL numbers and number of sub-steps of the last step
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
    !$acc update self(clip_vol,fmin_raw,fmax_raw) async(1)
    !$acc wait(1)
    call MPI_ALLREDUCE(MPI_IN_PLACE,vol,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(clip_vol,clip,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(fmin_raw,fmin,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(fmax_raw,fmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    cfl(:) = cfl_last(:)
    nsub = nsub_last
  end subroutine vof_split_diagnostics
end module mod_vof_split
