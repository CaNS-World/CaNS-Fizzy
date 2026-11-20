! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 The CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_solver_vc
#if !defined(_CONSTANT_COEFFS_POISSON)
  use mpi
  use mod_common_mpi, only: myid
  use mod_param     , only: hypre_solver_i,hypre_tol,hypre_maxiter, &
                            hypre_precond_i,hypre_precond_maxiter,hypre_gmres_k_dim, &
                            hypre_pfmg_relax_type,hypre_pfmg_num_pre_relax,hypre_pfmg_num_post_relax, &
                            hypre_pfmg_max_levels
  use mod_types
  use, intrinsic :: iso_c_binding, only: C_PTR,C_NULL_PTR
  implicit none
  integer :: ierr
  private
  public solver_vc,solver_vc_helmholtz,solver_vc_end
  integer, parameter :: HYPRESolverSMG       = 1, &
                        HYPRESolverPFMG      = 2, &
                        HYPRESolverGMRES     = 3, &
                        HYPRESolverBiCGSTAB  = 4, &
                        HYPRESolverFlexGMRES = 5
  integer, parameter :: HYPREPrecondSMG  = 0, &
                        HYPREPrecondPFMG = 1, &
                        HYPREPrecondNone = 9
  integer, parameter :: HYPREProblemPressure  = 1, & ! entries of the `hypre_*` inputs
                        HYPREProblemDiffusion = 2
  type hypre_solver
    type(C_PTR) :: grid,stencil,precond,solver,mat,rhs,sol
    integer     :: stype,ptype,comm_hypre
    integer     :: ndim       ! number of directions of the problem (those with more than one grid cell)
    integer     :: idirs(3)   ! and which ones
    integer     :: nstencil
    integer     :: nlevels    ! number of multigrid levels down to a single grid cell
    logical     :: is_singular
  end type hypre_solver
  type(hypre_solver), save :: asolver_p,asolver_d
  logical           , save :: is_init_p = .false.,is_init_d = .false.
  contains
  subroutine solver_vc(ng,lo,hi,nh_p,nh_psi,nh_g,cbc,bc,dli,dzci,dzfi,is_bound,alpha12,psi,p,po)
    !
    ! solves the variable-coefficient Poisson equation for the pressure,
    !   div(1/rho grad p) = rhs, rho = alpha12(2) + (alpha12(1)-alpha12(2))*psi,
    ! with an iterative solver of the HYPRE library, chosen and tuned with the `hypre_*` inputs;
    ! `p` holds the right-hand side on input and the solution on output, and `po` the initial guess
    !
    ! the grid, the matrix and the vectors are created once; the coefficients change with `psi`, so
    ! that the matrix is filled and the solver is set up at every call
    !
    implicit none
    integer             , intent(in   ), dimension(3)                          :: ng,lo,hi
    integer             , intent(in   )                                        :: nh_p,nh_psi,nh_g
    character(len=1)    , intent(in   ), dimension(0:1,3)                      :: cbc
    real(rp)            , intent(in   ), dimension(0:1,3)                      :: bc
    real(rp)            , intent(in   ), dimension(3)                          :: dli
    real(rp)            , intent(in   ), dimension(lo(3)-nh_g:)                :: dzci,dzfi
    logical             , intent(in   ), dimension(0:1,3)                      :: is_bound
    real(rp)            , intent(in   ), dimension(2)                          :: alpha12
    real(rp)            , intent(in   ), dimension(lo(1)-nh_psi:,lo(2)-nh_psi:,lo(3)-nh_psi:) :: psi
    real(rp)            , intent(inout), dimension(lo(1)-nh_p:,lo(2)-nh_p:,lo(3)-nh_p:) :: p,po
    integer, parameter :: iq = HYPREProblemPressure
    !
    if(.not.is_init_p) then
      call init_matrix(ng,lo,hi,cbc,asolver_p)
      is_init_p = .true.
    end if
    call fill_matrix(ng,lo,hi,nh_p,nh_psi,nh_g,cbc,bc,dli,dzci,dzfi,is_bound,alpha12,psi,p,asolver_p)
    call create_solver(iq,asolver_p)
    call setup_solver(asolver_p)
    call solve_helmholtz(iq,asolver_p,ng,lo,hi,nh_p,nh_g,dzfi,p,po)
    call finalize_solver(asolver_p)
  end subroutine solver_vc
  !
  subroutine solver_vc_helmholtz(ng,lo,hi,cbc,coeffs,p)
    !
    ! solves a system whose matrix is given, e.g., that of the implicit diffusion (see `mod_impdiff`),
    ! with the solver of the second entry of the `hypre_*` inputs, and a zero initial guess;
    ! `coeffs(0,:,:,:)` is the diagonal, and `coeffs(2*idir-1,:,:,:)` and `coeffs(2*idir,:,:,:)` are the
    ! coefficients of the lower and upper neighbours along direction `idir`, with the boundary conditions
    ! included; `p` holds the right-hand side on input and the solution on output
    !
    implicit none
    integer         , intent(in   ), dimension(3)     :: ng,lo,hi
    character(len=1), intent(in   ), dimension(0:1,3) :: cbc
    real(rp)        , intent(in   ), dimension(0:,lo(1):,lo(2):,lo(3):) :: coeffs
    real(rp)        , intent(inout), dimension(lo(1):,lo(2):,lo(3):)    :: p
    real(rp), allocatable, dimension(:) :: matvalues
    real(rp), allocatable, dimension(:,:,:) :: po
    real(rp), dimension(lo(3):hi(3)) :: dzfi
    integer, parameter :: iq = HYPREProblemDiffusion
    logical, dimension(3) :: is_dir
    integer :: i,j,k,q,qq,idir,ndim,nstencil
    !
    if(.not.is_init_d) then
      call init_matrix(ng,lo,hi,cbc,asolver_d)
      asolver_d%is_singular = .false.
      is_init_d = .true.
    end if
    ndim = asolver_d%ndim; nstencil = asolver_d%nstencil
    is_dir(:) = .false.
    is_dir(asolver_d%idirs(1:ndim)) = .true.
    allocate(matvalues(product(hi(:)-lo(:)+1)*nstencil))
    q = 0
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          q = q + 1
          qq = (q-1)*nstencil
          matvalues(qq+1) = coeffs(0,i,j,k)
          do idir=1,3
            !
            ! single periodic cell: its neighbours are the cell itself
            !
            if(.not.is_dir(idir)) matvalues(qq+1) = matvalues(qq+1) + coeffs(2*idir-1,i,j,k) + coeffs(2*idir,i,j,k)
          end do
          do idir=1,ndim
            matvalues(qq+2*idir  ) = coeffs(2*asolver_d%idirs(idir)-1,i,j,k)
            matvalues(qq+2*idir+1) = coeffs(2*asolver_d%idirs(idir)  ,i,j,k)
          end do
        end do
      end do
    end do
    call HYPRE_StructMatrixSetBoxValues(asolver_d%mat,lo(asolver_d%idirs(1:ndim)),hi(asolver_d%idirs(1:ndim)),nstencil, &
                                        [(q,q=0,nstencil-1)],matvalues,ierr)
    call HYPRE_StructMatrixAssemble(asolver_d%mat,ierr)
    deallocate(matvalues)
    allocate(po(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    po(:,:,:) = 0._rp
    dzfi(:)   = 1._rp
    call create_solver(iq,asolver_d)
    call setup_solver(asolver_d)
    call solve_helmholtz(iq,asolver_d,ng,lo,hi,0,0,dzfi,p,po)
    call finalize_solver(asolver_d)
  end subroutine solver_vc_helmholtz
  !
  subroutine solver_vc_end
    implicit none
    if(is_init_p) then
      call finalize_matrix(asolver_p)
      is_init_p = .false.
    end if
    if(is_init_d) then
      call finalize_matrix(asolver_d)
      is_init_d = .false.
    end if
  end subroutine solver_vc_end
  !
  subroutine init_matrix(ng,lo,hi,cbc,asolver)
    !
    ! creates the grid, the stencil, the matrix and the vectors of the problem;
    ! a periodic direction with a single grid cell is left out, i.e., a 2D run is a 2D problem,
    ! since the multigrid solvers cannot coarsen such a direction
    !
    implicit none
    integer             , intent(in ), dimension(3)     :: ng,lo,hi
    character(len=1)    , intent(in ), dimension(0:1,3) :: cbc
    type(hypre_solver)  , intent(out)                   :: asolver
    integer, dimension(3,7) :: offsets
    integer, dimension(3) :: periods
    type(C_PTR) :: grid,stencil,mat,rhs,sol
    integer :: q,idir,ndim,nstencil,idirs(3),m
    integer :: comm_hypre
    !
    comm_hypre = MPI_COMM_WORLD
    !
    ndim = 0; idirs(:) = 0
    do idir=1,3
      if(.not.(ng(idir) == 1.and.cbc(0,idir)//cbc(1,idir) == 'PP')) then
        ndim = ndim + 1
        idirs(ndim) = idir
      end if
    end do
    if(ndim == 0) then
      ndim = 1; idirs(1) = 3
    end if
    nstencil = 1+2*ndim
    periods(:) = 0
    do q=1,ndim
      idir = idirs(q)
      if(cbc(0,idir)//cbc(1,idir) == 'PP') periods(q) = ng(idir)
    end do
    !
    ! create grid object
    !
    call HYPRE_StructGridCreate(comm_hypre,ndim,grid,ierr)
    call HYPRE_StructGridSetPeriodic(grid,periods(1:ndim),ierr)
    call HYPRE_StructGridSetExtents(grid,lo(idirs(1:ndim)),hi(idirs(1:ndim)),ierr)
    call HYPRE_StructGridAssemble(grid,ierr)
    !
    ! setup the finite-difference stencil: the cell, and its lower and upper neighbours along each direction
    !
    call HYPRE_StructStencilCreate(ndim,nstencil,stencil,ierr)
    offsets(:,:) = 0
    do q=1,ndim
      offsets(q,2*q  ) = -1
      offsets(q,2*q+1) =  1
    end do
    do q=1,nstencil
      call HYPRE_StructStencilSetElement(stencil,q-1,offsets(1:ndim,q),ierr)
    end do
    !
    ! create coefficient matrix, and solution & right-hand-side vectors
    !
    call HYPRE_StructMatrixCreate(comm_hypre,grid,stencil,mat,ierr)
    call HYPRE_StructMatrixInitialize(mat,ierr)
    call HYPRE_StructVectorCreate(comm_hypre,grid,sol,ierr)
    call HYPRE_StructVectorInitialize(sol,ierr)
    call HYPRE_StructVectorCreate(comm_hypre,grid,rhs,ierr)
    call HYPRE_StructVectorInitialize(rhs,ierr)
    asolver%grid       = grid
    asolver%stencil    = stencil
    asolver%mat        = mat
    asolver%rhs        = rhs
    asolver%sol        = sol
    asolver%precond    = C_NULL_PTR
    asolver%solver     = C_NULL_PTR
    asolver%comm_hypre = comm_hypre
    asolver%ndim       = ndim
    asolver%idirs(:)   = idirs(:)
    asolver%nstencil   = nstencil
    !
    ! without a Dirichlet boundary the solution is defined up to a constant
    !
    asolver%is_singular = all(cbc(:,:) /= 'D')
    !
    ! number of levels of a multigrid hierarchy that halves the grid down to a single cell
    !
    asolver%nlevels = 1
    do q=1,ndim
      m = ng(idirs(q))
      do while(m > 1)
        m = (m+1)/2
        asolver%nlevels = asolver%nlevels + 1
      end do
    end do
  end subroutine init_matrix
  !
  subroutine fill_matrix(ng,lo,hi,nh_p,nh_psi,nh_g,cbc,bc,dli,dzci,dzfi,is_bound,alpha12,psi,p,asolver)
    !
    ! computes the coefficients of the matrix, and adds the contributions of the boundary conditions
    ! to the right-hand side `p`;
    !
    ! each equation is multiplied by the height of its cell, 1/dzfi(k), which makes the matrix symmetric
    ! also on a stretched grid; the right-hand side leaves this routine scaled in the same way
    !
    implicit none
    integer             , intent(in   ), dimension(3)                          :: ng,lo,hi
    integer             , intent(in   )                                        :: nh_p,nh_psi,nh_g
    character(len=1)    , intent(in   ), dimension(0:1,3)                      :: cbc
    real(rp)            , intent(in   ), dimension(0:1,3)                      :: bc
    real(rp)            , intent(in   ), dimension(3)                          :: dli
    real(rp)            , intent(in   ), dimension(lo(3)-nh_g:)                :: dzci,dzfi
    logical             , intent(in   ), dimension(0:1,3)                      :: is_bound
    real(rp)            , intent(in   ), dimension(2)                          :: alpha12
    real(rp)            , intent(in   ), dimension(lo(1)-nh_psi:,lo(2)-nh_psi:,lo(3)-nh_psi:) :: psi
    real(rp)            , intent(inout), dimension(lo(1)-nh_p:,lo(2)-nh_p:,lo(3)-nh_p:) :: p
    type(hypre_solver)  , intent(in   )                                        :: asolver
    real(rp), allocatable, dimension(:) :: matvalues
    real(rp), dimension(0:1,3) :: sgn,factor
    real(rp), dimension(0:1,3) :: cs ! coefficients of the lower and upper neighbours along each direction
    logical , dimension(3) :: is_dir,is_per
    real(rp) :: cc,alpha,dalpha,vol,rhsb
    integer :: i,j,k,q,qq,idir,ib,ndim,nstencil
    integer :: ijk(3),ijkn(3)
    !
    ndim = asolver%ndim; nstencil = asolver%nstencil
    is_dir(:) = .false.
    is_dir(asolver%idirs(1:ndim)) = .true.
    is_per(:) = cbc(0,:)//cbc(1,:) == 'PP'
    sgn(:,:)    = 0._rp
    factor(:,:) = 0._rp
    do q=1,3
      do qq=0,1
        if(is_bound(qq,q)) then
          select case(cbc(qq,q))
          case('N')
            sgn(qq,q)    =  1._rp
            if(q == 3) then
              if(qq == 0) factor(qq,q) =  dzci(lo(3)-1)**(-1)*bc(qq,q)
              if(qq == 1) factor(qq,q) = -dzci(hi(3)  )**(-1)*bc(qq,q)
            else
              if(qq == 0) factor(qq,q) =  dli(q)**(-1)*bc(qq,q)
              if(qq == 1) factor(qq,q) = -dli(q)**(-1)*bc(qq,q)
            end if
          case('D')
            sgn(qq,q)    = -1._rp
            factor(qq,q) = -2._rp*bc(qq,q)
          end select
        end if
      end do
    end do
    !
    alpha = alpha12(2); dalpha = alpha12(1)-alpha12(2)
    allocate(matvalues(product(hi(:)-lo(:)+1)*nstencil))
    q = 0
    do k=lo(3),hi(3)
      vol = dzfi(k)**(-1)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          q = q + 1
          ijk(:) = [i,j,k]
          !
          ! 1/rho at the cell faces, times the metric terms
          !
          do idir=1,3
            do ib=0,1
              ijkn(:) = ijk(:); ijkn(idir) = ijkn(idir) + 2*ib-1
              cs(ib,idir) = (alpha + dalpha*0.5*(psi(ijkn(1),ijkn(2),ijkn(3))+psi(i,j,k)))**(-1)
            end do
          end do
          cs(:,1) = cs(:,1)*dli(1)*dli(1)*vol
          cs(:,2) = cs(:,2)*dli(2)*dli(2)*vol
          cs(0,3) = cs(0,3)*dzci(k-1)
          cs(1,3) = cs(1,3)*dzci(k  )
          !
          cc   = 0._rp
          rhsb = 0._rp
          do idir=1,3
            if(.not.is_dir(idir)) then
              cs(:,idir) = 0._rp ! single periodic cell: its neighbours are the cell itself
              cycle
            end if
            cc = cc - (cs(0,idir)+cs(1,idir))
            if(.not.is_per(idir)) then
              if(is_bound(0,idir).and.ijk(idir) == lo(idir)) then
                cc   = cc   + sgn(0,idir)*cs(0,idir)
                rhsb = rhsb + cs(0,idir)*factor(0,idir)
                cs(0,idir) = 0._rp
              end if
              if(is_bound(1,idir).and.ijk(idir) == hi(idir)) then
                cc   = cc   + sgn(1,idir)*cs(1,idir)
                rhsb = rhsb + cs(1,idir)*factor(1,idir)
                cs(1,idir) = 0._rp
              end if
            end if
          end do
          p(i,j,k) = p(i,j,k)*vol + rhsb
          qq = (q-1)*nstencil
          matvalues(qq+1) = cc
          do idir=1,ndim
            matvalues(qq+2*idir  ) = cs(0,asolver%idirs(idir))
            matvalues(qq+2*idir+1) = cs(1,asolver%idirs(idir))
          end do
        end do
      end do
    end do
    call HYPRE_StructMatrixSetBoxValues(asolver%mat,lo(asolver%idirs(1:ndim)),hi(asolver%idirs(1:ndim)),nstencil, &
                                        [(q,q=0,nstencil-1)],matvalues,ierr)
    call HYPRE_StructMatrixAssemble(asolver%mat,ierr)
    deallocate(matvalues)
  end subroutine fill_matrix
  !
  subroutine configure_pfmg(iq,asolver,solver,maxiter,maxerror)
    implicit none
    integer           , intent(in) :: iq
    type(hypre_solver), intent(in) :: asolver
    type(C_PTR)       , intent(in) :: solver
    integer           , intent(in) :: maxiter
    real(rp)          , intent(in) :: maxerror
    integer :: relax_type,max_levels
    call HYPRE_StructPFMGSetMaxIter(solver,maxiter,ierr)
    call HYPRE_StructPFMGSetTol(solver,maxerror,ierr)
    call HYPRE_StructPFMGSetRelChange(solver,1,ierr)
    !
    ! relaxation method:
    ! 0: Jacobi
    ! 1: weighted Jacobi
    ! 2: red/black Gauss-Seidel (symmetric: RB pre- and post-relaxation), the default here since the matrix is symmetric
    ! 3: red/black Gauss-Seidel (nonsymmetric: RB pre- and post-relaxation)
    !
    relax_type = hypre_pfmg_relax_type(iq)
    if(relax_type < 0) relax_type = 2
    call HYPRE_StructPFMGSetRelaxType(solver,relax_type,ierr)
    call HYPRE_StructPFMGSetNumPreRelax(solver,hypre_pfmg_num_pre_relax(iq),ierr)
    call HYPRE_StructPFMGSetNumPostRelax(solver,hypre_pfmg_num_post_relax(iq),ierr)
    !
    ! number of levels: the operator of a singular problem vanishes on a level with a single grid cell,
    ! which spoils the convergence; by default, such a problem stops one level short of it
    !
    max_levels = hypre_pfmg_max_levels(iq)
    if(max_levels < 0) max_levels = merge(asolver%nlevels-1,0,asolver%is_singular)
    if(max_levels > 0) call HYPRE_StructPFMGSetMaxLevels(solver,max_levels,ierr)
  end subroutine configure_pfmg
  !
  subroutine create_solver(iq,asolver)
    !
    ! creates the solver and its preconditioner from the entry `iq` of the `hypre_*` inputs
    !
    implicit none
    integer           ,         intent(in   ) :: iq
    type(hypre_solver), target, intent(inout) :: asolver
    type(C_PTR) :: solver,precond
    integer  :: stype,precond_id,maxiter
    real(rp) :: maxerror
    !
    ! setup solver
    !
    ! note: this part was taken from the Paris Simulator code
    !       freely available under a GPL license
    !       http://www.ida.upmc.fr/~zaleski/paris
    !
    stype    = hypre_solver_i(iq)
    maxiter  = hypre_maxiter(iq)
    maxerror = hypre_tol(iq)
    precond  = C_NULL_PTR
    precond_id = HYPREPrecondNone
    if      ( stype == HYPRESolverSMG ) then
      call HYPRE_StructSMGCreate(asolver%comm_hypre,solver,ierr)
      call HYPRE_StructSMGSetMaxIter(solver,maxiter,ierr)
      call HYPRE_StructSMGSetTol(solver,maxerror,ierr)
    else if ( stype == HYPRESolverPFMG ) then
      call HYPRE_StructPFMGCreate(asolver%comm_hypre,solver,ierr)
      call configure_pfmg(iq,asolver,solver,maxiter,maxerror)
    else if ( stype == HYPRESolverGMRES .or. &
              stype == HYPRESolverBiCGSTAB .or. &
              stype == HYPRESolverFlexGMRES ) then
      if      ( stype == HYPRESolverGMRES ) then
        call HYPRE_StructGMRESCreate(asolver%comm_hypre,solver,ierr)
        call HYPRE_StructGMRESSetMaxIter(solver,maxiter,ierr)
        call HYPRE_StructGMRESSetTol(solver,maxerror,ierr)
        if(hypre_gmres_k_dim(iq) > 0) call HYPRE_StructGMRESSetKDim(solver,hypre_gmres_k_dim(iq),ierr)
      else if ( stype == HYPRESolverBiCGSTAB ) then
        call HYPRE_StructBiCGSTABCreate(asolver%comm_hypre,solver,ierr)
        call HYPRE_StructBiCGSTABSetMaxIter(solver,maxiter,ierr)
        call HYPRE_StructBiCGSTABSetTol(solver,maxerror,ierr)
      else if ( stype == HYPRESolverFlexGMRES ) then
        call HYPRE_StructFGMRESCreate(asolver%comm_hypre,solver,ierr)
        call HYPRE_StructFGMRESSetMaxIter(solver,maxiter,ierr)
        call HYPRE_StructFGMRESSetTol(solver,maxerror,ierr)
        if(hypre_gmres_k_dim(iq) > 0) call HYPRE_StructFGMRESSetKDim(solver,hypre_gmres_k_dim(iq),ierr)
      end if
      !
      ! preconditioner
      !
      precond_id = hypre_precond_i(iq)
      if(precond_id == HYPREPrecondSMG) then
        call HYPRE_StructSMGCreate(asolver%comm_hypre,precond,ierr)
        call HYPRE_StructSMGSetMaxIter(precond,hypre_precond_maxiter(iq),ierr)
        call HYPRE_StructSMGSetTol(precond,0._rp,ierr)
        call HYPRE_StructSMGSetZeroGuess(precond,ierr)
      else if(precond_id == HYPREPrecondPFMG) then
        call HYPRE_StructPFMGCreate(asolver%comm_hypre,precond,ierr)
        call configure_pfmg(iq,asolver,precond,hypre_precond_maxiter(iq),0._rp)
        call HYPRE_StructPFMGSetZeroGuess(precond,ierr)
      else
        precond_id = HYPREPrecondNone
      end if
      if      ( stype == HYPRESolverGMRES ) then
        call HYPRE_StructGMRESSetPrecond(solver,precond_id,precond,ierr)
      else if ( stype == HYPRESolverBiCGSTAB ) then
        call HYPRE_StructBiCGSTABSetPrecond(solver,precond_id,precond,ierr)
      else if ( stype == HYPRESolverFlexGMRES ) then
        call HYPRE_StructFGMRESSetPrecond(solver,precond_id,precond,ierr)
      end if
    end if
    asolver%precond = precond
    asolver%solver  = solver
    asolver%stype   = stype
    asolver%ptype   = precond_id
  end subroutine create_solver
  !
  subroutine setup_solver(asolver)
    implicit none
    type(hypre_solver), target, intent(inout) :: asolver
    type(C_PTR), pointer :: solver,mat,rhs,sol
    integer    , pointer :: stype
    !
    solver => asolver%solver
    stype  => asolver%stype
    mat    => asolver%mat
    rhs    => asolver%rhs
    sol    => asolver%sol
    !
    ! setup solver
    !
    ! note: this part was taken from the Paris Simulator code
    !       freely available under a GPL license
    !       http://www.ida.upmc.fr/~zaleski/paris
    !
    if      ( stype == HYPRESolverSMG ) then
      call HYPRE_StructSMGSetup(solver,mat,rhs,sol,ierr)
    else if ( stype == HYPRESolverPFMG ) then
      call HYPRE_StructPFMGSetup(solver,mat,rhs,sol,ierr)
    else if ( stype == HYPRESolverGMRES ) then
      call HYPRE_StructGMRESSetup(solver,mat,rhs,sol,ierr)
    else if ( stype == HYPRESolverBiCGSTAB ) then
      call HYPRE_StructBiCGSTABSetup(solver,mat,rhs,sol,ierr)
    else if ( stype == HYPRESolverFlexGMRES ) then
      call HYPRE_StructFGMRESSetup(solver,mat,rhs,sol,ierr)
    end if
  end subroutine setup_solver
  !
  subroutine solve_helmholtz(iq,asolver,ng,lo,hi,nh_p,nh_g,dzfi,p,po)
    !
    ! solves the system for the (scaled) right-hand side `p` with initial guess `po`, and returns the
    ! solution in both; warns if the solver did not converge within its maximum number of iterations
    !
    implicit none
    integer           ,         intent(in   )               :: iq
    type(hypre_solver), target, intent(in   )               :: asolver
    integer           ,         intent(in   ), dimension(3) :: ng,lo,hi
    integer           ,         intent(in   )               :: nh_p,nh_g
    real(rp)          ,         intent(in   ), dimension(lo(3)-nh_g:) :: dzfi
    real(rp)          ,         intent(inout), dimension(lo(1)-nh_p:,lo(2)-nh_p:,lo(3)-nh_p:) :: p,po
    type(C_PTR), pointer :: solver,mat,rhs,sol
    integer    , pointer :: stype
    integer , dimension(asolver%ndim) :: lo_s,hi_s
    real(rp) :: pmean(2)
    integer  :: k,niter
    solver  => asolver%solver
    mat     => asolver%mat
    rhs     => asolver%rhs
    sol     => asolver%sol
    stype   => asolver%stype
    lo_s(:) = lo(asolver%idirs(1:asolver%ndim))
    hi_s(:) = hi(asolver%idirs(1:asolver%ndim))
    !
    if(asolver%is_singular) then
      !
      ! the right-hand side of a singular problem must sum to zero: remove what round-off leaves
      !
      pmean(1) = sum(p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      call MPI_ALLREDUCE(MPI_IN_PLACE,pmean(1),1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      pmean(1) = pmean(1)/(1._rp*ng(1)*ng(2)*ng(3))
      p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) - pmean(1)
    end if
    call HYPRE_StructVectorSetBoxValues(rhs,lo_s,hi_s,p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),ierr)
    call HYPRE_StructVectorAssemble(rhs,ierr)
    !
    ! create soluction vector
    !
    call HYPRE_StructVectorSetBoxValues(sol,lo_s,hi_s,po(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),ierr)
    call HYPRE_StructVectorAssemble(sol,ierr)
    !
    ! setup solver, and solve
    !
    ! note: this part was based on the the Paris Simulator code
    !       freely available under a GPL license; see:
    !       http://www.ida.upmc.fr/~zaleski/paris/
    !
    niter = 0
    if      ( stype == HYPRESolverSMG ) then
      call HYPRE_StructSMGSolve(solver,mat,rhs,sol,ierr)
      call HYPRE_StructSMGGetNumIterations(solver,niter,ierr)
    else if ( stype == HYPRESolverPFMG ) then
      call HYPRE_StructPFMGSolve(solver,mat,rhs,sol,ierr)
      call HYPRE_StructPFMGGetNumIteration(solver,niter,ierr)
    else if ( stype == HYPRESolverGMRES ) then
      call HYPRE_StructGMRESSolve(solver,mat,rhs,sol,ierr)
      call HYPRE_StructGMRESGetNumIteratio(solver,niter,ierr)
    else if ( stype == HYPRESolverBiCGSTAB ) then
      call HYPRE_StructBiCGSTABSolve(solver,mat,rhs,sol,ierr)
      call HYPRE_StructBiCGSTABGetNumItera(solver,niter,ierr)
    else if ( stype == HYPRESolverFlexGMRES ) then
      call HYPRE_StructFGMRESSolve(solver,mat,rhs,sol,ierr)
      call HYPRE_StructFGMRESGetNumIter(solver,niter,ierr)
    end if ! stype
    !
    ! end of part based on the Paris Simulator code
    !
    if(niter >= hypre_maxiter(iq)) then
      if(myid == 0) print*, 'WARNING: HYPRE solver not converged after ',niter,' iterations (maximum: `hypre_maxiter`).'
    end if
    !
    ! fecth results and update old solution
    !
    call HYPRE_StructVectorGetBoxValues(sol,lo_s,hi_s,p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),ierr)
    if(asolver%is_singular) then
      !
      ! fix the constant of the solution: zero mean
      !
      pmean(:) = 0._rp
      do k=lo(3),hi(3)
        pmean(1) = pmean(1) + sum(p(lo(1):hi(1),lo(2):hi(2),k))/dzfi(k)
        pmean(2) = pmean(2) + (hi(1)-lo(1)+1)*(hi(2)-lo(2)+1)/dzfi(k)
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,pmean,2,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      pmean(1) = pmean(1)/pmean(2)
      p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) - pmean(1)
    end if
    po(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  end subroutine solve_helmholtz
  !
  subroutine finalize_solver(asolver)
    implicit none
    type(hypre_solver), target, intent(in) :: asolver
    type(C_PTR), pointer :: precond,solver
    integer    , pointer :: stype,ptype
    !
    precond => asolver%precond
    solver  => asolver%solver
    stype   => asolver%stype
    ptype   => asolver%ptype
    !
    ! note: this part was based on the the Paris Simulator code
    !       freely available under a GPL license; see:
    !       http://www.ida.upmc.fr/~zaleski/paris/
    !
    if      ( stype == HYPRESolverSMG ) then
      call HYPRE_StructSMGDestroy(solver,ierr)
    else if ( stype == HYPRESolverPFMG ) then
      call HYPRE_StructPFMGDestroy(solver,ierr)
    else if ( stype == HYPRESolverGMRES ) then
      call HYPRE_StructGMRESDestroy(solver,ierr)
    else if ( stype == HYPRESolverBiCGSTAB ) then
      call HYPRE_StructBiCGSTABDestroy(solver,ierr)
    else if ( stype == HYPRESolverFlexGMRES ) then
      call HYPRE_StructFGMRESDestroy(solver,ierr)
    end if
    if      ( ptype == HYPREPrecondSMG ) then
      call HYPRE_StructSMGDestroy(precond,ierr)
    else if ( ptype == HYPREPrecondPFMG ) then
      call HYPRE_StructPFMGDestroy(precond,ierr)
    end if
  end subroutine finalize_solver
  !
  subroutine finalize_matrix(asolver)
    implicit none
    type(hypre_solver), target, intent(in) :: asolver
    type(C_PTR), pointer :: grid,stencil,mat,rhs,sol
    !
    grid    => asolver%grid
    stencil => asolver%stencil
    mat     => asolver%mat
    rhs     => asolver%rhs
    sol     => asolver%sol
    !
    call HYPRE_StructGridDestroy(grid,ierr)
    call HYPRE_StructStencilDestroy(stencil,ierr)
    call HYPRE_StructMatrixDestroy(mat,ierr)
    call HYPRE_StructVectorDestroy(rhs,ierr)
    call HYPRE_StructVectorDestroy(sol,ierr)
  end subroutine finalize_matrix
#endif
end module mod_solver_vc
