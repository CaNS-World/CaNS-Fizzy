! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
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
                            hypre_pfmg_max_levels,hypre_is_symmetric,hypre_print_level, &
                            is_impdiff,is_impdiff_1d
  use mod_types
  use, intrinsic :: iso_c_binding, only: C_PTR,C_NULL_PTR,C_INT,C_SIZE_T,C_ASSOCIATED
#if defined(_OPENACC) || defined(_OPENMP)
  use cudafor       , only: cudaDeviceSynchronize
#endif
  implicit none
  integer :: ierr
  private
  public solver_vc,solver_vc_helmholtz,solver_vc_init,solver_vc_end
  integer, parameter :: HYPRESolverSMG       = 1, &
                        HYPRESolverPFMG      = 2, &
                        HYPRESolverGMRES     = 3, &
                        HYPRESolverBiCGSTAB  = 4, &
                        HYPRESolverFlexGMRES = 5, &
                        HYPRESolverPCG       = 6
  integer, parameter :: HYPREPrecondSMG  = 0, &
                        HYPREPrecondPFMG = 1, &
                        HYPREPrecondNone = 9
  integer, parameter :: HYPREProblemPressure  = 1, & ! entries of the `hypre_*` inputs
                        HYPREProblemDiffusion = 2
  type hypre_solver
    type(C_PTR) :: grid,stencil,precond,solver,mat,rhs,sol
    integer     :: stype,ptype,comm_hypre
    integer     :: ndim         ! number of directions of the problem (those with more than one grid cell)
    integer     :: idirs(3)     ! and which ones
    integer     :: nstencil
    logical     :: is_symmetric ! symmetric storage: HYPRE stores the diagonal and the lower neighbours only
    integer     :: nvalues      ! number of values per grid cell passed to HYPRE in `matvalues`
    integer     :: ientries(7)  ! and their stencil entries
    integer     :: ival(0:1,3)  ! position among them of the lower (0) and upper (1) neighbour along each
                                ! direction (0: not passed, i.e., direction left out, or upper neighbour
                                ! with symmetric storage)
    integer     :: nlevels      ! number of multigrid levels down to a single grid cell
    logical     :: is_singular
  end type hypre_solver
  type(hypre_solver), save :: asolver_p,asolver_d
  logical           , save :: is_init_p = .false.,is_init_d = .false.
  !
  ! values of the matrix and of the vectors passed to and from HYPRE, in the memory of HYPRE
  ! (the device memory in GPU builds); shared by both problems, which have the same grid
  !
  real(rp), allocatable, dimension(:)    , save :: matvalues,rhsvalues,solvalues
  real(rp), allocatable, dimension(:,:,:), save :: po_d ! (zero) initial guess of the implicit diffusion
#if defined(_OPENACC) || defined(_OPENMP)
  integer, parameter :: HYPREMemoryDevice = 1, & ! `HYPRE_MEMORY_DEVICE` and `HYPRE_EXEC_DEVICE` of HYPREf.h
                        HYPREExecDevice   = 1
  !
  ! HYPRE has no Fortran wrappers for these functions
  !
  interface
    integer(C_INT) function HYPRE_DeviceInitialize() bind(C,name='HYPRE_DeviceInitialize')
      import :: C_INT
    end function HYPRE_DeviceInitialize
#if defined(_HYPRE_UMPIRE)
    integer(C_INT) function HYPRE_SetUmpireDevicePoolSize(nbytes) bind(C,name='HYPRE_SetUmpireDevicePoolSize')
      import :: C_INT,C_SIZE_T
      integer(C_SIZE_T), value :: nbytes
    end function HYPRE_SetUmpireDevicePoolSize
#endif
  end interface
#endif
  contains
  subroutine solver_vc(ng,lo,hi,nh_p,nh_psi,nh_g,cbc,bc,dli,dzci,dzf,wzc,is_bound,alpha12,psi,p,po)
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
    real(rp)            , intent(in   ), dimension(lo(3)-nh_g:)                :: dzci,dzf,wzc
    logical             , intent(in   ), dimension(0:1,3)                      :: is_bound
    real(rp)            , intent(in   ), dimension(2)                          :: alpha12
    real(rp)            , intent(in   ), dimension(lo(1)-nh_psi:,lo(2)-nh_psi:,lo(3)-nh_psi:) :: psi
    real(rp)            , intent(inout), dimension(lo(1)-nh_p:,lo(2)-nh_p:,lo(3)-nh_p:) :: p,po
    integer, parameter :: iq = HYPREProblemPressure
    !
    if(.not.is_init_p) then
      call init_matrix(ng,lo,hi,cbc,hypre_is_symmetric(iq),asolver_p)
      is_init_p = .true.
    end if
    call fill_matrix(ng,lo,hi,nh_p,nh_psi,nh_g,cbc,bc,dli,dzci,dzf,wzc,is_bound,alpha12,psi,p,asolver_p)
    call create_solver(iq,asolver_p)
    call setup_solver(asolver_p)
    call solve_helmholtz(iq,asolver_p,ng,lo,hi,nh_p,nh_g,dzf,p,po)
    !
    ! the solver is kept for the next step, unless HYPRE solves other systems (implicit diffusion): its work
    ! vectors would then stay in the memory pool of HYPRE in between, which keeps the pool from reusing its
    ! free blocks, and its footprint grows
    !
    call finalize_solver(asolver_p,is_impdiff .and. .not.is_impdiff_1d)
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
    real(rp)        , intent(in   ), dimension(0:,lo(1):,lo(2):,lo(3):), contiguous :: coeffs
    real(rp)        , intent(inout), dimension(lo(1):,lo(2):,lo(3):)                :: p
    real(rp), dimension(lo(3):hi(3)) :: dzf
    integer, parameter :: iq = HYPREProblemDiffusion
    real(rp) :: cc
    integer  :: mxm,mxp,mym,myp,mzm,mzp ! positions of the coefficients among the values of a grid cell
    integer  :: i,j,k,q,qq,ndim,nstencil,nvalues,nx,ny
    !
    if(.not.is_init_d) then
      call init_matrix(ng,lo,hi,cbc,hypre_is_symmetric(iq),asolver_d)
      asolver_d%is_singular = .false.
      is_init_d = .true.
      allocate(po_d(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      !$acc        enter data create(   po_d)
      !$omp target enter data map(alloc:po_d)
    end if
    ndim = asolver_d%ndim; nstencil = asolver_d%nstencil; nvalues = asolver_d%nvalues
    mxm = asolver_d%ival(0,1); mxp = asolver_d%ival(1,1)
    mym = asolver_d%ival(0,2); myp = asolver_d%ival(1,2)
    mzm = asolver_d%ival(0,3); mzp = asolver_d%ival(1,3)
    nx = hi(1)-lo(1)+1; ny = hi(2)-lo(2)+1
    if(ndim == 3) then
      !
      ! the entries of `coeffs` are those of the stencil of HYPRE, in the same order
      ! (with symmetric storage HYPRE ignores the upper neighbours)
      !
      !$acc parallel     loop collapse(3) default(present) async(1)
      !$omp target teams loop collapse(3)
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            po_d(i,j,k) = 0._rp
          end do
        end do
      end do
      !$acc wait(1)
      !$acc   host_data use_device(     coeffs)
      !$omp target data use_device_addr(coeffs)
      call HYPRE_StructMatrixSetBoxValues(asolver_d%mat,lo,hi,nstencil,[(q,q=0,nstencil-1)],coeffs,ierr)
      !$omp end target data
      !$acc end   host_data
    else
      call alloc_matvalues(product(hi(:)-lo(:)+1)*nvalues)
      !$acc parallel     loop collapse(3) default(present) private(q,qq,cc) async(1)
      !$omp target teams loop collapse(3)                  private(q,qq,cc)
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            q = 1 + (i-lo(1)) + nx*((j-lo(2)) + ny*(k-lo(3)))
            qq = (q-1)*nvalues
            !
            ! a direction left out of the problem (no lower neighbour, e.g. mxm = 0) has a single periodic cell,
            ! whose neighbours are the cell itself
            !
            cc = coeffs(0,i,j,k)
            if(mxm == 0) cc = cc + coeffs(1,i,j,k) + coeffs(2,i,j,k)
            if(mym == 0) cc = cc + coeffs(3,i,j,k) + coeffs(4,i,j,k)
            if(mzm == 0) cc = cc + coeffs(5,i,j,k) + coeffs(6,i,j,k)
            matvalues(qq+1) = cc
            if(mxm > 0) matvalues(qq+mxm) = coeffs(1,i,j,k)
            if(mxp > 0) matvalues(qq+mxp) = coeffs(2,i,j,k)
            if(mym > 0) matvalues(qq+mym) = coeffs(3,i,j,k)
            if(myp > 0) matvalues(qq+myp) = coeffs(4,i,j,k)
            if(mzm > 0) matvalues(qq+mzm) = coeffs(5,i,j,k)
            if(mzp > 0) matvalues(qq+mzp) = coeffs(6,i,j,k)
            po_d(i,j,k) = 0._rp
          end do
        end do
      end do
      !$acc wait(1)
      !$acc   host_data use_device(     matvalues)
      !$omp target data use_device_addr(matvalues)
      call HYPRE_StructMatrixSetBoxValues(asolver_d%mat,lo(asolver_d%idirs(1:ndim)),hi(asolver_d%idirs(1:ndim)),nvalues, &
                                          asolver_d%ientries(1:nvalues),matvalues,ierr)
      !$omp end target data
      !$acc end   host_data
    end if
    call HYPRE_StructMatrixAssemble(asolver_d%mat,ierr)
    dzf(:) = 1._rp
    call create_solver(iq,asolver_d)
    call setup_solver(asolver_d)
    call solve_helmholtz(iq,asolver_d,ng,lo,hi,0,0,dzf,p,po_d)
    call finalize_solver(asolver_d,.true.)
  end subroutine solver_vc_helmholtz
  !
  subroutine solver_vc_init
    !
    ! initializes HYPRE, before any other call to it; a GPU build calls this after `initmpi` has set
    ! the device of the task
    !
    implicit none
    !
    ! `HYPRE_Init` is the Fortran name of `HYPRE_Initialize` in all versions of HYPRE (the latter
    ! has a Fortran wrapper only since version 2.29)
    !
    call HYPRE_Init(ierr)
#if defined(_OPENACC) || defined(_OPENMP)
#if defined(_HYPRE_UMPIRE)
    !
    ! the device memory pool of HYPRE starts at this size (4 GiB by default) and grows on demand;
    ! a small initial size leaves room for several tasks sharing one GPU
    !
    ierr = HYPRE_SetUmpireDevicePoolSize(512_C_SIZE_T*1024**2)
#endif
    !
    ! HYPRE takes the device of the task only here, as in its test drivers; otherwise its memory pool
    ! stays on device 0, and a task on another device fails
    !
    ierr = HYPRE_DeviceInitialize()
    !
    ! objects in device memory, and computations on the device
    !
    call HYPRE_SetMemoryLocation(HYPREMemoryDevice,ierr)
    call HYPRE_SetExecutionPolicy(HYPREExecDevice,ierr)
#endif
  end subroutine solver_vc_init
  !
  subroutine solver_vc_end
    !
    ! destroys the HYPRE objects, and finalizes HYPRE
    !
    implicit none
    if(is_init_p) then
      call finalize_solver(asolver_p,.true.)
      call finalize_matrix(asolver_p)
      is_init_p = .false.
    end if
    if(is_init_d) then
      call finalize_matrix(asolver_d)
      is_init_d = .false.
      !$acc        exit data delete(    po_d)
      !$omp target exit data map(delete:po_d)
      deallocate(po_d)
    end if
    if(allocated(matvalues)) then
      !$acc        exit data delete(    matvalues)
      !$omp target exit data map(delete:matvalues)
      deallocate(matvalues)
    end if
    if(allocated(rhsvalues)) then
      !$acc        exit data delete(    rhsvalues,solvalues)
      !$omp target exit data map(delete:rhsvalues,solvalues)
      deallocate(rhsvalues,solvalues)
    end if
    call HYPRE_Finalize(ierr)
  end subroutine solver_vc_end
  !
  subroutine init_matrix(ng,lo,hi,cbc,is_symmetric,asolver)
    !
    ! creates the grid, the stencil, the matrix and the vectors of the problem;
    ! a periodic direction with a single grid cell is left out, i.e., a 2D run is a 2D problem,
    ! since the multigrid solvers cannot coarsen such a direction
    !
    ! with symmetric storage (`is_symmetric`; both matrices of this module are symmetric, see `fill_matrix`
    ! and `mod_impdiff`), HYPRE stores the diagonal and the coefficients of the lower neighbours only, and
    ! takes the coefficient of an upper neighbour from the lower one of the next cell; values given for
    ! the upper neighbours are ignored
    !
    implicit none
    integer             , intent(in ), dimension(3)     :: ng,lo,hi
    character(len=1)    , intent(in ), dimension(0:1,3) :: cbc
    logical             , intent(in )                   :: is_symmetric
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
      if(.not.(ng(idir) == 1 .and. cbc(0,idir)//cbc(1,idir) == 'PP')) then
        ndim = ndim + 1
        idirs(ndim) = idir
      end if
    end do
    if(ndim == 0) then
      ndim = 1; idirs(1) = 3
    end if
    nstencil = 1+2*ndim
    !
    ! no symmetric storage for a problem with a direction with a single grid cell, whose couplings all
    ! vanish: PFMG has failed on such a problem with symmetric storage (HYPRE issue #1510)
    !
    asolver%is_symmetric = is_symmetric
    if(is_symmetric .and. any(ng(idirs(1:ndim)) == 1)) then
      asolver%is_symmetric = .false.
      if(myid == 0) print*, 'WARNING: symmetric storage of a HYPRE matrix not used: a direction has a single grid cell.'
    end if
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
    if(asolver%is_symmetric) call HYPRE_StructMatrixSetSymmetric(mat,1,ierr)
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
    ! values per grid cell passed to HYPRE: the diagonal, and the lower and upper neighbours along each
    ! direction of the problem (stencil entries 0, 2*q-1 and 2*q), or the lower ones only (symmetric storage)
    !
    asolver%ival(:,:) = 0
    if(.not.asolver%is_symmetric) then
      asolver%nvalues = nstencil
      asolver%ientries(1:nstencil) = [(q,q=0,nstencil-1)]
      do q=1,ndim
        asolver%ival(0,idirs(q)) = 2*q
        asolver%ival(1,idirs(q)) = 2*q+1
      end do
    else
      asolver%nvalues = 1+ndim
      asolver%ientries(1:1+ndim) = [0,(2*q-1,q=1,ndim)]
      do q=1,ndim
        asolver%ival(0,idirs(q)) = 1+q
      end do
    end if
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
    !
    ! work vectors of the values passed to HYPRE (see `alloc_matvalues` for the matrix)
    !
    if(.not.allocated(rhsvalues)) then
      allocate(rhsvalues(product(hi(:)-lo(:)+1)), &
               solvalues(product(hi(:)-lo(:)+1)))
      !$acc        enter data create(   rhsvalues,solvalues)
      !$omp target enter data map(alloc:rhsvalues,solvalues)
    end if
  end subroutine init_matrix
  !
  subroutine alloc_matvalues(n)
    !
    ! allocates `matvalues` with at least `n` entries; the two problems may pass different numbers of
    ! values per grid cell, and a 3D diffusion problem passes its matrix directly
    !
    implicit none
    integer, intent(in) :: n
    if(allocated(matvalues)) then
      if(size(matvalues) >= n) return
      !$acc        exit data delete(    matvalues)
      !$omp target exit data map(delete:matvalues)
      deallocate(matvalues)
    end if
    allocate(matvalues(n))
    !$acc        enter data create(   matvalues)
    !$omp target enter data map(alloc:matvalues)
  end subroutine alloc_matvalues
  !
  subroutine fill_matrix(ng,lo,hi,nh_p,nh_psi,nh_g,cbc,bc,dli,dzci,dzf,wzc,is_bound,alpha12,psi,p,asolver)
    !
    ! computes the coefficients of the matrix, and adds the contributions of the boundary conditions
    ! to the right-hand side `p`;
    !
    ! each equation is multiplied by the height of its cell, dzf(k), which makes the matrix symmetric
    ! also on a stretched grid; the right-hand side leaves this routine scaled in the same way
    !
    implicit none
    integer             , intent(in   ), dimension(3)                          :: ng,lo,hi
    integer             , intent(in   )                                        :: nh_p,nh_psi,nh_g
    character(len=1)    , intent(in   ), dimension(0:1,3)                      :: cbc
    real(rp)            , intent(in   ), dimension(0:1,3)                      :: bc
    real(rp)            , intent(in   ), dimension(3)                          :: dli
    real(rp)            , intent(in   ), dimension(lo(3)-nh_g:)                :: dzci,dzf,wzc
    logical             , intent(in   ), dimension(0:1,3)                      :: is_bound
    real(rp)            , intent(in   ), dimension(2)                          :: alpha12
    real(rp)            , intent(in   ), dimension(lo(1)-nh_psi:,lo(2)-nh_psi:,lo(3)-nh_psi:) :: psi
    real(rp)            , intent(inout), dimension(lo(1)-nh_p:,lo(2)-nh_p:,lo(3)-nh_p:) :: p
    type(hypre_solver)  , intent(in   )                                        :: asolver
    real(rp), dimension(0:1,3) :: sgn,factor
    logical , dimension(3) :: is_per
    real(rp) :: cc,alpha,dalpha,vol,rhsb,sgn_op,dxi,dyi
    real(rp) :: cxm,cxp,cym,cyp,czm,czp ! coefficients of the lower and upper neighbours along each direction
    real(rp) :: sxm,sxp,sym,syp,szm,szp,fxm,fxp,fym,fyp,fzm,fzp
    logical  :: is_xm,is_xp,is_ym,is_yp,is_zm,is_zp
    integer  :: mxm,mxp,mym,myp,mzm,mzp ! positions of the coefficients among the values of a grid cell
    integer  :: i,j,k,q,qq,ndim,nvalues,nx,ny
    !
    ndim = asolver%ndim; nvalues = asolver%nvalues
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
    ! PCG needs a positive-definite matrix: with it, the system is assembled and solved with the opposite sign
    !
    sgn_op = merge(-1._rp,1._rp,hypre_solver_i(HYPREProblemPressure) == HYPRESolverPCG)
    !
    alpha = alpha12(2); dalpha = alpha12(1)-alpha12(2)
    nx = hi(1)-lo(1)+1; ny = hi(2)-lo(2)+1
    dxi = dli(1); dyi = dli(2)
    mxm = asolver%ival(0,1); mxp = asolver%ival(1,1)
    mym = asolver%ival(0,2); myp = asolver%ival(1,2)
    mzm = asolver%ival(0,3); mzp = asolver%ival(1,3)
    call alloc_matvalues(product(hi(:)-lo(:)+1)*nvalues)
    !
    ! boundary condition of each face of the domain on this task, as scalars
    !
    is_xm = is_bound(0,1) .and. .not.is_per(1); is_xp = is_bound(1,1) .and. .not.is_per(1)
    is_ym = is_bound(0,2) .and. .not.is_per(2); is_yp = is_bound(1,2) .and. .not.is_per(2)
    is_zm = is_bound(0,3) .and. .not.is_per(3); is_zp = is_bound(1,3) .and. .not.is_per(3)
    sxm = sgn(0,1); sxp = sgn(1,1); fxm = factor(0,1); fxp = factor(1,1)
    sym = sgn(0,2); syp = sgn(1,2); fym = factor(0,2); fyp = factor(1,2)
    szm = sgn(0,3); szp = sgn(1,3); fzm = factor(0,3); fzp = factor(1,3)
    !$acc parallel     loop collapse(3) default(present) private(cc,vol,rhsb,q,qq,cxm,cxp,cym,cyp,czm,czp) async(1)
    !$omp target teams loop collapse(3)                  private(cc,vol,rhsb,q,qq,cxm,cxp,cym,cyp,czm,czp)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          q = 1 + (i-lo(1)) + nx*((j-lo(2)) + ny*(k-lo(3)))
          vol = dzf(k)
          !
          ! 1/rho at the cell faces (volume average along z), times the metric terms
          !
          cxm = (alpha + dalpha*0.5*(psi(i-1,j,k)+psi(i,j,k)))**(-1)
          cxp = (alpha + dalpha*0.5*(psi(i+1,j,k)+psi(i,j,k)))**(-1)
          cym = (alpha + dalpha*0.5*(psi(i,j-1,k)+psi(i,j,k)))**(-1)
          cyp = (alpha + dalpha*0.5*(psi(i,j+1,k)+psi(i,j,k)))**(-1)
          czm = (alpha + dalpha*(wzc(k-1)*psi(i,j,k-1)+(1.-wzc(k-1))*psi(i,j,k  )))**(-1)
          czp = (alpha + dalpha*(wzc(k  )*psi(i,j,k  )+(1.-wzc(k  ))*psi(i,j,k+1)))**(-1)
          cxm = cxm*dxi*dxi*vol; cxp = cxp*dxi*dxi*vol
          cym = cym*dyi*dyi*vol; cyp = cyp*dyi*dyi*vol
          czm = czm*dzci(k-1)
          czp = czp*dzci(k  )
          !
          ! a direction left out of the problem (no lower neighbour, e.g. mxm = 0) has a single periodic cell,
          ! whose neighbours are the cell itself
          !
          cc   = 0._rp
          rhsb = 0._rp
          if(mxm > 0) then
            cc = cc - (cxm+cxp)
            call bc_diag_rhs(is_xm .and. i == lo(1),sxm,fxm,cxm,cc,rhsb)
            call bc_diag_rhs(is_xp .and. i == hi(1),sxp,fxp,cxp,cc,rhsb)
          end if
          if(mym > 0) then
            cc = cc - (cym+cyp)
            call bc_diag_rhs(is_ym .and. j == lo(2),sym,fym,cym,cc,rhsb)
            call bc_diag_rhs(is_yp .and. j == hi(2),syp,fyp,cyp,cc,rhsb)
          end if
          if(mzm > 0) then
            cc = cc - (czm+czp)
            call bc_diag_rhs(is_zm .and. k == lo(3),szm,fzm,czm,cc,rhsb)
            call bc_diag_rhs(is_zp .and. k == hi(3),szp,fzp,czp,cc,rhsb)
          end if
          p(i,j,k) = sgn_op*(p(i,j,k)*vol + rhsb)
          qq = (q-1)*nvalues
          matvalues(qq+1) = sgn_op*cc
          if(mxm > 0) matvalues(qq+mxm) = sgn_op*cxm
          if(mxp > 0) matvalues(qq+mxp) = sgn_op*cxp
          if(mym > 0) matvalues(qq+mym) = sgn_op*cym
          if(myp > 0) matvalues(qq+myp) = sgn_op*cyp
          if(mzm > 0) matvalues(qq+mzm) = sgn_op*czm
          if(mzp > 0) matvalues(qq+mzp) = sgn_op*czp
        end do
      end do
    end do
    !$acc wait(1)
    !$acc   host_data use_device(     matvalues)
    !$omp target data use_device_addr(matvalues)
    call HYPRE_StructMatrixSetBoxValues(asolver%mat,lo(asolver%idirs(1:ndim)),hi(asolver%idirs(1:ndim)),nvalues, &
                                        asolver%ientries(1:nvalues),matvalues,ierr)
    !$omp end target data
    !$acc end   host_data
    call HYPRE_StructMatrixAssemble(asolver%mat,ierr)
  end subroutine fill_matrix
  !
  pure subroutine bc_diag_rhs(is_at,sgn,factor,cf,cc,rhsb)
    !
    ! boundary condition on a face of the domain, for a cell next to it (`is_at`): the coefficient `cf` of
    ! the neighbour across the face goes to the diagonal `cc` with the sign `sgn`, and to the right-hand
    ! side `rhsb` with the factor of the boundary value `factor` (see `fill_matrix`)
    !
    !$acc routine seq
    !$omp declare target
    implicit none
    logical , intent(in   ) :: is_at
    real(rp), intent(in   ) :: sgn,factor
    real(rp), intent(inout) :: cf,cc,rhsb
    if(is_at) then
      cc   = cc   + sgn*cf
      rhsb = rhsb + cf*factor
      cf   = 0._rp
    end if
  end subroutine bc_diag_rhs
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
    ! 2: red/black Gauss-Seidel (symmetric: RB pre- and post-relaxation), the default here since the matrix is symmetric,
    !    except on GPUs, where HYPRE's implementation is several times slower than weighted Jacobi, and with PCG,
    !    which it can stall (e.g., with a periodic direction of a few grid cells)
    ! 3: red/black Gauss-Seidel (nonsymmetric: RB pre- and post-relaxation)
    !
    relax_type = hypre_pfmg_relax_type(iq)
    if(relax_type < 0) then
#if !(defined(_OPENACC) || defined(_OPENMP))
      relax_type = 2
      if(hypre_solver_i(iq) == HYPRESolverPCG) relax_type = 1
#else
      relax_type = 1
#endif
    end if
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
    ! creates the solver and its preconditioner from the entry `iq` of the `hypre_*` inputs;
    ! a Krylov solver kept from a previous solve (see `finalize_solver`) is used again, with its work
    ! vectors, while its multigrid preconditioner is created anew, since HYPRE cannot set up a multigrid
    ! solver twice
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
    stype    = hypre_solver_i(iq)
    maxiter  = hypre_maxiter(iq)
    maxerror = hypre_tol(iq)
    solver   = asolver%solver
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
              stype == HYPRESolverFlexGMRES .or. &
              stype == HYPRESolverPCG ) then
      if      ( C_ASSOCIATED(solver) ) then
        continue ! (kept from a previous solve)
      else if ( stype == HYPRESolverGMRES ) then
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
      else if ( stype == HYPRESolverPCG ) then
        call HYPRE_StructPCGCreate(asolver%comm_hypre,solver,ierr)
        call HYPRE_StructPCGSetMaxIter(solver,maxiter,ierr)
        call HYPRE_StructPCGSetTol(solver,maxerror,ierr)
        call HYPRE_StructPCGSetTwoNorm(solver,1,ierr)
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
      else if ( stype == HYPRESolverPCG ) then
        call HYPRE_StructPCGSetPrecond(solver,precond_id,precond,ierr)
      end if
    end if
    if(hypre_print_level(iq) > 0) call set_print_level(stype,solver,hypre_print_level(iq))
    asolver%precond = precond
    asolver%solver  = solver
    asolver%stype   = stype
    asolver%ptype   = precond_id
  end subroutine create_solver
  !
  subroutine set_print_level(stype,solver,print_level)
    !
    ! printout of HYPRE for a solver (on the first task), e.g., to follow its convergence:
    !   1: SMG and PFMG print the residual norm of each iteration at the end of a solve, and the Krylov
    !      solvers a summary of the solve
    !   2: the Krylov solvers print the residual norm at each iteration
    ! the logging of the solver is turned on, which SMG and PFMG need for their printout; the preconditioner
    ! stays silent
    !
    implicit none
    integer    , intent(in) :: stype,print_level
    type(C_PTR), intent(in) :: solver
    select case(stype)
    case(HYPRESolverSMG)
      call HYPRE_StructSMGSetPrintLevel(solver,print_level,ierr)
      call HYPRE_StructSMGSetLogging(solver,1,ierr)
    case(HYPRESolverPFMG)
      call HYPRE_StructPFMGSetPrintLevel(solver,print_level,ierr)
      call HYPRE_StructPFMGSetLogging(solver,1,ierr)
    case(HYPRESolverGMRES)
      call HYPRE_StructGMRESSetPrintLevel(solver,print_level,ierr)
      call HYPRE_StructGMRESSetLogging(solver,1,ierr)
    case(HYPRESolverBiCGSTAB)
      call HYPRE_StructBiCGSTABSetPrintLev(solver,print_level,ierr)
      call HYPRE_StructBiCGSTABSetLogging(solver,1,ierr)
    case(HYPRESolverFlexGMRES)
      call HYPRE_StructFGMRESSetPrintLevel(solver,print_level,ierr)
      call HYPRE_StructFGMRESSetLogging(solver,1,ierr)
    case(HYPRESolverPCG)
      call HYPRE_StructPCGSetPrintLevel(solver,print_level,ierr)
      call HYPRE_StructPCGSetLogging(solver,1,ierr)
    end select
  end subroutine set_print_level
  !
  subroutine print_solve(iq,stype,solver,niter)
    !
    ! one line per solve with `hypre_print_level` > 0, the same for all solvers and versions of HYPRE
    ! (HYPRE 2 keeps the residual norms of SMG and PFMG, but does not print them): the problem, the
    ! number of iterations, and the final relative residual
    !
    implicit none
    integer    , intent(in) :: iq,stype,niter
    type(C_PTR), intent(in) :: solver
    real(rp) :: res
    res = 0._rp
    select case(stype)
    case(HYPRESolverSMG)
      call HYPRE_StructSMGGetFinalRelative(solver,res,ierr)
    case(HYPRESolverPFMG)
      call HYPRE_StructPFMGGetFinalRelativ(solver,res,ierr)
    case(HYPRESolverGMRES)
      call HYPRE_StructGMRESGetFinalRelati(solver,res,ierr)
    case(HYPRESolverBiCGSTAB)
      call HYPRE_StructBiCGSTABGetFinalRel(solver,res,ierr)
    case(HYPRESolverFlexGMRES)
      call HYPRE_StructFGMRESGetFinalRel(solver,res,ierr)
    case(HYPRESolverPCG)
      call HYPRE_StructPCGGetFinalRelative(solver,res,ierr)
    end select
    if(myid == 0) print '(A,I0,A,ES10.3)', ' HYPRE '//trim(merge('pressure ','diffusion',iq == HYPREProblemPressure))// &
                                          ' solve: ',niter,' iterations, final relative residual ',res
  end subroutine print_solve
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
    else if ( stype == HYPRESolverPCG ) then
      call HYPRE_StructPCGSetup(solver,mat,rhs,sol,ierr)
    end if
  end subroutine setup_solver
  !
  subroutine solve_helmholtz(iq,asolver,ng,lo,hi,nh_p,nh_g,dzf,p,po)
    !
    ! solves the system for the (scaled) right-hand side `p` with initial guess `po`, and returns the
    ! solution in both; warns if the solver did not converge within its maximum number of iterations
    !
    implicit none
    integer           ,         intent(in   )               :: iq
    type(hypre_solver), target, intent(in   )               :: asolver
    integer           ,         intent(in   ), dimension(3) :: ng,lo,hi
    integer           ,         intent(in   )               :: nh_p,nh_g
    real(rp)          ,         intent(in   ), dimension(lo(3)-nh_g:) :: dzf
    real(rp)          ,         intent(inout), dimension(lo(1)-nh_p:,lo(2)-nh_p:,lo(3)-nh_p:) :: p,po
    type(C_PTR), pointer :: solver,mat,rhs,sol
    integer    , pointer :: stype
    integer , dimension(asolver%ndim) :: lo_s,hi_s
    real(rp) :: pmean(2),pmean_1,psum
    integer  :: i,j,k,q,nx,ny,niter
    solver  => asolver%solver
    mat     => asolver%mat
    rhs     => asolver%rhs
    sol     => asolver%sol
    stype   => asolver%stype
    lo_s(:) = lo(asolver%idirs(1:asolver%ndim))
    hi_s(:) = hi(asolver%idirs(1:asolver%ndim))
    nx = hi(1)-lo(1)+1; ny = hi(2)-lo(2)+1
    !
    pmean_1 = 0._rp
    if(asolver%is_singular) then
      !
      ! the right-hand side of a singular problem must sum to zero: remove what round-off leaves
      !
      !$acc parallel     loop collapse(3) default(present) reduction(+:pmean_1) async(1)
      !$omp target teams loop collapse(3)                  reduction(+:pmean_1)
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            pmean_1 = pmean_1 + p(i,j,k)
          end do
        end do
      end do
      !$acc wait(1)
      call MPI_ALLREDUCE(MPI_IN_PLACE,pmean_1,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      pmean_1 = pmean_1/(1._rp*ng(1)*ng(2)*ng(3))
    end if
    !
    ! right-hand side and initial guess
    !
    !$acc parallel     loop collapse(3) default(present) private(q) async(1)
    !$omp target teams loop collapse(3)                  private(q)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          q = 1 + (i-lo(1)) + nx*((j-lo(2)) + ny*(k-lo(3)))
          rhsvalues(q) = p(i,j,k) - pmean_1
          solvalues(q) = po(i,j,k)
        end do
      end do
    end do
    !$acc wait(1)
    !$acc   host_data use_device(     rhsvalues,solvalues)
    !$omp target data use_device_addr(rhsvalues,solvalues)
    call HYPRE_StructVectorSetBoxValues(rhs,lo_s,hi_s,rhsvalues,ierr)
    call HYPRE_StructVectorAssemble(rhs,ierr)
    call HYPRE_StructVectorSetBoxValues(sol,lo_s,hi_s,solvalues,ierr)
    call HYPRE_StructVectorAssemble(sol,ierr)
    !$omp end target data
    !$acc end   host_data
    !
    ! solve
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
    else if ( stype == HYPRESolverPCG ) then
      call HYPRE_StructPCGSolve(solver,mat,rhs,sol,ierr)
      call HYPRE_StructPCGGetNumIterations(solver,niter,ierr)
    end if ! stype
    if(hypre_print_level(iq) > 0) call print_solve(iq,stype,solver,niter)
    if(niter >= hypre_maxiter(iq)) then
      if(myid == 0) print*, 'WARNING: HYPRE solver not converged after ',niter,' iterations (maximum: `hypre_maxiter`).'
    end if
    !
    ! fetch results and update old solution
    !
    !$acc   host_data use_device(     solvalues)
    !$omp target data use_device_addr(solvalues)
    call HYPRE_StructVectorGetBoxValues(sol,lo_s,hi_s,solvalues,ierr)
    !$omp end target data
    !$acc end   host_data
    call sync_hypre()
    pmean_1 = 0._rp
    if(asolver%is_singular) then
      !
      ! fix the constant of the solution: zero mean
      !
      pmean(:) = 0._rp
      !$acc parallel     loop gang default(present) private(psum) reduction(+:pmean_1) async(1)
      !$omp target teams loop                       private(psum) reduction(+:pmean_1)
      do k=lo(3),hi(3)
        psum = 0._rp
        !$acc loop collapse(2) private(q) reduction(+:psum)
        !$omp loop collapse(2) private(q) reduction(+:psum)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            q = 1 + (i-lo(1)) + nx*((j-lo(2)) + ny*(k-lo(3)))
            psum = psum + solvalues(q)
          end do
        end do
        pmean_1 = pmean_1 + psum*dzf(k)
      end do
      do k=lo(3),hi(3)
        pmean(2) = pmean(2) + (hi(1)-lo(1)+1)*(hi(2)-lo(2)+1)*dzf(k)
      end do
      !$acc wait(1)
      pmean(1) = pmean_1
      call MPI_ALLREDUCE(MPI_IN_PLACE,pmean,2,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      pmean_1 = pmean(1)/pmean(2)
    end if
    !$acc parallel     loop collapse(3) default(present) private(q) async(1)
    !$omp target teams loop collapse(3)                  private(q)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          q = 1 + (i-lo(1)) + nx*((j-lo(2)) + ny*(k-lo(3)))
          p( i,j,k) = solvalues(q) - pmean_1
          po(i,j,k) = p(i,j,k)
        end do
      end do
    end do
  end subroutine solve_helmholtz
  !
  subroutine sync_hypre()
    !
    ! waits for the work of HYPRE on the device, which uses its own stream(s), before the OpenACC queues
    ! touch the same data
    !
    implicit none
#if defined(_OPENACC) || defined(_OPENMP)
    integer :: istat
    istat = cudaDeviceSynchronize()
#endif
  end subroutine sync_hypre
  !
  subroutine finalize_solver(asolver,is_all)
    !
    ! destroys what `create_solver` creates at every solve: a multigrid solver, or the preconditioner
    ! of a Krylov solver; and the Krylov solver itself with `is_all`, otherwise kept for the next solve
    !
    implicit none
    type(hypre_solver), target, intent(inout) :: asolver
    logical           ,         intent(in   ) :: is_all
    type(C_PTR), pointer :: precond,solver
    integer    , pointer :: stype,ptype
    !
    precond => asolver%precond
    solver  => asolver%solver
    stype   => asolver%stype
    ptype   => asolver%ptype
    if      ( ptype == HYPREPrecondSMG ) then
      call HYPRE_StructSMGDestroy(precond,ierr)
    else if ( ptype == HYPREPrecondPFMG ) then
      call HYPRE_StructPFMGDestroy(precond,ierr)
    end if
    precond = C_NULL_PTR
    ptype   = HYPREPrecondNone
    if(.not.C_ASSOCIATED(solver)) return
    if      ( stype == HYPRESolverSMG ) then
      call HYPRE_StructSMGDestroy(solver,ierr)
    else if ( stype == HYPRESolverPFMG ) then
      call HYPRE_StructPFMGDestroy(solver,ierr)
    else if ( .not.is_all ) then
      return
    else if ( stype == HYPRESolverGMRES ) then
      call HYPRE_StructGMRESDestroy(solver,ierr)
    else if ( stype == HYPRESolverBiCGSTAB ) then
      call HYPRE_StructBiCGSTABDestroy(solver,ierr)
    else if ( stype == HYPRESolverFlexGMRES ) then
      call HYPRE_StructFGMRESDestroy(solver,ierr)
    else if ( stype == HYPRESolverPCG ) then
      call HYPRE_StructPCGDestroy(solver,ierr)
    end if
    solver = C_NULL_PTR
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
