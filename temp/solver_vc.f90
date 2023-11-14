module mod_solver
  use mpi_f08
  use mod_common_mpi, only: ierr
  use mod_types
  use, intrinsic :: iso_c_binding, only: C_PTR
  implicit none
  private
  public init_bc_rhs,init_matrix_3d,create_solver,setup_solver, &
         add_constant_to_diagonal,solve_helmholtz,finalize_matrix,finalize_solver, &
         hypre_solver, add_constant_to_boundary, &
         HYPRESolverSMG,HYPRESolverPFMG,HYPRESolverGMRES,HYPRESolverBiCGSTAB
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  public init_fft_reduction,init_n_2d_matrices,create_n_solvers,setup_n_solvers, &
         add_constant_to_n_diagonals,solve_n_helmholtz_2d, &
         finalize_n_matrices,finalize_n_solvers,init_comm_slab, &
         init_transpose_slab_uneven,transpose_slab,alltoallw, &
         init_n_3d_matrices,solve_n_helmholtz_3d
  type alltoallw
    integer            :: counts
    integer            :: disps
    type(MPI_DATATYPE) :: types
  end type alltoallw
#endif
!  ! parameterized derived types like below would make more sense but
!  ! they are very well supported yet
!  type alltoallw(n) ! Parameterized derived types are not working so well yet!
!    integer, len :: n
!    integer           , dimension(n) :: counts
!    integer           , dimension(n) :: disps
!    type(MPI_DATATYPE), dimension(n) :: types
!  end type alltoallw
  integer, parameter :: HYPRESolverSMG      = 1, &
                        HYPRESolverPFMG     = 2, &
                        HYPRESolverGMRES    = 3, &
                        HYPRESolverBiCGSTAB = 4
  type hypre_solver
    type(C_PTR) :: grid,stencil,precond,solver,mat,rhs,sol
    integer     :: stype,comm_hypre
  end type hypre_solver
  contains
  subroutine init_bc_rhs(cbc,bc,dl,is_bound,is_centered,lo,hi,periods, &
                         dx1,dx2,dy1,dy2,dz1,dz2,rhsx,rhsy,rhsz,bcx,bcy,bcz)
    !
    ! description
    !
    implicit none
    character(len=1)  , intent(in ), dimension(0:1,3) :: cbc
    real(rp)          , intent(in ), dimension(0:1,3) ::  bc
    real(rp)          , intent(in ), dimension(0:1,3) ::  dl
    logical           , intent(in ), dimension(0:1,3) ::  is_bound
    logical           , intent(in ), dimension(    3) ::  is_centered
    integer           , intent(in ), dimension(    3) :: lo,hi,periods
    real(rp)          , intent(in ), target, dimension(lo(1)-1:) :: dx1,dx2
    real(rp)          , intent(in ), target, dimension(lo(2)-1:) :: dy1,dy2
    real(rp)          , intent(in ), target, dimension(lo(3)-1:) :: dz1,dz2
    real(rp)          , intent(out), dimension(lo(2):,lo(3):,0:)    :: rhsx
    real(rp)          , intent(out), dimension(lo(1):,lo(3):,0:)    :: rhsy
    real(rp)          , intent(out), dimension(lo(1):,lo(2):,0:)    :: rhsz
    real(rp)          , intent(in ), dimension(lo(2):,lo(3):,0:), optional :: bcx
    real(rp)          , intent(in ), dimension(lo(1):,lo(3):,0:), optional :: bcy
    real(rp)          , intent(in ), dimension(lo(1):,lo(2):,0:), optional :: bcz
    integer, dimension(3) :: qqq
    real(rp), dimension(0:1,3) :: factor,sgn
    integer :: i,j,k,q,qq,idir,ib
    integer, dimension(3,3) :: eye
    real(rp) :: rhs
    !
    qqq(:) = 0
    where(.not.is_centered(:)) qqq(:) = 1
    factor(:,:) = 0._rp
    sgn(   :,:) = 0._rp
    do q=1,3
      do qq=0,1
        if(is_bound(qq,q)) then
          select case(cbc(qq,q))
          case('N')
            factor(qq,q) = 1._rp*dl(qq,q)*bc(qq,q)
            if(qq == 1) factor(qq,q) = -factor(qq,q)
            sgn(   qq,q) = 1._rp
          case('D')
            if(is_centered(q)) then
              factor(qq,q) = -2._rp*bc(qq,q)
              sgn(   qq,q) = -1._rp
            else
              factor(qq,q) = -1._rp*bc(qq,q)
              sgn(   qq,q) =  0._rp
            end if
          end select
        end if
      end do
    end do
    !
    eye(:,:) = 0
    do idir=1,3
      eye(idir,idir) = 1
    end do
    rhsx(:,:,:) = 0._rp
    rhsy(:,:,:) = 0._rp
    rhsz(:,:,:) = 0._rp
    do idir=1,3
      do k=lo(3),hi(3),max((hi(3)-lo(3))*eye(idir,3),1)
        do j=lo(2),hi(2),max((hi(2)-lo(1))*eye(idir,2),1)
          do i=lo(1),hi(1),max((hi(1)-lo(1))*eye(idir,1),1)
            if(periods(idir) == 0) then
              select case(idir)
              case(1)
                if(     i == lo(idir)) then
                  ib = 0
                else if(i == hi(idir)) then
                  ib = 1
                end if
                if(is_bound(ib,idir)) then
                  rhs = factor(ib,idir)/(dx1(i-(1-ib)+qqq(idir))*dx2(i))
                  if(present(bcx).and.bc(ib,idir)/=0._rp) rhs = rhs*bcx(j,k,ib)/bc(ib,idir)
                  rhsx(j,k,ib) = rhs
                end if
              case(2)
                if(     j == lo(idir)) then
                  ib = 0
                else if(j == hi(idir)) then
                  ib = 1
                end if
                if(is_bound(ib,idir)) then
                  rhs = factor(ib,idir)/(dy1(j-(1-ib)+qqq(idir))*dy2(j))
                  if(present(bcy).and.bc(ib,idir)/=0._rp) rhs = rhs*bcy(i,k,ib)/bc(ib,idir)
                  rhsy(i,k,ib) = rhs
                end if
              case(3)
                if(     k == lo(idir)) then
                  ib = 0
                else if(k == hi(idir)) then
                  ib = 1
                end if
                if(is_bound(ib,idir)) then
                  rhs = factor(ib,idir)/(dz1(k-(1-ib)+qqq(idir))*dz2(k))
                  if(present(bcz).and.bc(ib,idir)/=0._rp) rhs = rhs*bcz(i,j,ib)/bc(ib,idir)
                  rhsz(i,j,ib) = rhs
                end if
              end select
            end if
          end do
        end do
      end do
    end do
  end subroutine init_bc_rhs
  subroutine init_matrix_3d(cbc,bc,dl,is_uniform_grid,is_bound,is_centered,lo,hi,periods, &
                            dx1,dx2,dy1,dy2,dz1,dz2,alpha,alpha_bc,asolver,lambda,is_bound_outflow)
    !
    ! description
    !
    implicit none
    integer, parameter :: nstencil = 7
    character(len=1)  , intent(in ), dimension(0:1,3) :: cbc
    real(rp)          , intent(in ), dimension(0:1,3) ::  bc
    real(rp)          , intent(in ), dimension(0:1,3) ::  dl
    logical           , intent(in ), dimension(0:1,3) ::  is_bound
    logical           , intent(in )                   ::  is_uniform_grid
    logical           , intent(in ), dimension(    3) ::  is_centered
    integer           , intent(in ), dimension(    3) :: lo,hi,periods
    real(rp)          , intent(in ), target, dimension(lo(1)-1:) :: dx1,dx2
    real(rp)          , intent(in ), target, dimension(lo(2)-1:) :: dy1,dy2
    real(rp)          , intent(in ), target, dimension(lo(3)-1:) :: dz1,dz2
    real(rp)          , intent(in )                              :: alpha
    real(rp)          , intent(in ), dimension(0:1,3)            :: alpha_bc
    type(hypre_solver), intent(out)                              :: asolver
    real(rp)          , intent(in ), optional, dimension(:)      :: lambda
    logical           , intent(in ), optional, dimension(0:1,3)  :: is_bound_outflow
    integer, dimension(3         ) :: qqq
    integer, dimension(3,nstencil) :: offsets
    real(rp), dimension(product(hi(:)-lo(:)+1)*nstencil) :: matvalues
    real(rp), dimension(0:1,3) :: factor,sgn
    type(C_PTR) :: grid,stencil,mat,rhs,sol
    integer :: i,j,k,q,qq
    real(rp) :: cc,cxm,cxp,cym,cyp,czm,czp
    integer            :: comm_hypre
    logical, dimension(0:1,3) :: is_bound_outflow_aux
    !
    comm_hypre = MPI_COMM_WORLD%MPI_VAL
    qqq(:) = 0
    where(.not.is_centered(:)) qqq(:) = 1
    factor(:,:) = 0._rp
    sgn(   :,:) = 0._rp
    do q=1,3
      do qq=0,1
        if(is_bound(qq,q)) then
          select case(cbc(qq,q))
          case('N')
            factor(qq,q) = 1._rp*dl(qq,q)*bc(qq,q)
            if(qq == 1) factor(qq,q) = -factor(qq,q)
            sgn(   qq,q) = 1._rp
          case('D')
            if(is_centered(q)) then
              factor(qq,q) = -2._rp*bc(qq,q)
              sgn(   qq,q) = -1._rp
            else
              factor(qq,q) = -1._rp*bc(qq,q)
              sgn(   qq,q) =  0._rp
            end if
          end select
        end if
      end do
    end do
    if(present(is_bound_outflow)) then
      is_bound_outflow_aux(:,:) = is_bound_outflow(:,:)
    else
      is_bound_outflow_aux(:,:) = .false.
    end if
    !
    ! create 3D grid object
    !
    call HYPRE_StructGridCreate(comm_hypre,3,grid,ierr)
    call HYPRE_StructGridSetPeriodic(grid,periods,ierr)
    call HYPRE_StructGridSetExtents(grid,lo,hi,ierr)
    call HYPRE_StructGridAssemble(grid,ierr)
    !
    ! setup the finite-difference stencil
    !
    call HYPRE_StructStencilCreate(3,nstencil,stencil,ierr)
    offsets = reshape([ 0, 0, 0, &
                       -1, 0, 0, &
                        1, 0, 0, &
                        0,-1, 0, &
                        0, 1, 0, &
                        0, 0,-1, &
                        0, 0, 1 ],shape(offsets))
    do q=1,nstencil
      call HYPRE_StructStencilSetElement(stencil,q-1,offsets(:,q),ierr)
    end do
    !
    ! create coefficient matrix, and solution & right-hand-side vectors
    !
    call HYPRE_StructMatrixCreate(comm_hypre,grid,stencil,mat,ierr)
    if(is_uniform_grid) call HYPRE_StructMatrixSetSymmetric(mat,1,ierr)
    call HYPRE_StructMatrixInitialize(mat,ierr)
    call HYPRE_StructVectorCreate(comm_hypre,grid,sol,ierr)
    call HYPRE_StructVectorInitialize(sol,ierr)
    call HYPRE_StructVectorCreate(comm_hypre,grid,rhs,ierr)
    call HYPRE_StructVectorInitialize(rhs,ierr)
    q = 0
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          q = q + 1
          cxm = 1._rp/(dx1(i-1+qqq(1))*dx2(i))
          cxp = 1._rp/(dx1(i  +qqq(1))*dx2(i))
          cym = 1._rp/(dy1(j-1+qqq(2))*dy2(j))
          cyp = 1._rp/(dy1(j  +qqq(2))*dy2(j))
          czm = 1._rp/(dz1(k-1+qqq(3))*dz2(k))
          czp = 1._rp/(dz1(k  +qqq(3))*dz2(k))
          !
          alphaxp = (0.5d0*(alpha(i+1,j,k)+alpha(i,j,k)))**(-1)
          alphaxm = (0.5d0*(alpha(i-1,j,k)+alpha(i,j,k)))**(-1)
          alphayp = (0.5d0*(alpha(i,j+1,k)+alpha(i,j,k)))**(-1)
          alphaym = (0.5d0*(alpha(i,j-1,k)+alpha(i,j,k)))**(-1)
          alphazp = (0.5d0*(alpha(i,j,k+1)+alpha(i,j,k)))**(-1)
          alphazm = (0.5d0*(alpha(i,j,k-1)+alpha(i,j,k)))**(-1)
          !
          cxm = alphaxm*dli(1)**2
          cxp = alphaxp*dli(1)**2
          cym = alphaym*dli(2)**2
          cyp = alphayp*dli(2)**2
          czm = alphazm*dli(3)**2 ! change to dzci
          czp = alphazp*dli(3)**2 ! change to dzci
          cc  = -(cxm+cxp+cym+cyp+czm+czp)
          if(periods(1) == 0) then
            if(is_bound(0,1).and.i == lo(1)) then
              cc = cc + sgn(0,1)*cxm + alpha_bc(0,1)
              cxm = 0._rp
              if( is_bound_outflow_aux(0,1) ) then
                cc = cc + sgn(0,1)*cxp
                cxp = 0._rp
              end if
            end if
            if(is_bound(1,1).and.i == hi(1)) then
              cc = cc + sgn(1,1)*cxp + alpha_bc(1,1)
              cxp = 0._rp
              if( is_bound_outflow_aux(1,1) ) then
                cc = cc + sgn(1,1)*cxm
                cxm = 0._rp
              end if
            end if
          end if
          if(periods(2) == 0) then
            if(is_bound(0,2).and.j == lo(2)) then
              cc = cc + sgn(0,2)*cym + alpha_bc(0,2)
              cym = 0._rp
              if( is_bound_outflow_aux(0,2) ) then
                cc = cc + sgn(0,2)*cyp
                cyp = 0._rp
              end if
            end if
            if(is_bound(1,2).and.j == hi(2)) then
              cc = cc + sgn(1,2)*cyp + alpha_bc(1,2)
              cyp = 0._rp
              if( is_bound_outflow_aux(1,2) ) then
                cc = cc + sgn(1,2)*cym
                cym = 0._rp
              end if
            end if
          end if
          if(periods(3) == 0) then
            if(is_bound(0,3).and.k == lo(3)) then
              cc = cc + sgn(0,3)*czm + alpha_bc(0,3)
              czm = 0._rp
              if( is_bound_outflow_aux(0,3) ) then
                cc = cc + sgn(0,3)*czp
                czp = 0._rp
              end if
            end if
            if(is_bound(1,3).and.k == hi(3)) then
              cc = cc + sgn(1,3)*czp + alpha_bc(1,3)
              czp = 0._rp
              if( is_bound_outflow_aux(1,3) ) then
                cc = cc + sgn(1,3)*czm
                czm = 0._rp
              end if
            end if
          end if
          qq = (q-1)*nstencil
          matvalues(qq+1) = cc
          matvalues(qq+2) = cxm
          matvalues(qq+3) = cxp
          matvalues(qq+4) = cym
          matvalues(qq+5) = cyp
          matvalues(qq+6) = czm
          matvalues(qq+7) = czp
          if(periods(1).eq.0) then
            if(    ii.eq.  1) then
              rhs = rhs + cxm*factor(0,1)
              cc = cc + sgn(0,1)*cxm
              cxm = 0.d0
            elseif(ii.eq.ng(1)) then
              rhs = rhs + cxp*factor(1,1)
              cc = cc + sgn(1,1)*cxp
              cxp = 0.d0
            endif
          endif
          if(periods(2).eq.0) then
            if(    jj.eq.  1 ) then
              rhs = rhs + cym*factor(0,2)
              cc = cc + sgn(0,2)*cym
              cym = 0.d0
            elseif(jj.eq.ng(2)) then
              rhs = rhs + cyp*factor(1,2)
              cc = cc + sgn(1,2)*cyp
              cyp = 0.d0
            endif
          endif
          if(periods(3).eq.0) then
            if(    kk.eq.  1) then
              rhs = rhs + czm*factor(0,3)
              cc = cc + sgn(0,3)*czm
              czm = 0.d0
            elseif(kk.eq.ng(3)) then
              rhs = rhs + czp*factor(1,3)
              cc = cc + sgn(1,3)*czp
              czp = 0.d0
            endif
          endif
        end do
      end do
    end do
    call HYPRE_StructMatrixSetBoxValues(mat,lo,hi,nstencil, &
                                        [0,1,2,3,4,5,6],matvalues,ierr)
    call HYPRE_StructMatrixAssemble(mat,ierr)
    asolver%grid       = grid
    asolver%stencil    = stencil
    asolver%mat        = mat
    asolver%rhs        = rhs
    asolver%sol        = sol
    asolver%comm_hypre = comm_hypre
  end subroutine init_matrix_3d
  subroutine create_solver(maxiter,maxerror,stype,asolver)
    implicit none
    integer           ,         intent(   in) :: maxiter
    real(rp)          ,         intent(   in) :: maxerror
    integer           ,         intent(   in) :: stype
    type(hypre_solver), target, intent(inout) :: asolver
    type(C_PTR) :: solver,precond
    integer :: precond_id
    !
    ! setup solver
    !
    ! note: this part was taken from the Paris Simulator code
    !       freely available under a GPL license
    !       http://www.ida.upmc.fr/~zaleski/paris
    !
    if      ( stype == HYPRESolverSMG ) then
      call HYPRE_StructSMGCreate(asolver%comm_hypre,solver,ierr)
      call HYPRE_StructSMGSetMaxIter(solver,maxiter,ierr)
      call HYPRE_StructSMGSetTol(solver,maxerror,ierr)
      call hypre_structSMGsetLogging(solver,1,ierr)
      call HYPRE_StructSMGSetPrintLevel(solver,1,ierr)
    else if ( stype == HYPRESolverPFMG ) then
      call HYPRE_StructPFMGCreate(asolver%comm_hypre,solver,ierr)
      call HYPRE_StructPFMGSetMaxIter(solver,maxiter,ierr)
      call HYPRE_StructPFMGSetTol(solver,maxerror,ierr)
      call HYPRE_structPFMGsetLogging(solver,1,ierr)
      call HYPRE_StructPFMGSetPrintLevel(solver,1,ierr)
      call HYPRE_StructPFMGSetRelChange(solver,1,ierr)
      ! Relaxiation Method: 2 is the fastest if symm matrix
      ! 0: Jacobi
      ! 1: Weighted Jacobi (default)
      ! 2: Red/Black Gauss-Seidel (symmetric: RB pre- and post-relaxation)
      ! 3: Red/Black Gauss-Seidel (nonsymmetric: RB pre- and post-relaxation)
      call HYPRE_StructPFMGSetRelaxType(solver,1,ierr)
      call HYPRE_StructPFMGSetNumPreRelax(solver,1,ierr)
      call HYPRE_StructPFMGSetNumPostRelax(solver,1,ierr)
    else if ( stype == HYPRESolverGMRES .or. &
              stype == HYPRESolverBiCGSTAB   ) then
      if      ( stype == HYPRESolverGMRES ) then
        call HYPRE_StructGMRESCreate(asolver%comm_hypre,solver,ierr)
        call HYPRE_StructGMRESSetMaxIter(solver,maxiter,ierr)
        call HYPRE_StructGMRESSetTol(solver,maxerror,ierr)
        !call HYPRE_StructGMRESSetLogging(solver, 1 ,ierr)
      else if ( stype == HYPRESolverBiCGSTAB ) then
        call HYPRE_StructBiCGSTABCreate(asolver%comm_hypre,solver,ierr)
        call HYPRE_StructBiCGSTABSetMaxIter(solver,maxiter,ierr)
        call HYPRE_StructBiCGSTABSetTol(solver,maxerror,ierr)
      end if
      ! Use PFMG as preconditioner
      call HYPRE_StructPFMGCreate(asolver%comm_hypre,precond,ierr)
      call HYPRE_StructPFMGSetMaxIter(precond,10,ierr)
      call HYPRE_StructPFMGSetTol(precond,0._rp,ierr)
      call HYPRE_StructPFMGSetZeroGuess(precond,ierr)
      call HYPRE_StructPFMGSetRelChange(precond,1,ierr)
      call HYPRE_StructPFMGSetRelaxType(precond,2,ierr)
      precond_id = 1   ! Set PFMG as preconditioner
      if      ( stype == HYPRESolverGMRES ) then
        call HYPRE_StructGMRESSetPrecond(solver,precond_id,precond,ierr)
      else if ( stype == HYPRESolverBiCGSTAB ) then
        call HYPRE_StructBiCGSTABSetPrecond(solver,precond_id,precond,ierr)
      end if
      asolver%precond = precond
    end if
    asolver%solver  = solver
    asolver%stype   = stype
  end subroutine create_solver
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
    else if ( stype == HYPRESolverGMRES .or. &
             stype == HYPRESolverBiCGSTAB   ) then
      if      ( stype == HYPRESolverGMRES ) then
        call HYPRE_StructGMRESSetup(solver,mat,rhs,sol,ierr)
      else if ( stype == HYPRESolverBiCGSTAB ) then
        call HYPRE_StructBiCGSTABSetup(solver,mat,rhs,sol,ierr)
      end if
    end if
  end subroutine setup_solver
  subroutine add_constant_to_diagonal(lo,hi,alpha,mat)
    implicit none
    integer    , intent(in   ), dimension(3) :: lo,hi
    real(rp)   , intent(in   ) :: alpha
    type(C_PTR), intent(inout) :: mat
    real(rp), dimension(product(hi(:)-lo(:)+1)) :: matvalues
    integer :: i,j,k,q
    q = 0
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          q=q+1
          matvalues(q) = alpha
        end do
      end do
    end do
    call HYPRE_StructMatrixAddToBoxValues(mat,lo,hi,1,[0],matvalues,ierr)
  end subroutine add_constant_to_diagonal
  subroutine solve_helmholtz(asolver,lo,hi,p,po)
    implicit none
    type(hypre_solver), target, intent(in   )               :: asolver
    integer           ,         intent(in   ), dimension(3) :: lo,hi
    real(rp)          ,         intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p,po
    type(C_PTR), pointer :: solver,mat,rhs,sol
    integer    , pointer :: stype
    solver  => asolver%solver
    mat     => asolver%mat
    rhs     => asolver%rhs
    sol     => asolver%sol
    stype   => asolver%stype
    !
    call HYPRE_StructVectorSetBoxValues(rhs,lo,hi,p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),ierr)
    call HYPRE_StructVectorAssemble(rhs,ierr)
    !
    ! create soluction vector
    !
    call HYPRE_StructVectorSetBoxValues(sol,lo,hi,po(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),ierr)
    call HYPRE_StructVectorAssemble(sol,ierr)
    !
    ! setup solver, and solve
    !
    ! note: this part was based on the the Paris Simulator code
    !       freely available under a GPL license; see:
    !       http://www.ida.upmc.fr/~zaleski/paris/
    !
    if      ( stype == HYPRESolverSMG ) then
      call HYPRE_StructSMGSolve(solver,mat,rhs,sol,ierr)
      !call HYPRE_StructSMGGetNumIterations(solver,num_iterations,ierr)
    else if ( stype == HYPRESolverPFMG ) then
      call HYPRE_StructPFMGSolve(solver,mat,rhs,sol,ierr)
      !call HYPRE_StructPFMGGetNumIteration(solver,num_iterations,ierr)
    else if ( stype == HYPRESolverGMRES ) then
      call HYPRE_StructGMRESSolve(solver,mat,rhs,sol,ierr)
      !call HYPRE_StructGMRESGetNumIteratio(solver, num_iterations,ierr)
    else if ( stype == HYPRESolverBiCGSTAB ) then
      call HYPRE_StructBiCGSTABSolve(solver,mat,rhs,sol,ierr)
      !call HYPRE_StructBiCGSTABGetNumItera(solver, num_iterations,ierr)
    end if ! stype
    !
    ! end of part based on the Paris Simulator code
    !
    ! fecth results
    !
    call HYPRE_StructVectorGetBoxValues(sol,lo,hi,p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),ierr)
    !$OMP PARALLEL WORKSHARE
    po(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    !$OMP END PARALLEL WORKSHARE
  end subroutine solve_helmholtz
  subroutine finalize_solver(asolver)
    implicit none
    type(hypre_solver), target, intent(in) :: asolver
    type(C_PTR), pointer :: precond,solver
    integer    , pointer :: stype
    !
    precond => asolver%precond
    solver  => asolver%solver
    stype   => asolver%stype
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
      call HYPRE_StructPFMGDestroy(precond,ierr)
    else if ( stype == HYPRESolverBiCGSTAB ) then
      call HYPRE_StructBiCGSTABDestroy(solver,ierr)
      call HYPRE_StructPFMGDestroy(precond,ierr)
    end if
  end subroutine finalize_solver
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
  subroutine add_constant_to_boundary(lo,hi,is_bound,alpha_bc,mat)
    implicit none
    integer    , intent(in ), dimension(    3) :: lo,hi
    logical    , intent(in ), dimension(0:1,3) :: is_bound
    real(rp)   , intent(in ), dimension(0:1,3) :: alpha_bc
    type(C_PTR), intent(inout) :: mat
    real(rp), dimension(product(hi(:)-lo(:)+1)) :: matvalues
    integer :: i,j,k,q
    real(rp) :: cc
    q = 0
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          q=q+1
          cc = 0._rp
          if(is_bound(0,1).and.i == lo(1)) then
            cc = cc + alpha_bc(0,1)
          end if
          if(is_bound(1,1).and.i == hi(1)) then
            cc = cc + alpha_bc(1,1)
          end if
          if(is_bound(0,2).and.j == lo(2)) then
            cc = cc + alpha_bc(0,2)
          end if
          if(is_bound(1,2).and.j == hi(2)) then
            cc = cc + alpha_bc(1,2)
          end if
          if(is_bound(0,3).and.k == lo(3)) then
            cc = cc + alpha_bc(0,3)
          end if
          if(is_bound(1,3).and.k == hi(3)) then
            cc = cc + alpha_bc(1,3)
          end if
          matvalues(q) = cc
        end do
      end do
    end do
    call HYPRE_StructMatrixAddToBoxValues(mat,lo,hi,1,[0],matvalues,ierr)
  end subroutine add_constant_to_boundary
  subroutine add_constant_to_n_boundaries(n,lo,hi,is_bound,alpha_bc,mat)
    implicit none
    integer    , intent(in   )               :: n
    integer    , intent(in   ), dimension(2) :: lo,hi
    logical    , intent(in   ), dimension(0:1,2) :: is_bound
    real(rp)   , intent(in   ), dimension(0:1,2) :: alpha_bc
    type(C_PTR), intent(inout), dimension(:) :: mat
    integer :: q
    logical , dimension(0:1,3) :: is_bound_aux
    real(rp), dimension(0:1,3) :: alpha_bc_aux
    alpha_bc_aux(:,1:2) = alpha_bc(:,:)
    alpha_bc_aux(:,  3) = 0._rp
    is_bound_aux(:,1:2) = is_bound(:,:)
    is_bound_aux(:,  3) = .false.
    do q=1,n
      call add_constant_to_boundary([lo(1),lo(2),1],[hi(1),hi(2),1],is_bound_aux,alpha_bc_aux,mat(q))
    end do
  end subroutine add_constant_to_n_boundaries
  subroutine add_constant_to_n_3d_boundaries(nslices,lo_sp,hi_sp,is_bound,alpha_bc,mat)
    implicit none
    integer, intent(in)                       :: nslices
    integer, intent(in), dimension(3,nslices) :: lo_sp,hi_sp
    logical    , intent(in   ), dimension(0:1,3) :: is_bound
    real(rp)   , intent(in   ), dimension(0:1,3) :: alpha_bc
    type(C_PTR), intent(inout), dimension(:) :: mat
    integer :: q
    do q=1,nslices
      call add_constant_to_boundary(lo_sp(:,q),hi_sp(:,q),is_bound,alpha_bc,mat(q))
    end do
  end subroutine add_constant_to_n_3d_boundaries
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  subroutine init_matrix_2d(cbc,bc,dl,is_uniform_grid,is_bound,is_centered,lo,hi,periods, &
                            dl1_1,dl1_2,dl2_1,dl2_2,alpha,alpha_bc,comm,asolver,is_bound_outflow)
    !
    ! description
    !
    implicit none
    integer, parameter :: nstencil = 5
    character(len=1)  , intent(in ), dimension(0:1,2) :: cbc
    real(rp)          , intent(in ), dimension(0:1,2) ::  bc
    real(rp)          , intent(in ), dimension(0:1,2) ::  dl
    logical           , intent(in )                   ::  is_uniform_grid
    logical           , intent(in ), dimension(0:1,2) ::  is_bound
    logical           , intent(in ), dimension(    2) ::  is_centered
    integer           , intent(in ), dimension(    2) :: lo,hi,periods
    real(rp)          , intent(in ), target, dimension(lo(1)-1:) :: dl1_1,dl1_2
    real(rp)          , intent(in ), target, dimension(lo(2)-1:) :: dl2_1,dl2_2
    real(rp)          , intent(in )                              :: alpha
    real(rp)          , intent(in ), dimension(0:1,2)            :: alpha_bc
    type(MPI_COMM)    , intent(in )                   :: comm
    type(hypre_solver), intent(out)                              :: asolver
    logical           , intent(in ), optional, dimension(0:1,2)  :: is_bound_outflow
    integer, dimension(2         ) :: qqq
    integer, dimension(2,nstencil) :: offsets
    real(rp), dimension(product(hi(:)-lo(:)+1)*nstencil) :: matvalues
    real(rp), dimension(0:1,2) :: factor,sgn
    type(C_PTR) :: grid,stencil,mat,rhs,sol
    integer :: i1,i2,q,qq
    real(rp) :: cc,c1m,c1p,c2m,c2p
    integer            :: comm_hypre
    logical, dimension(0:1,2) :: is_bound_outflow_aux
    !
    comm_hypre = comm%MPI_VAL
    !
    qqq(:) = 0
    where(.not.is_centered(:)) qqq(:) = 1
    factor(:,:) = 0._rp
    sgn(   :,:) = 0._rp
    do q=1,2
      do qq=0,1
        if(is_bound(qq,q)) then
          select case(cbc(qq,q))
          case('N')
            factor(qq,q) = 1._rp*dl(qq,q)*bc(qq,q)
            if(qq == 1) factor(qq,q) = -factor(qq,q)
            sgn(   qq,q) = 1._rp
          case('D')
            if(is_centered(q)) then
              factor(qq,q) = -2._rp*bc(qq,q)
              sgn(   qq,q) = -1._rp
            else
              factor(qq,q) = -1._rp*bc(qq,q)
              sgn(   qq,q) =  0._rp
            end if
          end select
        end if
      end do
    end do
    if(present(is_bound_outflow)) then
      is_bound_outflow_aux(:,:) = is_bound_outflow(:,:)
    else
      is_bound_outflow_aux(:,:) = .false.
    end if
    !
    ! create 2D grid object
    !
    call HYPRE_StructGridCreate(comm_hypre,2,grid,ierr)
    call HYPRE_StructGridSetPeriodic(grid,periods,ierr)
    call HYPRE_StructGridSetExtents(grid,lo,hi,ierr)
    call HYPRE_StructGridAssemble(grid,ierr)
    !
    ! setup the finite-difference stencil
    !
    call HYPRE_StructStencilCreate(2,nstencil,stencil,ierr)
    offsets = reshape([ 0, 0, &
                       -1, 0, &
                        1, 0, &
                        0,-1, &
                        0, 1 ],shape(offsets))
    do q=1,nstencil
      call HYPRE_StructStencilSetElement(stencil,q-1,offsets(:,q),ierr)
    end do
    !
    ! create coefficient matrix, and solution & right-hand-side vectors
    !
    call HYPRE_StructMatrixCreate(comm_hypre,grid,stencil,mat,ierr)
    if(is_uniform_grid) call HYPRE_StructMatrixSetSymmetric(mat,1,ierr)
    call HYPRE_StructMatrixInitialize(mat,ierr)
    call HYPRE_StructVectorCreate(comm_hypre,grid,sol,ierr)
    call HYPRE_StructVectorInitialize(sol,ierr)
    call HYPRE_StructVectorCreate(comm_hypre,grid,rhs,ierr)
    call HYPRE_StructVectorInitialize(rhs,ierr)
    q = 0
    do i2=lo(2),hi(2)
      do i1=lo(1),hi(1)
        c1m = 1._rp/(dl1_1(i1-1+qqq(1))*dl1_2(i1))
        c1p = 1._rp/(dl1_1(i1  +qqq(1))*dl1_2(i1))
        c2m = 1._rp/(dl2_1(i2-1+qqq(2))*dl2_2(i2))
        c2p = 1._rp/(dl2_1(i2  +qqq(2))*dl2_2(i2))
        cc  = -(c1m+c1p+c2m+c2p) + alpha
        if(periods(1) == 0) then
          if(is_bound(0,1).and.i1 == lo(1)) then
            cc = cc + sgn(0,1)*c1m + alpha_bc(0,1)
            c1m = 0._rp
            if( is_bound_outflow_aux(0,1) ) then
              cc = cc + sgn(0,1)*c1p
              c1p = 0._rp
            end if
          end if
          if(is_bound(1,1).and.i1 == hi(1)) then
            cc = cc + sgn(1,1)*c1p + alpha_bc(1,1)
            c1p = 0._rp
            if( is_bound_outflow_aux(1,1) ) then
              cc = cc + sgn(1,1)*c1m
              c1m = 0._rp
            end if
          end if
        end if
        if(periods(2) == 0) then
          if(is_bound(0,2).and.i2 == lo(2)) then
            cc = cc + sgn(0,2)*c2m + alpha_bc(0,2)
            c2m = 0._rp
            if( is_bound_outflow_aux(0,2) ) then
              cc = cc + sgn(0,2)*c2p
              c2p = 0._rp
            end if
          end if
          if(is_bound(1,2).and.i2 == hi(2)) then
            cc = cc + sgn(1,2)*c2p + alpha_bc(1,2)
            c2p = 0._rp
            if( is_bound_outflow_aux(1,2) ) then
              cc = cc + sgn(1,2)*c2m
              c2m = 0._rp
            end if
          end if
        end if
        q  = q + 1
        qq = (q-1)*nstencil
        matvalues(qq+1) = cc
        matvalues(qq+2) = c1m
        matvalues(qq+3) = c1p
        matvalues(qq+4) = c2m
        matvalues(qq+5) = c2p
      end do
    end do
    call HYPRE_StructMatrixSetBoxValues(mat,lo,hi,nstencil, &
                                        [0,1,2,3,4],matvalues,ierr)
    call HYPRE_StructMatrixAssemble(mat,ierr)
    asolver%grid       = grid
    asolver%stencil    = stencil
    asolver%mat        = mat
    asolver%rhs        = rhs
    asolver%sol        = sol
    asolver%comm_hypre = comm_hypre
  end subroutine init_matrix_2d
  subroutine solve_n_helmholtz_2d(asolver,lo_out,hi_out,nh,lo,hi,p,po)
    implicit none
    type(hypre_solver), target, intent(inout), dimension(:) :: asolver
    integer           ,         intent(in   )               :: lo_out,hi_out,nh
    integer           ,         intent(in   ), dimension(2) :: lo,hi
#ifdef _FFT_Z
    real(rp)          ,         intent(inout), dimension(lo(1)-nh:,lo(2)-nh:,lo_out-nh:) :: p,po
#elif  _FFT_Y
    real(rp)          ,         intent(inout), dimension(lo(1)-nh:,lo_out-nh:,lo(2)-nh:) :: p,po
#elif  _FFT_X
    real(rp)          ,         intent(inout), dimension(lo_out-nh:,lo(1)-nh:,lo(2)-nh:) :: p,po
#endif
    type(C_PTR), pointer :: solver,mat,rhs,sol
    integer    , pointer :: stype
    integer :: i_out,n
    do i_out=lo_out,hi_out
      n = i_out-lo_out+1
      solver  => asolver(n)%solver
      mat     => asolver(n)%mat
      rhs     => asolver(n)%rhs
      sol     => asolver(n)%sol
      stype   => asolver(n)%stype
      !
      ! setup solution and rhs vectors
      !
#ifdef _FFT_Z
      call HYPRE_StructVectorSetBoxValues(rhs,lo,hi, p(lo(1):hi(1),lo(2):hi(2),i_out),ierr)
      call HYPRE_StructVectorSetBoxValues(sol,lo,hi,po(lo(1):hi(1),lo(2):hi(2),i_out),ierr)
#elif  _FFT_Y
      call HYPRE_StructVectorSetBoxValues(rhs,lo,hi, p(lo(1):hi(1),i_out,lo(2):hi(2)),ierr)
      call HYPRE_StructVectorSetBoxValues(sol,lo,hi,po(lo(1):hi(1),i_out,lo(2):hi(2)),ierr)
#elif  _FFT_X
      call HYPRE_StructVectorSetBoxValues(rhs,lo,hi, p(i_out,lo(1):hi(1),lo(2):hi(2)),ierr)
      call HYPRE_StructVectorSetBoxValues(sol,lo,hi,po(i_out,lo(1):hi(1),lo(2):hi(2)),ierr)
#endif
      call HYPRE_StructVectorAssemble(rhs,ierr)
      call HYPRE_StructVectorAssemble(sol,ierr)
      !
      ! setup solver, and solve
      !
      ! note: this part was based on the the Paris Simulator code
      !       freely available under a GPL license; see:
      !       http://www.ida.upmc.fr/~zaleski/paris/
      !
      if      ( stype == HYPRESolverSMG ) then
        call HYPRE_StructSMGSolve(solver,mat,rhs,sol,ierr)
        !call HYPRE_StructSMGGetNumIterations(solver,num_iterations,ierr)
      else if ( stype == HYPRESolverPFMG ) then
        call HYPRE_StructPFMGSolve(solver,mat,rhs,sol,ierr)
        !call HYPRE_StructPFMGGetNumIteration(solver,num_iterations,ierr)
      else if ( stype == HYPRESolverGMRES ) then
        call HYPRE_StructGMRESSolve(solver,mat,rhs,sol,ierr)
        !call HYPRE_StructGMRESGetNumIteratio(solver, num_iterations,ierr)
      else if ( stype == HYPRESolverBiCGSTAB ) then
        call HYPRE_StructBiCGSTABSolve(solver,mat,rhs,sol,ierr)
        !call HYPRE_StructBiCGSTABGetNumItera(solver, num_iterations,ierr)
      end if ! stype
      !
      ! end of part based on the Paris Simulator code
      !
      ! fecth results
      !
#ifdef _FFT_Z
      call HYPRE_StructVectorGetBoxValues(sol,lo,hi,p(lo(1):hi(1),lo(2):hi(2),i_out),ierr)
      !$OMP PARALLEL WORKSHARE
      po(lo(1):hi(1),lo(2):hi(2),i_out) = p(lo(1):hi(1),lo(2):hi(2),i_out)
      !$OMP END PARALLEL WORKSHARE
#elif  _FFT_Y
      call HYPRE_StructVectorGetBoxValues(sol,lo,hi,p(lo(1):hi(1),i_out,lo(2):hi(2)),ierr)
      !$OMP PARALLEL WORKSHARE
      po(lo(1):hi(1),i_out,lo(2):hi(2)) = p(lo(1):hi(1),i_out,lo(2):hi(2))
      !$OMP END PARALLEL WORKSHARE
#elif  _FFT_X
      call HYPRE_StructVectorGetBoxValues(sol,lo,hi,p(i_out,lo(1):hi(1),lo(2):hi(2)),ierr)
      !$OMP PARALLEL WORKSHARE
      po(i_out,lo(1):hi(1),lo(2):hi(2)) = p(i_out,lo(1):hi(1),lo(2):hi(2))
      !$OMP END PARALLEL WORKSHARE
#endif
    end do
  end subroutine solve_n_helmholtz_2d
  subroutine solve_n_helmholtz_3d(asolver,nslices,lo_sp,hi_sp,nh,lo,hi,p,po)
    implicit none
    type(hypre_solver), target, intent(inout), dimension(:) :: asolver
    integer           ,         intent(in   )               :: nslices
    integer           ,         intent(in   ), dimension(3,nslices) :: lo_sp,hi_sp
    integer           ,         intent(in   )               :: nh
    integer           ,         intent(in   ), dimension(3) :: lo,hi
    real(rp)          ,         intent(inout), dimension(lo(1)-nh:,lo(2)-nh:,lo(3)-nh:) :: p,po
    type(C_PTR), pointer :: solver,mat,rhs,sol
    integer    , pointer :: stype
    integer, dimension(3) :: lo_s,hi_s
    integer :: q
    do q=1,nslices
      solver  => asolver(q)%solver
      mat     => asolver(q)%mat
      rhs     => asolver(q)%rhs
      sol     => asolver(q)%sol
      stype   => asolver(q)%stype
      lo_s(:) = lo_sp(:,q)
      hi_s(:) = hi_sp(:,q)
      !
      ! setup soluction and rhs vectors
      !
      call HYPRE_StructVectorSetBoxValues(rhs,lo_s,hi_s, p(lo_s(1):hi_s(1),lo_s(2):hi_s(2),lo_s(3):hi_s(3)),ierr)
      call HYPRE_StructVectorSetBoxValues(sol,lo_s,hi_s,po(lo_s(1):hi_s(1),lo_s(2):hi_s(2),lo_s(3):hi_s(3)),ierr)
      call HYPRE_StructVectorAssemble(rhs,ierr)
      call HYPRE_StructVectorAssemble(sol,ierr)
      !
      ! setup solver, and solve
      !
      ! note: this part was based on the the Paris Simulator code
      !       freely available under a GPL license; see:
      !       http://www.ida.upmc.fr/~zaleski/paris/
      !
      if      ( stype == HYPRESolverSMG ) then
        call HYPRE_StructSMGSolve(solver,mat,rhs,sol,ierr)
        !call HYPRE_StructSMGGetNumIterations(solver,num_iterations,ierr)
      else if ( stype == HYPRESolverPFMG ) then
        call HYPRE_StructPFMGSolve(solver,mat,rhs,sol,ierr)
        !call HYPRE_StructPFMGGetNumIteration(solver,num_iterations,ierr)
      else if ( stype == HYPRESolverGMRES ) then
        call HYPRE_StructGMRESSolve(solver,mat,rhs,sol,ierr)
        !call HYPRE_StructGMRESGetNumIteratio(solver, num_iterations,ierr)
      else if ( stype == HYPRESolverBiCGSTAB ) then
        call HYPRE_StructBiCGSTABSolve(solver,mat,rhs,sol,ierr)
        !call HYPRE_StructBiCGSTABGetNumItera(solver, num_iterations,ierr)
      end if ! stype
      !
      ! end of part based on the Paris Simulator code
      !
      ! fecth results
      !
      call HYPRE_StructVectorGetBoxValues(sol,lo_s,hi_s,p(lo_s(1):hi_s(1),lo_s(2):hi_s(2),lo_s(3):hi_s(3)),ierr)
      !$OMP PARALLEL WORKSHARE
      po(lo_s(1):hi_s(1),lo_s(2):hi_s(2),lo_s(3):hi_s(3)) = p(lo_s(1):hi_s(1),lo_s(2):hi_S(2),lo_s(3):hi_s(3))
      !$OMP END PARALLEL WORKSHARE
    end do
  end subroutine solve_n_helmholtz_3d
  subroutine init_n_2d_matrices(cbc,bc,dl,is_uniform_grid,is_bound,is_centered,lo_out,hi_out,lo,hi,periods, &
                                dl1_1,dl1_2,dl2_1,dl2_2,alpha,alpha_bc,lambda,comm,asolver,is_bound_outflow)
    character(len=1)  , intent(in   ), dimension(0:1,2) :: cbc
    real(rp)          , intent(in   ), dimension(0:1,2) ::  bc
    real(rp)          , intent(in   ), dimension(0:1,2) ::  dl
    logical           , intent(in   )                   ::  is_uniform_grid
    logical           , intent(in   ), dimension(0:1,2) ::  is_bound
    logical           , intent(in   ), dimension(    2) ::  is_centered
    integer           , intent(in   )                   :: lo_out,hi_out
    integer           , intent(in   ), dimension(    2) :: lo,hi,periods
    real(rp)          , intent(in   ), target, dimension(lo(1)-1:) :: dl1_1,dl1_2
    real(rp)          , intent(in   ), target, dimension(lo(2)-1:) :: dl2_1,dl2_2
    real(rp)          , intent(in   )                              :: alpha
    real(rp)          , intent(in   ), dimension(0:1,2)            :: alpha_bc
    real(rp)          , intent(in   ), dimension(:) :: lambda
    type(MPI_COMM)    , intent(in   ), dimension(:) :: comm
    type(hypre_solver), intent(inout), dimension(:) :: asolver
    logical           , intent(in   ), optional, dimension(0:1,2)  :: is_bound_outflow
    type(hypre_solver) :: asolver_aux
    integer :: i_out,q
    logical, dimension(0:1,2) :: is_bound_outflow_aux
    is_bound_outflow_aux(:,:) = .false.
    if( present(is_bound_outflow) ) is_bound_outflow_aux(:,:) = is_bound_outflow(:,:)
    do i_out=lo_out,hi_out
      q = i_out-lo_out+1
      asolver_aux = asolver(q)
      call init_matrix_2d(cbc,bc,dl,is_uniform_grid,is_bound,is_centered,lo,hi,periods, &
                          dl1_1,dl1_2,dl2_1,dl2_2,alpha,alpha_bc,comm(q),asolver_aux,is_bound_outflow_aux)
      call add_constant_to_diagonal([lo(1),lo(2),1],[hi(1),hi(2),1],lambda(q),asolver_aux%mat)
      asolver(q) = asolver_aux
    end do
  end subroutine init_n_2d_matrices
  subroutine init_n_3d_matrices(idir,nslice,cbc,bc,dl,is_uniform_grid,is_bound,is_centered,lo,periods, &
                                lo_sp,hi_sp,dl1_1,dl1_2,dl2_1,dl2_2,dl3_1,dl3_2,alpha,alpha_bc,lambda,asolver, &
                                is_bound_outflow)
    integer           , intent(in   )                   :: idir,nslice
    character(len=1)  , intent(in   ), dimension(0:1,3) :: cbc
    real(rp)          , intent(in   ), dimension(0:1,3) ::  bc
    real(rp)          , intent(in   ), dimension(0:1,3) ::  dl
    logical           , intent(in   )                   ::  is_uniform_grid
    logical           , intent(in   ), dimension(0:1,3) ::  is_bound
    logical           , intent(in   ), dimension(    3) ::  is_centered
    integer           , intent(in   ), dimension(    3) :: lo,periods
    integer           , intent(in   ), dimension(:,:  ) :: lo_sp,hi_sp
    real(rp)          , intent(in   ), target, dimension(lo(1)-1:) :: dl1_1,dl1_2
    real(rp)          , intent(in   ), target, dimension(lo(2)-1:) :: dl2_1,dl2_2
    real(rp)          , intent(in   ), target, dimension(lo(3)-1:) :: dl3_1,dl3_2
    real(rp)          , intent(in   )                              :: alpha
    real(rp)          , intent(in   ), dimension(0:1,3)            :: alpha_bc
    real(rp)          , intent(in   ), dimension(:) :: lambda
    type(hypre_solver), intent(inout), dimension(:) :: asolver
    logical           , intent(in   ), optional, dimension(0:1,3)  :: is_bound_outflow
    type(hypre_solver) :: asolver_aux
    integer :: q
    logical, dimension(0:1,3) :: is_bound_outflow_aux
    is_bound_outflow_aux(:,:) = .false.
    if( present(is_bound_outflow) ) is_bound_outflow_aux(:,:) = is_bound_outflow(:,:)
    do q=1,nslice
      asolver_aux = asolver(q)
      call init_matrix_3d(cbc,bc,dl,is_uniform_grid,is_bound,is_centered,lo_sp(:,q),hi_sp(:,q),periods, &
                          dl1_1(lo_sp(1,q)-1:hi_sp(1,q)+1),dl1_2(lo_sp(1,q)-1:hi_sp(1,q)+1), &
                          dl2_1(lo_sp(2,q)-1:hi_sp(2,q)+1),dl2_2(lo_sp(2,q)-1:hi_sp(2,q)+1), &
                          dl3_1(lo_sp(3,q)-1:hi_sp(3,q)+1),dl3_2(lo_sp(3,q)-1:hi_sp(3,q)+1), &
                          alpha,alpha_bc, & 
                          asolver_aux,lambda(lo_sp(idir,q)-lo(idir)+1:hi_sp(idir,q)-lo(idir)+1), &
                          is_bound_outflow_aux)
      asolver(q) = asolver_aux
    end do
  end subroutine init_n_3d_matrices
  subroutine create_n_solvers(n,maxiter,maxerror,stype,asolver)
    integer           , intent(   in) :: n
    integer           , intent(   in) :: maxiter
    real(rp)          , intent(   in) :: maxerror
    integer           , intent(   in) :: stype
    type(hypre_solver), intent(inout), dimension(:) :: asolver
    integer :: q
    do q=1,n
      call create_solver(maxiter,maxerror,stype,asolver(q))
    end do
  end subroutine create_n_solvers
  subroutine setup_n_solvers(n,asolver)
    integer           , intent(   in) :: n
    type(hypre_solver), intent(inout), dimension(:) :: asolver
    integer :: q
    do q=1,n
      call setup_solver(asolver(q))
    end do
  end subroutine setup_n_solvers
  subroutine add_constant_to_n_diagonals(n,lo,hi,alpha,mat)
    integer    , intent(in   )               :: n
    integer    , intent(in   ), dimension(2) :: lo,hi
    real(rp)   , intent(in   )               :: alpha
    type(C_PTR), intent(inout), dimension(:) :: mat
    integer :: q
    do q=1,n
      call add_constant_to_diagonal([lo(1),lo(2),1],[hi(1),hi(2),1],alpha,mat(q))
    end do
  end subroutine add_constant_to_n_diagonals
  subroutine finalize_n_matrices(n,asolver)
    integer           , intent(   in) :: n
    type(hypre_solver), intent(inout), dimension(:) :: asolver
    integer :: q
    do q=1,n
      call finalize_matrix(asolver(q))
    end do
  end subroutine finalize_n_matrices
  subroutine finalize_n_solvers(n,asolver)
    integer           , intent(   in) :: n
    type(hypre_solver), intent(inout), dimension(:) :: asolver
    integer :: q
    do q=1,n
      call finalize_solver(asolver(q))
    end do
  end subroutine finalize_n_solvers
  subroutine init_fft_reduction(idir,n,bc,is_centered,dl,arrplan,normfft,lambda)
    use mod_fft       , only: fftini,eigenvalues
    implicit none
    integer         , intent(in )                 :: idir
    integer         , intent(in ), dimension(3  ) :: n
    character(len=1), intent(in ), dimension(0:1) :: bc
    logical         , intent(in )                 :: is_centered
    real(rp)        , intent(in )                 :: dl
    type(C_PTR)     , intent(out), dimension(2  ) :: arrplan
    real(rp)        , intent(out)                 :: normfft
    real(rp)        , intent(out), dimension(:  ) :: lambda
    call fftini(idir,n,bc,is_centered,arrplan,normfft)
    call eigenvalues(n(idir),bc,is_centered,lambda)
    lambda(:) = lambda(:)/dl**2
  end subroutine init_fft_reduction
  subroutine solve_n_helmholtz_2d_old(asolver,maxiter,maxerror,lo_out,hi_out,lo,hi,lambda,alpha_old,p,po)
    implicit none
    type(hypre_solver), target, intent(inout)               :: asolver
    integer           ,         intent(in   )               :: maxiter
    real(rp)          ,         intent(in   )               :: maxerror
    integer           ,         intent(in   )               :: lo_out,hi_out
    integer           ,         intent(in   ), dimension(2) :: lo,hi
    real(rp)          ,         intent(inout)               :: alpha_old
    real(rp)          ,         intent(in   ), dimension(:) :: lambda
#ifdef _FFT_Z
    real(rp)          ,         intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo_out-1:) :: p,po
#elif  _FFT_Y
    real(rp)          ,         intent(inout), dimension(lo(1)-1:,lo_out-1:,lo(2)-1:) :: p,po
#elif  _FFT_X
    real(rp)          ,         intent(inout), dimension(lo_out-1:,lo(1)-1:,lo(2)-1:) :: p,po
#endif
    type(C_PTR), pointer :: solver,mat,rhs,sol
    integer    , pointer :: stype
    real(rp), dimension(product(hi(:)-lo(:)+1)) :: solvalues,rhsvalues
    integer :: i_out,ii,i1,i2,q
    real(rp) :: alpha
    solver  => asolver%solver
    mat     => asolver%mat
    rhs     => asolver%rhs
    sol     => asolver%sol
    stype   => asolver%stype
    do i_out=lo_out,hi_out
      q = 0
      do i2=lo(2),hi(2)
        do i1=lo(1),hi(1)
          q = q + 1
#ifdef _FFT_Z
          rhsvalues(q) = p( i1,i2,i_out)
          solvalues(q) = po(i1,i2,i_out)
#elif  _FFT_Y
          rhsvalues(q) = p( i1,i_out,i2)
          solvalues(q) = po(i1,i_out,i2)
#elif  _FFT_X
          rhsvalues(q) = p( i_out,i1,i2)
          solvalues(q) = po(i_out,i1,i2)
#endif
        end do
      end do
      ii = i_out-lo_out+1
      alpha = lambda(ii) - alpha_old
      call add_constant_to_diagonal([lo(1),lo(2),1],[hi(1),hi(2),1],alpha,mat)
      alpha_old = lambda(ii)
      !
      call HYPRE_StructVectorSetBoxValues(rhs,lo,hi, &
                                          rhsvalues,ierr)
      call HYPRE_StructVectorAssemble(rhs,ierr)
      !
      ! create soluction vector
      !
      call HYPRE_StructVectorSetBoxValues(sol,lo,hi, &
                                          solvalues,ierr)
      call HYPRE_StructVectorAssemble(sol,ierr)
      !
      ! setup solver, and solve
      !
      ! note: this part was based on the the Paris Simulator code
      !       freely available under a GPL license; see:
      !       http://www.ida.upmc.fr/~zaleski/paris/
      !
      call create_solver(maxiter,maxerror,stype,asolver)
      call setup_solver(asolver)
      if      ( stype == HYPRESolverSMG ) then
        call HYPRE_StructSMGSolve(solver,mat,rhs,sol,ierr)
        !call HYPRE_StructSMGGetNumIterations(solver,num_iterations,ierr)
      else if ( stype == HYPRESolverPFMG ) then
        call HYPRE_StructPFMGSolve(solver,mat,rhs,sol,ierr)
        !call HYPRE_StructPFMGGetNumIteration(solver,num_iterations,ierr)
      else if ( stype == HYPRESolverGMRES ) then
        call HYPRE_StructGMRESSolve(solver,mat,rhs,sol,ierr)
        !call HYPRE_StructGMRESGetNumIteratio(solver, num_iterations,ierr)
      else if ( stype == HYPRESolverBiCGSTAB ) then
        call HYPRE_StructBiCGSTABSolve(solver,mat,rhs,sol,ierr)
        !call HYPRE_StructBiCGSTABGetNumItera(solver, num_iterations,ierr)
      end if ! stype
      !
      ! end of part based on the Paris Simulator code
      !
      ! fecth results
      !
      call HYPRE_StructVectorGetBoxValues(sol,lo,hi,solvalues,ierr)
      q = 0
      do i2=lo(2),hi(2)
        do i1=lo(1),hi(1)
          q = q + 1
#ifdef _FFT_Z
          p( i1,i2,i_out) = solvalues(q)
          po(i1,i2,i_out) = p(i1,i2,i_out)
#elif  _FFT_Y
          p( i1,i_out,i2) = solvalues(q)
          po(i1,i_out,i2) = p(i1,i_out,i2)
#elif  _FFT_X
          p( i_out,i1,i2) = solvalues(q)
          po(i_out,i1,i2) = p(i_out,i1,i2)
#endif
        end do
      end do
      call finalize_solver(asolver)
    end do
  end subroutine solve_n_helmholtz_2d_old
  subroutine init_transpose_slab_uneven(idir,nhi,nho,dims,lo_p,lo_s,n_p,n_s,comm_block,transpose_params)
    implicit none
    integer, intent(in)                            :: idir,nhi,nho
    integer, intent(in), dimension(3)              :: dims,lo_p,lo_s,n_p,n_s ! n.b. lower bounds wrt [0,0,0]
    type(MPI_COMM)     , intent(in)                :: comm_block
    type(alltoallw), intent(out), dimension(:,:) :: transpose_params
    integer, dimension(3) :: n_i
    type(MPI_DATATYPE) :: type_sub_p,type_sub_s
    integer(MPI_ADDRESS_KIND) :: lb,ext_rp
    integer, dimension(3,3) :: eye
    integer :: i,j,k,ii,jj,kk,irank
    integer, dimension(3,product(dims)) :: lo_p_a,lo_s_a,n_p_a,n_s_a,n_i_a
    !
    eye(:,:) = 0
    do i=1,3
      eye(i,i) = 1
    end do
    call MPI_ALLGATHER(lo_p,3,MPI_INTEGER,lo_p_a,3,MPI_INTEGER,comm_block)
    call MPI_ALLGATHER(lo_s,3,MPI_INTEGER,lo_s_a,3,MPI_INTEGER,comm_block)
    call MPI_ALLGATHER(n_p ,3,MPI_INTEGER,n_p_a ,3,MPI_INTEGER,comm_block)
    call MPI_ALLGATHER(n_s ,3,MPI_INTEGER,n_s_a ,3,MPI_INTEGER,comm_block)
    transpose_params(:,:)%counts = 1
    call MPI_TYPE_GET_EXTENT(MPI_REAL_RP,lb,ext_rp)
    irank = 0
    do k=0,dims(3)-1
      do j=0,dims(2)-1
        do i=0,dims(1)-1
          ii = lo_s_a(1,irank+1)*eye(1,idir) + 1
          jj = lo_s_a(2,irank+1)*eye(2,idir) + 1
          kk = lo_s_a(3,irank+1)*eye(3,idir) + 1
          transpose_params(irank+1,1)%disps = extent_f(ii,jj,kk,n_p(:)+2*nhi)*ext_rp
          ii = lo_p_a(1,irank+1)*(1-eye(1,idir)) + 1
          jj = lo_p_a(2,irank+1)*(1-eye(2,idir)) + 1
          kk = lo_p_a(3,irank+1)*(1-eye(3,idir)) + 1
          transpose_params(irank+1,2)%disps = extent_f(ii,jj,kk,n_s(:)+2*nho)*ext_rp
          !
          n_i(:) = n_s_a(:,irank+1)*(eye(:,idir)) + n_p(:)*(1-eye(:,idir))
          call MPI_TYPE_CREATE_SUBARRAY(3,n_p(:)+2*nhi,n_i(:),[0,0,0],MPI_ORDER_FORTRAN,MPI_REAL_RP,type_sub_p)
          n_i(:) = n_s(:)*(eye(:,idir)) + n_p_a(:,irank+1)*(1-eye(:,idir))
          call MPI_TYPE_CREATE_SUBARRAY(3,n_s(:)+2*nho,n_i(:),[0,0,0],MPI_ORDER_FORTRAN,MPI_REAL_RP,type_sub_s)
          call MPI_TYPE_COMMIT(type_sub_p)
          call MPI_TYPE_COMMIT(type_sub_s)
          transpose_params(irank+1,1)%types = type_sub_p
          transpose_params(irank+1,2)%types = type_sub_s
          irank = irank + 1
        end do
      end do
    end do
  end subroutine init_transpose_slab_uneven
  subroutine init_transpose_slab(idir,nhi,nho,dims,n_p,n_s,transpose_params)
    implicit none
    integer, intent(in)                            :: idir,nhi,nho
    integer, intent(in), dimension(3)              :: dims,n_p,n_s
    type(alltoallw), intent(out), dimension(:,:)   :: transpose_params
    integer, dimension(3) :: n_i
    type(MPI_DATATYPE) :: type_sub_p,type_sub_s
    integer(MPI_ADDRESS_KIND) :: lb,ext_rp
    integer, dimension(3,3) :: eye
    integer :: i,j,k,ii,jj,kk,irank
    !
    eye(:,:) = 0
    do i=1,3
      eye(i,i) = 1
    end do
    select case(idir)
      case(1)
        n_i(:) = [n_s(1),n_p(2),n_p(3)]
      case(2)
        n_i(:) = [n_p(1),n_s(2),n_p(3)]
      case(3)
        n_i(:) = [n_p(1),n_p(2),n_s(3)]
    end select
    transpose_params(:,:)%counts = 1
    call MPI_TYPE_GET_EXTENT(MPI_REAL_RP,lb,ext_rp)
    irank = 0
    do k=0,dims(3)-1
      do j=0,dims(2)-1
        do i=0,dims(1)-1
          ii = irank*n_s(1)*eye(1,idir) + 1
          jj = irank*n_s(2)*eye(2,idir) + 1
          kk = irank*n_s(3)*eye(3,idir) + 1
          transpose_params(irank+1,1)%disps = extent_f(ii,jj,kk,n_p+2*nhi)*ext_rp
          ii = i*n_p(1)*(1-eye(1,idir)) + 1
          jj = j*n_p(2)*(1-eye(2,idir)) + 1
          kk = k*n_p(3)*(1-eye(3,idir)) + 1
          transpose_params(irank+1,2)%disps = extent_f(ii,jj,kk,n_s+2*nho)*ext_rp
          irank = irank + 1
        end do
      end do
    end do
    call MPI_TYPE_CREATE_SUBARRAY(3,n_p+2*nhi,n_i,[0,0,0],MPI_ORDER_FORTRAN,MPI_REAL_RP,type_sub_p)
    call MPI_TYPE_CREATE_SUBARRAY(3,n_s+2*nho,n_i,[0,0,0],MPI_ORDER_FORTRAN,MPI_REAL_RP,type_sub_s)
    call MPI_TYPE_COMMIT(type_sub_p)
    call MPI_TYPE_COMMIT(type_sub_s)
    transpose_params(:,1)%types = type_sub_p
    transpose_params(:,2)%types = type_sub_s
  end subroutine init_transpose_slab
  pure integer function extent_f(i,j,k,n) ! memory position of array with point i,j,k assuming Fortran ordering
     implicit none
     integer, intent(in) :: i,j,k,n(3)
     extent_f = (i-1) + (j-1)*n(1) + (k-1)*n(1)*n(2)
  end function extent_f
  subroutine transpose_slab(nhi,nho,t_param,comm_block,arr_in,arr_out)
    implicit none
    integer, intent(in)                            :: nhi,nho
    type(alltoallw), dimension(:,:), intent(in)    :: t_param
    type(MPI_COMM), intent(in)                     :: comm_block
    real(rp), intent(in ), dimension(1-nhi:,1-nhi:,1-nhi:) :: arr_in
    real(rp), intent(out), dimension(1-nho:,1-nho:,1-nho:) :: arr_out
    call MPI_ALLTOALLW(arr_in( 1,1,1),t_param(:,1)%counts,t_param(:,1)%disps,t_param(:,1)%types, &
                       arr_out(1,1,1),t_param(:,2)%counts,t_param(:,2)%disps,t_param(:,2)%types, &
                       comm_block)
  end subroutine transpose_slab
  subroutine init_comm_slab(lo_idir,hi_idir,lo_s_idir,hi_s_idir,myid,comms_slab)
    implicit none
    integer, intent(in) :: lo_idir,hi_idir,lo_s_idir,hi_s_idir,myid
    type(MPI_COMM), dimension(hi_s_idir-lo_s_idir+1), intent(out) :: comms_slab
    type(MPI_COMM) :: comm
    integer :: i,n,icolor
    do i=lo_idir,hi_idir
      n = i-lo_s_idir+1
      icolor = MPI_UNDEFINED
      if( i >= lo_s_idir .and. i <= hi_s_idir) icolor=i
      ! n.b. -- order of ranks could be set by e.g. memory layout (i-1 + (j-1)*n(1) + (k-1)*n(1)*n(2)*k)
      call MPI_COMM_SPLIT(MPI_COMM_WORLD,icolor,myid,comm)
      if(i >= lo_s_idir .and. i <= hi_s_idir) comms_slab(n) = comm
    end do
  end subroutine init_comm_slab
#endif
end module mod_solver
