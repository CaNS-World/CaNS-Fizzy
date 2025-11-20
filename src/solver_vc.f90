! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 The CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_solver_vc
#if !defined(_CONSTANT_COEFFS_POISSON)
  use mpi
  use mod_types
  use, intrinsic :: iso_c_binding, only: C_PTR
  implicit none
  integer :: ierr
  private
  public solver_vc
  integer, parameter, public :: HYPRESolverSMG      = 1, &
                                HYPRESolverPFMG     = 2, &
                                HYPRESolverGMRES    = 3, &
                                HYPRESolverBiCGSTAB = 4
  type hypre_solver_t
    type(C_PTR) :: grid,stencil,precond,solver,mat,rhs,sol
    integer     :: stype,comm_hypre
  end type hypre_solver_t
  contains
  subroutine solver_vc(ng,lo,hi,cbc,bc,dli,dzci,dzfi,is_bound,alpha12,psi,p,po)
    !
    ! hypre solver parameter hard-coded for now
    !
    integer , parameter  :: maxiter  = 1000
    real(rp), parameter  :: maxerror = 1.e-9
    integer , parameter  :: stype    = HYPRESolverSMG
    type(hypre_solver_t) :: asolver
    !
    integer             , intent(in   ), dimension(3)                          :: ng,lo,hi
    character(len=1)    , intent(in   ), dimension(0:1,3)                      :: cbc
    real(rp)            , intent(in   ), dimension(0:1,3)                      :: bc
    real(rp)            , intent(in   ), dimension(3)                          :: dli
    real(rp)            , intent(in   ), dimension(lo(3)-1:)                   :: dzci,dzfi
    logical             , intent(in   ), dimension(0:1,3)                      :: is_bound
    real(rp)            , intent(in   ), dimension(2)                          :: alpha12
    real(rp)            , intent(in   ), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: psi
    real(rp)            , intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p,po ! updt_rhs_b done here!
    !
    ! iterative solver created and destroyed every time step for now
    !
    call init_matrix_3d(ng,lo,hi,cbc,bc,dli,dzci,dzfi,is_bound,alpha12,psi,p,asolver)
    call create_solver(maxiter,maxerror,stype,asolver)
    call setup_solver(asolver)
    call solve_helmholtz(asolver,lo,hi,p,po)
    call finalize_matrix(asolver)
    call finalize_solver(asolver)
  end subroutine solver_vc
  !
  subroutine init_matrix_3d(ng,lo,hi,cbc,bc,dli,dzci,dzfi,is_bound,alpha12,psi,p,asolver)
    !
    implicit none
    integer, parameter :: nstencil = 7
    integer             , intent(in   ), dimension(3)                          :: ng,lo,hi
    character(len=1)    , intent(in   ), dimension(0:1,3)                      :: cbc
    real(rp)            , intent(in   ), dimension(0:1,3)                      :: bc
    real(rp)            , intent(in   ), dimension(3)                          :: dli
    real(rp)            , intent(in   ), dimension(lo(3)-1:)                   :: dzci,dzfi
    logical             , intent(in   ), dimension(0:1,3)                      :: is_bound
    real(rp)            , intent(in   ), dimension(2)                          :: alpha12
    real(rp)            , intent(in   ), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: psi
    real(rp)            , intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p ! updt_rhs_b done here!
    type(hypre_solver_t), intent(out  )                                        :: asolver
    integer, dimension(3,nstencil) :: offsets
    real(rp), dimension(product(hi(:)-lo(:)+1)*nstencil) :: matvalues
    real(rp), dimension(0:1,3) :: sgn,factor
    integer , dimension(3) :: periods
    type(C_PTR) :: grid,stencil,mat,rhs,sol
    real(rp) :: alphaxp,alphaxm,alphayp,alphaym,alphazp,alphazm
    real(rp) :: cc,cxm,cxp,cym,cyp,czm,czp
    real(rp) :: alpha,dalpha
    integer :: i,j,k,q,qq
    integer :: comm_hypre
    !
    comm_hypre = MPI_COMM_WORLD
    !
    sgn(:,:) = 0._rp
    do q=1,3
      do qq=0,1
        if(is_bound(qq,q)) then
          select case(cbc(qq,q))
          case('N')
            sgn(qq,q)    =  1._rp
            factor(qq,q) = 1.d0/dli(q)*bc(qq,q) ! N.B.: only valid for constant grid spacing
          case('D')
            sgn(qq,q)    = -1._rp
            factor(qq,q) = -2.d0*bc(qq,q)
          end select
        end if
      end do
    end do
    periods(:) = 0
    where(cbc(0,:)//cbc(1,:) == 'PP') periods(:) = ng(:)
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
    call HYPRE_StructMatrixInitialize(mat,ierr)
    call HYPRE_StructVectorCreate(comm_hypre,grid,sol,ierr)
    call HYPRE_StructVectorInitialize(sol,ierr)
    call HYPRE_StructVectorCreate(comm_hypre,grid,rhs,ierr)
    call HYPRE_StructVectorInitialize(rhs,ierr)
    !
    alpha = alpha12(2); dalpha = alpha12(1)-alpha12(2)
    q = 0
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          q = q + 1
          alphaxp = (alpha + dalpha*0.5*(psi(i+1,j,k)+psi(i,j,k)))**(-1)
          alphaxm = (alpha + dalpha*0.5*(psi(i-1,j,k)+psi(i,j,k)))**(-1)
          alphayp = (alpha + dalpha*0.5*(psi(i,j+1,k)+psi(i,j,k)))**(-1)
          alphaym = (alpha + dalpha*0.5*(psi(i,j-1,k)+psi(i,j,k)))**(-1)
          alphazp = (alpha + dalpha*0.5*(psi(i,j,k+1)+psi(i,j,k)))**(-1)
          alphazm = (alpha + dalpha*0.5*(psi(i,j,k-1)+psi(i,j,k)))**(-1)
          !
          cxm = dli(1)*dli(1)
          cxp = dli(1)*dli(1)
          cym = dli(2)*dli(2)
          cyp = dli(2)*dli(2)
          czm = dzci(k-1)*dzfi(k)
          czp = dzci(k  )*dzfi(k)
          !
          cxm = cxm*alphaxm
          cxp = cxp*alphaxp
          cym = cym*alphaym
          cyp = cyp*alphayp
          czm = czm*alphazm
          czp = czp*alphazp
          cc  = -(cxm+cxp+cym+cyp+czm+czp)
          if(periods(1) == 0) then
            if(is_bound(0,1).and.i == lo(1)) then
              cc = cc + sgn(0,1)*cxm
              p(i,j,k) = p(i,j,k) + cxm*factor(0,1)
              cxm = 0.
            end if
            if(is_bound(1,1).and.i == hi(1)) then
              cc = cc + sgn(1,1)*cxp
              p(i,j,k) = p(i,j,k) + cxp*factor(1,1)
              cxp = 0.
            end if
          end if
          if(periods(2) == 0) then
            if(is_bound(0,2).and.j == lo(2)) then
              cc = cc + sgn(0,2)*cym
              p(i,j,k) = p(i,j,k) + cym*factor(0,2)
              cym = 0.
            end if
            if(is_bound(1,2).and.j == hi(2)) then
              cc = cc + sgn(1,2)*cyp
              p(i,j,k) = p(i,j,k) + cyp*factor(1,2)
              cyp = 0.
            end if
          end if
          if(periods(3) == 0) then
            if(is_bound(0,3).and.k == lo(3)) then
              cc = cc + sgn(0,3)*czm
              p(i,j,k) = p(i,j,k) + czm*factor(0,3)
              czm = 0.
            end if
            if(is_bound(1,3).and.k == hi(3)) then
              cc = cc + sgn(1,3)*czp
              p(i,j,k) = p(i,j,k) + czp*factor(1,3)
              czp = 0.
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
  !
  subroutine create_solver(maxiter,maxerror,stype,asolver)
    implicit none
    integer             ,         intent(   in) :: maxiter
    real(rp)            ,         intent(   in) :: maxerror
    integer             ,         intent(   in) :: stype
    type(hypre_solver_t), target, intent(inout) :: asolver
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
      !call hypre_structSMGsetLogging(solver,1,ierr)
      !call HYPRE_StructSMGSetPrintLevel(solver,1,ierr)
    else if ( stype == HYPRESolverPFMG ) then
      call HYPRE_StructPFMGCreate(asolver%comm_hypre,solver,ierr)
      call HYPRE_StructPFMGSetMaxIter(solver,maxiter,ierr)
      call HYPRE_StructPFMGSetTol(solver,maxerror,ierr)
      !call HYPRE_structPFMGsetLogging(solver,1,ierr)
      !call HYPRE_StructPFMGSetPrintLevel(solver,1,ierr)
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
  !
  subroutine setup_solver(asolver)
    implicit none
    type(hypre_solver_t), target, intent(inout) :: asolver
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
  !
  subroutine solve_helmholtz(asolver,lo,hi,p,po)
    implicit none
    type(hypre_solver_t), target, intent(in   )               :: asolver
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
    ! fecth results and update old solution
    !
    call HYPRE_StructVectorGetBoxValues(sol,lo,hi,p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),ierr)
    po(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  end subroutine solve_helmholtz
  !
  subroutine finalize_solver(asolver)
    implicit none
    type(hypre_solver_t), target, intent(in) :: asolver
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
  !
  subroutine finalize_matrix(asolver)
    implicit none
    type(hypre_solver_t), target, intent(in) :: asolver
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
