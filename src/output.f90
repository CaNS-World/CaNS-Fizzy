! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_output
  use mpi
  use decomp_2d_io
  use mod_common_mpi, only:ierr,myid
  use mod_load      , only: io_write_subset
  use mod_types
  use mod_posix     , only: link_file,rename_file,remove_file
  implicit none
  private
  public is_output_due, &
         out0d,gen_alias,out1d,out1d_chan,out2d,out3d,write_log_output,write_visu_2d,write_visu_3d,write_visu_subset, &
         velstats,cmpt_total_area,cmpt_mean_pos_and_vel,cmpt_mean_mass,cmpt_mean_energy
  character(len=*), parameter :: fmt_dp = '(*(es24.16e3,1x))', &
                                 fmt_sp = '(*(es15.8e2,1x))'
#if !defined(_SINGLE_PRECISION)
  character(len=*), parameter :: fmt_rp = fmt_dp
#else
  character(len=*), parameter :: fmt_rp = fmt_sp
#endif
  contains
  pure logical function is_output_due(istep,interval,period,time_start,time_end,restart) result(is_due)
    !
    ! whether an output is due at the end of a time step: with `period > 0`, if a multiple of `period`
    ! was crossed (at startup, with `restart` present: for a new run only); otherwise, if `istep` is a
    ! multiple of `interval`
    !
    implicit none
    integer , intent(in) :: istep,interval
    real(rp), intent(in) :: period,time_start,time_end
    logical , intent(in), optional :: restart
    real(rp) :: ratio_start,ratio_end
    !
    is_due = .false.
    if(period > 0.) then
      if(present(restart)) then
        is_due = .not.restart
      else
        ratio_start = time_start/period
        ratio_end   = time_end/period
        is_due = floor(ratio_end  +tol(istep  ,ratio_end  ),kind=i8) > &
                 floor(ratio_start+tol(istep-1,ratio_start),kind=i8)
      end if
    else if(interval > 0) then
      is_due = mod(istep,interval) == 0
    end if
  contains
    pure real(rp) function tol(k,ratio)
      !
      ! round-off tolerance of the time after k steps, which accumulates errors of up to ~k*epsilon*time;
      ! capped at sqrt(epsilon) of a period, and evaluated with the step count of each end of a time step,
      ! so that consecutive steps agree
      !
      implicit none
      integer , intent(in) :: k
      real(rp), intent(in) :: ratio
      tol = min(4*epsilon(1._rp)*max(1,k)*max(1._rp,abs(ratio)),sqrt(epsilon(1._rp)))
    end function tol
  end function is_output_due
  !
  subroutine out0d(fname,n,var)
    !
    ! appends the first n entries of an array
    ! var to a file
    ! fname -> name of the file
    ! n     -> number of entries
    ! var   -> input array of real values
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in) :: n
    real(rp), intent(in), dimension(:) :: var
    integer :: iunit
    !
    if(myid == 0) then
      open(newunit=iunit,file=fname,position='append')
      write(iunit,fmt_rp) var(1:n)
      close(iunit)
    end if
  end subroutine out0d
  !
  subroutine gen_alias(myid,datadir,fname,fname_alias)
    !
    ! this subroutine generates a symlink with name `datadir//fname_alias`, pointing to
    ! file `fname` (relative to `datadir`), created with a temporary name and then renamed,
    ! which replaces an existing alias atomically; it is called by task `myid`
    !
    integer, intent(in) :: myid
    character(len=*), intent(in) :: datadir,fname,fname_alias
    character(len=:), allocatable :: flink,ftmp
    integer :: istat
    if(myid /= 0) return
    flink = trim(datadir)//trim(fname_alias)
    ftmp  = flink//'.tmp'
    istat = remove_file(ftmp) ! possible leftover of an interrupted call
    istat = link_file(fname,ftmp)
    if(istat == 0) istat = rename_file(ftmp,flink)
    if(istat /= 0) print*, 'WARNING: could not generate the alias ', flink, ' to ', trim(fname), '.'
  end subroutine gen_alias
  !
  subroutine out1d(fname,ng,nh,nh_g,lo,hi,idir,l,dl,z_g,dz,p)
    !
    ! writes the profile of a variable averaged
    ! over two domain directions
    !
    ! fname -> name of the file
    ! ng    -> global domain sizes
    ! lo,hi -> upper and lower extents of the input array
    ! idir  -> direction of the profile
    ! dl,l  -> uniform grid spacing and length arrays
    ! z_g   -> global z coordinate array (grid is non-uniform in z)
    ! dz    -> local z grid spacing array (should work also with the global one)
    ! p     -> 3D input scalar field
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in)               :: nh,nh_g
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(1-nh_g:) :: z_g
    real(rp), intent(in), dimension(lo(3)-nh_g:) :: dz
    real(rp), intent(in), dimension(lo(1)-nh:,lo(2)-nh:,lo(3)-nh:) :: p
    real(rp), allocatable, dimension(:) :: p1d
    integer :: i,j,k
    integer :: iunit
    real(rp) :: grid_area_ratio,p1d_s
    !
    allocate(p1d(ng(idir)))
    !$acc enter data create(p1d)
    !$acc parallel loop default(present)
    do k=1,size(p1d)
      p1d(k) = 0._rp
    end do
    select case(idir)
    case(3)
      grid_area_ratio = dl(1)*dl(2)/(l(1)*l(2))
      !$acc parallel loop gang default(present) private(p1d_s)
      do k=lo(3),hi(3)
        p1d_s = 0._rp
        !$acc loop collapse(2) reduction(+:p1d_s)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            p1d_s = p1d_s + p(i,j,k)*grid_area_ratio
          end do
        end do
        p1d(k) = p1d_s
      end do
      !$acc exit data copyout(p1d)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do k=1,ng(3)
          write(iunit,fmt_rp) z_g(k),p1d(k)
        end do
        close(iunit)
      end if
    case(2)
      grid_area_ratio = dl(1)/(l(1)*l(3))
      !$acc parallel loop gang default(present) private(p1d_s)
      do j=lo(2),hi(2)
        p1d_s = 0._rp
        !$acc loop collapse(2) reduction(+:p1d_s)
        do k=lo(3),hi(3)
          do i=lo(1),hi(1)
            p1d_s = p1d_s + p(i,j,k)*dz(k)*grid_area_ratio
          end do
        end do
        p1d(j) = p1d_s
      end do
      !$acc exit data copyout(p1d)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1),ng(2),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do j=1,ng(2)
          write(iunit,fmt_rp) (j-.5)*dl(2),p1d(j)
        end do
        close(iunit)
      end if
    case(1)
      grid_area_ratio = dl(2)/(l(2)*l(3))
      !$acc parallel loop gang default(present) private(p1d_s)
      do i=lo(1),hi(1)
        p1d_s = 0._rp
        !$acc loop collapse(2) reduction(+:p1d_s)
        do k=lo(3),hi(3)
          do j=lo(2),hi(2)
            p1d_s = p1d_s + p(i,j,k)*dz(k)*grid_area_ratio
          end do
        end do
        p1d(i) = p1d_s
      end do
      !$acc exit data copyout(p1d)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1),ng(1),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do i=1,ng(1)
          write(iunit,fmt_rp) (i-.5)*dl(1),p1d(i)
        end do
        close(iunit)
      end if
    end select
  end subroutine out1d
  !
  subroutine out2d(fname,ng,nh,nh_g,lo,hi,inorm,islice,p,varname,time,istep,x_g,y_g,z_g,is_pack)
    !
    ! saves a planar slice of a scalar field using the structured subset engine
    !
    ! fname  -> name of the output file
    ! ng     -> global domain sizes
    ! nh     -> local halo extents
    ! nh_g   -> halo width of the global coordinate arrays
    ! lo,hi  -> local interior ownership in global indices
    ! inorm  -> plane is perpendicular to direction inorm (1,2,3)
    ! islice -> plane is of constant index islice in direction inorm
    ! p      -> 3D input scalar field
    ! varname,time,istep,x_g,y_g,z_g -> subset metadata
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,nh,lo,hi
    integer , intent(in)               :: nh_g
    integer , intent(in) :: inorm,islice
    real(rp), intent(in), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: p
    character(len=*), intent(in) :: varname
    real(rp), intent(in) :: time
    integer , intent(in) :: istep
    real(rp), intent(in), dimension(1-nh_g:) :: x_g
    real(rp), intent(in), dimension(1-nh_g:) :: y_g
    real(rp), intent(in), dimension(1-nh_g:) :: z_g
    logical , intent(in), optional :: is_pack
    integer, dimension(3) :: lo_out,hi_out
    !
    select case(inorm)
    case(1)
      lo_out(:) = [islice,1     ,1   ]
      hi_out(:) = [islice,ng(2),ng(3)]
    case(2)
      lo_out(:) = [1    ,islice,1   ]
      hi_out(:) = [ng(1),islice,ng(3)]
    case(3)
      lo_out(:) = [1    ,1    ,islice]
      hi_out(:) = [ng(1),ng(2),islice]
    end select
    call io_write_subset(fname,varname,MPI_COMM_WORLD,ng,nh,nh_g,lo,hi,lo_out,hi_out,[1,1,1],p,time,istep,x_g,y_g,z_g,is_pack)
  end subroutine out2d
  !
  subroutine out3d(fname,ng,nh,nh_g,lo,hi,lo_out,hi_out,nskip,p,varname,time,istep,x_g,y_g,z_g,is_pack)
    !
    ! saves a structured subset of a scalar field using the subset engine
    !
    ! fname  -> name of the output file
    ! ng     -> global domain sizes
    ! nh     -> local halo extents
    ! nh_g   -> halo width of the global coordinate arrays
    ! lo,hi  -> local interior ownership in global indices
    ! lo_out -> first global point written in each direction
    ! hi_out -> last  global point written in each direction
    ! nskip  -> output stride in each direction
    ! p      -> 3D input scalar field
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,nh,lo,hi,lo_out,hi_out
    integer , intent(in)               :: nh_g
    integer , intent(in), dimension(3) :: nskip
    real(rp), intent(in), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: p
    character(len=*), intent(in) :: varname
    real(rp), intent(in) :: time
    integer , intent(in) :: istep
    real(rp), intent(in), dimension(1-nh_g:) :: x_g
    real(rp), intent(in), dimension(1-nh_g:) :: y_g
    real(rp), intent(in), dimension(1-nh_g:) :: z_g
    logical , intent(in), optional :: is_pack
    !
    call io_write_subset(fname,varname,MPI_COMM_WORLD,ng,nh,nh_g,lo,hi,lo_out,hi_out,nskip,p,time,istep,x_g,y_g,z_g,is_pack)
  end subroutine out3d
  !
  subroutine write_log_output(fname,fname_fld,varname,lo_out,hi_out,nskip,time,istep)
    !
    ! appends information about a saved binary file to a log file
    ! this file is used to generate a xdmf file for visualization of field data
    !
    ! fname     -> name of the output log file
    ! fname_fld -> name of the saved binary file (excluding the directory)
    ! varname   -> name of the variable that is saved
    ! lo_out    -> first element of the field that is saved in each direction, e.g. [1,1,1]
    ! hi_out    -> last  element of the field that is saved in each direction, e.g. [ng(1),ng(2),ng(3)]
    ! nskip     -> step size between nmin and nmax, e.g. [1,1,1] if the whole array is saved
    ! time      -> physical time
    ! istep     -> time step number
    !
    implicit none
    character(len=*), intent(in) :: fname,fname_fld,varname
    integer , intent(in), dimension(3) :: lo_out,hi_out,nskip
    real(rp), intent(in)               :: time
    integer , intent(in)               :: istep
    character(len=100) :: cfmt
    integer :: iunit
    !
    write(cfmt, '(A)') '(A,A,A,9i5,E16.7e3,i7)'
    if(myid  ==  0) then
      open(newunit=iunit,file=fname,position='append')
      write(iunit,trim(cfmt)) trim(fname_fld),' ',trim(varname),lo_out,hi_out,nskip,time,istep
      close(iunit)
    end if
  end subroutine write_log_output
  !
  subroutine write_visu_subset(datadir,fname_bin,fname_log,varname,ng,nh,nh_g,lo,hi,lo_out,hi_out,nskip,time,istep,p, &
                               x_g,y_g,z_g,is_pack)
    !
    ! wraps subset data output and subset metadata logging
    !
    implicit none
    character(len=*), intent(in)          :: datadir,fname_bin,fname_log,varname
    integer , intent(in), dimension(3)    :: ng,nh,lo,hi,lo_out,hi_out,nskip
    integer , intent(in)                  :: nh_g
    real(rp), intent(in)                  :: time
    integer , intent(in)                  :: istep
    real(rp), intent(in), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: p
    real(rp), intent(in), dimension(1-nh_g:)                                :: x_g
    real(rp), intent(in), dimension(1-nh_g:)                                :: y_g
    real(rp), intent(in), dimension(1-nh_g:)                                :: z_g
    logical , intent(in), optional        :: is_pack
    !
    call out3d(trim(datadir)//trim(fname_bin),ng,nh,nh_g,lo,hi,lo_out,hi_out,nskip,p,varname,time,istep,x_g,y_g,z_g,is_pack)
    call write_log_output(trim(datadir)//trim(fname_log),trim(fname_bin),trim(varname),lo_out,hi_out,nskip,time,istep)
  end subroutine write_visu_subset
  !
  subroutine write_visu_3d(datadir,fname_bin,fname_log,varname,ng,nh,nh_g,lo,hi,lo_out,hi_out,nskip,time,istep,p, &
                           x_g,y_g,z_g,is_pack)
    !
    ! wraps the calls of out3d and write_log_output into the same subroutine
    !
    implicit none
    character(len=*), intent(in)          :: datadir,fname_bin,fname_log,varname
    integer , intent(in), dimension(3)    :: ng,nh,lo,hi,lo_out,hi_out,nskip
    integer , intent(in)                  :: nh_g
    real(rp), intent(in)                  :: time
    integer , intent(in)                  :: istep
    real(rp), intent(in), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: p
    real(rp), intent(in), dimension(1-nh_g:)                                :: x_g
    real(rp), intent(in), dimension(1-nh_g:)                                :: y_g
    real(rp), intent(in), dimension(1-nh_g:)                                :: z_g
    logical , intent(in), optional        :: is_pack
    !
    call write_visu_subset(datadir,fname_bin,fname_log,varname,ng,nh,nh_g,lo,hi,lo_out,hi_out,nskip,time,istep,p, &
                           x_g,y_g,z_g,is_pack)
  end subroutine write_visu_3d
  !
  subroutine write_visu_2d(datadir,fname_bin,fname_log,varname,ng,nh,nh_g,lo,hi,inorm,nslice,time,istep,p,x_g,y_g,z_g,is_pack)
    !
    ! wraps the calls of out2d and write-log_output into the same subroutine
    !
    implicit none
    character(len=*), intent(in)          :: datadir,fname_bin,fname_log,varname
    integer , intent(in), dimension(3)    :: ng,nh,lo,hi
    integer , intent(in)                  :: nh_g
    integer , intent(in)                  :: inorm,nslice
    real(rp), intent(in)                  :: time
    integer , intent(in)                  :: istep
    real(rp), intent(in), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: p
    real(rp), intent(in), dimension(1-nh_g:)                                :: x_g
    real(rp), intent(in), dimension(1-nh_g:)                                :: y_g
    real(rp), intent(in), dimension(1-nh_g:)                                :: z_g
    logical , intent(in), optional        :: is_pack
    integer , dimension(3) :: lo_out,hi_out
    !
    select case(inorm)
    case(1)
      lo_out(:) = [nslice,1     ,1   ]
      hi_out(:) = [nslice,ng(2),ng(3)]
    case(2)
      lo_out(:) = [1    ,nslice,1    ]
      hi_out(:) = [ng(1),nslice,ng(3)]
    case(3)
      lo_out(:) = [1    ,1    ,nslice]
      hi_out(:) = [ng(1),ng(2),nslice]
    end select
    call write_visu_subset(datadir,fname_bin,fname_log,varname,ng,nh,nh_g,lo,hi,lo_out,hi_out,[1,1,1],time,istep,p, &
                           x_g,y_g,z_g,is_pack)
  end subroutine write_visu_2d
  !
  subroutine out1d_chan(fname,ng,nh,nh_g,lo,hi,idir,l,dl,z_g,u,v,w) ! e.g. for a channel with streamwise dir in x
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in)               :: nh,nh_g
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(1-nh_g:) :: z_g
    real(rp), intent(in), dimension(lo(1)-nh:,lo(2)-nh:,lo(3)-nh:) :: u,v,w
    real(rp), allocatable, dimension(:) :: um,vm,wm,u2,v2,w2,uw
    integer :: i,j,k
    integer :: iunit
    integer :: q
    real(rp) :: grid_area_ratio
    !
    q = ng(idir)
    select case(idir)
    case(3)
      grid_area_ratio = dl(1)*dl(2)/(l(1)*l(2))
      allocate(um(0:q+1),vm(0:q+1),wm(0:q+1),u2(0:q+1),v2(0:q+1),w2(0:q+1),uw(0:q+1))
      um(:) = 0.
      vm(:) = 0.
      wm(:) = 0.
      u2(:) = 0.
      v2(:) = 0.
      w2(:) = 0.
      uw(:) = 0.
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            um(k) = um(k) + u(i,j,k)
            vm(k) = vm(k) + v(i,j,k)
            wm(k) = wm(k) + 0.50*(w(i,j,k-1) + w(i,j,k))
            u2(k) = u2(k) + u(i,j,k)**2
            v2(k) = v2(k) + v(i,j,k)**2
            w2(k) = w2(k) + 0.50*(w(i,j,k)**2+w(i,j,k-1)**2)
            uw(k) = uw(k) + 0.25*(u(i-1,j,k) + u(i,j,k))* &
                                 (w(i,j,k-1) + w(i,j,k))
          end do
        end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,um(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vm(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,wm(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,u2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,v2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,w2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,uw(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      um(:) = um(:)*grid_area_ratio
      vm(:) = vm(:)*grid_area_ratio
      wm(:) = wm(:)*grid_area_ratio
      u2(:) = u2(:)*grid_area_ratio - um(:)**2
      v2(:) = v2(:)*grid_area_ratio - vm(:)**2
      w2(:) = w2(:)*grid_area_ratio - wm(:)**2
      uw(:) = uw(:)*grid_area_ratio - um(:)*wm(:)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do k=1,ng(3)
          write(iunit,fmt_rp) z_g(k),um(k),vm(k),wm(k), &
                                     u2(k),v2(k),w2(k), &
                                     uw(k)
        end do
        close(iunit)
      end if
    case(2)
    case(1)
    end select
  end subroutine out1d_chan
  !
  subroutine out2d_duct(fname,ng,nh,nh_g,lo,hi,idir,l,dl,z_g,u,v,w) ! e.g. for a duct
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in)               :: nh,nh_g
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(1-nh_g:) :: z_g
    real(rp), intent(in), dimension(lo(1)-nh:,lo(2)-nh:,lo(3)-nh:) :: u,v,w
    real(rp), allocatable, dimension(:,:) :: um,vm,wm,u2,v2,w2,uv,uw,vw
    integer :: i,j,k
    integer :: iunit
    integer :: p,q
    real(rp) :: x_g,y_g,grid_area_ratio
    !
    select case(idir) ! streamwise direction
    case(3)
    case(2)
      grid_area_ratio = dl(2)/l(2)
      p = ng(1)
      q = ng(3)
      allocate(um(p,q),vm(p,q),wm(p,q),u2(p,q),v2(p,q),w2(p,q),uv(p,q),vw(p,q))
      !
      um(:,:) = 0.
      vm(:,:) = 0.
      wm(:,:) = 0.
      u2(:,:) = 0.
      v2(:,:) = 0.
      w2(:,:) = 0.
      uv(:,:) = 0.
      vw(:,:) = 0.
      do k=lo(3),hi(3)
        do i=lo(1),hi(1)
          do j=lo(2),hi(2)
            um(i,k) = um(i,k) + 0.5*(u(i-1,j,k)+u(i,j,k))
            vm(i,k) = vm(i,k) + v(i,j,k)
            wm(i,k) = wm(i,k) + 0.5*(w(i,j,k-1)+w(i,j,k))
            u2(i,k) = u2(i,k) + 0.5*(u(i-1,j,k)**2+u(i,j,k)**2)
            v2(i,k) = v2(i,k) + v(i,j,k)**2
            w2(i,k) = w2(i,k) + 0.5*(w(i,j,k-1)**2+w(i,j,k)**2)
            vw(i,k) = vw(i,k) + 0.25*(v(i,j-1,k) + v(i,j,k))* &
                                     (w(i,j,k-1) + w(i,j,k))
            uv(i,k) = uv(i,k) + 0.25*(u(i-1,j,k) + u(i,j,k))* &
                                     (v(i,j-1,k) + v(i,j,k))
          end do
        end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,um(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vm(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,wm(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,u2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,v2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,w2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vw(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,uv(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      um(:,:) = um(:,:)*grid_area_ratio
      vm(:,:) = vm(:,:)*grid_area_ratio
      wm(:,:) = wm(:,:)*grid_area_ratio
      u2(:,:) = u2(:,:)*grid_area_ratio - um(:,:)**2
      v2(:,:) = v2(:,:)*grid_area_ratio - vm(:,:)**2
      w2(:,:) = w2(:,:)*grid_area_ratio - wm(:,:)**2
      vw(:,:) = vw(:,:)*grid_area_ratio - vm(:,:)*wm(:,:)
      uv(:,:) = uv(:,:)*grid_area_ratio - um(:,:)*vm(:,:)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do k=1,ng(3)
          do i=1,ng(1)
            x_g = (i-.5)*dl(1)
            write(iunit,fmt_rp) x_g,z_g(k),um(i,k),vm(i,k),wm(i,k), &
                                           u2(i,k),v2(i,k),w2(i,k), &
                                           vw(i,k),uv(i,k)
          end do
        end do
        close(iunit)
      end if
    case(1)
      grid_area_ratio = dl(1)/l(1)
      p = ng(2)
      q = ng(3)
      allocate(um(p,q),vm(p,q),wm(p,q),u2(p,q),v2(p,q),w2(p,q),uv(p,q),uw(p,q))
      !
      um(:,:) = 0.
      vm(:,:) = 0.
      wm(:,:) = 0.
      u2(:,:) = 0.
      v2(:,:) = 0.
      w2(:,:) = 0.
      uv(:,:) = 0.
      uw(:,:) = 0.
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            um(j,k) = um(j,k) + u(i,j,k)
            vm(j,k) = vm(j,k) + 0.5*(v(i,j-1,k)+v(i,j,k))
            wm(j,k) = wm(j,k) + 0.5*(w(i,j,k-1)+w(i,j,k))
            u2(j,k) = u2(j,k) + u(i,j,k)**2
            v2(j,k) = v2(j,k) + 0.5*(v(i,j-1,k)**2+v(i,j,k)**2)
            w2(j,k) = w2(j,k) + 0.5*(w(i,j,k-1)**2+w(i,j,k)**2)
            uv(j,k) = uv(j,k) + 0.25*(u(i-1,j,k) + u(i,j,k))* &
                                     (v(i,j-1,k) + v(i,j,k))
            uw(j,k) = uw(j,k) + 0.25*(u(i-1,j,k) + u(i,j,k))* &
                                     (w(i,j,k-1) + w(i,j,k))
          end do
        end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,um(1,1),ng(2)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vm(1,1),ng(2)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,wm(1,1),ng(2)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,u2(1,1),ng(2)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,v2(1,1),ng(2)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,w2(1,1),ng(2)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,uv(1,1),ng(2)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,uw(1,1),ng(2)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      um(:,:) = um(:,:)*grid_area_ratio
      vm(:,:) = vm(:,:)*grid_area_ratio
      wm(:,:) = wm(:,:)*grid_area_ratio
      u2(:,:) = u2(:,:)*grid_area_ratio - um(:,:)**2
      v2(:,:) = v2(:,:)*grid_area_ratio - vm(:,:)**2
      w2(:,:) = w2(:,:)*grid_area_ratio - wm(:,:)**2
      uv(:,:) = uv(:,:)*grid_area_ratio - um(:,:)*vm(:,:)
      uw(:,:) = uw(:,:)*grid_area_ratio - um(:,:)*wm(:,:)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do k=1,ng(3)
          do j=1,ng(2)
            y_g = (j-.5)*dl(2)
            write(iunit,fmt_rp) y_g,z_g(k),um(j,k),vm(j,k),wm(j,k), &
                                           u2(j,k),v2(j,k),w2(j,k), &
                                           uv(j,k),uw(j,k)
          end do
        end do
        close(iunit)
      end if
    end select
  end subroutine out2d_duct
  !
  ! two-fluid diagnostics
  !
  subroutine velstats(fname,n,nh_v,nh_p,nh_s,nh_psi,nh_g,dl,l,zc,u,v,w,p,s,psi)
    !
    ! computes a few volume-averaged profiles
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer, intent(in), dimension(3) :: n
    integer , intent(in) :: nh_v,nh_p,nh_s,nh_psi,nh_g
    real(rp), intent(in), dimension(3) :: dl,l
    real(rp), intent(in), dimension(1-nh_g:) :: zc
    real(rp), intent(in), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(in), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    real(rp), intent(in), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s
    real(rp), intent(in), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), allocatable, dimension(:) :: u1_1,u1_2,v1_1,v1_2,w1_1,w1_2,p1_1,p1_2,s1_1,s1_2,c1_1,c1_2, &
                                           u2_1,u2_2,v2_1,v2_2,w2_1,w2_2,p2_1,p2_2,s2_1,s2_2,c2_1,c2_2
    integer :: i,j,k,ii,jj
    integer :: iunit
    integer, dimension(3) :: ng
    integer :: nn
    real(rp) :: grid_area_ratio
    !
    nn = ng(3)
    allocate(u1_1(nn),v1_1(nn),w1_1(nn),p1_1(nn),s1_1(nn),c1_1(nn), &
             u1_2(nn),v1_2(nn),w1_2(nn),p1_2(nn),s1_2(nn),c1_2(nn), &
             u2_1(nn),v2_1(nn),w2_1(nn),p2_1(nn),s2_1(nn),c2_1(nn), &
             u2_2(nn),v2_2(nn),w2_2(nn),p2_2(nn),s2_2(nn),c2_2(nn))
    grid_area_ratio = dl(1)*dl(2)/l(1)*l(2)
    do k=1,n(3)
      u1_1(k)  = 0.
      v1_1(k)  = 0.
      w1_1(k)  = 0.
      p1_1(k)  = 0.
      s1_1(k)  = 0.
      c1_1(k)  = 0.
      u1_2(k)  = 0.
      v1_2(k)  = 0.
      w1_2(k)  = 0.
      p1_2(k)  = 0.
      s1_2(k)  = 0.
      c1_2(k)  = 0.
      u2_1(k)  = 0.
      v2_1(k)  = 0.
      w2_1(k)  = 0.
      p2_1(k)  = 0.
      s2_1(k)  = 0.
      c2_1(k)  = 0.
      u2_2(k)  = 0.
      v2_2(k)  = 0.
      w2_2(k)  = 0.
      p2_2(k)  = 0.
      s2_2(k)  = 0.
      c2_2(k)  = 0.
      do j=1,n(2)
        do i=1,n(1)
          u1_1(k)  = u1_1(k) + grid_area_ratio*(   psi(i,j,k))*(0.5*(u(i,j,k)+u(i-1,j,k)))
          u1_2(k)  = u1_2(k) + grid_area_ratio*(   psi(i,j,k))*(0.5*(u(i,j,k)+u(i-1,j,k)))**2
          v1_1(k)  = v1_1(k) + grid_area_ratio*(   psi(i,j,k))*(0.5*(v(i,j,k)+v(i,j-1,k)))
          v1_2(k)  = v1_2(k) + grid_area_ratio*(   psi(i,j,k))*(0.5*(v(i,j,k)+v(i,j-1,k)))**2
          w1_1(k)  = w1_1(k) + grid_area_ratio*(   psi(i,j,k))*(0.5*(w(i,j,k)+w(i,j,k-1)))
          w1_2(k)  = w1_2(k) + grid_area_ratio*(   psi(i,j,k))*(0.5*(w(i,j,k)+w(i,j,k-1)))**2
          p1_1(k)  = p1_1(k) + grid_area_ratio*(   psi(i,j,k))*p(i,j,k)
          p1_2(k)  = p1_2(k) + grid_area_ratio*(   psi(i,j,k))*p(i,j,k)**2
          s1_1(k)  = s1_1(k) + grid_area_ratio*(   psi(i,j,k))*s(i,j,k)
          s1_2(k)  = s1_2(k) + grid_area_ratio*(   psi(i,j,k))*s(i,j,k)**2
          c1_1(k)  = c1_1(k) + grid_area_ratio*(   psi(i,j,k))
          c1_2(k)  = c1_2(k) + grid_area_ratio*(   psi(i,j,k))**2
          u2_1(k)  = u2_1(k) + grid_area_ratio*(1.-psi(i,j,k))*(0.5*(u(i,j,k)+u(i-1,j,k)))
          u2_2(k)  = u2_2(k) + grid_area_ratio*(1.-psi(i,j,k))*(0.5*(u(i,j,k)+u(i-1,j,k)))**2
          v2_1(k)  = v2_1(k) + grid_area_ratio*(1.-psi(i,j,k))*(0.5*(v(i,j,k)+v(i,j-1,k)))
          v2_2(k)  = v2_2(k) + grid_area_ratio*(1.-psi(i,j,k))*(0.5*(v(i,j,k)+v(i,j-1,k)))**2
          w2_1(k)  = w2_1(k) + grid_area_ratio*(1.-psi(i,j,k))*(0.5*(w(i,j,k)+w(i,j,k-1)))
          w2_2(k)  = w2_2(k) + grid_area_ratio*(1.-psi(i,j,k))*(0.5*(w(i,j,k)+w(i,j,k-1)))**2
          p2_1(k)  = p2_1(k) + grid_area_ratio*(1.-psi(i,j,k))*p(i,j,k)
          p2_2(k)  = p2_2(k) + grid_area_ratio*(1.-psi(i,j,k))*p(i,j,k)**2
          s2_1(k)  = s2_1(k) + grid_area_ratio*(1.-psi(i,j,k))*s(i,j,k)
          s2_2(k)  = s2_2(k) + grid_area_ratio*(1.-psi(i,j,k))*s(i,j,k)**2
          c2_1(k)  = c2_1(k) + grid_area_ratio*(1.-psi(i,j,k))
          c2_2(k)  = c2_2(k) + grid_area_ratio*(1.-psi(i,j,k))**2
        end do
      end do
    end do
    call mpi_allreduce(MPI_IN_PLACE,u1_1,nn,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,v1_1,nn,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,w1_1,nn,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,p1_1,nn,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,c1_1,nn,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,u1_2,nn,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,v1_2,nn,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,w1_2,nn,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,p1_2,nn,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,c1_2,nn,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,u2_1,nn,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,v2_1,nn,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,w2_1,nn,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,p2_1,nn,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,c2_1,nn,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,u2_2,nn,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,v2_2,nn,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,w2_2,nn,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,p2_2,nn,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,c2_2,nn,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    if(myid == 0) then
      open(newunit=iunit,file=fname)
      do k=1,nn
        write(iunit,fmt_rp) zc(k), &
                            u1_1(k),v1_1(k),w1_1(k),p1_1(k),s1_1(k),c1_1(k), &
                            u1_2(k),v1_2(k),w1_2(k),p1_2(k),s1_2(k),c1_2(k), &
                            u2_1(k),v2_1(k),w2_1(k),p2_1(k),s2_1(k),c2_1(k), &
                            u2_2(k),v2_2(k),w2_2(k),p2_2(k),s2_2(k),c2_2(k)
      end do
      close(iunit)
    end if
  end subroutine velstats
  !
  subroutine cmpt_total_area(n,nh_psi,nh_g,dli,dzci,dzfi,psi,area)
    !
    ! computes some volume-averaged profiles along z (assumes constant grid spacing)
    !
    implicit none
    integer , intent(in ), dimension(3) :: n
    integer , intent(in) :: nh_psi,nh_g
    real(rp), intent(in ), dimension(3) :: dli
    real(rp), intent(in ), dimension(1-nh_g:) :: dzci,dzfi
    real(rp), intent(in ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(out) :: area
    real(rp), dimension(8) :: mx,my,mz
    real(rp) :: dpsidx,dpsidy,dpsidz
    integer  :: i,j,k
    area = 0._rp
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          mx(1) = 0.25*((psi(i+1,j  ,k  )+psi(i+1,j+1,k  )+psi(i+1,j  ,k+1)+psi(i+1,j+1,k+1)) - &
                        (psi(i  ,j  ,k  )+psi(i  ,j+1,k  )+psi(i  ,j  ,k+1)+psi(i  ,j+1,k+1)))*dli(1)
          mx(2) = 0.25*((psi(i+1,j  ,k  )+psi(i+1,j-1,k  )+psi(i+1,j  ,k+1)+psi(i+1,j-1,k+1)) - &
                        (psi(i  ,j  ,k  )+psi(i  ,j-1,k  )+psi(i  ,j  ,k+1)+psi(i  ,j-1,k+1)))*dli(1)
          mx(3) = 0.25*((psi(i+1,j  ,k  )+psi(i+1,j+1,k  )+psi(i+1,j  ,k-1)+psi(i+1,j+1,k-1)) - &
                        (psi(i  ,j  ,k  )+psi(i  ,j+1,k  )+psi(i  ,j  ,k-1)+psi(i  ,j+1,k-1)))*dli(1)
          mx(4) = 0.25*((psi(i+1,j  ,k  )+psi(i+1,j-1,k  )+psi(i+1,j  ,k-1)+psi(i+1,j-1,k-1)) - &
                        (psi(i  ,j  ,k  )+psi(i  ,j-1,k  )+psi(i  ,j  ,k-1)+psi(i  ,j-1,k-1)))*dli(1)
          mx(5) = 0.25*((psi(i  ,j  ,k  )+psi(i  ,j+1,k  )+psi(i  ,j  ,k+1)+psi(i  ,j+1,k+1)) - &
                        (psi(i-1,j  ,k  )+psi(i-1,j+1,k  )+psi(i-1,j  ,k+1)+psi(i-1,j+1,k+1)))*dli(1)
          mx(6) = 0.25*((psi(i  ,j  ,k  )+psi(i  ,j-1,k  )+psi(i  ,j  ,k+1)+psi(i  ,j-1,k+1)) - &
                        (psi(i-1,j  ,k  )+psi(i-1,j-1,k  )+psi(i-1,j  ,k+1)+psi(i-1,j-1,k+1)))*dli(1)
          mx(7) = 0.25*((psi(i  ,j  ,k  )+psi(i  ,j+1,k  )+psi(i  ,j  ,k-1)+psi(i  ,j+1,k-1)) - &
                        (psi(i-1,j  ,k  )+psi(i-1,j+1,k  )+psi(i-1,j  ,k-1)+psi(i-1,j+1,k-1)))*dli(1)
          mx(8) = 0.25*((psi(i  ,j  ,k  )+psi(i  ,j-1,k  )+psi(i  ,j  ,k-1)+psi(i  ,j-1,k-1)) - &
                        (psi(i-1,j  ,k  )+psi(i-1,j-1,k  )+psi(i-1,j  ,k-1)+psi(i-1,j-1,k-1)))*dli(1)
          !
          my(1) = 0.25*((psi(i  ,j+1,k  )+psi(i+1,j+1,k  )+psi(i  ,j+1,k+1)+psi(i+1,j+1,k+1)) - &
                        (psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j  ,k+1)+psi(i+1,j  ,k+1)))*dli(2)
          my(2) = 0.25*((psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j  ,k+1)+psi(i+1,j  ,k+1)) - &
                        (psi(i  ,j-1,k  )+psi(i+1,j-1,k  )+psi(i  ,j-1,k+1)+psi(i+1,j-1,k+1)))*dli(2)
          my(3) = 0.25*((psi(i  ,j+1,k  )+psi(i+1,j+1,k  )+psi(i  ,j+1,k-1)+psi(i+1,j+1,k-1)) - &
                        (psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j  ,k-1)+psi(i+1,j  ,k-1)))*dli(2)
          my(4) = 0.25*((psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j  ,k-1)+psi(i+1,j  ,k-1)) - &
                        (psi(i  ,j-1,k  )+psi(i+1,j-1,k  )+psi(i  ,j-1,k-1)+psi(i+1,j-1,k-1)))*dli(2)
          my(5) = 0.25*((psi(i  ,j+1,k  )+psi(i-1,j+1,k  )+psi(i  ,j+1,k+1)+psi(i-1,j+1,k+1)) - &
                        (psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j  ,k+1)+psi(i-1,j  ,k+1)))*dli(2)
          my(6) = 0.25*((psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j  ,k+1)+psi(i-1,j  ,k+1)) - &
                        (psi(i  ,j-1,k  )+psi(i-1,j-1,k  )+psi(i  ,j-1,k+1)+psi(i-1,j-1,k+1)))*dli(2)
          my(7) = 0.25*((psi(i  ,j+1,k  )+psi(i-1,j+1,k  )+psi(i  ,j+1,k-1)+psi(i-1,j+1,k-1)) - &
                        (psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j  ,k-1)+psi(i-1,j  ,k-1)))*dli(2)
          my(8) = 0.25*((psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j  ,k-1)+psi(i-1,j  ,k-1)) - &
                        (psi(i  ,j-1,k  )+psi(i-1,j-1,k  )+psi(i  ,j-1,k-1)+psi(i-1,j-1,k-1)))*dli(2)
          !
          mz(1) = 0.25*((psi(i  ,j  ,k+1)+psi(i+1,j  ,k+1)+psi(i  ,j+1,k+1)+psi(i+1,j+1,k+1)) - &
                        (psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j+1,k  )+psi(i+1,j+1,k  )))*dzci(k  )
          mz(2) = 0.25*((psi(i  ,j  ,k+1)+psi(i+1,j  ,k+1)+psi(i  ,j-1,k+1)+psi(i+1,j-1,k+1)) - &
                        (psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j-1,k  )+psi(i+1,j-1,k  )))*dzci(k  )
          mz(3) = 0.25*((psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j+1,k  )+psi(i+1,j+1,k  )) - &
                        (psi(i  ,j  ,k-1)+psi(i+1,j  ,k-1)+psi(i  ,j+1,k-1)+psi(i+1,j+1,k-1)))*dzci(k  )
          mz(4) = 0.25*((psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j-1,k  )+psi(i+1,j-1,k  )) - &
                        (psi(i  ,j  ,k-1)+psi(i+1,j  ,k-1)+psi(i  ,j-1,k-1)+psi(i+1,j-1,k-1)))*dzci(k  )
          mz(5) = 0.25*((psi(i  ,j  ,k+1)+psi(i-1,j  ,k+1)+psi(i  ,j+1,k+1)+psi(i-1,j+1,k+1)) - &
                        (psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j+1,k  )+psi(i-1,j+1,k  )))*dzci(k+1)
          mz(6) = 0.25*((psi(i  ,j  ,k+1)+psi(i-1,j  ,k+1)+psi(i  ,j-1,k+1)+psi(i-1,j-1,k+1)) - &
                        (psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j-1,k  )+psi(i-1,j-1,k  )))*dzci(k+1)
          mz(7) = 0.25*((psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j+1,k  )+psi(i-1,j+1,k  )) - &
                        (psi(i  ,j  ,k-1)+psi(i-1,j  ,k-1)+psi(i  ,j+1,k-1)+psi(i-1,j+1,k-1)))*dzci(k  )
          mz(8) = 0.25*((psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j-1,k  )+psi(i-1,j-1,k  )) - &
                        (psi(i  ,j  ,k-1)+psi(i-1,j  ,k-1)+psi(i  ,j-1,k-1)+psi(i-1,j-1,k-1)))*dzci(k  )
          !
          dpsidx = .125*(mx(1)+mx(2)+mx(3)+mx(4)+mx(5)+mx(6)+mx(7)+mx(8))
          dpsidy = .125*(my(1)+my(2)+my(3)+my(4)+my(5)+my(6)+my(7)+my(8))
          dpsidz = .125*(mz(1)+mz(2)+mz(3)+mz(4)+mz(5)+mz(6)+mz(7)+mz(8))
          area = area + sqrt(dpsidx**2+dpsidy**2+dpsidz**2)/(dli(1)*dli(2)*dzfi(k))
        enddo
      enddo
    enddo
    call MPI_ALLREDUCE(MPI_IN_PLACE,area,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  end subroutine cmpt_total_area
  !
  subroutine cmpt_mean_pos_and_vel(iphase,lo,hi,nh_v,nh_psi,nh_g,dl,dzf_g,zc_g,psi,u,v,w,pos,vel)
    !
    ! computes the mean position and velocity of one of the phases (e.g., useful for single-droplet or bubble benchmarks)
    !
    ! not too difficult to port to GPU, if later needed
    !
    implicit none
    integer , intent(in )               :: iphase ! phase 1 or 2?
    integer , intent(in ), dimension(3) :: lo,hi
    integer , intent(in) :: nh_v,nh_psi,nh_g
    real(rp), intent(in ), dimension(3) :: dl
    real(rp), intent(in ), dimension(lo(3)-nh_g:) :: dzf_g,zc_g
    real(rp), intent(in ), dimension(lo(1)-nh_psi:,lo(2)-nh_psi:,lo(3)-nh_psi:) :: psi
    real(rp), intent(in ), dimension(lo(1)-nh_v:,lo(2)-nh_v:,lo(3)-nh_v:) :: u,v,w
    real(rp), intent(out), dimension(3) :: pos,vel
    real(rp) :: aux,sgn,xc,yc,zc,vol_loc,vol
    integer  :: i,j,k
    !
    aux = 0.
    sgn = 1.
    if(      iphase == 1 ) then
      aux = aux
      sgn = sgn
    else if( iphase == 2 ) then
      aux =  1.
      sgn = -1.
    end if
    !
    pos(:) = 0._rp
    vel(:) = 0._rp
    vol    = 0.
    do k=lo(3),hi(3)
      zc = zc_g(k)
      do j=lo(2),hi(2)
        yc = (j-0.5)*dl(2)
        do i=lo(1),hi(1)
          xc = (i-0.5)*dl(1)
          !
          vol_loc = (aux+sgn*psi(i,j,k))*dl(1)*dl(2)*dzf_g(k)
          vol     = vol + vol_loc
          !
          pos(1) = pos(1) + xc*vol_loc
          pos(2) = pos(2) + yc*vol_loc
          pos(3) = pos(3) + zc*vol_loc
          vel(1) = vel(1) + 0.5*(u(i,j,k)+u(i-1,j,k))*vol_loc
          vel(2) = vel(2) + 0.5*(v(i,j,k)+v(i,j-1,k))*vol_loc
          vel(3) = vel(3) + 0.5*(w(i,j,k)+w(i,j,k-1))*vol_loc
        end do
      end do
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE,vol,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    pos(:) = pos(:)/vol
    vel(:) = vel(:)/vol
    call MPI_ALLREDUCE(MPI_IN_PLACE,pos,3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,vel,3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  end subroutine cmpt_mean_pos_and_vel
  !
  subroutine cmpt_mean_mass(n,nh_psi,nh_g,l,dl,dzf,rho12,psi,mass12)
    !
    ! computes mean mass of the system
    !
    implicit none
    integer , intent(in ), dimension(3)        :: n
    integer , intent(in) :: nh_psi,nh_g
    real(rp), intent(in ), dimension(3)        :: l,dl
    real(rp), intent(in ), dimension(1-nh_g:)       :: dzf
    real(rp), intent(in ), dimension(2)        :: rho12
    real(rp), intent(in ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(out), dimension(2)        :: mass12
    real(rp) :: vcell,mass1,mass2
    integer :: i,j,k
    !
    mass1 = 0._rp
    mass2 = 0._rp
    !$acc data copy(mass1,mass2) async(1)
    !$acc parallel loop collapse(3) default(present) reduction(+:mass1,mass2) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          vcell = dl(1)*dl(2)*dzf(k)
          mass1 = mass1 + rho12(1)*psi(i,j,k)*vcell
          mass2 = mass2 + rho12(2)*(1.-psi(i,j,k))*vcell
        end do
      end do
    end do
    !$acc end data
    !$acc wait(1)
    mass12(:) = [mass1,mass2]/product(l(:))
    call MPI_ALLREDUCE(MPI_IN_PLACE,mass12,2,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  end subroutine cmpt_mean_mass
  !
  subroutine cmpt_mean_energy(n,nh_v,nh_psi,nh_g,l,dl,dzf,rho12,psi,u,v,w,en12)
    !
    ! computes mean kinetic energy of the system
    !
    implicit none
    integer , intent(in ), dimension(3)        :: n
    integer , intent(in) :: nh_v,nh_psi,nh_g
    real(rp), intent(in ), dimension(3)        :: l,dl
    real(rp), intent(in ), dimension(1-nh_g:)       :: dzf
    real(rp), intent(in ), dimension(2)        :: rho12
    real(rp), intent(in ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(out), dimension(2)        :: en12
    real(rp) :: vcell,ecell,en1,en2
    integer :: i,j,k
    !
    en1 = 0._rp
    en2 = 0._rp
    !$acc data copy(en1,en2) async(1)
    !$acc parallel loop collapse(3) default(present) reduction(+:en1,en2) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          vcell = dl(1)*dl(2)*dzf(k)
          ecell = 0.5*((0.5*(u(i,j,k)+u(i-1,j,k)))**2+(0.5*(v(i,j,k)+v(i,j-1,k)))**2+(0.5*(w(i,j,k)+w(i,j,k-1)))**2)
          en1 = en1 + rho12(1)*psi(i,j,k)*vcell*ecell
          en2 = en2 + rho12(2)*(1.-psi(i,j,k))*vcell*ecell
        end do
      end do
    end do
    !$acc end data
    !$acc wait(1)
    en12(:) = [en1,en2]/product(l(:))
    call MPI_ALLREDUCE(MPI_IN_PLACE,en12,2,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  end subroutine cmpt_mean_energy
end module mod_output
