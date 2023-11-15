! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_output
  use mpi
  use decomp_2d_io
  use mod_common_mpi, only:ierr,myid
  use mod_types
  implicit none
  private
  public out0d,gen_alias,out1d,out1d_chan,out2d,out3d,write_log_output,write_visu_2d,write_visu_3d
  contains
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
    if (myid  ==  0) then
      open(newunit=iunit,file=fname,position='append')
      write(iunit,'(*(E16.7e3))') var(1:n)
      close(iunit)
    end if
  end subroutine out0d
  !
  subroutine gen_alias(myid,datadir,fname,fname_alias)
    !
    ! this subroutine generates a symlink with name `fname_alias`, pointing to
    ! file `datadir//fname` using the `execute_command_line` intrinsic;
    ! it is called by task `myid`
    !
    integer, intent(in) :: myid
    character(len=*), intent(in) :: datadir,fname,fname_alias
    if(myid == 0) call execute_command_line('ln -sf '//trim(fname)//' '//trim(datadir)//fname_alias)
  end subroutine gen_alias
  !
  subroutine out1d(fname,ng,lo,hi,idir,l,dl,z_g,dz,p)
    !
    ! writes the profile of a variable averaged
    ! over two domain directions
    !
    ! fname -> name of the file
    ! ng    -> global domain sizes
    ! lo,hi -> upper and lower extents of the input array
    ! idir  -> direction of the profile
    ! l,dl  -> uniform grid spacing and length arrays
    ! z_g   -> global z coordinate array (grid is non-uniform in z)
    ! dz    -> local z grid spacing array (should work also with the global one)
    ! p     -> 3D input scalar field
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: z_g,dz
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    real(rp), allocatable, dimension(:) :: p1d
    integer :: i,j,k
    integer :: iunit
    real(rp) :: grid_area_ratio,p1d_s
    !
    allocate(p1d(ng(idir)))
    !$acc enter data create(p1d)
    !$acc kernels default(present)
    p1d(:) = 0._rp
    !$acc end kernels
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
          write(iunit,'(2E16.7e3)') z_g(k),p1d(k)
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
          write(iunit,'(2E16.7e3)') (j-.5)*dl(2),p1d(j)
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
          write(iunit,'(2E16.7e3)') (i-.5)*dl(1),p1d(i)
        end do
        close(iunit)
      end if
    end select
  end subroutine out1d
  !
  subroutine out2d(fname,inorm,islice,p)
    use mod_common_mpi, only: ipencil => ipencil_axis
    !
    ! saves a planar slice of a scalar field into a binary file
    !
    ! fname  -> name of the output file
    ! inorm  -> plane is perpendicular to direction
    !           inorm (1,2,3)
    ! islice -> plane is of constant index islice
    !           in direction inorm
    ! p      -> 3D input scalar field
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in) :: inorm,islice
    real(rp), intent(in), dimension(:,:,:) :: p
    !
    select case(inorm)
    case(1) !normal to x --> yz plane
       call decomp_2d_write_plane(ipencil,p,inorm,islice,'.',fname,'dummy')
    case(2) !normal to y --> zx plane
       call decomp_2d_write_plane(ipencil,p,inorm,islice,'.',fname,'dummy')
    case(3) !normal to z --> xy plane
       call decomp_2d_write_plane(ipencil,p,inorm,islice,'.',fname,'dummy')
    end select
  end subroutine out2d
  !
  subroutine out3d(fname,nskip,p)
    use mod_common_mpi, only: ipencil => ipencil_axis
    !
    ! saves a 3D scalar field into a binary file
    !
    ! fname  -> name of the output file
    ! nskip  -> array with the step size for which the
    !           field is written; i.e.: [1,1,1]
    !           writes the full field
    ! p      -> 3D input scalar field
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: nskip
    real(rp), intent(in), dimension(:,:,:) :: p
    integer :: fh
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp
    !
    call MPI_FILE_OPEN(MPI_COMM_WORLD, fname, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)
    disp = 0_MPI_OFFSET_KIND
#if !defined(_OPENACC)
    call decomp_2d_write_every(ipencil,p,nskip(1),nskip(2),nskip(3),fname,.true.)
#else
    !
    ! alternative way of writing a full 3D field that was needed in the past
    !
    block
      use decomp_2d
      use mod_load, only: io_field
      integer, dimension(3) :: ng,lo,hi
      ng(:) = [nx_global,ny_global,nz_global]
      select case(ipencil)
      case(1)
        lo(:) = xstart(:)
        hi(:) = xend(:)
      case(2)
        lo(:) = ystart(:)
        hi(:) = yend(:)
      case(3)
        lo(:) = zstart(:)
        hi(:) = zend(:)
      end select
      if(any(nskip /= 1) .and. myid == 0) &
        print*, 'Warning: `nskip` should be `[1,1,1]` if `io_field()` is used to output 3D field data'
      call io_field('w',fh,ng,[0,0,0],lo,hi,disp,p)
    end block
#endif
    call MPI_FILE_CLOSE(fh,ierr)
  end subroutine out3d
  !
  subroutine write_log_output(fname,fname_fld,varname,nmin,nmax,nskip,time,istep)
    !
    ! appends information about a saved binary file to a log file
    ! this file is used to generate a xdmf file for visualization of field data
    !
    ! fname     -> name of the output log file
    ! fname_fld -> name of the saved binary file (excluding the directory)
    ! varname   -> name of the variable that is saved
    ! nmin      -> first element of the field that is saved in each direction, e.g. [1,1,1]
    ! nmax      -> last  element of the field that is saved in each direction, e.g. [ng(1),ng(2),ng(3)]
    ! nskip     -> step size between nmin and nmax, e.g. [1,1,1] if the whole array is saved
    ! time      -> physical time
    ! istep     -> time step number
    !
    implicit none
    character(len=*), intent(in) :: fname,fname_fld,varname
    integer , intent(in), dimension(3) :: nmin,nmax,nskip
    real(rp), intent(in)               :: time
    integer , intent(in)               :: istep
    character(len=100) :: cfmt
    integer :: iunit
    !
    write(cfmt, '(A)') '(A,A,A,9i5,E16.7e3,i7)'
    if (myid  ==  0) then
      open(newunit=iunit,file=fname,position='append')
      write(iunit,trim(cfmt)) trim(fname_fld),' ',trim(varname),nmin,nmax,nskip,time,istep
      close(iunit)
    end if
  end subroutine write_log_output
  !
  subroutine write_visu_3d(datadir,fname_bin,fname_log,varname,nmin,nmax,nskip,time,istep,p)
    !
    ! wraps the calls of out3d and write_log_output into the same subroutine
    !
    implicit none
    character(len=*), intent(in)          :: datadir,fname_bin,fname_log,varname
    integer , intent(in), dimension(3)    :: nmin,nmax,nskip
    real(rp), intent(in)                  :: time
    integer , intent(in)                  :: istep
    real(rp), intent(in), dimension(:,:,:) :: p
    !
    call out3d(trim(datadir)//trim(fname_bin),nskip,p)
    call write_log_output(trim(datadir)//trim(fname_log),trim(fname_bin),trim(varname),nmin,nmax,nskip,time,istep)
  end subroutine write_visu_3d
  !
  subroutine write_visu_2d(datadir,fname_bin,fname_log,varname,inorm,nslice,ng,time,istep,p)
    !
    ! wraps the calls of out2d and write-log_output into the same subroutine
    !
    implicit none
    character(len=*), intent(in)          :: datadir,fname_bin,fname_log,varname
    integer , intent(in)                  :: inorm,nslice
    integer , intent(in), dimension(3)    :: ng
    real(rp), intent(in)                  :: time
    integer , intent(in)                  :: istep
    real(rp), intent(in), dimension(:,:,:) :: p
    integer , dimension(3) :: nmin_2d,nmax_2d
    !
    call out2d(trim(datadir)//trim(fname_bin),inorm,nslice,p)
    select case(inorm)
    case(1)
      nmin_2d(:) = [nslice,1    ,1    ]
      nmax_2d(:) = [nslice,ng(2),ng(3)]
    case(2)
      nmin_2d(:) = [1    ,nslice,1    ]
      nmax_2d(:) = [ng(1),nslice,ng(3)]
    case(3)
      nmin_2d(:) = [1    ,1    ,nslice]
      nmax_2d(:) = [ng(1),ng(2),nslice]
    end select
    call write_log_output(trim(datadir)//trim(fname_log),trim(fname_bin),trim(varname),nmin_2d,nmax_2d,[1,1,1],time,istep)
  end subroutine write_visu_2d
  !
  subroutine out1d_chan(fname,ng,lo,hi,idir,l,dl,z_g,u,v,w) ! e.g. for a channel with streamwise dir in x
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: z_g
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w
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
          write(iunit,'(8E16.7e3)') z_g(k),um(k),vm(k),wm(k), &
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
  subroutine out2d_duct(fname,ng,lo,hi,idir,l,dl,z_g,u,v,w) ! e.g. for a duct
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: z_g
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w
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
      um(:,:) =      um(:,:)*grid_area_ratio
      vm(:,:) =      vm(:,:)*grid_area_ratio
      wm(:,:) =      wm(:,:)*grid_area_ratio
      u2(:,:) = sqrt(u2(:,:)*grid_area_ratio - um(:,:)**2)
      v2(:,:) = sqrt(v2(:,:)*grid_area_ratio - vm(:,:)**2)
      w2(:,:) = sqrt(w2(:,:)*grid_area_ratio - wm(:,:)**2)
      vw(:,:) =      vw(:,:)*grid_area_ratio - vm(:,:)*wm(:,:)
      uv(:,:) =      uv(:,:)*grid_area_ratio - um(:,:)*vm(:,:)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do k=1,ng(3)
          do i=1,ng(1)
            x_g = (i-.5)*dl(1)
            write(iunit,'(10E16.7e3)') x_g,z_g(k),um(i,k),vm(i,k),wm(i,k), &
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
      call MPI_ALLREDUCE(MPI_IN_PLACE,um(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vm(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,wm(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,u2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,v2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,w2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,uv(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,uw(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      um(:,:) =      um(:,:)*grid_area_ratio
      vm(:,:) =      vm(:,:)*grid_area_ratio
      wm(:,:) =      wm(:,:)*grid_area_ratio
      u2(:,:) = sqrt(u2(:,:)*grid_area_ratio - um(:,:)**2)
      v2(:,:) = sqrt(v2(:,:)*grid_area_ratio - vm(:,:)**2)
      w2(:,:) = sqrt(w2(:,:)*grid_area_ratio - wm(:,:)**2)
      uv(:,:) =      uv(:,:)*grid_area_ratio - um(:,:)*vm(:,:)
      uw(:,:) =      uw(:,:)*grid_area_ratio - um(:,:)*wm(:,:)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do k=1,ng(3)
          do j=1,ng(2)
            y_g = (j-.5)*dl(2)
            write(iunit,'(10E16.7e3)') y_g,z_g(k),um(j,k),vm(j,k),wm(j,k), &
                                                  u2(j,k),v2(j,k),w2(j,k), &
                                                  uv(j,k),uw(j,k)
          end do
        end do
        close(iunit)
      end if
    end select
  end subroutine out2d_duct
  !
  subroutine velstats(fname,n,dl,l,zc,u,v,w,p,s,psi)
    !
    ! computes a few volume-averaged profiles
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer, intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: dl,l
    real(rp), intent(in), dimension(0:) :: zc
    real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w,p,s,psi
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
        write(iunit,'(25E16.7e3)') zc(k), &
                                 u1_1(k),v1_1(k),w1_1(k),p1_1(k),s1_1(k),c1_1(k), &
                                 u1_2(k),v1_2(k),w1_2(k),p1_2(k),s1_2(k),c1_2(k), &
                                 u2_1(k),v2_1(k),w2_1(k),p2_1(k),s2_1(k),c2_1(k), &
                                 u2_2(k),v2_2(k),w2_2(k),p2_2(k),s2_2(k),c2_2(k)
      end do
      close(iunit)
    end if
  end subroutine velstats
  !
  subroutine cmpt_total_area(n,dli,dzci,dzfi,psi,area)
    !
    ! computes some volume-averaged profiles along z (assumes constant grid spacing)
    !
    implicit none
    integer , intent(in ), dimension(3) :: n
    real(rp), intent(in ), dimension(3) :: dli
    real(rp), intent(in ), dimension(0:) :: dzci,dzfi
    real(rp), intent(in ), dimension(0:,0:,0:) :: psi
    real(rp), intent(out) :: area
    real(rp), dimension(8) :: mx,my,mz
    real(rp) :: dpsidx,dpsidy,dpsidz
    integer  :: i,j,k
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
end module mod_output
