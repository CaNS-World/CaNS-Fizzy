! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_output_acdi
  use mod_output, only: out0d
  use mod_param, only: datadir
  use mpi
  use mod_common_mpi, only:ierr,myid
  use mod_types
  implicit none
  private
  public out0d_acdi,out1d_acdi_channel
  contains
  !
  subroutine out0d_acdi(n,dl,dzf,psi,time)
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: dzf
    real(rp), intent(in), dimension(0:,0:,0:) :: psi
    real(rp), intent(in) :: time
    real(rp) :: vcell
    real(rp) :: vpsi
    integer :: i,j,k
    real(rp), dimension(2) :: var
    !
    ! Calculate the total volume of the phase field psi
    !
    vpsi = 0._rp
    !$acc data copy(vpsi) async(1)
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          vcell = dl(1)*dl(2)*dzf(k)
          vpsi = vpsi + psi(i,j,k)*vcell
        end do
      end do
    end do
    !$acc end data
    !$acc wait(1)
    call MPI_ALLREDUCE(MPI_IN_PLACE,vpsi,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    var(1) = time
    var(2) = vpsi
    call out0d(trim(datadir)//'acdi_volume.out',2,var)
  end subroutine out0d_acdi
  !
  subroutine out1d_acdi_channel(fname,ng,lo,hi,l,dl,dzc_g,dzf_g,zc_g,zf_g,u,v,w,p,psi)
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,lo,hi
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: dzc_g,dzf_g,zc_g,zf_g
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w,p,psi
    real(rp), allocatable, dimension(:,:) :: buf
    real(rp) :: tmp_x,tmp_y,tmp_z
    integer :: i,j,k
    integer :: iunit
    integer :: nvars
    character(len=30) cfmt
    real(rp) :: grid_area_ratio
    real(rp) :: buf01,buf02,buf03,buf04,buf05,buf06,buf07,buf08,buf09,buf10, &
                buf11,buf12,buf13,buf14,buf15,buf16,buf17,buf18,buf19,buf20, &
                buf21,buf22,buf23,buf24,buf25,buf26,buf27,buf28,buf29,buf30, &
                buf31,buf32,buf33,buf34,buf35,buf36,buf37,buf38,buf39,buf40, &
                buf41,buf42,buf43,buf44,buf45,buf46
    !
    nvars = 36
    grid_area_ratio = dl(1)*dl(2)/(l(1)*l(2))
    allocate(buf(nvars,ng(3)))
    !$acc enter data create(buf) copyin(dzc_g,dzf_g,zc_g,zf_g)
    !$acc kernels default(present)
    buf(:,:) = 0._rp
    !$acc end kernels
    !$acc parallel loop gang default(present)
    do k=lo(3),hi(3)
      buf01 = 0._rp
      buf02 = 0._rp
      buf03 = 0._rp
      buf04 = 0._rp
      buf05 = 0._rp
      buf06 = 0._rp
      buf07 = 0._rp
      buf08 = 0._rp
      buf09 = 0._rp
      buf10 = 0._rp
      buf11 = 0._rp
      buf12 = 0._rp
      buf13 = 0._rp
      buf14 = 0._rp
      buf15 = 0._rp
      buf16 = 0._rp
      buf17 = 0._rp
      buf18 = 0._rp
      buf19 = 0._rp
      buf20 = 0._rp
      buf21 = 0._rp
      buf22 = 0._rp
      buf23 = 0._rp
      buf24 = 0._rp
      buf25 = 0._rp
      buf26 = 0._rp
      buf27 = 0._rp
      buf28 = 0._rp
      buf29 = 0._rp
      buf30 = 0._rp
      buf31 = 0._rp
      buf32 = 0._rp
      buf33 = 0._rp
      buf34 = 0._rp
      buf35 = 0._rp
      buf36 = 0._rp
      !buf37 = 0._rp
      !buf38 = 0._rp
      !buf39 = 0._rp
      !buf40 = 0._rp
      !buf41 = 0._rp
      !buf42 = 0._rp
      !buf43 = 0._rp
      !buf44 = 0._rp
      !buf45 = 0._rp
      !buf46 = 0._rp
      !$acc loop vector collapse(2)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          !
          ! phase field
          !
          buf01 = buf01 + psi(i,j,k)
          ! 
          buf02 = buf02 + psi(i,j,k)**2
          !
          ! velocity
          !
          buf03 = buf03  + u(i,j,k)
          buf04 = buf04  + v(i,j,k)
          buf05 = buf05  + 0.5*(w(i,j,k)+w(i,j,k-1))
          buf06 = buf06  + 0.5*(psi(i,j,k)+psi(i+1,j,k))*u(i,j,k)
          buf07 = buf07  + 0.5*(psi(i,j,k)+psi(i,j+1,k))*v(i,j,k)
          buf08 = buf08  + 0.5*psi(i,j,k)*(w(i,j,k)+w(i,j,k-1))
          !            
          buf09 = buf09  + u(i,j,k)**2
          buf10 = buf10  + v(i,j,k)**2
          buf11 = buf11  + 0.25*(w(i,j,k)+w(i,j,k-1))**2
          buf12 = buf12  + 0.5*(psi(i,j,k)+psi(i+1,j,k))*u(i,j,k)**2
          buf13 = buf13  + 0.5*(psi(i,j,k)+psi(i,j+1,k))*v(i,j,k)**2
          buf14 = buf14  + 0.25*psi(i,j,k)*(w(i,j,k)+w(i,j,k-1))**2
          buf15 = buf15  + 0.25*(u(i-1,j,k)+u(i,j,k))*(w(i,j,k-1)+w(i,j,k))
          buf16 = buf16  + 0.25*psi(i,j,k)*(u(i-1,j,k)+u(i,j,k))*(w(i,j,k)+w(i,j,k-1))
          !            
          buf17 = buf17 + u(i,j,k)**3
          buf18 = buf18 + v(i,j,k)**3
          buf19 = buf19 + 0.125*(w(i,j,k)+w(i,j,k-1))**3
          buf20 = buf20 + 0.5*(psi(i,j,k)+psi(i+1,j,k))*u(i,j,k)**3
          buf21 = buf21 + 0.5*(psi(i,j,k)+psi(i,j+1,k))*v(i,j,k)**3
          buf22 = buf22 + 0.125*psi(i,j,k)*(w(i,j,k)+w(i,j,k-1))**3
          !
          ! pressure
          !
          buf23 = buf23 + p(i,j,k)
          !
          buf24 = buf24 + p(i,j,k)**2
          !
          ! vorticity
          !
          tmp_x = (w(i,j+1,k)-w(i,j,k))/dl(2) - (v(i,j,k+1)-v(i,j,k))/dzc_g(k)
          tmp_y = (u(i,j,k+1)-u(i,j,k))/dzc_g(k) - (w(i+1,j,k)-w(i,j,k))/dl(1)
          tmp_z = (v(i+1,j,k)-v(i,j,k))/dl(1) - (u(i,j+1,k)-u(i,j,k))/dl(2)
          !
          buf25 = buf25 + tmp_x
          buf26 = buf26 + tmp_y
          buf27 = buf27 + tmp_z
          buf28 = buf28 + 0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i,j+1,k+1)+psi(i,j,k+1))*tmp_x
          buf29 = buf29 + 0.25*(psi(i,j,k)+psi(i+1,j,k)+psi(i+1,j,k+1)+psi(i,j,k+1))*tmp_y
          buf30 = buf30 + 0.25*(psi(i,j,k)+psi(i+1,j,k)+psi(i+1,j+1,k)+psi(i,j+1,k))*tmp_z
          !            
          buf31 = buf31 + tmp_x**2
          buf32 = buf32 + tmp_y**2
          buf33 = buf33 + tmp_z**2
          buf34 = buf34 + 0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i,j+1,k+1)+psi(i,j,k+1))*tmp_x**2
          buf35 = buf35 + 0.25*(psi(i,j,k)+psi(i+1,j,k)+psi(i+1,j,k+1)+psi(i,j,k+1))*tmp_y**2
          buf36 = buf36 + 0.25*(psi(i,j,k)+psi(i+1,j,k)+psi(i+1,j+1,k)+psi(i,j+1,k))*tmp_z**2
          !
          ! scalar
          !
          !buf37 = buf37 + s(i,j,k)
          !buf38 = buf38 + psi(i,j,k)*s(i,j,k)
          !!            
          !buf39 = buf39 + s(i,j,k)**2
          !buf40 = buf40 + psi(i,j,k)*s(i,j,k)**2
          !buf41 = buf41 + s(i,j,k)*0.5*(u(i,j,k)+u(i-1,j,k))
          !buf42 = buf42 + s(i,j,k)*0.5*(v(i,j,k)+v(i,j-1,k))
          !buf43 = buf43 + s(i,j,k)*0.5*(w(i,j,k)+w(i,j,k-1))
          !buf44 = buf44 + psi(i,j,k)*s(i,j,k)*0.5*(u(i,j,k)+u(i-1,j,k))
          !buf45 = buf45 + psi(i,j,k)*s(i,j,k)*0.5*(v(i,j,k)+v(i,j-1,k))
          !buf46 = buf46 + psi(i,j,k)*s(i,j,k)*0.5*(w(i,j,k)+w(i,j,k-1))
          !
        end do
      end do
      buf( 1,k) = buf01*grid_area_ratio
      buf( 2,k) = buf02*grid_area_ratio
      buf( 3,k) = buf03*grid_area_ratio
      buf( 4,k) = buf04*grid_area_ratio
      buf( 5,k) = buf05*grid_area_ratio
      buf( 6,k) = buf06*grid_area_ratio
      buf( 7,k) = buf07*grid_area_ratio
      buf( 8,k) = buf08*grid_area_ratio
      buf( 9,k) = buf09*grid_area_ratio
      buf(10,k) = buf10*grid_area_ratio
      buf(11,k) = buf11*grid_area_ratio
      buf(12,k) = buf12*grid_area_ratio
      buf(13,k) = buf13*grid_area_ratio
      buf(14,k) = buf14*grid_area_ratio
      buf(15,k) = buf15*grid_area_ratio
      buf(16,k) = buf16*grid_area_ratio
      buf(17,k) = buf17*grid_area_ratio
      buf(18,k) = buf18*grid_area_ratio
      buf(19,k) = buf19*grid_area_ratio
      buf(20,k) = buf20*grid_area_ratio
      buf(21,k) = buf21*grid_area_ratio
      buf(22,k) = buf22*grid_area_ratio
      buf(23,k) = buf23*grid_area_ratio
      buf(24,k) = buf24*grid_area_ratio
      buf(25,k) = buf25*grid_area_ratio
      buf(26,k) = buf26*grid_area_ratio
      buf(27,k) = buf27*grid_area_ratio
      buf(28,k) = buf28*grid_area_ratio
      buf(29,k) = buf29*grid_area_ratio
      buf(30,k) = buf30*grid_area_ratio
      buf(31,k) = buf31*grid_area_ratio
      buf(32,k) = buf32*grid_area_ratio
      buf(33,k) = buf33*grid_area_ratio
      buf(34,k) = buf34*grid_area_ratio
      buf(35,k) = buf35*grid_area_ratio
      buf(36,k) = buf36*grid_area_ratio
      !buf(37,k) = buf37*grid_area_ratio
      !buf(38,k) = buf38*grid_area_ratio
      !buf(39,k) = buf39*grid_area_ratio
      !buf(40,k) = buf40*grid_area_ratio
      !buf(41,k) = buf41*grid_area_ratio
      !buf(42,k) = buf42*grid_area_ratio
      !buf(43,k) = buf43*grid_area_ratio
      !buf(44,k) = buf44*grid_area_ratio
      !buf(45,k) = buf45*grid_area_ratio
      !buf(46,k) = buf46*grid_area_ratio
    end do
    !$acc update self(buf)
    call MPI_ALLREDUCE(MPI_IN_PLACE,buf(1,1),size(buf),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    if(myid == 0) then
      nvars = 36
      write(cfmt,'(A,I3,A)') '(',nvars+2+2,'ES26.18)'
      open(newunit=iunit,file=fname//'.out')
      do k=1,ng(3)
        write(iunit,trim(cfmt)) zc_g(k),zf_g(k),(buf(i,k),i=1,nvars),dzc_g(k),dzf_g(k)
      end do
      close(iunit)
    end if
    !$acc exit data delete(buf)
    !$acc wait
  end subroutine out1d_acdi_channel
end module mod_output_acdi
