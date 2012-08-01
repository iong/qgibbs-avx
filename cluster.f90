program cluster
    use vgw
    use vgwref, only:vgwrefinit, vgwref0, vgwrefcleanup
    use xyz
    implicit none
    integer, parameter :: Npts=20
    integer :: i, N
    double precision :: U(Npts, 2)=0d0, kT=.1, BL=100.0, deBoer=0.078, beta
    double precision, allocatable :: r(:,:), r0(:,:)
    double precision :: LNP, ENRG, sf, endtime(2), begtime(2)
    character(256) :: arg
    

    call get_command_argument(1, arg)
    call load_xyz(arg, r0)
    N = size(r0, 2)

    allocate(r(3,N))
    
    call get_command_argument(2, arg)
    read(arg, *) deBoer

    call get_command_argument(3, arg)
    read(arg, *) kT
    beta = 1d0/kT

   
    call vgwinit(N, 'LJ', massx=1d0/deBoer**2)
    call cpu_time(begtime(1))
    do i = 1,Npts
        sf = 0.95 + (1.5-0.95)*real(i-1)/real(Npts-1)
        call vgw0(r0 * sf, bl, beta, U(i, 1))
    end do
    call cpu_time(endtime(1))
    call vgwcleanup()

!   call vgwrefinit(N, 'LJ', massx=1d0/deBoer**2)
!   call cpu_time(begtime(2))
!   do i = 1,Npts
!       sf = 0.95 + (1.5-0.95)*real(i-1)/real(Npts-1)
!       r = r0*sf
!       call vgwref0(r, bl, beta, U(i, 2))
!   end do
!   call cpu_time(endtime(2))
!   call vgwrefcleanup()

    write (*,*) endtime - begtime
    open(33,file='Uref.dat')
    write(33, '(3F16.8)') (0.95 + (1.5-0.95)*real(i-1)/real(Npts-1), U(i,:), i=1,Npts)
    close(33)
    
end program cluster
