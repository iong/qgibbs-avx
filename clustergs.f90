program clustergs
!    use vgw
    use vgw
    use vgwspfm
    use vgwfm
    implicit none

    character(256) :: arg, coords

    integer, parameter :: Npts=10
    integer :: i, Natom=180

    double precision :: U(4, Npts), endtime(4), begtime(4)
    double precision :: deBoer, deBoerMin=0.01, deBoerMax=0.35, kT=0.1d0, &
        RC=100d0, BL=100d0

    double precision, allocatable :: r(:,:), r0(:,:)

    if (command_argument_count() < 2) then
        write (2,*) 'The program needs two argumens: <N> <coords.dat>'
    end if

    call get_command_argument(1, arg)
    read(arg, *) Natom
    
    allocate(r0(3,Natom))

    call get_command_argument(2, coords)
    open(33, file=trim(coords))
    read(33, *) r0
    close(33)
    
    !if (command_argument_count() >= 3) then
    !    call get_command_argument(3, arg)
    !    read(arg, *) deBoer
    !    call vgwspfminit(natom, 'LJ', massx=1/deBoer**2)
    !end if

    !call vgwinit(natom, 'pH2-4g', RCUTOFF=RC)


    U = 0

    call vgwfminit(natom, 'LJ')
    U(1,:) = classical_Utot(r0(:,1:2))
    call vgwfmcleanup()


    call cpu_time(begtime(1))
    do i = 1,Npts
        deBoer = deBoerMin + (deBoerMax-deBoerMin)*real(i-1)/real(Npts-1)

        write (*,*) i, deBoer

        call vgwinit(natom, 'LJ', massx=1/deBoer**2)
        call vgw0(r0, BL, 1d0/kT, U(2, i))
        call vgwcleanup()
    end do
    call cpu_time(endtime(1))

    call cpu_time(begtime(2))
    do i = 1,Npts
        deBoer = deBoerMin + (deBoerMax-deBoerMin)*real(i-1)/real(Npts-1)

        write (*,*) i, deBoer

        call vgwfminit(natom, 'LJ', massx=1d0/deBoer**2)
        call vgw0fm(r0, BL, 1d0/kT, U(3, i))
        call vgwfmcleanup()
    end do
    call cpu_time(endtime(2))

    call cpu_time(begtime(3))
    do i = 1,Npts
        deBoer = deBoerMin + (deBoerMax-deBoerMin)*real(i-1)/real(Npts-1)

        write (*,*) i, deBoer

        call vgwspfminit(natom, 'LJ', massx=1/deBoer**2)
        call vgw0spfm(r0, BL, 1d0/kT, U(4, i))
        call vgwspfmcleanup()
    end do
    call cpu_time(endtime(3))

    
    write (*,*) endtime - begtime

    arg=coords(1:len_trim(coords)-4) // '_groundstate.dat'
    open(33,file=trim(arg),status='REPLACE')
    write (33, '(5F15.6)') (&
        deBoerMin + (deBoerMax-deBoerMin)*real(i-1)/real(Npts-1), U(:,i), &
        i=1,Npts)
    close(33)
end program clustergs
