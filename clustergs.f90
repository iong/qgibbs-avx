program clustergs
    use vgw
    use vgwspfm
    use vgwfm
    use utils, only: M_PI
    implicit none

    character(256) :: arg, coords, waste

    integer :: i, Natom=180, Npts=20

    double precision :: endtime(3), begtime(3), mass, taustop
    double precision ::  deBoerMin=0.025, deBoerMax=0.35, &
        kT=0.001d0, RC=100d0, BL=100d0, sigma0 = 2.749, epsilon0 = 35.6d0

    double precision, allocatable :: r(:,:), r0(:,:), deBoer(:), U(:,:), Uinf(:)

!    if (command_argument_count() < 2) then
!        write (2,*) 'The program needs two argumens: <N> <coords.dat>'
!    end if
!
!    call get_command_argument(1, arg)
!    read(arg, *) Natom
!
!    allocate(r0(3,Natom))
!
!    call get_command_argument(2, coords)
!    open(33, file=trim(coords))
!    read(33, *) r0
!    close(33)

    call get_command_argument(1, coords)

    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read(arg, *) deBoerMin
        deBoerMax = deBoerMin
        Npts = 1
    end if

    open(33, file=trim(coords))
    read (33, *) Natom


    allocate(r0(3,Natom), Uinf(Npts), U(4,Npts), deBoer(Npts))
    U = 0

    read(33,*) U(1,1) !deBoer(1), U(1,1)
    read(33, *) (waste, r0(:,i), i=1,Natom)
    close(33)
    r0 = r0 / sigma0



    !call vgwfminit(natom, 'LJ')
    !U(1,:) = classical_Utot(r0, BL)
    !call vgwfmcleanup()
    U(1,:) = U(1,1)

    if (Npts>1) then
	    do i = 1,Npts
	        deBoer(i) = exp(log(deBoerMin) &
	            + log(deBoerMax/deBoerMin)*real(i-1)/real(Npts-1))
	        !deBoer(i) = deBoerMin + (deBoerMax-deBoerMin)*real(i-1)/real(Npts-1)
	        Uinf(i) = 1.5*Natom*kT
	    end do
	else
        deBoer = deBoerMin
    end if


    call cpu_time(begtime(1))
    do i = 1,Npts
        taustop = (deBoer(Npts)/deBoer(i))/kT
        mass = 1/deBoer(i)**2
        write (*,*) taustop
        call vgwinit(natom, 'LJ', massx=mass)
        call vgw0(r0, BL, taustop, U(2, i))
        call vgwcleanup()

        call vgwfminit(natom, 'LJ', massx=mass)
        call vgw0fm(r0, BL, taustop, U(3, i))
        call vgwfmcleanup()

        call vgwspfminit(natom, 'LJ', RCUTOFF=100d0, massx=mass)
        call vgw0spfm(r0, BL, taustop, U(4, i))
        call vgwspfmcleanup()
    end do
    call cpu_time(endtime(1))


    U(2,:) = U(2,:) - Uinf!+ U(1,1) - U(2,1)
    U(3,:) = U(3,:) - Uinf!+ U(1,1) - U(3,1)
    U(4,:) = U(4,:) - Uinf!+ U(1,1) - U(4,1)

    U( (/ 2, 3, 4 /), :) = U( (/ 2, 3, 4 /), :) * epsilon0

    write (*,*) endtime - begtime

    !if (Npts==1) then
    !    write (*, '(5F15.6)') deBoer(1), U(:,1)
    !else
	    arg=coords(1:len_trim(coords)-4) // '_groundstate.dat'
	    open(33,file=trim(arg),status='REPLACE')
	    write (33, '(5F15.6)') (deBoer(i), U(:,i), i=1,Npts)
	    close(33)
    !end if

end program clustergs
