program clustergs
    use vgw
    use vgwspfm
    use vgwfm
    use utils, only: M_PI
    implicit none
!    include 'mkl_service.fi'

    character(256) :: arg, coords, waste

    integer :: i, Natom

    double precision :: endtime(3), begtime(3), mass, taustop
    double precision ::  deBoer=0.1, kT=0.01d0, RC=100d0, BL=100d0, sigma0 = 2.749, epsilon0 = 35.6d0

    double precision, allocatable :: r0(:,:)
    double precision :: U(4), Uinf(4), rt(3)

    call get_command_argument(1, coords)

    open(33, file=trim(coords))
    read (33, *) Natom


    allocate(r0(3,Natom))
    U = 0

    read(33,*) deBoer, U(1)

    read(33, *) (waste, r0(:,i), i=1,Natom)
    close(33)

    r0 = r0*1.05d0

    !r0 = r0 / sigma0

    !write (*,*) 'max threads =', mkl_domain_get_max_threads( MKL_BLAS )

    taustop = 1d0/kT
    mass = 1/deBoer**2
    call vgwinit(natom, 'LJ', RCUTOFF=100d0, massx=mass)
!        call vgw0(r0, BL, taustop, U(2, i), rt(1, i))
    call vgwcleanup()

    call vgwfminit(natom, 'LJ', massx=mass)
!    call vgw0fm(r0, BL, taustop, U(3), rt(2))
    call vgwfmcleanup()

    arg=coords(1:len_trim(coords)-4) // '_gs_tau.dat'
    open(60,file=trim(arg),status='REPLACE')
    call vgwspfminit(natom, 'LJ', RCUTOFF=100d0, massx=mass)
    arg=coords(1:len_trim(coords)-4) // '_final.xyz'
    call vgw0spfm(r0, BL, taustop, U(4), rt(3), 60, arg)
    call vgwspfmcleanup()
    close(60)

    U(2:4) = U(2:4) * epsilon0

    arg=coords(1:len_trim(coords)-4) // '_gs.dat'
    open(33,file=trim(arg),status='REPLACE')

    arg=coords(1:len_trim(coords)-4) // '_rt.dat'
    open(34,file=trim(arg),status='REPLACE')

    write (33, '(I8, 4F15.6)') Natom, U
    write (34, '(I8, 3F15.6)') Natom, rt

    close(33)
    close(34)

end program clustergs
