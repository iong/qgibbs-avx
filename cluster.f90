program cluster
!    use vgw
    use vgwspfm
    use vgwm
    implicit none
    integer, parameter :: Npts=20, natom=180
    integer :: i
    double precision :: U(4, Npts), rmin=2.8, rmax=15, mass=2.0, kT=14.0, RC=100d0, BL=19.7122
    double precision, allocatable :: r(:,:), Y(:), FX(:,:), InvMeff(:,:,:), &
        SqrtMeff(:,:,:), r0(:,:)
    double precision :: LNP, ENRG, sf, endtime(2), begtime(2)
    
    allocate(Y(1+21*Natom), FX(3,Natom), InvMeff(3,3,Natom), &
        SqrtMeff(3,3,Natom), r(3,Natom), r0(3,Natom))
    
    !call vgwinit(natom, 'pH2-4g', RCUTOFF=RC)
    call vgwspfminit(natom, 'pH2-4g')
    

    open(33, file='pH2180.dat')
    read(33, *) r0
    close(33)

    U = 0

    call cpu_time(begtime(1))
    do i = 1,Npts
        sf = 0.95 + (1.5-0.95)*real(i-1)/real(Npts-1)

        write (*,*) i, sf
        r = r0 * sf
        !call vgw0v(r, 100d0, (/1d0/kT/), 0d0, U(1:1))

        !call vgwquenchspb(Natom,MASS,r,FX,LNP,U(2),&
        !              ENRG,1d0/kT,1d-4,100d0,1d-4,RC,Y, InvMeff,SqrtMeff)

        !call INITIALIZE_VGWM(Natom,1d-4,1d-4,sf*BL,8d0)
        !call VGWMQUENCH(r, 1d0/kT, U(3, i))
        !call CLEANUP_VGWM()

        write (*,*) U(3,i)
    end do
    call cpu_time(endtime(1))
     
    call cpu_time(begtime(2))
    do i = 1,Npts
       write (*,*) i, '/', Npts
       sf = 0.95 + (1.5-0.95)*real(i-1)/real(Npts-1)
       r = r0 * sf

       call vgw0spfm(r, sf*BL, 1d0/kT, U(4, i))
        write (*,*) U(3:4,i)
    end do
    call cpu_time(endtime(2))
    
    write (*,*) endtime - begtime
    
    open(33,file='FC-IonutFC_pH2_180.dat',status='REPLACE')
        write (33, '(5F15.6)') (0.95 + (1.5-0.95)*real(i-1)/real(Npts-1), &
                U(:,i), i=1,Npts)
    close(33)
end program cluster
