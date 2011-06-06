program mccluster
    use vgw
    use vgwspfm
    use vgwm
    implicit none
    integer, parameter :: Npts=20, natom=55, NMC=100
    integer :: i
    double precision :: U0(4), Unew(4), rmin=2.8, rmax=15, mass=2.0, kT=20.0, RC=100.0
    double precision, allocatable :: r(:,:), Y(:), FX(:,:), InvMeff(:,:,:), &
        SqrtMeff(:,:,:)
    double precision :: LNP, ENRG, sf, et(2), bt(2), xstep=0.75, rjold(3), dr(3), rn, r10(3), p
    
    allocate(Y(1+21*Natom), FX(3,Natom), InvMeff(3,3,Natom), &
        SqrtMeff(3,3,Natom), r(3,Natom))
    
    call vgwinit(natom, 'pH2-4g', RCUTOFF=RC)
    call vgwspfminit(natom, 'pH2-4g', RCUTOFF=RC)
    !call INITIALIZE_VGWM(Natom,1d-4,1d-4,0d0,0d0, 'fmcore55')
    
    open(33, file='H2_55_Lowest.dat')
    read(33, *) r
    close(33)

    U0 = 0
    
    call INITIALIZE_VGWM(Natom,1d-4,1d-4,0d0,0d0)
    call VGWMQUENCH(r, 1d0/kT, U0(2))
    call CLEANUP_VGWM()
    
    !call vgw0spfm(r, 100d0, 1d0/kT, Unew(3))
    
    call INITIALIZE_VGWM(Natom,1d-4,1d-4,0d0,0d0, 'fmcore55')
    call VGWMQUENCH(r, 1d0/kT, U0(4))
    call CLEANUP_VGWM()
    
    r10 = r(:,1)
    open(33,file='VGW_FM_SP_55core_MCWALK.dat',status='REPLACE')
    do i = 1,NMC
        write (*,*) i
        call random_number(dr)
        
        rjold = r(:,1)
        r(:,1) = r(:,1) + xstep*(dr-0.5)
        if (sum((r(:,1) - r10)**2) > 3.99d0**2) then
            cycle
        end if

        !call vgwquenchspb(Natom,MASS,r,FX,LNP,U(2),&
        !              ENRG,1d0/kT,1d-4,100d0,1d-4,RC,Y, InvMeff,SqrtMeff)

        call INITIALIZE_VGWM(Natom,1d-4,1d-4,0d0,0d0)
        call VGWMQUENCH(r, 1d0/kT, Unew(2))
        call CLEANUP_VGWM()
        
        !call vgw0spfm(r, 100d0, 1d0/kT, Unew(3))
        
        call INITIALIZE_VGWM(Natom,1d-4,1d-4,0d0,0d0, 'fmcore55')
        call VGWMQUENCH(r, 1d0/kT, Unew(4))
        call CLEANUP_VGWM()
        
        p = exp( - ( Unew(2) - U0(2) ) / kT)
        call random_number(rn)
        
        write (33, '(I10,7F15.6)') i, Unew - U0, xstep, rn, p
        
        
        if ( p > rn) then
            U0 = Unew
        else
            r(:,1) = rjold
        end if
    end do
    close(33)
end program mccluster
