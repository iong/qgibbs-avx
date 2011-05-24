program dimer
    use vgw
    use vgwm
    use vgwfm
    implicit none
    integer :: i, Npts=150, natom=2
    double precision :: r(3,3), U(4), rmin=2.8, rmax=15, mass=2.0, kT=20.0, RC=30.0
    double precision, allocatable :: Y(:), FX(:,:), InvMeff(:,:,:), &
        SqrtMeff(:,:,:), Ucorr(:)
    double precision :: LNP, ENRG
    
    allocate(Y(1+21*Natom), FX(3,Natom), InvMeff(3,3,Natom), &
        SqrtMeff(3,3,Natom), Ucorr(Npts))
    
    call vgwinit(size(r,2), 'pH2-4g', RCUTOFF=RC)
    write (*,*) 'NOWW =========='
    call vgwfminit(natom, 'pH2-4g', RCUTOFF=RC)

    call INITIALIZE_VGWM(Natom,1d-4,1d-5,0d0,0d0)
    
    open(33,file='VGW_FM_SP.dat')
    
    r=0; U=0
    
    do i = 1,Npts
        r(1,2) = rmin + real(i-1)*(rmax - rmin)/real(Npts-1)
        call vgw0(r(:,1:2), 100d0, 1d0/kT, 0d0, U(1))

        call vgwquenchspb(natom,MASS,r(:,1:2),FX(:,1:2),LNP,U(2),&
                      ENRG,1d0/kT,1d-4,100d0,1d-4,RC,Y, InvMeff,SqrtMeff)

        call VGWMQUENCH(r(:,1:2),1d0/kT,U(3))

        call vgw0fm(r(:,1:2), 100d0, 1d0/kT, 0d0, U(4))
         

        Ucorr(i) = U(3) - U(2)
        write (33, '(5F12.6)') r(1,2), U
    end do
    close(33)
    call CLEANUP_VGWM()
    call vgwfmcleanup()
    
    Natom = 3
    call INITIALIZE_VGWM(Natom,1d-4,1d-5,0d0,0d0)
    call vgwfminit(natom, 'pH2-4g', RCUTOFF=RC)
    open(33,file='VGW_FM_SP_3b.dat')
    do i = 1,Npts
        r(1,2) = rmin + real(i-1)*(rmax - rmin)/real(Npts-1)
        r(1,3) = 0.5d0*r(1,2)
        r(2,3) = sqrt(3d0) * r(1,3)
        
        call vgw0(r, 100d0, 1d0/kT, 0d0, U(1))
        
        call vgwquenchspb(Natom,MASS,r,FX,LNP,U(2),&
                      ENRG,1d0/kT,1d-4,100d0,1d-4,RC,Y, InvMeff,SqrtMeff)

        call VGWMQUENCH(r,1d0/kT,U(3))
        
        call vgw0fm(r, 100d0, 1d0/kT, 0d0, U(4))
        
        write (33, '(6F12.6)') r(1,2), U, U(2) + 3.0*Ucorr(i)
    end do
    close(33)
    
end program dimer
