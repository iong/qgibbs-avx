program cluster
    use vgw
    use vgwm
    implicit none
    integer :: i, Npts=20, natom=54
    double precision :: U(4), rmin=2.8, rmax=15, mass=2.0, kT=20.0, RC=100.0
    double precision, allocatable :: r(:,:), Y(:), FX(:,:), InvMeff(:,:,:), &
        SqrtMeff(:,:,:), r0(:,:)
    double precision :: LNP, ENRG, sf
    
    allocate(Y(1+21*Natom), FX(3,Natom), InvMeff(3,3,Natom), &
        SqrtMeff(3,3,Natom), r(3,Natom), r0(3,Natom))
    
    call vgwinit(natom, 'pH2-4g', RCUTOFF=RC)

    open(33, file='H2_54surf_Lowest.dat')
    read(33, *) r0
    close(33)

    call INITIALIZE_VGWM(Natom,1d-4,1d-5,0d0,0d0)
    open(33,file='VGW_FM_SP_54_surf.dat',status='REPLACE')
    do i = 1,Npts
        sf = 0.8 + (1.5-0.8)*real(i-1)/real(Npts-1)
        r = r0 * sf
        call vgw0v(r, 100d0, (/1d0/kT/), 0d0, U(1:1))
        
        call vgwquenchspb(Natom,MASS,r,FX,LNP,U(2),&
                      ENRG,1d0/kT,1d-4,100d0,1d-4,RC,Y, InvMeff,SqrtMeff)
        
        call VGWMQUENCH(r, 1d0/kT, U(3))
        
        write (33, '(5F15.6)') sf, U
    end do
    close(33)
end program cluster
