program pimc
    implicit none
    double precision, allocatable :: r(:,:,:), Utot(:), Umat(:,:,:), &
        Ueff(:,:), rho(:), rhocum(:), Ueffavg(:), ringdist(:), &
        Uring(:)
    
    double precision :: beta, betaeff, kT, K, xstep, rstep, dr, Z, &
        Zcum, logw, A, rn, mass = 2.0
    double precision, parameter :: rmin=2.0, rmax=10.0, &
        rw=6.0, kw=0.0

    integer :: Natoms=2, Nbeads, Npts=100
    integer :: nxacc=0, nxtrials=0, nracc=0, nrtrials=0, NMC=46*10**6, &
        mcblen, nmcblocks=2
    
    character(256) :: arg
    
    integer :: i, j, kk
    
    call get_command_argument(1, arg)
    read (arg, *) kT
    call get_command_argument(2, arg)
    read (arg, *) Nbeads
    
    allocate(r(3,Natoms,Nbeads), Utot(Nbeads), ringdist(Nbeads), &
        rho(Npts), rhocum(Npts), Ueff(Npts,nmcblocks), Ueffavg(Npts), &
        Uring(Nbeads))

    mcblen = NMC/nmcblocks
    beta = 1.0/kT    
    betaeff = beta/Nbeads
    mass = mass * 0.0206148690370024d0
    K = mass / betaeff**2
    
    
    
    xstep = sqrt( -2.0 * log(0.5)/(K*betaeff) )
    rstep = xstep
    dr = (rmax-rmin)/(Npts-1)
    
    r = 0.0
    call random_number(rn)
    r(1,2,1) = rw + rn*5.0
    
    call initialize_ring(r, xstep)

    do j=1,Nbeads
        Utot(j) = SilveraGoldman(r(:,:,j))
        Uring(j) = 0.5*K*sum((r(:,:,j) - r(:,:,ring_index(j+1)))**2)
    end do

    logw = -beta*kw*(r(1,2,1) - rw)**2
    logw = -betaeff*kw*sum((sqrt(sum((r(:,2,:) - r(:,1,:))**2,1)) - rw)**2)

    do i=1,1000000
       call mc_disp()
    end do

    rho = 0.0
    rhocum=0.0
    Z = 0.0
    Zcum = 0
    do i=1,NMC
        if (mod(i,Nbeads*100) == 0) then
            call mc_disp_ring()
        else
            call mc_disp()
        end if
        ringdist = sqrt( sum((r(:,2,:) - r(:,1,:))**2, 1) )
        do j = 1,Nbeads
            kk = int( (ringdist(j) - rmin + 0.5*dr)/dr ) + 1
            kk =max(kk, 1)
            
            if (kk > Npts) cycle

            rho(kk) = rho(kk) + exp(-logw) ! 1.0/w = exp(-logw)
        end do
        !write(*,*) nxacc, nxtrials, xstep

        Z = Z + exp(-logw)
        

        if (mod(i,100000) == 0) then
            write (*,*) rstep
            !rewind(69)
            !write(69,'(6F12.6)') r
            !flush(69)
        end if
        
        
        
        if (mod(i,mcblen) == 0) then
            rhocum = rhocum + rho
            Zcum = Zcum + Z
            rho = rho / (Z*Nbeads)
            Ueff(:,(i-1)/mcblen+1) = -1.0/beta * log(rho+1d-10)
            rho = 0.0
            Z = 0
        end if
    end do
    
    write(arg,'("Ueff_PIMC_kT=",F5.2,"_N=",I2,"_all.dat")') kT, Nbeads
    open(33,file=arg,status='REPLACE')
    write (33, '(3F12.6)') (rmin + dr*(i-1), Ueff(i,:), i=1,Npts)
    close(33)
    
    A = 0
    do j=1,Npts
        A = A + rhocum(j)*(rmin + dr*(j-1))**2
    end do
    A = A * dr * 4*3.141592653
    rhocum = rhocum/A
    Ueffavg = -1.0/beta * log(rhocum+1d-10)
    
    write(arg,'("Ueff_PIMC_kT=",F5.2,"_N=",I2,".dat")') kT, Nbeads
    open(33,file=arg,status='REPLACE')
    write (33, '(3F12.6)') (rmin + dr*(i-1), Ueffavg(i), rhocum(i), i=1,Npts)
    close(33)
contains

pure function ring_index(idx)
    integer, intent(in) :: idx
    integer :: ring_index

    ring_index = idx
    if (idx == 0) then
        ring_index = Nbeads
    else if (idx == Nbeads + 1) then
        ring_index = 1
    end if
end function

pure function SilveraGoldman(rv)
    double precision, intent(in) :: rv(3,2)
    double precision :: SilveraGoldman
    
    double precision, parameter :: palpha = 14.3757840987608, &
        pbeta = 2.961819, &
        pgamma = 0.035470953, &
        pC6 = 84105.6402716801, &
        pC8 = 417373.713330885, &
        pC9 = 146845.504557468, &
        pC10 = 2613703.27624221, &
        prc = 4.4021164021164

    double precision :: rsq, r, LJ, fc
    
    rsq = sum((rv(:,2) - rv(:,1))**2)
    r = sqrt(rsq)
    
    LJ = (pC6*rsq + pC8) / rsq**4 - (pC9 * r  - pC10)/rsq**5

    fc = 1.0
    if ( r <= prc ) then
        fc = exp( - (prc / r - 1.0)**2 )
    end if
    
    SilveraGoldman = exp(palpha - pbeta*r - pgamma*rsq) - LJ*fc
end function


subroutine initialize_ring(ring, xstep)
    double precision :: ring(:,:,:), xstep
    double precision :: phi
    integer :: i
    
    
    do i=2,size(ring,3)
        phi = 2.0 * (i-1) * acos(-1.0)/size(ring,3)
        ring(1,:,i) = ring(1,:,1)
        ring(2,:,i) = ring(2,:,1) + xstep*cos(phi)
        ring(3,:,i) = ring(3,:,1) + xstep*sin(phi)
    end do
end subroutine

subroutine mc_disp()
    implicit none
    double precision :: rn, rold(3), Unew, dr(3), logp, logwn, Uringnew(2)
    integer :: j, p, ibead
    integer, parameter :: nstepcheck = 200
    
    nxtrials = nxtrials + 1
    call random_number(rn)
    j = int(rn * (Nbeads*Natoms - 1) ) + 2
    p = (j-1)/Nbeads + 1
    ibead = mod(j, Nbeads) + 1
    
    call random_number(dr)
    if (ibead == 1 .and. p == 2) then
        dr(2:3) = 0.5d0
    end if
    rold = r(:,p,ibead)
    r(:,p,ibead) = r(:,p,ibead) + xstep * (dr - 0.5)

    logwn = -betaeff*kw*sum((sqrt(sum((r(:,2,:) - r(:,1,:))**2,1)) - rw)**2)

    if (r(1,2,1) <= 1.75*rmax) then
        Unew = SilveraGoldman(r(:,:,ibead))
        Uringnew(1) = 0.5*K*sum((r(:,:,ibead) - r(:,:,ring_index(ibead-1)))**2)
        Uringnew(2) = 0.5*K*sum((r(:,:,ibead) - r(:,:,ring_index(ibead+1)))**2)
        logp = -betaeff*(Unew + sum(Uringnew) -  Utot(ibead) &
                        - Uring(ring_index(ibead - 1)) - Uring(ibead) ) 
        call random_number(rn)
        if ((logp + logwn - logw ) > log(rn)) then
            Utot(ibead) = Unew
            Uring(ring_index(ibead - 1)) = Uringnew(1)
            Uring(ibead) = Uringnew(2)
            nxacc = nxacc + 1
            logw = logwn
        else
            r(:,p,ibead) = rold
        end if
    else
        r(:,p,ibead) = rold
    end if
    
    if (mod(nxtrials, nstepcheck) == 0) then
        if (nxacc > (0.4 * nstepcheck)) then
            xstep = xstep*1.08351d0
        elseif (nxacc < (0.2 * nstepcheck)) then
            xstep = xstep/1.04792d0
        end if
        !write (*,*) 'xstep =', xstep
        nxacc = 0
    end if
end subroutine


subroutine mc_disp_ring()
    implicit none
    double precision :: rn, rold(3,2,Nbeads), Unew(Nbeads), logp, logwn
    integer :: j
    integer, parameter :: nstepcheck = 200
    
    nrtrials = nrtrials + 1
    call random_number(rn)
    rold = r
    r(1,2,:) = r(1,2,:) + rstep*(rn-0.5)

    logwn = -betaeff*kw*sum((sqrt(sum((r(:,2,:) - r(:,1,:))**2,1)) - rw)**2)

    if (r(1,2,1) <= 1.25*rmax) then
        do j=1,Nbeads
            Unew(j) = SilveraGoldman(r(:,:,j))
        end do
        logp = -betaeff*(sum(Unew)  -  sum(Utot))
        call random_number(rn)
        if ((logp + logwn - logw ) > log(rn)) then
            Utot = Unew
            nracc = nracc + 1
            logw = logwn
        else
            r = rold
        end if
    else
        r = rold
    end if
    
    if (mod(nrtrials, nstepcheck) == 0) then
        if (nracc > (0.4 * nstepcheck)) then
            rstep = rstep*1.08351d0
        elseif (nracc < (0.2 * nstepcheck)) then
            rstep = rstep/1.04792d0
        end if
        rstep = min(rstep, 0.25*(rmax-rmin))
        !write (*,*) 'xstep =', xstep
        nracc = 0
    end if
end subroutine

end program pimc
