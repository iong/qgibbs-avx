program pimcsc
    use utils
    implicit none
    double precision, allocatable :: r(:,:,:), str(:,:,:), stK(:), &
        Ueff(:,:), rho(:), rhocum(:), Ueffavg(:), &
        ringdist(:), rold(:,:)
    
    double precision :: beta, betaeff, kT, K, xstep, rstep, dr, Z, &
        Zcum, A, rn, Utot, mass = 2.0, w, deBoer=0.3
    double precision, parameter :: rmin=.8, rmax=5.0

    integer :: Natoms=2, Nbeads, Npts=100
    integer :: nxacc=0, nxtrials=0, nracc=0, nrtrials=0, NMC=4*10**5, &
        mcblen, nmcblocks=2
    
    character(256) :: arg
    
    integer :: i, j, kk

    dr=0.05
    Npts = int((rmax - rmin) / dr) + 1
    
    call get_command_argument(1, arg)
    read (arg, *) kT
    call get_command_argument(2, arg)
    read (arg, *) Nbeads
    
    allocate(r(3,Natoms,Nbeads), str(3,Natoms,Nbeads), stK(Nbeads), &
        ringdist(Nbeads), rold(3,Nbeads),&
        rho(Npts), rhocum(Npts), Ueff(Npts,nmcblocks), Ueffavg(Npts))

    mcblen = NMC/nmcblocks
    beta = 1.0/kT    
    !mass = mass * 0.0206148690370024d0

    mass = 1/deBoer**2
    
    
    str=0d0
    stK(1) = 0.0

    call random_number(rn)
    str(1,2,1) = (rmax -rmin)*rn + rmin

    do j=2,Nbeads
        stK(j) = mass* kT**2 *real(j*Nbeads)/real(j-1)
        do i=1,3
            str(i,1,j) = gaussran(1d0/sqrt(beta*stK(j)), 0d0)
            str(i,2,j) = gaussran(1d0/sqrt(beta*stK(j)), 0d0)
        end do
    end do
    write(*,*) stK
    
    xstep = 1.0/sqrt(beta*stK(Nbeads))
    rstep = xstep

    call staging_to_real(str, r)
    Utot = stUtot(r)

    do i=1,1000000
        if (mod(i,2*Nbeads-1) == 0) then
            call mc_disp_ring()
        else
            call mc_disp()
        end if
    end do

    rho = 0.0
    rhocum=0.0
    Z = 0.0
    Zcum = 0
    
    do i=1,NMC*Nbeads*Natoms
        if (mod(i,2*Nbeads-1) == 0) then
            call mc_disp_ring()
        else
            call mc_disp()
        end if
        !call staging_to_real(str, r)
        ringdist = sqrt( sum((r(:,2,:) - r(:,1,:))**2, 1) )
        w = ringdist(1)**2
        do j = 1,Nbeads
            kk = int( (ringdist(j) - rmin + 0.5*dr)/dr ) + 1
            if (kk < 1) cycle
            if (kk > Npts) cycle

            rho(kk) = rho(kk) + w ! 1.0/w = exp(-logw)
        end do
        !write(*,*) nxacc, nxtrials, xstep

        Z = Z + w
        

        if (mod(i,100000) == 0) then
            !write (*,*) xstep, rstep, i
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
    
    !write(arg,'("Ueff_PIMC_kT=",F5.2,"_N=",I2,"_all.dat")') kT, Nbeads
    !open(33,file=arg,status='REPLACE')
    !write (33, '(3F12.6)') (rmin + dr*(i-1), Ueff(i,:), i=1,Npts)
    !close(33)
    
    A = 0
    A = sum(rhocum) * dr * 4*3.141592653
    do j=1,Npts
        !A = A + 4*3.141592653*dr*rhocum(j)*(rmin + dr*(j-1))**2
        rhocum(j) = rhocum(j) /(rmin + dr*(j-1))**2 
    end do
    rhocum = rhocum/A
    Ueffavg = -kT * log(rhocum)
    Ueffavg = Ueffavg - Ueffavg(Npts)
    
    write(arg,'("Ueff_PIMC_kT=",F5.2,"_N=",I2,".dat")') kT, Nbeads
    open(33,file=arg,status='REPLACE')
    do i=1,Npts
        if (rhocum(i) == 0.0) cycle
        write (33, '(3F12.6)') rmin + dr*(i-1), Ueffavg(i), rhocum(i)
    end do
    close(33)
contains

subroutine staging_to_real(q, x)
    implicit none
    double precision, intent(in) :: q(:,:,:)
    double precision, intent(out) :: x(:,:,:)
    integer :: Nb, i, j

    Nb = size(q, 3)

    x = spread(q(:,:,1), 3, Nb)
    do i=2,Nb
        do j=i,Nb
            x(:,:,i) = x(:,:,i) + real(i-1)/real(j-1)*q(:,:,j)
        end do
    end do
end subroutine

subroutine staging_to_real_p(dq, jq, r)
    double precision, intent(in) :: dq(3)
    double precision, intent(out) :: r(:,:)
    integer, intent(in) :: jq

    integer :: i

    if (jq == 1) then
        do i=1,size(r, 2)
            r(:,i) = r(:,i) + dq
        end do
    else
        do i=2,jq
            r(:,i) = r(:,i) + real(i-1)/real(jq-1)*dq
        end do
    end if
end subroutine


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

function LennardJones(rv)
    double precision, intent(in) :: rv(3,2)
    double precision :: LennardJones

    double precision :: rsq, y6

    rsq = sum((rv(:,2) - rv(:,1))**2)

    y6 = 1d0/rsq**3

    LennardJones = 4d0*(y6**2 - y6)
end function


function stUtot(r)
    double precision, intent(in) :: r(:,:,:)
    double precision :: stUtot
    integer :: j
   
    stUtot = 0.0
    do j=1,Nbeads
        stUtot = stUtot + LennardJones(r(:,:,j))
    end do
end function


subroutine mc_disp()
    implicit none
    double precision :: rn, strold(3), Unew, dq(3), logp
    integer :: j, p, ibead
    integer, parameter :: nstepcheck = 200
    
    nxtrials = nxtrials + 1
    call random_number(rn)
    j = int(rn * ((Nbeads-1)*Natoms) )
    p = j/(Nbeads-1) + 1
    ibead = mod(j, Nbeads-1) + 2
    
    call random_number(dq)
    dq = xstep * (dq  - 0.5)
    strold = str(:,p,ibead)
    str(:,p,ibead) = str(:,p,ibead) + dq

    rold = r(:,p,:)
    call staging_to_real_p(dq, ibead, r(:,p,:))

    if (minval(sum((r(:,2,:) - r(:,1,:))**2, 1)) <= 1.5*rmax**2) then
        Unew = stUtot(r)
        logp = -beta*( (Unew - Utot)/real(Nbeads) &
            + 0.5*stK(ibead)*sum(str(:,p,ibead)**2 - strold**2) )

        call random_number(rn)
        if (logp > log(rn)) then
            Utot = Unew
            nxacc = nxacc + 1
        else
            str(:,p,ibead) = strold
            r(:,p,:) = rold
        end if
    else
        str(:,p,ibead) = strold
        r(:,p,:) = rold
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
    double precision :: dq, rn, qold, Unew, logp
    integer :: j
    integer, parameter :: nstepcheck = 200
    
    nrtrials = nrtrials + 1
    call random_number(dq)
    dq = rstep*(dq - 0.5)

    qold = str(1,2,1)
    rold = r(:,2,:)

    do j=1,Nbeads
        r(1,2,j) = r(1,2,j) + dq
    end do

    if (minval(sum((r(:,2,:) - r(:,1,:))**2, 1)) <= 1.5*rmax**2) then
        Unew = stUtot(r)

        logp = -beta*(Unew - Utot)/real(Nbeads)

        call random_number(rn)
        if (logp > log(rn)) then
            Utot = Unew
            nracc = nracc + 1
        else
            str(1,2,1) = qold
            r(:,2,:) = rold
        end if
    else
        str(1,2,1) = qold
        r(:,2,:) = rold
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
end program pimcsc
