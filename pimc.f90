program pimc
    use pairint
    implicit none
    double precision, allocatable :: r(:,:,:), Utot(:), Umat(:,:,:), Ueff(:,:), &
        rho(:), rhocum(:), Ueffavg(:), ringdist(:), Ucacheline(:)
    
    double precision :: beta, betaeff, kT, K, xstep, rstep, Uring, dr
    double precision, parameter :: rmin=2.8, rmax=10.0, mass = 2.0
    
    integer :: Natoms=2, Nbeads, Npts=50
    integer :: nxacc=0, nxtrials=0, nracc=0, nrtrials=0, NMC=200000000, mcblen, nmcblocks=2
    
    character(256) :: arg
    
    integer :: i, j, kk
    
    class(pair_interactions), pointer :: pi
    type(SilveraGoldman), allocatable, target :: sg
    
    call get_command_argument(1, arg)
    read (arg, *) kT
    call get_command_argument(2, arg)
    read (arg, *) Nbeads
    
    allocate(r(3,Natoms,Nbeads), Utot(Nbeads), Umat(Natoms,Natoms,Nbeads), &
        Ueff(Npts,nmcblocks), rho(Npts), rhocum(Npts), Ueffavg(Npts), ringdist(Npts), &
        Ucacheline(Natoms), sg)

    pi => sg

    mcblen = NMC/nmcblocks
    beta = 1.0/kT    
    betaeff = beta/Nbeads
    K = mass / betaeff**2
    
    
    
    xstep = sqrt( -2.0 * log(0.5)/(K*betaeff) )
    rstep = 1.1
    dr = (rmax-rmin)/(Npts-1)
    
    r = 0.0
    r(1,2,1) = 0.75*(rmax-rmin)*rmin
    
    call initialize_ring(r, xstep)
    call init_interaction_matrix()
    rho = 0.0
    do i=1,NMC

        call mc_disp()
        ringdist = sqrt( sum((r(:,2,:) - r(:,1,:))**2, 1) )
        do j = 1,Nbeads
            kk = int( (ringdist(j) - rmin + 0.5*dr)/dr )
            kk = min(max(kk, 1), Npts)
            rho(kk) = rho(kk) + 1.0
        end do
        
        
        
        if (mod(i,mcblen) == 0) then
            rhocum = rhocum + rho
            rho = rho / (mcblen*Nbeads)
            Ueff(:,(i-1)/mcblen+1) = -1.0/beta * log(rho+1d-10)
            rho = 0.0
        end if
    end do
    
    write(arg,'("Ueff_PIMC_kT=",F5.2,"_all.dat")') kT
    open(33,file=arg,status='REPLACE')
    write (33, '(3F12.6)') (rmin + dr*(i-1), Ueff(i,:), i=1,Npts)
    close(33)
    
    rhocum = rhocum/(NMC*Nbeads)
    Ueffavg = -1.0/beta * log(rhocum+1d-10)
    
    !Ueffavg = sum(Ueff, 2)/nmcblocks
    !do i=1,Npts
    !    Ueffstd = Ueffstd + (Ueff(:,i) - Ueffavg)**2
    !end do
    !Ueffstd = sqrt(Ueffstd/nmcblocks)
    
    write(arg,'("Ueff_PIMC_kT=",F5.2,".dat")') kT
    open(33,file=arg,status='REPLACE')
    write (33, '(2F12.6)') (rmin + dr*(i-1), Ueffavg(i), i=1,Npts)
    close(33)
contains

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
    double precision :: rn, rj(3), dUring, dUpot, dr(3), logp
    integer :: j, p, ibead
    integer, parameter :: nstepcheck = 200
    
    nxtrials = nxtrials + 1
    call random_number(rn)
    j = int(rn * (Nbeads*Natoms - 1) ) + 2
    p = (j-1)/Nbeads + 1
    ibead = mod(j, Nbeads) + 1
    
    call random_number(dr)
    rj = r(:,p,ibead) + xstep * (dr - 0.5)
    if (ibead == 1 .and. p == 2) then
        rj = 0.0
        rj(1) = r(1,p,ibead) + xstep * (dr(1) - 0.5)
    end if
    
    if (ibead /= 1 .or. p /=2 .or. rj(1) <= rmax) then
        call dU_bead_disp(p, ibead, rj, dUpot, dUring)
        logp = -beta*(dUpot + dUring)
        call random_number(rn)
        if (logp > log(rn)) then
            call accept_disp(p, ibead, rj, dUpot, dUring)
            nxacc = nxacc + 1
        end if
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

subroutine mc_disp_ring(p)
    implicit none
    integer, intent(in) :: p
    double precision :: dr(3), rold(3,Nbeads), Umatold(Natoms,Natoms,Nbeads), &
        Utotold(Nbeads)
    integer :: i
    integer, parameter :: nstepcheck = 50

    nrtrials = nrtrials + 1
    rold = r(:,p,:)
    call random_number(dr)
    do i=1,Nbeads
        r(:,p,i) = r(:,p,i) + Nbeads*xstep*(dr-0.5)
    end do
    Utotold = Utot
    Umatold = Umat
    
    call init_interaction_matrix()
    call random_number(dr)
    if (-betaeff*sum(Utot-Utotold) > log(dr(1)) ) then
        nracc = nracc + 1
        !write (*,'(G12.5,8F8.3," xxx ",8F8.3)') rstep,Utotold, Utot
    else
        Utot = Utotold
        Umat = Umatold
        r(:,p,:) = rold
    end if
    
    if (mod(nrtrials, nstepcheck) == 0) then
        if (nracc > (0.4 * nstepcheck)) then
            rstep = rstep*1.08351d0
        elseif (nracc < (0.2 * nstepcheck)) then
            rstep = rstep/1.04792d0
        end if
        write (*,*) real(nracc)/real(nstepcheck), xstep*Nbeads
        nracc = 0
    end if
end subroutine
    
subroutine init_interaction_matrix()
    implicit none
    integer :: i, j

    do i=1,Nbeads
        do j=1,Natoms
            call pi%Up(r(:, j, i)/100d0, r(:,j+1:Natoms, i)/100d0, 100d0, &
                    Umat(j+1:Natoms, j, i))
        end do
    end do

    Utot = 0d0
    do i=1,Nbeads
        do j=1,Natoms-1
            Utot(i) = Utot(i) + sum(Umat(j+1:Natoms, j, i))
        end do
    end do

    do i=1,Nbeads
        do j=2,Natoms
            Umat(1:j-1,j, i) = Umat(j, 1:j-1, i)
        end do
    end do

    Uring =  0.5*K*sum( (r(:,:, 1) - r(:,:,Nbeads))**2 )
    do i=1,Nbeads-1
        Uring = Uring + 0.5*K*sum( (r(:,:, i+1) - r(:,:,i))**2 )
    end do
end subroutine

subroutine dU_bead_disp(j, ibead, rj, dUpot, dUring)
    implicit none
    integer, intent(in) :: j, ibead
    real*8, intent(in) :: rj(3)
    real*8, intent(out) :: dUpot, dUring

    Ucacheline = 0d0
    call pi%Up(rj/100d0, r(:,:, ibead)/100d0, 100d0, Ucacheline)
    Ucacheline(j) = 0d0

    dUpot =  sum(Ucacheline) - sum(Umat(: , j, ibead))
    
    if (ibead==Nbeads) then
        dUring = 0.5*K* ( sum( (rj - r(:,j,1))**2 ) &
                + sum( (rj - r(:,j,Nbeads-1))**2 ) &
                - sum( (r(:,j,ibead) - r(:,j,1))**2 ) &
                - sum( (r(:,j,ibead) - r(:,j,Nbeads-1))**2) )
    else if (ibead == 1) then
        dUring = 0.5*K* ( sum( (rj - r(:,j,Nbeads))**2 ) &
                + sum( (rj - r(:,j,2))**2 ) &
                - sum( (r(:,j,ibead) - r(:,j,Nbeads))**2 ) &
                - sum( (r(:,j,ibead) - r(:,j,2))**2) )
    else
        dUring = 0.5*K* ( sum( (rj - r(:,j,ibead-1))**2 ) &
                + sum( (rj - r(:,j, ibead+1))**2 ) &
                - sum( (r(:,j,ibead) - r(:,j, ibead-1))**2 ) &
                - sum( (r(:,j,ibead) - r(:,j, ibead+1))**2) )
    end if
end subroutine

subroutine accept_disp(p, ibead, rj, dUpot, dUring)
    implicit none
    double precision, intent(in) :: rj(3), dUpot, dUring
    integer, intent(in) :: p, ibead
    integer :: i, ibox

    Umat(:,p, ibead) = Ucacheline
    Umat(p,:, ibead) = Ucacheline
    Utot(ibead) = Utot(ibead) + dUpot
    Uring = Uring + dUring
    r(:,p,ibead) = rj
end subroutine

end program pimc