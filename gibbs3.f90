program gibbs3
    use utils
    use pairint
    implicit none
    real*8 :: rhcore
    real*8 :: Vtot, V(2), bl(2), kT, beta, U0(4), xstep(2), &
            logVstep, Vstep
    real*8 :: Z, Um(2), Vm(2), Nm(2), pm(2), rho0(2), rhom(2), mum(2)
    real*8, allocatable :: rs(:,:, :), mcblock(:,:)
    real*8, allocatable, target :: Umat(:,:,:), Umatnew(:,:,:)
    real*8, allocatable :: Ucacheline(:,:)
    integer :: Ntot, N(2), Nswap=0, Nvol
    integer :: npickatt(2), npickattmax, nswapacc = 0, nmum(2)
    integer :: nxtrials(2)=0, nxacc(2)=0, nvoltrials=0, nvolacc=0
    integer :: imc, NMC=10000000, jmc, Nburn=5000000, mcblen=1000000, ib
    integer :: logfd=31
    character(LEN=256) :: arg

    class(pair_interactions), pointer :: pi
    type(SilveraGoldman), allocatable, target :: sg
    type(LennardJones), allocatable, target :: lj

    real*8 :: rn

    call get_command_argument(1, arg)
    read (arg, *) N(1)
    call get_command_argument(2, arg)
    read (arg, *) N(2)
    call get_command_argument(3, arg)
    read (arg, *) rho0(1)
    call get_command_argument(4, arg)
    read (arg, *) rho0(2)
    call get_command_argument(5, arg)
    read (arg, *) kT

    if (command_argument_count() >= 6) then
        call get_command_argument(6, arg)
        read(arg, *) Nswap
    end if

    write(arg, '("kT=",F5.2,".log")') kT
    open(logfd,file=trim(arg))

    Ntot = sum(N)
    allocate(Umat(Ntot, Ntot, 4), Umatnew(Ntot, Ntot, 4), rs(3, Ntot, 2), &
        Ucacheline(Ntot,2), y6cache(Ntot), mcblock(14, NMC/mcblen))

    V = N/rho0
    Vtot = sum(V)
    Vstep=0.01*minval(V)

    bl = V**(1d0/3d0)

    if (minval(bl) < 5.0) then
        write (*,'("Box #",I1," is too small. L = ",F8.4, &
            " < 5*sigma = ",F8.4)') minloc(bl), minval(bl), 5.0
        write (*,*) 'Consider increasing the number of particles or decreasing the density'
    end if

    Nvol = 1!minval(N)/20
    write (*,*) 'Nvol =', Nvol

    beta = 1.0/kT

    if (kT>= 20.0) then
        allocate(sg)
        pi => sg
        rhcore = 2.4
        xstep = 3.0/bl
    else
        allocate(lj)
        pi => lj
        rhcore = 1.0/(0.5 + 0.5*sqrt(1 - kT*log10(1d-15)))**(1d0/6d0)
        xstep = 0.3/bl
    end if




    call populate_cube2(rs(:,1:N(1), 1))
    call populate_cube2(rs(:,1:N(2), 2))
    call init_interaction_matrix(bl, U0)

    do imc=1,1000000
        call random_number(rn)
        jmc=int(rn*Ntot) + 1
        call mc_move(jmc)
        call cumulate()
        if (mod(imc,100000) == 0) then
            write (logfd,'("NMC =",I10,", N =",2I5,", V =",2G12.5,", rho =",2F9.6,", U =",2G12.5,", p =",2F8.4)') &
                imc, N, V, N/V, U0(1:2), pm/Z
        end if
    end do

    Z = 0; Um = 0; Nm = 0; pm = 0; Vm = 0; rhom = 0; mum = 0; nmum = 0
    do imc=1,NMC+Nburn
        call random_number(rn)
        jmc=1+int(rn*(Ntot + Nswap + Nvol))
        if (jmc<= Ntot) then
            call mc_move(jmc)
        elseif (jmc <= Ntot + Nswap) then
            call mc_swap()
        elseif (jmc <= Ntot + Nswap + Nvol) then
            call mc_vol()
        endif
        if (mod(imc,100000) == 0) then
            write (logfd,'("NMC =",I10,", N =",2I5,", V =",2G12.5,", rho =",&
                2F9.6,", U =",2G12.5,", p =",2F8.4,", mu =",2F8.2,&
                ", swapacc = ",F8.5)') &
                imc, N, V, N/V, U0(1:2), pm/Z, -kT*log(mum/nmum), &
                real(nswapacc*(Nswap+Nvol+Ntot))/(Z*real(Ntot))
            call dump_xyz(imc)
        end if
        call cumulate()
        if (imc == Nburn) then
            Z = 0; Um = 0; Nm = 0; pm = 0; Vm = 0; rhom = 0; mum = 0; nmum = 0
        end if
        if (imc <= Nburn) cycle

        if (mod(imc-Nburn,mcblen) == 0) then
            ib = (imc - Nburn) / mcblen
            mcblock(1, ib) = kT
            mcblock(2:3, ib) = rhom/Z
            mcblock(4:5, ib) = -kT*log(mum/real(nmum))
            mcblock(6:7, ib) = pm/Z
            mcblock(8:9, ib) = Um/Z
            mcblock(10:11, ib) = Nm/Z
            mcblock(12:13, ib) = Vm/Z
            mcblock(14, ib) = real(nswapacc*(Ntot+Nvol+Nswap))/(Z*Ntot)

            call dump_avg(ib)
            Z = 0; Um = 0; Nm = 0; pm = 0; Vm = 0; rhom = 0; mum = 0; nmum = 0
            nswapacc = 0
        end if
    enddo

    close(logfd)
contains

subroutine populate_cube2(rs)
    real*8, intent(out) :: rs(:,:)
    logical, allocatable :: occupied(:)
    integer :: np, nu, i, j, lattpos
    np = size(rs, 2)
    nu = ceiling(real(np)**(1.0/3.0))

    allocate(occupied(nu*nu*nu))
    occupied = .FALSE.

    do i=1,np
        do
            call random_number(rn)
            lattpos = int(rn*nu**3)
            if (.not.occupied(lattpos+1)) exit
        end do
        occupied(lattpos+1) = .TRUE.

        do j=1,3
            rs(j, i) = real(mod(lattpos, nu))
            lattpos = lattpos / nu
        end do
    end do

    rs = (rs+0.5)/real(nu)

    deallocate(occupied)
end subroutine

pure function too_close(rsj, rs, blj)
    real*8, intent(in) :: rsj(3), rs(:,:), blj
    logical :: too_close
    real*8 :: drsq, rshc_sq, drs(3)
    integer :: i

    too_close = .FALSE.
    rshc_sq = (rhcore / blj)**2
    do i=1,size(rs, 2)
        drs = rsj - rs(:,i)
        drsq = sum(min_image(drs)**2)
        if (drsq < rshc_sq) then
            too_close = .TRUE.
            return
        end if
    enddo
end function


subroutine mc_move(jin)
    implicit none
    integer, intent(in) :: jin
    integer, parameter :: nstepcheck=200
    real*8 :: rsj(3), dU, dr(3), logp, dvir
    integer :: ibox, j, i

    ibox = 1
    j = jin
    if (j > N(1)) then
        ibox = 2
        j = j - N(1)
    end if

    nxtrials(ibox) = nxtrials(ibox) + 1

    call random_number(dr)
    rsj = rs(:, j, ibox) + xstep(ibox) * (dr-0.5)
    rsj = rsj - floor(rsj)

    call dU_particle_move(j, ibox, rsj, dU)
    logp = -beta * dU
    call random_number(rn)
    if (logp>log(rn)) then
        nxacc(ibox) = nxacc(ibox)  + 1
        rs(:, j, ibox) = rsj
        call dU_accept_move(j, ibox, dvir)
        U0(ibox) = U0(ibox) + dU
        U0(ibox+2) = U0(ibox+2) + dvir
    end if

    if (mod(nxtrials(ibox), nstepcheck) == 0) then
        if (nxacc(ibox) > (0.4 * nstepcheck)) then
            xstep(ibox) = xstep(ibox)*1.08351d0
        elseif (nxacc(ibox) < (0.2 * nstepcheck)) then
            xstep(ibox) = xstep(ibox)/1.04792d0
        end if
        xstep(ibox) = min(xstep(ibox), 0.25)
        xstep(ibox) = max(xstep(ibox), 0.1/bl(ibox))
        nxacc(ibox) = 0
    end if
end subroutine

subroutine mc_swap()
    implicit none
    integer :: isrc, idest, jkill
    real*8 :: pacc, pacc0, rsn(3), dvirsrc, dvirdest, dUdest, dUsrc, rc

    call random_number(rn)
    isrc = 1 + int(2.0*rn)
    idest = 3 - isrc

    if (N(isrc) == 0) return

    call random_number(rsn)
    if (too_close(rsn, rs(:,1:N(idest), idest), bl(idest)) ) return

    call random_number(rn)
    jkill = 1 + int(rn*N(isrc))

    call dU_test_particle(idest, rsn, dUdest)
    mum(idest) = mum(idest) + V(idest)*exp(-beta*dUdest)/(N(idest)+1)
    nmum(idest) = nmum(idest) + 1

    pacc0 = ( V(idest) * real(N(isrc)) ) / ( V(isrc) * real(N(idest) + 1) )


    dUsrc = sum(Umat(:, jkill, isrc)) + pi%lrcorr(N(isrc), V(isrc), 0.5*bl(isrc)) &
             - pi%lrcorr(N(isrc)-1, V(isrc), 0.5*bl(isrc))

    pacc = pacc0 * exp( - beta * (dUdest - dUsrc) )

    call random_number(rn)
    if (pacc > rn) then
        !write (*,*) 'swap'
        nswapacc = nswapacc + 1
        dvirsrc = sum(Umat(:, jkill, isrc+2))
        call dU_accept_particle_swap(jkill, isrc, idest, dvirdest)

        rs(:, jkill, isrc) = rs(:, N(isrc), isrc)
        N(isrc) = N(isrc) - 1

        rs(:, N(idest)+1, idest) = rsn
        N(idest) = N(idest) + 1

        U0(isrc) = U0(isrc) - dUsrc
        U0(idest) = U0(idest) + dUdest
        U0(isrc+2) = U0(isrc+2) - dvirsrc
        U0(idest+2) = U0(idest+2) + dvirdest
    end if
end subroutine


subroutine mc_vol()
    implicit none
    integer, parameter :: nstepcheck=50
    real*8 :: U0new(4), Vn(2), bln(2), logp

    nvoltrials = nvoltrials + 1

    call random_number(rn)
    Vn(1) = V(1) + Vstep*(rn-0.5)
    Vn(2) = Vtot - Vn(1)
    bln = Vn**(1.0/3.0)

    ! ensure bl > 5\sigma
    if (minval(bln) < 5.0) return

    U0new = Utot_vol_scale(bln)
    !write (*,'(4ES16.6)') U0new(1:2) - U0(1:2), N*log(Vn/V)
    logp = sum ( -beta * (U0new(1:2) - U0(1:2)) + real(N)*(log(Vn/V)) )
    call random_number(rn)
    if (logp>log(rn)) then
        nvolacc = nvolacc + 1
        V = Vn
        bl = bln
        U0 = U0new
        call dU_accept_vol_scale()
        !write (*,'("v",2I5,2F9.2,2F9.6,2F15.6,2F8.4)') N, V, N/V, U0(1:2), xstep, logVstep
    end if

    if (mod(nvoltrials, nstepcheck) == 0) then
        if (nvolacc > 0.55 * nstepcheck) then
            Vstep = Vstep*1.18351d0
        elseif (nvolacc < 0.45 * nstepcheck) then
            Vstep = Vstep/1.14792d0
        end if
        !Vstep = max(Vstep, 0.001*Vtot)
        nvolacc = 0
    end if
end subroutine


subroutine update_virial_cache(Nbox)
    integer, intent(in) :: Nbox
    Ucacheline(1:Nbox, 2) = 12.0*Ucacheline(1:Nbox, 1) &
                                + 24.0*y6cache(1:Nbox)
end subroutine

subroutine init_interaction_matrix(bln, Utot, Udestmat)
    implicit none
    real*8, intent(in) :: bln(2)
    real*8, intent(out) :: Utot(4)
    real*8, optional, target, intent(out) :: Udestmat(:,:,:)
    real*8, pointer :: Umatp(:,:,:)
    integer :: i, ibox, j

    Umatp => Umat
    if (present(Udestmat)) then
        Umatp => Udestmat
    end if

    Umatp = 0d0
    do ibox=1,2
        do i=1, N(ibox)-1
            call pi%Up(rs(:, i, ibox), rs(:,i+1:N(ibox), ibox), bln(ibox), &
                    Umatp(i+1:N(ibox), i, ibox), &
                    Umatp(i+1:N(ibox), i, ibox + 2))
        end do
    end do

    Utot = 0d0
    do i=1,4
        ibox = 1 + mod(i-1,2)
        do j=1,N(ibox)-1
            Utot(i) = Utot(i) + sum(Umatp(j+1:N(ibox), j, i))
        end do
    end do

    do i=1,4
        ibox = 1 + mod(i-1,2)
        do j=2,N(ibox)
            Umatp(1:j-1,j, i) = Umatp(j, 1:j-1, i)
        end do
    end do

    Utot(1) = Utot(1) + pi%lrcorr(N(1), bln(1)**3, 0.5*bln(1))
    Utot(2) = Utot(2) + pi%lrcorr(N(2), bln(2)**3, 0.5*bln(2))
end subroutine

subroutine dU_particle_move(j, ibox, rsj, du)
    implicit none
    integer, intent(in) :: j, ibox
    real*8, intent(in) :: rsj(3)
    real*8, intent(out) :: dU

    Ucacheline = 0d0
    call pi%Up(rsj, rs(:,1:N(ibox), ibox), bl(ibox), &
            Ucacheline(1:N(ibox), 1))
    Ucacheline(j,1) = 0d0

    dU =  sum(Ucacheline( 1:N(ibox), 1)) - sum(Umat( 1:N(ibox) , j, ibox))
end subroutine 

subroutine dU_accept_move(p, pbox, dvir)
    implicit none
    integer, intent(in) :: p, pbox
    real*8, intent(out) :: dvir
    integer :: i, ibox

    call update_virial_cache(N(pbox))
    Ucacheline(p,2) = 0.0
    dvir =  sum( Ucacheline( 1:N(pbox), 2) ) &
                                    - sum(Umat( 1:N(pbox) , p, pbox + 2))

    do i=1,2
        ibox = pbox + (i-1)*2
        Umat(1:N(pbox),p, ibox) = Ucacheline(1:N(pbox), i)
        Umat(p,1:N(pbox), ibox) = Ucacheline(1:N(pbox), i)
    end do
end subroutine

function Utot_vol_scale(bln) result(Unew)
    implicit none
    real*8, intent(in) :: bln(2)
    real*8 :: Unew(4)

    call init_interaction_matrix(bln, Unew, Umatnew)
end function

subroutine dU_accept_vol_scale()
    implicit none
    Umat = Umatnew
end subroutine

subroutine dU_test_particle(idest, rsj, dU)
    implicit none
    integer, intent(in) :: idest
    real*8, intent(in) :: rsj(3)
    real*8, intent(out) :: dU

    dU = 0d0
    if (N(idest) == 0) return

    Ucacheline = 0d0
    call pi%Up(rsj, rs(:,1:N(idest), idest), bl(idest), &
            Ucacheline(1:N(idest), 1))
    dU = sum(Ucacheline(1:N(idest),1 ))
    dU = dU + pi%lrcorr(N(idest)+1, V(idest), 0.5*bl(idest)) &
             - pi%lrcorr(N(idest), V(idest), 0.5*bl(idest))
end subroutine

subroutine dU_accept_particle_swap(p, psrc, pdest, dvirdest)
    implicit none
    integer, intent(in) :: p, psrc, pdest
    real*8, intent(out) :: dvirdest
    integer :: i, isrc, idest


    call update_virial_cache(N(pdest))
    dvirdest =  sum( Ucacheline( 1:N(pdest), 2) ) &
            - sum(Umat( 1:N(pdest) , p, pdest + 2))


    do i=1,2
        isrc = psrc + (i-1)*2
        idest = pdest + (i-1)*2

        Umat(1:N(psrc), p, isrc) = Umat(1:N(psrc),N(psrc), isrc) 
        Umat(p, 1:N(psrc), isrc) = Umat(1:N(psrc),N(psrc), isrc) 
        Umat(p, p, isrc) = 0d0
        Umat(1:N(psrc), N(psrc), isrc) = 0d0
        Umat(N(psrc), 1:N(psrc), isrc) = 0d0

        Umat(1:N(pdest), N(pdest)+1, idest) = Ucacheline(1:N(pdest), i)
        Umat(N(pdest)+1, 1:N(pdest), idest) = Ucacheline(1:N(pdest), i)
        Umat(N(pdest)+1, N(pdest)+1, idest) = 0d0
    end do
end subroutine

subroutine cumulate()
    integer :: i
    Z = Z + 1.0
    do i=1,2
        if (N(i) ==0) cycle
        Um(i) = Um(i) + U0(i)/N(i)
        Nm(i) = Nm(i) + N(i)
        Vm(i) = Vm(i) + V(i)
        rhom(i) = rhom(i) + N(i)/V(i)
        pm(i) = pm(i) + (N(i)*kT + U0(i+2)/3.0)/V(i)
    end do
end subroutine


subroutine dump_avg(nmcb)
    implicit none
    integer, intent(in) :: nmcb
    integer :: i
    character(256) :: fname

    write(fname, '("Um_",F5.2,"_NMC.dat")') kT
    open(30,file=trim(fname))
    write(30,'(I10, 14G16.7)') (i*mcblen+Nburn, mcblock(:,i), i=1,nmcb)
    close(30)
end subroutine

subroutine dump_xyz(nmcnow)
    implicit none
    integer, intent(in) :: nmcnow
    character(len=256) :: fname
    integer :: ibox, j

    do ibox=1,2
        write(fname, '("gibbs_box",I1,"_",I10,".xyz")') ibox, nmcnow
        call replace_char(fname, ' ', '0')
        open(34,file=trim(fname),status='REPLACE')
        write(34,*) N(ibox)
        write(34,*) bl(ibox)
        write(34,'(("X ", 3F13.8))') rs(:,1:N(ibox),ibox)
        close(34)
    end do
end subroutine
end program gibbs3
