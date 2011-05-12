program qgibbs
    use vgw
    use utils, only: min_image, replace_char
    implicit none
    real*8, parameter :: M_PI = 3.141592653d0
    real*8 :: rhcore
    real*8 :: Vtot, V(2), bl(2), kT, beta, U0(4), xstep(2), Vstep
    real*8 :: Z, Um(2), Vm(2), Nm(2), pm(2), rho0(2), rhom(2), mum(2)
    real*8, allocatable :: rs(:,:, :), mcblock(:,:)
    integer :: Ntot, N(2), Nswap=0, Nvol
    integer :: nswapacc = 0, nmum(2)
    integer :: nxtrials(2)=0, nxacc(2)=0, nvoltrials=0, nvolacc=0
    integer :: imc, NMC=15000000, jmc, Nequil=1000000, mcblen=1000000, ib, NMCstart
    integer :: logfd=31
    character(LEN=256) :: arg


    real*8 :: rn

    if (command_argument_count() == 3) then
        call get_command_argument(1, arg)
        call load_xyz(arg)
        read(arg, '(7X,I10)') NMCstart
        Nequil = 0
        NMCstart = NMCstart + 1
        
        call get_command_argument(2, arg)
        read (arg, *) kT
        
        call get_command_argument(3, arg)
        read(arg, *) Nswap
        
        write (*,*) NMCstart, kT
    else
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
        call get_command_argument(6, arg)
        read(arg, *) Nswap
        
        Ntot = sum(N)
        allocate(rs(3, Ntot, 2))
        
        V = N/rho0
        Vtot = sum(V)
        bl = V**(1d0/3d0)
        
        if (minval(bl) < 5.0) then
            write (*,'("Box #",I1," is too small. L = ",F8.4, &
               & " < 5*sigma = ",F8.4)') minloc(bl), minval(bl), 5.0
            write (*,*) 'Consider increasing the number of particles or decreasing the density'
        end if
        
        call populate_cube2(rs(:,1:N(1), 1))
        call populate_cube2(rs(:,1:N(2), 2))
        
        NMCstart = Nequil + 1
    end if
    
    write (*,*) N, V, bl


    write(arg, '("qgibbs_kT=",F5.2,".log")') kT
    open(logfd,file=trim(arg), status='REPLACE')
    
    Vstep=0.01*minval(V)
    xstep = 3.0/bl

    allocate(mcblock(16, NMC/mcblen))


    Nvol = 1!minval(N)/20

    beta = 1.0/kT

    rhcore = 2.4

    call vgwinit(Ntot, 'pH2')
    U0 = total_energy(bl)

    do imc=1,Nequil
        call random_number(rn)
        jmc=int(rn*Ntot) + 1
        call mc_move(jmc)
        write (*,*) imc
        call cumulate()
        if (mod(imc,100000) == 0) then
            write (logfd,'("NMC =",I10,", N =",2I5,", V =",2F9.2,", rho =",&
                & 2F9.6,", U =",2F15.6,", p =",2F8.4)') &
                imc, N, V, N/V, U0(1:2), pm/Z
            call dump_xyz(imc)
        end if
    end do

    Z = 0; Um = 0; Nm = 0; pm = 0; Vm = 0; rhom = 0; mum = 0; nmum = 0
    do imc=NMCstart,NMC
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
            write (logfd,'("NMC =",I10,", N =",2I5,", V =",2ES12.5,", rho =",&
                & 2F9.6,", U =",2ES12.5,", p =",2F8.4,", mu =",2F8.4,&
                & ", swapacc = ",F8.5)') &
                imc, N, V, N/V, U0(1:2), pm/Z, mum/nmum, &
                real(nswapacc*(Nswap+Nvol+Ntot))/(Z*real(Ntot))
            call dump_xyz(imc)
        end if
        call cumulate()
        if (mod(imc,mcblen) == 0) then
            ib = imc / mcblen
            mcblock(1, ib) = kT
            mcblock(2:3, ib) = rhom/Z
            mcblock(4:5, ib) = -kT*log(mum/real(nmum))
            mcblock(6:7, ib) = pm/Z
            mcblock(8:9, ib) = Um/Z
            mcblock(10:11, ib) = Nm/Z
            mcblock(12:13, ib) = Vm/Z
            mcblock(14, ib) = real(nswapacc*(Ntot+Nvol+Nswap))/(Z*Ntot)
            mcblock(15:16, ib) = nmum

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
    real*8 :: rso(3), rsn(3), dr(3), p, U0new(4)
    integer :: ibox, j

    ibox = 1
    j = jin
    if (j > N(1)) then
        ibox = 2
        j = j - N(1)
    end if

    nxtrials(ibox) = nxtrials(ibox) + 1


    call random_number(dr)
    rso = rs(:,j,ibox)
    rsn = rso + xstep(ibox) * (dr-0.5)
    rs(:,j,ibox) = rsn - floor(rsn)

    U0new = total_energy(bl, ibox)
    p = exp(-beta * sum(U0new(1:2) - U0(1:2)) )
    call random_number(rn)
    if (p>rn) then
        nxacc(ibox) = nxacc(ibox)  + 1
        U0 = U0new
    else
        rs(:, j, ibox) = rso
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
    real*8 :: pacc, pacc0, rsn(3), rso(3), U0new(4)

    call random_number(rn)
    isrc = 1 + int(2.0*rn)
    idest = 3 - isrc

    if (N(isrc) == 0) return

    call random_number(rsn)
    if (too_close(rsn, rs(:,1:N(idest), idest), bl(idest)) ) return

    call random_number(rn)
    jkill = 1 + int(rn*N(isrc))

    rso = rs(:,jkill, isrc)
    rs(:,jkill,isrc) = rs(:,N(isrc),isrc)
    rs(:,N(idest) + 1,idest) = rsn

    pacc0 = ( V(idest) * real(N(isrc)) ) / ( V(isrc) * real(N(idest) + 1) )

    N(isrc) = N(isrc) - 1
    N(idest) = N(idest) + 1

    U0new = total_energy(bl)

    pacc = pacc0 * exp( - beta * sum(U0new(1:2) - U0(1:2)) )
    mum(idest) = mum(idest) + V(idest)/N(idest)*exp(-beta*(U0new(idest) - U0(idest) ))
    nmum(idest) = nmum(idest) + 1

    call random_number(rn)
    if (pacc > rn) then
        write (*,'(4G15.6,F6.4)') U0new(1:2), U0(1:2), pacc
        nswapacc = nswapacc + 1
        U0 = U0new
    else
        N(isrc) = N(isrc) + 1
        N(idest) = N(idest) - 1
        rs(:,N(isrc),isrc) = rso
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

    U0new = total_energy(bln)
    logp = sum ( -beta * (U0new(1:2) - U0(1:2)) + real(N)*(log(Vn/V)) )
    call random_number(rn)
    if (logp>log(rn)) then
        nvolacc = nvolacc + 1
        V = Vn
        bl = bln
        U0 = U0new
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

function total_energy(bln, ibox) result(utot)
    implicit none
    real*8, intent(in) :: bln(2)
    integer, intent(in), optional :: ibox
    real*8 :: utot(4)
    real*8, allocatable ::  y(:)

    allocate(y(1+9*Ntot))
    
    if (present(ibox)) then
        Utot = U0
        call vgw0(rs(:,1:N(ibox),ibox)*bln(ibox), bln(ibox), beta, &
            0d0, Utot(ibox), y)
    else
        call vgw0(rs(:,1:N(1),1)*bln(1), bln(1), beta, 0d0, Utot(1), y)
        call vgw0(rs(:,1:N(2),2)*bln(2), bln(2), beta, 0d0, Utot(2), y)
    end if
    Utot(3:4) = 0.0
    
    deallocate(y)
end function

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

    write(fname, '("qgibbs_Um_",F5.2,"_NMC.dat")') kT
    open(30,file=trim(fname))
    write(30,'(I10, 16ES16.7)') (i*mcblen, mcblock(:,i), i=1,nmcb)
    close(30)
end subroutine

subroutine dump_xyz(nmcnow)
    implicit none
    integer, intent(in) :: nmcnow
    character(len=256) :: fname
    integer :: ibox

    write(fname, '("qgibbs_",I10,".xyz")') nmcnow
    call replace_char(fname, ' ', '0')
    open(33,file=trim(fname),status='REPLACE')
    do ibox=1,2
        write(33,*) N(ibox)
        write(33,*) bl(ibox)
        write(33,'(("X ", 3F13.8))') rs(:,1:N(ibox),ibox)
    end do
    close(33)
end subroutine

subroutine load_xyz(fname)
    character(*), intent(in) :: fname
    character(16) :: aname
    real*8, allocatable :: rs1(:,:), rs2(:,:)
    integer :: i
    
    open(33, file=trim(fname),action='READ')
    read(33, *) N(1)
    read(33, *) bl(1)
    allocate( rs1( 3, N(1) ) )
    read(33, *) (aname, rs1(:,i), i=1,N(1))
    
    read(33, *) N(2)
    read(33, *) bl(2)
    allocate( rs2( 3, N(2) ) )
    read(33, *) (aname, rs2(:,i), i=1,N(2))
    close(33)
    
    Ntot = sum(N)
    V = bl**3
    Vtot = sum(V)
    
    allocate(rs(3,Ntot,2))
    rs(:,1:N(1),1) = rs1
    rs(:,1:N(2),2) = rs2
    deallocate(rs1, rs2)
end subroutine load_xyz
end program qgibbs
