program qgibbs
    use utils, only: min_image, replace_char, M_PI
    use vgw
    implicit none
    real*8 :: rhcore
    real*8 :: Vtot, V(2), bl(2), kT, beta, U0(4), xstep(2), Vstep
    double precision :: VdeBroglie, deBoer
    real*8 :: Z, Um(2), Vm(2), Nm(2), pm(2), rho0(2), rhom(2), mum(2)
    real*8, allocatable :: rs(:,:, :), mcblock(:,:)
    integer :: Ntot, N(2), Nswap=0, Nvol
    integer :: nswapacc = 0, nmum(2)
    integer :: nxtrials(2)=0, nxaccsum(2) = 0, nxacc(2)=0, nvoltrials=0, nvolacc=0
    integer :: imc=1, NMC=15000000, jmc, Nequil=100000, mcblen=10000, ib, NMCstart
    integer :: logfd=31, swapfd=32, argshift,histfd=35
    character(LEN=256) :: arg, datadir, logfile, method='vgw', histfile
    logical :: restart = .FALSE., oldlog
    real*8 :: rn
    namelist/input_parameters/N, rho0, kT, deBoer, Nswap, method
    namelist/restart_parameters/imc, N, Vtot, V, bl, nxtrials, nxacc, xstep, nswapacc, nvoltrials, nvolacc, Vstep, rs

    interface
        double precision function energy(rs, bl)
            double precision, intent(in) :: rs(:,:), bl
        end function
    end interface
    procedure(energy), pointer :: total_energy


    if (command_argument_count() == 2) then
        restart = .TRUE.
        call get_command_argument(1, arg)
        open(33,file=trim(arg))
        read(33,NML=input_parameters)
        close(33)

        Ntot = sum(N)
        allocate(rs(3, Ntot, 2))

        call get_command_argument(2, arg)
        open(33,file=trim(arg))
        read(33,NML=restart_parameters)
        close(33)

        imc = imc + 1
    else
        argshift = 0
        call get_command_argument(1, arg)
        if (arg == '-m') then
            call get_command_argument(1, method)
            argshift = 2
        endif
        call get_command_argument(argshift+1, arg)
        read (arg, *) N(1)
        call get_command_argument(argshift+2, arg)
        read (arg, *) N(2)
        call get_command_argument(argshift+3, arg)
        read (arg, *) rho0(1)
        call get_command_argument(argshift+4, arg)
        read (arg, *) rho0(2)
        call get_command_argument(argshift+5, arg)
        read (arg, *) kT
        call get_command_argument(argshift+6, arg)
        read (arg, *) deBoer
        call get_command_argument(argshift+7, arg)
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

    total_energy =>total_energy_vgw
    !if (method == 'vgwspfm') then
    !    total_energy => total_energy_vgwspfm
    !end if


    write(datadir, '("qgibbs_kT=",F6.3)') kT
    call replace_char(datadir, ' ', '0')
    call system('mkdir -p '//trim(datadir))

    open(39,file=trim(datadir)//"/input.h")
    write(39,NML=input_parameters)
    close(39)

    logfile = trim(datadir)// '/block_average.dat'
    histfile = trim(datadir)// '/histogram_data.bin'
    if (restart) then
        inquire(file=trim(logfile), EXIST=oldlog)
        if (oldlog) then
            open(logfd,file=trim(logfile), status='OLD', position='APPEND')
            open(histfd,file=trim(histfile), status='OLD', position='APPEND', &
                form='UNFORMATTED', access='stream')
        else
            open(logfd,file=trim(logfile), status='NEW')
            call dump_block_avg(just_header=.TRUE.)

            open(histfd,file=trim(histfile), status='NEW', form='UNFORMATTED', &
                access='stream')
            call dump_histogram_header()
        end if
    else
        open(logfd,file=trim(logfile), status='REPLACE')
        call dump_block_avg(just_header=.TRUE.)

        open(histfd,file=trim(histfile), status='REPLACE', form='UNFORMATTED', &
            access='stream')
        call dump_histogram_header()
        Vstep=0.01*minval(V)
        xstep = 3.0/bl
    end if

    open(swapfd, file=trim(datadir)// '/swap.dat')

    Nvol = 1!minval(N)/20

    beta = 1.0/kT

    rhcore = 0.8 ! exp(-VLJ/kT) < 1e-9 for kT=1.2
    VdeBroglie = (2*M_PI*deBoer/sqrt(kT))**3


    U0(1) = total_energy(rs(:,:N(1), 1), bl(1))
    U0(2) = total_energy(rs(:,:N(2), 2), bl(2))

    call reset_averages()
    do
        if (imc > Nequil) exit

        call random_number(rn)
        jmc=int(rn*Ntot) + 1
        call mc_move(jmc)
        call cumulate()

        if (mod(imc,mcblen) == 0) then
            call dump_block_avg()
            call checkpoint()
        end if
        if (mod(imc,10*mcblen) == 0) then
            call dump_xyz()
        end if
        imc = imc + 1
    end do

    call reset_averages()
    do
        if (imc > NMC) exit

        call random_number(rn)
        jmc=1+int(rn*(Ntot + Nswap + Nvol))
        if (jmc<= Ntot) then
            call mc_move(jmc)
        elseif (jmc <= Ntot + Nswap) then
            call mc_swap()
        elseif (jmc <= Ntot + Nswap + Nvol) then
            call mc_vol()
        endif

        call cumulate()

        write (histfd) N(1),real(V(1),4), real(U0,4)

        if (mod(imc,mcblen) == 0) then
            call dump_block_avg()
            call checkpoint()
        end if
        if (mod(imc,10*mcblen) == 0) then
            call dump_xyz()
        end if

        imc = imc + 1
    enddo

    close(logfd)
    close(swapfd)
    close(histfd)

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
    end subroutine populate_cube2

    function too_close(rsj, ibox, jskip)
        real*8, intent(in) :: rsj(3)
        integer, intent(in) :: ibox, jskip
        logical :: too_close
        real*8 :: drsq, rshc_sq, drs(3)
        integer :: i

        too_close = .FALSE.
        rshc_sq = (rhcore / bl(ibox))**2
        do i=1,N(ibox)
            if ( i == jskip ) cycle

            drs = rsj - rs(:,i, ibox)
            drsq = sum(min_image(drs)**2)

            if (drsq < rshc_sq) then
                too_close = .TRUE.
                return
            end if
        enddo
    end function too_close


    subroutine mc_move(jin)
        implicit none
        integer, intent(in) :: jin
        integer, parameter :: nstepcheck=200
        real*8 :: rso(3), rsn(3), dr(3), p, U0new
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
        rsn = rso + xstep(ibox) * (dr-0.5d0)
        rsn = rsn - floor(rsn)

        if (too_close(rsn, ibox, j)) then
            p = -1.0
        else
            rs(:,j,ibox) = rsn 

            U0new = total_energy(rs(:,:N(ibox), ibox), bl(ibox))

            p = exp(-beta * (U0new - U0(ibox)) )
        end if

        call random_number(rn)
        if (p>rn) then
            nxacc(ibox) = nxacc(ibox)  + 1
            nxaccsum(ibox) = nxaccsum(ibox) + 1
            U0(ibox) = U0new
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
    end subroutine mc_move

    subroutine mc_swap()
        implicit none
        integer :: isrc, idest, jkill
        real*8 :: pacc, pacc0, rsn(3), rso(3), U0new(4)

        call random_number(rn)
        isrc = 1 + int(2.0*rn)
        idest = 3 - isrc

        if (N(isrc) == 0) return

        call random_number(rsn)
        if (too_close(rsn, idest, -1) ) return

        call random_number(rn)
        jkill = 1 + int(rn*N(isrc))

        rso = rs(:,jkill, isrc)
        rs(:,jkill,isrc) = rs(:,N(isrc),isrc)
        rs(:,N(idest) + 1,idest) = rsn

        pacc0 = ( V(idest) * real(N(isrc)) ) / ( V(isrc) * real(N(idest) + 1) )

        N(isrc) = N(isrc) - 1
        N(idest) = N(idest) + 1

        U0new(1) = total_energy(rs(:,:N(1), 1), bl(1))
        U0new(2) = total_energy(rs(:,:N(2), 2), bl(2))

        pacc = pacc0 * exp( - beta * sum(U0new(1:2) - U0(1:2)) )
        mum(idest) = mum(idest) + V(idest)/(VdeBroglie * N(idest)) * &
             exp(-beta*( U0new(idest) - U0(idest) ) )
        nmum(idest) = nmum(idest) + 1

        call random_number(rn)
        if (pacc > rn) then
            !write (*,'(4G15.6,F6.4)') U0new(1:2), U0(1:2), pacc
            nswapacc = nswapacc + 1
            U0 = U0new
            write (swapfd, '(3I10)') imc, N
        else
            N(isrc) = N(isrc) + 1
            N(idest) = N(idest) - 1
            rs(:,N(isrc),isrc) = rs(:,jkill,isrc)
            rs(:,jkill,isrc) = rso
        end if
    end subroutine mc_swap


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

        U0new(1) = total_energy(rs(:,:N(1), 1), bln(1))
        U0new(2) = total_energy(rs(:,:N(2), 2), bln(2))

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
    end subroutine mc_vol



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
    end subroutine cumulate

    subroutine reset_averages()
        rhom = 0
        mum = 0
        nmum = 0
        pm = 0
        Um = 0
        Nm = 0
        Vm = 0
        nswapacc = 0
        nxaccsum = 0
        Z = 0
    end subroutine reset_averages

    subroutine dump_block_avg(just_header)
        implicit none
        logical, optional :: just_header
        integer :: i
        character(256) :: fname
        double precision :: mu(2) = 0d0

        if (present(just_header)) then
            if (just_header) then
                write(logfd, '("# NMC kT rho1 rho2 mu1 mu2 p1 p2 U1 U2 N1 N2 V1 V2 swapacc1 swapacc2 mum1 mum2 nmum1 nmum2")')
                flush(logfd)
                return
            end if
        end if

        if (product(nmum) /= 0 .and. product(mum) /= 0) then
            mu = -kT*log(mum/real(nmum))
        end if

        write(logfd,'(I10, 19(" ", G18.6))') &
             imc, kT, rhom/Z, mu, pm/Z, Um/Z, Nm/Z, Vm/Z,&
             real(nswapacc)/real(nxaccsum), mum, nmum

        flush(logfd)

        call reset_averages()
    end subroutine dump_block_avg

    subroutine dump_histogram_header()
        integer :: i, j
        integer :: header_length = 256

        write(histfd) &
            "###This header is ", header_length, "bytes long."//new_line("x")//&
            "### Start reading the binary data from bye 257." //new_line('x')//&
            "[('N1', 'i4'), ('V1', 'f4'), ('E1', 'f4'), ('E2', 'f4')]"//&
            new_line('x')

        inquire(histfd, pos=i)
        write(histfd) (/ ('#', j=i,header_length - 1) /)
        write(histfd) new_line('x')
    end subroutine

    subroutine dump_xyz()
        implicit none
        character(len=256) :: fname
        integer :: ibox

        write(fname, '(I10,".xyz")') imc
        call replace_char(fname, ' ', '0')
        open(33,file=trim(datadir)//'/'//trim(fname),status='REPLACE')
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

    subroutine checkpoint()
        open(33,file=trim(datadir)//'/checkpoint.h',status='REPLACE')
        write(33,NML=restart_parameters)
        close(33)
    end subroutine checkpoint


    double precision function total_energy_vgw(rs, bl)
        implicit none
        double precision, intent(in) :: rs(:,:), bl

        call vgwinit('LJ')
        total_energy_vgw = vgw0(rs*bl, bl, beta, deBoer)
    end function
end program qgibbs
