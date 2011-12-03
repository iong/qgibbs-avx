    function total_energy(N, rs, bl, deBoer) result(utot)
        use vgw
        implicit none
        double precision, intent(in) :: rs(3,N), bl, deBoer
        double precision :: utot

        call vgwinit(N, 'LJ', massx=1d0/deBoer**2)

        call vgw0(rs, bl, beta, Utot, scale=.TRUE.)

        call vgwcleanup()

    end function total_energy
