    function total_energy(N, rs, bl, deBoer) result(utot)
        use vgwspfm
        implicit none
        double precision, intent(in) :: rs(3,N), bl, deBoer
        double precision :: utot

        call vgwspfminit(N, 'LJ', massx=1d0/deBoer**2)

        call vgw0spfm(rs, bl, beta, Utot, scale=.TRUE.)

        call vgwspfmcleanup()

    end function total_energy
