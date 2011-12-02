double precision function total_energy(rs, bl, beta, deBoer)
    use vgwspfm
    implicit none 
    double precision, intent(in) :: rs(:,:), bl, beta, deBoer
    integer N

    N = size(rs, 2)

    call vgwspfminit(N, 'LJ', massx=1d0/deBoer**2)

    call vgw0spfm(rs, bl, beta, total_energy, scale2bl=.TRUE.)

    call vgwspfmcleanup()

end function total_energy
