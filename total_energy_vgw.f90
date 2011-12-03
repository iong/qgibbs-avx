double precision function total_energy(rs, bl, beta, deBoer)
    use vgw
    implicit none
    double precision, intent(in) :: rs(:,:), bl, beta, deBoer
    integer N

    N = size(rs, 2)

    call vgwinit(N, 'LJ', massx=1d0/deBoer**2)

    call vgw0(rs*bl, bl, beta, total_energy)

    call vgwcleanup()

end function total_energy
