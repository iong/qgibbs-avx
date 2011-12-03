double precision function total_energy(rs, bl, beta, deBoer)
    use vgwspfm
    implicit none 
    double precision, intent(in) :: rs(:,:), bl, beta, deBoer
    integer N

    N = size(rs, 2)
    write (*,*) 'N =', N

    call vgwspfminit(N, 'LJ', massx=1d0/deBoer**2)

    call vgw0spfm(rs*bl, bl, beta, total_energy)

    call vgwspfmcleanup()

end function total_energy
