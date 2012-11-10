    vgw_atol = 1d-4
    dt0 = 0d0
    dtmax = 0d0
    dtmin = 0d0


    Ulrc = 0d0
    UXXlrc = 0d0

    if (species=='pH2-4g') then
        NGAUSS=4
        LJA(1:4) = (/ 1.038252215127D0, 0.5974039109464D0, 0.196476572277834D0, &
                    0.06668611771781D0 /)
        LJC(1:4) = (/ 96609.488289873d0, 14584.62075507514d0, -365.460614956589d0, &
                    -19.5534697800036d0 /)
        mass = 2.0*0.020614788876D0
        rc = 8.0
        dt0=5d-4
        dtmax = 2d-3
        dtmin = 1d-5
        vgw_atol = (/ 1d-3, 1d-4, 1d0 /)
    else if (species == 'LJ') then
        NGAUSS = 3
        LJA(1:3) = (/ 6.65, 0.79, 2.6 /)
        LJC(1:3) = (/ 1840d0, -1.48d0, -23.2d0 /)
        mass = 1.0
        rc = 2.75
        rfullmatsq = 1.8d0**2
        vgw_atol = (/ 1d-5, 1d-7, .1d0 /)
        dtmax = .25d0
        dt0=1d-5
        dtmin=1d-7
    else if (species == 'TT:Ne-Ne') then
        NGAUSS = 3
        LJA(1:3) = (/ 10.2645, 2.25997, 0.611381 /)
        LJC(1:3) = (/ 7065.59, -10.6659, -0.239751 /)
        mass = 1.0
        rc = 2.75
        rfullmatsq = 1.8d0**2
        vgw_atol = (/ 1d-5, 1d-7, .1d0 /)
        dtmax = .25d0
        dt0=1d-5
        dtmin=1d-7
    end if

    if (present(M)) then
        mass = M
    end if
    if (present(rcutoff)) then
        rc = rcutoff
    end if


    if (species == 'LJ') then
        Ulrc = 8d0*M_PI/9d0*(1d0/rc**9 - 3d0/rc**3)
        UXXlrc = 32d0*M_PI*(2d0/rc**11 - 1d0/rc**5)
    else if (species == 'TT:Ne-Ne') then
        C6 = 1.2007
        C8 = 0.4983
        C10 = 0.2484

        Ulrc = -4d0*M_PI*(C10/(7.0*rc**7) + C8/(5.0*rc**5) + C6/(3.0*rc**3))
        UXXlrc = -8.0*M_PI*(5.0 * C10 + 4.0*C8*rc**2 + 3.0*C6*rc**4) / &
            (3.0 * rc**9)
    end if
