SUBROUTINE vgw0fm(Q0, BL_, TAUMAX,TAUI, W)
    use omp_lib
    IMPLICIT NONE
    REAL*8, intent(in) :: Q0(:,:), TAUMAX, TAUI, BL_
    REAL*8, intent(out) :: W
    real*8 :: GDET
    real*8 :: DT, next_stop
    integer :: j, info
    logical :: mm = .FALSE.


    Natom = size(Q0, 2)
    BL = BL_

    if (TAUI <= 0.0d0) then
        T = TAUMIN
    else
        T = 0.5*TAUI
    endif

    if (TAUI <= 0.0d0) then
        call init_gaussians(Q0, TAUMIN)
    end if
    next_stop = 0.5d0*TAUMAX
    do
        DT = 1d-5!1d-2*sqrt(T)
        if (T+DT > next_stop) then
            DT = next_stop - T
            T = next_stop
        else
            T = T + DT
        end if
        call RHSSfm(DT, mm)
        if (T == next_stop) exit
    end do

    GU = G
    call dpotrf('U', 3*Natom, GU, 3*Natom, info)
    GDET=1.0
    DO j=1,3*Natom
        GDET = GDET * GU(j,j)
    ENDDO
    GDET = GDET**2

    W=-(1/TAUMAX)*(2.0*gama - 0.5*log(GDET))! - 3.0*Natom*log(2.0*sqrt(M_PI)))
    !write (*,*) gama
END SUBROUTINE