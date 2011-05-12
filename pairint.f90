module pairint
    use utils
    implicit none
    real*8, allocatable :: y6cache(:)

   type,abstract:: pair_interactions
    contains
        procedure(UXXp),deferred :: Up
        procedure(pi_lrcorr),deferred :: lrcorr
    end type pair_interactions

    abstract interface
        subroutine UXXp(this, rsj, rs, blj, U, V)
            import pair_interactions
            class(pair_interactions), intent(in) :: this
            real*8,intent(in) :: rsj(3), rs(:,:), blj
            real*8, intent(out) :: U(:)
            real*8, intent(out), optional :: V(:)
        end subroutine
        real*8 function pi_lrcorr(this, N, V, rc)
            import pair_interactions
            class(pair_interactions), intent(in) :: this
            real*8,intent(in) :: V, rc
            integer,intent(in) :: N
        end function
    end interface

    type,extends(pair_interactions) :: LennardJones
    contains
        procedure:: Up => ULJp
        procedure:: lrcorr => lj_lrcorr
    end type LennardJones

    type, extends(pair_interactions) :: SilveraGoldman
        real*8 :: A = 14.376, B = 2.9618, G = 0.035471, &
            C6 = 84106d0, C8 = 4.1737d5, C9 = 1.4685d5, C10 = 2.6137d6, &
            rcrit = 4.4021
    contains
        procedure :: Up => USGp
        procedure:: lrcorr => sg_lrcorr
    end type SilveraGoldman

contains
subroutine ULJp(this, rsj, rs, blj, U, V)
    implicit none
    class(LennardJones), intent(in) :: this
    real*8, intent(in) :: rsj(3), rs(:,:), blj
    real*8, intent(out) :: U(:)
    real*8, intent(out), optional :: V(:)
    real*8 :: drs(3, size(rs, 2)), drs0(3)
    integer :: Nl, i

    Nl = size(rs, 2)

    do i=1,Nl
        drs0 = rsj - rs(:,i)
        drs(:,i) = min_image(drs0)
    enddo
    y6cache(1:Nl) = 1.0 / sum(drs**2, 1)**3
    forall (i=1:Nl, y6cache(i) < 64.0)
        y6cache(i) = 0.0
    end forall
    y6cache(1:Nl) = 1.0/(blj**6) * y6cache(1:Nl)
    U = 4d0*(y6cache(1:Nl)**2 - y6cache(1:Nl))
    if (present(V)) then
        V(1:Nl) = 12.0*U(1:Nl) + 24.0*y6cache(1:Nl)
    end if
end subroutine

real*8 function lj_lrcorr(this, N, V, rc)
    class(LennardJones), intent(in) :: this
    real*8,intent(in) :: V, rc
    integer,intent(in) :: N

    lj_lrcorr = real(N**2)/V*8.0/9.0*M_PI*(1.0/rc**9 - 3.0/rc**3)
end function

subroutine USGp(this, rsj, rs, blj, U, V)
    implicit none
    class(SilveraGoldman), intent(in) :: this
    real*8,intent(in) :: rsj(3), rs(:,:), blj
    real*8, intent(out) :: U(:)
    real*8, intent(out), optional :: V(:)
    real*8 :: blsq, drs(3), drsq
    real*8, dimension(size(rs,2)) :: VLJ, r1, r2, r6, r8, Ul
    integer :: rcutidx(size(rs,2))
    integer :: N, i, j, nrcut


    N = size(rs, 2)
    blsq = blj**2
    nrcut = 0
    do i=1,N
        drs = rsj - rs(:,i)
        drsq = sum(min_image(drs)**2)
        if (drsq < 0.25) then
            nrcut = nrcut + 1
            r2(nrcut) = drsq
            rcutidx(nrcut) = i
        end if
    end do

    r2(1:nrcut) = r2(1:nrcut) * blj**2
    r1(1:nrcut) = sqrt(r2(1:nrcut))
    r6(1:nrcut) = r2(1:nrcut)**3
    r8(1:nrcut) = r6(1:nrcut)*r2(1:nrcut)

    VLJ(1:nrcut) = this%C6/r6(1:nrcut) + this%C8/r8(1:nrcut) &
            - this%C9/(r8(1:nrcut) * r1(1:nrcut) ) + this%C10/(r8(1:nrcut) * r2(1:nrcut))
    !firi = 6.0*C6 / r6 + 8.0*C8 / r8 - 9.0*C9 / (r8*r1) + 10.0*C10 / (r8*r2)

    Ul(1:nrcut) = exp(this%A - this%B*r1(1:nrcut) - this%G*r2(1:nrcut))! -VLJ*fc!-SGcutoff
    do i=1,nrcut
        if (r1(i) <= this%rcrit) then
            Ul(i) = Ul(i) - VLJ(i)*exp( - (this%rcrit / r1(i) - 1.0)**2 )
        else
            Ul(i) = Ul(i) - VLJ(i)
        endif
    end do
    U = 0.0
    U(rcutidx(1:nrcut)) = Ul(1:nrcut)
    if (present(V)) then
        V=0.0
    end if
end subroutine

real*8 function sg_lrcorr(this, N, V, rc)
    class(SilveraGoldman), intent(in) :: this
    real*8,intent(in) :: V, rc
    integer,intent(in) :: N

    sg_lrcorr = 2.0*M_PI*N**2/V * ( this%C6/(3.0 * rc**3) + this%C8/(5.0 * rc**5) &
            - this%C9/(6.0 * rc**6) + this%C10/(7.0 * rc**7) )
end function
end module
