subroutine interaction_lists(r)
    implicit none
    real*8, intent(in) :: r(3,Natom)
    integer :: I,J, NN
    real*8 :: rsq,rcsq, dr(3), bl2
    

    rcsq=rc**2
    bl2=0.5*bl
!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(J, RSQ, dr, NN)
    do I=1,Natom-1
        NN = 0
        do J=I+1,Natom
            dr = r(:,j) - r(:,i)

            where (abs(dr) > bl2)
                dr  = dr - sign(bl, dr)
            end where

            rsq = sum(dr**2)
 
            if(rsq <= rcsq) then
                NN = NN + 1
                NBIDX(NN, I) = J
            endif
        enddo
        NNB(i) = NN
    enddo
!$OMP ENDDO

!$omp master
    NNB(Natom) = 0
    nnbmax = maxval(nnb)
    nnbmax = (nnbmax/8 + 1) * 8
!$omp end master
end subroutine interaction_lists
