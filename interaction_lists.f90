subroutine interaction_lists(x, y, z)
    implicit none
    real*8, intent(in) :: x(:), y(:), z(:)
    integer :: I,J, NN, N
    real*8 :: rsq(size(x)),rc2

    N = size(x)
    rc2=rc**2
!$omp parallel
!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(I, J, RSQ,NN)
    do I=1,N-1
        NN = 0
        rsq(i+1:) = min_image(x(i) - x(i+1:), bl)**2 &
                  + min_image(y(i) - y(i+1:), bl)**2 &
                  + min_image(z(i) - z(i+1:), bl)**2
 
        do J=I+1,N
            if(rsq(j) <= rc2) then
                NN = NN + 1
                NBIDX(NN, I) = J
            endif
        enddo
        NNB(i) = NN
    enddo
!$OMP ENDDO
!$omp end parallel

    NNB(N) = 0
    nnbmax = maxval(nnb)
    nnbmax = (nnbmax/8 + 1) * 8
end subroutine interaction_lists
