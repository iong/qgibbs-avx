! GDATA is a 4D matrix:
!     1st idx: NGAUSS
!     2nd idx: LJC
!     3rd idx: LJA
!     4th irdx: intra vs inter interaction

!MOLDATA
! 1st: molecule number
! 2nd: molecule type
! 3rd: offset
! 4th: number of atoms in the molecule (NAM)

!NAM number of atoms in the molecule


MODULE VGWM
    IMPLICIT NONE
    private
    INTEGER, SAVE :: N_ATOM,NEQ,MCNT,MSIZE,N_MOL,LRGMSZE,NPAIRS,N3ATOM,IARRAY(7)
    INTEGER, ALLOCATABLE, SAVE :: MOLDATA(:,:),IWORK(:)
    INTEGER, PARAMETER :: NUMOBSERV=3,MAXGAUSS=20,MAXTYPE=20,ISIZE=5
    DOUBLE PRECISION, SAVE :: ATOL,T,TAUI,ENOT,RCUT,RCUTSQ,BL,BL2
    DOUBLE PRECISION, PARAMETER :: ATOMICMASS=0.020614788876D0,PI=3.1415926D0,RTOL=0.0D0
    DOUBLE PRECISION, PARAMETER :: ANGTOCM=1.0D-24
    DOUBLE PRECISION, ALLOCATABLE, SAVE :: MASSARRY(:),GMAT(:),Y(:),RWORK(:),PARRAY(:,:)
    DOUBLE PRECISION, SAVE :: GDATA(2*MAXGAUSS+1,MAXTYPE,MAXTYPE,2),E_ZERO
    CHARACTER*100 :: CNST
    LOGICAL :: PBC,UPDATED,EXTENDED
    EXTERNAL DLSODE
    public :: INITIALIZE_VGWM, VGWMQUENCH, CLEANUP_VGWM
contains

    SUBROUTINE CLEANUP_VGWM()
        DEALLOCATE(MASSARRY,MOLDATA,PARRAY,GMAT,Y,IWORK,RWORK)
    END SUBROUTINE CLEANUP_VGWM

    SUBROUTINE JAC()
    END SUBROUTINE JAC

    SUBROUTINE INITIALIZE_VGWM(N,TI,AT,BLX,RCT, vgwdata)
        IMPLICIT NONE
        character(*), optional, intent(in) :: vgwdata
        INTEGER :: I,J,K,N,NSETS,ITYPE,MSIZE1,MSIZE2,NG,CP,CNT,ACNT,NAM,IND,LSZE,PCNT
        DOUBLE PRECISION :: BLX,COE,AT,ALPHA,TI,RCT
        CHARACTER*100 :: DUMMY

        N_ATOM=N
        N3ATOM=3*N_ATOM
        ATOL=AT
        TAUI=TI!MIN(TI,(0.5D0/TMIN)-TI)
        RCUT=RCT                                                  ! Cut off radius   
        RCUTSQ=RCUT**2       
        NPAIRS=0.5D0*N_ATOM*(N_ATOM-1)                            ! Maximum number of pairs
        PCNT=0                                                    ! Pair count

        IF(BLX.GT.0) THEN    ! Periodic boundary conditions ?
            PBC=.TRUE.
            BL=BLX
            BL2=BL/2.0D0
        ELSE
            PBC=.FALSE.
            BL=-1.0D0
        ENDIF

        ALLOCATE(MASSARRY(N3ATOM),MOLDATA(4,N_ATOM),PARRAY(2,NPAIRS))   ! Data matrix for observables, see VGWFUPDATE

        DO I=1,N_ATOM-1 ! Inititalize pair array
            DO J=I+1,N_ATOM
                PCNT=PCNT+1
                PARRAY(1,PCNT)=I
                PARRAY(2,PCNT)=J
            ENDDO
        ENDDO

        if (present(vgwdata)) then
            OPEN(UNIT=7,FILE=trim(vgwdata),STATUS='OLD')
        else
            OPEN(UNIT=7,FILE='vgwdata',STATUS='OLD') ! Read molecule data (mass, molecule, atom type)
        end if
        DO I=1,N_ATOM
            READ(7,*) MASSARRY(3*(I-1)+1),MOLDATA(1,I),MOLDATA(2,I),DUMMY
            MASSARRY(3*(I-1)+1)=1/(ATOMICMASS*MASSARRY(3*(I-1)+1))
            MASSARRY(3*(I-1)+2)=MASSARRY(3*(I-1)+1)
            MASSARRY(3*(I-1)+3)=MASSARRY(3*(I-1)+1)
        ENDDO
        CLOSE(7)

        OPEN(UNIT=7,FILE='gauss.dat',STATUS='OLD')

        READ(7,*) NSETS

        DO I=1,NSETS                                     ! Read Gaussian parameters from file into data array
            READ(7,*) ITYPE,MSIZE1,MSIZE2,NG              ! NG = number of gaussians
            GDATA(1,MSIZE1,MSIZE2,ITYPE)=NG
            GDATA(1,MSIZE2,MSIZE1,ITYPE)=NG
            DO J=1,NG
                READ(7,*) COE, ALPHA
                GDATA(2*J,MSIZE1,MSIZE2,ITYPE)=COE
                GDATA(2*J+1,MSIZE1,MSIZE2,ITYPE)=ALPHA
                IF(MSIZE1.NE.MSIZE2) THEN
                    GDATA(2*J,MSIZE2,MSIZE1,ITYPE)=COE
                    GDATA(2*J+1,MSIZE2,MSIZE1,ITYPE)=ALPHA
                ENDIF
            ENDDO
        ENDDO

        READ(7,*) CNST
        READ(7,*) E_ZERO

        CLOSE(7)

        NAM=0
        NEQ=0
        IND=1
        LRGMSZE=1
        N_MOL=0

        MSIZE1=MOLDATA(1,1)

        DO I=1,N_ATOM                                    ! Determine size of G matrix
            MSIZE2=MOLDATA(1,I)
            IF(MSIZE1.EQ.MSIZE2) THEN
                NAM=NAM+1
            ELSE
                N_MOL=N_MOL+1
                NEQ=NEQ+3*NAM*(3*NAM+1)/2
                DO J=1,NAM
                    MOLDATA(1,IND)=N_MOL                    ! Molecule number
                    MOLDATA(3,IND)=3*(J-1)                  ! Block offset for atom in molecule
                    MOLDATA(4,IND)=NAM                      ! Molecule Size
                    IF(NAM.GT.LRGMSZE) LRGMSZE=NAM          ! Largest molecule?
                    IND=IND+1
                ENDDO
                NAM=1
            ENDIF
            MSIZE1=MSIZE2
        ENDDO

        NEQ=NEQ+3*NAM*(3*NAM+1)/2
        N_MOL=N_MOL+1

        DO J=1,NAM
            MOLDATA(1,IND)=N_MOL
            MOLDATA(3,IND)=3*(J-1)
            MOLDATA(4,IND)=NAM
            IF(NAM.GT.LRGMSZE) LRGMSZE=NAM
            IND=IND+1
        ENDDO

        ALLOCATE(GMAT(NEQ))                               ! G matrix
        NEQ=NEQ+N3ATOM+1                                ! Number of equations plus gamma plus coords

        IARRAY(5)=20+16*NEQ                                     ! LRW variable (dlsode rwork array size)
        IARRAY(6)=20                                            ! LIW variable (dlsode iwork array size)

        ALLOCATE(Y(NEQ),RWORK(IARRAY(5)),IWORK(IARRAY(6)))

        CNT=1
        IND=1
        ACNT=1

        DO I=1,N_MOL                                 ! Initialize G matrix
            NAM=MOLDATA(4,IND)
            IND=IND+NAM
            DO J=1,3*NAM
                GMAT(CNT)=MASSARRY(ACNT)
                ACNT=ACNT+1
                CNT=CNT+1
                DO K=J+1,3*NAM
                    GMAT(CNT)=0.0D0
                    CNT=CNT+1
                ENDDO
            ENDDO
        ENDDO

        LSZE=LRGMSZE                       ! Return largest molecule size

    END SUBROUTINE INITIALIZE_VGWM

    SUBROUTINE VGWMQUENCH(Q,BETAMAX,EFFPOT)
        IMPLICIT NONE
        DOUBLE PRECISION, intent(in) :: Q(3,N_atom), BETAMAX
        DOUBLE PRECISION, intent(out) :: EFFPOT
        DOUBLE PRECISION :: QT(N3ATOM),UGAUSS
        INTEGER :: I

        Y=0.0D0

        CALL INITDLSODE

        IF(PBC) THEN
            QT=BL*reshape(Q, (/ N3atom /) )
        ELSE
            QT=reshape(Q, (/ N3atom /) )
        ENDIF

        IF(RCUT.GT.0) CALL PAIRS(QT)

        DO I=1,N3ATOM
            Y(I+1)=QT(I)                           ! COPY INITITAL COORDS TO Y VECTOR
        ENDDO

        DO I=1,NEQ-N3ATOM-1
            Y(N3ATOM+1+I)=GMAT(I)*TAUI
        ENDDO
        
        CALL GAUSSENERGY(Y(2),UGAUSS)
        Y(1)=-TAUI*UGAUSS
        
        

        T=TAUI

        CALL DLSODE(RHSM,NEQ,Y,T,0.5D0*BETAMAX,IARRAY(1),RTOL,ATOL,IARRAY(2),IARRAY(4),IARRAY(3),RWORK,IARRAY(5),IWORK,IARRAY(6),JAC,IARRAY(7))

        CALL LNPS(EFFPOT)

        EFFPOT=-EFFPOT/BETAMAX
        write (*,*) Y(1)
    END SUBROUTINE VGWMQUENCH

    SUBROUTINE PAIRS(Q)

        IMPLICIT NONE
        INTEGER :: I,J,K
        DOUBLE PRECISION :: Q(3,N_ATOM),RSQ,QIJ(3)

        NPAIRS=0

        DO I=1,N_ATOM-1
            DO J=I+1,N_ATOM
                RSQ=0.0D0
                DO K=1,3
                    QIJ(K)=Q(K,I)-Q(K,J)
                    IF(PBC) THEN                                            ! Periodic boundary conditions
                        IF(QIJ(K).GT.BL2) QIJ(K)=QIJ(K)-BL
                        IF(QIJ(K).LT.-BL2) QIJ(K)=QIJ(K)+BL
                    ENDIF
                    RSQ=RSQ+QIJ(K)**2
                ENDDO
                IF(RSQ.LE.RCUTSQ) THEN
                    NPAIRS=NPAIRS+1
                    PARRAY(1,NPAIRS)=I
                    PARRAY(2,NPAIRS)=J
                ENDIF
            ENDDO
        ENDDO

    END SUBROUTINE PAIRS

    SUBROUTINE INITDLSODE

        IMPLICIT NONE

        IWORK=0.0D0
        RWORK=0.0D0

        IARRAY(1)=1                    ! ITOL
        IARRAY(2)=1                    ! ITASK
        IARRAY(3)=1                    ! IOPT
        IARRAY(4)=1                    ! ISTATE
        IARRAY(5)=20+16*NEQ            ! LRW
        IARRAY(6)=20                   ! LIW
        IARRAY(7)=10                   ! MF

        IWORK(6)=100000
        IWORK(5)=4
        IWORK(7)=0
        IWORK(8)=0
        IWORK(9)=0
        IWORK(10)=0
        RWORK(5)=0.0D0
        RWORK(6)=0.0D0
        RWORK(7)=0.0D0
        RWORK(8)=0.0D0
        RWORK(9)=0.0D0
        RWORK(10)=0.0D0

    END SUBROUTINE INITDLSODE

    SUBROUTINE LNPS(LOGZ)

        IMPLICIT NONE
        INTEGER :: I,J,K,CNT,DIM,INFO,IND
        DOUBLE PRECISION :: C(3*LRGMSZE,3*LRGMSZE),LOGZ,GAMMA,DET
        DOUBLE PRECISION :: DETL
        CHARACTER*1 :: UPLO="U"

        GAMMA=Y(1)
        CNT=2+N3ATOM
        C=0.0D0
        IND=1
        DET=0.0D0

        DO I=1,N_MOL
            DIM=3*MOLDATA(4,IND)
            IND=IND+MOLDATA(4,IND)
            DO J=1,DIM
                DO K=J,DIM
                    C(J,K)=Y(CNT)
                    CNT=CNT+1
                ENDDO
            ENDDO

            DETL=1.0D0
            CALL DPOTRF(UPLO,DIM,C,DIM,INFO)

            DO J=1,DIM
                DETL=DETL*C(J,J)
            ENDDO

            DET=DET+2.0D0*LOG(DETL)

        ENDDO

        LOGZ=2.0D0*GAMMA-0.5D0*DET

    END SUBROUTINE LNPS

    SUBROUTINE INVDET(A,M,DET)
        IMPLICIT NONE
        DOUBLE PRECISION :: DET,A(3,3),M(3,3)

        M(1,1) = A(2,2)*A(3,3)-A(2,3)**2
        M(2,1) = -A(1,2)*A(3,3)+A(1,3)*A(2,3)
        M(3,1) = A(1,2)*A(2,3)-A(1,3)*A(2,2)

        M(1,2) = M(2,1)
        M(2,2) = A(1,1)*A(3,3)-A(1,3)**2
        M(3,2) = -A(1,1)*A(2,3)+A(1,3)*A(1,2)

        M(1,3) = M(3,1)
        M(2,3) = M(3,2)
        M(3,3) = A(1,1)*A(2,2)-A(1,2)**2

        DET = M(1,1)*A(1,1)+M(1,2)*A(1,2)+M(1,3)*A(1,3)

    END SUBROUTINE INVDET

    SUBROUTINE RHSM(NEQM,TT,YM,YPRIME)

        IMPLICIT NONE
        INTEGER :: I,J,K,PCNT,L,NEQM,CNT,NDIM,NAM,IND,MASSCNT
        DOUBLE PRECISION :: QIJ(3),U,TRUXX,GUX,GUG,GU(3*LRGMSZE,3*LRGMSZE),TT 
        DOUBLE PRECISION :: UPV(N3ATOM),UPM(3*LRGMSZE,3*LRGMSZE,N_MOL)
        DOUBLE PRECISION :: YM(NEQM),YPRIME(NEQM),G(3*LRGMSZE,3*LRGMSZE,N_MOL)

        CNT=2+N3ATOM
        NDIM=3*LRGMSZE
        IND=1

        DO I=1,N_MOL
            NAM=MOLDATA(4,IND)                       ! Inititalize G matrix
            IND=IND+NAM
            DO J=1,3*NAM
                G(J,J,I)=YM(CNT)
                CNT=CNT+1
                DO K=J+1,3*NAM
                    G(K,J,I)=YM(CNT)                                 
                    G(J,K,I)=YM(CNT)
                    CNT=CNT+1
                ENDDO
            ENDDO
        ENDDO

        U=0.0D0
        UPV=0.0D0
        UPM=0.0D0
        TRUXX=0.0D0

        DO PCNT=1,NPAIRS  
            I=PARRAY(1,PCNT)
            J=PARRAY(2,PCNT)
            DO K=1,3
                QIJ(K)=YM(3*(I-1)+K+1)-YM(3*(J-1)+K+1)
                IF(PBC) THEN                                            ! Periodic boundary conditions
                    IF(QIJ(K).GT.BL2) QIJ(K)=QIJ(K)-BL
                    IF(QIJ(K).LT.-BL2) QIJ(K)=QIJ(K)+BL
                ENDIF
            ENDDO
            IF(MOLDATA(1,I).EQ.MOLDATA(1,J)) THEN                                             ! Same molecule ?
                CALL RHSFC(I,J,QIJ,U,UPV,G,UPM,MOLDATA(1,I),MOLDATA(2,I),MOLDATA(2,J))
            ELSE
                CALL RHSSP(I,J,QIJ,U,UPV,G,UPM,MOLDATA(1,I),MOLDATA(1,J),MOLDATA(2,I),MOLDATA(2,J))     ! Different molecule.
            ENDIF
        ENDDO

        CALL FILLBLOCKS(UPM)                  ! Fill lower diagonal of G matrix

        DO I=1,N_MOL
            DO J=1,NDIM
                TRUXX=TRUXX+UPM(J,J,I)*G(J,J,I)
                DO K=J+1,NDIM
                    TRUXX=TRUXX+2.0D0*UPM(J,K,I)*G(K,J,I)
                ENDDO
            ENDDO
        ENDDO

        YPRIME(1)=-0.25D0*TRUXX-(U+E_ZERO)
        CNT=2
        IND=1

        DO I=1,N_MOL
            NDIM=3*MOLDATA(4,IND)
            DO J=1,NDIM
                GUX=0.0D0
                DO K=1,NDIM
                    GUX=GUX-G(J,K,I)*UPV(3*(IND-1)+K)
                ENDDO
                YPRIME(CNT)=GUX
                CNT=CNT+1
            ENDDO
            IND=IND+MOLDATA(4,IND)
        ENDDO

        IND=1
        MASSCNT=1

        !write (36,*) YPRIME(2+N3ATOM:1+N3ATOM)

        DO I=1,N_MOL
            GU=MATMUL(UPM(:,:,I),G(:,:,I))
            NDIM=3*MOLDATA(4,IND)
            IND=IND+MOLDATA(4,IND)
            DO J=1,NDIM
                DO K=J,NDIM
                    GUG=0.0D0
                    DO L=1,NDIM
                        GUG=GUG-G(J,L,I)*GU(L,K)
                    ENDDO
                    IF(J.EQ.K) THEN
                        GUG=GUG+MASSARRY(MASSCNT)
                        MASSCNT=MASSCNT+1
                    ENDIF
                    YPRIME(CNT)=GUG
                    CNT=CNT+1
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE RHSM

    SUBROUTINE RHSFC(IIND,JIND,QIJ,U,UPV,G,UPM,MNUM,TYPEI,TYPEJ)

        IMPLICIT NONE
        INTEGER :: I,J,K,IIND,JIND,IOF,JOF,MNUM,NG,TYPEI,TYPEJ
        DOUBLE PRECISION :: QIJ(3),UPV(3,N_ATOM),G(3*LRGMSZE,3*LRGMSZE,N_MOL),DETA,DETAG,ZQ(3),U,QZQ
        DOUBLE PRECISION :: UPM(3*LRGMSZE,3*LRGMSZE,N_MOL),A(3,3),M(3,3),AG(3,3),Z(3,3),EXPF,UX,UXX,AK,CK

        IOF=MOLDATA(3,IIND)
        JOF=MOLDATA(3,JIND)

        DO I=1,3
            DO J=I,3
                A(I,J)=G(IOF+I,IOF+J,MNUM)+G(JOF+I,JOF+J,MNUM)-G(IOF+I,JOF+J,MNUM)-G(IOF+J,JOF+I,MNUM)
            ENDDO
        ENDDO

        CALL INVDET(A,M,DETA)

        DETA=1.0D0/DETA
        A=M*DETA

        NG=INT(GDATA(1,TYPEI,TYPEJ,1))

        DO I=1,NG                                     ! Sum over intra-molecular Gaussians
            AG=A
            AK=GDATA(2*I+1,TYPEI,TYPEJ,1)
            CK=GDATA(2*I,TYPEI,TYPEJ,1)

            DO J=1,3
                AG(J,J)=AG(J,J)+AK
            ENDDO

            CALL INVDET(AG,M,DETAG)

            Z=-(AK**2/DETAG)*M

            DO J=1,3
                Z(J,J)=Z(J,J)+AK
            ENDDO

            Z(3,1)=Z(1,3)
            Z(2,1)=Z(1,2)
            Z(3,2)=Z(2,3)

            ZQ=MATMUL(Z,QIJ)
            QZQ=DOT_PRODUCT(QIJ,ZQ)
            ZQ=2.0D0*ZQ

            EXPF=SQRT(DETA/DETAG)*EXP(-QZQ)*CK
            U=U+EXPF

            DO J=1,3
                UX=-EXPF*ZQ(J)
                UPV(J,IIND)=UPV(J,IIND)+UX
                UPV(J,JIND)=UPV(J,JIND)-UX
                DO K=J,3
                    UXX=EXPF*(ZQ(J)*ZQ(K)-2.0D0*Z(J,K))
                    UPM(IOF+J,IOF+K,MNUM)=UPM(IOF+J,IOF+K,MNUM)+UXX             ! Fill upper triangle of blocks
                    UPM(JOF+J,JOF+K,MNUM)=UPM(JOF+J,JOF+K,MNUM)+UXX
                    UPM(IOF+J,JOF+K,MNUM)=UPM(IOF+J,JOF+K,MNUM)-UXX             ! Fill off diagonal blocks
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE RHSFC

    SUBROUTINE RHSSP(IIND,JIND,QIJ,U,UPV,G,UPM,MNUM1,MNUM2,TYPEI,TYPEJ)

        IMPLICIT NONE
        INTEGER :: I,J,K,IIND,JIND,IOF,JOF,MNUM1,MNUM2,NG,TYPEI,TYPEJ
        DOUBLE PRECISION :: QIJ(3),UPV(3,N_ATOM),G(3*LRGMSZE,3*LRGMSZE,N_MOL),DETA,DETAG,ZQ(3),U,QZQ
        DOUBLE PRECISION :: UPM(3*LRGMSZE,3*LRGMSZE,N_MOL),A(3,3),M(3,3),AG(3,3),Z(3,3),EXPF,UX,UXX,AK,CK

        IOF=MOLDATA(3,IIND)
        JOF=MOLDATA(3,JIND)

        DO I=1,3
            DO J=I,3
                A(I,J)=G(IOF+I,IOF+J,MNUM1)+G(JOF+I,JOF+J,MNUM2)
            ENDDO
        ENDDO

        CALL INVDET(A,M,DETA)

        DETA=1.0D0/DETA
        A=M*DETA

        NG=INT(GDATA(1,TYPEI,TYPEJ,2))

        DO I=1,NG                                     ! Sum over intra-molecular Gaussians
            AG=A
            AK=GDATA(2*I+1,TYPEI,TYPEJ,2)
            CK=GDATA(2*I,TYPEI,TYPEJ,2)
            DO J=1,3
                AG(J,J)=AG(J,J)+AK
            ENDDO

            CALL INVDET(AG,M,DETAG)

            Z=-(AK**2/DETAG)*M

            DO J=1,3
                Z(J,J)=Z(J,J)+AK
            ENDDO

            Z(3,1)=Z(1,3)
            Z(2,1)=Z(1,2)
            Z(3,2)=Z(2,3)

            ZQ=MATMUL(Z,QIJ)
            QZQ=DOT_PRODUCT(QIJ,ZQ)
            ZQ=2.0D0*ZQ

            EXPF=SQRT(DETA/DETAG)*EXP(-QZQ)*CK
            U=U+EXPF

            DO J=1,3
                UX=-EXPF*ZQ(J)
                UPV(J,IIND)=UPV(J,IIND)+UX
                UPV(J,JIND)=UPV(J,JIND)-UX
                DO K=J,3
                    UXX=EXPF*(ZQ(J)*ZQ(K)-2.0D0*Z(J,K))
                    UPM(IOF+J,IOF+K,MNUM1)=UPM(IOF+J,IOF+K,MNUM1)+UXX             ! Fill upper triangle of blocks
                    UPM(JOF+J,JOF+K,MNUM2)=UPM(JOF+J,JOF+K,MNUM2)+UXX
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE RHSSP

    SUBROUTINE FILLBLOCKS(UPM)
        IMPLICIT NONE
        INTEGER :: I,J,K,MSZE
        DOUBLE PRECISION :: UPM(3*LRGMSZE,3*LRGMSZE,N_MOL)

        MSZE=3*LRGMSZE

        DO I=1,N_MOL
            DO J=1,LRGMSZE
                DO K=J+1,LRGMSZE
                    UPM(3*(J-1)+2,3*(K-1)+1,I)=UPM(3*(J-1)+1,3*(K-1)+2,I)           ! Fill off diagonal blocks
                    UPM(3*(J-1)+3,3*(K-1)+1,I)=UPM(3*(J-1)+1,3*(K-1)+3,I)
                    UPM(3*(J-1)+3,3*(K-1)+2,I)=UPM(3*(J-1)+2,3*(K-1)+3,I)
                ENDDO
            ENDDO

            DO J=1,MSZE                  ! Fill lower triangle of molecular blocks       
                DO K=J+1,MSZE
                    UPM(K,J,I)=UPM(J,K,I)
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE FILLBLOCKS

    SUBROUTINE GAUSSENERGY(Q,UGAUSS)
        IMPLICIT NONE
        DOUBLE PRECISION :: Q(3,N_ATOM),UGAUSS,RSQ,AK,CK,QIJ(3)
        INTEGER :: I,J,K,PCNT,NG,TYPEI,TYPEJ

        UGAUSS=0.0D0

        DO PCNT=1,NPAIRS
            I=PARRAY(1,PCNT)
            J=PARRAY(2,PCNT)
            RSQ=0.0D0
            DO K=1,3
                QIJ(K)=Q(K,I)-Q(K,J)
                IF(PBC) THEN                                            ! Periodic boundary conditions
                    IF(QIJ(K).GT.BL2) QIJ(K)=QIJ(K)-BL
                    IF(QIJ(K).LT.-BL2) QIJ(K)=QIJ(K)+BL
                ENDIF
                RSQ=RSQ+QIJ(K)**2
            ENDDO
            TYPEI=MOLDATA(2,I)
            TYPEJ=MOLDATA(2,J)
            IF(MOLDATA(1,I).EQ.MOLDATA(1,J)) THEN
                NG=INT(GDATA(1,TYPEI,TYPEJ,1))
                DO K=1,NG
                    AK=GDATA(2*K+1,TYPEI,TYPEJ,1)
                    CK=GDATA(2*K,TYPEI,TYPEJ,1)
                    UGAUSS=UGAUSS+CK*EXP(-AK*RSQ)
                ENDDO
            ELSE
                NG=INT(GDATA(1,TYPEI,TYPEJ,2))
                DO K=1,NG
                    AK=GDATA(2*K+1,TYPEI,TYPEJ,2)
                    CK=GDATA(2*K,TYPEI,TYPEJ,2)
                    UGAUSS=UGAUSS+CK*EXP(-AK*RSQ)
                ENDDO
            ENDIF
        ENDDO

        UGAUSS=UGAUSS+E_ZERO

    END SUBROUTINE GAUSSENERGY
END MODULE VGWM
