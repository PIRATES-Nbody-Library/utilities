
C
C-------------------------------------------------------------------------
C                   Main program
C-------------------------------------------------------------------------
C
        PROGRAM MULTI_POSVEL

        include 'swift.inc'

        REAL*8, PARAMETER :: BINR=1.0d0
        REAL*8, PARAMETER :: RMAXBIN=3000.0d0
        INTEGER*4, PARAMETER :: NMAXBIN = NINT(RMAXBIN/BINR)
        INTEGER*4, PARAMETER :: NSMAX = 20
        INTEGER*4 NBIN(NMAXBIN)
        REAL*4, DIMENSION(:,:), ALLOCATABLE :: DATA4
        REAL*4, DIMENSION(-NPLMAX:NTPMAX) :: TA,TINC,TE,TR,TEX,TEY,TEZ,
     &                   TX,TY,TZ,TCAPOM,TOMEGA,TCAPM,TEPX,TEPY,TEPZ,
     &                   TVX,TVY,TVZ,TXH,TYH,TZH

        INTEGER*4 NT,IPL,K,IBIN,ISNAPTOT,NSNAP,IREAD,SIZ,ISNAP
        REAL*8 XJ(NPLMAX),YJ(NPLMAX),ZJ(NPLMAX),MASS(NPLMAX)
        REAL*8 VXJ(NPLMAX),VYJ(NPLMAX),VZJ(NPLMAX)
        REAL*8 VXH(NPLMAX),VYH(NPLMAX),VZH(NPLMAX)
        REAL*8 XH(NPLMAX),YH(NPLMAX),ZH(NPLMAX)
        REAL*8 VXB(NPLMAX),VYB(NPLMAX),VZB(NPLMAX)
        REAL*8 XB(NPLMAX),YB(NPLMAX),ZB(NPLMAX)
        REAL*8 EX(NPLMAX),EY(NPLMAX),EZ(NPLMAX),AINIT(NTPMAX*NSMAX)
        REAL*8 EPX(NPLMAX),EPY(NPLMAX),EPZ(NPLMAX)
        REAL*8 MAT(NPLMAX,NPLMAX),UMAT(NPLMAX,NPLMAX)
        real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
        real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)
        real*8 matp(NPLMAX,NTPMAX),umatp(NPLMAX,NTPMAX)
        REAL*8 ETA(NPLMAX),MU(NPLMAX),WIMPMFP(NTPMAX)
        REAL*8 EPXT,EPYT,EPZT,EXT,EYT,EZT,MPAS(3,3)
        INTEGER ISTAT(NTPMAX,NSTAT),OLOC(NPLMAX,NPLMAX)
        REAL*8 RSTAT(NTPMAX,NSTATR),etatp(NTPMAX)
        REAL*8 DTSNAP,RR,TIM,XT,YT,ZT,VXT,VYT,VZT,TDEB,q
        INTEGER NBOD,NTP,IERR,IFOL,ISTEP,KBIN,OLOCT(NPLMAX,NTPMAX)
        INTEGER IFLGCHK,IU,NLEFT,I,ID,J,ISIZ,ISNAPST,SNAPST,nref
        INTEGER IO_READ_HDR,IO_READ_LINE,NPO,NBOD0,NTP0(NSMAX)
        INTEGER IO_READ_HDR_R,IO_READ_LINE_R,II,ialpha,NS,NBSIMU
        INTEGER SNAPEXTR

        REAL*8 T0,TSTOP,DT,DTOUT,DTDUMP,DR,CX,CY,CZ,MATR(3,3),GM
        REAL*8 T,TMAX,EMAX(NTPMAX),VARPI,KR,KT,RCOLL,RX,RY,R,U

        REAL*8 RMIN,RMAX,RMAXU,QMIN,RPLSQ(NPLMAX),VX,VY,VZ,MSYS
        LOGICAL ERROR,WIMPS,EVAP,COLL
        LOGICAL*2 LCLOSE,OK
        REAL*8 A,E,INC,CAPOM,OMEGA,CAPM,CXB,CYB,CZB,CVXB,CVYB,CVZB
         REAL*8 ELH,ELK,ELP,ELQ,Z0,ALPHEVAP,REVAP,RADMIN
        REAL*8 SUMX,SUMY,SUMZ,SUMVX,SUMVY,SUMVZ,ETAL,j2rp2,j4rp4

        CHARACTER*256 OUTFILE(NSMAX),INPARFILE,INPLFILE,INTPFILE(NSMAX),
     &           FOPENSTAT,DIRO(NSMAX),DIRS,GNAME(NSMAX),DEV,EXTRFILE

 199    FORMAT(1X,6(1D14.6,1X))

        NPO=200
        print*,'ntpmax',ntpmax
        KBIN = NMAXBIN/10
        COLL=.FALSE.

        WRITE(*,*)'Enter number of individual simulations : '
        READ(*,*) NBSIMU
        write(*,*) nbsimu


        WRITE(*,*) ' '
        WRITE(*,*) 'Enter name of planet data file : '
        READ(*,'(A)')INPLFILE

        CALL IO_INIT_PL(INPLFILE,LCLOSE,IFLGCHK,NBOD0,
     &       MASS,XH,YH,ZH,VXH,VYH,VZH,RPLSQ,j2rp2,j4rp4)

        WRITE(*,*) ' '

        DO NS = 1,NBSIMU
          write(*,*) 'Enter name of test particle data file : '
          read(*,'(a)') intpfile(ns)
          call io_init_tp(intpfile(ns),ntp0(ns),xht,yht,zht,vxht,vyht,
     &     vzht,istat,rstat)
        END DO


c Stockage positions et vitesses des corps et particules pour un snapshot donne


          open(22,file='posvel_ECrvSngl_planets.dat',status='unknown')
          DO I=1,NBOD0
             write(22,199)xh(i),yh(i),zh(i),vxh(i),vyh(i),vzh(i)
          END DO
          close(22)

          do ns=1,nbsimu
          open(23,file='posvel_ECrvSngl_particles.dat',status='unknown')
          DO I=1,NTP0(ns)
             write(23,199)xht(i),yht(i),zht(i),vxht(i),vyht(i),vzht(i)
          END DO
          end do
          close(23)

c------------------------------------

        contains

C
C--------------------------------------------------------------------
C               Keplerian orbit stuff
C--------------------------------------------------------------------
C
        SUBROUTINE KEPLER_STUFF(GM,A,E,INC,CAPOM,OMEGA,CAPM,
     &                   R,CX,CY,CZ,VX,VY,VZ,EX,EY,EZ,EPX,EPY,EPZ)


        include 'swift.inc'

        REAL*8 GM,A,E,INC,CAPOM,OMEGA,CAPM
        REAL*8 R,RX,RY,CX,CY,CZ,VRX,VRY,VX,VY,VZ
        REAL*8 EX,EY,EZ,EPX,EPY,EPZ,U,F,F1,F2,F3,DELA
        LOGICAL OK
        INTEGER I

        EX = COS(OMEGA)*COS(CAPOM)-SIN(OMEGA)*COS(INC)*SIN(CAPOM)
        EY = COS(OMEGA)*SIN(CAPOM)+SIN(OMEGA)*COS(INC)*COS(CAPOM)
        EZ = SIN(OMEGA)*SIN(INC)
        EPX = -SIN(OMEGA)*COS(CAPOM)-COS(OMEGA)*COS(INC)*SIN(CAPOM)
        EPY = -SIN(OMEGA)*SIN(CAPOM)+COS(OMEGA)*COS(INC)*COS(CAPOM)
        EPZ = COS(OMEGA)*SIN(INC)
        U = CAPM
        OK = .FALSE.
        I=0
        DO WHILE (.NOT.OK)
          F = U-E*SIN(U)-CAPM
          F1 = 1.-E*COS(U)
          F2 = U-CAPM-F
          F3 = 1.-F1
          DELA = -F/F1
          DELA = -F/(F1+0.5*DELA*F2)
          DELA = -F/(F1+0.5*DELA*F2+DELA*DELA*F3/6.)
          U = U+DELA
          I = I+1
          OK = ((ABS(DELA).LT.1e-12).OR.(I.GE.500))
        END DO
        R = A*(1.-E*COS(U))
        RX = A*(COS(U)-E)
        RY = A*SQRT(1.-E*E)*SIN(U)
        VRX = -SQRT(GM/A)*A/R*SIN(U)
        VRY = SQRT(GM/A)*A/R*COS(U)
        CX = EX*RX+EPX*RY
        CY = EY*RX+EPY*RY
        CZ = EZ*RX+EPZ*RY
        VX = EX*VRX+EPX*VRY
        VY = EY*VRX+EPY*VRY
        VZ = EZ*VRX+EPZ*VRY
        END subroutine

C--------------------------------------
        SUBROUTINE KEPLER(CAPM,E,U)

        IMPLICIT NONE

        REAL*8 U,DELA,CAPM,F1,F2,F3,F,E
        LOGICAL OK
        INTEGER*4 I

        U = CAPM
        OK = .FALSE.
        I=0
        DO WHILE (.NOT.OK)
          F = U-E*SIN(U)-CAPM
          F1 = 1.-E*COS(U)
          F2 = U-CAPM-F
          F3 = 1.-F1
          DELA = -F/F1
          DELA = -F/(F1+0.5*DELA*F2)
          DELA = -F/(F1+0.5*DELA*F2+DELA*DELA*F3/6.)
          U = U+DELA
          I = I+1
          OK = ((ABS(DELA).LT.1e-12).OR.(I.GE.500))
        END DO
        END subroutine

C--------------------------------------------------------------------
C               For invariable plane
C--------------------------------------------------------------------
C
        SUBROUTINE INVAR(NBOD,MASS,X,Y,Z,VX,VY,VZ,A)

        include 'swift.inc'

        INTEGER I,NBOD

        REAL*8 MASS(NBOD),X(NBOD),Y(NBOD),Z(NBOD),
     &                    VX(NBOD),VY(NBOD),VZ(NBOD),
c     &                    XT(NPLMAX),YT(NPLMAX),ZT(NPLMAX),
c     &                    VXT(NPLMAX),VYT(NPLMAX),VZT(NPLMAX),
     &                    C1,C2,C3,C,BOT,A(3,3)

        C1 = SUM(MASS(1:NBOD)*
     &           (Y(1:NBOD)*VZ(1:NBOD)-Z(1:NBOD)*VY(1:NBOD)))
        C2 = SUM(MASS(1:NBOD)*
     &           (Z(1:NBOD)*VX(1:NBOD)-X(1:NBOD)*VZ(1:NBOD)))
        C3 = SUM(MASS(1:NBOD)*
     &           (X(1:NBOD)*VY(1:NBOD)-Y(1:NBOD)*VX(1:NBOD)))

        C = SQRT( C1*C1 + C2*C2 + C3*C3 )
        C1 = C1/C
        C2 = C2/C
        C3 = C3/C
        BOT = 1.0D0/(1.0D0 + C3)

        A(1,1) = 1.0D0 - C1*C1*BOT
        A(1,2) = -1.0D0*C1*C2*BOT
        A(1,3) = -1.0D0*C1

        A(2,1) = -1.0D0*C1*C2*BOT
        A(2,2) = 1.0D0 - C2*C2*BOT
        A(2,3) = -1.0D0*C2

        A(3,1) = C1
        A(3,2) = C2
        A(3,3) = C3


c       DO I=1,NBOD
c          XT(I) = A(1,1)*X(I) + A(1,2)*Y(I) + A(1,3)*Z(I)
c          YT(I) = A(2,1)*X(I) + A(2,2)*Y(I) + A(2,3)*Z(I)
c          ZT(I) = A(3,1)*X(I) + A(3,2)*Y(I) + A(3,3)*Z(I)
c          VXT(I) = A(1,1)*VX(I) + A(1,2)*VY(I) + A(1,3)*VZ(I)
c          VYT(I) = A(2,1)*VX(I) + A(2,2)*VY(I) + A(2,3)*VZ(I)
c          VZT(I) = A(3,1)*VX(I) + A(3,2)*VY(I) + A(3,3)*VZ(I)
c        END DO
c
C..    Check
c
c        C1 = SUM(MASS(1:NBOD)*
c     &           (YT(1:NBOD)*VZT(1:NBOD)-ZT(1:NBOD)*VYT(1:NBOD)))
c        C2 = SUM(MASS(1:NBOD)*
c     &           (ZT(1:NBOD)*VXT(1:NBOD)-XT(1:NBOD)*VZT(1:NBOD)))
c        C3 = SUM(MASS(1:NBOD)*
c     &           (XT(1:NBOD)*VYT(1:NBOD)-YT(1:NBOD)*VXT(1:NBOD)))
c        WRITE(*,*) ' C1,C2,C3 ',C1,C2,C3


        END subroutine
c-------------------------------------------------------------------------

      end program multi_posvel
