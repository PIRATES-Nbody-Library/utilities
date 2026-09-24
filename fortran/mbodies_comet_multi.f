C
C-------------------------------------------------------------------------
C                   Main program
C-------------------------------------------------------------------------
C
        PROGRAM MULTI_COMET

        include 'swift.inc'

        REAL*8, PARAMETER :: BINR=1.0d0
        REAL*8, PARAMETER :: RMAXBIN=3000.0d0
        INTEGER*4, PARAMETER :: NMAXBIN = NINT(RMAXBIN/BINR)
        INTEGER*4, PARAMETER :: NSMAX = 20
        INTEGER*4 NBIN(NMAXBIN)
        REAL*4, DIMENSION(:,:), ALLOCATABLE :: DATA4,comet
        integer*4, dimension(:), allocatable :: idtot,idcomet
        REAL*4, DIMENSION(-NPLMAX:NTPMAX) :: TA,TINC,TE,TR,TEX,TEY,TEZ,
     &                   TX,TY,TZ,TCAPOM,TOMEGA,TCAPM,TEPX,TEPY,TEPZ

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
        INTEGER*4 ISTAT(NTPMAX,NSTAT),OLOC(NPLMAX,NPLMAX),icom,ncom
        REAL*8 RSTAT(NTPMAX,NSTATR),etatp(NTPMAX),qmax
        REAL*8 DTSNAP,RR,TIM,XT,YT,ZT,VXT,VYT,VZT,TDEB,q
        INTEGER*4 NBOD,NTP,IERR,IFOL,ISTEP,KBIN,OLOCT(NPLMAX,NTPMAX)
        INTEGER*4 IFLGCHK,IU,NLEFT,I,ID,J,ISIZ,ISNAPST,SNAPST,nref
        INTEGER*4 IO_READ_HDR,IO_READ_LINE,NPO,NBOD0,NTP0(NSMAX)
        INTEGER*4 IO_READ_HDR_R,IO_READ_LINE_R,II,ialpha,NS,NBSIMU

        REAL*8 T0,TSTOP,DT,DTOUT,DTDUMP,DR,CX,CY,CZ,MATR(3,3),GM
        REAL*8 T,TMAX,EMAX(NTPMAX),VARPI,KR,KT,RCOLL,RX,RY,R,U

        REAL*8 RMIN,RMAX,RMAXU,QMIN,RPLSQ(NPLMAX),VX,VY,VZ,MSYS
        LOGICAL ERROR,WIMPS,EVAP,COLL
        LOGICAL*2 LCLOSE,OK
        REAL*8 A,E,INC,CAPOM,OMEGA,CAPM,CXB,CYB,CZB,CVXB,CVYB,CVZB
         REAL*8 ELH,ELK,ELP,ELQ,Z0,ALPHEVAP,REVAP,RADMIN
        REAL*8 SUMX,SUMY,SUMZ,SUMVX,SUMVY,SUMVZ,ETAL,j2rp2,j4rp4

        CHARACTER*256 OUTFILE(NSMAX),INPARFILE,INPLFILE,INTPFILE(NSMAX),
     &           FOPENSTAT,DIRO(NSMAX),DIRS,GNAME(NSMAX),DEV,EXTRFILE,
     &	         cometfile

 99     FORMAT(1X,1D14.6,1X)
 199    format(1X,6(1D14.6,1X))
 999    FORMAT(1X,7(1D14.6,1X))

        NPO=200
        print*,' '
        KBIN = NMAXBIN/10
        COLL=.FALSE.

        write(*,*) 'Enter name of outfile:'
        read(*,'(a)') cometfile
	    print*,cometfile
	    write(*,*) 'Enter value for periastron investigated:'
        read(*,*)qmax
        print*,qmax

        WRITE(*,*)'Enter number of individual simulations : '
        READ(*,*) NBSIMU
        DO I=1,NBSIMU
          WRITE(*,*)'Enter name of parameter data file : ',I
          READ(*,'(a)')inparfile
          call io_init_param(inparfile,t0,tstop,dt,dtout,dtdump,
     &         iflgchk,rmin,rmax,rmaxu,qmin,lclose,diro(i),dirs,
     &         gname(i),outfile(i),fopenstat)
        END DO

        WRITE(*,*) 'Enter beginning time (estimated) :'
        READ(*,*)TDEB
        IF (TDEB.GT.TSTOP) THEN
          PRINT*,'Time too large'
          STOP
        END IF

        WRITE(*,*) 'Enter final time (estimated) :'
        READ(*,*)TIM
        IF (TIM.GT.TSTOP) THEN
          PRINT*,'Time too large'
          STOP
        END IF
        WRITE(*,*) 'Enter timestep between snapshots :'
        READ(*,*)DTSNAP
        NSNAP = NINT(DTSNAP/DTOUT)
        SNAPST = FLOOR(TDEB/DTSNAP)
        ISNAPTOT = FLOOR(TIM/DTSNAP)+1-SNAPST
        print*,isnaptot
        WRITE(*,*) ' '
        WRITE(*,*) 'Enter name of planet data file : '
        READ(*,'(A)')INPLFILE

        CALL IO_INIT_PL(INPLFILE,LCLOSE,IFLGCHK,NBOD0,
     &       MASS,XH,YH,ZH,VXH,VYH,VZH,RPLSQ,j2rp2,j4rp4)

c        WRITE(*,*) ' '
c        WRITE(*,*) 'Enter name of extract file: '
c        READ(*,'(A)') EXTRFILE
        IU = 20
        DR = 180.0/PI
c        DEV=TRIM(DIRO(1))//'/'//TRIM(EXTRFILE)//'.dat'

        IF (BTEST(IFLGCHK,0)) THEN
           WRITE(*,*) ' Reading INTEGER*2 binary files '
        ELSE IF(BTEST(IFLGCHK,1)) THEN
           WRITE(*,*) ' Reading REAL*4 binary files '
        ELSE
           WRITE(*,*) ' Error: no binary file format specified '
           WRITE(*,*) '        in param file '
           STOP
        ENDIF
        OLOC=0
        DO J=2,nbod0
            OLOC(J,1)=-1
            OLOC(J,J)=1
        END DO

        DO NS = 1,NBSIMU
c          print*,ns
          write(*,*) 'Enter name of test particle data file : '
          read(*,'(a)') intpfile(ns)
          call io_init_tp(intpfile(ns),ntp0(ns),xht,yht,zht,vxht,vyht,
     &     vzht,istat,rstat)
          do i=1,ntp0(ns)
            j=sum(ntp0(1:ns-1))
            call orbel_xv2aeq(xht(i),yht(i),zht(i),
     &        vxht(i),vyht(i),vzht(i),mass(1),ialpha,ainit(j+i),e,q)
          end do
        END DO



        allocate(idtot(NTPMAX))
        allocate(idcomet(NTPMAX))
        allocate(comet(NTPMAX*isnaptot,7))

c        allocate(idcomet(sum(ntp0(1:ns))*isnaptot))
c        allocate(comet(sum(ntp0(1:ns))*isnaptot,5))

        icom=0
        idtot=0
        idcomet=0
        comet=0

        SIZ = MAX(NBOD0,3)+NBOD0

        DO NS = 1,NBSIMU
          OPEN(UNIT=IU, FILE=trim(diro(ns))//'/'//OUTFILE(NS),
     &                    STATUS='OLD',FORM='UNFORMATTED')
          mpas =0.0d0
          IREAD = 0
          ISNAP = 0
          ISNAPST = 0
          NTP = 0
          T = -1.0d0
          IPL = 2
c          OLOC=0
          DO WHILE(T.LT.TIM)
            IF (BTEST(IFLGCHK,0))  THEN        ! Bit 0 i set = Integer*2
              IERR = IO_READ_HDR(IU,T,NBOD,NLEFT)
c              print*,'flag 0',T,nbod,nleft
            ELSE
              IERR = IO_READ_HDR_R(IU,T,NBOD,NLEFT)
c              print*,'flag non 0',t,nbod,nleft
            ENDIF
            IF (IERR.NE.0) EXIT
c            IF (MOD(IREAD,50).EQ.0)
c     &
            print*,'T=',T,'Simu',NS,' NLEFT=',nleft,siz
            OK=(MOD(IREAD,NSNAP).EQ.0)
            IF (OK.AND.(ISNAPST.LT.SNAPST)) THEN
              ISNAPST = ISNAPST+1
              OK = .FALSE.
            END IF
            DO I=2,NBOD
              IF (BTEST(IFLGCHK,0))  THEN    ! Bit 0 i set = Integer*2
                IERR = IO_READ_LINE(IU,ID,A,E,INC,CAPOM,OMEGA,CAPM)
              ELSE
                IERR = IO_READ_LINE_R(IU,ID,A,E,INC,CAPOM,OMEGA,CAPM)
              ENDIF
              IF(IERR.NE.0) EXIT
            END DO
            NTP=0
            DO I=1,NLEFT
              IF(BTEST(IFLGCHK,0))  THEN ! BIT 0 IS SET
                IERR = IO_READ_LINE(IU,ID,A,E,INC,CAPOM,OMEGA,CAPM)
              ELSE
                IERR = IO_READ_LINE_R(IU,ID,A,E,INC,CAPOM,OMEGA,CAPM)
              END IF
              IF(IERR.NE.0) EXIT
              IF (OK) THEN
                IF ((E.LT.1.0d0).and.(e.ge.0.0d0)) THEN
                  NTP=NTP+1
                END IF
              END IF
            END DO
            IF (OK) THEN
              IF (NS.EQ.1) THEN
c              IF (NS.EQ.nref) THEN
                SIZ = SIZ+NBOD+NTP+1
              ELSE
                SIZ = SIZ+NTP
              END IF
            END IF
            IREAD=IREAD+1
          END DO
          CLOSE(IU)
        END DO
        ALLOCATE(DATA4(SIZ,11))
        DATA4 = SNGL(0.0d0)
        ISIZ = 0

        DO I=1,NBOD0
          DATA4(I,1) = SNGL(MASS(I))
        END DO
        DO I=2,NBOD0
          DATA4(I,2) = SNGL(MASS(1))
        END DO
        DO I=2,NBOD0
          DATA4(I,3) = SNGL(MASS(I))
        END DO
        DO I=1,3
          DO J=1,3
            DATA4(J,I+3) = SNGL(MPAS(J,I))
          END DO
        END DO
        DATA4(1,7) = REAL(ISNAPTOT)
        print*,nbod0,sum(ntp0(1:nbsimu))
        DATA4(2,7) = REAL(NBOD0)
        DATA4(3,7) = REAL(SUM(NTP0(1:NBSIMU)))
        DO I=1,NBOD0
          DO J=1,NBOD0
            DATA4(MAX(NBOD0,3)+J,I) = REAL(OLOC(J,I))
          END DO
        END DO
        ISIZ = MAX(NBOD0,3)+NBOD0

        DO NS = 1,NBSIMU
          OPEN(UNIT=IU+NS, FILE=TRIM(DIRO(NS))//'/'//OUTFILE(NS),
     &                    STATUS='OLD',FORM='UNFORMATTED')
        END DO
        IREAD = 0
        ISNAP = 0
        ISNAPST = 0
        NTP = 0
        T = -1.0d0

c        print*,oloc(1:nbod0,1:nbod0)

        DO WHILE(T.LT.TIM)
          OK = (MOD(IREAD,NSNAP).EQ.0)
          IF (OK.AND.(ISNAPST.LT.SNAPST)) THEN
            ISNAPST = ISNAPST+1
            OK = .FALSE.
          END IF

          NTP = 0
          DO NS = 1,NBSIMU
            IF (BTEST(IFLGCHK,0))  THEN        ! Bit 0 i set = Integer*2
              IERR = IO_READ_HDR(IU+NS,T,NBOD,NLEFT)
            ELSE
              IERR = IO_READ_HDR_R(IU+NS,T,NBOD,NLEFT)
            ENDIF
            IF (IERR.NE.0) EXIT
            IF (MOD(IREAD,20).EQ.0)
     &              print*,'T=',T,'Simu',NS,' NLEFT=',nleft
            XH(1) = 0.0d0
            YH(1) = 0.0d0
            ZH(1) = 0.0d0
            VXH(1) = 0.0d0
            VYH(1) = 0.0d0
            VZH(1) = 0.0d0
            DO I=2,NBOD
              IF (BTEST(IFLGCHK,0))  THEN    ! Bit 0 i set = Integer*2
                IERR = IO_READ_LINE(IU+NS,ID,A,E,INC,CAPOM,OMEGA,CAPM)
              ELSE
                IERR = IO_READ_LINE_R(IU+NS,ID,A,E,INC,CAPOM,OMEGA,CAPM)
              ENDIF
              IF(IERR.NE.0) EXIT
              IF ((OK).AND.(NS.EQ.1)) THEN
                TA(ID)=SNGL(A)
                TE(ID)=SNGL(E)
                TINC(ID)=SNGL(INC*DR)
                TCAPOM(ID)=SNGL(CAPOM*DR)
                TOMEGA(ID)=SNGL(OMEGA*DR)
                TCAPM(ID)=SNGL(CAPM*DR)
                GM=MASS(-ID)+MASS(1)
                CALL KEPLER_STUFF(GM,A,E,INC,CAPOM,OMEGA,CAPM,
     &    R,XH(-ID),YH(-ID),ZH(-ID),VXH(-ID),VYH(-ID),VZH(-ID),
     &       EX(-ID),EY(-ID),EZ(-ID),EPX(-ID),EPY(-ID),EPZ(-ID))
              END IF
            END DO
            IF ((OK).AND.(NS.EQ.1)) THEN
              call coord_h2b(nbod,mass,xh,yh,zh,vxh,vyh,vzh,
     &      xb,yb,zb,vxb,vyb,vzb,msys)
              DO I=1,NBOD
                TR(-I) = SNGL(SQRT(XB(I)*XB(I)+YB(I)*YB(I)+ZB(I)*ZB(I)))
                TX(-I) = SNGL(XB(I))
                TY(-I) = SNGL(YB(I))
                TZ(-I) = SNGL(ZB(I))
                TEX(-I) = SNGL(EX(I))
                TEY(-I) = SNGL(EY(I))
                TEZ(-I) = SNGL(EZ(I))
                TEPX(-I) = SNGL(EPX(I))
                TEPY(-I) = SNGL(EPY(I))
                TEPZ(-I) = SNGL(EPZ(I))
              END DO
            END IF

            DO I=1,NLEFT
              IF(BTEST(IFLGCHK,0))  THEN ! BIT 0 IS SET
                IERR = IO_READ_LINE(IU+NS,ID,A,E,INC,CAPOM,OMEGA,CAPM)
              ELSE
                IERR = IO_READ_LINE_R(IU+NS,ID,A,E,INC,CAPOM,OMEGA,CAPM)
              ENDIF
              ID = ID + SUM(NTP0(1:NS-1))
              IF(IERR.NE.0) EXIT
              IF (OK) THEN
                IF ((E.LT.1.0d0).and.(e.ge.0.0d0)) THEN
                  NTP=NTP+1
                  idtot(ntp)=id
                  TA(NTP)=SNGL(A)
                  TE(NTP)=SNGL(E)
                  TINC(NTP)=SNGL(INC*DR)
                  TCAPOM(NTP)=SNGL(CAPOM*DR)
                  TOMEGA(NTP)=SNGL(OMEGA*DR)
                  TCAPM(NTP)=SNGL(CAPM*DR)
                  GM = MASS(1)
                  CALL KEPLER_STUFF(GM,A,E,INC,CAPOM,OMEGA,CAPM,
     &               R,CX,CY,CZ,VX,VY,VZ,EXT,EYT,EZT,EPXT,EPYT,EPZT)
                  call coord_h2b_tp(1,cx,cy,cz,vx,vy,vzh,
     &               xb(1),yb(1),zb(1),vxb(1),vyb(1),vzb(1),
     &               cxb,cyb,czb,cvxb,cvyb,cvzb)
                  TR(NTP) = SNGL(AINIT(ID))
                  TX(NTP) = SNGL(CXB)
                  TY(NTP) = SNGL(CYB)
                  TZ(NTP) = SNGL(CZB)
                END IF
              END IF
            END DO
          END DO

          IF (OK) THEN
            ISNAP = ISNAP+1

            IF (ISNAP.le.ISNAPTOT) THEN

              DATA4(ISIZ+1,1) = SNGL(T)
              DATA4(ISIZ+1,2) = REAL(NBOD)
              DATA4(ISIZ+1,3) = REAL(NTP)
              DATA4(ISIZ+1,4) = REAL(ISNAP)
              DATA4(ISIZ+1,5) = REAL(ISNAPTOT)
              ISIZ = ISIZ+1
                DO I=-1,-NBOD,-1
                  DATA4(ISIZ-I,1)=TA(I)
                  DATA4(ISIZ-I,2)=TE(I)
                  DATA4(ISIZ-I,3)=TEX(I)
                  DATA4(ISIZ-I,4)=TEY(I)
                  DATA4(ISIZ-I,5)=TEZ(I)
                  DATA4(ISIZ-I,6)=TEPX(I)
                  DATA4(ISIZ-I,7)=TEPY(I)
                  DATA4(ISIZ-I,8)=TEPZ(I)
                  DATA4(ISIZ-I,9)=TX(I)
                  DATA4(ISIZ-I,10)=TY(I)
                  DATA4(ISIZ-I,11)=TZ(I)
                END DO
                DO I=1,NTP
                  DATA4(ISIZ+I+NBOD,1)=TA(I)
                  DATA4(ISIZ+I+NBOD,2)=TE(I)
                  DATA4(ISIZ+I+NBOD,3)=TINC(I)
                  DATA4(ISIZ+I+NBOD,4)=TCAPOM(I)
                  DATA4(ISIZ+I+NBOD,5)=TOMEGA(I)
                  DATA4(ISIZ+I+NBOD,6)=TCAPM(I)
                  DATA4(ISIZ+I+NBOD,7)=TR(I)
                  DATA4(ISIZ+I+NBOD,9)=TX(I)
                  DATA4(ISIZ+I+NBOD,10)=TY(I)
                  DATA4(ISIZ+I+NBOD,11)=TZ(I)

c Here we will store information on particles which periastron is inferior to the value given from the start

                IF (TA(I)*(1.D0-TE(I)).le.qmax) THEN
                  icom=icom+1
                  comet(icom,1)=dble(isnap)
                  comet(icom,2)=dble(idtot(i))
                  comet(icom,3)=ainit(idtot(i))
                  comet(icom,4)=ta(i)*(1-te(i))
                  comet(icom,5)=te(i)
                  comet(icom,6)=tinc(i)
                  comet(icom,7)=tomega(i)
                  idcomet(idtot(i))=1
                END IF
              END DO

            END IF

c            print*,isnap
c            print*,comet(1,1)

            ISIZ = ISIZ+NBOD+NTP
c print*,'Creating snapshot at T=',SNGL(T),' npts=',NBOD+NTP
          END IF
          IREAD=IREAD+1
        END DO
        DO NS = 1,NBSIMU
          CLOSE(IU+NS)
        END DO
        print*,siz,isnaptot

	    ncom=icom

        OPEN(18,FILE=cometfile,status='unknown')
        DO I=1,ncom
           if (comet(i,7).gt.360.0D0) then
              comet(i,7)=comet(i,5)-360.0D0
           end if
           WRITE(18,999)comet(i,1),comet(i,2),comet(i,3),
     &       comet(i,4),comet(i,5),comet(i,6),comet(i,7)
        END DO
        CLOSE(18)

        DEALLOCATE(DATA4)
        deallocate(idtot)
        deallocate(comet)
        deallocate(idcomet)


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
        INTEGER*4 I

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
        END SUBROUTINE

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
        END SUBROUTINE

C--------------------------------------------------------------------
C               For invariable plane
C--------------------------------------------------------------------
C
        SUBROUTINE INVAR(NBOD,MASS,X,Y,Z,VX,VY,VZ,A)

        include 'swift.inc'

        INTEGER*4 I,NBOD

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


        END SUBROUTINE
c-------------------------------------------------------------------------

        END PROGRAM MULTI_COMET
