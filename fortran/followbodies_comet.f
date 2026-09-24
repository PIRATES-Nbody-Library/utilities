
C
C-------------------------------------------------------------------------
C                   Main program
C-------------------------------------------------------------------------
C
        PROGRAM FOLLOWBODIES_FOMC

        include 'swift.inc'

        REAL*4, DIMENSION(:,:,:), ALLOCATABLE :: DATA4
        INTEGER*4, PARAMETER :: NFOLMAX = 20
        REAL*4, DIMENSION(500000,-NPLMAX:NFOLMAX) :: TT,TA,TINC,TE,
     &                                         TCAPOM,TOMEGA,TCAPM

        INTEGER*4 NB(-NPLMAX:NFOLMAX),NBMAX,NT
        INTEGER*4 NBC(-NPLMAX:NFOLMAX),NCG,J,JFOL
	REAL*8 XJ(NPLMAX),YJ(NPLMAX),ZJ(NPLMAX),MASS(NPLMAX)
	REAL*8 VXJ(NPLMAX),VYJ(NPLMAX),VZJ(NPLMAX),TTINC(NTPMAX)
        REAL*8 TTX(NTPMAX),TTY(NTPMAX),TTZ(NTPMAX),TTH(NTPMAX)
        REAL*8 TTR(NTPMAX)
	INTEGER ISTAT(NTPMAX,NSTAT),NPFOL,IDFOL(NFOLMAX)
        REAL*8 RSTAT(NTPMAX,NSTATR),RR,TTA(NTPMAX),TTE(NTPMAX)
	REAL*8 MOYM,XX,YY,ZZ,XP,DTT,TEV,EH,TEVII,CEPS,SEPS,EPS
	INTEGER*4 NBOD,NTP,IERR,IFOL,ISTEP
	INTEGER*4 IFLGCHK,IU,NLEFT,I,ID,OLOC(NPLMAX,NPLMAX)
        INTEGER*4 IO_READ_HDR,IO_READ_LINE,IREAD
        INTEGER*4 IO_READ_HDR_R,IO_READ_LINE_R,II

	REAL*8 T0,TSTOP,DT,DTOUT,DTDUMP,DR,TTCAPOM(NTPMAX),IMOY
	REAL*8 T,TMAX,EMAX(NTPMAX),VARPI,KR,KT,RCOLL,AMOY,EMOY
        REAL*8 MAT(NPLMAX,NPLMAX),UMAT(NPLMAX,NPLMAX),OR2,IR2
        REAL*8 ETA(NPLMAX),MU(NPLMAX),IRMOY,OMOY,KXMOY,KYMOY,KZMOY
        REAL*8 PPP(3),SIGP(3),CHI2,DQA,DQAA,RHA,SIGA

	REAL*8 RMIN,RMAX,RMAXU,QMIN,RPLSQ(NPLMAX),KX,KY,KZ,KXR,KYR
        LOGICAL ERROR,WIMPS,EVAP,COLL
        LOGICAL*2 LCLOSE,OK
        REAL*8 A,E,INC,CAPOM,OMEGA,CAPM,J2RP2,J4RP4,KZR,SIGI,SIGO
        REAL*8 ELH,ELK,ELP,ELQ,Z0,ALPHEVAP,REVAP,RADMIN,ORMOY,RHO

	CHARACTER*256 OUTFILE,INPARFILE,INPLFILE,INTPFILE,FOPENSTAT,
     &               DIRO,DIRS,GNAME,DEV,folfile,followfile

 199    format(1X,4(1D14.6,1X))

        COLL=.FALSE.
	WRITE(*,*) 'Enter name of parameter data file : '
	READ(*,'(a)') inparfile
        call io_init_param_multirun(inparfile,t0,tstop,dt,dtout,dtdump,
     &         iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile,diro,dirs,
     &         gname,fopenstat)

        WRITE(*,*) 'Enter final time (estimated) :'
        READ(*,*)TSTOP
        NT=NINT(TSTOP/DTOUT)
        print*,nt
        WRITE(*,*) 'One data over :'
        READ(*,*)NCG
	WRITE(*,*) ' '
	WRITE(*,*) 'Enter name of planet data file : '
	READ(*,'(A)') INPLFILE
        call io_init_pl(inplfile,lclose,iflgchk,nbod,mass,xj,yj,zj,
     &     vxj,vyj,vzj,rplsq,j2rp2,j4rp4)
        IU = 20
        NB = 0
        NBC = 0
        DR = 180.0D0/PI
        T=-1.0d0

        NPFOL = 0
	WRITE(*,*) ' '
        write(*,*) 'Enter name of followbodies file:'
        read(*,'(a)')followfile
c        write(*,*) 'Enter name of file with particles to follow'
c        read(*,'(a)')folfile
c        open(unit=123,file=folfile,status='unknown')
	WRITE(*,*) 'Enter number of particles to follow : '
c        READ(123,*)NPFOL
        read(*,*)npfol
        print*,npfol
        IF (NPFOL.GT.NFOLMAX) THEN
          WRITE(*,*)'Maximum number of particles allowed = ',NFOLMAX
          STOP
        END IF
        IF (NPFOL.GT.0) THEN
  	  WRITE(*,*) 'Enter ID''s of particles to follow : '
          READ(*,*)IDFOL(1:NPFOL)
c          do i=1,npfol
c             read(123,*)idfol(i)
c             print*,idfol(i)
c          end do
        END IF


        IF (BTEST(IFLGCHK,0)) THEN
           WRITE(*,*) ' Reading an INTEGER*2 binary file '
        ELSE IF(BTEST(IFLGCHK,1)) THEN
           WRITE(*,*) ' Reading a REAL*4 binary file '
        ELSE
           WRITE(*,*) ' Error: no binary file format specified '
           WRITE(*,*) '        in param file '
           STOP
        ENDIF
        NTP=0

        OPEN(UNIT=IU, FILE=TRIM(DIRO)//'/'//OUTFILE,
     &                    STATUS='OLD',FORM='UNFORMATTED')
        open(unit=10,file='toto.dat',status='UNKNOWN')
        IREAD = 0
        MAIN_LOOP: DO WHILE (T.LT.TSTOP)
          IF (BTEST(IFLGCHK,0))  THEN        ! Bit 0 i set = Integer*2
            IERR = IO_READ_HDR(IU,T,NBOD,NLEFT)
          ELSE
            IERR = IO_READ_HDR_R(IU,T,NBOD,NLEFT)
          ENDIF
          IF (IERR.NE.0) EXIT MAIN_LOOP
          if (mod(iread,50).eq.0) print*,'T = ',t,' NBOD=',nbod
          iread = iread+1
          write(10,*)sngl(T),nleft
          DO I=2,NBOD
            IF (BTEST(IFLGCHK,0))  THEN    ! Bit 0 i set = Integer*2
              IERR = IO_READ_LINE(IU,ID,A,E,INC,CAPOM,OMEGA,CAPM)
            ELSE
              IERR = IO_READ_LINE_R(IU,ID,A,E,INC,CAPOM,OMEGA,CAPM)
            ENDIF
            IF(IERR.NE.0) EXIT MAIN_LOOP
            IF (NBC(ID).EQ.0) THEN
              NB(ID)=NB(ID)+1
              TT(NB(ID),ID)=SNGL(T)
              TA(NB(ID),ID)=SNGL(A)
              TE(NB(ID),ID)=SNGL(E)
              TINC(NB(ID),ID)=SNGL(INC*DR)
              TCAPOM(NB(ID),ID)=SNGL(CAPOM*DR)
              TOMEGA(NB(ID),ID)=SNGL(OMEGA*DR)
              TCAPM(NB(ID),ID)=SNGL(CAPM*DR)
            END IF
            NBC(ID)=NBC(ID)+1
            IF (NBC(ID).EQ.NCG) NBC(ID)=0
          END DO
          NTP=0
          DO I=1,NLEFT
            IF(BTEST(IFLGCHK,0))  THEN ! BIT 0 IS SET
              IERR = IO_READ_LINE(IU,ID,A,E,INC,CAPOM,OMEGA,CAPM)
            ELSE
              IERR = IO_READ_LINE_R(IU,ID,A,E,INC,CAPOM,OMEGA,CAPM)
            ENDIF
            IF(IERR.NE.0) EXIT MAIN_LOOP
            IF ((E.LT.1.0d0).and.(e.gt.0.0d0)) THEN
              NTP=NTP+1
              TTA(NTP)=A
              TTE(NTP)=E
              TTINC(NTP)=INC
              TTCAPOM(NTP)=CAPOM
              CALL KEPLER_STUFF(A,E,INC,CAPOM,OMEGA,CAPM,
     &                          TTX(NTP),TTY(NTP),TTZ(NTP))
            END IF
            OK = .FALSE.
            DO J=1,NPFOL
              IF (ID.EQ.IDFOL(J)) THEN
                JFOL = J
                OK = .TRUE.
              END IF
            END DO
            IF (OK) THEN
              IF (NBC(ID).EQ.0) THEN
c                print*,jfol,nb(jfol),t
                NB(JFOL)=NB(JFOL)+1
                TT(NB(JFOL),JFOL)=SNGL(T)
                TA(NB(JFOL),JFOL)=SNGL(A)
                TE(NB(JFOL),JFOL)=SNGL(E)
                TINC(NB(JFOL),JFOL)=SNGL(INC*DR)
                TCAPOM(NB(JFOL),JFOL)=SNGL(CAPOM*DR)
                TOMEGA(NB(JFOL),JFOL)=SNGL(OMEGA*DR)
                TCAPM(NB(JFOL),JFOL)=SNGL(CAPM*DR)
              END IF
              NBC(JFOL)=NBC(JFOL)+1
              IF (NBC(JFOL).EQ.NCG) NBC(JFOL)=0
            END IF
          END DO

          IF ((NTP.GT.0).AND.(NBC(0).EQ.0)) THEN
            NB(0)=NB(0)+1
            TT(NB(0),0)=SNGL(T)
            AMOY = 0.0D0
            EMOY = 0.0D0
            KXMOY = 0.0d0
            KYMOY = 0.0d0
            KZMOY = 0.0d0
            OMOY = 0.0D0
            DO I = 1,NTP
              AMOY = AMOY+TTA(I)
              EMOY = EMOY+TTE(I)
              KXMOY = KXMOY+SIN(TTINC(I))*SIN(TTCAPOM(I))
              KYMOY = KYMOY-SIN(TTINC(I))*COS(TTCAPOM(I))
              KZMOY = KZMOY+COS(TTINC(I))
            END DO
            AMOY = AMOY/NTP
            EMOY = EMOY/NTP
            KXMOY = KXMOY/NTP
            KYMOY = KYMOY/NTP
            KZMOY = KZMOY/NTP
            IMOY = ACOS(KZMOY)
            OMOY = ATAN2(KXMOY,-KYMOY)
            IRMOY = 0.0D0
            ORMOY = 0.0D0
            DO I = 1,NTP
C              IRMOY = IRMOY+ACOS(
C     &          SIN(TTINC(I))*SIN(IMOY)*COS(TTCAPOM(I)-OMOY)
C     &                     +COS(TTINC(I))*COS(IMOY))
              KX = SIN(TTINC(I))*SIN(TTCAPOM(I))
              KY = -SIN(TTINC(I))*COS(TTCAPOM(I))
              KZ = COS(TTINC(I))
              KXR = KX+KX*KZMOY-KX*KXMOY*KXMOY-KXMOY*KYMOY*KY
     &                          -KXMOY*KZ-KXMOY*KZ*KZMOY
              KYR = KXMOY*KYMOY*KX-KY-KY*KZMOY+KY*KYMOY*KYMOY
     &                          +KYMOY*KZ+KYMOY*KZ*KZMOY
              KZR = KX*KXMOY+KY*KYMOY+KZ*KZMOY
              IRMOY = IRMOY+ACOS(KZR)
              ORMOY = ORMOY+ATAN2(KXR,KYR)
c              KXR = (TTX(I)+TTX(I)*KZMOY-TTX(I)*KXMOY*KXMOY
c     &              -KXMOY*KYMOY*TTY(I)-KXMOY*TTZ(I)
c     &              -KXMOY*TTZ(I)*KZMOY)/(1.0d0+KZMOY)
c              KYR = -(KXMOY*KYMOY*TTX(I)+TTY(I)+TTY(I)*KZMOY
c     &              -TTY(I)*KYMOY*KYMOY-KYMOY*TTZ(I)
c     &              -KYMOY*TTZ(I)*KZMOY)/(1.0d0+KZMOY)
c              KZR = TTX(I)*KXMOY+TTY(I)*KYMOY+TTZ(I)*KZMOY
c              TTR(I) = SQRT(KXR*KXR+KYR*KYR)
              TTX(I) = ACOS(KZR)
              TTY(I) = ATAN2(KXR,KYR)
c              TTZ(I) = KZR/TTR(I)
            END DO
            IRMOY = IRMOY/NTP
            ORMOY = ORMOY/NTP
            SIGI = 0.0D0
            SIGO = 0.0D0
            SIGA = 0.0d0
            DO I = 1,NTP
              KX = SIN(TTINC(I))*SIN(TTCAPOM(I))
              KY = -SIN(TTINC(I))*COS(TTCAPOM(I))
              KZ = COS(TTINC(I))
              KXR = KX+KX*KZMOY-KX*KXMOY*KXMOY-KXMOY*KYMOY*KY
     &                          -KXMOY*KZ-KXMOY*KZ*KZMOY
              KYR = KXMOY*KYMOY*KX-KY-KY*KZMOY+KY*KYMOY*KYMOY
     &                          +KYMOY*KZ+KYMOY*KZ*KZMOY
              KZR = KX*KXMOY+KY*KYMOY+KZ*KZMOY
              IR2 = ACOS(KZR)-IRMOY
              SIGI = SIGI+IR2*IR2
              OR2 = ATAN2(KXR,KYR)-ORMOY
              SIGO = SIGO+OR2*OR2
            END DO
            SIGI = SQRT(SIGI/NTP)
            SIGO = SQRT(SIGO/NTP)
            RHO = 0.0D0
            RHA = 0.d0
            DO I = 1,NTP
              KX = SIN(TTINC(I))*SIN(TTCAPOM(I))
              KY = -SIN(TTINC(I))*COS(TTCAPOM(I))
              KZ = COS(TTINC(I))
              KXR = KX+KX*KZMOY-KX*KXMOY*KXMOY-KXMOY*KYMOY*KY
     &                          -KXMOY*KZ-KXMOY*KZ*KZMOY
              KYR = KXMOY*KYMOY*KX-KY-KY*KZMOY+KY*KYMOY*KYMOY
     &                          +KYMOY*KZ+KYMOY*KZ*KZMOY
              KZR = KX*KXMOY+KY*KYMOY+KZ*KZMOY
              RHO = RHO+((ACOS(KZR)-IRMOY)/SIGI)
     &                       *((ATAN2(KXR,KYR)-ORMOY)/SIGO)
            END DO
            RHO = RHO/NTP

c            PPP(1)=4.0d0
c            ppp(2)=1d-9
c            ppp(3)=1d-9
c            OK = .FALSE.
c            do i=-10,60
c              CALL SOLXY(NTP,TTA,TTY,TTX,0.2*i,PPP(2),PPP(3),
c     &                                CHI2,DQA,DQAA)
c              print*,0.2*i,chi2
c            end do
c            print*,'---'
c            CALL QQMIN(NTP,TTA,TTY,TTX,PPP(1),PPP(2),PPP(3),CHI2)


c            I=0
c
c            DO WHILE(.NOT.OK)
c              CALL SOLXY(NTP,TTA,TTY,TTX,PPP(1),PPP(2),PPP(3),
c     &                                    CHI2,DQA,DQAA)
c              I=I+1
c              CALL SOLXY(NTP,TTA,TTY,TTX,PPP(1)+1d-6,PPP(2),PPP(3),
c     &                                    kxr,kyr,kzr)
c              OK = (((ABS(DQA/DQAA/PPP(1)).LT.1d-9)
c     &                .and.(abs(dqa).lt.1d-4)).OR.(i.gt.100))
c              print*,'chi2',sngl(CHI2),'dqa',sngl(dqa),
c     &              'alpha',sngl(ppp(1)),sngl((kxr-chi2)/1d-6)
cc     & sngl(sqrt(PPP(2)*PPP(2)+PPP(3)*PPP(3)))
c              IF (DQA.NE.0.0d0) PPP(1)=PPP(1)-DQA/DQAA
c              IF (PPP(1).lt.-1.) PPP(1)=-1.
c            end do
c            print*,sngl(ppp(1)),
c     &       sngl(sqrt(PPP(2)*PPP(2)+PPP(3)*PPP(3))*
c     &                        (250.0d0**(PPP(1))))
C
            TA(NB(0),0)=SNGL(RHO)
            TE(NB(0),0)=SNGL(EMOY)
            TINC(NB(0),0)=SNGL(IMOY*DR)
            TCAPOM(NB(0),0)=SNGL(OMOY*DR)
            TOMEGA(NB(0),0)=SNGL(IRMOY*DR)
            TCAPM(NB(0),0)=REAL(NTP)
          END IF
          NBC(0)=NBC(0)+1
          IF (NBC(0).EQ.NCG) NBC(0)=0
        END DO MAIN_LOOP
        close(10)
        CLOSE(IU)
        NBMAX=MAXVAL(NB)
        print*,'nbmax=',nbmax,nb
        ALLOCATE(DATA4(NBMAX,7,NBOD+NPFOL))
        DATA4=0.
        DO I=-2,-NBOD,-1
          DATA4(1:NB(I),1,-I-1)=TT(1:NB(I),I)
          DATA4(1:NB(I),2,-I-1)=TA(1:NB(I),I)
          DATA4(1:NB(I),3,-I-1)=TE(1:NB(I),I)
          DATA4(1:NB(I),4,-I-1)=TINC(1:NB(I),I)
          DATA4(1:NB(I),5,-I-1)=TCAPOM(1:NB(I),I)
          DATA4(1:NB(I),6,-I-1)=TOMEGA(1:NB(I),I)
          DATA4(1:NB(I),7,-I-1)=TCAPM(1:NB(I),I)
        END DO
        DATA4(1:NB(0),1,NBOD)=TT(1:NB(0),0)
        DATA4(1:NB(0),2,NBOD)=TA(1:NB(0),0)
        DATA4(1:NB(0),3,NBOD)=TE(1:NB(0),0)
        DATA4(1:NB(0),4,NBOD)=TINC(1:NB(0),0)
        DATA4(1:NB(0),5,NBOD)=TCAPOM(1:NB(0),0)
        DATA4(1:NB(0),6,NBOD)=TOMEGA(1:NB(0),0)
        DATA4(1:NB(0),7,NBOD)=TCAPM(1:NB(0),0)
        DO I=1,NPFOL
c          print*,I,sngl(ta(1:nb(i),i))
          DATA4(1:NB(I),1,NBOD+I)=TT(1:NB(I),I)
          DATA4(1:NB(I),2,NBOD+I)=TA(1:NB(I),I)
          DATA4(1:NB(I),3,NBOD+I)=TE(1:NB(I),I)
          DATA4(1:NB(I),4,NBOD+I)=TINC(1:NB(I),I)
          DATA4(1:NB(I),5,NBOD+I)=TCAPOM(1:NB(I),I)
          DATA4(1:NB(I),6,NBOD+I)=TOMEGA(1:NB(I),I)
          DATA4(1:NB(I),7,NBOD+I)=TCAPM(1:NB(I),I)
        END DO

c        print*,DATA4(1:NB(0),7,NBOD)
        DEV = TRIM(DIRO)//'/'//followfile

c        CALL IMWRITE(DEV,DATA4,NBMAX,DBLE(NBOD),0.,0.,7,
c     &                    0.,0.,0.,NBOD,0.,0.,0.,
c     &                    SNGL(-10000.),SNGL(0.001),0.,0.,0.,ERROR)
        OPEN(18,FILE=DEV,STATUS='UNKNOWN')
        WRITE(18,*)NBMAX,7,NBOD+NPFOL
        WRITE(18,*)DATA4
        CLOSE(18)

c        OPEN(19,file='fomb_meca1.dat',status='unknown')
c        print*,npfol
c        DO I=1,NPFOL
c           print*,nb(i)
c           do j=1,nb(i)
c         write(19,199)DATA4(j,1,NBOD+I),DATA4(j,2,NBOD+I),
c     &   DATA4(j,3,NBOD+I),DATA4(j,5,NBOD+I)+DATA4(j,6,NBOD+I)
c           end do
c        END DO
c        close(19)


        DEALLOCATE(DATA4)


        END



C--------------------------------------------------------------------
C               Keplerian orbit stuff
C--------------------------------------------------------------------
C
c       SUBROUTINE KEPLER_STUFF(GM,A,E,INC,CAPOM,OMEGA,CAPM,
c    &                   R,CX,CY,CZ,VX,VY,VZ,EX,EY,EZ,EPX,EPY,EPZ)
        SUBROUTINE KEPLER_STUFF(A,E,INC,CAPOM,OMEGA,CAPM,CX,CY,CZ)

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
c        VRX = -SQRT(GM/A)*A/R*SIN(U)
c        VRY = SQRT(GM/A)*A/R*COS(U)
        CX = EX*RX+EPX*RY
        CY = EY*RX+EPY*RY
        CZ = EZ*RX+EPZ*RY
c        VX = EX*VRX+EPX*VRY
c        VY = EY*VRX+EPY*VRY
c        VZ = EZ*VRX+EPZ*VRY
        END

C
C -----------------------------------------------------------------------------
C       XY linear solution
C -----------------------------------------------------------------------------
C

        SUBROUTINE SOLXYB(NDAT,A,OMEGA,INC,ALPHA,X0,Y0,Q,DQA,DQAA)

        IMPLICIT NONE

        INTEGER*4       NDAT            ! Dimension of parameter space

        REAL*8          A(NDAT),OMEGA(NDAT),INC(NDAT)  ! DATA
        REAL*8          EA,SA2,SASOCISI,SACOCISI,SA2SI2CO2,SA2SI2COSO
        REAL*8          CO,SO,CI,SI,LA,RA,X0,Y0,ALPHA
        REAL*8          SA2SI2,SLNA2,SLNACOCISI,SLNASOCISI,SLNA2SI2
        REAL*8          SLNA2SI2CO2,SLNA2SI2COSO
        REAL*8          DQA,DQAA,Q,DENO,DX0,DY0,DQAX,DQAY
        INTEGER*4       I

        SA2=0.
        SASOCISI=0.
        SACOCISI=0.
        SA2SI2CO2=0.
        SA2SI2COSO=0.
        SA2SI2=0.
        SLNA2=0.
        SLNACOCISI=0.
        SLNASOCISI=0.
        SLNA2SI2=0.
        SLNA2SI2CO2=0.
        SLNA2SI2COSO=0.
c        EA=EXP(ALPHA)
        EA = ALPHA
        DO I=1,NDAT
          RA=A(I)**EA
          SO=SIN(OMEGA(I))
          CO=COS(OMEGA(I))
          SI=SIN(INC(I))
          CI=COS(INC(I))
          LA=LOG(A(I))
          SA2=SA2+RA*RA
          SASOCISI=SASOCISI+RA*SO*CI*SI
          SACOCISI=SACOCISI+RA*CO*CI*SI
          SA2SI2CO2=SA2SI2CO2+RA*RA*SI*SI*CO*CO
          SA2SI2COSO=SA2SI2COSO+RA*RA*SI*SI*CO*SO
          SA2SI2=SA2SI2+RA*RA*SI*SI
          SLNA2=SLNA2+2.*RA*RA*LA
          SLNACOCISI=SLNACOCISI+RA*LA*CO*CI*SI
          SLNASOCISI=SLNASOCISI+RA*LA*SO*CI*SI
          SLNA2SI2=SLNA2SI2+2.*RA*RA*LA*SI*SI
          SLNA2SI2CO2=SLNA2SI2CO2+2.*RA*RA*LA*SI*SI*CO*CO
          SLNA2SI2COSO=SLNA2SI2COSO+2.*RA*RA*LA*SI*SI*CO*SO
        END DO
        DENO=-SA2*SA2+SA2*SA2SI2-SA2SI2CO2*SA2SI2
     &                  +SA2SI2CO2*SA2SI2CO2+SA2SI2COSO*SA2SI2COSO
        Y0 =-(SA2*SACOCISI-SA2SI2*SACOCISI+SA2SI2CO2*SACOCISI
     &                  +SASOCISI*SA2SI2COSO)/DENO
        X0=-(-SA2*SASOCISI+SA2SI2CO2*SASOCISI
     &                  -SACOCISI*SA2SI2COSO)/DENO
        DY0=-(SLNA2*SACOCISI+SA2*SLNACOCISI-SLNA2SI2*SACOCISI
     &        -SA2SI2*SLNACOCISI+SLNA2SI2CO2*SACOCISI
     &        +SA2SI2CO2*SLNACOCISI+SLNASOCISI*SA2SI2COSO
     &        +SASOCISI*SLNA2SI2COSO)/DENO+(SA2*SACOCISI
     &        -SA2SI2*SACOCISI+SA2SI2CO2*SACOCISI+SASOCISI*SA2SI2COSO)*
     &        (-2.*SA2*SLNA2+SLNA2*SA2SI2+SA2*SLNA2SI2
     &        -SLNA2SI2CO2*SA2SI2-SA2SI2CO2*SLNA2SI2
     &        +2.*SA2SI2CO2*SLNA2SI2CO2+2.*SA2SI2COSO+SLNA2SI2COSO)
     &                                                     /DENO/DENO
        DX0=-(-SLNA2*SASOCISI-SA2*SLNASOCISI+SLNA2SI2CO2*SASOCISI
     &        +SA2SI2CO2*SLNASOCISI-SLNACOCISI*SA2SI2COSO
     &        -SACOCISI*SLNA2SI2COSO)/DENO+(-SA2*SASOCISI
     &        +SA2SI2CO2*SASOCISI-SACOCISI*SA2SI2COSO)*
     &        (-2.*SA2*SLNA2+SLNA2*SA2SI2+SA2*SLNA2SI2
     &        -SLNA2SI2CO2*SA2SI2-SA2SI2CO2*SLNA2SI2
     &        +2.*SA2SI2CO2*SLNA2SI2CO2+2.*SA2SI2COSO+SLNA2SI2COSO)
     &                                                     /DENO/DENO
        DQA=0.
        DQAA=0.
        DQAX=0.
        DQAY=0.
        Q=0.
        DO I=1,NDAT
          SO=SIN(OMEGA(I))
          CO=COS(OMEGA(I))
          SI=SIN(INC(I))
          CI=COS(INC(I))
          LA=LOG(A(I))
          RA=A(I)**EA
          Q=Q+0.5*RA*RA*(X0*X0+Y0*Y0)+0.5*SI*SI*(1.-RA*RA*Y0*Y0*CO*CO
     &       -RA*RA*X0*X0*SO*SO+2.*RA*RA*X0*Y0*CO*SO)
     &       +RA*CI*SI*(-Y0*CO+X0*SO)
          DQA=DQA+RA*RA*LA*(X0*X0+Y0*Y0)+SI*SI*LA*RA*RA*(-Y0*Y0*CO*CO
     &       -X0*X0*SO*SO+2.*X0*Y0*CO*SO)+RA*LA*CI*SI*(-Y0*CO+X0*SO)
          DQAX=DQAX+2.*RA*RA*LA*X0+SI*SI*LA*RA*RA*(
     &       -2.*X0*SO*SO+2.*Y0*CO*SO)+RA*LA*CI*SI*SO
          DQAY=DQAY+2.*RA*RA*LA*Y0+SI*SI*LA*RA*RA*(-2.*Y0*CO*CO
     &       +2.*Y0*CO*SO)-RA*LA*CI*SI*CO
          DQAA=DQAA+2.*RA*RA*LA*LA*(X0*X0+Y0*Y0)+2.*SI*SI*LA*LA*RA*RA*(
     &       -Y0*Y0*CO*CO-X0*X0*SO*SO+2.*X0*Y0*CO*SO)
     &       +RA*LA*LA*CI*SI*(-Y0*CO+X0*SO)
        END DO
c        print*,dqaa,dqax*dx0,dqay*dy0
c        DQAA=(DQAA+DQAX*DX0+DQAY*DY0)*EA*EA+DQA*EA
c        DQA=DQA*EA
        DQAA=(DQAA+DQAX*DX0+DQAY*DY0)
        END


C
C------------------------------------------------------------------------------
C     Procedure qui, partant d'un point P et d'un vecteur u,
C     Recherche de scalaire l (ici x(0)) minimisant QDIFF(P+l*u)
C     (Minimisation de QDIFF a partir de P dans la direction de u )
C     Remplace ensuite P par l'endroit correspondant au minimum
C------------------------------------------------------------------------------
C
        SUBROUTINE QQMIN(NDAT,A,OMEGA,INC,XTRY,X0,Y0,YTRY)

        IMPLICIT NONE

        INTEGER*4       NDAT
	REAL*8		XTRY,YTRY,	! Point d'essai
     &                  A(NDAT),OMEGA(NDAT),INC(NDAT),
     &                  D1,D2,DENO,     ! Intermediaires
     &                  DD,             ! Difference
     &                  QDIFF,          ! Operateur a minimiser
     &                  NOR,            ! Norme
     &                  PI,X0,Y0,       ! Pi
     &                  X(-1:1),Y(-1:1),! Valeurs aux 3 points d'essai
     &                  EPS             ! Precision
        INTEGER*4       ITS,            ! Nombre d'iterations
     &                  I,J,            !
     &                  DI              ! Sens du deplacement
        LOGICAL         OK
        character*1     car
        PARAMETER(      PI=3.14159265358979323846,EPS=1d-7)

C                         Initialisation

        DO I=-1,1
          X(I)=(DFLOAT(I)+1.0d0)*0.1
          CALL SOLXY(NDAT,A,OMEGA,INC,X(I),X0,Y0,Y(I))
        END DO
c            print*,'deb...',sngl(x(-1)),sngl(x(0)),sngl(x(1))
c            print*,'----..',sngl(y(-1)-y(0)),sngl(y(0)),sngl(y(1)-y(0))
c        read(5,'(a)')car
C                          Recherche de l'encadrement initial du minimum

        OK=.FALSE.
        ITS=0.
        DO WHILE(.NOT.OK)
          OK=.TRUE.
          IF (Y(1).LT.Y(0)) THEN
            DI=1
            OK=.FALSE.
          END IF
          IF (Y(-1).LT.Y(0)) THEN
            DI=-1
            OK=.FALSE.
          END IF
          IF (.NOT.OK) THEN
            X(-DI)=X(0)
            Y(-DI)=Y(0)
            X(0)=X(DI)
            Y(0)=Y(DI)
            X(DI)=X(DI)+DI*0.1
c	    DO J=1,N
c              V(J)=P(J)+X(DI)*U(J)*NOR*0.1*EXP(DFLOAT(ITS))
c            END DO
            CALL SOLXY(NDAT,A,OMEGA,INC,X(DI),X0,Y0,Y(DI))
c            print*,'depl..',sngl(x(-1)),sngl(x(0)),sngl(x(1))
c            print*,'----..',sngl(y(-1)-y(0)),sngl(y(0)),sngl(y(1)-y(0))
          END IF
          ITS=ITS+1
        END DO
c            print*,'depl fini'
c            read(5,'(a)')car
C                Recherche du minimum par paraboles inverses
        ITS=0
        OK=.FALSE.
        DO WHILE(.NOT.OK)
          D1=X(0)-X(-1)
          D2=X(0)-X(1)
          IF ((ABS(D1).GT.EPS).AND.(ABS(D2).GT.EPS)) THEN
            D1=D1*(Y(0)-Y(1))
            D2=D2*(Y(0)-Y(-1))
            DENO=2.0*(D1-D2)
            D1=D1*(X(0)-X(-1))
            D2=D2*(X(0)-X(1))
            IF (DENO.NE.0) THEN
              XTRY=X(0)-(D1-D2)/DENO
c              print*,'par !'
            ELSE
              XTRY=X(-1)+X(1)-X(0)
c              print*,'lin !'
            END IF
          ELSE
            XTRY=X(-1)+X(1)-X(0)
c            print*,'lin !'
          END IF
          DD=ABS(X(0)-XTRY)
c          DO J=1,N
c            V(J)=P(J)+XTRY*U(J)*NOR*0.1
c          END DO
c          YTRY=QDIFF(V(1),V(2),V(3),V(4),V(5))
          CALL SOLXY(NDAT,A,OMEGA,INC,XTRY,X0,Y0,YTRY)
          DI=0
          IF (XTRY.LT.X(0)) DI=1
          IF (XTRY.GT.X(0)) DI=-1
          IF (YTRY.LT.Y(0)) THEN
            X(DI)=X(0)
            Y(DI)=Y(0)
            X(0)=XTRY
            Y(0)=YTRY
          ELSE
            X(-DI)=XTRY
            Y(-DI)=YTRY
          END IF
          ITS=ITS+1
          OK=((DD.LT.EPS*ABS(X(0))).OR.(ITS.GT.200))
c           print*,sngl(x(-1)),sngl(x(0)),sngl(x(1)),
c     &         sngl(y(-1)-y(0)),sngl(y(0)),sngl(y(1)-y(0))
        END DO
c            print*,'par fini'
c            read(5,'(a)')car
C                      Remplacement de P par l'endroit trouve
        XTRY=X(0)
        CALL SOLXY(NDAT,A,OMEGA,INC,XTRY,X0,Y0,YTRY)
c        print*,(sngl(p(j)),j=1,N)
c        print*,'        ----> ',y(0),' <----'
        END
C
C -----------------------------------------------------------------------------
C       XY linear solution
C -----------------------------------------------------------------------------
C

        SUBROUTINE SOLXY(NDAT,A,OMEGA,INC,ALPHA,X0,Y0,Q)

        IMPLICIT NONE

        INTEGER*4       NDAT            ! Dimension of parameter space

        REAL*8          A(NDAT),OMEGA(NDAT),INC(NDAT)  ! DATA
        REAL*8          EA,SA2,SASOCISI,SACOCISI,SA2SI2CO2,SA2SI2COSO
        REAL*8          CO,SO,CI,SI,LA,RA,X0,Y0,ALPHA
        REAL*8          SA2SI2,A2
        REAL*8          Q,DENO,DX0,DY0,DQAX,DQAY
        INTEGER*4       I

        SA2=0.
        SASOCISI=0.
        SACOCISI=0.
        SA2SI2CO2=0.
        SA2SI2COSO=0.
        SA2SI2=0.
c        EA=EXP(ALPHA)
        EA = ALPHA
        DO I=1,NDAT
          A2=A(I)*A(I)
          RA=A(I)**EA
          SO=SIN(OMEGA(I))
          CO=COS(OMEGA(I))
          SI=SIN(INC(I))
          CI=COS(INC(I))
          SA2=SA2+RA*RA*A2
          SASOCISI=SASOCISI+RA*SO*CI*SI*A2
          SACOCISI=SACOCISI+RA*CO*CI*SI*A2
          SA2SI2CO2=SA2SI2CO2+RA*RA*SI*SI*CO*CO*A2
          SA2SI2COSO=SA2SI2COSO+RA*RA*SI*SI*CO*SO*A2
          SA2SI2=SA2SI2+RA*RA*SI*SI*A2
        END DO
        DENO=-SA2*SA2+SA2*SA2SI2-SA2SI2CO2*SA2SI2
     &                  +SA2SI2CO2*SA2SI2CO2+SA2SI2COSO*SA2SI2COSO
        Y0 =-(SA2*SACOCISI-SA2SI2*SACOCISI+SA2SI2CO2*SACOCISI
     &                  +SASOCISI*SA2SI2COSO)/DENO
        X0=-(-SA2*SASOCISI+SA2SI2CO2*SASOCISI
     &                  -SACOCISI*SA2SI2COSO)/DENO

        Q=0.
        DO I=1,NDAT
          SO=SIN(OMEGA(I))
          CO=COS(OMEGA(I))
          SI=SIN(INC(I))
          CI=COS(INC(I))
          RA=A(I)**EA
          A2=A(I)*A(I)
          Q=Q+A2*(0.5*RA*RA*(X0*X0+Y0*Y0)+0.5*SI*SI*(1.
     &        -RA*RA*Y0*Y0*CO*CO
     &       -RA*RA*X0*X0*SO*SO+2.*RA*RA*X0*Y0*CO*SO)
     &       +RA*CI*SI*(-Y0*CO+X0*SO))
        END DO
        END
