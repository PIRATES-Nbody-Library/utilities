C
C-------------------------------------------------------------------------
C                   Main program
C-------------------------------------------------------------------------
C
        PROGRAM GEN_MULTI_REPRISE

        include 'swift.inc'

        INTEGER*4 NFILE,I

        CHARACTER*256 GOFILE,MVSFILE,PARAMFILE,PLFILE
        CHARACTER*256 SIMUNAME,TPFILE

5000    FORMAT('mvs_reprise_',I2.2,'.in')
6000    FORMAT('go_rmvs3_reprise_',I2.2,'.sh')
7001    FORMAT('_',I2.2,'/dump_param.dat')
7002    FORMAT('_',I2.2,'/dump_pl.dat')
7003    FORMAT('_',I2.2,'/dump_tp.dat')


c ........entree des parametres

        READ (*,'(a)')SIMUNAME
        READ(*,*) NFILE

        DO I=1,NFILE

           WRITE(MVSFILE,5000)I
           WRITE(PARAMFILE,7001)I
           WRITE(PLFILE,7002)I
           WRITE(TPFILE,7003)I
           OPEN(31,FILE=MVSFILE,STATUS='UNKNOWN')
           WRITE(31,'(a)')'__TMP_DIR__/run_'
     &     //TRIM(SIMUNAME)//TRIM(PARAMFILE)
           WRITE(31,'(a)')'__TMP_DIR__/run_'
     &     //TRIM(SIMUNAME)//TRIM(PLFILE)
           WRITE(31,'(a)')'__TMP_DIR__/run_'
     &     //TRIM(SIMUNAME)//TRIM(TPFILE)
           WRITE(31,'(a)')'1d-8'
           CLOSE(31)

           WRITE(GOFILE,6000)I
           OPEN(31,FILE=GOFILE,STATUS='UNKNOWN')
           WRITE(31,'(a)')'#!/bin/bash'
           WRITE(31,'(a)')'#SBATCH --account=__SBATCH_ACCOUNT__'
           WRITE(31,'(a)')'#SBATCH --job-name=go_rmvs'
           WRITE(31,'(a)')'#SBATCH --mail-type=ALL'
           WRITE(31,'(a)')'#SBATCH --mail-user=__SBATCH_EMAIL__'
     &
           WRITE(31,'(a)')'#SBATCH --partition=__SBATCH_PARTITION__'
           WRITE(31,'(a)')'#SBATCH --ntasks=1'
           WRITE(31,'(a)')'#SBATCH --nodes=1'
           WRITE(31,'(a)')'#SBATCH --mem=__SBATCH_MEM__'
           WRITE(31,'(a)')'#SBATCH --time=__SBATCH_TIME__'
           WRITE(31,'(a)')' '
           WRITE(31,'(a)')'__SWIFT_DIR__/main'
     &     //'/swift_rmvs3 < __MULTI_DIR__/simu_'
     &     //TRIM(SIMUNAME)//'/'//TRIM(MVSFILE)
           CLOSE(31)
           CALL SYSTEM('chmod ogu+x '//TRIM(GOFILE))

            END DO

        end program gen_multi_reprise
