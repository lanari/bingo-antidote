!-----Version: 15.03.2018
!               *************
!               * THERIAK D *
!               *************
!
!	  THERIAKD calls THERIAK(Ref.1) as a function with two line arguments.
!     First line argument is an input file, which contains information about
!     1. database, 2. temperature, 3. pressure and 4. command.
!     Second line argument is the file of bulk composition(e.g. THERIN).
!     version: 15.03.2018
!     ========================================
!     Example - Get total density of solids:
!	   
!	  infile (input file):
!     ----------
!     JUN92.bs
!     400
!     4000
!     DENSOL
!     ----------
!
!     INPUT unix: theriakd infile THERIN
!
!     OUTPUT unix: 2.81
!     =========================================
!
!     Add-on THERIAK D 
!     written by Erik Duesterhoeft & Christian de Capitani
!
!     Any questions relevant to the addon should be sent to:
!          Erik Duesterhoeft
!          Christian-Albrechts-Universität zu Kiel
!          Institut für Geowissenschaften
!          Ludewig-Meyn-Str. 10
!          24118 Kiel
!          GERMANY
!
!
!     Ref.1 :
!
!     THERIAK
!     Program written by Christian de Capitani
!     at the Department of Geological Sciences,
!     University of British Columbia, Vancouver, B.C. Canada
!     (May 1984 - Sept 1987)
!
!     revision: April 1987
!     minor changes: December 1987
!     major revision: July 1993
!     revisions: October 2002, July 2004, February 2005, May 2006, March 2007,
!                December 2007, May 2008, August 2009
!
!     for details of algorithm see:
!     de Capitani C. and Brown T.H. : The computation of chemical
!     equilibrium in complex systems containing non-ideal solutions.
!     Geochim. Cosmochim. Acta 51(1987):2639-2652
!
!     Any suggestions, complaints or comments are greatly appreciated
!     by the author and should be sent to:
!          Christian de Capitani
!          Mineralogisch-Petrographisches Institut
!          Universitaet Basel
!          Bernoullistrasse 30
!          CH-4056 BASEL
!
!
!
!
      IMPLICIT NONE
      INCLUDE 'theriak.cmn'
      include 'files.cmn'
!*****
!
!-----END OF COMMON VARIABLES
      LOGICAL*4 LPFILE,SIM,LOLO, IPR, CONCH
      LOGICAL*4 ISASOLID(COMAX)
      COMMON /LOSOLID/ ISASOLID
      INTEGER*4 I001,I002,I,I1,COMAY,NPTS,EINS,ALLES,ierr,j,BIN1, &
      LARG, status_read, LPP, NSS, NSP, II
      REAL*8 FF, R001, WASSER, RHOH2O, M_RHO, MOL, XMOL, WPRO, VPRO
      REAL*8 KDGRTCPX, KDGRTOPX, KDOPXCPX, KDOLGRT, KDOLOPX, &
      KDOLCPX, KDCTDCAR, KDGRTBT, KDCTDCHL
      REAL*8 M_K, M_H, M_O, M_SI, M_AL, M_FE, M_MG, M_CA, M_NA, &
      M_WASSER
      REAL*8 HTOT, STOT, CPTOT, KT, KS, EXPAN
      CHARACTER*16 CH16
      CHARACTER*80 KOMMENTAR,INFILE
      CHARACTER*500 CH001,CH002,SYREC,CHIN(2),ZEITSTRING


      progname='THERIAK_D'
      vers='15.03.2018'
      task='"Computation of equilibrium assemblages at given PT"'
      EINS=1
      ierr=9 					! is set to 9 to suppress output
      
      !call initialize('$THERIAK-FILES',ierr)
!*******************************************************************

        ierr=9
        IF (ierr.EQ.9) IPR=.FALSE.
!-----  first find working- and program-directory
!
        CALL GetEnvVar(tpath, itpath, ierr, dir, ext(job), os)
        ierr=9
        CALL getwork(wpath, iwpath, ierr, dir)
!------------------
!       eduester:
!       FAKE: Open UNIT=ini (= theriak.ini)     
!------------------

        clear=0

        ghelp=0

        batch=1

        LO1MAX=300
        EQUALX=1.0D-2
        TEST=1.0D-9
        DXMIN=1.0D-9 
        DXSCAN=1.0D0
        DXSTAR=1.0D-4
        STPSTA=100
        STPMAX=25
        GCMAX=500
        LPP=30
!*******************************************************************
! Get input- und dat-file
      if(ierr.ne.0) STOP
      REFAS=.FALSE.
!
!-----
      LARG=0
      IERR=0 
      DO I=1,5
      CALL GetLineArgs (I,LARGUM(I),IERR)
      IF(IERR.NE.0.OR.LARGUM(I).EQ.' ')THEN
      GOTO 399
      ELSE
      LARG=LARG+1
      END IF
      END DO
  399 CONTINUE
      IFNR=5
!-----
      IF (LARG.GT.0) THEN
      INFILE=LARGUM(1)
      CALL LABLA(INFILE,I1)
      IFNR=39
      OPEN (UNIT=IFNR,FILE=INFILE(1:I1),STATUS='UNKNOWN')
      END IF
!-----
      IF (LARG.GT.1) THEN
      filename(dat)=LARGUM(2)
      CALL LABLA(filename(dat),fnl(dat))
      END IF

      CHIN(1)=' '
      CHIN(2)='no'

         
!*******************************************************************
! Read database-name out of infile (e.g. JUN92.bs)
   
      READ (IFNR,FMT='(A500)') CH001
      CH002=CH001
      CALL TAXI(CH002,DBNAME)
      CALL LABLA(DBNAME,I001)
      !===== for theriak3
      filename(dbs)=DBNAME

      
!*******************************************************************
!
!------------------
!     open UNIT=dbs
!------------------
      j=dbs
      line=DBNAME
      path=wpath
      akzess=' '
      state='old'
      call openfile(j,ierr)
      if(ierr.ne.0) STOP
!------------------
!     open UNIT=dat
!------------------
      j=dat
      line=filename(j)(1:fnl(j))//ext(j)
      path=wpath
      akzess=' '
      state='old'
      call openfile(j,ierr)
      if(ierr.ne.0) STOP
!*****

      COMAY=COMAX
      CALL PROREAD(SYREC)
      
!-- pick nach proread!
      NPICK=0
      DO I=1,COMAX
       CALL TAXI(CH002,CH16)
       IF (CH16.NE.' ') THEN
       NPICK=NPICK+1
       PICK(NPICK)=CH16
       ELSE
       GOTO 5
       END IF
      END DO
    5 CONTINUE
    
!*****
!
!*******************************************************************
! Define type of Calculations
! (here the value is fixed to 'no') 

      CH001='no'

!*******************************************************************
!---
!---
      SIM=.FALSE.
      LOLO=.FALSE.
      CONCH=.FALSE.
      LPFILE=.FALSE.           
!-----
!-----READ  PRTCOD, FORMUL AND USE FROM SYREC
!-----SET UP FIRST NUN COLUMNS OF MATRIX
      CALL GELI(SYREC,FF)

      PRTCOD=0  ! is set to 0 to suppress long output
      
      IF (PRTCOD.EQ.0) TEST=DABS(TEST)
      DO 650,I=1,11
  650 PRTLOG(I)=.FALSE.
  

      CALL FIBLA(SYREC,I1)
      IF (I1.EQ.0) I1=1
      CH002=' '//SYREC(I1:)
      
!-----------------------------------------------------------------
! read dat-file (THERIN)
      CALL TAXI(SYREC,FORMUL)
      NBUL=1
      BULINE(1)=FORMUL
      CALL TAXI(SYREC,USE)	! important read USE-command
      CALL TAXI(SYREC,KOMMENTAR)
!---- define first three COMINS
      CALL LABLA(CHIN(1),I002)
      NCOMIN=3
      COMINS(1)=KOMMENTAR
      COMINS(2)='theriak version: '//vers(1:10)
      COMINS(3)='database: '//CHIN(1)(1:I002)
!+++++
      CALL CHEMIE(COMAY,NC,OXYDE,OXANZ,FORMUL,CHEM)
      CALL LABLA(USE,LUSE)
!-----------------------------------------------------------------
! read database-file
      CALL DBREAD
   
!---
      IF (PEXCL.GT.0) THEN
      DO I=1,NPHA
       IF (NAME(I)(1:1).EQ.'$') THEN
       NCOMIN=NCOMIN+1
       COMINS(NCOMIN)='excl.: '//NAME(I)(2:)
       IF (NCOMIN.GE.50) GOTO 691
       END IF
      END DO
      END IF
               
      
      IF (SEXCL.GT.0) THEN
      DO I=1,NSOL

       IF (EXSOL(I)) THEN
       NCOMIN=NCOMIN+1
       COMINS(NCOMIN)='excl.: '//SOLNAM(I)
       IF (NCOMIN.GE.50) GOTO 691
       END IF

      END DO
      END IF
  691 CONTINUE
!---
!*******************************************************************
!Read Pressure and Temperature from Infile

IF (IFNR.EQ.39) THEN
		READ (UNIT=IFNR,FMT='(F9.0)') R001
  		TC = R001
  		READ (UNIT=IFNR,FMT='(F9.0)') R001
  		P = R001
END IF
!write(*,*) P,TC
!*******************************************************************
! START CALCULATION

      CALL NURVONPT	! important
      CALL CALSTR	! important
      ALLES=NSOL
!--- eduester
      IF (LPFILE.OR.PRTCOD.EQ.0) CALL GIBBS_D(EINS,ALLES) ! important for complex solution
        	
!---      	                  
      IF (.NOT.LPFILE) THEN
      	      	      	            	         	
        CALL THERIA	! important 
      	CLOSE (dat)	! important      
!---
!loop-file 
      ELSE
      	CLOSE (dat)
      	CH001=CHIN(2)
      	CALL MINILOOP(CH001)
      	GOTO 999
!-----------------------------------------------------------------
  999 CONTINUE
!-----------------------------------------------------------------

      END IF

!*******************************************************************     
!                       OUTPUT
!*******************************************************************

!Read command from Infile
READ (UNIT=IFNR,FMT='(A500)', iostat = status_read) CH001

HTOT=0
STOT=0
CPTOT=0
KT=0
KS=0
EXPAN=0

!==========

WASSER=0
RHOH2O=0

!==========
!KD (Fe-Mg)
KDGRTCPX=0
KDGRTOPX=0
KDOPXCPX=0
KDOLGRT=0
KDOLOPX=0
KDOLCPX=0
KDCTDCAR=0
KDGRTBT=0
KDCTDCHL=0

!==========
XMOL=0
MOL=0
WPRO=0
VPRO=0

!==========
M_WASSER=0
M_RHO=0

!!=============
!MINERAL-PHASEN
!X(I,J) I=1:NUN2(PHASE) J=1:NUN (ELEMENTS) --> MENGE
!CHNAME(J) -Element Name (K,Mg,Fe ...)
!
!Vorgehensweise: Bsp. K in FSP:
! Erst K suchen in CHNAME(J) um J bekommen
! Dann such ich in STHP(I) nach FSP um I zu bekommen
! Dann hab X(FSP,K)

!CODE: *XGRTSI

IF (CH001(1:1).EQ.'*') THEN
	DO I=1,NUN2
	!write(*,*) 'first',  STPHNAM(I), NN(I)
	    IF (EMCODE(I) .EQ. 0) THEN !check for solutions, 0=no; 1= yes, there are
	    
		    !IF (STPHNAM(I)(1:3).EQ.CH001(3:5)) EXIT
		    IF (SHORTNAM(I)(1:4).EQ.CH001(3:6)) EXIT
	
		ELSE
		!wenn es eine solution ist
		    !write(*,*) SOLNAM(I)
		    IF (SOLNAM(EMCODE(I))(1:4).EQ.CH001(3:6)) EXIT
	
		END IF
	
	END DO
	
	IF (CH001(1:2).EQ.'*X') THEN	
		DO J=1, NUN
			IF (CHNAME(J)(1:2).EQ.CH001(7:8)) THEN
				XMOL=X(I,J)*NN(I) 
				EXIT
			END IF
		END DO
	END IF
	
	!DO I=1,NUN2
	!	DO J=1, NUN
	!		IF (CHNAME(J).EQ.'SI') EXIT
	!	END DO
	!	IF (STPHNAM(I)(1:6).EQ.'GRT') EXIT
	!END DO	
	
	
	MOL=NN(I)
	WPRO=WTPH(I)/WTSOL*100.0D0
	VPRO=VOLPH(I)/VOLSOL*100.0D0
END IF



! DENS=WTPH(I)/VOLPH(I)

! molares Volumen = VOLM(I)
! molares Gewicht = WTM(I)

!Mol-Percent
!MOL=X(I,II)/SUMM(I)*100.0D0
!=============
!VPRO=VOLPH(I)/VOLSOL*100.0D0



!===============
! WRITE OUTPUT
	IF (status_read .EQ. 0) THEN  
		 SELECT CASE (CH001)
				CASE ('DENSOL')
				! total density of solids
					write(scr,"(F20.6)") WTSOL/VOLSOL
				CASE ('VOLSOL')
				! total volume of solids
					write(scr,"(F12.6)") VOLSOL 
					!WRITE (UNIT=scr,FMT='(''volsol= '',1PE15.8)') VOLSOL
				CASE ('WTSOL')
				! total molar weight of solids
					write(scr,"(F12.6)") WTSOL 
				CASE ('GTOT')
				! Free Gibbs Energy (System)
					write(scr,"(F20.6)") GGTOT
				CASE ('HTOT')
				! Enthalpy
					CALL DSDVTEST5(HTOT,STOT,CPTOT)
					write(scr,"(F20.6)") HTOT
				CASE ('STOT')
				! Entropy
					CALL DSDVTEST5(HTOT,STOT,CPTOT)
					write(scr,"(F20.6)") STOT
				CASE ('CPTOT')
				! Heat capacity
					CALL DSDVTEST5(HTOT,STOT,CPTOT)
					write(scr,"(F20.6)") CPTOT
				CASE ('CPRHO')
				! Heat capacity * density
					CALL DSDVTEST5(HTOT,STOT,CPTOT)
					write(scr,"(F20.6)") CPTOT*WTSOL/VOLSOL
				CASE ('KT')
				! isothermal bulk modulus (GPa)
					CALL BULKEXPAN(KT,KS,EXPAN)
					write(scr,"(F20.4)") KT*0.0001
				CASE ('KS')
				! adiabatic bulk modulus (GPa)
					CALL BULKEXPAN(KT,KS,EXPAN)
					write(scr,"(F20.4)") KS*0.0001
				CASE ('EXPAN')
				! thermal expansivity
					CALL BULKEXPAN(KT,KS,EXPAN)
					write(scr,"(F20.10)") EXPAN
				CASE ('SHEAR')
				! shear modulus (GPa)
					write(scr,"(F20.4)") SSVELO*SSVELO &
					*WTSOL/VOLSOL*1000.0D0/1000000000.0D0
				CASE ('POIS')
				! poissons ratio
					CALL BULKEXPAN(KT,KS,EXPAN)
					write(scr,"(F20.6)") (3.0D0*KS*100000.0D0 &
					-2.0D0*(SSVELO*SSVELO &
					*WTSOL/VOLSOL*1000.0D0))/(6.0D0*KS*100000.0D0 &
					+2.0D0*(SSVELO*SSVELO*WTSOL/VOLSOL*1000.0D0))	
				CASE ('SWV')
				! S-wave velocity of solids [m/s] 
					write(scr,"(F12.6)") SSVELO	
			    CASE ('PWV')
				! P-wave velocity of solids [m/s] 
					write(scr,"(F12.6)") PPVELO	
					
							 
				CASE ('H')
				! H- Hydrogen in free fluid water [mol]
					!Searching H- Hydrogen in water [mol]
					DO i=1, NUN2
					  IF (STPHNAM(I).EQ. 'STEAM') THEN 		! Stable Phase Name
						WASSER = NN(I)*2 + WASSER			! Number of Moles
						RHOH2O = WTPH(I)/VOLPH(I)			! Density of free fluid water
					  END IF
					  IF (STPHNAM(I).EQ. 'water.fluid') THEN ! Stable Phase Name
					    WASSER = NN(I)*2 + WASSER			! Number of Moles
					    RHOH2O = WTPH(I)/VOLPH(I)			! Density of free fluid water
					  END IF
	
					!IF (STPHNAM(I).EQ. 'LIQtc_h2oL') THEN 	! H2O melt component
					!WASSER = NN(I)*2 + WASSER				! Number of Moles
					!END IF
					END DO		
							 
					write(scr,"(F12.6)") WASSER
					
				CASE ('DENH2O')
				! density of free fluid H2O [g/cm3] 
					!Searching H- Hydrogen in water [mol]
					DO i=1, NUN2
					  IF (STPHNAM(I).EQ. 'STEAM') THEN 		! Stable Phase Name
						WASSER = NN(I)*2 + WASSER			! Number of Moles
						RHOH2O = WTPH(I)/VOLPH(I)			! Density of free fluid water
					  END IF
					  IF (STPHNAM(I).EQ. 'water.fluid') THEN ! Stable Phase Name
					    WASSER = NN(I)*2 + WASSER			! Number of Moles
					    RHOH2O = WTPH(I)/VOLPH(I)			! Density of free fluid water
					  END IF
					END DO	
						
					write(scr,"(F12.6)") RHOH2O
					
				CASE ('H2OSOL')
				! H2O in solids [mol] 
					write(scr,"(F12.6)") H2OSOL
				CASE ('WH2OSOL')
				! H2O in solids [g] 
					write(scr,"(F12.6)") WH2OSOL
				CASE ('WPH2OSOL')
				! H2O in solids [wt%] 
					write(scr,"(F12.6)") WH2OSOL/WTSOL*100.0D0
				
				
				!CASE ('GRT')
				! *SI ??? ! Number of moles of a phases [mol] 
				!	write(scr,"(F12.6)") MOL
				!CASE ('WGRT')
				! Weight percent of a phase
				!	write(scr,"(F12.6)") WPRO
				!CASE ('VGRT')
				! Weight percent of a phase
				!	write(scr,"(F12.6)") VPRO
										
					
				CASE ('DENMELT')
				! density of melt [g/cm3] 

					!============
					! OLD !!!
					!Searching number of moles of melt [mol] --> its now molar amount of melt					
					DO i=1, NUN2	
						IF (STPHNAM(I)(1:3).EQ. 'LIQ') THEN 	! H2O melt component
						 M_WASSER = NN(I) + M_WASSER			! Number of Moles
						 M_RHO = WTPH(I)/VOLPH(I) 				! Density of melt
						END IF
											
						!IF (STPHNAM(I).EQ. 'LIQtc_h2oL') THEN 	! H2O melt component
						!M_WASSER = NN(I) + M_WASSER				! Number of Moles
						!M_RHO = WTPH(I)/VOLPH(I) 				! Density of melt
						!END IF
						
						!IF (STPHNAM(I).EQ. 'LIQtc_foL8') THEN 	! H2O melt component
						!M_WASSER = NN(I) + M_WASSER				! Number of Moles
						!END IF
						
						!IF (STPHNAM(I).EQ. 'LIQtc_faL8') THEN 	! H2O melt component
						!M_WASSER = NN(I) + M_WASSER				! Number of Moles
						!END IF
						
					END DO
					write(scr,"(F12.6)") M_RHO
					
				CASE ('MMELT')
				! MELT-Number of moles (Si Al Fe Mg Ca Na K H O)
								
					!!==================================
					! Liquids melt CaNa-KFMASH
					M_SI=0
					DO I=1,NUN2
					DO J=1, NUN
						IF (CHNAME(J).EQ.'SI') EXIT
					END DO
						IF (STPHNAM(I)(1:5).EQ.'LIQtc') EXIT
						!IF (STPHNAM(I).EQ.'LIQtc_h2oL') EXIT
						!IF (STPHNAM(I).EQ. 'LIQtc_foL8') EXIT
						!IF (STPHNAM(I).EQ. 'LIQtc_faL8') EXIT 	
					END DO	
					M_SI=X(I,J)*NN(I)
					
					
					
					M_AL=0
					DO I=1,NUN2
					DO J=1, NUN
						IF (CHNAME(J).EQ.'AL') EXIT
					END DO
						IF (STPHNAM(I)(1:5).EQ.'LIQtc') EXIT
					END DO	
					M_AL=X(I,J)*NN(I)
					

					M_FE=0
					DO I=1,NUN2
					DO J=1, NUN
						IF (CHNAME(J).EQ.'FE') EXIT
					END DO
						IF (STPHNAM(I)(1:5).EQ.'LIQtc') EXIT
					END DO	
					M_FE=X(I,J)*NN(I)
					
										
					M_MG=0
					DO I=1,NUN2
					DO J=1, NUN
						IF (CHNAME(J).EQ.'MG') EXIT
					END DO
						IF (STPHNAM(I)(1:5).EQ.'LIQtc') EXIT
					END DO	
					M_MG=X(I,J)*NN(I)
					
					
					M_CA=0
					DO I=1,NUN2
					DO J=1, NUN
						IF (CHNAME(J).EQ.'CA') EXIT
					END DO
						IF (STPHNAM(I)(1:5).EQ.'LIQtc') EXIT
					END DO	
					M_CA=X(I,J)*NN(I)
					
					
					M_NA=0
					DO I=1,NUN2
					DO J=1, NUN
						IF (CHNAME(J).EQ.'NA') EXIT
					END DO
						IF (STPHNAM(I)(1:5).EQ.'LIQtc') EXIT
					END DO	
					M_NA=X(I,J)*NN(I)
					
					
					M_K=0
					DO I=1,NUN2
					DO J=1, NUN
						IF (CHNAME(J).EQ.'K') EXIT
					END DO
						IF (STPHNAM(I)(1:5).EQ.'LIQtc') EXIT 
					END DO	
					M_K=X(I,J)*NN(I)
					
					
					M_H=0
					DO I=1,NUN2
					DO J=1, NUN
						IF (CHNAME(J).EQ.'H') EXIT
					END DO
						IF (STPHNAM(I)(1:5).EQ.'LIQtc') EXIT 
					END DO	
					M_H=X(I,J)*NN(I)
					
					
					M_O=0
					DO I=1,NUN2
					DO J=1, NUN
						IF (CHNAME(J).EQ.'O') EXIT
					END DO
						IF (STPHNAM(I)(1:5).EQ.'LIQtc') EXIT 
					END DO	
					M_O=X(I,J)*NN(I)
					
					write(scr,"(F12.6)") M_SI, M_AL, M_FE, M_MG, &
					 M_CA, M_NA, M_K, M_H, M_O
						
					 
				CASE ('BINGO')
				! XMAPTools (Si Al Fe Mg Ca Na K H O)

				    ! NNS =number of stable solid solution
				    ! NPS =number of stable solid solution
				    NSS=0
				    NSP=0
				    
				    DO I=1,NUN2	
				    IF (ISASOLID(I)) THEN
				    NSP=NSP+1
				    IF (EMCODE(I) .NE. 0) NSS=NSS+1 
                    END IF
				    END DO	
				    
					write(scr,*) NSP	
					write(scr,*) NUN
					write(CH001,*) '(',NUN,'A6)'	
					write(scr,CH001) CHNAME(1:NUN)
					
					write(CH001,*) '(',NUN,'F12.6)'
								
					DO I=1,NUN2		
					IF (ISASOLID(I)) THEN
					write(scr,CH001) X(I,1:NUN)			
					END IF
					END DO	
					

					!write(CH001,*) '(',NUN2,'2X,A16)'
					write(CH001,*) '(',NSP,'A16)'
					WRITE (scr,CH001) STPHNAM(1:NSP)
										
					write(CH001,*) '(',NSP,'F12.6)'
					write(scr,CH001) VOLPH(1:NSP)/VOLSOL*100.0D0	
                    write(scr,CH001) WTPH(1:NSP)/VOLPH(1:NSP)
                    write(scr,*) NSS
                    !write(scr,*) NSOL, NPHA, PHMAX, NMAX, NUN3

																	
				CASE ('KDCTDCAR')
				! KD-value(Fe-Mg) Chloritoid - Carpholite
					DO I=1, NUN2	
				   	  IF ((STPHNAM(I)(1:3).EQ.'CTD'.OR.&
				   	  &STPHNAM(I)(1:4).EQ.'CHTD').OR.&
				   	  &(STPHNAM(I)(1:5).EQ.'CTOID')) THEN
					  !IF (INDEX(STPHNAM(I),'CTD').NE.0) THEN 
						DO J=1,NUN2
						  IF ((STPHNAM(J)(1:3).EQ.'CAR'.OR.&
						  &STPHNAM(J)(1:4).EQ.'CRPH')) THEN
						  !IF (INDEX(STPHNAM(I),'CAR').NE.0) THEN 
          					KDCTDCAR=(MGFE(J)*(1.0D0-MGFE(I)))/&
          					(MGFE(I)*(1.0D0-MGFE(J)))
          				  END IF
						END DO
					  END IF
					END DO
					write(scr,"(F12.6)") KDCTDCAR

				CASE ('KDGRTCPX')
				! KD-value(Fe-Mg) Garnet - Clinopyroxene
					DO I=1, NUN2	
					  IF (STPHNAM(I)(1:6).EQ.'GARNET') THEN 
					    DO J=1,NUN2
						  IF (STPHNAM(J)(1:4).EQ. 'OMPH') THEN 
          					KDGRTCPX=(MGFE(J)*(1.0D0-MGFE(I)))/&
          					(MGFE(I)*(1.0D0-MGFE(J)))
          				  END IF
					    END DO
					  END IF
					END DO
					write(scr,"(F12.6)") KDGRTCPX
					
				CASE ('KDGRTOPX')
				! KD-value(Fe-Mg) Garnet- Orthopyroxene
					DO I=1, NUN2	
					  IF (STPHNAM(I)(1:6).EQ.'GARNET') THEN 
					    DO J=1,NUN2
						  IF (STPHNAM(J)(1:3).EQ. 'OPX') THEN 
          				    KDGRTOPX=(MGFE(J)*(1.0D0-MGFE(I)))/&
          				    (MGFE(I)*(1.0D0-MGFE(J)))
          			      END IF
					    END DO
					  END IF
					END DO
					write(scr,"(F12.6)") KDGRTOPX
					
				CASE ('KDOPXCPX')
				! KD-value(Fe-Mg) Orthopyroxene - Clinopyroxene
					DO I=1, NUN2	
				  	  IF (STPHNAM(I)(1:3).EQ.'OPX') THEN 
						DO J=1,NUN2
					  	  IF (STPHNAM(J)(1:4).EQ. 'OMPH') THEN 
          					KDOPXCPX=(MGFE(J)*(1.0D0-MGFE(I)))/&
          					(MGFE(I)*(1.0D0-MGFE(J)))
          			  	  END IF
						END DO
				  	  END IF
					END DO
					write(scr,"(F12.6)") KDOPXCPX
					
				CASE ('KDOLGRT')
				! KD-value(Fe-Mg) Olivine - Garnet
					DO I=1, NUN2	
					  IF (STPHNAM(I)(1:7).EQ.'OLIVINE') THEN 
						DO J=1,NUN2
						  IF (STPHNAM(J)(1:6).EQ. 'GARNET') THEN 
          					KDOLGRT=(MGFE(J)*(1.0D0-MGFE(I)))/&
          					(MGFE(I)*(1.0D0-MGFE(J)))
          				  END IF
						END DO
					  END IF
					END DO
					write(scr,"(F12.6)") KDOLGRT
					
				CASE ('KDOLOPX')
				! KD-value(Fe-Mg) Olivine - Orthopyroxene
					DO I=1, NUN2	
					  IF (STPHNAM(I)(1:7).EQ.'OLIVINE') THEN 
						DO J=1,NUN2
						  IF (STPHNAM(J)(1:3).EQ.'OPX') THEN 
          					KDOLOPX=(MGFE(J)*(1.0D0-MGFE(I)))/&
          					(MGFE(I)*(1.0D0-MGFE(J)))
          				  END IF
						END DO
					  END IF
					END DO
					write(scr,"(F12.6)") KDOLOPX
					
				CASE ('KDOLCPX')
				! KD-value(Fe-Mg) Olivine - Clinopyroxene
					DO I=1, NUN2	
					  IF (STPHNAM(I)(1:7).EQ.'OLIVINE') THEN 
						DO J=1,NUN2
						  IF (STPHNAM(J)(1:4).EQ. 'OMPH') THEN 
          					KDOLCPX=(MGFE(J)*(1.0D0-MGFE(I)))/&
          					(MGFE(I)*(1.0D0-MGFE(J)))
          				  END IF
						END DO
					  END IF
					END DO
					write(scr,"(F12.6)") KDOLCPX	
												
				CASE ('KDGRTBT')
				! KD-value(Fe-Mg) Garnet - Biotite
					DO I=1, NUN2	
					  IF (STPHNAM(I)(1:6).EQ.'GARNET') THEN 
						DO J=1,NUN2
						  IF (STPHNAM(J)(1:3).EQ. 'BIO') THEN 
          					KDGRTBT=(MGFE(J)*(1.0D0-MGFE(I)))/&
          					(MGFE(I)*(1.0D0-MGFE(J)))
          				  END IF
						END DO
				   	  END IF
					END DO
					write(scr,"(F12.6)") KDGRTBT
					
				CASE ('KDCTDCHL')
				! KD-value(Fe-Mg) Chloritoid - Chlorite
					DO I=1, NUN2	
					  IF ((STPHNAM(I)(1:3).EQ.'CTD'.OR.&
					  &STPHNAM(I)(1:4).EQ.'CHTD').OR.&
					  &(STPHNAM(I)(1:5).EQ.'CTOID')) THEN
						DO J=1,NUN2
						  IF (STPHNAM(J)(1:3).EQ.'CHL') THEN
          					KDCTDCHL=(MGFE(J)*(1.0D0-MGFE(I)))/&
          					(MGFE(I)*(1.0D0-MGFE(J)))
          				  END IF
						END DO
				   	  END IF
					END DO
					write(scr,"(F12.6)") KDCTDCHL
	

		 END SELECT
		 
		 SELECT CASE (CH001(1:2))
		 		CASE ('*X')
				! Number of moles of a element within a phases [mol] 
					write(scr,"(F12.6)") XMOL
				CASE ('*M')
				! Number of moles of a phases [mol] 
					write(scr,"(F12.6)") MOL
				CASE ('*W')
				! Weight percent of a phase
					write(scr,"(F12.6)") WPRO
				CASE ('*V')
				! Volume percent of a phase
					write(scr,"(F12.6)") VPRO
		END SELECT
    ELSE    	
        write(scr,*) 'ERROR-Define output'
    END IF  
!*******************************************************************
     
      END 
	  ! end of program 

      
!===================================================================
!                        SUBROUTINES
!===================================================================
!******************************
!-----test for x=0 -> a=0
!-----test for x=1 -> a=1
!-----test if solution is complex
!-----test for consistency with Gibbs-Duhem equation
!-----test only if MODELL not equal 'I'
!-----muss nach GRECAL (wegen T,P,RT,GG,GGK) aufgerufen werden
!----- eduester = GIBBSTEST without print
      SUBROUTINE GIBBS_D(E1,E2)
      IMPLICIT NONE
      INCLUDE 'theriak.cmn'
      include 'files.cmn'
!-----END OF COMMON VARIABLES
      INTEGER*4 I,IGD,IS,CODE,ITEST,I001,II,E1,E2,IE, &
      IC,NROW,NCOL,RANG,IKOMM
      REAL*8 FF,FFX,XX0(EMAX),XX1(EMAX),XX2(EMAX),F1, &
      GGD0,GGD1,GGD2,XGD,AOFG(EMAX),AOFM(EMAX),GDTEST, &
      XXX(EMAX),AAA(EMAX),MAT(250,250),MUE(EMAX)
!-----
      XGD=0.0001D0
      GDTEST=1D-1
      DO 500,IS=E1,E2
      NNEGEM(IS)=0
      NSTREM(IS)=0
      IKOMM=0

      IF (PRTLOG(2)) IKOMM=1
      IF (NEND(IS).LE.1) GOTO 500
!-----test solution model (only ideal part of site-mixing)
!-----
!-----test x=0 -> a=0
      IF (MODELL(IS).NE.'I') THEN
      FF=1.0D0/DBLE(NEND(IS)-1.0D0)
      DO 510,I=1,NEND(IS)
  510 XXX(I)=FF
      I001=0
      DO 520,I=1,NEND(IS)
      XXX(I)=0.0D0
      CALL ACTIVI(IS,XXX,AAA)
!dC      IF (AAA(I).NE.0.0D0) THEN
      IF (AAA(I).GT.2.0D-10) THEN
      I001=I001+1
      NNEGEM(IS)=I001
      NEGEM(IS,I001)=I
      END IF
  520 XXX(I)=FF
      IF (I001.GT.0) THEN
      IKOMM=IKOMM+1
!
      IF (.NOT.DRU) RETURN
!
      IF (IKOMM.EQ.1) THEN
!      WRITE (UNIT=scr,FMT=1000)
!      WRITE (UNIT=out,FMT=1000)
! 1000 FORMAT (/,' ----------------')
      END IF
!      WRITE (UNIT=scr,FMT=1001) SOLNAM(IS)
!      WRITE (UNIT=out,FMT=1001) SOLNAM(IS)
! 1001 FORMAT (/ &
!      1X,A16,' : ideal part of solution model has ', &
!      'non-zero activity for zero concentration')
!--
      IF (PRTLOG(2).OR.PRTLOG(1)) THEN
!      WRITE (UNIT=scr,FMT='('' '')')
!      WRITE (UNIT=out,FMT='('' '')')
!      WRITE (UNIT=scr,FMT=1010) (ABK(EM(IS,I)),I=1,NEND(IS))
!      WRITE (UNIT=out,FMT=1010) (ABK(EM(IS,I)),I=1,NEND(IS))
! 1010 FORMAT (29X,10(2X,A8))
      DO 530,I=1,NEND(IS)
      XXX(I)=0.0D0
      CALL ACTIVI(IS,XXX,AAA)
!      WRITE (UNIT=scr,FMT=1020) NAME(EM(IS,I)),(AAA(II),II=1,NEND(IS))
!      WRITE (UNIT=out,FMT=1020) NAME(EM(IS,I)),(AAA(II),II=1,NEND(IS))
! 1020 FORMAT (1X,A16,' = 0 -> A:',10F10.6)
  530 XXX(I)=FF
      END IF
!--
      END IF
!-----
      END IF
      IF (NNEGEM(IS).GT.0) THEN
!      WRITE (scr,1022) (NAME(EM(IS,NEGEM(IS,IE))),IE=1,NNEGEM(IS))
!      WRITE (out,1022) (NAME(EM(IS,NEGEM(IS,IE))),IE=1,NNEGEM(IS))
! 1022 FORMAT (/' the following endmembers may have negative', &
!      ' concentrations: ',20(1X,A16))
      END IF
!-----end of solution test 1
!-----
!-----test x=1 -> a=1
      IF (MODELL(IS).NE.'I') THEN
      DO 550,I=1,NEND(IS)
  550 XXX(I)=0.0D0
      I001=0
      DO 555,I=1,NEND(IS)
      XXX(I)=1.0D0
      CALL ACTIVI(IS,XXX,AAA)
      IF (DABS(AAA(I)-1.0D0).GT.1.0D-10) I001=1
  555 XXX(I)=0.0D0
      IF (I001.EQ.1) THEN
      IKOMM=IKOMM+1
      IF (IKOMM.EQ.1) THEN
!      WRITE (UNIT=scr,FMT=1030)
!      WRITE (UNIT=out,FMT=1030)
! 1030 FORMAT (/,' ----------------')
      END IF
!      WRITE (UNIT=scr,FMT=1031) SOLNAM(IS)
!      WRITE (UNIT=out,FMT=1031) SOLNAM(IS)
! 1031 FORMAT (/ &
!      1X,A16,' : ', &
!      'activity is not one for pure endmember'/' ')
!      WRITE (UNIT=scr,FMT=1040) (ABK(EM(IS,I)),I=1,NEND(IS))
!      WRITE (UNIT=out,FMT=1040) (ABK(EM(IS,I)),I=1,NEND(IS))
! 1040 FORMAT (29X,10(2X,A8))
      DO 560,I=1,NEND(IS)
      XXX(I)=1.0D0
      CALL ACTIVI(IS,XXX,AAA)
!      WRITE (UNIT=scr,FMT=1050) NAME(EM(IS,I)),(AAA(II),II=1,NEND(IS))
!      WRITE (UNIT=out,FMT=1050) NAME(EM(IS,I)),(AAA(II),II=1,NEND(IS))
! 1050 FORMAT (1X,A16,' = 1 -> A:',10F10.6)
  560 XXX(I)=0.0D0
!     STOP
      END IF
!-----
      END IF
!-----end of solution test 2
!-----
!-----test if solution is complex
      IF (MODELL(IS).NE.'I') THEN
      NROW=NUN
      NCOL=NEND(IS)
      DO 600,IC=1,NUN
      DO 600,IE=1,NEND(IS)
  600 MAT(IC,IE)=XX(EM(IS,IE),IC)
      RANG=0
      CALL KOMPLEX(NROW,NCOL,MAT,RANG)
      IF (RANG.LT.NEND(IS)) THEN
      IKOMM=IKOMM+1
      IF (IKOMM.EQ.1) THEN
!      WRITE (UNIT=scr,FMT=1060)
!      WRITE (UNIT=out,FMT=1060)
! 1060 FORMAT (/,' ----------------')
      END IF
!      WRITE (UNIT=scr,FMT=1061) SOLNAM(IS)
!      WRITE (UNIT=out,FMT=1061) SOLNAM(IS)
! 1061 FORMAT (/ &
!      1X,A16,' : ', &
!      'solution is complex, do not use "ideal" site mixing.')
      END IF
      END IF
!-----end of complex-test
!-----
!-----test Gibbs-Duhem relation
      DO 800,ITEST=1,NEND(IS)
      CODE=0
      F1=0.0D0
!     FF=1.0D0/DBLE(NEND(IS)+1)
      FF=0.5D0/DBLE(NEND(IS))
      DO 700,I=1,NEND(IS)
      F1=F1+FF
  700 XX0(I)=FF
      XX0(ITEST)=1.0D0-F1+FF
      CALL GNONID(IS,XX0,GGD0)
      CALL MUECAL(IS,XX0,MUE)
      DO 710,IGD=1,NEND(IS)
      DO 720,I=1,NEND(IS)
      XX1(I)=XX0(I)
      XX2(I)=XX0(I)
  720 CONTINUE
      XX1(IGD)=XX0(IGD)-XGD
      XX2(IGD)=XX0(IGD)+XGD
      DO 730,I=1,NEND(IS)
      XX1(I)=XX1(I)/(1.0D0-XGD)
      XX2(I)=XX2(I)/(1.0D0+XGD)
  730 CONTINUE
      FFX=XX2(IGD)-XX1(IGD)
      CALL GNONID(IS,XX1,GGD1)
      CALL GNONID(IS,XX2,GGD2)
      AOFG(IGD)=GGD0+(1.0D0-XX0(IGD))*(GGD2-GGD1)/FFX
      AOFG(IGD)=AOFG(IGD)-GG(EM(IS,IGD))
      AOFM(IGD)=MUE(IGD)-GG(EM(IS,IGD))
      IF (DABS(AOFG(IGD)-AOFM(IGD)).GT.GDTEST) CODE=1
  710 CONTINUE
!-----if inconsistent or PRTLOG(1):
      IF (CODE.EQ.1.OR.PRTLOG(1)) THEN
      IKOMM=IKOMM+1
      IF (CODE.EQ.1) THEN
      IF (IKOMM.EQ.1) THEN
!      WRITE (UNIT=scr,FMT=1070)
!      WRITE (UNIT=out,FMT=1070)
! 1070 FORMAT (/,' ----------------')
      END IF
!      WRITE (UNIT=scr,FMT=1071) SOLNAM(IS)
!      WRITE (UNIT=out,FMT=1071) SOLNAM(IS)
! 1071 FORMAT (/ &
!      1X,A16,' : Gibbs-Duhem test failed: ')
      ELSE
      IF (IKOMM.EQ.1) THEN
!      WRITE (UNIT=scr,FMT=1070)
!      WRITE (UNIT=out,FMT=1070)
      END IF
!      WRITE (UNIT=scr,FMT=2000) SOLNAM(IS)
!      WRITE (UNIT=out,FMT=2000) SOLNAM(IS)
! 2000 FORMAT (/1X,A16,' : Gibbs-Duhem test: ')
      END IF
!---
!      WRITE (UNIT=scr,FMT=2005) (NAME(EM(IS,I)),I=1,NEND(IS))
!      WRITE (UNIT=out,FMT=2005) (NAME(EM(IS,I)),I=1,NEND(IS))
! 2005 FORMAT (17X,10A16)
!      WRITE (UNIT=scr,FMT=2010) (XX0(I),I=1,NEND(IS))
!      WRITE (UNIT=out,FMT=2010) (XX0(I),I=1,NEND(IS))
! 2010 FORMAT (' X         ',10F16.8)
!      WRITE (UNIT=scr,FMT=2015) (AOFM(I),I=1,NEND(IS))
!      WRITE (UNIT=out,FMT=2015) (AOFM(I),I=1,NEND(IS))
! 2015 FORMAT (' MUE(calc) ',10F16.4)
! 2015 FORMAT (' MUE(calc) ',10(2X,1PE14.7))
!      WRITE (UNIT=scr,FMT=2020) (AOFG(I),I=1,NEND(IS))
!      WRITE (UNIT=out,FMT=2020) (AOFG(I),I=1,NEND(IS))
! 2020 FORMAT (' MUE(deriv)',10F16.4)
! 2020 FORMAT (' MUE(deriv)',10(2X,1PE14.7))
      DO 810,I=1,NEND(IS)
      AOFG(I)=AOFG(I)/RT
      IF (AOFG(I).LT.-150) AOFG(I)=-150
      IF (AOFG(I).GT.150) AOFG(I)=150
      AOFG(I)=DEXP(AOFG(I))
      AOFM(I)=AOFM(I)/RT
      IF (AOFM(I).LT.-150) AOFM(I)=-150
      IF (AOFM(I).GT.150) AOFM(I)=150
      AOFM(I)=DEXP(AOFM(I))
  810 CONTINUE
!      WRITE (UNIT=scr,FMT=2025) (AOFM(I),I=1,NEND(IS))
!      WRITE (UNIT=out,FMT=2025) (AOFM(I),I=1,NEND(IS))
! 2025 FORMAT (' A(calc)   ',10F16.8)
! 2025 FORMAT (' A(calc)   ',10(2X,1PE14.7))
!      WRITE (UNIT=scr,FMT=2030) (AOFG(I),I=1,NEND(IS))
!      WRITE (UNIT=out,FMT=2030) (AOFG(I),I=1,NEND(IS))
! 2030 FORMAT (' A(deriv)  ',10F16.8)
! 2030 FORMAT (' A(deriv)  ',10(2X,1PE14.7))
      END IF
  800 CONTINUE
!-----end of Gibbs-Duhem-test
!-----analyze Solution model
      I001=IS
!cdc      IF (NNEGEM(IS).GT.0) CALL ANALSOL(I001)
      IF (NSITE(IS).GT.0) CALL ANALSOL(I001)
      IF (NNEGEM(IS).GT.0.AND.NONEG(IS)) THEN
!      WRITE (UNIT=scr,FMT=2040) SOLNAM(IS)
!      WRITE (UNIT=out,FMT=2040) SOLNAM(IS)
! 2040 FORMAT (/ &
!      1X,A16,' : ', &
!      'all concentrations will be forced to be positive'/' ')
      NNEGEM(IS)=0
      END IF
!-----END analyze Solution model
!++++++
  500 CONTINUE
      RETURN
      END
!-----
!******************************
! by eduester -bulk modulus and thermal expansivity for solids
      SUBROUTINE BULKEXPAN(KT,KS,EXPAN)
      IMPLICIT NONE
      INCLUDE 'theriak.cmn'
      include 'files.cmn'
!-----
      REAL*8 DG(-2:2,-2:2), KT, KS, EXPAN, VG(-2:2,-2:2)
      REAL*8 DG0,DG1,DG2,DGM1,DGM2,TVOR,PVOR,DXT,DXP,SS0,SS1,SSM1, &
      TE0,TE1,TE2,TEM1,TEM2,HH0,HH1,HHM1,LECP, &
      DV0, DV1,DVM1,DVDPX,DS0,DSP1,DSPM1,DSDP, &
      DST1,DSTM1, DSDT, DVDTX, DVT1, DVTM1
      INTEGER*4 I,II
!-----
      ! It's set to 0.1 because of lower sensivity, if DX is smaller
      DXP=1!P/1D5
      DXT=0.1!T/1D5
      TVOR=TC
      !write(*,*) TC, P, T
      PVOR=P
      !Initialize
      DO I=-2,2
        DO II=-2,2
          DG(I,II)=0
        END DO
      END DO
      
!-----
      CALL RECHGSO(DG(0,0))
      TE0=T
!-----
      TC=TVOR+DXT
      CALL SORTOFNURVONPT
      CALL RECHGSO(DG(1,0))
      TE1=T
!-----
      TC=TVOR+2.0D0*DXT
      CALL SORTOFNURVONPT
      CALL RECHGSO(DG(2,0))
      TE2=T
!-----
      TC=TVOR-DXT
      CALL SORTOFNURVONPT
      CALL RECHGSO(DG(-1,0))
      TEM1=T
!-----
      TC=TVOR-2.0D0*DXT
      CALL SORTOFNURVONPT
      CALL RECHGSO(DG(-2,0))
      TEM2=T
      
! Oben
!-----
      P=PVOR+DXP
      TC=TVOR+DXT
      CALL SORTOFNURVONPT
      CALL RECHGSO(DG(1,1))
 !-----
      P=PVOR+DXP
      TC=TVOR
      CALL SORTOFNURVONPT
      CALL RECHGSO(DG(0,1))
!-----
      P=PVOR+DXP
      TC=TVOR-DXT
      CALL SORTOFNURVONPT
      CALL RECHGSO(DG(-1,1))
! Oben-Oben   
!-----
      P=PVOR+2.0D0*DXP
      TC=TVOR
      CALL SORTOFNURVONPT
      CALL RECHGSO(DG(0,2))   
      
! Unten
!-----
      P=PVOR-DXP
      TC=TVOR+DXT
      CALL SORTOFNURVONPT
      CALL RECHGSO(DG(1,-1))
 !-----
      P=PVOR-DXP
      TC=TVOR
      CALL SORTOFNURVONPT
      CALL RECHGSO(DG(0,-1))
!-----
      P=PVOR-DXP
      TC=TVOR-DXT
      CALL SORTOFNURVONPT
      CALL RECHGSO(DG(-1,-1))
! Unten-Unten   
!-----
      P=PVOR-2.0D0*DXP
      TC=TVOR
      CALL SORTOFNURVONPT
      CALL RECHGSO(DG(0,-2))   
!=====
       !DG(T-index,P-index)
       DV0=(DG(0,1)-DG(0,-1))/(2.0D0*DXP)
       DV1=(DG(0,2)-DG(0,0))/(2.0D0*DXP)
       DVM1=(DG(0,0)-DG(0,-2))/(2.0D0*DXP)
       
       DVT1=(DG(1,1)-DG(1,-1))/(2.0D0*DXP)
       DVTM1=(DG(-1,1)-DG(-1,-1))/(2.0D0*DXP)
       DVDTX=(DVT1-DVTM1)/(2.0D0*DXT)
       !DVDPX=(DV0-DVM1)/(1.0D0*DXP)
       !DVDPX=(DV0-DVM1)/(2.0D0*DXP)
       DVDPX=(DV1-DVM1)/(2.0D0*DXP)
       DS0=-(DG(1,0)-DG(-1,0))/(2.0D0*DXT)
       DSP1=-(DG(1,1)-DG(-1,1))/(2.0D0*DXT)
       DSPM1=-(DG(1,-1)-DG(-1,-1))/(2.0D0*DXT)
       DSDP=(DSP1-DSPM1)/(2.0D0*DXP)
       DST1=-(DG(2,0)-DG(0,0))/(2.0D0*DXT)
       DSTM1=-(DG(0,0)-DG(-2,0))/(2.0D0*DXT)
       DSDT=(DST1-DSTM1)/(2.0D0*DXT)      
!=====
  KS=DVDPX+(DSDP*DSDP)/DSDT
  KS=-DV0/KS
  KT=-DV0/DVDPX
  !KT=-4.366D0/DVDPX
  EXPAN=(1/DV0)*DVDTX
  !EXPAN=(1/4.366D0)*DVDTX
            
!=====
      !WRITE (UNIT=6,FMT='(''DG(0,0) = '',F20.8)') DG(0,0)
      !WRITE (UNIT=6,FMT='(''DG(1,0) = '',F20.8)') DG(1,0)
      !WRITE (UNIT=6,FMT='(''DG(2,0) = '',F20.8)') DG(2,0)
      !WRITE (UNIT=6,FMT='(''DG(-1,0) = '',F20.8)') DG(-1,0)
      !WRITE (UNIT=6,FMT='(''DG(-2,0) = '',F20.8)') DG(-2,0)

      !WRITE (UNIT=6,FMT='(''dv0 = '',F20.8)') DV0
      !WRITE (UNIT=6,FMT='(''dv1 = '',F20.8)') DV1
      !WRITE (UNIT=6,FMT='(''dvm1= '',F20.8)') DVM1
      !WRITE (UNIT=6,FMT='(''dvdp = '',F20.8)') DVDPX
      !WRITE (UNIT=6,FMT='(''ds0 = '',F20.8)') DS0
      !WRITE (UNIT=6,FMT='(''dsp1= '',F20.8)') DSP1    
      !WRITE (UNIT=6,FMT='(''dspm1= '',F20.8)') DSPM1             
      !WRITE (UNIT=6,FMT='(''dsdp= '',F20.8)') DSDP                
      !WRITE (UNIT=6,FMT='(''dst1= '',F20.8)') DST1          
      !WRITE (UNIT=6,FMT='(''dstm1= '',F20.8)') DSTM1            
      !WRITE (UNIT=6,FMT='(''dsdt= '',F20.8)') DSDT

      !write(*,*) 'Expansivity', EXPAN   
      !write(*,*) 'Bulk modulus', KT, KS, 1/KT
!=====
      TC=TVOR
      P=P+DXP
      CALL SORTOFNURVONPT
!=====
!==============================
      TC=TVOR
      P=PVOR
      CALL SORTOFNURVONPT
!*****
      RETURN
      END
!-----
!******************************
      SUBROUTINE DSDVTEST5(HH0,SS0,LECP)
      IMPLICIT NONE
      INCLUDE 'theriak.cmn'
      include 'files.cmn'
!-----
      REAL*8 DG0,DG1,DG2,DGM1,DGM2,TVOR,PVOR,DXT,DXP,SS0,SS1,SSM1, &
      TE0,TE1,TE2,TEM1,TEM2,HH0,HH1,HHM1,LECP,KUH,OU,PW,VTOTX
      REAL*8 WTMSOL, WTMTOT
      INTEGER*4 I
      LOGICAL*4 ISASOLID(COMAX)
      COMMON /LOSOLID/ ISASOLID
!-----
      DXP=P/1D5
      DXT=TC/1D5
      TVOR=TC
      PVOR=P
!-----
      CALL RECHG(DG0)
      TE0=T
!-----
      TC=TVOR+DXT
      CALL SORTOFNURVONPT
      CALL RECHG(DG1)
      TE1=T
!-----
      TC=TVOR+2.0D0*DXT
      CALL SORTOFNURVONPT
      CALL RECHG(DG2)
      TE2=T
!-----
      TC=TVOR-DXT
      CALL SORTOFNURVONPT
      CALL RECHG(DGM1)
      TEM1=T
!-----
      TC=TVOR-2.0D0*DXT
      CALL SORTOFNURVONPT
      CALL RECHG(DGM2)
      TEM2=T
!=====
      !WRITE (UNIT=6,FMT='(/''using 5 points'')')
      !WRITE (UNIT=6,FMT='(''g0 = '',F20.8)') DG0
      SS0=-(DG1-DGM1)/(2.0D0*DXT)
      !WRITE (UNIT=6,FMT='(''ss0 = '',F20.8)') SS0
      SS1=-(DG2-DG0)/(2.0D0*DXT)
      !WRITE (UNIT=6,FMT='(''ss1 = '',F20.8)') SS1
      SSM1=-(DG0-DGM2)/(2.0D0*DXT)
      !WRITE (UNIT=6,FMT='(''ssm1= '',F20.8)') SSM1
      HH0=DG0+T*SS0
      !WRITE (UNIT=6,FMT='(''hh0 = '',F20.8)') HH0
      HH1=DG1+TE1*SS1
      !WRITE (UNIT=6,FMT='(''hh1 = '',F20.8)') HH1
      HHM1=DGM1+TEM1*SSM1
      !WRITE (UNIT=6,FMT='(''hhm1= '',F20.8)') HHM1
      LECP=(HH1-HHM1)/(2.0D0*DXT)
      !WRITE (UNIT=6,FMT='(''lecp= '',F20.8)') LECP
      
      ! 'TS_tot'      
      KUH=SS0*TE0
      !WRITE (UNIT=6,FMT='(''TS_tot = '',F20.8)') KUH
!=====
	  ! heat capacity [J/(g*K)]
      DO I=1,NUN2  ! NUN2 = number of stable phases        
          WTMTOT=WTMTOT+WTM(I) ! molar weight of all stable phase

          IF (ISASOLID(I)) THEN
            WTMSOL=WTMSOL+WTM(I) ! molar weight of solids
          ENDIF 
      END DO   
      !write(*,*) WTMTOT   
	  LECP=LECP/WTMTOT
!===== calculate V_tot
      TC=TVOR
      PVOR=P
!-----
      P=PVOR+DXP
      CALL SORTOFNURVONPT
      CALL RECHG(DG1)
!-----
      P=PVOR-DXP
      CALL SORTOFNURVONPT
      CALL RECHG(DGM1)
!-----
      ! 'V_tot'
      VTOTX=(DG1-DGM1)/(2.0D0*DXP)*10.0D0
      !WRITE (UNIT=6,FMT='(''V_tot = '',F20.8)') VTOTX
      
      ! 'PV-_tot'
      PW=-VTOTX/10.0D0*P
      !WRITE (UNIT=6,FMT='(''PV_tot = '',F20.8)') PW
      
      ! 'U_tot'
      OU=DG0+KUH-VTOTX/10.0D0*P
      !WRITE (UNIT=6,FMT='(''U_tot = '',F20.8)') OU
!==============================
      TC=TVOR
      P=PVOR
      CALL SORTOFNURVONPT
!*****
      RETURN
      END
!-----
!******************************
! eduester: new subroutine, same as RECHG(DGX), but only for solids
      SUBROUTINE RECHGSO(DGX)
      IMPLICIT NONE
      INCLUDE 'theriak.cmn'
!
      LOGICAL*4 ISASOLID(COMAX)
      COMMON /LOSOLID/ ISASOLID
!-----
      REAL*8 DGX,XXX(EMAX),GGG
      INTEGER*4 I,II,IS
!-----
      DGX=0.0D0
      DO I=1,NUN2
!++++
      IF (ISASOLID(I)) THEN
!++++
       IF (NUMMER(I).GT.0) THEN
        DGX=DGX+NN(I)*GGK(NUMMER(I))
       ELSE
        IS=EMCODE(I)
        DO II=1,NEND(IS)
         XXX(II)=XEM(I,II)
        END DO
        CALL GNONID(IS,XXX,GGG)
        DO II=1,NEND(IS)
        GGG=GGG-XXX(II)*GG(EM(IS,II))+XXX(II)*GGK(EM(IS,II))
        END DO
        DGX=DGX+NN(I)*GGG
!++++
       END IF
!++++
       END IF
      END DO
!*****
      RETURN
      END
!-----
!******************************
      SUBROUTINE MANLOOP
      IMPLICIT NONE
      INCLUDE 'theriak.cmn'
      include 'files.cmn'
!-----END OF COMMON VARIABLES
      CHARACTER*500 SYREC
!----
   10 CONTINUE
      SYREC=' '
    3 WRITE (scr,110)
      WRITE (out,110)
  110 FORMAT (/ &
      ' -------------------------------'/ &
      ' define Temperature and Pressure'/ &
      ' -------------------------------')
      WRITE (scr,1000)
 1000 FORMAT(' Enter [ "?" | CR | "end" | T(C)  P(bar) ]: ')
      READ (UNIT=kbd,FMT='(A500)',END=99) SYREC
      IF (VERGL(SYREC,'end').OR.SYREC.eq.' ') RETURN
      IF (SYREC.EQ.'?') THEN
         CALL helpme('$THK-LOOP')
         GOTO 10
      END IF
      CALL GELI(SYREC,TC)
      CALL GELI(SYREC,P)
      CALL NURVONPT
      CALL CALSTR
      CALL PRININ
      CALL THERIA
      GOTO 10
!----
   99 CONTINUE
      RETURN
      END
!-----
!********************************
      SUBROUTINE SIMPL(NPTS,BIN1)
      IMPLICIT NONE
      INCLUDE 'theriak.cmn'
      include 'files.cmn'
!
!-----END OF COMMON VARIABLES
      INTEGER*4 I,I1,I2,I3,II,COMAY,EINS,ALLES,JX,NPTS,BIN1
      REAL*8 FF,CHEXX(2,COMAX),YWERT(0:1000),XXMIN,XXMAX,YYMIN,YYMAX, &
      FA,FB,BREIT,HOCH,Y0,Y1,NEUWERT(0:1000),SIZLAB,ANGE,XVAL(0:1000), &
      F1,F2,F3,F4
      CHARACTER*500 SYREC,BULKLINE(2)
      CHARACTER*16 XVARI,YVARI,BLANK
      CHARACTER*80 XTXT,YTXT,DATI,ZEITSTRING
!
      integer ierr,j
!----
      PRTLOG(9)=.TRUE.
      PRTLOG(2)=.TRUE.
      NVARTBL=0
      NROWTBL=0
      COMAY=COMAX
      EINS=1
      BLANK=' '
!     DO 500,I=1,8
! 500 PRTLOG(I)=.FALSE.
!---
      DO 510,I=1,2
      READ (UNIT=dat,FMT='(A500)') SYREC
      CALL GELI(SYREC,FF)
      CALL TAXI(SYREC,FORMUL)
      BULKLINE(I)=FORMUL
      CALL LABLA(FORMUL,I2)
      WRITE (UNIT=out,FMT='(170A1)') (FORMUL(I3:I3),I3=1,I2)
      CALL CHEMIE(COMAY,NC,OXYDE,OXANZ,FORMUL,CHE)
      DO 512,II=1,NC
  512 CHEXX(I,II)=CHE(II)
  510 CONTINUE
      NBUL=3
      BULINE(2)=BULKLINE(1)
      BULINE(3)=BULKLINE(2)
      XXMIN=0.0D0
      XXMAX=1.0D0
      YYMIN=0.0D0
      YYMAX=1.0D0
      BREIT=15.0D0
      HOCH=15.0D0
      XVARI='binary'
      YVARI='G'
!-----
      DO 520,II=1,NC
      CHEM(II)=0.0D0
      DO 520,I=1,2
  520 CHEM(II)=CHEM(II)+CHEXX(I,II)
      CALL DBREAD
!     WRITE (scr,2100) CHNAME(1),CHEM(CHMCOD(1))
      WRITE (out,2100) CHNAME(1),CHEM(CHMCOD(1))
 2100 FORMAT (' Bulk composition: ',A8,1X,F11.6)
      DO 522,I=2,NUN
!     WRITE (scr,1000) CHNAME(I),CHEM(CHMCOD(I))
      WRITE (out,1000) CHNAME(I),CHEM(CHMCOD(I))
 1000 FORMAT (19X,A8,1X,F11.6)
  522 CONTINUE
      CALL NURVONPT
      CALL CALSTR
      ALLES=NSOL
      CALL GIBBSTEST(EINS,NSOL)
!----
      FF=DBLE(NPTS)
      DO 550,JX=0,NPTS
      IF (BIN1.NE.0.AND.BIN1.NE.JX) GOTO 550
      FB=DBLE(JX)/FF
      FA=1.0D0-FB
      XVAL(JX)=FB
      DO 552,I=1,NUN
  552 BULK(I)=0.0D0
      DO 554,I=1,NUN
      II=CHMCOD(I)
  554 BULK(I)=FA*CHEXX(1,II)+FB*CHEXX(2,II)
!----
      CALL TRENNE(130)
      CALL NURVONPT
      CALL CALSTR
      CALL PRININ
      WRITE (UNIT=scr,FMT=1050) JX
      WRITE (UNIT=out,FMT=1050) JX
 1050 FORMAT (/,' THERIAK-loop = ',I6)
      CALL THERIA
      YWERT(JX)=GGTOT
  550 CONTINUE
!-----
!!!      CALL PRTTBL
!
      NVARTBL=NVARTBL+1
      IF (NVARTBL.GT.MAXVARTBL) THEN
      NVARTBL=NVARTBL-1
      END IF
      VARTBL(NVARTBL)='G(rel)'
!
      Y0=YWERT(0)
      Y1=YWERT(NPTS)
      YYMAX=-1D34
      YYMIN=1D34
      FF=DBLE(NPTS)
      DO 600,I=0,NPTS
      NEUWERT(I)=YWERT(I)-Y0-(Y1-Y0)*DBLE(I)/FF
      IF (NEUWERT(I).GT.YYMAX) YYMAX=NEUWERT(I)
      IF (NEUWERT(I).LT.YYMIN) YYMIN=NEUWERT(I)
!
      OUTTBL(I+1,NVARTBL)=NEUWERT(I)
!
      WRITE (*,1005) I,YWERT(I),NEUWERT(I)
      WRITE (out,1005) I,YWERT(I),NEUWERT(I)
 1005 FORMAT (1X,I4,2X,1PE20.12,2X,1PE20.12)
  600 CONTINUE
!
      CALL PRTTBL
!----
      CALL CPUTIME(ZEITSTRING)
      DATI=sdate
      CALL LABLA(DATI,J)
      DATI=DATI(1:J)//'  '//ZEITSTRING
!------------------
!     open UNIT=inf
!------------------
      j=inf
      line=filename(j)(1:fnl(j))//ext(j)
      path=wpath
      akzess=' '
      state='unknown'
      call openfile(j,ierr)
      if(ierr.ne.0) STOP
!-----
!      WRITE (inf,FMT='(A)') DATI
      CALL PUST(inf,DATI)
      WRITE (inf,FMT='(I4)') NBUL
      DO II=1,NBUL
!       WRITE (inf,FMT='(A)') BULINE(II)
       CALL PUST(inf,BULINE(II))
      END DO
      WRITE (inf,FMT='(I4)') NCOMIN
      DO II=1,NCOMIN
!       WRITE (inf,FMT='(A)') COMINS(II)
       CALL PUST(inf,COMINS(II))
      END DO
      CLOSE (UNIT=inf)
!----
!------------------
!     open UNIT=bin
!------------------
      j=bin
      line=filename(j)(1:fnl(j))//ext(j)
      path=wpath
      akzess=' '
      state=' '
      call openfile(j,ierr)
      if(ierr.ne.0) STOP
!-----
      IF (YYMAX-YYMIN.LT.1D-8) THEN
      WRITE (UNIT=*,FMT=2005) YYMIN
 2005 FORMAT (/,' Minimum and maximum are both = ',1PE10.3)
      YYMAX=YYMAX+DABS(YYMAX*0.05D0+1.0D0)
      YYMIN=YYMIN-DABS(YYMIN*0.05D0+1.0D0)
      END IF
!----
      F1=BREIT+1.0D0
      WRITE (bin,2010) F1,F1
 2010 FORMAT ('NULLPT     5   3'/ &
      'FONT      Helvetica'/ &
      'FAT   0.02'/ &
      'PUNKTE   97  0.5   ',F7.4,'  0.25  999  999 '/ &
      'PUNKTE   98  0.5   ',F7.4,'  0.25  999  999 '/ &
      'FAT   0.03')
      F1=BREIT
      F2=-2.0D0
      F3=0.2D0
      WRITE (bin,2011) DATI,F1,F2,F3
 2011 FORMAT ('TEXTB  ',A,/,3(2X,F9.4),'  0  0  0  0')
      F1=0.0D0
      F3=0.2D0
      DO II=1,NBUL
      F2=HOCH+0.8D0+DBLE(II-1)*0.35D0
      CALL LABLA(BULINE(II),I1)
      WRITE (bin,2012) II,BULINE(II)(1:I1),F1,F2,F3
 2012 FORMAT ('TEXT  ','Bulk(',I1,')= ',A,3(2X,F9.4),'  0  0  0  0')
      END DO
      F1=BREIT+1.5D0
      F3=0.2D0
      DO II=1,NCOMIN
      F2=HOCH-DBLE(II-1)*0.35D0
      CALL LABLA(COMINS(II),I1)
      IF (I1.NE.0) THEN
      WRITE (bin,2014) COMINS(II)(1:I1),F1,F2,F3
 2014 FORMAT ('TEXT  ',A,3(2X,F9.4),'  0  0  0  0')
      END IF
 !
      END DO
      WRITE (UNIT=XTXT,FMT=2015) TC,P
 2015 FORMAT ('x,T=',F8.2,'[C],P=',F8.1,'[Bar]')
      CALL COLLAPS(XTXT,I1)
!      XTXT='x'
      YTXT='G(rel) [J]'
      WRITE (bin,FMT='(A/A/A)') 'NPLOG2',XTXT(1:I1),YTXT
      WRITE (bin,2016) XXMIN,XXMAX,YYMIN,YYMAX,BREIT,HOCH
 2016 FORMAT(4(1PE20.12),0P,F10.3,F10.3)
!-----
      WRITE (bin,FMT='(''    0'',I5,''    0    0    0    0'')') &
      NPTS+1
      WRITE (bin,2018) (XVAL(I),NEUWERT(I),3-MIN0(1,I),I=0,NPTS)
 2018 FORMAT (7(2(1PE20.12),I2))
!      SIZLAB=0.35D0
!      ANGE=0.0D0
!      CALL LABLA(YVARI,I2)
!--
!      FORMUL=BULKLINE(1)
!      CALL LABLA(FORMUL,I1)
!      WRITE (bin,2025) XVAL(0),YYMAX,SIZLAB,ANGE,I1,FORMUL
!      FORMUL=BULKLINE(2)
!      CALL LABLA(FORMUL,I1)
!!--
!      WRITE (bin,2025) XVAL(NPTS),YYMAX,SIZLAB,ANGE,I1,FORMUL
! 2025 FORMAT (2(1PE20.12),0PF10.7,F10.4,I5,A170)
      CLOSE (UNIT=bin)
!
!----
!
!
!
      RETURN
      END
!-----
!********************************
      SUBROUTINE ADDTOBUL(TOADD)
      IMPLICIT NONE
      INCLUDE 'theriak.cmn'
!-----END OF COMMON VARIABLES
      REAL*8 FF,TOADD(COMAX),FAK
      INTEGER*4 I,IC,I0
      LOGICAL*4 MORE
!*****
      FAK=1.0D0
      I0=0
      DO 500,IC=1,NC
      IF (TOADD(IC).LT.0.0D0) THEN
      FF=-CHEM(IC)/TOADD(IC)
      IF (FF.LT.FAK) THEN
      FAK=FF
      I0=IC
      END IF
      END IF
  500 CONTINUE
      IF (FAK.LE.0.0D0) RETURN
      DO 700,IC=1,NC
      CHE(IC)=CHEM(IC)+TOADD(IC)*FAK
  700 CONTINUE
      IF (I0.NE.0) CHE(I0)=0.0D0
      MORE=.FALSE.
      DO 701,I=1,NC
      IF (CHE(I).EQ.0.0D0.NEQV.CHEM(I).EQ.0.0D0) MORE=.TRUE.
      CHEM(I)=CHE(I)
  701 CONTINUE
      IF (MORE) THEN
      CALL DBREAD
      ELSE
      DO 702,I=1,NUN
  702 BULK(I)=CHE(CHMCOD(I))
      END IF
!-----
      RETURN
      END
!-----
!********************************
      SUBROUTINE REMOPH(TOADD)
      IMPLICIT NONE
      INCLUDE 'theriak.cmn'
!-----END OF COMMON VARIABLES
      REAL*8 TOADD(COMAX)
      INTEGER*4 I,IC
      LOGICAL*4 MORE
!*****
      DO 700,IC=1,NC
      CHE(IC)=CHEM(IC)+TOADD(IC)
      IF (CHE(IC).LT.0.0D0) CHE(IC)=0.0D0
  700 CONTINUE
      MORE=.FALSE.
      DO 701,I=1,NC
      IF (CHE(I).EQ.0.0D0.NEQV.CHEM(I).EQ.0.0D0) MORE=.TRUE.
      CHEM(I)=CHE(I)
  701 CONTINUE
      IF (MORE) THEN
      CALL DBREAD
      ELSE
      DO 702,I=1,NUN
  702 BULK(I)=CHE(CHMCOD(I))
      END IF
!-----
      RETURN
      END
!-----
!********************************
      SUBROUTINE MINILOOP(CH001)
      IMPLICIT NONE
      INCLUDE 'theriak.cmn'
      include 'files.cmn'
!
!-----END OF COMMON VARIABLES
      INTEGER*4 NPTS,ierr,j,II
      REAL*8 FF,TPSTEP,T1,T2,P1,P2,DET,DEP,TOADD(COMAX),REFRA(50), &
      INFLOW(COMAX),OUTFLOW(COMAX),F1,NNREM(50),X16
      INTEGER*4 I,IC,IR,IP,I001,COMAY,NREP,KOUNT,SUBPH(50)
      CHARACTER*500 CH001,CH002,SYREC
      CHARACTER*32 REPHASE(50),CH16,TEXT
      CHARACTER*8 LOOPCODE
      CHARACTER*80 DATI,ZEITSTRING
      LOGICAL*4 SUBD
!----
!     DO 500,I=1,10
! 500 PRTLOG(I)=.FALSE.
      PRTLOG(9)=.TRUE.
      PRTLOG(2)=.TRUE.
      PRTLOG(6)=.TRUE.
      DO I=1,50
       NNREM(I)=0.0D0
      END DO
      CH16=' '
      NVARTBL=0
      NROWTBL=0
      NREP=0
      KOUNT=0
      COMAY=COMAX
      SUBD=.FALSE.
      CALL LABLA(CH001,I001)
      IF (NCOMIN.LT.50) THEN
      NCOMIN=NCOMIN+1
      COMINS(NCOMIN)='drv-file: '//CH001(1:I001)
      END IF
!------------------
!     open UNIT=drv
!------------------
      j=drv
      line=CH001
      path=wpath
      akzess=' '
      state='old'
      call openfile(j,ierr)
      if(ierr.ne.0) STOP
!-----
   10 CONTINUE
      READ (UNIT=drv,FMT='(A500)',END=888) SYREC
      IF (NCOMIN.LT.50) THEN
      NCOMIN=NCOMIN+1
      COMINS(NCOMIN)=SYREC(1:80)
      END IF
      CH002='special input: '//SYREC
      CALL PUST(scr,CH002)
!      CALL PUST(out,CH002)
      LOOPCODE='BLAH'
      IF (.NOT.SUBD) CALL TAXI(SYREC,LOOPCODE)
!*****
      IF (LOOPCODE(1:1).EQ.'!') GOTO 10
      IF (VERGL(LOOPCODE,'END')) GOTO 888
!*****
!*****
      IF (VERGL(LOOPCODE,'SUB')) THEN
      CALL TAXI(SYREC,CH001)
      CALL TAXI(SYREC,CH002)
!------------------------
!     open UNIT=41 and 42
!------------------------
      CALL LABLA(CH001,I001)
      OPEN (UNIT=41,FILE=CH001(1:I001))
      CALL LABLA(CH002,I001)
      OPEN (UNIT=42,FILE=CH002(1:I001))
!-----
      SUBD=.TRUE.
      GOTO 10
      END IF
!*****
!***** (subduction, not updated)
    1 IF (SUBD) THEN
      DO 110,I=1,NC
      INFLOW(I)=0.0D0
  110 OUTFLOW(I)=0.0D0
      CALL GELI(SYREC,TC)
      CALL GELI(SYREC,P)
      READ (UNIT=41,FMT='(100(2X,F12.6))',END=111) (INFLOW(I),I=1,NC)
  111 CONTINUE
      CALL ADDTOBUL(INFLOW)
      CALL TRENNE(130)
      CALL NURVONPT
      CALL CALSTR
      IF (PRTLOG(2).OR.PRTLOG(3).OR.PRTLOG(4)) CALL PRININ
      KOUNT=KOUNT+1
      WRITE (UNIT=scr,FMT=1050) KOUNT
      WRITE (UNIT=out,FMT=1050) KOUNT
 1050 FORMAT (/,' THERIAK-loop = ',I6)
      CALL THERIA
      CALL PRTTBL
!+++++
      IF (NREP.GT.0) THEN
      DO 150,IP=1,NUN2
      IF (NUMMER(IP).LE.0) THEN
      TEXT=SOLNAM(EMCODE(IP))
      ELSE
      TEXT=NAME(NUMMER(IP))
      END IF
      DO 160,IR=1,NREP
      IF (REPHASE(IR).EQ.TEXT) THEN
!---- phase to remove
      DO 140,IC=1,NC
  140 TOADD(IC)=0.0D0
      DO 170,IC=1,NUN
      TOADD(CHMCOD(IC))= &
      TOADD(CHMCOD(IC))-X(IP,IC)*NN(IP)*REFRA(IR)
  170 CONTINUE
      CALL REMOPH(TOADD)
      IF (SUBPH(IR).EQ.1) THEN
      DO 115,I=1,NC
  115 OUTFLOW(I)=OUTFLOW(I)-TOADD(I)
      END IF
      END IF
!---- end phase to remove
  160 CONTINUE
  150 CONTINUE
      END IF
!---
      WRITE (UNIT=42,FMT='(100(2X,F12.6))') (OUTFLOW(I),I=1,NC)
      GOTO 10
      END IF
!*****
!*****
      IF (VERGL(LOOPCODE,'COMP')) THEN
      CALL TAXI(SYREC,FORMUL)
      CALL TAXI(SYREC,USE)
      CALL CHEMIE(COMAY,NC,OXYDE,OXANZ,FORMUL,CHEM)
      CALL LABLA(USE,LUSE)
      CALL DBREAD
!     IF (PRTLOG(1)) STOP
      GOTO 10
      END IF
!*****
      IF (VERGL(LOOPCODE,'ADD')) THEN
      CALL TAXI(SYREC,FORMUL)
      CALL CHEMIE(COMAY,NC,OXYDE,OXANZ,FORMUL,CHE)
      DO 700,IC=1,NC
      TOADD(IC)=CHE(IC)
  700 CONTINUE
      CALL ADDTOBUL(TOADD)
      GOTO 10
      END IF
!*****
      IF (VERGL(LOOPCODE,'REMOVE')) THEN
      CALL TAXI(SYREC,CH16)
      CALL GELI(SYREC,FF)
      CALL GELI(SYREC,F1)
      IR=0
      DO 600,I=1,NREP
      IF (CH16.EQ.REPHASE(I)) THEN
      REFRA(I)=FF/100.0D0
      IR=I
      SUBPH(IR)=IDINT(F1)
      GOTO 601
      END IF
  600 CONTINUE
  601 CONTINUE
      IF (IR.EQ.0) THEN
      NREP=NREP+1
      IF (NREP.GT.50) THEN
      WRITE (UNIT=*,FMT='('' NREP too big'')')
      STOP
      END IF
      REPHASE(NREP)=CH16
      REFRA(NREP)=FF/100.0D0
      SUBPH(NREP)=IDINT(F1)
      END IF
      GOTO 10
      END IF
!*****
      IF (VERGL(LOOPCODE,'TP').OR.VERGL(LOOPCODE,'REF')) THEN
      T1=TC
      P1=P
      CALL GELI(SYREC,T2)
      CALL GELI(SYREC,P2)
      CALL GELI(SYREC,TPSTEP)
      IF (TPSTEP.LT.1.0D0) TPSTEP=1.0D0
      NPTS=IDINT(TPSTEP)
      DET=(T2-T1)/TPSTEP
      DEP=(P2-P1)/TPSTEP
      TC=T1
      P=P1
!-----
      DO 800,I=1,NPTS
      CALL TRENNE(130)
      TC=TC+DET
      P=P+DEP
      CALL NURVONPT
      CALL CALSTR
      IF (PRTLOG(2).OR.PRTLOG(3).OR.PRTLOG(4)) CALL PRININ
      KOUNT=KOUNT+1
      WRITE (UNIT=scr,FMT=1000) KOUNT
      WRITE (UNIT=out,FMT=1000) KOUNT
 1000 FORMAT (/,' THERIAK-loop = ',I6)
      CALL THERIA
      IF (REFAS) CALL SHOWREF
!+++++
      IF (NREP.GT.0) THEN
       DO IP=1,NUN2
        IF (NUMMER(IP).LE.0) THEN
         TEXT=SOLNAM(EMCODE(IP))
        ELSE
         TEXT=NAME(NUMMER(IP))
        END IF
        DO IR=1,NREP
         IF (REPHASE(IR).EQ.TEXT) THEN
          NNREM(IR)=NNREM(IR)+NN(IP)*(REFRA(IR))
         END IF
          CH16='nsum_'//REPHASE(IR)
          X16=NNREM(IR)
          CALL SETTBL(CH16,X16)
        END DO
       END DO
      END IF
!+++++
      CALL PRTTBL
!+++++
      IF (VERGL(LOOPCODE,'REF')) CALL MACHREF
!+++++
      IF (NREP.GT.0) THEN
      DO 650,IP=1,NUN2
      IF (NUMMER(IP).LE.0) THEN
      TEXT=SOLNAM(EMCODE(IP))
      ELSE
      TEXT=NAME(NUMMER(IP))
      END IF
      DO 660,IR=1,NREP
      IF (REPHASE(IR).EQ.TEXT) THEN
!---- phase to remove
      DO 640,IC=1,NC
  640 TOADD(IC)=0.0D0
      DO 670,IC=1,NUN
      TOADD(CHMCOD(IC))= &
      TOADD(CHMCOD(IC))-X(IP,IC)*NN(IP)*REFRA(IR)
  670 CONTINUE
      CALL REMOPH(TOADD)
      END IF
!---- end phase to remove
  660 CONTINUE
  650 CONTINUE
      END IF
!+++++
  800 CONTINUE
      END IF
!*****
      GOTO 10
!---
  888 CONTINUE
      CLOSE (UNIT=drv)
      CALL PRTTBL
!----
      CALL CPUTIME(ZEITSTRING)
      DATI=sdate
      CALL LABLA(DATI,J)
      DATI=DATI(1:J)//'  '//ZEITSTRING
!------------------
!     open UNIT=inf
!------------------
      j=inf
      line=filename(j)(1:fnl(j))//ext(j)
      path=wpath
      akzess=' '
      state='unknown'
      call openfile(j,ierr)
      if(ierr.ne.0) STOP
!-----
!      WRITE (inf,FMT='(A)') DATI
      CALL PUST(inf,DATI)
      WRITE (inf,FMT='(I4)') NBUL
      DO II=1,NBUL
!       WRITE (inf,FMT='(A)') BULINE(II)
       CALL PUST(inf,BULINE(II))
      END DO
      WRITE (inf,FMT='(I4)') NCOMIN
      DO II=1,NCOMIN
!       WRITE (inf,FMT='(A)') COMINS(II)
       CALL PUST(inf,COMINS(II))
      END DO
      CLOSE (UNIT=inf)
!----
      RETURN
      END
!-----
!********************************
      SUBROUTINE MACHREF
      IMPLICIT NONE
      INCLUDE 'theriak.cmn'
      include 'files.cmn'
!
!-----END OF COMMON VARIABLES
      INTEGER*4 I,II,IS
!-----
      REFAS=.TRUE.
      RENUN2=NUN2
      DO I=1,NUN2
       RENUM(I)=NUMMER(I)
       REEMC(I)=EMCODE(I)
       RENN(I)=NN(I)
       IF (NUMMER(I).EQ.0) THEN
        IS=EMCODE(I)
        DO II=1,NEND(IS)
         REXEM(I,II)=XEM(I,II)
        END DO
       END IF
      END DO
!-----
      RETURN
      END
!-----
!********************************
      SUBROUTINE SHOWREF
      IMPLICIT NONE
      INCLUDE 'theriak.cmn'
      include 'files.cmn'
!
!-----END OF COMMON VARIABLES
      INTEGER*4 I,II,IS,IGAR,ISTA,IAND
      REAL*8 SUMME,XXSC(EMAX),GSC
      CHARACTER*32 CH16
      REAL*8 X16
!-----
      SUMME=0.0D0
      DO I=1,RENUN2
       IF (RENUM(I).EQ.0) THEN
        IS=REEMC(I)
        DO II=1,NEND(IS)
         XXSC(II)=REXEM(I,II)
        END DO
        CALL GNONID(IS,XXSC,GSC)
!      WRITE (6,1000) I,RENN(I),GSC
!      WRITE (out,1000) I,RENN(I),GSC
! 1000 FORMAT ('sol: I,NN(I),G(I)=',I2,2(2X,1PE15.8))
       ELSE
        GSC=GG(RENUM(I))
!      WRITE (6,1002) I,RENN(I),GSC
!      WRITE (out,1002) I,RENN(I),GSC
! 1002 FORMAT ('non: I,NN(I),G(I)=',I2,2(2X,1PE15.8))
       END IF
       SUMME=SUMME+RENN(I)*GSC
      END DO
!-----
!      IGAR=0
!      ISTA=0
!      IAND=0
!      DO I=1,NUN2
!       IF (NUMMER(I).EQ.0) THEN
!        IF (SOLNAM(EMCODE(I)).EQ.'Grt') IGAR=I
!        IF (SOLNAM(EMCODE(I)).EQ.'St') ISTA=I
!       ELSE
!        IF (NAME(NUMMER(I)).EQ.'andalusite') IAND=I
!       END IF
!      END DO
!-----
      WRITE (6,2000) SUMME
      WRITE (out,2000) SUMME
 2000 FORMAT (' Difference to reference: ',1PE15.8)
      CH16='G_overstep'
      X16=SUMME
      CALL SETTBL(CH16,X16)
!-----
!      IF (IGAR.NE.0.AND.TC.GT.544.6) THEN
!      WRITE (6,2010) NN(IGAR)
!      WRITE (out,2010) NN(IGAR)
! 2010 FORMAT (' amount of garnet: ',1PE15.8)
!      X16=SUMME/NN(IGAR)/12.0D0
!      CH16='Gref_Grt'
!      CALL SETTBL(CH16,X16)
!      END IF
!-----
!      IF (ISTA.NE.0.AND.TC.GT.555.3) THEN
!      WRITE (6,2012) NN(ISTA)
!      WRITE (out,2012) NN(ISTA)
! 2012 FORMAT (' amount of staurolite: ',1PE15.8)
!      X16=SUMME/NN(ISTA)/48.0D0/(TC-555.3)
!      CH16='Gref_St'
!      CALL SETTBL(CH16,X16)
!      END IF
!-----
!      IF (IAND.NE.0.AND.TC.GT.558.9) THEN
!      WRITE (6,2014) NN(IAND)
!      WRITE (out,2014) NN(IAND)
! 2014 FORMAT (' amount of andalusite: ',1PE15.8)
!      X16=SUMME/NN(IAND)/5.0D0/(TC-558.9)
!      CH16='Gref_And'
!      CALL SETTBL(CH16,X16)
!      END IF
!-----
      RETURN
      END
!-----
