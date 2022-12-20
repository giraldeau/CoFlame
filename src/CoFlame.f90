!*****************************************************************
! CoFlame Code: Version 1.7, April 1st, 2015  
! A parallel CFD code for modeling laminar coflow diffusion flames (with soot included)
! This parallel code was developed by N. Eaves, Q. Zhang, S.B. Dworkin, F. Lui, H. Guo, and M. J. Thomson at the 
! University of Toronto, National Research Council, and Ryerson University.
! When using this code, please cite:
! M. Saffaripour, M. Kholghy, S.B. Dworkin, M.J. Thomson, Proc. Combust. Inst. 34 (2013) 1057-1065.
! N.A. Eaves, A. Veshkini, C. Riese, Q. Zhang, S.B. Dworkin, M.J. Thomson, Combust. Flame 159 (2012) 3179-3190
! Q. Zhang, H. Guo, F. Liu, G.J. Smallwood, M.J. Thomson, Proc. Combust. Inst. 32 (2009) 697-705
!*****************************************************************
! RDB modifications
! Version 05-05-19
!*****************************************************************
PROGRAM CoFlame
    ! *********************************************************************************************
    !This program controls the overall code. It sets up the grid, sets initial conditons, then directs the
    !iterative solution procedure.
    !**********************************************************************************************
    !Global Varaible Declaractions
    implicit none
    
    INCLUDE 'mpif.h'
    INTEGER MYID, NUMP, RC, IERR
    
!

!    INTEGER, PARAMETER :: LENIWK=73000,LENRWK=48000000,LENCWK=800,LENSYM=16,LENITP=100000,LENTP=200000,LINKCK=25,LINKTP=35,LOUT=6
    INTEGER, PARAMETER :: LENIWK=200000,LENRWK=20000000,LENCWK=8000,LENSYM=16
    INTEGER, PARAMETER :: LENITP=100000,LENTP=1022000,LINKCK=25,LINKTP=35,LOUT=6
!    INTEGER, PARAMETER :: LENIWK=23000,LENRWK=21000,LENCWK=300,LENSYM=16
!    INTEGER, PARAMETER :: LENITP=1000,LENTP=1022000,LINKCK=25,LINKTP=35,LOUT=6
    DOUBLE PRECISION ::  WORK(LENRWK),TPWRK(LENTP)
    integer :: IWORK(LENIWK) ,ITPWRK(LENITP) ! RDB originalmente declarado como double
    CHARACTER :: CWORK(LENCWK)*(LENSYM)
    INTEGER :: LENICK, LENRCK, LENCCK, LENIMC, LENRMC

    !CHEMKIN FILES
    ! N°1 -- DLR 2009 Kinetic Mechanism - N.A. Slavinskaya and P. Frank C&F 156 (2009)
    !CHARACTER(LEN=12), PARAMETER :: CHEMfile = 'chem2009.bin' , TRANfile = 'tran2009.bin'
    !INTEGER, PARAMETER :: NumNuc = 3 !# PAH for nucleation process
    ! N°2 -- DLR 2011 Kinetic Mechanism - S. Dworkin et al. C&F 158 (2011)
    CHARACTER(LEN=12), PARAMETER :: CHEMfile = 'chem.bin' , TRANfile = 'tran.bin'
!    CHARACTER(LEN=17), PARAMETER :: CHEMfile = 'chem2011toto1.bin' , TRANfile = 'tran2011toto1.bin'
    INTEGER, PARAMETER :: NumNuc = 3 !# PAH for nucleation process
    ! N°3 -- DLR 2014 Kinetic Mechanism - V. Chernov et al. C&F 161 (2014)
    !CHARACTER(LEN=12), PARAMETER :: CHEMfile = 'chem2014.bin' , TRANfile = 'tran2014.bin'
    !INTEGER, PARAMETER :: NumNuc = 3 !# PAH for nucleation process
    ! N°4 -- ABF Kinetic Mechanism - J. Appel, H. Bockhorn and M. Frenklach C&F 121 (2000).
!    CHARACTER(LEN=12), PARAMETER :: CHEMfile = 'chemABF1.bin' , TRANfile = 'tranABF1.bin'
!    INTEGER, PARAMETER :: NumNuc = 1
    ! N°5 -- MARINOV Kinetic Mechanism - Marinov, Combustion and Flame 114 (1998).
    !CHARACTER(LEN=12), PARAMETER :: CHEMfile = 'chemMari.bin' , TRANfile = 'tranMari.bin'
    !INTEGER, PARAMETER :: NumNuc = 1
     ! N°6 -- Combined_DLR_Chung_For_Gasoline. Kinetic mechanism provided by F.Liu. Used in C&F 180(2017)
!    CHARACTER(LEN=16), PARAMETER :: CHEMfile = 'chemDLRchung.bin' , TRANfile = 'tranDLRchung.bin'
!    INTEGER, PARAMETER :: NumNuc = 3 !# PAH for nucleation process
!    CHARACTER(LEN=14), PARAMETER :: CHEMfile = 'chemKM2red.bin' , TRANfile = 'tranKM2red.bin'
!    CHARACTER(LEN=11), PARAMETER :: CHEMfile = 'chemKM2.bin' , TRANfile = 'tranKM2.bin'
!    CHARACTER(LEN=16), PARAMETER :: CHEMfile = 'chemKM2toto1.bin' , TRANfile = 'tranKM2toto1.bin'
!    INTEGER, PARAMETER :: NumNuc = 3 !# PAH for nucleation process
!    CHARACTER(LEN=14), PARAMETER :: CHEMfile = 'chemExtRaj.bin' , TRANfile = 'tranExtRaj.bin'
!    INTEGER, PARAMETER :: NumNuc = 3 !# PAH for nucleation process

    integer :: L3, L2, L1, L0, LM1, LM2, M3, M2, M1, KK, MSection, I, J, K, K1, K2, indp, FIRST
    
    DOUBLE PRECISION :: DT, DTmax, PRESSURE, ALPHA_SURF, FS, minslopetar, slopetar, ABSO, TUBEEND
    
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: STATUS

    integer :: mechanism, sootmodel, last, iter, miniter, dtiter, printinitial, saveint, outbound, sootdetail
    
    logical :: lrun
    
    integer, dimension (6) :: LSOLVE
    
    DOUBLE PRECISION, DIMENSION(6) :: RELAX
    
    INTEGER, DIMENSION(5) :: NTIMES
    
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WT, XS
    
    LOGICAL,DIMENSION(:,:),ALLOCATABLE :: SOLIDCVP, SOLIDCVU, SOLIDCVV
    
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: X, XU, XDIF, XCV, XCVS, XCVI, XCVIP, FX, FXM
    
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: Y, YV, YDIF, YCV, YCVS, YCVR, YCVRS, ARX, ARXJ, ARXJP, R, RMN, FV, FVP, FY, FYM
    
    DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: U, V, P, PC, T, Rho, DOM_QR4Print
    
    DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE :: SpeciesMF, SootSec, SpeciesMFO, SootSecO 
    
    DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: FlowN, FlowE
    
    DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: DU, DV 
    
    DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: DTDX, DTDY, VCX, VCY
    
    DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE :: VKX, VKY, VTX, VTY, VSX, VSY, VSTX, VSTY
    
    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: SootGasRates
    
    DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: SootHacaRates, SootMolConc

    DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: SootNucRates
    
    DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: SootCoagRates
    
    DOUBLE PRECISION, PARAMETER :: AV = 6.022D+23, AMU = 1.d0/AV, C_MW = 12.011d0, C_MASS = C_MW*AMU, MaxTemp = 3000.d0
    DOUBLE PRECISION, PARAMETER :: DensityP = 1.9d0 , BOLTZMANN = 1.3807D-16
    
    DOUBLE PRECISION, PARAMETER :: SOLIDVIS = 1.0D+100, smallnum = 1.0D-100, CONDSOLID = 1.626D+06, CPSOLID = 5.10D+06
    DOUBLE PRECISION, PARAMETER :: SOLIDDENSITY = 7.85d0, PI = 3.141592d0

    DOUBLE PRECISION, PARAMETER :: GX = -980.d0, GY = 0.d0, PCKIN = 1013250.d0
    
    DOUBLE PRECISION, PARAMETER :: DFRCT = 1.8d0, AKF = 1.37d0, FVOL = 1.43d0  
    
    INTEGER, PARAMETER ::  NumSootReaction = 9
    
    DOUBLE PRECISION, PARAMETER :: C_aggl = 3.0d0, NumCIncep = 40.0d0, c_min_mono = 90000.d0
    
    DOUBLE PRECISION, PARAMETER :: Rad_OH = 2.750d-8/2.d0 , Rad_BAPYR=8.087D-8/2.d0, Rad_BGHIF=7.559D-8/2.d0
    DOUBLE PRECISION, PARAMETER :: Rad_BAPYRS=8.087D-8/2.d0 , Rad_C6H6 = 5.290D-8/2.d0 , Rad_A4 = 7.240D-8/2.d0
    
    INTEGER :: IC2H2, IBAPYR, IBAPYRS, IBGHIF, IO2, IH2O, IH, IH2, IC2H6, IN2, ICH4, IO, IOH, ICO
    INTEGER :: ICO2, IC2H4, IC6H6, IHE, IAR, IC3H6, IC3H9, IC4H6, IC4H8, IA4
    INTEGER :: IC7H16, IC7H8, IDME,IC3H8
    
    CHARACTER (1800) :: species_header,species_header2
    DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Rate         ! RDB tasas de formacion
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: RateLimit      ! RDB tasas de formacion
    DOUBLE PRECISION, PARAMETER :: totoratelimit = 0.0d0            ! RDB totoratelimit
    INTEGER, PARAMETER ::  reacdetail = 1     ! RDB tasas de formacion   write a detailed reaction rate results file
    INTEGER, PARAMETER ::  IREACP = 54        ! tasas de formacion   # reacciones para presentar el campo
    INTEGER :: IREAC, iR(IREACP)
    DOUBLE PRECISION :: DTDXtoto, DTDYtoto

    !End variable declaration

    ! ==========================

    !*********************PROGRAM START****************************************************

    CALL MPI_INIT(IERR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)   !give each CPU a rank
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMP,IERR)   !determine total number of CPUs
    
    lrun = .true. !set the program run flag to true
    
    if(myid==0)WRITE(6,*)
    if(myid==0)WRITE(6,*)
    if(myid==0)WRITE(6,*)'This calculation will be run on ',NUMP,' CPUs.'    
    ! ---------------------------------------------------------------------------
       
    CALL Prepare          !reads in the mechanism/input.dat file,other prep work
    call MPI_BARRIER(MPI_COMM_WORLD,IERR)
    
    ! RDB numero de reacciones del macanismo
!    IREAC = 723
!    IREAC = 960
     IREAC = 733 

    !set up array/loop parameters
    L2=L3+1
    L1=L3+2
    L0=L3+3
    LM1=L3+4
    LM2=L3+5
    
    M2=M3+1
    M1=M3+2
    
    !Set up equation relaxation factors, and how many TMDA iterations are performed for each equation each iteration
    RELAX(1) = 0.5d0 !0.2d0    !0.5d0 !U
    RELAX(2) = 0.5d0 !0.2d0    !0.5d0 !V
    RELAX(3) = 0.5d0    !0.75d0 !PC   0.2d0 !
    RELAX(4) = 0.5d0 !T
    RELAX(5) = 0.7d0 !RHO
    RELAX(6) = 0.5d0 !soot, when using 2-eq model
    
    NTIMES(1) = 12  !6 !U
    NTIMES(2) = 12  !6 !V
    NTIMES(3) = 40 !40 !PC
    NTIMES(4) = 10  !6  !T
    NTIMES(5) = 40 !soot, when using 2-eq model
    
    !allocated arrays based on input.dat file
    allocate (WT(KK))
    allocate (XS(MSection))
    allocate (SOLIDCVP(0:LM1,M1))
    allocate (SOLIDCVU(0:LM1,M1))
    allocate (SOLIDCVV(0:LM1,M1))
    
    allocate (U(0:LM1,M1))
    allocate (V(0:LM1,M1))
    allocate (P(0:LM1,M1))
    allocate (PC(0:LM1,M1))
    allocate (T(0:LM1,M1))
    allocate (Rho(0:LM1,M1))
    allocate (SpeciesMF(0:LM1,M1,KK))
    allocate (SpeciesMFO(0:LM1,M1,KK))
    
    allocate (SootSec(0:LM1,M1,MSection*2))
    allocate (SootSecO(0:LM1,M1,MSection*2))    
    
    allocate (X(0:L0))
    allocate (XU(0:LM1))
    allocate (XDIF(L0))
    allocate (XCV(0:L0))
    allocate (XCVS(L0))
    allocate (XCVI(0:L0))
    allocate (XCVIP(0:L0))
    allocate (FX(L0))
    allocate (FXM(L0))
    
    allocate (Y(M1))
    allocate (YV(M1))
    allocate (YDIF(M1))
    allocate (YCV(M1))
    allocate (YCVS(M1))
    allocate (YCVR(M1))
    allocate (YCVRS(M1))
    allocate (ARX(M1))
    allocate (ARXJ(M1))
    allocate (ARXJP(M1))
    allocate (R(M1))
    allocate (RMN(M1))
    allocate (FV(M1))
    allocate (FVP(M1))
    allocate (FY(M1))
    allocate (FYM(M1))
    
    allocate (DU(0:LM1,M1))
    allocate (DV(0:LM1,M1))
    
    allocate (FLOWN(0:LM1,M1))
    allocate (FLOWE(0:LM1,M1))
    
    allocate (VSX(L1,M1,2*MSection))
    allocate (VSY(L1,M1,2*MSection))
    allocate (VSTX(L1,M1,2*MSection))
    allocate (VSTY(L1,M1,2*MSection))
    
    allocate (VKX(L1,M1,KK))
    allocate (VKY(L1,M1,KK))
    allocate (VTX(L1,M1,KK))
    allocate (VTY(L1,M1,KK))
    
    allocate (DTDX(L1,M1))
    allocate (DTDY(L1,M1))

    allocate (VCX(L1,M1))
    allocate (VCY(L1,M1))
    
    allocate (SootGasRates(L1,M1,MSection,NumSootReaction)) 
    allocate (SootHacaRates(L1,M1,13),SootMolConc(L1,M1,10))
    allocate (SootNucRates(L1,M1,2*NumNuc))   
    allocate(SootCoagRates(L1,M1,2*MSection))   
    allocate(DOM_QR4Print(L3*NUMP+2,M1))

    allocate (Rate(L1,M1,IREACP))       ! RDB tasas de formación
    allocate (RateLimit(2,IREAC))    ! RDB tasas de formación

    !done allocating arrays
    
    call MPI_BARRIER(MPI_COMM_WORLD,IERR)
    
    ! indices reacciones interesantes
    iR=(/13,14,15,17,18,19,20,37,38,39,64,65,67,68,70,74,78,84,85,88,90,91,96,100,104,105,111,112,113,117,&
    118,125,126,130,131,133,148,150,151,152,154,155,156,163,174,175,177,192,205,246,247,251,252,253/)

    !initialize main solution arrays
    PC=0.0d0
    P=0.0d0
    RHO=1.0d0
    U=0.0d0
    V=0.0d0
    SpeciesMF = 0.0d0
    SpeciesMFO = 0.0d0
    SootSec = 0.0d0
    SootSecO = 0.0d0
    T = 300.d0
    VKY = 0.0d0
    VKX = 0.0d0
    VTX = 0.0d0
    VTY = 0.0d0
    VSX = 0.0d0
    VSY = 0.0d0
    VSTX = 0.0d0
    VSTY = 0.0d0
    VCX = 0.0d0
    VCY = 0.0d0
    DTDY = 0.0d0
    DTDX = 0.0d0
    SootGasRates=0.0d0
    SootHacaRates=0.0d0
    SootMolConc=0.0d0
    SootNucRates=0.0d0
    FLOWN = 0.0d0
    FLOWE = 0.0d0
    
    FIRST = 0  
            
    ! -----Get species molecular weight--------------------------    
    CALL CKWT(IWORK,WORK,WT)
    
    !set parameters for soot model
    if(SootModel.EQ.1) then    
        
        !calculate sectional spacing
        XS(1) = DLOG(C_MASS*NumCIncep)

        DO I=2,MSection
            XS(I)=XS(I-1)+DLOG(FS)
        ENDDO
        !done calculating sectional spacing
    endif   
    
    CALL ConstructMesh     !determines locations of nodes, and grid points
    call MPI_BARRIER(MPI_COMM_WORLD,IERR)
    
    CALL SetInitial       !sets the initial conditions for the calculation, and flags CVs as solid (if doing CHT) 
    call MPI_BARRIER(MPI_COMM_WORLD,IERR)
    
    CALL UpdateDensity  !updates the density field
    call MPI_BARRIER(MPI_COMM_WORLD,IERR)   
    
    !if printinitial equals 1, write a solution and tecplot file of the initial solution
    if(printinitial.eq.1) then 
        call OutputResultsf
        call MPI_BARRIER(MPI_COMM_WORLD,IERR) 
    endif 
    
    CALL OutputResultss !writes out the residuals, and increases timestep if suitable. also stops code if divergence occurs 
    call MPI_BARRIER(MPI_COMM_WORLD,IERR)  

    Do while (LRUN) !main program loop
    
        ITER=ITER+1   !increase iteration counter                
        
        CALL SetBoundary    !sets up the boundary conditions
        call MPI_BARRIER(MPI_COMM_WORLD,IERR) 

        !begin solving the conservation equations
        if(LSOLVE(1).EQ.1) CALL SolveU
        call MPI_BARRIER(MPI_COMM_WORLD,IERR)
        
        if(LSOLVE(2).EQ.1) CALL SolveV
        call MPI_BARRIER(MPI_COMM_WORLD,IERR)
        
        if(LSOLVE(3).EQ.1) CALL SolvePC
        call MPI_BARRIER(MPI_COMM_WORLD,IERR)       
!OBS1
        if((LSOLVE(4).EQ.1).OR.(LSOLVE(5).EQ.1).OR.(LSOLVE(6).EQ.1)) then
            CALL CVFaceFlow
            call MPI_BARRIER(MPI_COMM_WORLD,IERR)
            
            DO J=1,M1
                DO I=1, L1                
                    !CALL CVGradients(T,DTDX(I,J),DTDY(I,J),0)
                    CALL CVGradients(T,DTDXtoto,DTDYtoto,0)     !RDB
                    DTDX(I,J)=DTDXtoto      !/2.0d0 !RDB
                    DTDY(I,J)=DTDYtoto      !/2.0d0 !RDB
                ENDDO      
            ENDDO 
            
            call MPI_BARRIER(MPI_COMM_WORLD,IERR)            
            CALL SootDiff
            call MPI_BARRIER(MPI_COMM_WORLD,IERR)
!OBS2            
            if((LSOLVE(4).EQ.1).OR.(LSOLVE(6).EQ.1)) then
                CALL SpeciesDiff
                call MPI_BARRIER(MPI_COMM_WORLD,IERR)         
                Call CorDiffVel
                call MPI_BARRIER(MPI_COMM_WORLD,IERR)
            endif
        endif      
        
        if(LSOLVE(4).EQ.1) then 
            
             ! test SpeciesMF
             ! RDB
!            if(myid.eq.0) then
!            do I= 0,LM1
!            do J= 1,M1
!            do K= 1,KK
!                if(SpeciesMF(I,J,K).NE.SpeciesMF(I,J,K)) print*, 'toto1', myid,I,J,K,SpeciesMF(I,J,K), 'testSpecies'
!            enddo
!            enddo
!            enddo
!            endif

            CALL SolveSpecies

             ! test SpeciesMF
             ! RDB
!            if(myid.eq.0) then
!            do I= 0,LM1
!            do J= 1,M1
!            do K= 1,KK
!                if(SpeciesMF(I,J,K).NE.SpeciesMF(I,J,K)) print*, 'toto2', myid,I,J,K, 'testSpecies'
!            enddo
!            enddo
!            enddo
!            endif

            call MPI_BARRIER(MPI_COMM_WORLD,IERR)
            CALL UpdateDensity  !updates the density field          !RDB
            call MPI_BARRIER(MPI_COMM_WORLD,IERR)                   !RDB
        endif     
        
        if(LSOLVE(5).EQ.1) then 
            
            CALL SolveSoot
            call MPI_BARRIER(MPI_COMM_WORLD,IERR)         
                    
        endif      
        
        if(LSOLVE(6).EQ.1) then
            if((LSOLVE(5).EQ.1).OR.(LSOLVE(4).EQ.1)) then !OBS3
                if(LSOLVE(5).EQ.1) then
                    CALL SootDiff
                    call MPI_BARRIER(MPI_COMM_WORLD,IERR)
                endif
                if(LSOLVE(4).EQ.1) then
                    CALL SpeciesDiff
                    call MPI_BARRIER(MPI_COMM_WORLD,IERR)
                endif     
                Call CorDiffVel
                call MPI_BARRIER(MPI_COMM_WORLD,IERR)
            endif   
               
            CALL SolveT
            call MPI_BARRIER(MPI_COMM_WORLD,IERR)
            
        endif
        
        CALL UpdateDensity  !updates the density field
        call MPI_BARRIER(MPI_COMM_WORLD,IERR)         

        !new solution obtained!
        
        !output the results from the current iteration
        CALL OutputResultss !writes out the residuals, and increases timestep if suitable. also stops code if divergence occurs 
        call MPI_BARRIER(MPI_COMM_WORLD,IERR)
        if(iter.gt.0.and.mod(iter,saveint)==0) call OutputResultsf !writes out current solution restart file and two tecplot formatted files
        call MPI_BARRIER(MPI_COMM_WORLD,IERR)
        !done output
           
        IF(ITER .GE. LAST) THEN 
            LRUN=.FALSE.  !if we are done iterations, stop the program
            call OutputResultsf !writes out current solution restart file and two tecplot formatted files
            call MPI_BARRIER(MPI_COMM_WORLD,IERR)
        ENDIF
            
        IF(ITER.GE.LAST.AND.MYID.EQ.0)WRITE(6,*)'ITERATIONS ARE COMPLETE.'
        
        call MPI_BARRIER(MPI_COMM_WORLD,IERR)
        
    enddo     !done main program loop
    
    CALL MPI_FINALIZE(RC)
        
    CONTAINS
        
        !***************************************************************************************
        SUBROUTINE Prepare
            !This subroutines reads the input.dat file, and sets various parameters
            ! ===================
            
            !local variable declarations
            double precision UROUN, COMP
            double precision slopetarfactor
            !end variable declaration
            
            ! ====Initialize and Bcast the Chemkin and transport code=========!!
            OPEN(LINKCK, FORM='UNFORMATTED', STATUS='UNKNOWN',FILE=CHEMfile)
            OPEN(LINKTP, FORM='UNFORMATTED', STATUS='UNKNOWN',FILE=TRANfile)
            CALL CKLEN(LINKCK, LOUT, LENICK, LENRCK, LENCCK)

            ! RDB info mecanismo e inicializacion arreglos
            if (MYID .EQ. 0) print*,"LENICK",LENICK, "LENRCK",LENRCK, "LENCCK",LENCCK
            CALL MCLEN(LINKTP, LOUT, LENIMC, LENRMC)
            if (MYID .EQ. 0) print*, "LENIMC",LENIMC, "LENRMC",LENRMC

            CWORK = " "     ! agregada para inicializar vacio el arreglo
            IWORK = 0       ! agregada para inicializar  el arreglo
            WORK = 0.d0     ! agregada para inicializar  el arreglo
            CALL CKINIT(LENICK,LENRCK,LENCCK,LINKCK,LOUT,IWORK,WORK,CWORK)
!            if (MYID .EQ. 0) print*, (i,CWORK(i),i=1,LENCWK)
!            if (MYID .EQ. 0) print*, (IWORK(i),i=1,LENIWK)
!            if (MYID .EQ. 0) print*, (i,WORK(i),i=1,LENRWK)
            CLOSE(LINKCK)

            ITPWRK = 0       ! agregada para inicializar  el arreglo
            TPWRK = 0.d0     ! agregada para inicializar  el arreglo
            CALL MCINIT(LINKTP, LOUT, LENIMC, LENRMC, ITPWRK,TPWRK)
!            if (MYID .EQ. 0) print*, (i,ITPWRK(i),i=1,LENITP)
!            if (MYID .EQ. 0) print*, (i,TPWRK(i),i=1,LENTP)
            CLOSE(LINKTP)

            ! -----Calculate the perturbation number-----------------
            if (MYID .EQ. 0) THEN
                UROUN=1.d0
                1234    CONTINUE
                UROUN=UROUN*0.5d0
                COMP=1.d0+UROUN
                IF(COMP.NE.1.d0) GOTO 1234
!                ABSO=SQRT(2.d0*UROUN)
                ABSO=DSQRT(2.d0*UROUN)        !RDB double
!                print*,'ABSO', ABSO, ABSO**2, ABSO**5
            endif

            !broadcast perturbation number to all processors
            CALL MPI_BCAST(ABSO,1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)


            !   OPEN AND READ AN INPUT PARAMETER CARD
            IF(MYID.EQ.0)THEN
                OPEN(UNIT=214,STATUS='old',FORM='FORMATTED',FILE='input.dat')
                WRITE(6,*) 'Reading Values from the Input Card...'
                read(214,*)DT
                read(214,*)LSOLVE(1)
                read(214,*)LSOLVE(2)
                read(214,*)LSOLVE(3)
                read(214,*)LSOLVE(4)
                read(214,*)LSOLVE(5)
                read(214,*)LSOLVE(6)
                read(214,*)LAST
                read(214,*)saveint
                read(214,*)sootdetail
                read(214,*)MECHANISM
                read(214,*)SOOTMODEL
                read(214,*)ALPHA_SURF
                read(214,*)MSection
                read(214,*)FS
                read(214,*)PRESSURE
                read(214,*)L3
                read(214,*)M3
                read(214,*)DTmax
                read(214,*)Miniter
                read(214,*)MINSLOPETAR
                read(214,*)SLOPETARFACTOR
                read(214,*)OUTBOUND
                Write(6,*) 'Values from the Input Card:'
                Write(6,*)  'DT =',DT
                Write(6,*)  'LSOLVE(U)       =',LSOLVE(1)
                Write(6,*)  'LSOLVE(V)       =',LSOLVE(2)
                Write(6,*)  'LSOLVE(PC)      =',LSOLVE(3)
                Write(6,*)  'LSOLVE(Species) =',LSOLVE(4)
                Write(6,*)  'LSOLVE(Soot)    =',LSOLVE(5)
                Write(6,*)  'LSOLVE(Temp)    =',LSOLVE(6)
                Write(6,*)  'LAST =',LAST

                Write(6,*)

                if(SootModel.eq.3) then
                    write(6,*)'You are using the 2-equation soot model with benzene inception'
                elseif(SootModel.eq.2) then
                    write(6,*)'You are using the 2-equation soot model without benzene inception'
                elseif(SootModel.eq.1) then
                    write(6,*)'You are using the sectional soot model'
                    if(Alpha_surf.EQ.5) then
                        write(6,*)'You are using a functional form of alpha_surf'
                    else
                        write(6,*)'You are using a constant value of alpha_surf, which is ',alpha_surf
                    endif
                endif
                if(Mechanism.LE.3) then
                    write(6,*)'Condensation and nucleation based in heaviest PAH (BAPYR/BAPYR*S/BGHIF)'
                else
                    write(6,*)'Condensation and nucleation based in Pyrene (A4)'
                endif
                if(MYID.EQ.0) write(6,*) ' '

                !fix this
!               !if you try to use the section code with a non-PAH mechanism
!               if(SootModel.EQ.1.AND.MECHANISM.    )then
!                   write(6,*)'This mechanism is not compatible with the section model. Please select a valid mechanism'
!                   close(214)
!                   STOP
!               endif
!
!               !if you try to use the 2-eq Wooley(Benzene) code with a mechanism without benzene
!               if(SootModel.EQ.3.AND.MECHANISM.   )then !fix this
!                   write(6,*)'This mechanism is not compatible with the Wooley 2-eq (Benzene Inception) model. Please select a valid mechanism'
!                   close(214)
!                   STOP
!               endif

                Write(6,*)
                close(214)

                if(SootModel.eq.2.OR.SootModel.EQ.3) MSection = 1 !if using 2-equation model, automatically set the number of sections to 1
            endif

            call MPI_BARRIER(MPI_COMM_WORLD,IERR)

            !   BROADCAST THE INPUT CARD VARIABLES TO THE OTHER CPUS
!
            Call MPI_BCAST(DT,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(DTmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(LAST,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(saveint,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(LSOLVE(1),6,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(MECHANISM,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(SOOTMODEL,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(MSection,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(FS,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(PRESSURE,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(ALPHA_SURF,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(L3,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(M3,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(Miniter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(MINSLOPETAR,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(SLOPETARFACTOR,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(OUTBOUND,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(SOOTDETAIL,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

            call MPI_BARRIER(MPI_COMM_WORLD,IERR)

            !Set up species indicies, and PAH indicies, based on the used mechanism
            include 'mechanisms.inc'

            if(MYID.EQ.0) write(6,*)'You are using ',CHEMfile,' chemical kinetic library'
            if(MYID.EQ.0) write(6,*) 'KK =', KK
            if(MYID.EQ.0) write(6,*) ' '

            !set up initial slope_T target for auto-timestepping
            slopetar=slopetarfactor/(10.d0**(dlog10(DT/1.d-6)/dlog10(2.0d0)))    !slope target for auto-timestepping
            if((slopetar).LE.minslopetar) then
                slopetar = minslopetar                        !if target is too low, set it to the minimum
            endif

            !set dtiter to zero
            dtiter = 0

        end subroutine Prepare
        !*****************************************************************************************************************

        !***************************************************************************************************************
        SUBROUTINE ConstructMesh
        !This subroutine constructs the grid

            !local variable declarations
            double precision DX1, SHRINKSTART, SHRINKEND, SHRINKFAC, XEXPANDSTART, XSTRETCHFAC, DX
            double precision DY1, YSTRETCHFAC, YEXPANDSTART
            DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: XUAll
            !end local variable declarations

            !allocate arrays
            allocate (XUAll(0:NUMP*L3+4))

            !read in the meshing parameters
            IF(MYID==0)THEN
                OPEN(UNIT=215,STATUS='old',FORM='FORMATTED',FILE='meshing.dat')
                WRITE(6,*) 'Reading the Meshing Parameters...'

                read(215,*)DX1
                read(215,*)SHRINKSTART
                read(215,*)SHRINKEND
                read(215,*)SHRINKFAC
                read(215,*)XEXPANDSTART
                read(215,*)XSTRETCHFAC
                read(215,*)DY1
                read(215,*)YSTRETCHFAC
                read(215,*)YEXPANDSTART
                read(215,*)TUBEEND
                close(215)

            ENDIF

            !broadcast values for the radial (or y) direction mesh
            Call MPI_BCAST(DY1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(YSTRETCHFAC,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(YEXPANDSTART,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(Tubeend,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

            call MPI_BARRIER(MPI_COMM_WORLD,IERR)

            !building the grid in the x (or axial) direction**************************************


            if(MYID==0)then
                !calculating U CV centers on CPU-0
!               XUall(0)=0.0
!               XUall(1)=0.0
!               XUall(2)=0.0
!                DX=DX1*3.5d0
!                DO I=3, NUMP*L3+2   !NUMP*L3+1 is the actual end of the domain
!                    if(SHRINKSTART.LT.0.0.AND.(XUall(I-1)+DX1).LT.XEXPANDSTART) then
!                        XUall(I)=XUall(I-1)+DX1
!                    elseif((XUall(I-1)+DX1).LT.SHRINKSTART) then
!                        DX = DX1*3.5d0
!                        XUall(I)=XUall(I-1)+DX
!                    elseif((XUall(I-1)+DX1).GT.SHRINKSTART.AND.(XUall(I-1)+DX1).LT.SHRINKEND) then
!                        DX = DX*SHRINKFAC
!                        XUall(I)=XUall(I-1)+DX  
!                    elseif((XUall(I-1)+DX1).GT.SHRINKEND.AND.(XUall(I-1)+DX1).LT.XEXPANDSTART) then
!                        DX = DX1
!                        XUall(I)=XUall(I-1)+DX
!                    elseif((XUall(I-1)+DX1).GT.XEXPANDSTART) then
!                        DX = DX*XSTRETCHFAC
!                        XUall(I)= XUall(I-1)+DX
!                    endif                 
!                ENDDO            
!                        
!                XUall(NUMP*L3+3) = XUall(NUMP*L3+2)
!                XUall(NUMP*L3+4) = XUall(NUMP*L3+2) 

                XUall(0)=0.0d0
                XUall(1)=0.0d0
                XUall(2)=0.0d0
                DX=DX1*3.5d0
                DO I=3, NUMP*L3+2   !NUMP*L3+1 is the actual end of the domain
                    if((XUall(I-1)+DX).LT.SHRINKSTART.AND.TUBEEND.GT.0.0d0) then
                        XUall(I)=XUall(I-1)+DX
                    elseif((XUall(I-1)+DX).GT.SHRINKSTART.AND.(XUall(I-1)+DX).LT.SHRINKEND.AND.TUBEEND.GT.0.0d0) then
                        DX = DX*SHRINKFAC
                        XUall(I)=XUall(I-1)+DX  
                    elseif(((XUall(I-1)+DX).GT.SHRINKEND.OR.TUBEEND.LE.0.0d0).AND.(XUall(I-1)+DX).LT.XEXPANDSTART) then
                        DX = DX1
                        XUall(I)=XUall(I-1)+DX
                    elseif((XUall(I-1)+DX).GT.XEXPANDSTART) then
                        DX = DX*XSTRETCHFAC
                        XUall(I)= XUall(I-1)+DX
                    endif                 
                ENDDO             
                        
                XUall(NUMP*L3+3) = XUall(NUMP*L3+2)
                XUall(NUMP*L3+4) = XUall(NUMP*L3+2) 
                
            endif 
            
            !distribute the resulting U CV centers to the individual CPUs    
            DO indp=Nump-1,1,-1
                
                IF (MYID.EQ.0) THEN 
                           
                    DO I=0,LM1         
                        XU(I) = XUall(L3*(indp)+I)
                    ENDDO
                    
                    CALL MPI_SEND(XU(0),LM2,MPI_DOUBLE_PRECISION,indp,22691120+indp,MPI_COMM_WORLD,IERR)
               
                ELSEIF(MYID.EQ.indp) THEN
                    
                    CALL MPI_RECV(XU(0),LM2,MPI_DOUBLE_PRECISION,0,22691120+indp,MPI_COMM_WORLD,STATUS,IERR)
                    
                ENDIF
                
            ENDDO
            
            !set the U CV centers for CPU-0
            IF (MYID.EQ.0) THEN    
                XU(0:LM1) = XUall(0:LM1)
            ENDIF     
            
            !calculate CV centers for the V, and P, CVs in the axial direction
            X(0:L0)=0.5d0*(XU(1:L0+1)+XU(0:L0))

            !write out the number of CVs in the x (r)-dir, and the length 
            !if(MYID==NUMP-1)then
            !    write(6,*)'In the x-dir there are ',(NUMP*L3),'CVs.'
            !    write(6,*)'The Grid Extends to  ',XU(L1)
            !    write(6,*)
            !endif

            !calculate the spacings for the CV centers  
            XDIF(1:L0)=X(1:L0)-X(0:L0-1)
            ! --------
            
            !calculate the height (or z-direction length) of each V,P CV
            XCV(0:L0)=XU(1:L0+1)-XU(0:L0)
            ! --------

            if(myid.eq.0) then

                !calculate the height (or z-direction length) of each U CV
                XCVS(3:L0)=XDIF(3:L0)
                XCVS(3)=XCVS(3)+XDIF(2)  !this is to account for the fact that the first CV on the 
                                         !staggered grid is actually 1.5*DX long.
                xcvs(2)=0.0d0
                xcvs(1)=0.0d0
                
                !calculate the distance from each node on the U grid to the upper and lower face
                !of the associated control volume
                XCVI(3:L0)=0.5d0*XCV(3:L0)   !upper face distance
                XCVIP(3:L0)=XCVI(3:L0)     !lower face distance
                XCVIP(2)=XCV(2)  !this is to account for the fact that the first CV on the 
                                  !staggered grid is actually 1.5*DX long, thus making the distance to
                                  !its lower face DX, and not 0.5*DX
                xcvi(2)=0.0d0
                xcvip(1)=0.0d0
                xcvi(1)=0.0d0
                xcvip(0)=0.0d0
                xcvi(0)=0.0d0

            else if(myid.eq. (NUMP-1)) then

                !calculate the height (or z-direction length) of each U CV
                XCVS(1:L2)=XDIF(1:L2)
                XCVS(L2)=XCVS(L2)+XDIF(L1) !this is to account for the fact that the last CV on the 
                                            !staggered grid is actually 1.5*DX long.
                xcvs(L1)=0.0d0
                xcvs(L0)=0.0d0
                
                !calculate the distance from each node on the U grid to the upper and lower face
                !of the associated control volume
                XCVI(0:L3)=0.5d0*XCV(0:L3)   !upper face distance
                XCVIP(0:L3)=XCVI(0:L3)     !lower face distance
                XCVI(L2)=XCV(L2)    !this is to account for the fact that the last CV on the 
                                  !staggered grid is actually 1.5*DX long, thus making the distance to
                                  !its upper face DX, and not 0.5*DX
                xcvip(L2)=0.0d0
                xcvi(L1)=0.0d0
                xcvip(L1)=0.0d0
                xcvi(L0)=0.0d0
                xcvip(L0)=0.0d0

            else
                !calculate the height (or z-direction length) of each U CV
                XCVS(1:L0)=XDIF(1:L0)

                !calculate the distance from each node on the U grid to the upper and lower face
                !of the associated control volume
                XCVI(0:L0)=0.5d0*XCV(0:L0)
                XCVIP(0:L0)=XCVI(0:L0)
            ENDif
            ! ---------------------------------------------------------
            
            !not sure
            if(myid .ne.0 .and. myid.ne. (NUMP-1)) then
                FX(1:L0)=0.5d0*XCV(0:L0-1)/XDIF(1:L0)
                FXM(1:L0)=1.0d0-FX(1:L0)

            else if(myid.eq.0) then
                FX(3:L0)=0.5d0*XCV(2:L0-1)/XDIF(3:L0)
                FXM(3:L0)=1.0d0-FX(3:L0)

                FX(2)=0.0d0
                FXM(2)=1.0d0
                fx(1)=0.0d0
                fxm(1)=0.0d0

            else
                FX(1:L2)=0.5d0*XCV(0:L2-1)/XDIF(1:L2)
                FXM(1:L2)=1.0d0-FX(1:L2)

                FX(L1)=1.0d0
                FXM(L1)=0.0d0
                fx(L0)=0.0d0
                fxm(L0)=0.0d0
            endif
            ! ! ! ! ! !-----------X (z-direction, axial direction) DONE!------------       
  
            !building grid in the y (r, or radial) direction
            YV(1)=0.0d0
            YV(2)=0.0d0
            R(1)=0.0d0        ! RDB correccion Pipe

            !construct V grid (up to node 50, no stretching)
            DO J=3,YEXPANDSTART
                IF(YV(J-1).GT.0.45d0.AND.YV(J-1).LT.0.66d0) THEN
                    YV(J)=YV(J-1)+DY1/1.3d0
                ELSE
                    YV(J)=YV(J-1)+DY1
                ENDIF
            ENDDO
            DO J=YEXPANDSTART+1,M1
                YV(J)=YV(J-1)+DY1*YSTRETCHFAC**(J-YEXPANDSTART)
            ENDDO

            !construct the U,P grid
            Y(1)=YV(2)
            
            Y(2:M2)=0.5d0*(YV(3:M2+1)+YV(2:M2))
            
            Y(M1)=YV(M1)

            !calculate grid spacings
            YDIF(2:M1)=Y(2:M1)-Y(1:M1-1)

            !calculate length (or r-direction distance) of each U,P CV
            YCV(2:M2)=YV(3:M2+1)-YV(2:M2)

            !calculate length (or r-direction distance) of each V CV
            YCVS(3:M2)=YDIF(3:M2)
            
            YCVS(3)=YCVS(3)+YDIF(2)
            YCVS(M2)=YCVS(M2)+YDIF(M1)

            !r-co-ordinates for the U,P nodes. R(1) = 0 was set earlier
            !these values should be the same as Y(J)
            DO J = 2,M1
                R(J)=R(J-1)+YDIF(J)
            ENDDO    
            
            !r-co-ordinates for the V nodes.
            !these values should be the same as YV(J)    
            RMN(2)=R(1)
            
            DO J=3,M2
                RMN(J)=RMN(J-1)+YCV(J-1)
            ENDDO      
            
            RMN(M1)=R(M1)
            
            !Calculation of r*deltar term for P,U control volumes (deltar is the CV length, r is node location)
            YCVR(2:M2)=R(2:M2)*YCV(2:M2)
            ARX(2:M2)=YCVR(2:M2)

            !Calculation of r*deltar term for V control volumes (deltar is the CV length, r is node location)
            YCVRS(4:M3)=0.5d0*(R(4:M3)+R(3:M3-1))*YDIF(4:M3)
            
            YCVRS(3)=0.5d0*(R(3)+R(1))*YCVS(3)
            YCVRS(M2)=0.5d0*(R(M1)+R(M3))*YCVS(M2)    
         
            ARXJ(3:M3)=0.25d0*(1.0d0+RMN(3:M3)/R(3:M3))*ARX(3:M3)
            ARXJP(3:M3)=ARX(3:M3)-ARXJ(3:M3) 

            ARXJP(2)=ARX(2) ! ALL
            ARXJ(M2)=ARX(M2)


            FV(3:M3)=ARXJP(3:M3)/ARX(3:M3)
            FVP(3:M3)=1.0d0-FV(3:M3)
            
            !
            FY(3:M2)=0.5d0*YCV(2:M2-1)/YDIF(3:M2)
            FYM(3:M2)=1.0d0-FY(3:M2)

            FY(2)=0.d0
            FYM(2)=1.0d0
            FY(M1)=1.0d0
            FYM(M1)=0.d0
            !-----Y DONE------------
            
            !write out the number of CVs in the x (r)-dir, and the length
            if(MYID==NUMP-1)then
                write(6,*)'In the x-dir there are ',(NUMP*L3),'CVs.'
                write(6,*)'The Grid Extends to  x=',XU(L1)
                write(6,*)
                write(6,*)'In the r-dir there are ',(M3),'CVs.'
                write(6,*)'The Grid Extends to  r=',YV(M1)
                write(6,*)
            endif

            DEallocate(XUALL)
      
        end subroutine ConstructMesh
        !**********************************************************************************************
        
        !**********************************************************************************************
        SUBROUTINE SetInitial
        !this subroutine sets the initial guess, and flags CVs in the solid reion
        
            !local variable declarations
            double precision :: uair, ufuel, tair, tfuel
            double precision :: taperstart, roinit, riinit, roend, riend, ropv, ripv, rou, riu
            integer :: irest
            integer :: fuel1, fuel2, fuel3, fuel4,fuel5,fuel6, ox1, ox2, ox3
            double precision fuel1frac, fuel2frac, fuel3frac,fuel4frac,fuel5frac,fuel6frac, ox1frac, ox2frac, ox3frac
            double precision, dimension(:), allocatable :: SMassF, SMoleF
            double precision, dimension(:), allocatable :: Xrestart, Yrestart
            double precision, dimension(:,:), allocatable :: UTrestart, VTrestart, Trestart, Prestart, PCrestart, RHOrestart
            double precision, dimension(:,:), allocatable :: U_tot, V_tot  
            double precision, dimension(:,:,:), allocatable :: SpeciesRestart, SootSecRestart
            double precision :: Xlimit      ! para modificar rest.dat
            integer :: itoto, jtoto
            !done local variable declarations
            
            !allocate arrays
            allocate(SMassF(KK))
            allocate(SMoleF(KK))
            allocate(Xrestart(L3*Nump+2))
            allocate(Yrestart(M1))
            allocate(UTrestart(L3*Nump+2,M1))
            allocate(V_tot(L3*Nump+3,M1))
            allocate(U_tot(0:L3*Nump+3,M1))
            allocate(VTrestart(L3*Nump+2,M1))
            allocate(Trestart(L3*Nump+2,M1))
            allocate(Prestart(L3*Nump+2,M1))
            allocate(PCrestart(L3*Nump+2,M1))
            allocate(RHOrestart(L3*Nump+2,M1))
            allocate(SpeciesRestart(L3*Nump+2,M1,KK))
            allocate(SootSecRestart(L3*Nump+2,M1,2*MSection)) 

            !read in values from the initialvalues.dat file
            IF(MYID==0)THEN
                OPEN(UNIT=216,STATUS='old',FORM='FORMATTED',FILE='initialvalues.dat')
                WRITE(6,*) 'Reading the Initial Values Parameters...'

                read(216,*)IREST
                read(216,*)TAPERSTART
                read(216,*)ROINIT
                read(216,*)RIINIT
                read(216,*)ROEND
                read(216,*)RIEND
                read(216,*)uair
                read(216,*)ufuel
                read(216,*)tair
                read(216,*)tfuel
                read(216,*)fuel1
                read(216,*)fuel1frac
                read(216,*)fuel2
                read(216,*)fuel2frac
                read(216,*)fuel3
                read(216,*)fuel3frac
                read(216,*)fuel4
                read(216,*)fuel4frac
                read(216,*)fuel5
                read(216,*)fuel5frac
                read(216,*)fuel6
                read(216,*)fuel6frac
                read(216,*)ox1
                read(216,*)ox1frac
                read(216,*)ox2
                read(216,*)ox2frac
                read(216,*)ox3
                read(216,*)ox3frac
                read(216,*)printinitial
                close(216)

                if(IREST==0)then
                    Write(6,*)'This calculation is starting from scratch.'
                elseif(IREST==1)then
                    Write(6,*)'This calculation is continuing from a restart file.'
                else
                    Write(6,*)'ERROR - IREST not recognized. Program will terminate'
                    close(214)
                    STOP
                endif

            ENDIF

            call MPI_BARRIER(MPI_COMM_WORLD,IERR)

            !broadcast values from the initial values input card
            Call MPI_BCAST(TAPERSTART,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(ROINIT,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(RIINIT,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(ROEND,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(RIEND,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(Tfuel,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(TAir,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(Ufuel,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(UAir,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(IREST,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(printinitial,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(ox1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(ox2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(ox3,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(ox1frac,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(ox2frac,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(ox3frac,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(fuel1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(fuel2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(fuel3,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(fuel4,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(fuel5,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(fuel6,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(fuel1frac,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(fuel2frac,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(fuel3frac,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(fuel4frac,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(fuel5frac,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(fuel6frac,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

            call MPI_BARRIER(MPI_COMM_WORLD,IERR)
            
            !initially set that all CVs are not in a solid region
            SOLIDCVP = .false.
            SOLIDCVU = .false.
            SOLIDCVV = .false.

            if(IREST.EQ.0) then  !if starting from scratch, set arbitrary initial guess
                DO J=1,M1
                    DO I=0,L0
                    
                        if(Y(J).LT.RIEND) then !if the radius is less than the inner fuel tube radius, set velocity to that of the fuel
                            U(I,J)= Ufuel   !(Ufuel+Uair)/2.d0   ! Ufuel
                        else
                            U(I,J)=Uair  !otherwise set velocity to that of the air inlet
                        endif

                        !set species mole fractions to that of the oxidizer everywhere
                       ! if(X(i).GT.TUBEEND.OR.Y(J).GT.ROEND) then
                        if(X(i).GT.TUBEEND.OR.Y(J).GT.RIEND) then       !RDB
                            SpeciesMF(I,J,ox1)=ox1frac
                            SpeciesMF(I,J,ox2)=ox2frac
                            SpeciesMF(I,J,ox3)=ox3frac
                        endif

                        !except for locations that are inside the fuel tube
                       ! if(X(i).LE.TUBEEND.AND.Y(J).LE.ROEND) then
                        if(X(i).LE.TUBEEND.AND.Y(J).LE.RIEND) then      !RDB
                            SpeciesMF(I,J,fuel1)=fuel1frac
                            SpeciesMF(I,J,fuel2)=fuel2frac 
                            SpeciesMF(I,J,fuel3)=fuel3frac
                            SpeciesMF(I,J,fuel4)=fuel4frac   
                            SpeciesMF(I,J,fuel5)=fuel5frac   
                            SpeciesMF(I,J,fuel6)=fuel6frac                          
                        endif
                        
                        !set Temperature in the domain
                        if(X(i).LT.TUBEEND) then
                            T(I,J)= Tair    !Tfuel
                            if(Y(j).GE.ROINIT) T(I,J)= Tair
                        else !OBS7
                            T(I,J)= 1900.d0
!                            T(I,J)= Tair
                            if(Y(J).GT.ROINIT) T(I,J) = Tair
                            if(X(I).GT.0.3d0)  T(I,J) = Tair

!                            if((X(I).LT.0.3d0).and.(Y(J).LT.ROINIT).and.(J.NE.1))  then
!                                T(I,J) = Tair
!                                SpeciesMF(I,J,:)= 0.0d0
!                                SpeciesMF(I,J,fuel1)=fuel1frac
!                                SpeciesMF(I,J,fuel2)=fuel2frac
!                                SpeciesMF(I,J,fuel3)=fuel3frac
!                                SpeciesMF(I,J,fuel4)=fuel4frac
!                                SpeciesMF(I,J,fuel5)=fuel5frac
!                                SpeciesMF(I,J,fuel6)=fuel6frac
!                             endif
                        endif

                        !convert mol fractions to mass fractions
                        SMoleF(:) = SpeciesMF(I,J,:)
                        CALL CKXTY(SMoleF,IWORK,WORK,SMassF)
                        SpeciesMF(I,J,:) = SMassF(:)
                    ENDDO
                ENDDO
                
                !set the first iteration to zero
                ITER = 0
                FIRST = 0
            endif   
            
            call MPI_BARRIER(MPI_COMM_WORLD,IERR)     

            if(IREST.EQ.1) then !if starting from a restart file and continuing with the same grid

                IF (MYID.EQ.0) THEN
                    !read in the solution from the restart file. Note the species data is in mole, not mass, fraction
                    OPEN(UNIT=15320,STATUS='unknown',FORM='FORMATTED',FILE='rest.dat')
                    read(15320,*) ITER
                    DO I=1,L3*NUMP+2    ! L3*NUMP+2 = 3*64+2 = 194
                        DO J=1,M1
                            READ (15320,*) Yrestart(J),Xrestart(I),UTrestart(I,J),VTrestart(I,J),Trestart(I,J),Prestart(I,J), &
                                           PCrestart(I,J),RHOrestart(I,J),(SpeciesRestart(i,j,k),k=1,kk),                     &
                                           (SootSecRestart(i,j,k),k=1,MSection*2)
                        END DO
                    END DO
                    CLOSE(15320)

                    327 FORMAT(1X,1P172E26.15)  !format of restart file. should match the format of the SOLN file format being used for restart
!                    327 FORMAT(1X,1P309E26.15)  !format of restart file. should match the format of the SOLN file format being used for restart
                    !done reading in the solution

                    ! sobreescribir secciones de hollin
!                    do k=20,MSection       !1,MSection
!                     DO I=1,L3*NUMP+2    ! L3*NUMP+2 = 3*64+2 = 194
!!                      Xlimit = 6.5d0    !5.3d0
!!                      !if(k .GE. 31) Xlimit = 5.3d0
!!                      if(Xrestart(i).GT.Xlimit)then
!                        DO J=1,M1
!                    SootSecRestart(i,j,k) = SootSecRestart(i,j,k) * 1.d-3   !20       ! # agregados
!                    SootSecRestart(i,j,k+MSection) = SootSecRestart(i,j,k+MSection) * 1.d-4 !20 !# particulas primarias
!                        ENDDO
!!                      endif
!                     ENDDO
!                    enddo

                    ! disminuir temperatura
!                    DO I=1,L3*NUMP+2
!                        DO J=1,M1
!                            if(Trestart(I,J).GT.2700.d0) Trestart(I,J)=2700.0d0
!                            !if(Prestart(I,J).GT.2.d0) Prestart(I,J)=2.0d0
!                        ENDDO
!                    ENDDO

                    ! reiniciar campo presion
                    !Prestart(:,:) = 0.0d0

                    ! especies borde quemador
!                    DO J = 1,M1
!                        if(Y(J).GE.RIEND .AND. Y(J).LE.ROEND) SpeciesRestart(1,j,:) = SpeciesRestart(1,33,:)
!                    enddo

                    !   encendido tardio
!                    Trestart(4,3:7) = 1000.0d0       !RDB

                    !determine the velocities at the center of the U CVs, not the center of the P CVs
                    DO J=1,M1
                        U_tot(1,J) = UTrestart(1,J)
                        U_tot(0,J) = U_tot(1,J)

                        DO I=1,L3*NUMP+1
                            U_tot(I+1,J) = 2.0d0*UTrestart(I,J)- U_tot(I,J)
                        END DO
                    END DO                  
                    
                    !determine the velocities at the center of the V CVs, not the center of the P CVs
                    !also convert species mole fractions to mass fractions
                    V_tot = 0.d0
                    DO I=1,L3*NUMP+2
                        DO J=1,M1
                            IF(J.NE.M1) THEN
                                V_tot(I,J+1) = 2.0d0*VTrestart(I,J)- V_tot(I,J)
                            ENDIF
                            SMoleF(:) = SpeciesRestart(I,J,:)
                            CALL CKXTY(SMoleF,IWORK,WORK,SMassF)
                            SpeciesRestart(I,J,:) = SMassF(:)
                        END DO
                    END DO

                ENDIF    ! IF (MYID.EQ.0)

                call MPI_BARRIER(MPI_COMM_WORLD,IERR)

  !          Prestart(:,:)=1.d0      ! sobreescribir la presion de la solucion anterior
  !          PCrestart(:,:)=0.d0     ! sobreescribir la correccion de la presion de la solucion anterior
                ! ---------- distribute solution to internal CPUs ------------------------------
                DO indp=1,(NUMP-1)  !index of the internal CPUs


                    IF (MYID.EQ.0) THEN ! CPU0 distributes the solution

                        if (indp.eq.(NUMP-1)) then  ! the last cpu
                            DO J=1,M1
                                DO I=0,L1
                                    U(I,J)=U_tot(L3*indp+I,J)
                                    V(I,J)=V_tot(L3*indp+I,J)

                                    T(I,J)=Trestart(L3*indp+I,J)
                                    RHO(I,J)=RHOrestart(L3*indp+I,J)
                                    P(I,J)=Prestart(L3*indp+I,J)
                                    PC(I,J)=PCrestart(L3*indp+I,J)

                                    SpeciesMF(I,J,:)= SpeciesRestart(L3*indp+I,J,:)
                                    SootSec(I,J,:)= SootSecRestart(L3*indp+I,J,:)

                                ENDDO
                                U(L0,J)= U_tot(L3*indp+L0,J)    ! L0=L3+3 =6
                                U(LM1,J)= U(L0,J)   ! LM1=L3+4=7

                                V(L0,J)= V_tot(L3*indp+L0,J)    ! L0=L3+3 =6
                                V(LM1,J)= V(L0,J)   ! LM1=L3+4=7

                                T(L0,J)= T(L1,J)    !T
                                T(LM1,J)= T(L0,J)

                                RHO(L0,J)= RHO(L1,J)    !density
                                RHO(LM1,J)= RHO(L0,J)

                                P(L0,J)= P(L1,J)    !Pressure
                                P(LM1,J)= P(L0,J)

                                PC(L0,J)= PC(L1,J)  !Pressure
                                PC(LM1,J)= PC(L0,J)

                                SpeciesMF(L0,J,:)= SpeciesMF(L1,J,:)
                                SpeciesMF(LM1,J,:)= SpeciesMF(L0,J,:)

                                SootSec(L0,J,:)= SootSec(L1,J,:)
                                SootSec(LM1,J,:)= SootSec(L0,J,:)

                            ENDDO

                        else ! other cpus (not the last one)

                            DO I=0,LM1  ! 30 May 2011, LM1 =L3+4=7, L3*indp+i-1
                                DO J=1,M1

                                    U(I,J)=U_tot(L3*indp+I,J)
                                    V(I,J)=V_tot(L3*indp+I,J)
                                    T(I,J)=Trestart(L3*indp+I,J)
                                    RHO(I,J)=RHOrestart(L3*indp+I,J)
                                    P(I,J)=Prestart(L3*indp+I,J)
                                    PC(I,J)=PCrestart(L3*indp+I,J)
                                    SpeciesMF(I,J,:)= SpeciesRestart(L3*indp+I,J,:)
                                    SootSec(I,J,:)= SootSecRestart(L3*indp+I,J,:)

                                ENDDO
                            ENDDO

                        endif ! end of if (indp.eq.(NUMP-1)) then

                        CALL MPI_SEND(U(0,1), LM2*M1,MPI_DOUBLE_PRECISION,indp,2269110+indp,MPI_COMM_WORLD,IERR)    ! changed 20 May
                        CALL MPI_SEND(V(0,1), LM2*M1,MPI_DOUBLE_PRECISION,indp,3279110+indp,MPI_COMM_WORLD,IERR)    ! changed 20 May
                        CALL MPI_SEND(RHO(0,1), LM2*M1,MPI_DOUBLE_PRECISION,indp,4289110+indp,MPI_COMM_WORLD,IERR)
                        CALL MPI_SEND(T(0,1), LM2*M1,MPI_DOUBLE_PRECISION,indp,4299110+indp,MPI_COMM_WORLD,IERR)
                        CALL MPI_SEND(SpeciesMF(0,1,1), LM2*M1*KK,MPI_DOUBLE_PRECISION,indp,4309110+indp,MPI_COMM_WORLD,IERR)
                        CALL MPI_SEND(SootSec(0,1,1), LM2*M1*MSection*2,MPI_DOUBLE_PRECISION,indp,4319110+indp,MPI_COMM_WORLD,IERR)
                        CALL MPI_SEND(P(0,1), LM2*M1,MPI_DOUBLE_PRECISION,indp,4329110+indp,MPI_COMM_WORLD,IERR)
                        CALL MPI_SEND(PC(0,1), LM2*M1,MPI_DOUBLE_PRECISION,indp,4339110+indp,MPI_COMM_WORLD,IERR)

                    ELSEIF(MYID.eq.indp) then

                        CALL MPI_RECV(U(0,1), LM2*M1,MPI_DOUBLE_PRECISION,0,2269110+indp,MPI_COMM_WORLD,STATUS,IERR)    ! changed 20 May
                        CALL MPI_RECV(V(0,1), LM2*M1,MPI_DOUBLE_PRECISION,0,3279110+indp,MPI_COMM_WORLD,STATUS,IERR)    ! changed 20 May

                        CALL MPI_RECV(RHO(0,1), LM2*M1,MPI_DOUBLE_PRECISION,0,4289110+indp,MPI_COMM_WORLD,STATUS,IERR)
                        CALL MPI_RECV(T(0,1), LM2*M1,MPI_DOUBLE_PRECISION,0,4299110+indp,MPI_COMM_WORLD,STATUS,IERR)
	CALL MPI_RECV(SpeciesMF(0,1,1), LM2*M1*KK,MPI_DOUBLE_PRECISION,0,4309110+indp,MPI_COMM_WORLD,STATUS,IERR)
	CALL MPI_RECV(SootSec(0,1,1), LM2*M1*MSection*2,MPI_DOUBLE_PRECISION,0,4319110+indp,MPI_COMM_WORLD,STATUS,IERR)
	CALL MPI_RECV(P(0,1), LM2*M1,MPI_DOUBLE_PRECISION,0,4329110+indp,MPI_COMM_WORLD,STATUS,IERR)
                        CALL MPI_RECV(PC(0,1), LM2*M1,MPI_DOUBLE_PRECISION,0,4339110+indp,MPI_COMM_WORLD,STATUS,IERR)

                    END IF
                ENDDO   !DO indp=1,(NUMP-2)  !index of process

                !--------- the first cpu: CPU0 ---------------------------
                IF (MYID.EQ.0) THEN
                    DO J=1,M1
                        U(0,J) = U_tot(0,J)
                        U(1,J) = U_tot(1,J)
                        U(2,J) = U_tot(2,J)

                        V(0,J) = 0.d0
                        V(1,J) = V_tot(1,J)
                        V(2,J) = V_tot(2,J)

                        T(1,J)=Trestart(1,J)
                        T(0,J)=T(1,J)

                        RHO(1,J)=RHOrestart(1,J)
                        RHO(0,J)=RHO(1,J)

                        P(1,J)=Prestart(1,J)
                        P(0,J)=P(1,J)

                        PC(1,J)=PCrestart(1,J)
                        PC(0,J)=PC(1,J)

                        SpeciesMF(1,J,:)= SpeciesRestart(1,J,:)
                        SpeciesMF(0,J,:)= SpeciesMF(1,J,:)

                        SootSec(1,J,:)= SootSecRestart(1,J,:)
                        SootSec(0,J,:)= SootSec(1,J,:)

                        DO I=2,LM1  ! LM1=7
                            U(I,J)= U_tot(I,J)
                            V(I,J)= V_tot(I,J)

                            T(I,J)=Trestart(I,J)
                            RHO(I,J)=RHOrestart(I,J)
                            P(I,J)=Prestart(I,J)
                            PC(I,J)=PCrestart(I,J)

                            SpeciesMF(I,J,:)= SpeciesRestart(I,J,:)
                            SootSec(I,J,:)= SootSecRestart(I,J,:)
                        ENDDO

                    ENDDO

                ENDIF   ! end of MYID.eq.0 -----------------------------q

                !make sure no locations of zero species mass fractions
                DO I=0,LM1
                    DO J=1,M1
                    
                        IF(SUM(SpeciesMF(I,J,:)).EQ.0d0.OR.SUM(SpeciesMF(I,J,:)).NE.SUM(SpeciesMF(I,J,:)).OR.&
                                  SpeciesMF(I,J,4).NE.SpeciesMF(I,J,4)) THEN

                            SpeciesMF(I,J,:) = 0.d0
                            IF(Y(J).LT.ROEND) THEN
                                SpeciesMF(I,J,fuel1)=fuel1frac
                                SpeciesMF(I,J,fuel2)=fuel2frac
                                SpeciesMF(I,J,fuel3)=fuel3frac
                                SpeciesMF(I,J,fuel4)=fuel4frac
                                SpeciesMF(I,J,fuel5)=fuel5frac
                                SpeciesMF(I,J,fuel6)=fuel6frac
                            ELSE
                                SpeciesMF(I,J,ox1)=ox1frac
                                SpeciesMF(I,J,ox2)=ox2frac
                                SpeciesMF(I,J,ox3)=ox3frac
                            ENDIF

                            SMoleF(:) = SpeciesMF(I,J,:)
                            CALL CKXTY(SMoleF,IWORK,WORK,SMassF)
                            SpeciesMF(I,J,:) = SMassF(:)

                        ENDIF

                    ENDDO
                ENDDO

                !send out iteration number to all processes
                CALL MPI_BCAST(ITER,1, MPI_integer,0,MPI_COMM_WORLD,IERR)

                !calculate the last iteration to be performed, and record the starting iteration in FIRST
                FIRST = ITER
                LAST=LAST+ITER

            ENDIF   ! IF (IREST.EQ.1)

            !flag CVs that are in the solid region, and set velocities, PC, P, species, and soot to zero if so
            DO J=1,M1
                DO I=0,L0

!                    if(myid.eq.1) THEN
!                        write(6,*)I,J,X(I),Y(J),Tubeend,ROINIT,ROEND,RIEND,RIINIT,TAPERSTART,SOLIDCVP(I,J)
!                    endif

                    !figure out what CVs are in the solid region
                    if(X(i)<TAPERSTART) then
                        ROPV = ROINIT
                        RIPV = RIINIT
                    elseif(X(i)<TUBEEND) then
                        ROPV = ROINIT - (ROINIT-ROEND)*(X(i)-TAPERSTART)/(TUBEEND - TAPERSTART)
                        RIPV = RIINIT + (RIEND-RIINIT)*(X(I)-TAPERSTART)/(TUBEEND - TAPERSTART)
                    else
                        ROPV = ROEND
                        RIPV = RIEND
                    endif
                    if(X(i).LT.TUBEEND.AND.Y(J).LT.ROPV.AND.Y(J).GT.RIPV) SOLIDCVP(I,J) = .true.
                    if(J.NE.1) THEN
                        if(X(i).LT.TUBEEND.AND.Y(J-1).LT.ROPV.AND.Y(J).GT.RIPV) SOLIDCVV(I,J) = .true.
                    endif

                    if(I.NE.0.d0) THEN
                        if(X(i-1)<TAPERSTART) then
                            ROU = ROINIT
                            RIU = RIINIT
                        elseif(X(i-1)<TUBEEND) then
                            ROU = ROINIT - (ROINIT-ROEND)*(X(i-1)-TAPERSTART)/(TUBEEND - TAPERSTART)
                            RIU = RIINIT + (RIEND-RIINIT)*(X(I-1)-TAPERSTART)/(TUBEEND - TAPERSTART)
                        else
                            ROU = ROEND
                            RIU = RIEND
                        endif
                        if(X(i-1).LT.TUBEEND.AND.Y(J).LT.ROU.AND.Y(J).GT.RIU) SOLIDCVU(I,J) = .true.
                    else
                        if(XU(i)<TAPERSTART) then
                            ROU = ROINIT
                            RIU = RIINIT
                        elseif(XU(i)<TUBEEND) then
                            ROU = ROINIT - (ROINIT-ROEND)*(XU(i)-TAPERSTART)/(TUBEEND - TAPERSTART)
                            RIU = RIINIT + (RIEND-RIINIT)*(XU(I)-TAPERSTART)/(TUBEEND - TAPERSTART)
                        else
                            ROU = ROEND
                            RIU = RIEND
                        endif
                        if(XU(i).LT.TUBEEND.AND.Y(J).LT.ROU.AND.Y(J).GT.RIU) SOLIDCVU(I,J) = .true.
                    endif

                    if(SolidCVP(I,J)) then
                        P(I,J) = 0.0d0
                        PC(I,J) = 0.0d0
                        SpeciesMF(I,J,:)=0.0d0
                        SootSec(I,J,:)=0.0d0
                    endif
                    if(SolidCVV(I,J)) V(I,J)=0.0d0
                    if(SolidCVU(I,J)) U(I,J)=0.0d0
                                           
                ENDDO
            ENDDO
            
            call MPI_BARRIER(MPI_COMM_WORLD,IERR)              
            
            ! INLET BOUNDARY
            if (MYID .EQ. 0) then
                !loop over all CVs in the radial direction
                DO J=1,M1
                    SpeciesMF(1,J,:)=0.d0
                    SpeciesMF(0,J,:)=0.d0
                    !for radial locations less than the fuel tube
                    IF(Y(J).LT.RIINIT) THEN
                        !U(2,J)=UFUEL
                        U(2,J)= 2.0d0*Ufuel*(1.-(Y(J)/RIINIT)**2.d0)     ! RDB Parabolic inlet velocity profile
                        U(1,J)=U(2,J)
                        U(0,J)=U(1,J)
                        V(2,J)=0.0d0
                        V(1,J)=0.0d0
                        V(0,J)=0.0d0
                        T(1,J)=TFUEL
                        T(0,J)=T(1,J)
                        SpeciesMF(1,J,fuel1)=fuel1frac
                        SpeciesMF(0,J,fuel1)=SpeciesMF(1,J,fuel1)
                        SpeciesMF(1,J,fuel2)=fuel2frac
                        SpeciesMF(0,J,fuel2)=SpeciesMF(1,J,fuel2) 
                        SpeciesMF(1,J,fuel3)=fuel3frac 
                        SpeciesMF(0,J,fuel3)=SpeciesMF(1,J,fuel3)
                        SpeciesMF(1,J,fuel4)=fuel4frac 
                        SpeciesMF(0,J,fuel4)=SpeciesMF(1,J,fuel4)
                        SpeciesMF(1,J,fuel5)=fuel5frac 
                        SpeciesMF(0,J,fuel5)=SpeciesMF(1,J,fuel5)
                        SpeciesMF(1,J,fuel6)=fuel6frac 
                        SpeciesMF(0,J,fuel6)=SpeciesMF(1,J,fuel6)
                        
                        !convert mole fractions to mass fractions
                        SMoleF(:) = SpeciesMF(0,J,:)
                        CALL CKXTY(SMoleF,IWORK,WORK,SMassF)
                        SpeciesMF(0,J,:) = SMassF(:)
                        SpeciesMF(1,J,:) = SMassF(:)
          
                    
                    !for radial locations that are the surface of the burner lip
                    ELSEIF(Y(J).GE.RIINIT.AND.Y(J).LT.ROINIT) THEN
                        !U(2,J)=0.0d0
                        U(2,J)= Uair/(ROINIT-RIINIT)*(Y(J)-RIINIT)      !RDB variacion lineal de la velocidad sobre el borde
                        U(1,J)=U(2,J)
                        U(0,J)=U(1,J)
                        V(2,J)=0.d0
                        V(1,J)=0.d0
                        V(0,J)=0.d0
                        T(1,J)=TFUEL
                        T(0,J)=T(1,J)
                        
                        !if fuel tube is not modelled, set mass fraction of fuel
                        if(TUBEEND.LE.0.0) then       
                            SpeciesMF(1,J,fuel1)=fuel1frac
                            SpeciesMF(0,J,fuel1)=SpeciesMF(1,J,fuel1)
                            SpeciesMF(1,J,fuel2)=fuel2frac
                            SpeciesMF(0,J,fuel2)=SpeciesMF(1,J,fuel2) 
                            SpeciesMF(1,J,fuel3)=fuel3frac 
                            SpeciesMF(0,J,fuel3)=SpeciesMF(1,J,fuel3)
                            SpeciesMF(1,J,fuel4)=fuel4frac 
                            SpeciesMF(0,J,fuel4)=SpeciesMF(1,J,fuel4)
                            SpeciesMF(1,J,fuel5)=fuel5frac 
                            SpeciesMF(0,J,fuel5)=SpeciesMF(1,J,fuel5)
                            SpeciesMF(1,J,fuel6)=fuel6frac 
                            SpeciesMF(0,J,fuel6)=SpeciesMF(1,J,fuel6)
                            SMoleF(:) = SpeciesMF(0,J,:)
                            CALL CKXTY(SMoleF,IWORK,WORK,SMassF)
                            SpeciesMF(0,J,:) = SMassF(:)
                            SpeciesMF(1,J,:) = SMassF(:)
                        endif                          
            
                    !for radial location larger than the burner tube
                    ELSE
                        U(2,J)=UAIR
                        U(1,J)=UAIR
                        U(0,J)=U(1,J)
                        V(2,J)=0.d0
                        V(1,J)=0.d0
                        V(0,J)=0.d0
                        T(1,J)=TAIR
                        T(0,J)=T(1,J)
                        SpeciesMF(1,J,ox1)=ox1frac
                        SpeciesMF(0,J,ox1)=SpeciesMF(1,J,ox1)
                        SpeciesMF(1,J,ox2)=ox2frac
                        SpeciesMF(0,J,ox2)=SpeciesMF(1,J,ox2)
                        SpeciesMF(1,J,ox3)=ox3frac
                        SpeciesMF(0,J,ox3)=SpeciesMF(1,J,ox3)
                        SMoleF(:) = SpeciesMF(0,J,:)
                        CALL CKXTY(SMoleF,IWORK,WORK,SMassF)
                        SpeciesMF(0,J,:) = SMassF(:)
                        SpeciesMF(1,J,:) = SMassF(:)
                    ENDIF

                ENDDO
            endif
            
            DEallocate(SMassF)
            DEallocate(SMoleF)
            DEallocate(Xrestart)
            DEallocate(Yrestart)
            DEallocate(UTrestart)
            DEallocate(V_tot)
            DEallocate(U_tot)
            DEallocate(VTrestart)
            DEallocate(Trestart)
            DEallocate(Prestart)
            DEallocate(PCrestart)
            DEallocate(RHOrestart)
            DEallocate(SpeciesRestart)
            DEallocate(SootSecRestart)                           
                 
        end subroutine SetInitial
        !*************************************************************************************************************
        
        !*************************************************************************************************************              
        SUBROUTINE UpdateDensity
        !this subroutine calculates the new density field, and updates species mass fractions to account for soot  
        
            !local variable declaration
            double precision GLOB_soot, SM, RHOCK, RHO_GLOB, the_m_of_mono, the_m_of_soot, sootMF
            double precision, dimension(:), allocatable :: SMassF
            !end local variable declarations 
            
            !allocate arrays
            allocate(SMassF(KK))
            !done allocating arrays 
            
            ! Prevent physically impossible values for Soot
            do J=1,M1
                do I=0,L0      
                     
                    !only do this if not in a solid region
                    if(.NOT.SOLIDCVP(I,J)) then                    
                        
                        DO k = 1, MSection*2
                            IF (SootSec(i,j,k).LT.0.d0) SootSec(i,j,k)=0.d0
                        ENDDO                                                
                        
                        if(SootModel.Eq.1) then
                            !Make sure that number of primary particles is not less than the number of clusters
                            !and if clusters = 0, aggregates = 0
                            DO k=1,MSection
                                IF (SootSec(i,j,k+MSection).LT.SootSec(i,j,k)) SootSec(i,j,k+MSection)=SootSec(i,j,k)
                            enddo
                            
                            !if particle aggregate density (Na) is negligible, set Np and Na equal to zero
                           ! DO k=1,MSection
!                                GLOB_soot=SootSec(i,j,k)*rho(i,j)
!                                IF (GLOB_soot.LE.1.0d0) then
!                                    SootSec(i,j,k)=0.
!                                    SootSec(i,j,k+MSection)=SootSec(i,j,k)
!                                endif
!                            enddo
                            !if mass of primary particles is greater than the mass of aggregates in a section, set 
                            !# of primary particles in that section so the mass of primary particles equals the 
                            !mass of aggregates           
                            DO k=1,MSection 
                                if (SootSec(i,j,k+MSection).GT.SootSec(i,j,k)*DEXP(XS(k))/DEXP(XS(1))) then
                                    SootSec(i,j,k+MSection)=SootSec(i,j,k)*DEXP(XS(k))/DEXP(XS(1))
                                endif
                            enddo                               
                        else     
                            !if primary particle number density is 1 or less, set Ys and N to zero
                            GLOB_soot=SootSec(i,j,2)   !N  #/g-gas
                            IF (GLOB_soot.LE.1.0) then
                                SootSec(i,j,1)=0.d0
                                SootSec(i,j,2)=0.d0
                            endif

                            ! To avoid primary particle size being smaller than monomer
                            the_m_of_mono=c_min_mono*C_MW/Av    ! in g/#
                            the_m_of_soot=SootSec(i,j,1)/(SootSec(i,j,2)+smallnum) ! in g/#
                            if (the_m_of_soot .lt.the_m_of_mono) then
                                SootSec(i,j,2)=SootSec(i,j,1)/the_m_of_mono !change the number of particles so that each particle has the mass of a monomer
                            endif

                        endif
                    
                    endif
                
                enddo
            enddo
            
            !update the density field, and recalculate spieces mass fractions to account for soot
            do J=1,M1
                do I=0,L0
                                
                    !If termpeture is above MaxTemp, set to MaxTemp and write out a message
                    IF(T(I,J).GT.MaxTemp) then
                        T(I,J)=MaxTemp 
                    ENDIF      
                    
                    IF(T(I,J).LT.300.d0) T(I,J)=300.d0
                    
                    !If species mass fraction less than 0, set to 0
                    DO K=1,KK
                        IF(SpeciesMF(I,J,K).LT.0.0d0) SpeciesMF(I,J,K)=0.0d0
                    ENDDO
                    
                    !only do this if not in a solid region
                    if(.NOT.SOLIDCVP(I,J)) then
                        
                        !get soot mass fraction
                        if(SootModel.eq.1) then
                            call SootMassFrac(SootSec(i,j,1:MSection),XS,sootMF)
!                            if (sootMF.NE.sootMF) print*, myid, I,J,'soot Nan'
                            IF(sootMF.LT.0.0d0) sootMF = 0.0d0
                        else
                            sootMF = SootSec(I,J,1)
                        endif       
                        
                        !store current, unadjusted, mass fractions in SMassF
                        DO K=1,KK
                            SMassF(K) = SpeciesMF(I,J,K)
                        ENDDO
                        
                        !adjust mass fractions to account for soot
                        
                        !calculate summation of species and soot mass fractions 
                        SM=0.0d0
                        DO K=1,KK
!                             ! RDB Species NaN
!                            if(SpeciesMF(I,J,K).NE.SpeciesMF(I,J,K).and.myid.eq.0) then
!                                print*, myid, I,J,K,'species Nan'
!                            endif
                            SM=SM+SpeciesMF(I,J,K)
                        ENDDO
                        SM=SM+sootMF
                        
                        !alter all species mass fractions to account for soot
                        DO K=1,KK
                            SpeciesMF(I,J,K)=SpeciesMF(I,J,K)/(SM+smallnum)
                        ENDDO
                        
                        ! get mixture density accouting for soot
                        
                        !Use chemkin to calculate gaseous density (no soot), store in RHOCK
                        CALL CKRHOY(PRESSURE*PCKIN,T(I,J),SMassF,IWORK,WORK,RHOCK) !uses unmodified species mass fractions
!                        ! RDB test RHOCK NaN
!                        if (RHOCK.ne.RHOCK) then
!                        print*,myid, I,J,T(I,J),'RHO'
!                       ! do k=1,KK
!                       !     print*,k,SMassF(k)
!                       ! enddo
!                       ! pause
!                        endif
                        
                        if(iter.EQ.0.0) then !if starting from an arbitray guess, no need to account for soot       ! linea a revisar
                            RHO(I,J) = RHOCK
                        else                             
                            !figure out density including soot
                            RHO_GLOB = RHOCK+(sootMF*RHO(I,J))*(1.0d0-RHOCK/densityP)
                            if(ITER.NE.FIRST) THEN !if not the first iteration
                                !update density, using relaxation
                                RHO(I,J)=RHO_GLOB*RELAX(5)+RHO(I,J)*(1.0d0-RELAX(5))
                            ENDIF     
                        endif      
                    else !if in a solid region, set density to the solid density
                        RHO(I,J)= SolidDensity
                        SpeciesMF(I,J,:) = 0.d0
                        SootSec(I,J,:) = 0.d0 
                    endif

                enddo
            enddo  
            
            DEallocate(SMassF)                      

        
        end subroutine UpdateDensity                       
        !*******************************************************************************************************

        !******************************************************************************************************
        SUBROUTINE SetBoundary
        !This subroutine sets the boundary conditions at the top of the domain, and at the inner and outer radial locations
            
            !local variable declarations
            double precision FLOWIN, FLOWOUT, FACTOR, FL 
            !done local variable declaration
            
            
            !Calculate the total mass flow in the boundary
            if(MYID .eq. 0) then
                FLOWIN=0.d0
                DO J=2,M2
                    FLOWIN=FLOWIN+RHO(1,J)*U(2,J)*YCVR(J)
                ENDDO
            endif
            CALL MPI_BCAST(flowin,1, MPI_double_precision,0,MPI_COMM_WORLD,IERR)

            !imposing inner/outer radial boundary conditons on soot and species
            SpeciesMF(0:L0,1,:)=SpeciesMF(0:L0,2,:)
            SpeciesMF(0:L0,M1,:)=SpeciesMF(0:L0,M2,:)
            SootSec(0:L0,1,:)=SootSec(0:L0,2,:)
            SootSec(0:L0,M1,:)=SootSec(0:L0,M2,:)            
            
            !imposing inner/outer radial boundary conditons on U-velocity and GLOBerature
            if(OutBound.EQ.1) then !if outer radial boundary is a wall, set U velocity equal to zero
                U(0:L0,M1)=0.0d0
            else                   !if outer radial boundary is not a wall, use free-slip condition
                U(0:L0,M1)=U(0:L0,M2)
            endif  
            V(0:L0,M1)=0.0d0
            V(0:L0,1)=0.0d0
            U(0:L0,1)=U(0:L0,2)
            T(0:L0,1)=T(0:L0,2)
            T(0:L0,M1)=T(0:L0,M2)
            
            !Calculates mass flow leaving the boundary, and imposes the no gradient
            !boundary conditions at the top of the domain
            if(MYID .eq.(NUMP-1)) then
                FLOWOUT = 0.0d0
                !imposing conditon on species/soot
                SpeciesMF(L1,2:M2,:)=SpeciesMF(L2,2:M2,:)
                SpeciesMF(L0,2:M2,:)=SpeciesMF(L1,2:M2,:)
                SootSec(L1,2:M2,:)=SootSec(L2,2:M2,:)
                SootSec(L0,2:M2,:)=SootSec(L1,2:M2,:)

                !imposing condition on Temperature, and v-velocity
                T(L1,2:M2)=T(L2,2:M2)
                T(L0,2:M2)=T(L1,2:M2)
                V(L1,2:M2)=V(L2,2:M2)                
                
                DO J=2,M2
                    FLOWOUT=FLOWOUT+RHO(L1,J)*U(L2,J)*YCVR(J)
                ENDDO 
            endif

            !Adjust outlet velocity to acheive overall mass conservation
            if(myid.eq.(NUMP-1)) then
                FACTOR=FLOWIN/FLOWOUT   !Determine mass imbalance
                !adjust outlet flow accordingly
                U(L1,2:M2)=U(L2,2:M2)*FACTOR
                U(L0,2:M2)=U(L1,2:M2)
            endif
        
        end subroutine SetBoundary
        !*******************************************************************************************************

        !*******************************************************************************************************
        SUBROUTINE SolveU
        !This subroutine solves the u-momentum equation in a globally coupled manner. The coefficient matrix is formed
        !and solved using the TDMA for penta-diagonal matrices iterative algorithm.

            !local variable declarations
            integer :: IST, JST
            double precision :: FLOW, DIFF, ACOF
            double precision :: FL, FLP, FLM, GM, GMM, VOL, APT
            DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: GAM, AP, AIP, AIM, AJP, AJM, CON
            double precision, dimension (:), allocatable :: SMoleF
            double precision :: DUDZP, DUDZM, CON1, DRVDRP, DRVDRM, GAMNW, GAMNE, GAMN, GAMSW, GAMSE, GAMSS
            double precision :: DVDZN, DVDZS, CON3, CON2
            !end local variable declarations
            
            !allocate matrices
            allocate (GAM(LM1,M1))
            allocate (CON(LM1,M1))
            allocate (AP(LM1,M1))
            allocate (AIM(LM1,M1))
            allocate (AIP(LM1,M1))
            allocate (AJP(LM1,M1))
            allocate (AJM(LM1,M1))
            allocate(SMoleF(KK))
            !allocate matrices

            !initial matrix coefficients and GAM to zero
            AP = 0.0d0
            AIM = 0.0d0
            AIP = 0.0d0
            AJM = 0.0d0
            AJP = 0.0d0
            CON = 0.0d0
            GAM = 0.0d0

            !set up the starting index for DO loops that calculate matrix coefficients, and TDMA
            if(MYID .eq. 0) then
                IST=3
            else
                IST=2
            endif

            JST=2        
            
            !get the values of viscosity (GAM) and partially calculate constant terms (CON) in the u-mom equation
            do  J=1,M1
                do  I=1,L1                                    
                    
                    if(.NOT.SOLIDCVP(I,J)) then
                        if(J.NE.M1) then
                            CALL CKYTX(SpeciesMF(I,J,:),IWORK,WORK,SMoleF)
                            CALL MCAVIS(T(I,J),SMoleF,TPWRK,GAM(I,J))
                        else
                            GAM(I,J) = 0.d0
                        endif 
                    else
                        GAM(I,J) =  SOLIDVIS      
                    endif
                    
                enddo
            enddo
    
            if (MYID .eq. NUMP-1) then
                do J=1,M1
                    GAM(L1,J)=0.0d0
                enddo
            endif
    
            DO J=JST,M2
                DO I=IST,L2
                    CON(I,J) = RHO(I,J)*GX
                ENDDO
            ENDDO
            !done getting values of GAM and partially calculating CON

            !loop to calculate the values of the matrix coefficients for all CVs on a given CPU
            do J=JST,M2
                
                do I = IST, L2
                    
                    !calculate Asouth for the CV above and Anorth for the current CV
                    
                    !Special case when on the first CV on a CPU, must calculate Asouth for the current CV
                    if(I.EQ.IST) then !Asouth for the first CV on a CPU
                        if(myid.EQ.0) then   !For CPU-0, inlet velocity is prescribed at the south face
                            FLOW=ARX(J)*U(2,J)*RHO(1,J)
                            DIFF=ARX(J)*GAM(1,J)/XCV(2)
                        else  !other CPUs must average the flow in the current CV, with that below it, to get flow at the south face
                            FL=U(I-1,J)*(FX(I-1)*RHO(I-1,J)+FXM(I-1)*RHO(I-2,J))
                            FLP=U(I,J)*(FX(I)*RHO(I,J)+FXM(I)*RHO(I-1,J))
                            FLOW=ARX(J)*0.5d0*(FL+FLP)
                            DIFF=ARX(J)*GAM(I-1,J)/XCV(I-1)
                        endif                
                        CALL DIFLOW(FLOW,DIFF,ACOF)                            
                        
                        AIM(I,J)=ACOF+DMAX1(0.0D0,FLOW) !only need to set Asouth
                        
                    endif
                    
                    !calculate Asouth for the CV above and Anorth for the current CV for all other CVs
                    
                    IF(MYID .eq. (NUMP-1) .and. I .EQ. L2) GO TO 104 !special case when on the last CV of the last CPU
                    
                    !calcalate the flow at the upper face of the current CPU by averaging the flow at the current CV with that above it
                    FL=U(I,J)*(FX(I)*RHO(I,J)+FXM(I)*RHO(I-1,J))
                    FLP=U(I+1,J)*(FX(I+1)*RHO(I+1,J)+FXM(I+1)*RHO(I,J))
                    FLOW=ARX(J)*0.5d0*(FL+FLP)
                    DIFF=ARX(J)*GAM(I,J)/XCV(I)
                    GO TO 105
                104 FLOW=ARX(J)*U(L1,J)*RHO(L1,J) !calculate flow of the upper face for the last CV on the last CPU, which is prescribed, so no averaging needed
                    DIFF=ARX(J)*GAM(L1,J)/XCV(L2)
                    
                105 CALL DIFLOW(FLOW,DIFF,ACOF) 
                   
                    AIM(I+1,J)=ACOF+DMAX1(0.0D0,FLOW) !determine the Asouth for the CV above the current CV
                    AIP(I,J)=AIM(I+1,J)-FLOW   !determine the Anorth for the current CV
                    
                    !Asouth for the upper CV and Anorth for the current CV are calculated

                    !calculate Awest for the CV to the right, and Aeast for the current CV
                    IF (J .EQ. M2) GOTO 106 !special case if at the last CV in the domain
                    !due to the staggered grid, the flow through the east face of the CV has to be calculated in 2 parts, FL, FLM
                    FL=XCVI(I)*V(I,J+1)*(FY(J+1)*RHO(I,J+1)+FYM(J+1)*RHO(I,J))
                    FLM=XCVIP(I-1)*V(I-1,J+1)*(FY(J+1)*RHO(I-1,J+1)+FYM(J+1)* RHO(I-1,J))
                    GM=GAM(I,J)*GAM(I,J+1)/(YCV(J)*GAM(I,J+1)+YCV(J+1)*GAM(I,J)+ smallnum)*XCVI(I)
                    GMM=GAM(I-1,J)*GAM(I-1,J+1)/(YCV(J)*GAM(I-1,J+1)+YCV(J+1)* GAM(I-1,J)+smallnum)*XCVIP(I-1)
                    DIFF=RMN(J+1)*2.d0*(GM+GMM)
                    GO TO 107
                106 FL=XCVI(I)*V(I,M1)*RHO(I,M1)
                    FLM=XCVIP(I-1)*V(I-1,M1)*RHO(I-1,M1)
                    DIFF=R(M1)*(XCVI(I)*GAM(I,M1)+XCVIP(I-1)*GAM(I-1,M1))/YDIF(M1)
                107 FLOW=RMN(J+1)*(FL+FLM)
                    CALL DIFLOW(FLOW,DIFF,ACOF)
                    AJM(I,J+1)=ACOF+DMAX1(0.0D0,FLOW)
                    AJP(I,J)=AJM(I,J+1)-FLOW
                    !done calculating Awest and Aeast

                    !calculate Ap and the constant term, and coefficients (DU) for pressure correction
                    VOL=YCVR(J)*XCVS(I)
                    APT=(RHO(I,J)*XCVI(I)+RHO(I-1,J)*XCVIP(I-1)) /(XCVS(I)*DT)
                    AP(I,J)=AP(I,J)-APT
                    CON(I,J)=CON(I,J)+APT*U(I,J)
                    AP(I,J)=(-AP(I,J)*VOL+AIP(I,J)+AIM(I,J)+AJP(I,J)+AJM(I,J)) /RELAX(1)
                    CON(I,J)=CON(I,J)*VOL+(1.d0-RELAX(1))*AP(I,J)*U(I,J)
                    DU(I,J)=VOL/XDIF(I)
                    CON(I,J)=CON(I,J)+DU(I,J)*(P(I-1,J)-P(I,J))
                    DU(I,J)=DU(I,J)/AP(I,J)    
                              
                enddo
            enddo

            call MPI_BARRIER(MPI_COMM_WORLD,IERR)
        
            !Solve the coefficient matrix with TDMA
            Call TDMASolver(AP,AIP,AIM,AJP,AJM,CON,U,IST,JST,NTIMES(1))

            DEALLOCATE (GAM)
            DEALLOCATE (CON)
            DEALLOCATE (AP)
            DEALLOCATE (AIM)
            DEALLOCATE (AIP)
            DEALLOCATE (AJP)
            DEALLOCATE (AJM)
            DEALLOCATE(SMoleF)                                                  

        end subroutine SolveU
        !*************************************************************************************************************

        !**************************************************************************************************************
        SUBROUTINE SolveV
        !This subroutine solves the v-momentum equation in a globally coupled manner. The coefficient matrix is formed
        !and solved using the TDMA for penta-diagonal matrices iterative algorithm.

            !local variable declarations
            integer IST, JST
            double precision FLOW, DIFF, ACOF
            double precision FL, FLM, FLP, GM, GMM, AREA, APT, VOL
            DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: GAM, AP, AIP, AIM, AJP, AJM, CON
            double precision, dimension (:), allocatable :: SMoleF
            double precision :: DVDRP, DVDRM, CON1, DVDRVP, DVDRVM, CON2, DUDZP, DUDZM, CON3, GAMNE, GAMSE, GAME
            double precision :: GAMSW, GAMNW, GAMW, DUDRE, DUDRW, CON4, GAMV, DRVDR, UE, UW, DUDZ, CON7, CON6       
            !end local variable declarations
            
            !allocate matrices
            allocate (GAM(LM1,M1))
            allocate (CON(LM1,M1))
            allocate (AP(LM1,M1))
            allocate (AIM(LM1,M1))
            allocate (AIP(LM1,M1))
            allocate (AJP(LM1,M1))
            allocate (AJM(LM1,M1))
            allocate(SMoleF(KK))
            !allocate matrices

            !initial matrix coefficients and GAM to zero
            AP = 0.0d0
            AIM = 0.0d0
            AIP = 0.0d0
            AJM = 0.0d0
            AJP = 0.0d0
            CON = 0.0d0
            GAM = 0.0d0

            !set up the starting index for DO loops
            IST=2
            JST=3   
            
            !get the values of viscosity (GAM) and partially calculate AP, and constant terms (AP, CON) in the V-mom equation 
            do  J=1,M1
                do  I=1,L1                                    
                    
                    if(.NOT.SOLIDCVP(I,J)) then
                        if(J.NE.M1) then
                            CALL CKYTX(SpeciesMF(I,J,:),IWORK,WORK,SMoleF)
                            CALL MCAVIS(T(I,J),SMoleF,TPWRK,GAM(I,J))
                        else
                            GAM(I,J) = 0.d0
                        endif 
                    else
                        GAM(I,J) =  SOLIDVIS      
                    endif
                    
                enddo
            enddo
            
            DO I=IST,L2
                DO J=JST,M2
                    CON(I,J)= RHO(I,J)*GY
                    AP(I,J)=-1.d0*GAMV/(RMN(J)*RMN(J))
                ENDDO
            ENDDO
            !done gettings the values of viscosity (GAM) and partially calculating AP, and constant terms (AP, CON) in the V-mom equation
            
            !loop to calculate the values of the matrix coefficients for all CVs on a given CPU
            DO J=3,M2
                
                !calculate Asouth for the first CV on each CPU
                if(myid .eq.0) then !For CPU-0, inlet velocity is prescribed at the south face; however, due to staggered grid, the flow is a
                    FL=ARXJ(J)*U(2,J)*RHO(1,J)    !combination of FL and FLM
                    FLM=ARXJP(J-1)*U(2,J-1)*RHO(1,J-1)
                    FLOW=FL+FLM
                    DIFF=(ARXJ(J)*GAM(1,J)+ARXJP(J-1)*GAM(1,J-1))/XDIF(2)
                    CALL DIFLOW(FLOW,DIFF,ACOF)
                    AIM(2,J)=ACOF+DMAX1(0.0D0,FLOW)
                else
                    FL=ARXJ(J)*U(2,J)*(FX(2)*RHO(2,J)+FXM(2)*RHO(1,J))
                    FLM=ARXJP(J-1)*U(2,J-1)*(FX(2)*RHO(2,J-1)+FXM(2)* RHO(1,J-1))
                    GM=GAM(1,J)*GAM(2,J)/(XCV(1)*GAM(2,J)+XCV(2)*GAM(1,J)+ smallnum)*ARXJ(J)
                    GMM=GAM(1,J-1)*GAM(2,J-1)/(XCV(1)*GAM(2,J-1)+XCV(2)* GAM(1,J-1)+ smallnum)*ARXJP(J-1)
                    DIFF=2.d0*(GM+GMM)
                    FLOW=FL+FLM
                    CALL DIFLOW(FLOW,DIFF,ACOF)
                    AIM(2,J)=ACOF+DMAX1(0.0D0,FLOW)
                endif


                DO I=2,L2
                    !calculate Anorth for the current CV, and Asouth for the CV below
                    IF(myid .eq.(NUMP-1) .and. I .EQ. L2) then !special case if at the last CV in the domain, as flow is prescibed at that interface
                        FL=ARXJ(J)*U(L1,J)*RHO(L1,J)
                        FLM=ARXJP(J-1)*U(L1,J-1)*RHO(L1,J-1)
                        DIFF=(ARXJ(J)*GAM(L1,J)+ARXJP(J-1)*GAM(L1,J-1))/XDIF(L1)
                    else
                        FL=ARXJ(J)*U(I+1,J)*(FX(I+1)*RHO(I+1,J)+FXM(I+1)*RHO(I,J))
                        FLM=ARXJP(J-1)*U(I+1,J-1)*(FX(I+1)*RHO(I+1,J-1)+FXM(I+1)* RHO(I,J-1))
                        GM=GAM(I,J)*GAM(I+1,J)/(XCV(I)*GAM(I+1,J)+XCV(I+1)*GAM(I,J)+ smallnum)*ARXJ(J)
                        GMM=GAM(I,J-1)*GAM(I+1,J-1)/(XCV(I)*GAM(I+1,J-1)+XCV(I+1)*GAM(I,J-1)+ smallnum)*ARXJP(J-1)
                        DIFF=2.d0*(GM+GMM)
                    endif
                    FLOW=FL+FLM
                    CALL DIFLOW(FLOW,DIFF,ACOF)
                    AIM(I+1,J)=ACOF+DMAX1(0.0D0,FLOW)
                    AIP(I,J)=AIM(I+1,J)-FLOW
                    !Asouth and Anorth calculated
                    
                    !calculate Aeast of the current CV, and Awest of the CV to the right
                    IF(J .EQ. M2) GO TO 206 !special case of when at the last CV in the domain
                    AREA=R(J)*XCV(I)
                    FL=V(I,J)*(FY(J)*RHO(I,J)+FYM(J)*RHO(I,J-1))*RMN(J)
                    FLP=V(I,J+1)*(FY(J+1)*RHO(I,J+1)+FYM(J+1)*RHO(I,J))*RMN(J+1)
                    FLOW=(FV(J)*FL+FVP(J)*FLP)*XCV(I)
                    DIFF=AREA*GAM(I,J)/YCV(J)
                    GO TO 207
                206 AREA=R(M1)*XCV(I)
                    FLOW=AREA*V(I,M1)*RHO(I,M1)
                    DIFF=AREA*GAM(I,M1)/YCV(M2)
                207 CALL DIFLOW(FLOW,DIFF,ACOF)
                    AJM(I,J+1)=ACOF+DMAX1(0.0D0,FLOW)
                    AJP(I,J)=AJM(I,J+1)-FLOW
                    !Awest and Aeast calculated
                    
                    !calculate Ap, b, and pressure correction coefficients
                    VOL=YCVRS(J)*XCV(I)
                    APT=(ARXJ(J)*RHO(I,J)+ARXJP(J-1)*RHO(I,J-1))/(YCVRS(J)*DT)   
                    AP(I,J)=AP(I,J)-APT  !AP partially calculated in gammaANDsource
                    CON(I,J)=CON(I,J)+APT*V(I,J)  
                    AP(I,J)=(-AP(I,J)*VOL+AIP(I,J)+AIM(I,J)+AJP(I,J)+AJM(I,J)) /RELAX(2)
                    CON(I,J)=CON(I,J)*VOL+(1.d0-RELAX(2))*AP(I,J)*V(I,J)
                    DV(I,J)=VOL/YDIF(J)
                    CON(I,J)=CON(I,J)+DV(I,J)*(P(I,J-1)-P(I,J))
                    DV(I,J)=DV(I,J)/AP(I,J)  
                
                enddo
            enddo
            
            call MPI_BARRIER(MPI_COMM_WORLD,IERR)
            
            !Solve the coefficient matrix with TDMA
            Call TDMASolver(AP,AIP,AIM,AJP,AJM,CON,V,IST,JST,NTIMES(2))

            DEALLOCATE (GAM)
            DEALLOCATE (CON)
            DEALLOCATE (AP)
            DEALLOCATE (AIM)
            DEALLOCATE (AIP)
            DEALLOCATE (AJP)
            DEALLOCATE (AJM)
            DEALLOCATE(SMoleF)                                                              

        end subroutine SolveV
        !*************************************************************************************************************

        !**************************************************************************************************************
        SUBROUTINE SolvePC
        !This subroutine solves the pressure correction equation in a globally coupled manner. The coefficient matrix is formed
        !and solved using the TDMA for penta-diagonal matrices iterative algorithm.

            !local variable declarations
            integer IST, JST
            double precision FLOW, ARHO
            DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: AP, AIP, AIM, AJP, AJM, CON
            DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: DUtemp, DVtemp
            !end local variable declarations
            
            !allocate matrices
            allocate (CON(LM1,M1))
            allocate (AP(LM1,M1))
            allocate (AIM(LM1,M1))
            allocate (AIP(LM1,M1))
            allocate (AJP(LM1,M1))
            allocate (AJM(LM1,M1))
            allocate (DUtemp(M1))
            allocate (DVtemp(M1))
            !allocate matrices
            
            !initialize matrices
            CON=0.d0
            AP =0.d0
            AIM=0.d0
            AIP=0.d0
            AJP=0.d0
            AJM=0.d0
            
            !set up the starting index for DO loops
            IST=2
            JST=2                   
            
            !pass the values for DU and DV to 1 lower cpu rank
            !set up temporary arrays for passing DU and DV
            DUtemp(:) = DU(2,:)
            DVtemp(:) = DV(2,:)
            
            if(myid .ne. 0) then
                CALL MPI_SEND(DUtemp(1), 1*M1, MPI_DOUBLE_PRECISION,myid-1,11700+myid,MPI_COMM_WORLD,IERR)
                CALL MPI_SEND(DVtemp(1), 1*M1, MPI_DOUBLE_PRECISION,myid-1,11800+myid,MPI_COMM_WORLD,IERR)
            endif

            if(myid .ne. (NUMP-1)) then
                CALL MPI_RECV(DUtemp(1), 1*M1, MPI_DOUBLE_PRECISION, myid+1,11700+myid+1,MPI_COMM_WORLD,STATUS,IERR)
                CALL MPI_RECV(DVtemp(1), 1*M1, MPI_DOUBLE_PRECISION, myid+1,11800+myid+1,MPI_COMM_WORLD,STATUS,IERR)
            endif
            
            !store passed values in DU, and DV
            DU(L1,:) = DUtemp(:)
            DV(L1,:) = DVtemp(:)
            !done passing DU and DV to 1 lower cpu rank

            call MPI_BARRIER(MPI_COMM_WORLD,IERR)

            !pass DU and DV to 1 higher cpu rank
            !set up temporary arrays for passing DU and DV
            DUtemp(:) = DU(L2,:)
            DVtemp(:) = DV(L2,:)
            if(myid .ne. (NUMP-1)) then
                CALL MPI_SEND(DUtemp(1), 1*M1, MPI_DOUBLE_PRECISION,myid+1,11400+myid,MPI_COMM_WORLD,IERR)
                CALL MPI_SEND(DVtemp(1), 1*M1, MPI_DOUBLE_PRECISION,myid+1,11500+myid,MPI_COMM_WORLD,IERR)
            endif

            !receive DV and DU
            if(myid .ne. 0) then
                CALL MPI_RECV(DUtemp(1), 1*M1, MPI_DOUBLE_PRECISION, myid-1,11400+myid-1,MPI_COMM_WORLD,STATUS,IERR)
                CALL MPI_RECV(DVtemp(1), 1*M1, MPI_DOUBLE_PRECISION, myid-1,11500+myid-1,MPI_COMM_WORLD,STATUS,IERR)
            endif

            !store passed values in DU, and DV
            DU(1,:) = DUtemp(:)
            DV(1,:) = DVtemp(:)
            !done passing DU and DV to 1 higher cpu rank
            !done passing DV and DU

            call MPI_BARRIER(MPI_COMM_WORLD,IERR)


            !loop to calculate the values of the matrix coefficients for all CVs on a given CPU
            DO J=2,M2
                !set Asouth, and source term for bottom CV
                if(myid .eq. 0) then !if on CPU-0
                    ARHO=ARX(J)*RHO(1,J)
                    CON(2,J)=CON(2,J)+ARHO*U(2,J)
                    AIM(2,J)=0.d0
                !set Asouth and source term for other processes at bottom boundary
                else
                    ARHO=ARX(J)*(FX(2)*RHO(2,J)+FXM(2)*RHO(1,J))
                    FLOW=ARHO*U(2,J)
                    CON(2,J)=CON(2,J)+FLOW
                    AIM(2,J)=ARHO*DU(2,J)   
                endif               

                DO I=2,L2
                !calculated Anorth, Asouth and source term for interior nodes               
                    IF(myid .eq. (NUMP-1) .and. I .EQ. L2) then !special case if at last CV in the domain
                        ARHO=ARX(J)*RHO(L1,J)
                        CON(I,J)=CON(I,J)-ARHO*U(L1,J)
                        AIP(I,J)=0.0d0
                    else
                        ARHO=ARX(J)*(FX(I+1)*RHO(I+1,J)+FXM(I+1)*RHO(I,J))
                        FLOW=ARHO*U(I+1,J)
                        CON(I,J)=CON(I,J)-FLOW
                        CON(I+1,J)=CON(I+1,J)+FLOW
                        AIP(I,J)=ARHO*DU(I+1,J)
                        AIM(I+1,J)=AIP(I,J)  
                    endif
                    
                    !calculate Awest, Aeast, and source term for interior nodes
                    IF(J .EQ. M2) then !special case if at last CV in the domain
                        ARHO=RMN(M1)*XCV(I)*RHO(I,M1)
                        CON(I,J)=CON(I,J)-ARHO*V(I,M1)
                        AJP(I,J)=0.0d0
                    else
                        ARHO=RMN(J+1)*XCV(I)*(FY(J+1)*RHO(I,J+1)+FYM(J+1)*RHO(I,J))
                        FLOW=ARHO*V(I,J+1)
                        CON(I,J)=CON(I,J)-FLOW
                        CON(I,J+1)=CON(I,J+1)+FLOW
                        AJP(I,J)=ARHO*DV(I,J+1)
                        AJM(I,J+1)=AJP(I,J)
                    endif  
                    
                    !determine Ap
                    AP(I,J)=AIP(I,J)+AIM(I,J)+AJP(I,J)+AJM(I,J)
                    !initialize pressure correction to zero
                    PC(I,J)=0.0d0
                
                ENDDO
            ENDDO 
            
            call MPI_BARRIER(MPI_COMM_WORLD,IERR)
            
            !Solve the coefficient matrix with TDMA  
            Call TDMASolver(AP,AIP,AIM,AJP,AJM,CON,PC,IST,JST,NTIMES(3))                                                                
                    
            !Correct the pressure and velocity fields
    
            if(myid.eq.0) then
                DO  J=2,M2
                    DO  I=2,L1
                        P(I,J)=P(I,J)+PC(I,J)*RELAX(3)
                        IF(I .NE. 2) U(I,J)=U(I,J)+DU(I,J)*(PC(I-1,J)-PC(I,J))
                        IF(J .NE. 2) V(I,J)=V(I,J)+DV(I,J)*(PC(I,J-1)-PC(I,J))
                    enddo
                enddo

            else if(myid.eq.(NUMP-1)) then
                DO  J=2,M2
                    DO  I=1,L2
                        P(I,J)=P(I,J)+PC(I,J)*RELAX(3)
                        U(I,J)=U(I,J)+DU(I,J)*(PC(I-1,J)-PC(I,J))
                        IF(J .NE. 2) V(I,J)=V(I,J)+DV(I,J)*(PC(I,J-1)-PC(I,J))
                    enddo
                enddo

            else
                DO  J=2,M2
                    DO  I=1,L1
                        P(I,J)=P(I,J)+PC(I,J)*RELAX(3)
                        U(I,J)=U(I,J)+DU(I,J)*(PC(I-1,J)-PC(I,J))
                        IF(J .NE. 2) V(I,J)=V(I,J)+DV(I,J)*(PC(I,J-1)-PC(I,J))
                    enddo
                enddo
            endif

            DEALLOCATE (CON)
            DEALLOCATE (AP)
            DEALLOCATE (AIM)
            DEALLOCATE (AIP)
            DEALLOCATE (AJP)
            DEALLOCATE (AJM)

        end subroutine SolvePC
        !*************************************************************************************************************

        !*************************************************************************************************************
        SUBROUTINE CVFaceFlow
        !This subroutine calculates the flow across the North, and East, faces of each CV, to be used when calculating the
        !matrix coefficients when solving the species, soot, and GLOBerature equations

            !local variable declarations
            double precision :: AREA
            !end local variable declarations

            DO J = 2, M2
                !calculate the flow for the North face of the CV below the first CV on a CPU
                !or, calculate the flow for the South face of the first CV on a CPU
                IF(.NOT.SOLIDCVP(1,J)) THEN !if not in a solid region
                    if(myid.eq.0) then !special case for CPU-0, flow is prescribed at the face
                        FLOWN(1,J) = ARX(J)*U(2,J)*RHO(1,J)
                    else
                        FLOWN(1,J)=ARX(J)*U(2,J)*(FX(2)*RHO(2,J)+FXM(2)*RHO(1,J))
                    endif
                ELSE
                    FLOWN(1,J) = 0.d0
                ENDIF

                DO I=2,L2
                    !calculate the flow at the north face of all other CVs
                    IF(myid .eq. (NUMP-1) .and. I .EQ. L2) then !special case for last CV of domain, flow is prescribed at the north face
                        FLOWN(I,J)=ARX(J)*U(L1,J)*RHO(L1,J)
                    ELSE
                        if(SOLIDCVP(I,J)) then !if in a solid region
                            FLOWN(I,J)=0.d0
                        else
                            FLOWN(I,J)=ARX(J)*U(I+1,J)*(FX(I+1)*RHO(I+1,J)+FXM(I+1)*RHO(I,J))
                        endif
                    ENDIF

                    !calculate the flow at the east face for all CVs
                    AREA=RMN(J+1)*XCV(I)

                    IF(J .EQ. M2) then !special case is at the end of the domain, the flow is prescribed at the east face
                        FLOWE(I,J)=AREA*V(I,M1)*RHO(I,M1)
                    ELSE
                        if(SOLIDCVP(I,J).OR.SOLIDCVP(I,J+1)) then !if the CV is in the solid region, or the CV to the right is in the solid region
                            FLOWE(I,J)=0.d0
                        else
                            FLOWE(I,J)=AREA*V(I,J+1)*(FY(J+1)*RHO(I,J+1)+FYM(J+1)*RHO(I,J))
                        endif
                    ENDIF

                ENDDO
            ENDDO

        end subroutine CVFaceFlow
        !******************************************************************************************************************

        !*******************************************************************************************************************
        SUBROUTINE SolveSpecies
        !This subroutine solves the species equations in a locally coupled manner (species are coupled to other species in a given
        !CV, but not in a global sense). The resulting linear system at each CV is solved using Gaussian Elimination.

            !local variable declarations
            double precision DIFF, ACOF
            double precision AREA, VOL, APT
            DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE :: GAMS, APS, AIPS, AIMS, AJPS, AJMS, CONS
            DOUBLE PRECISION,DIMENSION(:,:,:,:),ALLOCATABLE :: DSDY
            DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: WDD, ScrubRate, ScrubRateOld, SMassF, WDDP, SMoleF, DIFFU
            DOUBLE PRECISION :: CONXP, CONXM, CONYP, CONYM, CONX_R_D, CONY_R_D, SOURCE, DENOM
            DOUBLE PRECISION :: OMEGA(IREAC) ! RDB tasas de formación
            INTEGER :: ll                   ! RDB tasas de formación
            INTEGER :: totoFlag,totoFlagP(NUMP)   !RDB totoratelimit
            DOUBLE PRECISION :: toto              !RDB totoratelimit
!            DOUBLE PRECISION, PARAMETER :: totoratelimit = 0.02d0     !RDB totoratelimit
            !end local variable declarations
            
            !allocate matrices
            allocate (GAMS(LM1,M1,KK))
            allocate (CONS(LM1,M1,KK))
            allocate (APS(LM1,M1,KK))
            allocate (AIMS(LM1,M1,KK))
            allocate (AIPS(LM1,M1,KK))
            allocate (AJPS(LM1,M1,KK))
            allocate (AJMS(LM1,M1,KK))
            allocate (DSDY(KK,KK,LM1,M1))
            allocate(WDD(KK))
            allocate(ScrubRate(KK))
            allocate(ScrubRateOld(KK))
            allocate(SMassF(KK))
            allocate(SMoleF(KK))
            allocate(DIFFU(KK))
            allocate(WDDP(KK))
            !allocate matrices

            if(reacdetail .eq. 1) RateLimit(:,:)=0.d0         ! RDB tasas de formación

            !Set the old species mass fractions to the current values
            SpeciesMFO = SpeciesMF

            !initialize GAMS, CONS, APS, AIMS, AIPS, AJPS, AJMS, DSDY
            GAMS = 0.0d0
            APS = 0.0d0
            CONS = 0.0d0
            AIMS = 0.0d0
            AIPS = 0.0d0
            AJPS = 0.0d0
            AJMS = 0.0d0
            DSDY = 0.0d0

            !get the values for GAMS, DSDY, and partially calculate CONS and APS
            !calculate GAMS at each CV
            DO J=1, M1
                DO I=1, L1

                  if (T(I,J).GT.3350.0) T(I,J) = 2950.0
                
                    if(.NOT.SOLIDCVP(I,J)) then !if not in a solid region        
                        CALL CKYTX(SpeciesMF(I,J,:),IWORK,WORK,SMoleF)  
                        CALL MCADIF(PRESSURE*PCKIN,T(I,J),SMoleF,TPWRK,DIFFU)
                        DO K=1, KK
                            !IF(DENOM.EQ.0) DENOM = smallnum            ! RDB
                            DENOM = 1.d0 - SpeciesMF(I,J,K) + smallnum  ! RDB
                            GAMS(I,J,K)= RHO(I,J)*(1.d0-SMoleF(K))/(DENOM)*DIFFU(K)
                        ENDDO                                  
                    else !if in a solid region, set GAMS to zero
                        GAMS(I,J,:) = 0.0d0
                    endif
                                   
                ENDDO
            ENDDO 
            
            !set GAMS to zero at boundary and CPU interface CVs
            do i=1,L1
                GAMS(I,1,:)=0.0d0
                GAMS(I,M1,:)=0.0d0
            enddo

            if(myid .eq.0) then
                DO J=1,M1
                    GAMS(1,J,:)=0.0d0
                enddo
            endif

            if(myid .eq.(NUMP-1)) then
                DO J=1,M1
                    GAMS(L1,J,:)=0.0d0
                enddo
            endif
            !done calculating GAMS
            
            totoFlag = 0        ! RDB totoratelimit

            !Calculate derivative of source wrt mass fraction (DSDY), and part of APS and CONS
            DO J=2, M2
                DO I=2,L2
                    
                    if(.NOT.SOLIDCVP(I,J)) then !if not in a solid region   
                        
                        !get the net production rate of each species using chemkin, store in WDD
                        CALL CKWYP(PRESSURE*PCKIN,T(I,J),SpeciesMF(I,J,:),IWORK,WORK,WDD)

                        ! RDB test WDD
                        do K=1,KK
                            if(WDD(K).NE.WDD(K))then
                                print*,"WDD NaN ",K,I,J,myid
                                WDD(K) = 0.0d0      !RDB
                            endif
                            if(WDD(K).GE.1.d50)then    !HUGE(0.d0)/
                                print*,"WDD Inf ",K,I,J,myid
                                WDD(K) = 1.d50      !RDB   HUGE(0.d0)/
                            endif
                            if(WDD(K).LE.-1.d50)then   !HUGE(0.d0)/
                                print*,"WDD -Inf",K,I,J,myid
                                WDD(K) = -1.d50     !RDB   HUGE(0.d0)/
                            endif
                        enddo
                        !
                        !---------------- RDB tasas de formación ---------------
                        ! bloque para imprimir las tasas de formacion de las diferentes reacciones
                      if(reacdetail .eq. 1) then
                        OMEGA(:)=0.d0           ! tasas de formacion
                        ! Rates of progress for the reactions
                        ! cgs units - moles/(cm**3*sec)
                        call CKQYP  (PRESSURE*PCKIN,T(I,J), SpeciesMF(I,J,:), IWORK, WORK, OMEGA)
                        ! la longitud de omega esta determinado por el numero de reacciones
                        ! reacciones seleccionadas para observar su campo
                        ! numero de reacciones definido por IREAC
                        do ll=1,size(iR)
                            Rate(I,J,ll)=OMEGA(iR(ll))
                        enddo
                        !
                        ! se almacenan los límites de las reacciones
                        do ll=1,IREAC
                            RateLimit(1,ll) = min(RateLimit(1,ll),OMEGA(ll))
                            RateLimit(2,ll) = max(RateLimit(2,ll),OMEGA(ll))
                        enddo
                      endif
                        !---------------- RDB tasas de formación ---------------

                        !Calculate the effect of soot on WDD
                        CALL SurfaceRates !get the species feedback rates, on a per section basis
                        
                        CALL NuclChem !get the nucleation rates
                            
                        CALL GasScrubbing(ScrubRate) !calculate the scrubbing rates

                        !update WDD to include the effect of gas phase scrubbing
!                        ScrubRate = 0.0d0   !RDB
                        WDD = WDD + ScrubRate

                        !set ScrubRateOld to ScrubRate
                        ScrubRateOld = ScrubRate
                        
                        !calculate the derivative of source with respect to mass fraction
                        DO K2=1,KK
                            !perturb the K2th species' mass fraction
                            !SpeciesMF(I,J,K2)=SpeciesMFO(I,J,K2)+(ABSO+ABSO*DABS(SpeciesMFO(I,J,K2)))
                            SpeciesMF(I,J,K2)=SpeciesMFO(I,J,K2)+(ABSO**2+ABSO*DABS(SpeciesMFO(I,J,K2)))        ! RDB
                            
!                            !calculate the new net rate of formation for all species with perturb value
                            CALL CKWYP(PRESSURE*PCKIN,T(I,J),SpeciesMF(I,J,:),IWORK,WORK,WDDP)
                            
                            ! RDB test WDD
                            do K=1,KK
                                if(WDDP(K).NE.WDDP(K))then
                                    !print*,"WDDP NaN",K,K2,I,J,myid
                                    WDDP(K) = 0.0d0      !RDB
                                endif
                                if(WDDP(K).GE.1.d50)then   !HUGE(0.d0)/
                                    !print*,"WDDP Inf",K,K2,I,J,myid
                                    WDDP(K) = 1.d50      !RDB  HUGE(0.d0)/
                                endif
                                if(WDD(K).LE.-1.d50)then   !HUGE(0.d0)/
                                    !print*,"WDDP -Inf",K,K2,I,J,myid
                                    WDDP(K) = -1.d50     !RDB  HUGE(0.d0)/
                                endif
                            enddo

                            if(ScrubRateOld(K2).NE.0.d0) then !If a species is involved with soot
                            
                                CALL SurfaceRates !get the species feedback rates, on a per section basis
                                
                                CALL NuclChem !get the nucleation rates
                            
                                CALL GasScrubbing(ScrubRate) !calculate the scrubbing rates
                                
                                !update WDDP to include the effect of gas phase scrubbing                            
                                WDDP = WDDP + ScrubRate
                                
                            else !just use the old scrubbing rates with WDDP
                            
                                WDDP = WDDP + ScrubRateOld
                                
                            endif  
!         
!                            !calculate DSDY                        
                            DO K1=1,KK
                                !DSDY(K1,K2,I,J)=-(WDD(K1)-WDDP(K1))*WT(K1)/(ABSO+ABSO*DABS(SpeciesMFO(I,J,K2)))
            DSDY(K1,K2,I,J)=-(WDD(K1)-WDDP(K1))*WT(K1)/(ABSO**2+ABSO*DABS(SpeciesMFO(I,J,K2)))      ! RDB

!                                ! RDB testDSDY
!                                if(DSDY(K1,K2,I,J).ne.DSDY(K1,K2,I,J))then
!                                print*,'testDSDY',WDD(K1),WDDP(K1),WT(K1),SpeciesMFO(I,J,K2),K1,K2,I,J
!                                DSDY(K1,K2,I,J) = 0.0d0     !RDB
!                                endif
                            ENDDO
                            
                            !set the K2th species' mass fraction back to non perturbed value
                            SpeciesMF(I,J,K2) = SpeciesMFO(I,J,K2)
                        ENDDO
                        !DSDY Calculated
                        
                        !calculate part of the source term (CONS) due to correction and thermal diffusion velocities
                        DO K=1,KK

                            !calculate diffusion across the top face of the CV
                            if(myid .eq.(NUMP-1) .and. I.EQ.L2)THEN
                                CONXP=RHO(I+1,J)*SpeciesMF(I+1,J,K)*(VCX(I+1,J)+VTX(I+1,J,K))
                            ELSEIF(SOLIDCVP(I,J)) then !if CV is in a solid region
                                CONXP = 0.d0
                            ELSE
     CONXP=(RHO(I+1,J)*SpeciesMF(I+1,J,K)*(VCX(I+1,J)+VTX(I+1,J,K))*XCV(I)+RHO(I,J)*SpeciesMF(I,J,K) &
                                *(VCX(I,J)+VTX(I,J,K))*XCV(I+1))/(2.d0*XDIF(I+1))
                            ENDIF
                            
                            !calculate diffusion across the bottom face of the CV
                            if(myid .eq. 0 .and. I.EQ.2)THEN
                                CONXM=RHO(I-1,J)*SpeciesMF(I-1,J,K)*(VCX(I-1,J)+VTX(I-1,J,K))
                            ELSEIF(SOLIDCVP(I-1,J)) THEN !if the CV below is in the solid region, no diffusion from the bottom face of the CV
                                CONXM=0.0d0
                            ELSE
     CONXM=(RHO(I-1,J)*SpeciesMF(I-1,J,K)*(VCX(I-1,J)+VTX(I-1,J,K))*XCV(I)+RHO(I,J)*SpeciesMF(I,J,K)  &
                                *(VCX(I,J)+VTX(I,J,K))*XCV(I-1))/(2.d0*XDIF(I))
                            ENDIF

                            !calculate diffusion across the east of right face of the CV
                            IF(J.EQ.M2)THEN
                                CONYP=RHO(I,J+1)*SpeciesMF(I,J+1,K)*(VCY(I,J+1)+VTY(I,J+1,K))
                            ELSEIF(SOLIDCVP(I,J+1).OR.SOLIDCVP(I,J)) THEN !if the CV to the right or the current CV is in the solid region, no diffusion across the east face
                                CONYP=0.0d0
                            ELSE
                CONYP=(RHO(I,J+1)*SpeciesMF(I,J+1,K)*(VCY(I,J+1)+VTY(I,J+1,K))*YCV(J)+RHO(I,J)  &
                                *SpeciesMF(I,J,K)*(VCY(I,J)+VTY(I,J,K))*YCV(J+1))/(2.d0*YDIF(J+1))
                            ENDIF

                            !calculate diffusion across the west face of the CV
                            IF(J.EQ.2)THEN
                                CONYM=RHO(I,J-1)*SpeciesMF(I,J-1,K)*(VCY(I,J-1)+VTY(I,J-1,K))
                            ELSEIF(SOLIDCVP(I,J-1).OR.SOLIDCVP(I,J)) THEN !if the CV to the left or the current CV  is in the solid region, no diffusion across the west face
                                CONYM=0.0d0
                            ELSE
                CONYM=(RHO(I,J-1)*SpeciesMF(I,J-1,K)*(VCY(I,J-1)+VTY(I,J-1,K))*YCV(J)+RHO(I,J)  &
                                *SpeciesMF(I,J,K)*(VCY(I,J)+VTY(I,J,K))*YCV(J-1))/(2.d0*YDIF(J))
                            ENDIF

                            CONX_R_D=-(CONXP-CONXM)/XCV(I)
                            CONY_R_D=-(RMN(J+1)*CONYP-RMN(J)*CONYM)/(YCV(J)*R(J))

                            !------ reaction rates limited --------- ! RDB
                            toto = WDD(K)*WT(K)
                           ! print*,'toto', toto
                            if(dabs(toto).GT.totoratelimit.and.totoratelimit.ne.0.0d0)then
                                if(toto.gt.0.0d0) then
                                toto = totoratelimit
                                else
                                toto = -totoratelimit
                                endif
                                totoFlag = 1
                            endif
                            !------ reaction rates limited --------- ! RDB
                            !SOURCE=CONX_R_D+CONY_R_D + +WDD(K)*WT(K)
                            SOURCE=CONX_R_D+CONY_R_D + toto         !RDB totoratelimit

                            CONS(I,J,K)=DMAX1(0.0D0,SOURCE)
                            APS(I,J,K)=-DMAX1(0.0D0,-SOURCE)/(SpeciesMF(I,J,K)+smallnum)

                        ENDDO

                    !if in a solid region, set APS, CONS, DSDY to zero
                    else
                         APS(I,J,:) = 0.0d0
                         CONS(I,J,:) = 0.0d0
                         DSDY(:,:,I,J)=0.0d0
                    endif 
                ENDDO
            ENDDO
            !done getting the values for GAMS and DSDY, and partially calculating APS and CONS

            ! --- indicador si se limito la tasa de produccion ---  !RDB totoratelimit
            call MPI_BARRIER(MPI_COMM_WORLD,IERR)
            if(myid .ne. 0) then
                !print*,'send',myid
                call MPI_SEND(totoFlag,1,MPI_INTEGER,0,4439110+myid,MPI_COMM_WORLD,IERR)
            else
                totoFlagP(1)=totoFlag
                do indp = 1, NUMP-1
                call MPI_RECV(totoFlagP(indp+1),1,MPI_INTEGER,indp,4439110+indp,MPI_COMM_WORLD,STATUS,IERR)
                enddo
                if(sum(totoFlagP).GE.1)then
                    print*,'totoFlag'
                endif
            endif
            !call MPI_BARRIER(MPI_COMM_WORLD,IERR)
             ! --- indicador si se limito la tasa de produccion ---

            !calculate the matrix coefficients
            DO K=1,KK
                include 'speciessootcoef.inc'
                
                        CONS(I,J,K)=CONS(I,J,K)+APT*SpeciesMF(I,J,K)
                        APS(I,J,K)=(-APS(I,J,K)*VOL+AIPS(I,J,K)+AIMS(I,J,K)+AJPS(I,J,K)+AJMS(I,J,K))
        
                    ENDDO
                ENDDO
            ENDDO
            
            call MPI_BARRIER(MPI_COMM_WORLD,IERR)
            
             ! RDB test NaN GaussSolver
!             if (myid.eq.0)then
!            do I= 1,LM1
!            do J= 1,M1
!            do K= 1,KK
!                if(CONS(I,J,K).NE.CONS(I,J,K)) print*, 'CONS', myid,I,J,K
!                if(APS(I,J,K) .NE.APS(I,J,K))  then
!                    print*, 'APS',  myid,I,J,K
!                    pause
!                endif
!                if(AIPS(I,J,K).NE.AIPS(I,J,K)) print*, 'AIPS', myid,I,J,K
!                if(AIMS(I,J,K).NE.AIMS(I,J,K)) print*, 'AIMS', myid,I,J,K
!                if(AJPS(I,J,K).NE.AJPS(I,J,K)) print*, 'AJPS', myid,I,J,K
!                if(AJMS(I,J,K).NE.AJMS(I,J,K)) print*, 'AJMS', myid,I,J,K
! !               do K1=1,KK
! !                   if(DSDY(K1,K,I,J).NE.DSDY(K1,K,I,J)) print*, 'DSDY', myid,I,J,K,K1,DSDY(K1,K,I,J)
! !               enddo
!            enddo
!            enddo
!            enddo
!             endif
            
            !solve the linear systems at each CV
            !DSDY=0.0d0      !RDB
            Call GaussSolver (CONS, APS, AIMS, AIPS, AJPS, AJMS, SpeciesMFO, SpeciesMF, DSDY, KK)


            DEALLOCATE (GAMS)
            DEALLOCATE (CONS)
            DEALLOCATE (APS)
            DEALLOCATE (AIMS)
            DEALLOCATE (AIPS)
            DEALLOCATE (AJPS)
            DEALLOCATE (AJMS)
            DEALLOCATE (DSDY)
            DEALLOCATE(WDD)
            DEALLOCATE(ScrubRate)
            DEALLOCATE(ScrubRateOld)
            DEALLOCATE(SMassF)
            DEALLOCATE(SMoleF)
            DEALLOCATE(DIFFU)
            DEALLOCATE(WDDP)
            
        end subroutine SolveSpecies
        !****************************************************************************************************

        !****************************************************************************************************
        SUBROUTINE SolveSoot
        !This subroutine solves the soot equations. When using the sectional code, the equations are solved in a locally coupled
        !manner (species are coupled to other species in a given CV, but not in a global sense). The resulting linear system
        !at each CV is solved using Gaussian Elimination. When using the 2-eq code, each equation is solved globally using the
        !TDMA for pentadiagonal matrices iterative algorithm.

            !local variable declarations
            double precision DIFF, ACOF
            double precision AREA, VOL, APT
            DOUBLE PRECISION, DIMENSION(:,:,:), Allocatable :: CS1, CNS1
            DOUBLE PRECISION, DIMENSION(:,:), Allocatable :: CS2, CNS2 
            DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE :: GAMS, APS, AIPS, AIMS, AJPS, AJMS, CONS
            DOUBLE PRECISION,DIMENSION(:,:,:,:),ALLOCATABLE :: DSDY
            DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: SOURCE, SOURCEOLD
            DOUBLE PRECISION :: SourceK, Pdiffusivity    
            !end local variable declarations
            
            !allocate matrices
            allocate(SOURCE(2*MSection))
            allocate(SOURCEOLD(2*MSection))
            allocate (GAMS(LM1,M1,2*MSection))
            allocate (CONS(LM1,M1,2*MSection))
            allocate (APS(LM1,M1,2*MSection))
            allocate (AIMS(LM1,M1,2*MSection))
            allocate (AIPS(LM1,M1,2*MSection))
            allocate (AJPS(LM1,M1,2*MSection))
            allocate (AJMS(LM1,M1,2*MSection))
            if(SootModel.EQ.1) then !only allocate DSDY if using sectional code
                allocate (DSDY(2*MSection,2*MSection,LM1,M1))
                allocate(CS1(MSection,MSection,MSection))
                allocate(CNS1(MSection,MSection,MSection))
                allocate(CS2(MSection,MSection))
                allocate(CNS2(MSection,MSection))
            endif     
            !allocate matrices                                  
            
            !Set the old species mass fractions to the current values, if using sectional code
            if(SootModel.Eq.1) SootSecO = SootSec

            !initialize GAMS, CONS, APS, AIMS, AIPS, AJPS, AJMS, DSDY
            GAMS = 0.0d0
            APS  = 0.0d0
            CONS = 0.0d0
            AIMS = 0.0d0
            AIPS = 0.0d0
            AJPS = 0.0d0
            AJMS = 0.0d0
            DSDY = 0.0d0

            !gets the values for GAMS, DSDY, and partially calculate CONS and APS
            
            !calculate GAMS at each CV
            DO J=1, M1
                DO I=1, L1
 
                  if (T(I,J).GT.3350.0) T(I,J) = 2950.0

                    if(.NOT.SOLIDCVP(I,J)) then !if not in a solid region        
                        
                        !set diffusion coefficient 
                        if(SootModel.EQ.1) then !for the section model, include ordinary diffusion
                            DO K=1, MSection
                                CALL Pdiff(K,Pdiffusivity)
                                GAMS(I,J,K) = Pdiffusivity*RHO(I,J)
                                GAMS(I,J,K+MSection) = GAMS(I,J,K)
                            ENDDO 
                        else !for the 2-eq code, do not include ordinary diffusion
                            GAMS(I,J,:) = 0.0d0
                        endif                                   
                        
                    else !if in a solid region, set diffusion velocities and GAMS to zero
                        GAMS(I,J,:) = 0.0d0
                    endif
                            
                ENDDO
            ENDDO
            
            !set GAMS to zero at boundary and CPU interface CVs
            if(myid==0)then
                GAMS(1,:,:)=0.0d0
            endif

            if(myid == (NUMP-1) ) then
                GAMS(L1,:,:)=0.0d0
            endif

            do i=1,L1
                GAMS(I,M1,:)=0.0d0
                GAMS(I,1,:) =0.0d0
            enddo
            !done calculating diffusion velocities and GAMS  
            
            !calculate DSDY (if using sectional code), and partially APS and CONS
            DO J=2, M2
                DO I=2, L2
                
                    CS1 = 0.0d0
                    CNS1 = 0.0d0
                    CS2 = 0.0d0
                    CNS2 = 0.0d0
                
                    if(.NOT.SOLIDCVP(I,J)) then !if not in a solid region 
                        
                        CALL SurfaceRates !get the species feedback rates, on a per section basis (only 1 section for 2-eq code)
                        
                        CALL NuclChem !get the nucleation rates, on a per species basis

                        !get the soot source term due to surface reactions/condensation/coagulation/sintering/obliteration/fragmentation for each section (stored in SOURCE)
                        SOURCE=0.d0
                        if(SootModel.EQ.1) then
                            CALL SootSecRate(SOURCE,.false.,CS1,CNS1,CS2,CNS2)
                        else
                            CALL Soot2eqRate(SOURCE) !get the source term for mass fraction and number density due to surface reaction (stored in SOURCE)
                        endif                                  
                        
                        !calculate the part of the source term due to thermothoresis and correctional diffusion velocities
                        DO K=1,2*MSection                                    
                            SourceK = Source(K) -(RHO(I+1,J)*(VSTX(I+1,J,K))*SootSec(I+1,J,K)- RHO(I-1,J)*(VSTX(I-1,J,K)) &
                            *SootSec(I-1,J,K))/(XDIF(I)+XDIF(I+1)) - (RHO(I,J+1)*(VSTY(I,J+1,K))*SootSec(I,J+1,K)*R(J+1)-   &
                            RHO(I,J-1)*(VSTY(I,J-1,K))*SootSec(I,J-1,K)*R(J-1))/((YDIF(J)+YDIF(J+1))*R(J)) 
                            CONS(I,J,K) = DMAX1(0.0D0,SourceK)
                            APS(I,J,K)=-DMAX1(0.0D0,-SourceK)/(SootSec(I,J,K)+smallnum) 
                        ENDDO
                        
                        !calculate DSDY, if using the sectional code 
                        if(SootModel.EQ.1) then
                        
                            SOURCEOLD = SOURCE !set SOURCEOLD to the unperturbed SOURCE values
                             
                            DO K2=1,2*MSection
                                !SootSec(I,J,K2)=SootSecO(I,J,K2)+(ABSO+ABSO*DABS(SootSecO(I,J,K2))) !perturb the current soot sectional value
                                SootSec(I,J,K2)=SootSecO(I,J,K2)+(ABSO**4+ABSO*DABS(SootSecO(I,J,K2)))      ! RDB
                                               
                                CALL SurfaceRates !get the species feedback rates, on a per section basis
                                
                                CALL SootSecRate(SOURCE,.true.,CS1,CNS1,CS2,CNS2) !get the soot source term due to surface reactions/condensation/sintering/obliteration/fragmentation for each section (use old Coag rates)          
                                
                                !calcualte DSDY
                                DO K1=1,2*MSection
                                    !DSDY(K1,K2,I,J) = -(SOURCEOLD(K1)-SOURCE(K1))/(ABSO+ABSO*DABS(SootSecO(I,J,K2)))
                                    DSDY(K1,K2,I,J) = -(SOURCEOLD(K1)-SOURCE(K1))/(ABSO**4+ABSO*DABS(SootSecO(I,J,K2)))     ! RDB se disminuyó la magnitud de la perturbación
                                ENDDO
                                
                                SootSec(I,J,K2) = SootSecO(I,J,K2) !set current soot sectional value back to the original value, as SootSecO never got modified
                            ENDDO         

                        endif
                        
                    else !if in a solid region
                        APS(I,J,:) = 0.0d0
                        CONS(I,J,:) = 0.0d0
                        if(SootModel.EQ.1) DSDY(:,:,I,J) = 0.0d0
                    endif
                
                ENDDO
            ENDDO
            !done calculating GAMS, DSDY, and partially calculating APS and CONS                 

            !calculate the matrix coefficients
            DO K=1,2*MSection
                include 'speciessootcoef.inc'
                        If(SootModel.Eq.1) then
                            CONS(I,J,K)=CONS(I,J,K)+APT*SootSec(I,J,K)             
                            APS(I,J,K)=(-APS(I,J,K)*VOL+AIPS(I,J,K)+AIMS(I,J,K)+AJPS(I,J,K)+AJMS(I,J,K))
                        else
                            CONS(I,J,K)=CONS(I,J,K)+APT*SootSec(I,J,K)
                            APS(I,J,K)=(-APS(I,J,K)*VOL+AIPS(I,J,K)+AIMS(I,J,K)+AJPS(I,J,K)+AJMS(I,J,K)) /RELAX(6)
                            CONS(I,J,K)=CONS(I,J,K)*VOL+(1.d0-RELAX(6))*APS(I,J,K)*SootSec(I,J,K)
                        endif      
                    ENDDO
                ENDDO
            ENDDO
            
            call MPI_BARRIER(MPI_COMM_WORLD,IERR)
            
            if(SootModel.EQ.1) then
                !solve the linear systems at each CV, using Gauss Elimination
                Call GaussSolver (CONS, APS, AIMS, AIPS, AJPS, AJMS, SootSecO, SootSec, DSDY, 2*MSection)
            else
                !solve the 2 soot equations in a globally coupled manner, using TDMA
                Call TDMASolver(APS(:,:,1),AIPS(:,:,1),AIMS(:,:,1),AJPS(:,:,1),AJMS(:,:,1),CONS(:,:,1),SootSec(:,:,1),2,2,NTIMES(5)) 
                Call TDMASolver(APS(:,:,2),AIPS(:,:,2),AIMS(:,:,2),AJPS(:,:,2),AJMS(:,:,2),CONS(:,:,2),SootSec(:,:,2),2,2,NTIMES(5))     
            endif
           
            DEALLOCATE (GAMS)
            DEALLOCATE (CONS)
            DEALLOCATE (APS)
            DEALLOCATE (AIMS)
            DEALLOCATE (AIPS)
            DEALLOCATE (AJPS)
            DEALLOCATE (AJMS)
            if(SootModel.EQ.1) then !only allocate DSDY if using sectional code
                DEALLOCATE (DSDY)
                DEALLOCATE (CS1)
                DEALLOCATE (CNS1)
                DEALLOCATE (CS2)
                DEALLOCATE (CNS2) 
            endif   
            DEALLOCATE(SOURCE)  
            DEALLOCATE(SOURCEOLD)             
            
        end subroutine SolveSoot
        !******************************************************************************************************
        
        !******************************************************************************************************
        SUBROUTINE SolveT
        !This subroutine solves the energy equation in a globally coupled manner. The coefficient matrix is formed
        !and solved using the TDMA for penta-diagonal matrices iterative algorithm.
                            
            !local variable declarations
            double precision DIFF, ACOF
            double precision AREA, VOL, APT
            DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE :: GAMS, APS, AIPS, AIMS, AJPS, AJMS, CONS
            DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE :: DOM_Inputs
            LOGICAL,DIMENSION(:,:),ALLOCATABLE :: DOM_SOLID_GLOB, DOM_SOLIDCV
          DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: COND, CPB, DOM_GLOB, DOM_SVF, DOM_XCO, DOM_XCO2, DOM_XH2O
            DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: DOM_QR, QR, QR_for_send
            DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: SMoleF, Xprint, Yprint, CP, enthalpy, WDD, ScrubRate
        DOUBLE PRECISION :: SootVF, SootMF, CPSOOT, HSOOT, S_X_DIFF, S_Y_DIFF, W_REAC, Soot_p_rate, Correct, Source
            DOUBLE PRECISION :: SootOxidation
            DOUBLE PRECISION ::  DCPBDX, DCPBDY
            DOUBLE PRECISION :: toto              !RDB totoratelimit
            integer :: ll
            DOUBLE PRECISION :: OMEGA2(IREAC) ! RDB tasas de formación
!            DOUBLE PRECISION, PARAMETER :: totoratelimit = 0.02d0     !RDB totoratelimit
            !end local variable declarations
            
            !allocate matrices
            allocate (GAMS(LM1,M1,1))
            allocate (CONS(LM1,M1,1))
            allocate (APS(LM1,M1,1))
            allocate (AIMS(LM1,M1,1))
            allocate (AIPS(LM1,M1,1))
            allocate (AJPS(LM1,M1,1))
            allocate (AJMS(LM1,M1,1))
            allocate(SMoleF(KK))
            allocate(CP(KK))
            allocate(WDD(KK))
            allocate(ScrubRate(KK))
            allocate(enthalpy(KK))
            allocate(COND(LM1,M1))
            allocate(CPB(0:LM1,M1)) !start array at zero for compatibility with CVGradients subroutine
            allocate(DOM_Inputs(5,M1,L1))
            allocate(DOM_SOLIDCV(M1,L1))         
            allocate(DOM_GLOB(L3*NUMP+2,M1))
            allocate(DOM_SVF(L3*NUMP+2,M1))
            allocate(DOM_XCO(L3*NUMP+2,M1))
            allocate(DOM_XCO2(L3*NUMP+2,M1))
            allocate(DOM_XH2O(L3*NUMP+2,M1))
            allocate(DOM_QR(L3*NUMP+2,M1))
            allocate(DOM_SOLID_GLOB(L3*NUMP+2,M1))
            allocate(Xprint(L3*NUMP+2))
            allocate(QR_for_send(M1,L2))
            allocate(QR(L2,M1))
            !allocate matrices
            
            !initialize arrays
            GAMS = 0.d0
            CONS = 0.d0
            APS = 0.d0
            AIMS = 0.d0
            AIPS = 0.d0
            AJPS = 0.d0
            AJMS = 0.d0
            
            !gets the values for GAMS, and partially calculate CONS and APS            
            !calculate GAMS
            DO I=1,L1
                DO J=1,M1                            

                  if (T(I,J).GT.3350.0) T(I,J) = 2950.0

                    if(.NOT.SOLIDCVP(I,J)) then !if not in a solid region
                        CALL CKYTX(SpeciesMF(I,J,:),IWORK,WORK,SMoleF) !get species mole fractions, store in SMoleF
                        CALL MCACON(T(I,J),SMoleF,TPWRK,COND(I,J)) !get mixture averaged conductivity ERG/CM*K*S
                        CALL CKCPBS(T(I,J),SpeciesMF(I,J,:),IWORK,WORK,CPB(I,J)) !get mixture averaged heat capacity ergs/(gm*K)
                        GAMS(I,J,1) = COND(I,J)/CPB(I,J) !assign GAMS
                    else !set conductivity and heat capacity to that of the solid
                        COND(I,J) = CONDSOLID
                        CPB(I,J) = CPSOLID
                        GAMS(I,J,1) = CONDSOLID/CPSOLID
                    endif
                    
                ENDDO
                
                GAMS(I,M1,1)= 0.0d0
                GAMS(I,1,1)=0.0d0
                
            ENDDO

            if(myid.eq.(NUMP-1)) GAMS(L1,:,1)=0.0d0
            !done calculating GAMS

            !calculate radiation terms using discrete ordinates method (DOM)  
            !make array of inputs for DOM subroutine on each CPU                                  
            DO J=1,M1
                DO I=1,L1  
                
                    DOM_Inputs(1,J,I) = T(I,J) !get GLOBerature 
                
                    if(.NOT.SOLIDCVP(I,J)) then !if not in a solid region                       
                        CALL CKYTX(SpeciesMF(I,J,:),IWORK,WORK,SMoleF) !get species mole fractions, store in SMoleF
                        DOM_Inputs(2,J,I) = SMoleF(ICO)
                        DOM_Inputs(3,J,I) = SMoleF(ICO2)
                        DOM_Inputs(4,J,I) = SMoleF(IH2O)
                        
                        !get soot volume fraction
                        if(SootModel.eq.1) then
                            call SootMassFrac(SootSec(i,j,1:MSection),XS,sootMF)
                            SootVF = sootMF*RHO(I,J)/DensityP
                        else
                            SootVF = SootSec(I,J,1)*Rho(I,J) / DensityP
                        endif
                        
                        DOM_Inputs(5,J,I) = SootVF
                    else
                        DOM_Inputs(2,J,I)=0.d0
                        DOM_Inputs(3,J,I)=0.d0
                        DOM_Inputs(4,J,I)=0.d0
                        DOM_Inputs(5,J,I)=0.d0
                    endif                               
                    
                    DOM_SOLIDCV(J,I) = SOLIDCVP(I,J)
                    
                ENDDO
            ENDDO  
            
            !send DOM input arrays to CPU-0, and assemble global DOM input array to pass to DOM subroutine
            if(myid.ne.0) then !if not on CPU-0, send DOM_Inputs to CPU-0
                CALL MPI_SEND(XU(2), L2, MPI_DOUBLE_PRECISION,0,100000+myid,MPI_COMM_WORLD,IERR)
                CALL MPI_SEND(DOM_Inputs(1,1,2), 5*L2*M1,MPI_DOUBLE_PRECISION,0,110000+myid,MPI_COMM_WORLD,IERR)
                CALL MPI_SEND(DOM_SOLIDCV(1,2), L2*M1,MPI_LOGICAL,0,120000+myid,MPI_COMM_WORLD,IERR)

            else !if you are CPU-0
            
                !assemble the part of the DOM input arrays from CPU-0
                do i=1,L1 !index of print in i direction
                    Xprint(i)=XU(i)
                    do j=1,M1
                        DOM_GLOB(i,j)=DOM_Inputs(1,J,I)
                        DOM_XCO(i,j)=DOM_Inputs(2,J,I)
                        DOM_XCO2(i,j)=DOM_Inputs(3,J,I)
                        DOM_XH2O(i,j)=DOM_Inputs(4,J,I)
                        DOM_SVF(i,j)=DOM_Inputs(5,J,I)
                        DOM_SOLID_GLOB(I,J) = DOM_SOLIDCV(J,I)
                    enddo
                enddo
!                
                do indp=1,(NUMP-1)  !receive the DOM_Inputs from the other processes into the DOM input arrays
                    CALL MPI_RECV(Xprint(L3*indp+2),L2,MPI_DOUBLE_PRECISION,indp,100000+indp,MPI_COMM_WORLD,STATUS,IERR)
                    CALL MPI_RECV(DOM_Inputs(1,1,2),5*L2*M1,MPI_DOUBLE_PRECISION,indp,110000+indp,MPI_COMM_WORLD,STATUS,IERR)
                    CALL MPI_RECV(DOM_SOLIDCV(1,2),L2*M1,MPI_LOGICAL,indp,120000+indp,MPI_COMM_WORLD,STATUS,IERR)
                    
                    !assemble DOM_Inputs from other CPUs into global arrays
                    do i=2,L1 !index of print in i direction
                        do j=1,M1
                            DOM_GLOB(L3*indp+i,j)=DOM_Inputs(1,J,I)
                            DOM_XCO(L3*indp+i,j)=DOM_Inputs(2,J,I)
                            DOM_XCO2(L3*indp+i,j)=DOM_Inputs(3,J,I)
                            DOM_XH2O(L3*indp+i,j)=DOM_Inputs(4,J,I)
                            DOM_SVF(L3*indp+i,j)=DOM_Inputs(5,J,I)
                            DOM_SOLID_GLOB(L3*indp+I,J) = DOM_SOLIDCV(J,I)
                        enddo
                    enddo               
                                        
                enddo
                
                ! call dom using the DOM global array
                CALL DOM(DOM_GLOB,DOM_XCO,DOM_XCO2,DOM_XH2O,DOM_SVF,DOM_QR,Xprint,YV,PRESSURE,DOM_SOLID_GLOB)             

            endif
            call MPI_BARRIER(MPI_COMM_WORLD,IERR)    
             
            !DOM radiation calculated on CPU-0
            
            !distribute the results of DOM from CPU-0 to all other CPUs
            if(myid.eq.0) then !if CPU-0, assign the relevant results of DOM_QR to QR
                DO J=2,M2
                    DO I=2,L2
                        QR(i,j)=DOM_QR(i,j)
                    enddo
                enddo
            endif

            do indp=1,(NUMP-1) !send DOM_QR to other CPU-0s
                if(myid.eq.0) then !if CPU-0, send the relevant part of DOM_QR via QR_for_send
                
                    DO I=2,L2
                        DO J=1,M1
                            QR_for_send(J,I) = DOM_QR(L3*indp+i,J)
                        ENDDO
                    ENDDO           

                    CALL MPI_SEND(QR_for_send(1,2), L3*M1, MPI_DOUBLE_PRECISION,indp,36000+indp,MPI_COMM_WORLD,IERR)

                else if(myid.eq.indp) then !receive the relevant part of DOM_QR into QR_for_send on CPU indp, and then store in QR
                
                    CALL MPI_RECV(QR_for_send(1,2), L3*M1,MPI_DOUBLE_PRECISION, 0,36000+indp,MPI_COMM_WORLD,STATUS,IERR)
                    
                    DO I=2,L2
                        DO J=1,M1
                            QR(I,J) = QR_for_send(J,I)
                        ENDDO
                    ENDDO                     
                               
                endif      
            enddo
            
            call MPI_BARRIER(MPI_COMM_WORLD,IERR)   
            !DOM results distributed 
            
            !partially calculate the source term (CONS), and APS
            DO I=2,L2
                DO J=2,M2
                    
                    if(.NOT.SOLIDCVP(I,J)) then !if not in a solid region
                        
                        CALL CKWYP(PRESSURE*PCKIN,T(I,J),SpeciesMF(I,J,:),IWORK,WORK,WDD) !get reaction rates moles/(cm**3*sec)
                        ! RDB test WDD
                        do K=1,KK
!                            if(WDD(K).NE.WDD(K))then
!                                print*,"WDD T NaN ",K,I,J,myid
!                                WDD(K) = 0.0d0      !RDB
!                            endif
                            if(WDD(K).GE.1.d50)then    !HUGE(0.d0)/
                                print*,"WDD T Inf ",K,I,J,myid
                                WDD(K) = 1.d50      !RDB   HUGE(0.d0)/
                            endif
                            if(WDD(K).LE.-1.d50)then   !HUGE(0.d0)/
                                print*,"WDD T -Inf",K,I,J,myid
                                WDD(K) = -1.d50     !RDB   HUGE(0.d0)/
                            endif
                        enddo
                        !
                        CALL SOOTHCP(T(I,J),CPSOOT,HSOOT) !get the enthalpy and specific heat capacity of soot
!                        print*,'HSOOT',HSOOT
                        
                        !Calculate the effect of soot on WDD, and get the soot production rate
                        
                        CALL SurfaceRates !get the species feedback rates, on a per section basis
                        
                        CALL NuclChem !get the nucleation rates
                            
                        CALL GasScrubbing(ScrubRate) !calculate the scrubbing rates         
                        
                        !update WDD to include the effect of gas phase scrubbing                            
                        WDD = WDD + ScrubRate
                        
                        !get the soot production rate, based on the values from SurfaceRates, and NuclChem
                        CALL SootProRate(Soot_p_rate,SootOxidation) !g-soot/cm^3/s
                        Soot_p_rate = Soot_p_rate - SootOxidation   !SootOxidation rate is returned as a positive value
                        DO K = 1, NumNuc
                            Soot_p_rate = Soot_p_rate + SootNucRates(I,J,NumNuc+K) !include contribution from nucleation
                        ENDDO     
                        !done calculating the effect of soot on WDD, and getting the soot production rate
                        
                        !calculate the flux of energy due to species diffusion, and energy release due to chemical reaction
                        !first, initialize diffusion energy flux, and enthalpy release to zero
                        S_X_DIFF = 0.0d0
                        S_Y_DIFF = 0.0d0
                        W_REAC  = 0.0d0
                        !get species heat capacity, and species enthalpy
                        CALL CKHML (T(I,J),IWORK,WORK,enthalpy)     !ergs/mole
                        CALL CKCPMS(T(I,J),IWORK,WORK,CP)           !ergs/(gm*K) Specific heats at constant pressure in mass units for the species
                        DO K=1,KK
                            S_X_DIFF = S_X_DIFF - RHO(I,J)*CP(K)*(VKX(I,J,K)+VCX(I,J)+VTX(I,J,K))*SpeciesMF(I,J,K)
                            S_Y_DIFF = S_Y_DIFF - RHO(I,J)*CP(K)*(VKY(I,J,K)+VCY(I,J)+VTY(I,J,K))*SpeciesMF(I,J,K)
                            !
                            !------ reaction rates limited --------- ! RDB
                            ! se utiliza el mismo limite que para las especies
                            toto = WDD(K)*WT(K)
                          !  print*,'toto'
                            if(dabs(toto).GT.totoratelimit.and.totoratelimit.ne.0.0d0)then
                                if(toto.gt.0.0d0) then
                                toto = totoratelimit
                                else
                                toto = -totoratelimit
                                endif
                                !totoFlag = 1
                            endif
                            !
!                if (myid.eq.0 .and. I.eq.2 .and. (J.eq.26 .or. J.eq.27 .or. J.eq.28) &
!                 .and. (K.eq.25 .or. K.eq.1 .or. K.eq.19)) toto=0.d0
!                 if ((K.eq.3 .or. K.eq.37 .or. K.eq.38)) toto=0.d0
                            !------ reaction rates limited --------- ! RDB
                            !W_REAC = W_REAC+ WDD(K)*Enthalpy(K)
                            W_REAC = W_REAC + toto/WT(K)*Enthalpy(K)       ! RDB totoratelimit
!                            if (myid.eq.0 .and. I.eq.2 .and. J.eq.27 .and. dabs(WDD(K)*Enthalpy(K)).gt.1.D8) &
!                            print*,'Enth',i,K,WDD(K)*Enthalpy(K),WDD(K),Enthalpy(K) !,WT(K)

                        ENDDO
!                        if (myid.eq.0 .and. I.eq.2 .and. J.eq.26) print*,'W_REAC',W_REAC

                        ! RDB test reacciones
!                        if (myid.eq.0 .and. I.eq.2 .and. J.eq.27 )then
!                        OMEGA2(:)=0.d0           ! tasas de formacion
!                        ! Rates of progress for the reactions
!                        ! cgs units - moles/(cm**3*sec)
!                        call CKQYP  (PRESSURE*PCKIN,T(I,J), SpeciesMF(I,J,:), IWORK, WORK, OMEGA2)
!                        do ll=1,size(OMEGA2)
!                        if (dabs(OMEGA2(ll)).gt.5.d-4) print*,'OMEGA2',ll,OMEGA2(ll)
!                        enddo
!                        endif

!                        if (myid.eq.0 .and. I.eq.2 .and. (J.eq.26 .or. J.eq.27 .or. J.eq.28)) then
!                        W_REAC= 0.0d0
!                       ! print*,'W_REAC',J
!                        endif
                             
                        !account for the effect of soot diffusion and production rate
                        !get soot mass fraction
                        do k=1, MSection
                            if(SootModel.eq.1) then
                                sootMF = DEXP(XS(K))*SootSec(i,j,K)
                                IF(sootMF.LT.0.0d0) sootMF = 0.0d0
                            else
                                sootMF = SootSec(I,J,1)
                            endif 
                            S_X_DIFF = S_X_DIFF - RHO(I,J)*CPSOOT*(VSX(I,J,K)+VSTX(I,J,K))*sootMF
                            S_Y_DIFF = S_Y_DIFF - RHO(I,J)*CPSOOT*(VSY(I,J,K)+VSTY(I,J,K))*sootMF
                        enddo      
                        W_REAC = W_REAC + (soot_p_rate/C_MW)*HSOOT   !HSOOT is ergs/mol carbons. must divide soot_p_rate by molecular weight of carbon to get
                                                                     !the number of mols of carbons in the soot.  
                    
                    else !if in a solid region
                        W_REAC = 0.0d0
                        S_X_DIFF = 0.0d0
                        S_Y_DIFF = 0.0d0
                    endif
                    
                    CALL CVGradients((1.d0/CPB),DCPBDX, DCPBDY,0) !get the gradient of heat capacity across the control volume 
                    
                    !calculate the final energy source term
                    CORRECT=-COND(I,J)*DTDX(I,J)*DCPBDX-COND(I,J)*DTDY(I,J)*DCPBDY
                    SOURCE=(-W_REAC+S_X_DIFF*DTDX(I,J)+S_Y_DIFF*DTDY(I,J)+QR(I,J)*10.0d0)/CPB(I,J)+CORRECT
                    !SOURCE=(-W_REAC+S_X_DIFF*DTDX(I,J)+S_Y_DIFF*DTDY(I,J)+QR(I,J)*10.0d0)/CPB(I,J)+CORRECT/2.d0 !RDB
!if (myid.eq.0 .and. I.eq.2 .and. J.eq.26) print*,'SOURCE',SOURCE,-W_REAC/CPB(I,J),S_X_DIFF*DTDX(I,J)/CPB(I,J),&
!S_Y_DIFF*DTDY(I,J)/CPB(I,J),QR(I,J)*10.0d0/CPB(I,J),CORRECT
!if (myid.eq.0 .and. I.eq.2 .and. J.eq.26) print*,'CORRECT',-COND(I,J),DTDX(I,J),DCPBDX,DTDY(I,J),DCPBDY,S_X_DIFF
                    
                    CONS(I,J,1)=DMAX1(0.0D0,SOURCE)
                    APS(I,J,1)=-DMAX1(0.0D0,-SOURCE)/(T(I,J)+smallnum)
                      
                ENDDO
            ENDDO
            !done calculating GAMS, and partially calculating APS, and CONS
            
            !calculate the matrix coefficients
            K=1 !using gamS, apS etc, due to similar structure for DIFF and FLOW to species/soot. Thus set K=1, and using same
                !matrix coefficient subroutine
            include 'speciessootcoef.inc'
                
                    CONS(I,J,K)=CONS(I,J,K)+APT*T(I,J)
                    APS(I,J,K)=(-APS(I,J,K)*VOL+AIPS(I,J,K)+AIMS(I,J,K)+AJPS(I,J,K)+AJMS(I,J,K)) /RELAX(4)
                    CONS(I,J,K)=CONS(I,J,K)*VOL+(1.d0-RELAX(4))*APS(I,J,K)*T(I,J)    
                
                ENDDO
            ENDDO 
            
            call MPI_BARRIER(MPI_COMM_WORLD,IERR)    
            
            !solve the energy equation
            Call TDMASolver(APS(:,:,1),AIPS(:,:,1),AIMS(:,:,1),AJPS(:,:,1),AJMS(:,:,1),CONS(:,:,1),T,2,2,NTIMES(4))
            
            !EXPORT DOM QR
            !DOM_QR4Print(L3*NUMP+2,M1)
            DO J=1,M1
                DO I=1,L3*NUMP+2
                        DOM_QR4Print(i,j)=DOM_QR(i,j)
                enddo
            enddo

!            if (myid.eq.0) print*,'T',T(2,23),T(2,24),T(2,25),T(2,26),T(2,27),T(2,28),T(2,29)

            DEALLOCATE(GAMS)
            DEALLOCATE(CONS)
            DEALLOCATE(APS)
            DEALLOCATE(AIMS)
            DEALLOCATE(AIPS)
            DEALLOCATE(AJPS)
            DEALLOCATE(AJMS)
            DEALLOCATE(SMoleF) 
            DEALLOCATE(CP)
            DEALLOCATE(WDD)
            DEALLOCATE(ScrubRate)
            DEALLOCATE(enthalpy)
            DEALLOCATE(COND)
            DEALLOCATE(CPB)
            DEALLOCATE(DOM_Inputs)
            DEALLOCATE(DOM_GLOB)
            DEALLOCATE(DOM_SVF)          
            DEALLOCATE(DOM_XCO)
            DEALLOCATE(DOM_XCO2)
            DEALLOCATE(DOM_XH2O)
            DEALLOCATE(DOM_QR)
            DEALLOCATE(DOM_SOLIDCV)       
            DEALLOCATE(DOM_SOLID_GLOB)                      
            DEALLOCATE(QR)
            DEALLOCATE(QR_for_send)
            DEALLOCATE(Xprint)  

        end subroutine SolveT
        !******************************************************************************************************                     
                                   
        !************************************************************************************************************
        SUBROUTINE OutputResultss
        !This subroutine outputs values pertaining to the current solution to the main output file (unit 6). It gives the
        !average values of U,V,P,PC,T, and Soot Volume Fraction over the whole domain. Additionally, it prints out the rate 
        !of change of these average values, weighted by the number of control volumes and the timestep size. It outputs
        !the maximum time a CPU spent on the current iteration, the average time a CPU spent on the current iteration, and the load balancing.
        !FInally, tt determines if the timestep should be increased due to the current rates of change of the solution variables,
        !or due to a change in the input.dat file
            
            !local variable declarations
            DOUBLE PRECISION :: T_to_sum, U_to_sum, V_to_sum, volsoot_to_sum, P_to_sum, PC_to_sum, A5_to_sum
            DOUBLE PRECISION :: T_av_sub, U_av_sub, V_av_sub, volsoot_av_sub, P_av_sub, PC_av_sub, A5_av_sub
            DOUBLE PRECISION :: T_av_whole, U_av_whole, V_av_whole, volsoot_av_whole, P_av_whole, PC_av_whole, A5_av_whole
            DOUBLE PRECISION :: slope_T, slope_fv, slope_U, slope_V, slope_P, slope_PC, slope_A5, TimeItn, Time_new
            DOUBLE PRECISION :: TimeItnMax, TimeItnMin, TimeItnSum, timebalance
            DOUBLE PRECISION :: maxslope, DTnew
            DOUBLE PRECISION :: sootMF, sootVF
            DOUBLE PRECISION, SAVE :: T_av_whole_old, volsoot_av_whole_old, U_av_whole_old, V_av_whole_old, P_av_whole_old
            DOUBLE PRECISION, SAVE :: PC_av_whole_old, A5_av_whole_old
            DOUBLE PRECISION, SAVE :: Time_old
            INTEGER :: counts, countreturn, countmax
            !end local variable declarations
            
            !initialize summation variables
            T_to_sum=0.d0
            U_to_sum=0.d0
            V_to_sum=0.d0
            volsoot_to_sum=0.d0
            P_to_sum = 0.d0
            PC_to_sum = 0.d0
            A5_to_sum = 0.d0
            
            !sum up the solution variables on the locally
            do i=2,L2
                do j=2,M2
                    T_to_sum=T_to_sum+T(I,J)
                    U_to_sum=U_to_sum+U(I,J)
                    V_to_sum=V_to_sum+V(I,J)
                    !get soot volume fraction
                    if(SootModel.eq.1) then
                        call SootMassFrac(SootSec(i,j,1:MSection),XS,sootMF)
                        SootVF = sootMF*RHO(I,J)/DensityP
                    else
                        SootVF = SootSec(I,J,1)*Rho(I,J) / DensityP
                    endif
                    volsoot_to_sum=volsoot_to_sum+SootVF
                    P_to_sum = P_to_sum + P(I,J) 
                    PC_to_sum = PC_to_sum + PC(I,J) 
                    A5_to_sum = A5_to_sum + SpeciesMF(I,J,IBAPYR)+ SpeciesMF(I,J,IBAPYRS)+ SpeciesMF(I,J,IBGHIF) 
                enddo
            enddo
            !done summing up the variables locally
            
            !get a local average value for the solution variables
            T_av_sub=T_to_sum/DFLOAT((L2-1)*(M2-1))
            U_av_sub=U_to_sum/DFLOAT((L2-1)*(M2-1))
            V_av_sub=V_to_sum/DFLOAT((L2-1)*(M2-1))
            volsoot_av_sub=volsoot_to_sum/DFLOAT((L2-1)*(M2-1))
            P_av_sub=P_to_sum/DFLOAT((L2-1)*(M2-1))
            PC_av_sub=PC_to_sum/DFLOAT((L2-1)*(M2-1)) 
            A5_av_sub=A5_to_sum/DFLOAT((L2-1)*(M2-1)) 
            !done getting local average values 
            
            !sum up all the local averages, and send them to CPU-0
            CALL MPI_REDUCE(T_av_sub,T_av_whole,1,MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,IERR)
            CALL MPI_REDUCE(U_av_sub,U_av_whole,1,MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,IERR)
            CALL MPI_REDUCE(V_av_sub,V_av_whole,1,MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,IERR)
            CALL MPI_REDUCE(volsoot_av_sub,volsoot_av_whole,1,MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,IERR)
            CALL MPI_REDUCE(P_av_sub,P_av_whole,1,MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,IERR)
            CALL MPI_REDUCE(PC_av_sub,PC_av_whole,1,MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,IERR) 
            CALL MPI_REDUCE(A5_av_sub,A5_av_whole,1,MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,IERR) 
            !done summing up all the local averages, and sending them to CPU-0   
            
            !divide the sums by the number of processes
            if (myid == 0) then
                T_av_whole=T_av_whole/DFLOAT(NUMP)
                U_av_whole=U_av_whole/DFLOAT(NUMP)
                V_av_whole=V_av_whole/DFLOAT(NUMP)
                volsoot_av_whole=volsoot_av_whole/DFLOAT(NUMP)
                P_av_whole=P_av_whole/DFLOAT(NUMP)
                PC_av_whole=PC_av_whole/DFLOAT(NUMP)
                A5_av_whole=A5_av_whole/DFLOAT(NUMP)
            endif      
            !divide the sums by the number of processes
            
            !broadcast these sums
            CALL MPI_BCAST(T_av_whole,1, MPI_double_precision,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(U_av_whole,1, MPI_double_precision,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(V_av_whole,1, MPI_double_precision,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(volsoot_av_whole,1, MPI_double_precision,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(P_av_whole,1, MPI_double_precision,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(PC_av_whole,1, MPI_double_precision,0,MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(A5_av_whole,1, MPI_double_precision,0,MPI_COMM_WORLD,IERR)
            !broadcast these sums
            
            !calculate the rate of change of these sums
            slope_U=(U_av_whole-U_av_whole_old)/(U_av_whole_old+smallnum)/DT
            slope_fv=(volsoot_av_whole-volsoot_av_whole_old)/(volsoot_av_whole_old+smallnum)/DT
            slope_V=(V_av_whole-V_av_whole_old)/(V_av_whole_old+smallnum)/DT
            slope_P=(P_av_whole-P_av_whole_old)/(P_av_whole_old+smallnum)/DT
            slope_PC=(PC_av_whole-PC_av_whole_old)/(PC_av_whole_old+smallnum)/DT
            slope_T=(T_av_whole-T_av_whole_old)/(T_av_whole_old+smallnum)/DT
            slope_A5=(A5_av_whole-A5_av_whole_old)/(A5_av_whole_old+smallnum)/DT
            !calculate the rate of change of these sums
            
            !determine how long this CPU took to get through the iteration
            call system_clock(countreturn,counts,countmax)
            Time_new = dble(countreturn)                   
            TimeItn = (Time_new - Time_old)/dble(counts)                
            
            !determine the min, max, and summation of iteration times
            Call MPI_ALLREDUCE(TimeItn,TimeItnMax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
            Call MPI_ALLREDUCE(TimeItn,TimeItnMin,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
            Call MPI_ALLREDUCE(TimeItn,TimeItnSum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

            !calculate load balancing and average time each CPU took
            timebalance=1.d0-(TimeItnMax-TimeItnMin)/TimeItnMin
            timebalance=100.d0*timebalance
            TimeItnSum=TimeItnSum/NUMP

            !output the reults to unit 6 (main output file)
            if (myid.eq.0)    then
                IF(mod(ITER,10).NE.0) GOTO 402
                WRITE(6,401)
                401 FORMAT(1X,'  ITER',5X,'U_avg',8X,'V_avg',7x,        &
                  'T_avg',7x,'FV_avg',7x,'P_avg',7x,'PC_avg',7x,'A5_avg',7x,    &
                  'slopeT',7x,'slopeFV',7x,'slopeU',7x,'slopeV',7x,'slopeP',7x, &
                   'slopePC',7x,'slopeA5',7x,'time',7x,'TimeAvg',7x,'TimeBal')

                402 WRITE(6,403) ITER,U_av_whole,V_av_whole,        &
                   T_av_whole,  &
                   volsoot_av_whole,P_av_whole,PC_av_whole,A5_av_whole,slope_T, &
                   slope_fv,slope_U,slope_V,slope_P,slope_PC,slope_A5,TimeItnMax, &
                   TimeItnSum,timebalance
                403 FORMAT(1X,I6,1P17E13.4)
           endif

            if(slope_T.NE.slope_T.OR.slope_U.NE.slope_U.OR.slope_V.NE.slope_V.OR. &
                slope_P.NE.slope_P.OR.slope_fv.NE.slope_fv.OR.slope_A5.NE.slope_A5) then
                LRUN = .FALSE.
                IF(MYID.EQ.0)WRITE(6,*)'Solution Diverged. Program Stopped' 
            endif

           !set all the _old values to the current values
            Time_old=Time_new
            U_av_whole_old=U_av_whole
            V_av_whole_old=V_av_whole
            T_av_whole_old = T_av_whole
            volsoot_av_whole_old=volsoot_av_whole
            P_av_whole_old=P_av_whole
            PC_av_whole_old=PC_av_whole 
            A5_av_whole_old=A5_av_whole
            
            DTITER = 1 + DTITER
            if(myid==0) then
                maxslope = slope_T
                if(DTITER.GE.miniter.AND.abs(maxslope).LE.slopetar) then
                    if(DT+DT.LE.DTmax) then
                        DT = DT+DT
                        WRITE(6,*)'Changing DT to ',DT,' ...'
                        DTITER = 0
                        if((slopetar/10).LE.minslopetar) then
                            slopetar = minslopetar
                        else
                            slopetar = slopetar/10.d0
                        endif
                    else
                        DT = DTmax
                        slopetar = 1.d-10
                        WRITE(6,*)'Changing DT to ',DT,' ...'
                    endif
                endif
            endif
            CALL MPI_BCAST(DT,1, MPI_double_precision,0,MPI_COMM_WORLD,IERR)
            
            !check if timestep has been changed in the input file, every saveint iterations
            if(mod(ITER,saveint)==0) then
                if(myid.eq.0) then
                    OPEN(UNIT=214,STATUS='old',FORM='FORMATTED',FILE='input.dat')
                    WRITE(6,*) 'Seeing if the timestep has been changed.'
                    read(214,*)DTnew
                    if(DTnew.GT.DT) then
                        DT = DTnew
                        WRITE(6,*)'Changing DT to ',DT,' ...'
                        DTITER = 0
                    else
                        WRITE(6,*)'Current MaxSlope = ',slope_T
                        WRITE(6,*)'Current DT = ',DT
                    endif
                    close(214)
                endif
            endif
            CALL MPI_BCAST(DT,1, MPI_double_precision,0,MPI_COMM_WORLD,IERR)
        
        end subroutine OutputResultss
        !************************************************************************************************                           
         
        !******************************************************************************************************
        SUBROUTINE OutputResultsf
        !This subroutine produces the restart file, two solution files (one for overall soot parameters, one for all other variables)  
        !Also, if the flag is set and the sectional model is being used, a file containing highly detailed sectional values will be 
        !written
        
            !local variable declarations
        DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: UT, VT
        DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE :: TotalSoot , TotalHACA, TotalMoleCon
        INTEGER, PARAMETER :: NumTotalSootVar = 16, TotalDetail = 6 , NumHacaData = 13, NumMoleCon = 10
        INTEGER, PARAMETER :: NumAxialVar = 27
        DOUBLE PRECISION :: SootMF, SootVol, TOTMASS, DNUM, TOTNUM, Diameter, TotDiffX, TotDiffY, PNX
        CHARACTER (24) :: FSPECNAME, FMASSNAME, FSOOTNAME, FSOLNAME, FHACANAME, FMCONNAME, FAXIALNAME
        DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE :: Xprint, Yaux
        DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE :: FsootCenter, TempCenter, CHradCent, FormCent, OxCent
        DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE :: NuclCent, CondCent, HACACent, OxO2Cent, OxOHCent
        DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE :: MaxSootLoc, FsootMax, TempMax, FormMax, OxMax
        DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE :: SootMax, NuclMax, CondMax, HACAMax, OxO2Max, OxOHMax
        DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE :: BetaSoot, AlphaNucl, AlphaCond, AlphaHACA
        DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE :: AlphaOxO2, AlphaOxOH, AlphaForm, AlphaOx
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: UTPrint, VTPrint, TPrint, RHOPrint,PPrint,PCPrint
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: PforPrint, PCforprint, TforPrint, RHOforPrint, QR4Print
      DOUBLE PRECISION,DIMENSION(:,:,:), ALLOCATABLE :: XSforPrint, SforPrint,SootforPrint,XSprint,Sprint,SootPrint
      DOUBLE PRECISION,DIMENSION(:,:,:), ALLOCATABLE :: TotSootPrint, TotHacaPrint, TotMolePrint, TotAxialPrint
      DOUBLE PRECISION,DIMENSION(:,:,:,:), ALLOCATABLE :: DetailSoot, DetailSootTotal
      CHARACTER (24) :: FRNAME      ! RDB tasas de formación
      DOUBLE PRECISION,DIMENSION(:,:,:), ALLOCATABLE :: RPrint, RforPrint,RLimitProcs       ! RDB tasas de formación
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: RLimitPrint                           ! RDB tasas de formación
      CHARACTER:: str4*4                ! tasas de formacion
      CHARACTER:: text_in*1000          ! tasas de formacion


            INTEGER :: Loop, ll, ii
            !end local variable declarations
            
            !allocate matrices
            allocate(Yaux(M1)) !A.Jerez C&F 2013 simulations
            !Centerline Vars
            allocate(FSootCenter(NUMP*L3+2))
            allocate(TempCenter(NUMP*L3+2))
            allocate(CHradCent(NUMP*L3+2))
            allocate(NuclCent(NUMP*L3+2))
            allocate(CondCent(NUMP*L3+2))
            allocate(HACACent(NUMP*L3+2))
            allocate(OxO2Cent(NUMP*L3+2))
            allocate(OxOHCent(NUMP*L3+2))
            allocate(FormCent(NUMP*L3+2))
            allocate(OxCent(NUMP*L3+2))
            !fSootMax flowline vars
            allocate(MaxSootLoc(NUMP*L3+2))
            allocate(FSootMax(NUMP*L3+2))
            allocate(TempMax(NUMP*L3+2))
            allocate(NuclMax(NUMP*L3+2))
            allocate(CondMax(NUMP*L3+2))
            allocate(HACAMax(NUMP*L3+2))
            allocate(OxO2Max(NUMP*L3+2))
            allocate(OxOHMax(NUMP*L3+2))
            allocate(FormMax(NUMP*L3+2))
            allocate(OxMax(NUMP*L3+2))
            !Integrates vars
            allocate(BetaSoot(NUMP*L3+2))
            allocate(AlphaNucl(NUMP*L3+2))
            allocate(AlphaCond(NUMP*L3+2))
            allocate(AlphaHACA(NUMP*L3+2))
            allocate(AlphaOxO2(NUMP*L3+2))
            allocate(AlphaOxOH(NUMP*L3+2))
            allocate(AlphaForm(NUMP*L3+2))
            allocate(AlphaOx(NUMP*L3+2))
            !Original Vars
            allocate(UT(M1,L1))
            allocate(VT(M1,L1))
            allocate(TotalSoot(NumTotalSootVar,M1,L1))
            allocate(TotalHACA(NumHacaData,M1,L1),TotalMoleCon(NumMoleCon,M1,L1))
            allocate(XSforPrint(KK,M1,L1))  !A.Jerez mod 2017
            allocate(SforPrint(KK,M1,L1))
            allocate(SootforPrint(2*MSection,M1,L1))
            allocate(PforPrint(M1,L1))
            allocate(PCforPrint(M1,L1))
            allocate(TforPrint(M1,L1))
            allocate(RHOforPrint(M1,L1))
            allocate(Xprint(NUMP*L3+2))
            allocate(UTPrint(M1,NUMP*L3+2))
            allocate(VTPrint(M1,NUMP*L3+2))
            allocate(TPrint(M1,NUMP*L3+2))
            allocate(QR4Print(M1,NUMP*L3+2))
            allocate(RHOPrint(M1,NUMP*L3+2))
            allocate(PPrint(M1,NUMP*L3+2))
            allocate(PCPrint(M1,NUMP*L3+2))
            allocate(SPrint(KK,M1,NUMP*L3+2))
            allocate(XSPrint(KK,M1,NUMP*L3+2))  !A.Jerez Mod 2017
            allocate(SootPrint(MSection*2,M1,NUMP*L3+2))
            allocate(TotSootPrint(NumTotalSootVar,M1,NUMP*L3+2))
            allocate(TotHacaPrint(NumHacaData,M1,NUMP*L3+2))
            allocate(TotMolePrint(NumMoleCon,M1,NUMP*L3+2))
            allocate(TotAxialPrint(NumAxialVar,M1,NUMP*L3+2))
            allocate(DetailSoot(TotalDetail,Msection,M1,L1))
            allocate(DetailSootTotal(TotalDetail,Msection,M1,NUMP*L3+2))

            if(reacdetail .eq. 1) then
            allocate(RforPrint(IREACP,M1,L1))           ! RDB tasas de formación
            allocate(RPrint(IREACP,M1,NUMP*L3+2))       ! RDB tasas de formación
            allocate(RLimitProcs(2,IREAC,NUMP+1))    ! RDB tasas de formación
            allocate(RlimitPrint(2,IREAC))           ! RDB tasas de formación
            endif
            !done allocating matrices
            

            DO J=1,M1
                Yaux(J) = Y(J) * 1.0d0
            END DO

            if(myid.eq.0) then
            DO J=1,M1
                DO I=1,L3*NUMP+2
                    QR4Print(j,I) = DOM_QR4Print(i,j)
                enddo
            enddo
            endif

            !populate the matrices taht will be sent to CPU-0 for printing to a file
            DO J=1,M1
                DO I=1,L1
                    !determine U and V velocities at the center of the P control volumes
                    UT(J,I)=(U(I+1,J)+U(I,J))/2.d0
                    if(J.NE.M1) then
                        VT(J,I)=(V(I,J+1)+V(I,J))/2.d0
                    else
                        VT(J,I) = V(I,M1)
                    endif
                    !if not in fuel tube region, convert mass fractions to mole fractions

                    !convert mass fractions to mole fractions and store in SforPrint
                    if(.NOT.SOLIDCVP(I,J)) then !if not in a solid region, convert mass fractions to mole fractions, and store in SforPrint
                        CALL CKYTX(SpeciesMF(I,J,:),IWORK,WORK,SforPrint(:,J,I))
                        XSforPrint(:,J,I) = SpeciesMF(I,J,:)
                    else !if in fuel tube region, mole fractions = 0.0
                        SforPrint(:,J,I) = 0.0d0
                        XSforPrint(:,J,I) = 0.0d0
                    endif 
                    
                    !set up the matrix that holds P, PC, T, and RHO
                    PforPrint(J,I) = P(I,J)
                    PCforPrint(J,I) = PC(I,J)
                    TforPrint(J,I) = T(I,J)
                    RHOforPrint(J,I) = RHO(I,J)
                    if(reacdetail .eq.1 ) RforPrint(:,J,I) = Rate(I,J,:)      ! RDB tasas de formación

                    DO K=1,2*MSection
                        SootforPrint(K,J,I)=SootSec(I,J,K)
                    ENDDO     
                ENDDO
            ENDDO     
            
            !deal with special cases for UT and VT
            if(myid.eq.0) then
                DO J=1,M1
                    UT(J,1)=U(1,J)
                ENDDO
            endif

            if(myid.eq.(NUMP-1)) then
                DO J=1,M1
                    UT(J,L1)=U(L1,J)
                ENDDO
            endif

            DO I=1,L1
                VT(1,I)=V(I,1)
            ENDDO
            !done calculating UT, VT, storing species mole fractions in SpeciesMFO

            !calculate soot vol fraction, average diameters and total #/cm^3 for primaries and aggregates, and total soot mass per cm^3.
            !Store these parameters and soot diffusional and thermophoretic velocities in TotalSoot
            if(SootModel.Eq.1) then  !if using the sectional model
                DO I=1,L1
                    DO J=1,M1

                        !first get soot mass fraction
                        call SootMassFrac(SootSec(i,j,1:MSection),XS,sootMF)
                        !convert to soot volume fraction
                        SootVol = SootMF*rho(i,j)/densityP
                        IF(SootVol.GT.1d-16) THEN
                            TotalSoot(1,J,I) = 1.0d6*SootVol !store soot volume fraction in ppm (parts per million)

                            !now get aggregate and primary particle values
                            TOTMASS = 0.d0
                            DNUM = 0.d0
                            TOTNUM = 0.d0

                            !calculate fractal equivalent aggregate diameters and total number and mass of aggregates
                            DO K=1, MSection
                                IF(MSection.EQ.1) then
                                    Diameter = (6.d0*DEXP(XS(K))/PI/densityP)**(1.d0/3.d0)   !fix this, the diameters of aggregates should not be based on solid spheres
                                ELSE
                                    PNX = SootSec(I,J,K+MSection)/(SootSec(I,J,K)+smallnum)  !number of primaries per aggregate
                                    IF(PNX.NE.0) THEN
                                        Diameter = (6.d0*(DEXP(XS(K))/PNX)/PI/densityP)**(1.d0/3.d0) !priamry particle diameter
                                    ELSE
                                        Diameter = 0.d0
                                    ENDIF
                                    Diameter = Diameter*(FVOL*PNX)**(1.d0/DFRCT) !fractal diameter
                                ENDIF
                                TOTMASS = TOTMASS + DEXP(XS(K))*SootSec(I,J,K)*rho(I,J) !g-soot/cm^3
                                DNUM = DNUM + Diameter*SootSec(I,J,K)*rho(I,J) !cm*(#aggregate/cm^3)
                                TOTNUM = TOTNUM + SootSec(I,J,K)*rho(I,J) !#aggregates/cm^3
                            ENDDO

                            TotalSoot(2,J,I) = TOTMASS !store total mass
                            TotalSoot(3,J,I) = TOTNUM !store total number of aggregates

                            TotalSoot(4,J,I) = 1.0d7*(6.d0*(TOTMASS/(TOTNUM+smallnum))/PI/densityP)**(1.d0/3.d0) !average fractal diameter of clusters on a mass basis, in nanometers
                            TotalSoot(5,J,I) = 1.0d7*DNUM / (TOTNUM+smallnum)   !average fractal diameter of aggregates on a number basis, in nanometers

                            !calculate volume equivalent aggregate diameters
                            DNUM = 0.d0
                            TOTNUM = 0.d0

                            DO K=1, MSection
                                Diameter = (6.d0*DEXP(XS(K))/PI/densityP)**(1.d0/3.d0)
                                TOTMASS = TOTMASS + DEXP(XS(K))*SootSec(I,J,K)*rho(I,J) !g-soot/cm^3
                                DNUM = DNUM + Diameter*SootSec(I,J,K)*rho(I,J) !cm*(#aggregate/cm^3)
                                !TOTNUM = TOTNUM + SootSec(I,J,K)*rho(I,J) !#aggregates/cm^3
                                TOTNUM = TOTNUM + SootSec(I,J,K)
                            ENDDO

                            TOTNUM = TOTNUM*RHO(I,J)

                            TotalSoot(6,J,I) = 1.0d7*(6.d0*(TOTMASS/(TOTNUM+smallnum))/PI/densityP)**(1.d0/3.d0) !average volume diameter of clusters on a mass basis, in nanometers
                            TotalSoot(7,J,I) = 1.0d7*DNUM / (TOTNUM+smallnum)   !average volume diameter of aggregates on a number basis, in nanometers

                            !calculate primary particle diameters and total number
                            DNUM = 0.d0
                            TOTNUM = 0.d0

                            DO K=1, MSection
                                Diameter = (6.d0*(SootSec(I,J,K)*DEXP(XS(K))/(SootSec(I,J,K+MSection) &
                                    +smallnum))/PI/densityP)**(1.d0/3.d0) !diameter of the primaries in section K
                                DNUM = DNUM + Diameter*SootSec(I,J,K+MSection)*rho(I,J)
                                !TOTNUM = TOTNUM + SootSec(I,J,K+MSection)*rho(I,J) !#primaries/cm^3
                                TOTNUM = TOTNUM + SootSec(I,J,K+MSection)
                            ENDDO

                            TOTNUM = TOTNUM*rho(I,J)

                            TotalSoot(8,J,I) = TOTNUM  !store total number of primaries (#primaries/cm^3)
                            TotalSoot(9,J,I) = 1.0d7*DNUM / (TOTNUM + smallnum)   !average diameter of primaries on a number basis, in nanometers
                            !done calculating number of primaries and average diameters

                            !calculate soot diffusional velocities
                            TotDiffX=0.d0
                            TotDiffY=0.d0
                            DO K=1, MSection
                                TotDiffX = TotDiffX + VSX(I,J,K)*SootSec(I,J,K)*DEXP(XS(K))/(sootMF+smallnum)
                                TotDiffY = TotDiffY + VSY(I,J,K)*SootSec(I,J,K)*DEXP(XS(K))/(sootMF+smallnum)
                            ENDDO
                            TotalSoot(10,J,I) = TotDiffX+VSTX(I,J,1) !store soot x-direction diffusional velocity
                            TotalSoot(11,J,I) = TotDiffY+VSTY(I,J,1) !store soot y-direction diffusional velocity
                        ELSE
                            TotalSoot(1:11,J,I) = 0.d0
                        ENDIF
                    ENDDO
                ENDDO
            else !if using the 2-eq Model

                DO I=1,L1
                    DO J=1,M1
                        TotalSoot(1,J,I) = 1.0d6*SootSec(I,J,1)/DensityP*RHO(I,J)  !store soot volume fraction (in ppm)
                        TotalSoot(2,J,I) = SootSec(I,J,1)*RHO(I,J)  !store g-soot/cm3
                        TotalSoot(3,J,I) = SootSec(I,J,2)*RHO(I,J) !store #/cm^3
                        TotalSoot(4,J,I) = ((6.0d0*SootSec(I,J,1))/(PI*densityP*SootSec(I,J,2)+smallnum))**(1.0d0/3.0d0)  !diameter of soot
                        TotalSoot(5,J,I) = ((6.0d0*SootSec(I,J,1))/(PI*densityP*SootSec(I,J,2)+smallnum))**(1.0d0/3.0d0)  !diameter of soot
                        TotalSoot(6,J,I) = ((6.0d0*SootSec(I,J,1))/(PI*densityP*SootSec(I,J,2)+smallnum))**(1.0d0/3.0d0)  !diameter of soot
                        TotalSoot(7,J,I) = ((6.0d0*SootSec(I,J,1))/(PI*densityP*SootSec(I,J,2)+smallnum))**(1.0d0/3.0d0)  !diameter of soot
                        TotalSoot(8,J,I) = SootSec(I,J,2)*RHO(I,J) !store #/cm^3
                        TotalSoot(9,J,I) = ((6.0d0*SootSec(I,J,1))/(PI*densityP*SootSec(I,J,2)+smallnum))**(1.0d0/3.0d0)  !diameter of soot
                        TotalSoot(10,J,I) = VSTX(I,J,1)  !soot thermal diffusion in the x-direction
                        TotalSoot(11,J,I) = VSTY(I,J,1)  !soot thermal diffusion in the y-direction
                    ENDDO
                ENDDO

            endif
            !done calculating soot parameters and storing in TotalSoot

            !calculate total Inception, HACA Surface Growth, Condensation, and OH/O2 Oxidation rates and store in TotalSoot in g-soot/cc/s
            if(SootModel.Eq.1) then  !if using the sectional model
                DO I=1,L1
                    DO J=1,M1
                        TotalSoot(12,J,I) = 0.d0
                        DO K=1,NumNuc
                            TotalSoot(12,J,I) = TotalSoot(12,J,I)+SootNucRates(I,J,NumNuc+K) !nucleation
                        ENDDO
                        DO K=1,13
                            TotalHACA(K,J,I) = SootHacaRates(I,J,K) !for print HACA reactions rates
                        ENDDO
                        DO k=1,10
                            TotalMoleCon(K,J,I) = SootMolConc(I,J,K) !for print soot mole concentration species
                        ENDDO
                        TotalSoot(13,J,I) = sum(SootGasRates(I,J,:,4))*2.d0*C_MW  !HACA C2H2 addition
                        if(Mechanism.LE.3) then
                            TotalSoot(14,J,I) = sum(SootGasRates(I,J,:,7))*20.d0*C_MW &
                              + sum(SootGasRates(I,J,:,8))*18.d0*C_MW + sum(SootGasRates(I,J,:,9))*20.d0*C_MW  !PAH condensation
                        else
                            TotalSoot(14,J,I) = sum(SootGasRates(I,J,:,7))*16.d0*C_MW   !PAH condensation A4=7
                        end if
                        TotalSoot(15,J,I) = sum(SootGasRates(I,J,:,5))*2.d0*C_MW  !O2 oxidation
                        TotalSoot(16,J,I) = sum(SootGasRates(I,J,:,6))*1.d0*C_MW  !OH oxidation
                    ENDDO
                ENDDO
            else !if using 2-eq code

                DO I=1,L1
                    DO J=1,M1
                        TotalSoot(12,J,I) = SootNucRates(I,J,1)*2.d0*C_MW + SootNucRates(I,J,2)*6.d0*C_MW !nucleation
                        TotalSoot(13,J,I) = SootGasRates(I,J,1,1)*2.d0*C_MW   !C2H2 addition
                        TotalSoot(14,J,I) = 0.d0                               !Condensation
                        TotalSoot(15,J,I) = SootGasRates(I,J,1,2)*2.d0*C_MW    !O2 Oxidation
                        TotalSoot(16,J,I) = SootGasRates(I,J,1,3)*1.d0*C_MW    !OH Oxidation
                    ENDDO
                ENDDO

            endif

            !done calculating the total rates of soot mass creation/destruction and storing in TotalSoot

            call MPI_BARRIER(MPI_COMM_WORLD,IERR)

            !start assembling the output matrices on CPU-0
            if(myid.ne.0) then !if not CPU-0
                    !send all required data to CPU-0
                    CALL MPI_SEND(x(2), L2, MPI_DOUBLE_PRECISION,0,12590+myid,MPI_COMM_WORLD,IERR)
                    CALL MPI_SEND(ut(1,2), L2*M1,MPI_DOUBLE_PRECISION,0,22690+myid,MPI_COMM_WORLD,IERR)
                    CALL MPI_SEND(vt(1,2), L2*M1,MPI_DOUBLE_PRECISION,0,32790+myid,MPI_COMM_WORLD,IERR)
                    CALL MPI_SEND(PforPrint(1,2),L2*M1,MPI_DOUBLE_PRECISION,0,42890+myid,MPI_COMM_WORLD,IERR)
                    CALL MPI_SEND(PCforPrint(1,2),L2*M1,MPI_DOUBLE_PRECISION,0,52890+myid,MPI_COMM_WORLD,IERR)
                    CALL MPI_SEND(TforPrint(1,2),L2*M1,MPI_DOUBLE_PRECISION,0,62890+myid,MPI_COMM_WORLD,IERR)
                    CALL MPI_SEND(RHOforPrint(1,2),L2*M1,MPI_DOUBLE_PRECISION,0,72890+myid,MPI_COMM_WORLD,IERR)
                    CALL MPI_SEND(SforPrint(1,1,2), L2*M1*KK,MPI_DOUBLE_PRECISION,0,82990+myid,MPI_COMM_WORLD,IERR)
                    CALL MPI_SEND(XSforPrint(1,1,2), L2*M1*KK,MPI_DOUBLE_PRECISION,0,92990+myid,MPI_COMM_WORLD,IERR)
                    CALL MPI_SEND(SootforPrint(1,1,2), L2*M1*MSection*2,MPI_DOUBLE_PRECISION,0,102990+myid,MPI_COMM_WORLD,IERR)
                    CALL MPI_SEND(TotalSoot(1,1,2), L2*M1*NumTotalSootVar,MPI_DOUBLE_PRECISION,0,112990+myid,MPI_COMM_WORLD,IERR)
                    CALL MPI_SEND(TotalHaca(1,1,2), L2*M1*NumHacaData,MPI_DOUBLE_PRECISION,0,122990+myid,MPI_COMM_WORLD,IERR)
                    CALL MPI_SEND(TotalMoleCon(1,1,2), L2*M1*NumMoleCon,MPI_DOUBLE_PRECISION,0,132990+myid,MPI_COMM_WORLD,IERR)
                    if(reacdetail .eq.1)then
                    CALL MPI_SEND(RforPrint(1,1,2), L2*M1*IREACP,MPI_DOUBLE_PRECISION,0,133990+myid,MPI_COMM_WORLD,IERR)    ! tasas de formación reacciones
                    CALL MPI_SEND(RateLimit(1,1), 2*IREAC,MPI_DOUBLE_PRECISION,0,134990+myid,MPI_COMM_WORLD,IERR)
                    endif       ! tasas de formación limites


            else !if you are CPU-0, assemble the data into output arrays

                !place the data from CPU-0 in the output arrays
                do i=1,L2 !index of print in i direction
                    Xprint(i)=x(i)

                    do j=1,M1
                        UTprint(j,i)=ut(j,i)
                        VTprint(j,i)=vt(j,i)
                        Tprint(j,i)=TforPrint(j,i)
                        RHOprint(j,i)=RHOforPrint(j,i)
                        Pprint(j,i)=PforPrint(j,i)
                        PCprint(j,i)=PCforPrint(j,i)
                        Sprint(:,j,i)=SforPrint(:,j,i)
                        XSprint(:,j,i)=XSforPrint(:,j,i)
                        Sootprint(:,j,i) = Sootforprint(:,j,i)
                        TotSootprint(:,j,i) = TotalSoot(:,j,i)
                        TotHacaprint(:,j,i) = TotalHaca(:,j,i)
                        TotMolePrint(:,j,i) = TotalMoleCon(:,j,i)
                        if(reacdetail .eq.1) Rprint(:,j,i) = RforPrint(:,j,i)        ! tasas de formación
                    enddo
                enddo
                if(reacdetail .eq.1)then
                RLimitProcs(1,:,1)=RateLimit(1,:)           ! tasas limites de formación: min
                RLimitProcs(2,:,1)=RateLimit(2,:)           ! tasas limites de formación: max
                endif

                !collect data from all other CPUs and place in output arrays
                do indp = 1, (NUMP-1)
                    CALL MPI_RECV(Xprint(L3*indp+2), L2, MPI_DOUBLE_PRECISION,indp,12590+indp,MPI_COMM_WORLD,STATUS,IERR)
                    CALL MPI_RECV(utprint(1,L3*indp+2), L2*M1,MPI_DOUBLE_PRECISION,indp,22690+indp,MPI_COMM_WORLD,STATUS,IERR)
                    CALL MPI_RECV(vtprint(1,L3*indp+2), L2*M1,MPI_DOUBLE_PRECISION,indp,32790+indp,MPI_COMM_WORLD,STATUS,IERR)
                    CALL MPI_RECV(Pprint(1,L3*indp+2), L2*M1,MPI_DOUBLE_PRECISION,indp,42890+indp,MPI_COMM_WORLD,STATUS,IERR)
                    CALL MPI_RECV(PCprint(1,L3*indp+2), L2*M1,MPI_DOUBLE_PRECISION,indp,52890+indp,MPI_COMM_WORLD,STATUS,IERR)
                    CALL MPI_RECV(Tprint(1,L3*indp+2), L2*M1,MPI_DOUBLE_PRECISION,indp,62890+indp,MPI_COMM_WORLD,STATUS,IERR)
                    CALL MPI_RECV(RHOprint(1,L3*indp+2), L2*M1,MPI_DOUBLE_PRECISION,indp,72890+indp,MPI_COMM_WORLD,STATUS,IERR)
                    CALL MPI_RECV(Sprint(1,1,L3*indp+2), L2*M1*KK,MPI_DOUBLE_PRECISION,indp,82990+indp,MPI_COMM_WORLD,STATUS,IERR)
                 CALL MPI_RECV(XSprint(1,1,L3*indp+2), L2*M1*KK,MPI_DOUBLE_PRECISION,indp,92990+indp,MPI_COMM_WORLD,STATUS,IERR)
         CALL MPI_RECV(Sootprint(1,1,L3*indp+2), L2*M1*MSection*2,MPI_DOUBLE_PRECISION,indp,102990+indp,MPI_COMM_WORLD,STATUS,IERR)
 CALL MPI_RECV(TotSootPrint(1,1,L3*indp+2), L2*M1*NumTotalSootVar,MPI_DOUBLE_PRECISION,indp,112990+indp,MPI_COMM_WORLD,STATUS,IERR)
 CALL MPI_RECV(TotHacaPrint(1,1,L3*indp+2), L2*M1*NumHacaData,MPI_DOUBLE_PRECISION,indp,122990+indp,MPI_COMM_WORLD,STATUS,IERR)
 CALL MPI_RECV(TotMolePrint(1,1,L3*indp+2), L2*M1*NumHacaData,MPI_DOUBLE_PRECISION,indp,132990+indp,MPI_COMM_WORLD,STATUS,IERR)
 if(reacdetail .eq.1)then
 CALL MPI_RECV(Rprint(1,1,L3*indp+2), L2*M1*IREACP,MPI_DOUBLE_PRECISION,indp,133990+indp,MPI_COMM_WORLD,STATUS,IERR)    ! tasas de formación
 CALL MPI_RECV(RLimitProcs(1,1,indp+1), 2*IREAC,MPI_DOUBLE_PRECISION,indp,134990+indp,MPI_COMM_WORLD,STATUS,IERR)    ! tasas de formación
 endif


                enddo
                !done collecting data

                ! ------- RDB tasas de formación -----
                ! recolección de los limites por cada procesador
                if(reacdetail .eq.1)then
                do ll = 1,IREAC
                RlimitPrint(1,ll) = RlimitProcs(1,ll,1)
                RlimitPrint(2,ll) = RlimitProcs(2,ll,1)
                do indp = 1, (NUMP-1)
                    RlimitPrint(1,ll) = min(RlimitPrint(1,ll),RlimitProcs(1,ll,indp+1))        ! RDB tasas de formación: min
                    RlimitPrint(2,ll) = max(RlimitPrint(2,ll),RlimitProcs(2,ll,indp+1))        ! RDB tasas de formación: max
                enddo
                enddo
                endif

                !#################################
                !## CENTERLINE VARS CALCULATION ##
                !##     A. JEREZ 2017 MOD       ##
                !#################################

                DO I=1,L3*NUMP+2
                    FSootCenter(I) =   TotSootPrint(1,2,I)
                    TempCenter(I)  =   Tprint(2,I)
                    CHradCent(I)   =   Sprint(25,2,I)
                    NuclCent(I)    =   TotSootPrint(12,2,I)
                    CondCent(I)    =   TotSootPrint(14,2,I)
                    HACACent(I)    =   TotSootPrint(13,2,I)
                    OxO2Cent(I)    =   TotSootPrint(15,2,I)
                    OxOHCent(I)    =   TotSootPrint(16,2,I)
                    FormCent(I)    =   NuclCent(I) + CondCent(I) + HACACent(I)
                    OxCent(I)      =   OxO2Cent(I) + OxOHCent(I)
                ENDDO

                !#################################
                !## fsootmax VARS CALCULATION   ##
                !##     A. JEREZ 2017 MOD       ##
                !#################################
                MaxSootLoc = maxloc(TotSootPrint(1,:,:),dim=1)
                !DO I=1,L3*NUMP+2
                !    write(6,*) 'I= ', I, '  ', MaxSootLoc(I)
                !ENDDO
                DO I=1,L3*NUMP+2
                    FSootMax(I)   =   TotSootPrint(1,MaxSootLoc(I),I)
                    TempMax(I)    =   Tprint(MaxSootLoc(I),I)
                    NuclMax(I)    =   TotSootPrint(12,MaxSootLoc(I),I)
                    CondMax(I)    =   TotSootPrint(14,MaxSootLoc(I),I)
                    HACAMax(I)    =   TotSootPrint(13,MaxSootLoc(I),I)
                    OxO2Max(I)    =   TotSootPrint(15,MaxSootLoc(I),I)
                    OxOHMax(I)    =   TotSootPrint(16,MaxSootLoc(I),I)
                    FormMax(I)    =    NuclMax(I) + CondMax(I) + HACAMax(I)
                    OxMax(I)      =    OxO2Max(I) + OxOHMax(I)
                ENDDO

                !#################################
                !## INTEGRATE VARS CALCULATION  ##
                !##     A. JEREZ 2017 MOD       ##
                !#################################

                DO I=1,L3*NUMP+2
                    DO J=1,M1-2
                        BetaSoot(I)= BetaSoot(I) + ((Yaux(J+1)-Yaux(J))*((Yaux(J)+Yaux(J+1))/2) &
                                    *TotSootPrint(1,J,I))*2.d0*3.14159d0
                        AlphaNucl(I)= AlphaNucl(I) + ((Yaux(J+1)-Yaux(J))*((Yaux(J)+Yaux(J+1))/2) &
                                    *TotSootPrint(12,J,I))*2.d0*3.14159d0
                        AlphaCond(I)= AlphaCond(I) + ((Yaux(J+1)-Yaux(J))*((Yaux(J)+Yaux(J+1))/2) &
                                    *TotSootPrint(14,J,I))*2.d0*3.14159d0
                        AlphaHACA(I)= AlphaHACA(I) + ((Yaux(J+1)-Yaux(J))*((Yaux(J)+Yaux(J+1))/2) &
                                    *TotSootPrint(13,J,I))*2.d0*3.14159d0
                        AlphaOxO2(I)= AlphaOxO2(I) + ((Yaux(J+1)-Yaux(J))*((Yaux(J)+Yaux(J+1))/2) &
                                    *TotSootPrint(15,J,I))*2.d0*3.14159d0
                        AlphaOxOH(I)= AlphaOxOH(I) + ((Yaux(J+1)-Yaux(J))*((Yaux(J)+Yaux(J+1))/2) &
                                    *TotSootPrint(16,J,I))*2.d0*3.14159d0
                        AlphaForm(I)= AlphaNucl(I) + AlphaCond(I) + AlphaHACA(I)
                        AlphaOx(I)  = AlphaOxO2(I) + AlphaOxOH(I)
                    ENDDO
                ENDDO




                !write the 2 main tecplot readable files (species data, and soot total data), and the solution file from which retarts can be done
                write(6,*) ' '
                !start with writing species data file
                write(6,*)'## Writing a tecplot file of the Species solution...'
                if(ITER.LE.9999999)WRITE(FSPECNAME,'("SPEC_",I7,".dat")')ITER
                if(ITER.LE.999999)WRITE(FSPECNAME,'("SPEC_0",I6,".dat")')ITER
                if(ITER.LE.99999)WRITE(FSPECNAME,'("SPEC_00",I5,".dat")')ITER
                if(ITER.LE.9999)WRITE(FSPECNAME,'("SPEC_000",I4,".dat")')ITER
                if(ITER.LE.999)WRITE(FSPECNAME,'("SPEC_0000",I3,".dat")')ITER
                if(ITER.LE.99)WRITE(FSPECNAME,'("SPEC_00000",I2,".dat")')ITER
                if(ITER.LE.9)WRITE(FSPECNAME,'("SPEC_000000",I1,".dat")')ITER

                OPEN (4,FILE=FSPECNAME,STATUS='UNKNOWN')

                write(4,*)'TITLE ="Species"'
                write(4,34) trim(species_header)!,trim(species_header2)
                write(4,*) ' ZONE T=', '"' ,'Species', '"' , ',', ' I =', M1,',' , ' J =', L3*NUMP+2,  ',' ,' F =POINT'

                DO I=1,L3*NUMP+2
                    DO J=1,M1
                        !A.Jerez C&F 2013 simulations
                        WRITE(4,33)Yaux(J),Xprint(I),UTprint(J,I),VTprint(J,I),Tprint(J,I), &
                                RHOprint(J,I),Pprint(J,I),PCprint(J,I),(Sprint(K,J,I),K=1,KK)
                    enddo
                enddo
                CLOSE(4)
                33  FORMAT(1X,1P150E13.4)   !format of SPEC file, should be KK + 8
   !             33  FORMAT(1X,1P250E13.4)   !format of SPEC file, should be KK + 8
                34  FORMAT(1X,A1600)        !SPEC header format. should be equal to species_header character variable
                !done writing species data fileq

                !start with writing Mass Fraction species data file
                write(6,*)'## Writing a tecplot file of Species Mass Fraction solution...'
                if(ITER.LE.9999999)WRITE(FMASSNAME,'("MASS_",I7,".dat")')ITER
                if(ITER.LE.999999)WRITE(FMASSNAME,'("MASS_0",I6,".dat")')ITER
                if(ITER.LE.99999)WRITE(FMASSNAME,'("MASS_00",I5,".dat")')ITER
                if(ITER.LE.9999)WRITE(FMASSNAME,'("MASS_000",I4,".dat")')ITER
                if(ITER.LE.999)WRITE(FMASSNAME,'("MASS_0000",I3,".dat")')ITER
                if(ITER.LE.99)WRITE(FMASSNAME,'("MASS_00000",I2,".dat")')ITER
                if(ITER.LE.9)WRITE(FMASSNAME,'("MASS_000000",I1,".dat")')ITER

                OPEN (846,FILE=FMASSNAME,STATUS='UNKNOWN')

                write(846,*)'TITLE ="MF_Species"'
                write(846,34) trim(species_header)!,trim(species_header2)
                write(846,*) ' ZONE T=', '"' ,'Species', '"' , ',', ' I =', M1,',' , ' J =', L3*NUMP+2,  ',' ,' F =POINT'

                DO I=1,L3*NUMP+2
                    DO J=1,M1
                        WRITE(846,998) Yaux(J),Xprint(I),UTprint(J,I),VTprint(J,I),Tprint(J,I), &
                       ! WRITE(846,33) Yaux(J),Xprint(I),UTprint(J,I),VTprint(J,I),Tprint(J,I), &
                               RHOprint(J,I),Pprint(J,I),PCprint(J,I),(XSprint(K,J,I),K=1,KK)
                    enddo
                enddo
                CLOSE(846)
                998  FORMAT(1X,1P250E13.4)   !format of SPEC file, should be KK + 8
                !334  FORMAT(1X,A1900)        !SPEC header format. should be equal to species_header character variable
                !done writing species data fileq

                !------------- tasas de formación ---------------
                if(reacdetail .eq.1)then
                !start with writing formation rate for each reaction data file
                write(6,*)'## Writing a tecplot file of formation rate for each reaction solution...'
                if(ITER.LE.9999999)WRITE(FRNAME,'("R_",I7,".dat")')ITER
                if(ITER.LE.999999)WRITE(FRNAME,'("R_0",I6,".dat")')ITER
                if(ITER.LE.99999)WRITE(FRNAME,'("R_00",I5,".dat")')ITER
                if(ITER.LE.9999)WRITE(FRNAME,'("R_000",I4,".dat")')ITER
                if(ITER.LE.999)WRITE(FRNAME,'("R_0000",I3,".dat")')ITER
                if(ITER.LE.99)WRITE(FRNAME,'("R_00000",I2,".dat")')ITER
                if(ITER.LE.9)WRITE(FRNAME,'("R_000000",I1,".dat")')ITER

                OPEN (847,FILE=FRNAME,STATUS='UNKNOWN')

                write(847,*)'TITLE ="formation rates for each reaction"'
                write(847,*) 'VARIABLES ="r (cm)","z (cm)"'
                text_in=' '
                do ll=1,size(iR)
                    write(str4,'(I4)') iR(ll)
                    text_in = trim(adjustl(text_in)) //',"R'// trim(adjustl(str4)) //'"'
                    if (mod(ll,40)==0) then
                        write(847,*) text_in
                        !print*,text_in
                        text_in=' '
                    endif
                enddo
                write(847,*) text_in
                !write(847,*) 'VARIABLES ="r (cm)","z (cm)","R10","R19","R47"'
                write(847,*) ' ZONE T=', '"' ,'Species', '"' , ',', ' I =', M1,',', ' J =', L3*NUMP+2, ',' ,' F =POINT'

                DO I=1,L3*NUMP+2
                    DO J=1,M1
                        WRITE(847,998)Yaux(J),Xprint(I),(Rprint(K,J,I),K=1,IREACP)
                    enddo
                enddo
                CLOSE(847)
!                998  FORMAT(1X,1P150E13.4)   !format of SPEC file, should be KK + 8
                !334  FORMAT(1X,A1900)        !SPEC header format. should be equal to species_header character variable
                !done writing formation rates data file

                !start with writing formation rate limits for each reaction data file
                write(6,*)'## Writing a tecplot file of formation rate for each reaction solution...'
                if(ITER.LE.9999999)WRITE(FRNAME,'("RL_",I7,".dat")')ITER
                if(ITER.LE.999999)WRITE(FRNAME,'("RL_0",I6,".dat")')ITER
                if(ITER.LE.99999)WRITE(FRNAME,'("RL_00",I5,".dat")')ITER
                if(ITER.LE.9999)WRITE(FRNAME,'("RL_000",I4,".dat")')ITER
                if(ITER.LE.999)WRITE(FRNAME,'("RL_0000",I3,".dat")')ITER
                if(ITER.LE.99)WRITE(FRNAME,'("RL_00000",I2,".dat")')ITER
                if(ITER.LE.9)WRITE(FRNAME,'("RL_000000",I1,".dat")')ITER

                OPEN (848,FILE=FRNAME,STATUS='UNKNOWN')

                write(848,*) 'TITLE ="formation rates limits for each reaction"'
                write(848,*) 'VARIABLES ="N","min","max"'
                write(848,*) 'ZONE T=', '"' ,'Limits', '"' , ',', ' I =', IREAC ,',',' F =POINT'
                do ll=1,IREAC
                    WRITE(848,999) ll,RlimitPrint(1,ll),RlimitPrint(2,ll)
                enddo
                write(848,*)' '
                do indp = 0, (NUMP-1)
                    write(848,*) '# indices of domain ',L3*indp+2,L3*(indp+1)+1
                    write(str4,'(I4)') indp+1
                    write(848,*) 'ZONE T=', '"D'//trim(adjustl(str4))//'"' , ',', ' I =', IREAC ,',',' F =POINT'
                    do ll=1,IREAC
                        WRITE(848,999) ll,RlimitProcs(1,ll,indp+1),RlimitProcs(2,ll,indp+1)
                    enddo
                enddo
                CLOSE(848)
                999 FORMAT(1X,I6,1X,1P2E15.4)
                endif
                !done writing formation rates limit data file
                !------------- tasas de formación ---------------

                !write soot file
                if(ITER.LE.9999999)WRITE(FSOOTNAME,'("soot_",I7,".dat")')ITER
                if(ITER.LE.999999)WRITE(FSOOTNAME,'("soot_0",I6,".dat")')ITER
                if(ITER.LE.99999)WRITE(FSOOTNAME,'("soot_00",I5,".dat")')ITER
                if(ITER.LE.9999)WRITE(FSOOTNAME,'("soot_000",I4,".dat")')ITER
                if(ITER.LE.999)WRITE(FSOOTNAME,'("soot_0000",I3,".dat")')ITER
                if(ITER.LE.99)WRITE(FSOOTNAME,'("soot_00000",I2,".dat")')ITER
                if(ITER.LE.9)WRITE(FSOOTNAME,'("soot_000000",I1,".dat")')ITER

                OPEN (UNIT=42,FILE=FSOOTNAME,STATUS='UNKNOWN')

                write(42,*)'TITLE ="sootoutp"'
                write(42,*)'VARIABLES ="r (cm)","z (cm)","vol_frac","Tot_Mass","Num_agg",'
                write(42,*) '"Df_av_mass", "Df_av_num","Dv_av_mass", "Dv_av_num","Num_prime",'
                write(42,*) '"Dp_av","Usoot","Vsoot","IncepRate","HACARate","CondRate",'
                write(42,*) '"O2_OxRate","OH_OxRate","QR_rad"'
                write(42,*) ' ZONE T=', '"' ,'sootoutp', '"' , ',', ' I =', M1,',' , ' J =', L3*NUMP+2,  ',' ,' F =POINT' 

                DO I=1,L3*NUMP+2
                    DO J=1,M1
                        WRITE(42,900) Yaux(J),Xprint(I),(TotSootPrint(K,J,I),k=1,NumTotalSootVar),QR4Print(j,I)
                    enddo
                enddo
                close(42)
                900 FORMAT(1X,1P19E15.4)    !soot file format
                !done writing soot file

                write(6,*)'## Writing a tecplot file of the HACA rates...'
                !write HACA soot file
                if(ITER.LE.9999999)WRITE(FHACANAME,'("HACA_",I7,".dat")')ITER
                if(ITER.LE.999999)WRITE(FHACANAME,'("HACA_0",I6,".dat")')ITER
                if(ITER.LE.99999)WRITE(FHACANAME,'("HACA_00",I5,".dat")')ITER
                if(ITER.LE.9999)WRITE(FHACANAME,'("HACA_000",I4,".dat")')ITER
                if(ITER.LE.999)WRITE(FHACANAME,'("HACA_0000",I3,".dat")')ITER
                if(ITER.LE.99)WRITE(FHACANAME,'("HACA_00000",I2,".dat")')ITER
                if(ITER.LE.9)WRITE(FHACANAME,'("HACA_000000",I1,".dat")')ITER

                OPEN (UNIT=42,FILE=FHACANAME,STATUS='UNKNOWN')

                write(42,*)'TITLE ="HACAoutp"'
                write(42,*)'VARIABLES ="r (cm)","z (cm)","k1(H)","k-1(H2)","k2(OH)",'
                write(42,*) '"k-2(H2O)", "k3(H)","k4(C2H2)", "k5(O2)","X_soot",'
                write(42,*) '"SR1","SR2","SR3","SR4","SR5"'
                write(42,*) ' ZONE T=', '"' ,'HACAoutp', '"' , ',', ' I =', M1,',' , ' J =', L3*NUMP+2,  ',' ,' F =POINT'

                DO I=1,L3*NUMP+2
                    DO J=1,M1
                        WRITE(42,930) Yaux(J),Xprint(I),(TotHacaPrint(K,J,I),k=1,NumHacaData)
                    enddo
                enddo
                close(42)
                930 FORMAT(1X,1P15E15.4)    !Haca file format
                !done writing Haca file

                 write(6,*)'## Writing a tecplot file of the Mole Concentration...'
                !write Mole Concentration file
                if(ITER.LE.9999999)WRITE(FMCONNAME,'("MCON_",I7,".dat")')ITER
                if(ITER.LE.999999)WRITE(FMCONNAME,'("MCON_0",I6,".dat")')ITER
                if(ITER.LE.99999)WRITE(FMCONNAME,'("MCON_00",I5,".dat")')ITER
                if(ITER.LE.9999)WRITE(FMCONNAME,'("MCON_000",I4,".dat")')ITER
                if(ITER.LE.999)WRITE(FMCONNAME,'("MCON_0000",I3,".dat")')ITER
                if(ITER.LE.99)WRITE(FMCONNAME,'("MCON_00000",I2,".dat")')ITER
                if(ITER.LE.9)WRITE(FMCONNAME,'("MCON_000000",I1,".dat")')ITER

                OPEN (UNIT=42,FILE=FMCONNAME,STATUS='UNKNOWN')

                write(42,*)'TITLE ="MCONoutp"'
                write(42,*)'VARIABLES ="r (cm)","z (cm)","H","H2","OH",'
                write(42,*) '"H2O", "C2H2","O2", "BAPYR","BGHIF",'
                write(42,*) '"BAPYRS","A4"'
                write(42,*) ' ZONE T=', '"' ,'MCONoutp', '"' , ',', ' I =', M1,',' , ' J =', L3*NUMP+2,  ',' ,' F =POINT'

                DO I=1,L3*NUMP+2
                    DO J=1,M1
                        WRITE(42,932) Yaux(J),Xprint(I),(TotMolePrint(K,J,I),k=1,NumMoleCon)
                    enddo
                enddo
                close(42)
                932 FORMAT(1X,1P13E15.4)    !Mole Concentration file format
                !done writing Mole Concentration file

                !start with writing Axial Results data file
                write(6,*)'## Writing a tecplot file of axial profile solutions...'
                if(ITER.LE.9999999)WRITE(FAXIALNAME,'("AXIA_",I7,".dat")')ITER
                if(ITER.LE.999999)WRITE(FAXIALNAME,'("AXIA_0",I6,".dat")')ITER
                if(ITER.LE.99999)WRITE(FAXIALNAME,'("AXIA_00",I5,".dat")')ITER
                if(ITER.LE.9999)WRITE(FAXIALNAME,'("AXIA_000",I4,".dat")')ITER
                if(ITER.LE.999)WRITE(FAXIALNAME,'("AXIA_0000",I3,".dat")')ITER
                if(ITER.LE.99)WRITE(FAXIALNAME,'("AXIA_00000",I2,".dat")')ITER
                if(ITER.LE.9)WRITE(FAXIALNAME,'("AXIA_000000",I1,".dat")')ITER

                OPEN (92,FILE=FAXIALNAME,STATUS='UNKNOWN')
                write(92,*)'TITLE ="AxialProfiles"'
                write(92,*)'VARIABLES ="z (cm)", "FsootCenter", "TempCenter", "CHradCent", "FormCent",'
                write(92,*) '"OxCent", "NuclCent", "CondCent","HACACent", "OxO2Cent", "OxOHCent",'
                write(92,*) '"FsootMax", "TempMax", "FormMax", "OxMax",'
                write(92,*) '"NuclMax", "CondMax", "HACAMax","OxO2Max", "OxOHMax",'
                write(92,*) '"BetaSoot", "AlphaNucl", "AlphaCond", "AlphaHACA",'
                write(92,*) '"AlphaOxO2", "AlphaOxOH", "AlphaForm", "AlphaOx"'
                write(92,*) ' ZONE T=', '"' ,'Axial-Results', '"' , ',', ' I =', L3*NUMP+2,',' ,' F =POINT'

                DO I=1,L3*NUMP+2
                        !A.Jerez MOD 2017
                        WRITE(92,942)Xprint(I),FsootCenter(I), TempCenter(I), CHradCent(I), &
                               FormCent(I), OxCent(I), NuclCent(I), CondCent(I),  &
                               HACACent(I), OxO2Cent(I), OxOHCent(I), &
                               FsootMax(I), TempMax(I), FormMax(I), OxMax(I), &
                               NuclMax(I), CondMax(I), HACAMax(I), &
                               OxO2Max(I), OxOHMax(I), &
                               BetaSoot(I), AlphaNucl(I), AlphaCond(I), AlphaHACA(I), &
                               AlphaOxO2(I), AlphaOxOH(I), AlphaForm(I), AlphaOx(I)
                enddo
                942  FORMAT(1X,1P30E15.4)   !format of Axial profile results
                CLOSE(92)
                !done writing axial results

                !write restart file
                write(6,*)'## Writing a restart file for the current solution...'
                if(ITER.LE.9999999)WRITE(FSOLNAME,'("SOLN_",I7,".dat")')ITER
                if(ITER.LE.999999)WRITE(FSOLNAME,'("SOLN_0",I6,".dat")')ITER
                if(ITER.LE.99999)WRITE(FSOLNAME,'("SOLN_00",I5,".dat")')ITER
                if(ITER.LE.9999)WRITE(FSOLNAME,'("SOLN_000",I4,".dat")')ITER
                if(ITER.LE.999)WRITE(FSOLNAME,'("SOLN_0000",I3,".dat")')ITER
                if(ITER.LE.99)WRITE(FSOLNAME,'("SOLN_00000",I2,".dat")')ITER
                if(ITER.LE.9)WRITE(FSOLNAME,'("SOLN_000000",I1,".dat")')ITER


                OPEN(2,FILE=FSOLNAME,STATUS='UNKNOWN')
                
                write(2,*) iter  
                
                DO I=1,L3*NUMP+2
                    DO J=1,M1
                        WRITE(2,35)Y(J),Xprint(I),UTprint(J,I),VTprint(J,I),Tprint(J,I), &
                                Pprint(J,I),PCprint(J,I),RHOprint(J,I),(Sprint(K,J,I),K=1,KK), &
                                (SootPrint(K,J,I),K=1,MSection*2)
                    enddo
                enddo
                close(2)
                35  FORMAT(1X,1P172E26.15)  !SOLN file format, should be equal to KK+MSectionTOTAL+8
    !            35  FORMAT(1X,1P309E26.15)  !SOLN file format, should be equal to KK+MSectionTOTAL+8
                !done writing solution file
                !done writing files

            endif !done writing files from CPU-0
            
            !if input flag is set, and the sectional model is used, output a file with per section information about soot.
            if((SootDetail.Eq.1).AND.(SootModel.EQ.1)) then
                
                !populate the arrays to be sent to CPU-0
                DO I = 2, L1
                    DO J = 1, M1
                            
                        DO K = 1, MSection
                            if(myid.NE.0) then
                                !number of aggregates in section k (#/cm^3)
                                DetailSoot(1,K,J,I) = SootSec(I,J,K)*rho(I,J)
                                !mass of aggregates in section k (g/cm^3)
                                DetailSoot(2,K,J,I) = DEXP(XS(K))*SootSec(I,J,K)*rho(I,J)
                                !number of primaries in section k (#/cm^3)
                                DetailSoot(3,K,J,I) = SootSec(I,J,K+MSection)*rho(I,J)
                                !diameter of primaries in section k (#/cm^3)
                                DetailSoot(4,K,J,I) = 1.0d7*(6.d0*(SootSec(I,J,K)*DEXP(XS(K)) &
                                        /(SootSec(I,J,K+MSection)+smallnum))/PI/densityP)**(1.d0/3.d0)
                                !fractal aggregate diameter in section k
                                DetailSoot(5,K,J,I) = DetailSoot(4,K,J,I)*(FVOL*(SootSec(I,J,K+MSection) &
                                      /(SootSec(I,J,K)+smallnum)))**(1.d0/DFRCT)
                                !volume equivalent aggregate diameters in section k
                                DetailSoot(6,K,J,I) = 1.0d7*(6.d0*DEXP(XS(K))/PI/densityP)**(1.d0/3.d0)
                            else
                                !number of aggregates in section k (#/cm^3)
                                DetailSootTotal(1,K,J,I) = SootSec(I,J,K)*rho(I,J)
                                !mass of aggregates in section k (g/cm^3)
                                DetailSootTotal(2,K,J,I) = DEXP(XS(K))*SootSec(I,J,K)*rho(I,J)
                                !number of primaries in section k (#/cm^3)
                                DetailSootTotal(3,K,J,I) = SootSec(I,J,K+MSection)*rho(I,J)
                                !diameter of primaries in section k (#/cm^3)
                                DetailSootTotal(4,K,J,I) = 1.0d7*(6.d0*(SootSec(I,J,K)*DEXP(XS(K)) &
                                     /(SootSec(I,J,K+MSection)+smallnum))/PI/densityP)**(1.d0/3.d0)
                                !fractal aggregate diameter in section k
                                DetailSootTotal(5,K,J,I) = DetailSootTotal(4,K,J,I) &
                                     *(FVOL*(SootSec(I,J,K+MSection)/(SootSec(I,J,K) &
                                     +smallnum)))**(1.d0/DFRCT)
                                !volume equivalent aggregate diameters in section k
                                DetailSootTotal(6,K,J,I) = 1.0d7*(6.d0*DEXP(XS(K))/PI/densityP)**(1.d0/3.d0)
                            endif                
                        ENDDO
                    
                    ENDDO
                ENDDO                                             
                !assemble the output matrix on CPU-0
                if(myid.ne.0) then !if not cpu-0, send data to CPU-0
                
                CALL MPI_SEND(DetailSoot(1,1,1,2), L2*M1*MSection*TotalDetail, &
                                MPI_DOUBLE_PRECISION,0,115990+myid,MPI_COMM_WORLD,IERR)
                
                else !if you are CPU-0, assemble the output arrays and write the file 
                
                    do indp = 1, (NUMP-1)
                        CALL MPI_RECV(DetailSootTotal(1,1,1,L3*indp+2), &
                                L2*M1*MSection*TotalDetail,MPI_DOUBLE_PRECISION, &
                                indp,115990+indp,MPI_COMM_WORLD,STATUS,IERR)
                    enddo

                    !write detailed soot file
                    write(6,*)'## Writing detailed sectional soot model ouput file....'
                    if(ITER.LE.9999999)WRITE(FSOOTNAME,'("sootdetail_",I7,".dat")')ITER
                    if(ITER.LE.999999)WRITE(FSOOTNAME,'("sootdetail_0",I6,".dat")')ITER
                    if(ITER.LE.99999)WRITE(FSOOTNAME,'("sootdetail_00",I5,".dat")')ITER
                    if(ITER.LE.9999)WRITE(FSOOTNAME,'("sootdetail_000",I4,".dat")')ITER
                    if(ITER.LE.999)WRITE(FSOOTNAME,'("sootdetail_0000",I3,".dat")')ITER
                    if(ITER.LE.99)WRITE(FSOOTNAME,'("sootdetail_00000",I2,".dat")')ITER
                    if(ITER.LE.9)WRITE(FSOOTNAME,'("sootdetail_000000",I1,".dat")')ITER

                    OPEN (UNIT=52,FILE=FSOOTNAME,STATUS='UNKNOWN')

                    if(MSection.eq.35) then
                        write(52,*)'TITLE ="detailed"'
                        write(52,*)'VARIABLES ="r (cm)","z (cm)","A1","A2","A3","A4","A5",'
                        write(52,*)'"A6","A7","A8","A9","A10","A11","A12","A13","A14","A15",'
                        write(52,*)'"A16","A17","A18","A19","A20","A21","A22","A23","A24",'
                        write(52,*)'"A25","A26","A27","A28","A29","A30","A31","A32","A33",'
                        write(52,*)'"A34","A35","M1","M2","M3","M4","M5","M6","M7","M8",'
                        write(52,*)'"M9","M10","M11","M12","M13","M14","M15","M16","M17",'
                        write(52,*)'"M18","M19","M20","M21","M22","M23","M24","M25","M26",'
                        write(52,*)'"M27","M28","M29","M30","M31","M32","M33","M34","M35",'
                        write(52,*)'"P1","P2","P3","P4","P5","P6","P7","P8",'
                        write(52,*)'"P9","P10","P11","P12","P13","P14","P15","P16","P17",'
                        write(52,*)'"P18","P19","P20","P21","P22","P23","P24","P25","P26",'
                        write(52,*)'"P27","P28","P29","P30","P31","P32","P33","P34","P35",'
                        write(52,*)'"Dp1","Dp2","Dp3","Dp4","Dp5","Dp6","Dp7","Dp8",'
                        write(52,*)'"Dp9","Dp10","Dp11","Dp12","Dp13","Dp14","Dp15","Dp16","Dp17",'
                        write(52,*)'"Dp18","Dp19","Dp20","Dp21","Dp22","Dp23","Dp24","Dp25","Dp26",'
                        write(52,*)'"Dp27","Dp28","Dp29","Dp30","Dp31","Dp32","Dp33","Dp34","Dp35",'
                        write(52,*)'"Df1","Df2","Df3","Df4","Df5","Df6","Df7","Df8",'
                        write(52,*)'"Df9","Df10","Df11","Df12","Df13","Df14","Df15","Df16","Df17",'
                        write(52,*)'"Df18","Df19","Df20","Df21","Df22","Df23","Df24","Df25","Df26",'
                        write(52,*)'"Df27","Df28","Df29","Df30","Df31","Df32","Df33","Df34","Df35",'
                        write(52,*)'"Dv1","Dv2","Dv3","Dv4","Dv5","Dv6","Dv7","Dv8",'
                        write(52,*)'"Dv9","Dv10","Dv11","Dv12","Dv13","Dv14","Dv15","Dv16","Dv17",'
                        write(52,*)'"Dv18","Dv19","Dv20","Dv21","Dv22","Dv23","Dv24","Dv25","Dv26",'
                        write(52,*)'"Dv27","Dv28","Dv29","Dv30","Dv31","Dv32","Dv33","Dv34","Dv35"'
                        write(52,*) ' ZONE T=', '"' ,'detailed', '"' , ',', &
                          ' I =', M1,',' , ' J =', L3*NUMP+2,  ',' ,' F =POINT'
                          
                    elseif(MSection.EQ.16) then
                        write(52,*)'TITLE ="detailed"'
                        write(52,*)'VARIABLES ="r (cm)","z (cm)","A1","A2","A3","A4","A5",'
                        write(52,*)'"A6","A7","A8","A9","A10","A11","A12","A13","A14","A15",'
                        write(52,*)'"A16","M1","M2","M3","M4","M5","M6","M7","M8",'
                        write(52,*)'"M9","M10","M11","M12","M13","M14","M15","M16",'
                        write(52,*)'"P1","P2","P3","P4","P5","P6","P7","P8",'
                        write(52,*)'"P9","P10","P11","P12","P13","P14","P15","P16",'
                        write(52,*)'"Dp1","Dp2","Dp3","Dp4","Dp5","Dp6","Dp7","Dp8",'
                        write(52,*)'"Dp9","Dp10","Dp11","Dp12","Dp13","Dp14","Dp15","Dp16",'
                        write(52,*)'"Df1","Df2","Df3","Df4","Df5","Df6","Df7","Df8",'
                        write(52,*)'"Df9","Df10","Df11","Df12","Df13","Df14","Df15","Df16",'
                        write(52,*)'"Dv1","Dv2","Dv3","Dv4","Dv5","Dv6","Dv7","Dv8",'
                        write(52,*)'"Dv9","Dv10","Dv11","Dv12","Dv13","Dv14","Dv15","Dv16",'
                        write(52,*) ' ZONE T=', '"' ,'detailed', '"' , ',', &
                              ' I =', M1,',' , ' J =', L3*NUMP+2,  ',' ,' F =POINT'
                    endif


                    DO I=1,L3*NUMP+2
                        DO J=1,M1
                            !WRITE(52,910) Y(J),Xprint(I),(DetailSootTotal(1,K,J,I),k=1,MSection), &
                            WRITE(52,910) Yaux(J),Xprint(I),(DetailSootTotal(1,K,J,I),k=1,MSection), &
                                                         (DetailSootTotal(2,K,J,I),k=1,MSection), &
                                                         (DetailSootTotal(3,K,J,I),k=1,MSection), &
                                                         (DetailSootTotal(4,K,J,I),k=1,MSection), &
                                                         (DetailSootTotal(5,K,J,I),k=1,MSection), &
                                                         (DetailSootTotal(6,K,J,I),k=1,MSection)
                        enddo
                    enddo
                    close(52)
                    910 FORMAT(1X,1P218E15.4)   !sootdetail file format
                    !done writing detailed soot file
                
                endif 
                
            endif       
            !done writing detailed soot file
            
            deallocate(UT)
            deallocate(VT)
            deallocate(TotalSoot)
            deallocate(Xprint)
            deallocate(UTPrint)
            deallocate(VTPrint)
            deallocate(TPrint)
            deallocate(RHOPrint)
            deallocate(PPrint)
            deallocate(PCPrint)
            deallocate(SPrint)
            deallocate(SootPrint)
            deallocate(TotSootPrint)  
            deallocate(TotalHACA)
            deallocate(TotalMoleCon)
            deallocate(PforPrint)
            deallocate(PCforPrint)
            deallocate(TforPrint)
            deallocate(RHOforPrint) 
            deallocate(Yaux)
            deallocate(RLimitProcs)    ! tasas de formación
            deallocate(RlimitPrint)    ! tasas de formación
       
        end subroutine OutputResultsf                                     
        !*************************************************************************************************************************
        
        !*****************************************************************************************************************************************************
        SUBROUTINE TDMASolver (PA, N, S, E, W, B, Solution,IST,JST,NumIters)
        !This subroutine takes in the coefficients of a penta-diagonal matrix from each CPU (AP,N,S,E,W,B), and the current solution, and compiles them into global 
        !matrices. The global penta-diagonal matrices are then solved using the iterative TDMA for penta-diagonal matrices algorithm. The solution is 
        !then distributed to each CPU
        
            !local variable declarations
            DOUBLE PRECISION, DIMENSION(:,:) :: PA, N, S, E, W, B
            DOUBLE PRECISION, DIMENSION(0:,:) :: Solution
            INTEGER :: IST, JST, NumIters
            DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: GLOBP, GLOBN, GLOBS, GLOBE, GLOBW, GLOBB, GLOBSOL
            DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: LOCP, LOCN, LOCS, LOCE, LOCW, LOCB, LOCSOL
            DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PT, QT
            INTEGER :: ISTF, JSTF, IT1, IT2, JT1, JT2, II, JJ, NTT
            DOUBLE PRECISION :: temp, denom
            !done local variable declarations
            
            !allocate matrices
            allocate(GLOBP(2:M2,L3*NUMP+2))
            allocate(GLOBN(2:M2,L3*NUMP+2))
            allocate(GLOBS(2:M2,L3*NUMP+2))
            allocate(GLOBE(2:M2,L3*NUMP+2))
            allocate(GLOBW(2:M2,L3*NUMP+2))
            allocate(GLOBB(2:M2,L3*NUMP+2))
            allocate(GLOBSOL(M1,L3*NUMP+2))
            allocate(PT(max(M1,L3*NUMP+2)))
            allocate(QT(max(M1,L3*NUMP+2)))
            allocate(LOCP(2:M2,L2))
            allocate(LOCN(2:M2,L2))
            allocate(LOCS(2:M2,L2))
            allocate(LOCE(2:M2,L2))
            allocate(LOCW(2:M2,L2))
            allocate(LOCB(2:M2,L2))
            allocate(LOCSOL(M1,0:L0))
            !done allocating matrices    
            
            !set up loop sizes for TDMA
            ISTF=IST-1
            JSTF=JST-1

            IT1=L2+(NUMP-1)*L3+IST
            IT2=L3+(NUMP-1)*L3+IST

            JT1=M2+JST
            JT2=M3+JST
            !done setting up loop sizes
            
            !populate "LOC" matrices, which will be used to send coefficients to CPU-0 
            if(myid.NE.0) then
                DO I = 2,L2
                    DO J=2,M2
                        LOCP(J,I) = PA(I,J)
                        LOCN(J,I) = N(I,J)
                        LOCS(J,I) = S(I,J)
                        LOCE(J,I) = E(I,J)
                        LOCW(J,I) = W(I,J)
                        LOCB(J,I) = B(I,J)
                    ENDDO
                ENDDO
                
                DO I = 1,L1
                    DO J=1,M1
                        LOCSOL(J,I) = SOLUTION(I,J)
                    ENDDO
                ENDDO
            endif      
            !done populating LOC matrices         
            
            !assemble the global matrices
            if (myid .ne. 0) then !if not CPU-0, send matrix to CPU-0
                CALL MPI_SEND(LOCP(2,2), L3*M3, MPI_DOUBLE_PRECISION,0,61000+myid,MPI_COMM_WORLD,IERR)
                CALL MPI_SEND(LOCB(2,2), L3*M3, MPI_DOUBLE_PRECISION,0,62000+myid,MPI_COMM_WORLD,IERR)
                CALL MPI_SEND(LOCS(2,2), L3*M3, MPI_DOUBLE_PRECISION,0,63000+myid,MPI_COMM_WORLD,IERR)
                CALL MPI_SEND(LOCW(2,2), L3*M3, MPI_DOUBLE_PRECISION,0,64000+myid,MPI_COMM_WORLD,IERR)
                CALL MPI_SEND(LOCN(2,2), L3*M3, MPI_DOUBLE_PRECISION,0,65000+myid,MPI_COMM_WORLD,IERR)
                CALL MPI_SEND(LOCE(2,2), L3*M3, MPI_DOUBLE_PRECISION,0,66000+myid,MPI_COMM_WORLD,IERR)
                CALL MPI_SEND(LOCSOL(1,2), L2*M1, MPI_DOUBLE_PRECISION,0,67000+myid,MPI_COMM_WORLD,IERR)

            else !if you are CPU-0

                !assemble the global arrays, starting with data from CPU-0
                do j=2,M2
                    do i=IST,L2
                        GLOBP(J,I)=PA(i,j)
                        GLOBB(J,I)=B(i,j)
                        GLOBS(J,I)=S(i,j)
                        GLOBW(J,I)=W(i,j)
                        GLOBN(J,I)=N(i,j)
                        GLOBE(J,I)=E(i,j)
                    enddo
                enddo

                do j=1,M1
                    do i=1,L2
                        GLOBSOL(J,I)=SOLUTION(i,j)
                    enddo
                enddo

                !receive data from other CPUs
                do indp=1,(NUMP-1)
                    CALL MPI_RECV(GLOBP(2,L3*indp+2), L3*M3,  MPI_DOUBLE_PRECISION,indp,61000+indp,MPI_COMM_WORLD,STATUS,IERR)
                    CALL MPI_RECV(GLOBB(2,L3*indp+2), L3*M3,  MPI_DOUBLE_PRECISION,indp,62000+indp,MPI_COMM_WORLD,STATUS,IERR)
                    CALL MPI_RECV(GLOBS(2,L3*indp+2), L3*M3,  MPI_DOUBLE_PRECISION,indp,63000+indp,MPI_COMM_WORLD,STATUS,IERR)
                    CALL MPI_RECV(GLOBW(2,L3*indp+2), L3*M3,  MPI_DOUBLE_PRECISION,indp,64000+indp,MPI_COMM_WORLD,STATUS,IERR)
                    CALL MPI_RECV(GLOBN(2,L3*indp+2), L3*M3,  MPI_DOUBLE_PRECISION,indp,65000+indp,MPI_COMM_WORLD,STATUS,IERR)
                    CALL MPI_RECV(GLOBE(2,L3*indp+2), L3*M3,  MPI_DOUBLE_PRECISION,indp,66000+indp,MPI_COMM_WORLD,STATUS,IERR)
                    CALL MPI_RECV(GLOBSOL(1,L3*indp+2), L2*M1, MPI_DOUBLE_PRECISION,indp,67000+indp,MPI_COMM_WORLD,STATUS,IERR) 
                    
                enddo
                !done receiving data, and assembling global arrays 

                !start the TDMA iterative solution process
                DO NTT=1,NumIters
                ! *-----------------------------------------------------------------------
                    DO J=JST,M2
                        PT(ISTF)=0.d0
                        QT(ISTF)=GLOBSOL(J,ISTF)
                        DO I=IST,L2+L3*(NUMP-1)
                            DENOM=GLOBP(J,I)-PT(I-1)*GLOBS(J,I)
                            PT(I)=GLOBN(J,I)/DENOM
                            TEMP=GLOBB(J,I)+GLOBE(J,I)*GLOBSOL(J+1,I)+ GLOBW(J,I)*GLOBSOL(J-1,I)
                            QT(I)=(TEMP+GLOBS(J,I)*QT(I-1))/DENOM
                        enddo
                        DO  II=IST,L2+L3*(NUMP-1)
                            I=IT1-II
                            GLOBSOL(J,I)=GLOBSOL(J,I+1)*PT(I)+QT(I)
                        enddo
                    enddo
                ! *-----------------------------------------------------------------------
                    DO  JJ=JST,M3
                        J=JT2-JJ
                        PT(ISTF)=0.d0
                        QT(ISTF)=GLOBSOL(J,ISTF)
                        DO  I=IST,L2+L3*(NUMP-1)
                            DENOM=GLOBP(J,I)-PT(I-1)*GLOBS(J,I)
                            PT(I)=GLOBN(J,I)/DENOM
                            TEMP=GLOBB(J,I)+GLOBE(J,I)*GLOBSOL(J+1,I)+ GLOBW(J,I)*GLOBSOL(J-1,I)
                            QT(I)=(TEMP+GLOBS(J,I)*QT(I-1))/DENOM
                        enddo
                        DO  II=IST,L2+L3*(NUMP-1)
                            I=IT1-II
                            GLOBSOL(J,I)=GLOBSOL(J,I+1)*PT(I)+QT(I)
                        enddo
                    enddo
                ! *-----------------------------------------------------------------------
                    DO  I=IST,L2+L3*(NUMP-1)
                        PT(JSTF)=0.d0
                        QT(JSTF)=GLOBSOL(JSTF,I)
                        DO  J=JST,M2
                            DENOM=GLOBP(J,I)-PT(J-1)*GLOBW(J,I)
                            PT(J)=GLOBE(J,I)/DENOM
                            TEMP=GLOBB(J,I)+GLOBN(J,I)*GLOBSOL(J,I+1)+ GLOBS(J,I)*GLOBSOL(J,I-1)
                            QT(J)=(TEMP+GLOBW(J,I)*QT(J-1))/DENOM
                        enddo
                        DO  JJ=JST,M2
                            J=JT1-JJ
                            GLOBSOL(J,I)=GLOBSOL(J+1,I)*PT(J)+QT(J)
                        enddo
                    enddo
                ! *-----------------------------------------------------------------------
                    DO II=IST,L3+L3*(NUMP-1)
                        I=IT2-II
                        PT(JSTF)=0.d0
                        QT(JSTF)=GLOBSOL(JSTF,I)
                        DO J=JST,M2
                            DENOM=GLOBP(J,I)-PT(J-1)*GLOBW(J,I)
                            PT(J)=GLOBE(J,I)/DENOM
                            TEMP=GLOBB(J,I)+GLOBN(J,I)*GLOBSOL(J,I+1)+ GLOBS(J,I)*GLOBSOL(J,I-1)
                            QT(J)=(TEMP+GLOBW(J,I)*QT(J-1))/DENOM
                        enddo
                        DO  JJ=JST,M2
                            J=JT1-JJ
                            GLOBSOL(J,I)=GLOBSOL(J+1,I)*PT(J)+QT(J)
                        enddo
                    enddo
                ! ************************************************************************
                !done 1 TDMA iteration          
                enddo
                !done TDMA iteration

            endif

            !distribute the new solution from CPU-0 to all other CPUs
            if(myid.eq.0) then !if CPU-0, take part of the GLOBSOL, then pass parts of GLOBSOL to all other CPUs
                do j=jst,M2
                    do i=IST,L0
                        Solution(i,j)=GLOBSOL(J,I)
                    enddo
                enddo

                !send data to each CPU, except for the last one
                do indp=1,(NUMP-2)
                    CALL MPI_SEND(GLOBSOL(1,L3*indp), LM1*M1, MPI_DOUBLE_PRECISION,indp,61000+indp,MPI_COMM_WORLD,IERR)
                enddo

                !send data to the last CPU
                CALL MPI_SEND(GLOBSOL(1,L3*(NUMP-1)), L1*M1, MPI_DOUBLE_PRECISION,NUMP-1,61000+NUMP-1,MPI_COMM_WORLD,IERR)

            else !if not CPU-0, receive the data sent from CPU-0

                if(myid.NE.(NUMP-1)) then
                    CALL MPI_RECV(LOCSOL(1,0), LM1*M1,  MPI_DOUBLE_PRECISION,0,61000+myid,MPI_COMM_WORLD,STATUS,IERR)
                else
                    CALL MPI_RECV(LOCSOL(1,0), L1*M1,  MPI_DOUBLE_PRECISION,0,61000+myid,MPI_COMM_WORLD,STATUS,IERR)
                endif

            endif
            !solution distributed

            !Move new solution from LOCSOL to SOLUTION
            IF(myid.NE.0.AND.myid.NE.(NUMP-1)) THEN
                DO I = 0, L0
                    DO J = JST, M1
                        SOLUTION(I,J) = LOCSOL(J,I)
                    ENDDO
                ENDDO
            ENDIF
            
            IF(myid.EQ.(NUMP-1)) THEN
                DO I = 0, L2 
                    DO J = JST, M1
                        SOLUTION(I,J) = LOCSOL(J,I)
                    ENDDO
                ENDDO
            ENDIF                
        
            DEALLOCATE(GLOBP)
            DEALLOCATE(GLOBN)
            DEALLOCATE(GLOBS)
            DEALLOCATE(GLOBE)
            DEALLOCATE(GLOBW)
            DEALLOCATE(GLOBB)
            DEALLOCATE(GLOBSOL)
            DEALLOCATE(PT)
            DEALLOCATE(QT)
            DEALLOCATE(LOCP)
            DEALLOCATE(LOCN)
            DEALLOCATE(LOCS)
            DEALLOCATE(LOCE)
            DEALLOCATE(LOCW)
            DEALLOCATE(LOCB)
            DEALLOCATE(LOCSOL) 
            
            call MPI_BARRIER(MPI_COMM_WORLD,IERR)     

        end subroutine TDMASolver
        !*********************************************************************************************************************************

        !*****************************************************************************************************************************************************
        SUBROUTINE GaussSolver (B, AP, S, N, E, W, SolOld, Sol, DSDY, NumVars)
        !This subroutine takes in the coefficients of a NumVars by NumVars matrix, the old/current solution, and the rate of change of 
        !the source term with respect to species or soot mass fraction (DSDY). The species equations are coupled at each CV (not globally), and solved using 
        !gaussian elimination. 
        !The ghost cells on each CPU are then updated
        
            !local variable declarations
            INTEGER, INTENT(IN) :: NumVars
!            DOUBLE PRECISION, DIMENSION(:,:,:) :: AP, N, S, E, W, B
!            DOUBLE PRECISION, DIMENSION(0:,:,:) :: SolOld, Sol
!            DOUBLE PRECISION, DIMENSION(:,:,:,:) :: DSDY
            DOUBLE PRECISION, DIMENSION(LM1,M1,NumVars) :: AP, N, S, E, W, B    ! RDB
            DOUBLE PRECISION, DIMENSION(0:LM1,M1,NumVars) :: SolOld, Sol        ! RDB
            DOUBLE PRECISION, DIMENSION(NumVars,NumVars,LM1,M1) :: DSDY         ! RDB
            DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: SOLTEMP,SOLTEMPSEND
            DOUBLE PRECISION :: Source
            DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: C_gauss, temp_yy_gauss, yy_gauss
            DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: B_gauss
            DOUBLE PRECISION :: BMAX, CTEMP, TIMER, SUMsp 
            INTEGER :: KS, IIMAX, KSJ, KI, KJ, KL, L
            DOUBLE PRECISION :: TOTO            !RDB
            !done local variable declarations
            
            !allocate matrices
            allocate(C_gauss(max0(KK,2*MSection)))
            allocate(temp_yy_gauss(max0(KK,2*MSection)))
            allocate(yy_gauss(max0(KK,2*MSection)))
            allocate(B_gauss(max0(KK,2*MSection),max0(KK,2*MSection)))
            allocate(SOLTEMP(NumVars,2:M2,0:L0))  
            allocate(SOLTEMPSEND(NumVars,2:M2,0:L0))           
            
            !done allocating matrices
            
            SOLTEMP = 0.0d0
            SOLTEMPSEND = 0.0d0
            B_gauss = 0.0d0  !RDB
            C_gauss = 0.0d0  !RDB

            !start solving the species equations
            DO I=2,L2
                DO J=2,M2
                
                    if(.NOT.SOLIDCVP(I,J)) then !if not in a solid region

                        DO K=1,NumVars
                            SOURCE=0.0d0
                            DO L=1,NumVars
                                ! RDB test DSDY Infinity
                                if(DSDY(K,L,I,J).GT.1.d50)then     !HUGE(0.0d0)/
                                    print*,'DSDY HUGE',myid,K,L,I,J,DSDY(K,L,I,J)!,HUGE(0.0d0)
                                    DSDY(K,L,I,J) = 1.d50      !RDB    HUGE(0.0d0)/
                                elseif(DSDY(K,L,I,J).LT.-1.d50)then        !HUGE(0.0d0)/
                                    print*,'DSDY-HUGE',myid,K,L,I,J,DSDY(K,L,I,J)!,-HUGE(0.0d0)
                                    DSDY(K,L,I,J) = -1.d50     !RDB    HUGE(0.0d0)/
                                endif
                                SOURCE=SOURCE+DSDY(K,L,I,J)*SolOld(I,J,L)
                            ENDDO
                            SOURCE=SOURCE*YCVR(J)*XCV(I)
   C_gauss(K)=N(I,J,K)*Sol(I+1,J,K)+S(I,J,K)*Sol(I-1,J,K)+E(I,J,K)*SOL(I,J+1,K)+W(I,J,K)*SOL(I,J-1,K)+ &
                            B(I,J,K)*YCVR(J)*XCV(I)-SOURCE
                            B_gauss(K,K)=AP(I,J,K)-DSDY(K,K,I,J)*YCVR(J)*XCV(I)

                            DO KJ=1,NumVars
                                IF(KJ.NE.K)then
                                    B_gauss(K,KJ)=-DSDY(K,KJ,I,J)*YCVR(J)*XCV(I)
                                ENDIF      
                            ENDDO

                        ENDDO

!                        ! test B_gauss y C_gauss
!                        do k=1,NumVars
!        if ((B_gauss(K,K).NE.B_gauss(K,K)).or.(DABS(B_gauss(K,K)).GT.HUGE(0.0d0))) then
!            print*,'testBB',myid,I,J,K,AP(I,J,K),DSDY(K,K,I,J)   !RDB
!        endif
!        if ((C_gauss(K).NE.C_gauss(K)).or.(DABS(C_gauss(K)).GT.HUGE(0.0d0))) then
!            print*,'testC',myid,I,J,K,C_gauss(K)   !RDB
!        endif
!                        do kJ=1,NumVars
!        if ((B_gauss(KJ,K).NE.B_gauss(KJ,K)).or.(DABS(B_gauss(KJ,K)).GT.HUGE(0.0d0))) then
!            print*, 'testB',myid,I,J,K,KJ,B_gauss(K,KJ) !RDB
!        endif
!                        enddo
!                        enddo

                        DO KI=1,NumVars-1
                            ! RDB loop normalizacion
                            ! normalización de cada ecuación por su mayor coeficiente
                            ! este loop permite mejorar considerablemente la estabilidad numérica, pero es lento
                            DO KS=KI,NumVars
                                BMAX=DABS(B_gauss(KS,KI))
                                DO KSJ=KI,NumVars
                                    IF(DABS(B_gauss(KS,KSJ)).GT.BMAX) then
                                        BMAX=DABS(B_gauss(KS,KSJ))
                                    ENDIF
                                ENDDO
                                IF(BMAX.NE.0.0d0)then
                                    DO KSJ=KI,NumVars   !KI,NumVars
                                        B_gauss(KS,KSJ)=B_gauss(KS,KSJ)/BMAX
                                    ENDDO
                                    C_gauss(KS)=C_gauss(KS)/BMAX
                                !ELSE
                                !    print*,'BMAX 0',KS,KI,myid
                                ENDIF
                            ENDDO
                            ! fin loop normalizacion

                            ! inicio bloque reemplazado RDB
                            IIMAX=KI
                            BMAX=DABS(B_gauss(KI,KI))
                            DO KS=KI+1,NumVars
                                IF(DABS(B_gauss(KS,KI)).GT.BMAX) THEN
                                    BMAX=DABS(B_gauss(KS,KI))
                                    IIMAX=KS
                                ENDIF
                            ENDDO
                            if(IIMAX.NE.KI)then
                                DO KSJ=1,NumVars   !KSJ=KI,NumVars
                                    temp_yy_gauss(KSJ)=B_gauss(KI,KSJ)
                                    B_gauss(KI,KSJ)=B_gauss(IIMAX,KSJ)
                                    B_gauss(IIMAX,KSJ)=temp_yy_gauss(KSJ)
                                ENDDO
                                CTEMP=C_gauss(KI)
                                C_gauss(KI)=C_gauss(IIMAX)
                                C_gauss(IIMAX)=CTEMP
                            endif
                            ! fin bloque reemplazado RDB

!                            IIMAX=KI
!                            BMAX=DABS(B_gauss(KI,KI))
!                            DO KS=KI+1,NumVars
!                                IF(DABS(B_gauss(KS,KI)).GT.BMAX) THEN
!                                    BMAX=DABS(B_gauss(KS,KI))
!                                    IIMAX=KS
!                                    DO KSJ=KI,NumVars
!                                        temp_yy_gauss(KSJ)=B_gauss(KI,KSJ)
!                                        B_gauss(KI,KSJ)=B_gauss(IIMAX,KSJ)
!                                        B_gauss(IIMAX,KSJ)=temp_yy_gauss(KSJ)
!                                    ENDDO
!                                    CTEMP=C_gauss(KI)
!                                    C_gauss(KI)=C_gauss(IIMAX)
!                                    C_gauss(IIMAX)=CTEMP
!                                ENDIF
!                            ENDDO

!                            ! test matriz singular
!                            if(BMAX.EQ.(0.0d0)) then
!                                print*,myid,I,J,KI,'matriz singular'
!!                                B_gauss(KI,KI) = smallnum       !RDB
!                            endif

                            DO KS=KI+1,NumVars
                                TIMER= B_gauss(KS,KI)/B_gauss(KI,KI)
                                DO KJ=KI,NumVars
                                    B_gauss(KS,KJ)=B_gauss(KS,KJ) - TIMER*B_gauss(KI,KJ)
                                ENDDO
                                C_gauss(KS)=C_gauss(KS) - TIMER*C_gauss(KI)
                            ENDDO

                        ENDDO


!                        ! test B_gauss2
!                        do KI=1,NumVars
  !     if ((B_gauss(KI,KI).NE.B_gauss(KI,KI)).or.(DABS(B_gauss(KI,KI)).GT.HUGE(0.0d0))) then
  !      print*,'testBB2',myid,I,J,KI,AP(I,J,KI),DSDY(KI,KI,I,J),B_gauss(KI,KI)   !RDB
  !     endif
!                        do KJ=KI,NumVars
!       if ((B_gauss(KI,KJ).NE.B_gauss(KI,KJ)).or.(DABS(B_gauss(KI,KJ)).GT.HUGE(0.0d0))) then
!        print*,'testB2',myid,I,J,KI,KJ,B_gauss(KI,KJ)   !RDB
!       endif
!                        enddo
!                        enddo
             
                        yy_gauss(NumVars)=C_gauss(NumVars)/B_gauss(NumVars,NumVars)
!    if ((yy_gauss(NumVars).NE.yy_gauss(NumVars)).or.(DABS(yy_gauss(NumVars)).GT.HUGE(0.0d0))) then         !RDB is NaN or Infinity
!        print*,'yy_gauss1',myid,NumVars,I,J,C_gauss(NumVars),B_gauss(NumVars,NumVars),yy_gauss(NumVars)    !RDB
!    endif
!                        if (yy_gauss(NumVars).LT.0.0d0) then
! !                               print*,'yy_g<0', myid,I,J,NumVars,yy_gauss(NumVars)
!                                yy_gauss(NumVars) = 0.0d0      !RDB
!                        endif
!                        if (yy_gauss(NumVars).GT.1.0d0) then
!                                print*,'yy_g>1', myid,I,J,NumVars,yy_gauss(NumVars)
!                                yy_gauss(NumVars) = 1.0d0      !RDB
!                        endif
                        DO KI=NumVars-1,1,-1
                            SUMsp=0.0D0
                            DO KJ=KI+1,NumVars
                                SUMsp=SUMsp+ B_gauss(KI,KJ)*yy_gauss(KJ)        !RDB
                            ENDDO
                            yy_gauss(KI)=(C_gauss(KI)-SUMsp)/B_gauss(KI,KI)

!     if ((yy_gauss(KI).NE.yy_gauss(KI)).or.(DABS(yy_gauss(KI)).GT.HUGE(0.0d0))) then
!        print*,'yy_gauss2',myid,KI,I,J,C_gauss(KI),SUMsp,B_gauss(KI,KI),yy_gauss(KI)    !RDB
!     endif
!                            if (yy_gauss(KI).LT.0.0d0) then
! !                               print*,'yy_g<0', myid,I,J,KI,yy_gauss(KI)
!                                yy_gauss(KI) = 0.0d0      !RDB
!                            endif
!                            if (yy_gauss(KI).GT.1.0d0) then
!                                print*,'yy_g>1', myid,I,J,KI,yy_gauss(KI)
!                                yy_gauss(KI) = 1.0d0      !RDB
!                            endif

!                            if (yy_gauss(KI).NE.yy_gauss(KI)) then
!                                print*,'prob',myid,I,J,KI,yy_gauss(KI)
!                            endif

                        ENDDO

             
                        DO K=1,NumVars
                            Sol(I,J,K)=yy_gauss(K)
                        ENDDO

                    else
                        Sol(I,J,:)=0.0d0
                    endif
                ENDDO       !J
            ENDDO       !I
            !done solving the species equation at each CV                    
            
            !update the ghost cells on each CPU
            
            !populate matrices to be sent to higher rank CPU
            DO I=L3,L2
                DO J=2,M2
                    DO K=1,NumVars
                        SOLTEMPSEND(K,J,I) = SOL(I,J,K)
                    ENDDO
                ENDDO
            ENDDO                 

            !pass ghost cell data to the higher rank CPU
            if(myid.ne.(NUMP-1)) then
                CALL MPI_SEND(SOLTEMPSEND(1,2,L3), 2*M3*NumVars,MPI_DOUBLE_PRECISION,myid+1,3770+myid,MPI_COMM_WORLD,IERR)
!                if(myid.eq.2) print*,'SEND', SOLTEMPSEND(1:10,2:5,L3:L2)
            endif 
            
            if(myid.ne.0) then
                CALL MPI_RECV(SOLTEMP(1,2,0), 2*M3*NumVars,MPI_DOUBLE_PRECISION, myid-1,3770+myid-1,MPI_COMM_WORLD,STATUS,IERR)
!                if(myid.eq.3) print*,'RECV',SOLTEMP(1:10,2:5,0:1)
            
                !place ghost cell data into SOL array
                DO I=0,1
                    DO J=2,M2
                        DO K=1,NumVars
                            SOL(I,J,K) = SOLTEMP(K,J,I)
                        ENDDO
                    ENDDO
                ENDDO    
            
            endif                
            
            call MPI_BARRIER(MPI_COMM_WORLD,IERR)
            !done passing ghost cell data to higher rank CPU
            
            !pass ghost cell data to lower rank CPU
            
            !populate matrices to be sent to lower rank CPU
            DO I=2,3
                DO J=2,M2
                    DO K=1,NumVars
                        SOLTEMPSEND(K,J,I) = SOL(I,J,K)
                    ENDDO
                ENDDO
            ENDDO
            
            if(myid.ne.0) then
                CALL MPI_SEND(SOLTEMPSEND(1,2,2), 2*M3*NumVars, MPI_DOUBLE_PRECISION,myid-1,4650+myid,MPI_COMM_WORLD,IERR)
            endif
            
            if(myid.ne.(NUMP-1)) then
                CALL MPI_RECV(SOLTEMP(1,2,L1), 2*M3*NumVars,MPI_DOUBLE_PRECISION, myid+1,4650+myid+1,MPI_COMM_WORLD,STATUS,IERR)
                       
                !place ghost cell data into SOL array
                DO I=L1,L0
                    DO J=2,M2
                        DO K=1,NumVars
                            SOL(I,J,K) = SOLTEMP(K,J,I)
                        ENDDO
                    ENDDO
                ENDDO  
            
            endif
            
            call MPI_BARRIER(MPI_COMM_WORLD,IERR)
            !done passing solution to ghost cells
          
            DEALLOCATE(C_gauss)
            DEALLOCATE(temp_yy_gauss)
            DEALLOCATE(yy_gauss)
            DEALLOCATE(B_gauss)
            DEALLOCATE(SOLTEMP)
            DEALLOCATE(SOLTEMPSEND)
            
            call MPI_BARRIER(MPI_COMM_WORLD,IERR)  
            
        end subroutine gausssolver
        !************************************************************************************************************************************************       

        !**************************************************************************************************************************************************
        SUBROUTINE CVGradients(Variable, XGrad, YGrad,Flag)
        !This subroutine takes in a variable array (such as temperature, species, soot), and returns the gradient of this variable across a CV, where the CV 
        !is defined by the current values of the global variables I and J. 
        
            !local variable declarations
            DOUBLE PRECISION, DIMENSION (0:,:) :: Variable
            DOUBLE PRECISION :: XGrad, YGrad
            DOUBLE PRECISION :: TXM, TXP, DX, TYP, TYM, DY
            INTEGER :: FLAG
            !done local variable declaration
                    
            !calculate x-direction gradient
            IF(myid .eq.0 .and. I.EQ.1) THEN
                DX=XDIF(2)
                TXM=Variable(I,J)
                TXP=Variable(I+1,J)
            ELSEIF(myid .eq. (NUMP-1) .and. I.EQ.L1) THEN
                DX=XDIF(L1)
                TXP=Variable(I,J)
                TXM=Variable(I-1,J)
            ELSEIF(myid .eq.0 .and. I.EQ.2)THEN
                DX=XCV(2) *2.5d0     !RDB
                TXP=(Variable(I,J)*XCV(I+1)+Variable(I+1,J)*XCV(I))/(2.d0*XDIF(I+1))
                TXM=Variable(I-1,J)
            ELSEIF(myid .eq.(NUMP-1) .and. I.EQ.L2)THEN
                DX=XCV(L2)
                TXP=Variable(I+1,J)
                TXM=(Variable(I,J)*XCV(I-1)+Variable(I-1,J)*XCV(I))/(2.d0*XDIF(I))
            ELSE
                DX=XCV(I)
                TXP=(Variable(I,J)*XCV(I+1)+Variable(I+1,J)*XCV(I))/(2.d0*XDIF(I+1))
                TXM=(Variable(I,J)*XCV(I-1)+Variable(I-1,J)*XCV(I))/(2.d0*XDIF(I))
            ENDIF
            
            IF(FLAG.EQ.1.AND.SOLIDCVP(I-1,J)) TXM = Variable(I,J)
            IF(FLAG.EQ.1.AND.SOLIDCVP(I+1,J)) TXP = Variable(I,J)
            
            XGrad = (TXP-TXM)/DX
           
            !done calculating x-direction gradient
            
            !calculate y-direction gradient
            IF(J.EQ.1) THEN
                DY=YDIF(2)
                TYM=Variable(I,J)
                TYP=Variable(I,J+1)
            ELSEIF(J.EQ.M1) THEN
                DY=YDIF(M1)
                TYP=Variable(I,J)
                TYM=Variable(I,J-1)
            ELSEIF(J.EQ.2) THEN
                DY=YCV(J)
                TYP=(Variable(I,J)*YCV(J+1)+Variable(I,J+1)*YCV(J))/(2.d0*YDIF(J+1))
                TYM=Variable(I,J-1)
            ELSEIF(J.EQ.M2)THEN
                DY=YCV(J)
                TYP=Variable(I,J+1)
                TYM=(Variable(I,J)*YCV(J-1)+Variable(I,J-1)*YCV(J))/(2.d0*YDIF(J))
            ELSE
                DY=YCV(J)
                TYP=(Variable(I,J)*YCV(J+1)+Variable(I,J+1)*YCV(J))/(2.d0*YDIF(J+1))
                TYM=(Variable(I,J)*YCV(J-1)+Variable(I,J-1)*YCV(J))/(2.d0*YDIF(J))
            ENDIF
            
            IF(J.NE.1) THEN
                IF(J.NE.M1) THEN
                    IF(FLAG.EQ.1.AND.SOLIDCVP(I,J-1)) TYM = Variable(I,J)
                    IF(FLAG.EQ.1.AND.SOLIDCVP(I,J+1)) TYP = Variable(I,J)
                ENDIF     
            ENDIF     
            
            YGrad = (TYP - TYM)/DY  
               
            !done calculating y-direction gradient

        end subroutine cvgradients
        !***********************************************************************************************************

        !**************************************************************************************************************************************************
        SUBROUTINE SpeciesDiff
        !This subroutine calculates the diffusion velocities for each species for each CV on a given CPU
        
            !local variable declarations
            DOUBLE PRECISION :: DYKDX, DYKDY, DENOM, DEFF
            DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: SMoleF, DIFFU, TDR
            !done local variable declarations
            
            !allocate matrices
            allocate(SMoleF(KK))
            allocate(DIFFU(KK))
            allocate(TDR(KK))
            !done allocating matrices 
            
            !loop over all CVs 
            DO J=2,M2
                DO I=1,L1 
                
                    if(.NOT.SOLIDCVP(I,J)) then !if not in a solid region
                    
                        CALL CKYTX(SpeciesMF(I,J,:),IWORK,WORK,SMoleF) !get species mole fractions
                        CALL MCADIF(PRESSURE*PCKIN,T(I,J),SMoleF,TPWRK,DIFFU)   !get species diffusion coefficients
                        CALL MCATDR(T(I,J),SMoleF,ITPWRK,TPWRK,TDR)          !get species thermal diffusion ratios
                        
                        DO K=1,KK  !loop over all species, and get the normal and thermal diffusional velocities
                            !IF(DENOM.eq.0.d0) DENOM = smallnum
                            DENOM = 1.d0 - SpeciesMF(I,J,K) + smallnum  !RDB
                            DEFF = DIFFU(K)*(1.d0-SmoleF(K))/DENOM
                            CALL CVGradients(SpeciesMF(:,:,K),DYKDX,DYKDY,1) !get the gradient of the current species across the current CV
                            VKX(I,J,K)=  -DEFF*DYKDX/(SpeciesMF(I,J,K)+smallnum)
                            VKY(I,J,K)=  -DEFF*DYKDY/(SpeciesMF(I,J,K)+smallnum)
                            VTX(I,J,K)= DIFFU(K)*TDR(K)*DTDX(I,J)/(SMoleF(K)*T(I,J)+smallnum)
                            VTY(I,J,K)= DIFFU(K)*TDR(K)*DTDY(I,J)/(SMoleF(K)*T(I,J)+smallnum)
                        ENDDO
                        !done getting diffusion velocities

                    else

                        VKX(I,J,:) = 0.d0
                        VKY(I,J,:) = 0.d0
                        VTX(I,J,:) = 0.d0
                        VTY(I,J,:) = 0.d0
                    endif
                ENDDO
            ENDDO

            !deal with special cases
            if(myid .eq.0) then
                VTX(1,:,:)=0.d0
                VTY(1,:,:)=0.d0
                VKX(1,:,:)=0.d0
                VKY(1,:,:)=0.d0
            endif

            if(myid .eq.(NUMP-1)) then
                VTX(L1,:,:)=0.d0
                VTY(L1,:,:)=0.d0
                VKX(L1,:,:)=0.d0
                VKY(L1,:,:)=0.d0
            endif
            
            VTX(:,1,:) = 0.d0
            VTY(:,1,:) = 0.d0
            VKX(:,1,:) = 0.d0
            VKY(:,1,:) = 0.d0
            
            VTX(:,M1,:) = 0.d0
            VTY(:,M1,:) = 0.d0
            VKX(:,M1,:) = 0.d0
            VKY(:,M1,:) = 0.d0
                    
            
            DEALLOCATE(SMoleF)
            DEALLOCATE(DIFFU)
            DEALLOCATE(TDR)
        
        end subroutine SpeciesDiff
        !*************************************************************************************************************************************
        
        !**************************************************************************************************************************************************
        SUBROUTINE SootDiff
        !This subroutine calculates the diffusion velocities for each soot section for each CV on a given CPU
        
            !local variable declarations
            DOUBLE PRECISION :: DADX, DADY, DPDX, DPDY, Pdiffusivity, DEFF, SootMF, VIS
            DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: SMoleF
            !done local variable declarations
            
            !allocate arrays
            allocate(SMoleF(KK))
            !done allocating arrays
        
            !loop over all CVs
            DO J=2,M2
                DO I=1,L1
                    
                    if(.NOT.SOLIDCVP(I,J)) then !if not in a solid region
                    
                        CALL CKYTX(SpeciesMF(I,J,:),IWORK,WORK,SMoleF)
                        CALL MCAVIS(T(I,J),SMoleF,TPWRK,VIS)
                    
                        DO K=1,MSection !loop over the sections
                            
                            !first calculate ordinary soot diffusion, which can only be accounted for with sectional model
                            if(SootModel.Eq.1) then
                                CALL Pdiff (k,Pdiffusivity)
                                DEFF = Pdiffusivity
                                CALL CVGradients(SootSec(:,:,K),DADX,DADY,1)
                                CALL CVGradients(SootSec(:,:,MSection+K),DPDX,DPDY,1)
                                SootMF = DEXP(XS(K))*SootSec(I,J,K)
                                VSX(I,J,K)=-DEFF*(DADX*DEXP(XS(K)))/(SootMF+smallnum)
                                VSY(I,J,K)=-DEFF*(DADY*DEXP(XS(K)))/(SootMF+smallnum)
                            endif

                        ENDDO

                       !calculate thermophoresis diffusion velocities
                        VSTX(I,J,:) = -0.55d0*VIS/RHO(I,J)*DTDX(I,J)/T(I,J) 
                        VSTY(I,J,:) = -0.55d0*VIS/RHO(I,J)*DTDY(I,J)/T(I,J)
                        
                        !if near a solid region, set thermophoretic velocity to zero
                        IF(J.GT.2.AND.I.GT.1.AND.J.LT.M2) THEN
                            if(SOLIDCVP(I-2,J).OR.SOLIDCVP(I,J+2).OR.SOLIDCVP(I,J-2).OR.SOLIDCVP(I+2,J)) then
                                VSTX(I,J,:) = 0.d0 
                                VSTY(I,J,:) = 0.d0
                            endif
                        ENDIF     
                        
                        IF(J.GT.1) THEN
                            if(SOLIDCVP(I-1,J).OR.SOLIDCVP(I,J+1).OR.SOLIDCVP(I,J-1).OR.SOLIDCVP(I+1,J)) then
                                VSTX(I,J,:) = 0.d0 
                                VSTY(I,J,:) = 0.d0
                            endif    
                        ENDIF           

                        !done calculating soot diffusion velocities at a given CV

                    else
                        VSTX(I,J,:) = 0.d0
                        VSTX(I,J,:) = 0.d0 
                        VSX(I,J,:) = 0.d0 
                        VSY(I,J,:) = 0.d0
                    endif               
                ENDDO
            ENDDO
            !done calculating diffusion at all CVs

            !deal with special cases
            if(myid .eq.0) then
                VSTX(1,:,:)=0.d0
                VSTY(1,:,:)=0.d0
                VSX(1,:,:)=0.d0
                VSY(1,:,:)=0.d0
            endif

            if(myid .eq.(NUMP-1)) then
                VSTX(L1,:,:)=0.d0
                VSTY(L1,:,:)=0.d0
                VSX(L1,:,:)=0.d0
                VSY(L1,:,:)=0.d0
            endif
            
            VSTX(:,1,:) = 0.d0
            VSTY(:,1,:) = 0.d0
            VSX(:,1,:) = 0.d0
            VSY(:,1,:) = 0.d0
            
            VSTX(:,M1,:) = 0.d0
            VSTY(:,M1,:) = 0.d0
            VSX(:,M1,:) = 0.d0
            VSY(:,M1,:) = 0.d0 
        
        end subroutine SootDiff
        !**************************************************************************************************************
        
        !**************************************************************************************************************************************************
        SUBROUTINE CorDiffVel
        !This subroutine calculates the correction diffusion velocities for each CV on each CPU
        
            !initialize correction velocities to zero
            VCX = 0.0d0
            VCY = 0.0d0
            
            DO J=2, M2
                DO I=1,L1
                    
                    DO K=1, KK   !calculate the contribution to correction diffusion velocity due to species
                        VCX(I,J) = VCX(I,J) - SpeciesMF(I,J,K)*(VKX(I,J,K)+VTX(I,J,K))
                        VCY(I,J) = VCY(I,J) - SpeciesMF(I,J,K)*(VKY(I,J,K)+VTY(I,J,K))
                    ENDDO
                    
                    !calculate the contribution to correction diffusion velocity due to soot
                    if(SootModel.Eq.1) then
                        DO K=1,Msection 
                            VCX(I,J) = VCX(I,J) - DEXP(XS(K))*SootSec(I,J,K)*(VSX(I,J,K)+VSTX(I,J,K))
                            VCY(I,J) = VCY(I,J) - DEXP(XS(K))*SootSec(I,J,K)*(VSY(I,J,K)+VSTY(I,J,K))
                        ENDDO 
                    else
                        VCX(I,J) = VCX(I,J) - SootSec(I,J,1)*(VSTX(I,J,K))
                        VCY(I,J) = VCY(I,J) - SootSec(I,J,1)*(VSTY(I,J,K)) 
                    endif
                    
                ENDDO
            ENDDO
            !done calculating correction diffusion velocities
        
        end subroutine CorDiffVel
        !***************************************************************************************************************************************************  
        
        !*****************************************************************************************************************************************************
        SUBROUTINE SurfaceRates
        !This subroutine calculates the reaction rates of soot with gaseous species for a given CV defined by the global I and J variables, and stores the results 
        !in SootGasRates. When using the sectional code the rates will be determined on a per section basis. Rates are stored as mol/cc/s
        
            !local variable declarations
            DOUBLE PRECISION, PARAMETER :: GasCon = 1.987D-3, CHI = 2.3D+15 
            DOUBLE PRECISION :: RT, ssRatio, PNX, denom
            DOUBLE PRECISION :: fr1, rr1, fr2, rr2, fr3, fr4, fr5, fr6
            DOUBLE PRECISION, DIMENSION(6) :: SR
            DOUBLE PRECISION :: Area, BETACOLL
            DOUBLE PRECISION :: Mass_OH, Mass_BAPYR, Mass_BGHIF, Mass_BAPYRS, Mass_C6H6, Mass_A4, DUMMY
            DOUBLE PRECISION :: H, H2, OH, H2O, C2H2, O2, BAPYR, BGHIF, BAPYRS, C6H6, A4
            DOUBLE PRECISION, PARAMETER :: CondEff = 1.0d0 !condensation efficiency 0.5d0 F. Liu original code (S. Dworkin 2001)
            DOUBLE PRECISION :: SOOTSECSUM                 !Condensation efficiency 1.0D0 N.A Eaves CPC (2016).
            DOUBLE PRECISION :: d_soot, a_soot
            INTEGER :: G
            !end local variable declarations

            !initialize reaction rates
            SR = 0.0d0
            SootGasRates(I,J,:,:) = 0.0d0
            Dummy = 10.d0
            
            if(SootModel.Eq.1) then !if using the sectional code
            
                SOOTSECSUM = sum(SootSec(I,J,:)) !determine if any soot is present in this CV
                
                if(SOOTSECSUM.GT.0.d0) then !if there is soot present
            
                    !============================================================================
                    !   From Berkeley Soot code, these are the reactions considered for the HACA mechanism
                    !   Surface growth is currently described by the following scheme:
                    !     1.  Csoot-H + H = Csoot* + H2         (fR1, rR1)
                    !     2.  Csoot-H + OH = Csoot* + H2O   (fR2, rR2)
                    !     3.  Csoot*  + H -> Csoot              (fR3)
                    !     4.  Csoot*  + C2H2 -> Csoot-H + H     (fR4)
                    !     5.  Csoot*  + O2 -> products (2CO)    (fR5)
                    !     6.  Csoot   + OH -> products (CO)          (fR6)
                    !*****************************************************************************
                    !   The following rate coefficients are from Appel et al. (2000)
                    
                    !first, get species molar concentrations
                    H = RHO(I,J)*SpeciesMF(I,J,IH)/WT(IH)
                    SootMolConc(I,J,1) = H
                    H2 = RHO(I,J)*SpeciesMF(I,J,IH2)/WT(IH2)
                    SootMolConc(I,J,2) = H2
                    OH = RHO(I,J)*SpeciesMF(I,J,IOH)/WT(IOH)
                    SootMolConc(I,J,3) = OH
                    H2O = RHO(I,J)*SpeciesMF(I,J,IH2O)/WT(IH2O)
                    SootMolConc(I,J,4) = H2O
                    C2H2 = RHO(I,J)*SpeciesMF(I,J,IC2H2)/WT(IC2H2)
                    SootMolConc(I,J,5) = C2H2
                    O2 = RHO(I,J)*SpeciesMF(I,J,IO2)/WT(IO2)
                    SootMolConc(I,J,6) = O2
                    BAPYR = RHO(I,J)*SpeciesMF(I,J,IBAPYR)/WT(IBAPYR)
                    SootMolConc(I,J,7) = BAPYR
                    BGHIF = RHO(I,J)*SpeciesMF(I,J,IBGHIF)/WT(IBGHIF)
                    SootMolConc(I,J,8) = BGHIF
                    BAPYRS = RHO(I,J)*SpeciesMF(I,J,IBAPYRS)/WT(IBAPYRS)
                    SootMolConc(I,J,9) = BAPYRS

                    A4 = RHO(I,J)*SpeciesMF(I,J,IA4)/WT(IA4)                !for PYRENE condensation
                    SootMolConc(I,J,10) = A4
                    !done getting species molar concentrations
                    
                    !get mass of OH, BAPYR, BGHIF, and BAPYR
                    Mass_OH = WT(IOH)*AMU
                    Mass_BAPYR = WT(IBAPYR)*AMU
                    Mass_BAPYRS = WT(IBAPYRS)*AMU
                    Mass_BGHIF = WT(IBGHIF)*AMU

                    Mass_A4 =  WT(IA4)*AMU           !for Pyrene condensation

                    !done getting masses  
                    
                    !calculate part of the reaction rates, neglecting soot surface surface area density             
                    RT = GasCon*T(I,J)
                    
                    fR1 = 4.2D+13 * DEXP(-13.0D0 / RT) * H
                    SootHacaRates(I,J,1) = fR1

                    rR1 = 3.9D+12  * DEXP(-11.0D0 / RT) * H2
                    SootHacaRates(I,J,2) = rR1

                    fR2 = 1.0D+10 * (T(I,J)**0.734D0)* DEXP(-1.43D0 / RT) * OH
                    SootHacaRates(I,J,3) = fR2

                    rR2 = 3.68D+08 * (T(I,J)**1.139D0) * DEXP(-17.1D0 / RT) * H2O
                    SootHacaRates(I,J,4) = rR2

                    fR3 = 2.0D+13   * H
                    SootHacaRates(I,J,5) = fR3

                    fR4 = 8.0D+7   * (T(I,J)**1.56D0) * DEXP(-3.8D0 / RT) * C2H2
                    SootHacaRates(I,J,6) = fR4

                    fR5 = 2.2D+12 * DEXP(-7.5D0 / RT) * O2
                    SootHacaRates(I,J,7) = fR5

                    fR6 = 0.13D0 * OH        ! gamma = 0.13 from Neoh et al.
                    
                    !if alpha_surf is equal to 5, use the functional form. If it is not 5, just use the value
                    !from the input.dat file   !FIX THIS!!
                    if(Alpha_Surf.EQ.5) CALL Alpha_function
                    
                    !determine the ratio of abstracted sites to non-abstracted sites (ssRatio)
                    denom = rR1 + rR2 + fR3 + fR4 + fR5 + fR1 + fR2
                    !if denominator is not zero, reactions will occur
                    if(denom.NE.0.d0) then
                        
                        ssRatio = (fR1 + fR2) / denom !determine the ssRatio
                        SootHacaRates(I,J,8) = ssRatio

                        !calculate the net rate for each reaction, considering reverse rates, but still neglecting soot surface area density
                        SR(1) = (fR1*(1.d0-ssratio) - rR1 * ssRatio) * CHI
                        SootHacaRates(I,J,9) = SR(1)
                        
                        SR(2) = (fR2*(1.d0-ssratio) - rR2 * ssRatio) * CHI
                        SootHacaRates(I,J,10) = SR(2)
                        
                        SR(3) = fR3 * ssRatio * CHI 
                        SootHacaRates(I,J,11) = SR(3)
                        
                        SR(4) = fR4 * ssRatio * CHI
                        SootHacaRates(I,J,12) = SR(4)
                        
                        SR(5) = fR5 * ssRatio * CHI  
                        SootHacaRates(I,J,13) = SR(5)
                        
                        SR(6) = fR6
                        
                        !determiine soot surface area density, and the resulting reaction rates. Store in SootGasRates
                        DO K=1,MSection
                        
                            !determine number of primary particles per aggregate
                            IF(SootSec(I,J,K).LE.smallnum.OR.SootSec(I,J,K+MSection).LE.smallnum) THEN !if number density of aggregates or primaries is low, assume 1 primary per aggregate
                                PNX=1.D0
                            ELSE
                                PNX=SootSec(I,J,K+MSection)/SootSec(I,J,K) !else, just Np/Na
                            ENDIF
                            
                            Area=PNX*4.D0*PI*(3.D0*(DEXP(XS(K))/PNX)/4.D0/PI/densityP)**(2.d0/3.d0) !calculate soot surface area per aggregate(cm^2/aggregate)
                            
                            DO G=1,5  !special case for reaction 6                          
                                SootGasRates(I,J,K,G) = Alpha_surf*Area*(SootSec(I,J,K)*RHO(I,J))*SR(G)/AV 
                            ENDDO
                            
                            if(SootSec(I,J,K).GT.0d0) then !if there is soot present in this section, calculation OH oxidation/condensation rates
                            
                                CALL BETA_FREE(Rad_OH, Mass_OH,0, DUMMY, PNX,K,BETACOLL)
                                SootGasRates(I,J,K,6) = BETACOLL*SR(6)*(SootSec(I,J,K)*RHO(I,J))
                                
                                if(Mechanism.LE.3) then
                                    CALL BETA_AEROSOL(Rad_BAPYR, Mass_BAPYR,0, DUMMY, PNX,K,BETACOLL)
                                    SootGasRates(I,J,K,7) = CondEff*BETACOLL*BAPYR*(SootSec(I,J,K)*RHO(I,J))

                                    CALL BETA_AEROSOL(Rad_BGHIF, Mass_BGHIF,0, DUMMY, PNX,K,BETACOLL)
                                    SootGasRates(I,J,K,8) = CondEff*BETACOLL*BGHIF*(SootSec(I,J,K)*RHO(I,J))

                                    CALL BETA_AEROSOL(Rad_BAPYRS, Mass_BAPYRS,0, DUMMY, PNX,K,BETACOLL)
                                    SootGasRates(I,J,K,9) = CondEff*BETACOLL*BAPYRS*(SootSec(I,J,K)*RHO(I,J))
                                else    !Pyrene condensation
                                    !CALL BETA_AEROSOL(Rad_C6H6, Mass_C6H6,0, DUMMY, PNX,K,BETACOLL)
                                    !SootGasRates(I,J,K,7) = CondEff*BETACOLL*C6H6*(SootSec(I,J,K)*RHO(I,J))
                                    CALL BETA_AEROSOL(Rad_A4, Mass_A4,0, DUMMY, PNX,K,BETACOLL)
                                    SootGasRates(I,J,K,7) = CondEff*BETACOLL*A4*(SootSec(I,J,K)*RHO(I,J))
                                    SootGasRates(I,J,K,8) = 0.0d0
                                    SootGasRates(I,J,K,9) = 0.0d0
                                end if
                            else
                                SootGasRates(I,J,K,6)=0.d0
                                SootGasRates(I,J,K,7)=0.d0
                                SootGasRates(I,J,K,8)=0.d0
                                SootGasRates(I,J,K,9)=0.d0
                            endif           
                        
                        ENDDO      
                                        
                    else !if denom = 0    
                        SootGasRates(I,J,:,:) = 0.0d0
                    endif         
                
                else !if no soot is present
                    SootGasRates(I,J,:,:) = 0.0d0
                endif
            
            else !if using 2-eq code                 
                    
                !get species molar concentrations
                OH = RHO(I,J)*SpeciesMF(I,J,IOH)/WT(IOH)
                C2H2 = RHO(I,J)*SpeciesMF(I,J,IC2H2)/WT(IC2H2)
                O2 = RHO(I,J)*SpeciesMF(I,J,IO2)/WT(IO2)
                                
                !calculating diameter of soot particle in [cm]. See Woolley et. al 2009 for reference
                d_soot =  ((6.0d0*SootSec(I,J,1))/(PI*densityP*SootSec(I,J,2)+smallnum))**(1.0d0/3.0d0)
                
                !calculate surface area per unit volume in [cm^2/cm^3]. See Woolley et. al 2009 for reference
                a_soot = PI*d_soot**(2.0d0)*RHO(I,J)*(SootSec(I,J,2)+smallnum)
                
                !Calculate rate of surface growth in [mol/cc/sec]. See Woolley et. al 2009 for reference
                !C2H2 + nC --> (n+2) C +H2
                !R_g = 2.0d0*750.0d0*100.0d0*exp(-12100.0d0/T)*C2H2*a_soot !multiplied by 100 to convert k rate from [m/s] to [cm/s]
                SootGasRates(I,J,1,1) = 0.7d0*50000.d0*dexp(-24000.d0/1.987d0/T(I,J))*a_soot*c2h2  !reaction rate of growth in mol/cc/sec !rate is from Fairweather 1992 model !ADJ GROWTH
                
                !calculate rate of oxidation for O2 in [mol/cc/s]. See Woolley et. al 2009 for reference
                !0.5O2 + C --> CO
                !R_o2 = 715.0d0*100.0d0*a_soot*T**0.50d0*exp(-19680.0d0/T)*O2 !multiplied by 100 to convert k rate from [m/s] to [cm/s]
                SootGasRates(I,J,1,2) = 1.78d4*100.d0*dsqrt(T(I,J))*dexp(-39000.d0/1.987d0/T(I,J))*a_soot*o2  !reaction rate of O2 oxidation in mol/cc/sec !rate is from Fairweather 1992 model

                !calculate rate of oxidation for OH in [mol/cc/s]. See Woolley et. al 2009 for reference
                !OH + C --> CO + H
                !R_oh = 0.36*100.0d0*a_soot*T**0.50d0*OH !multiplied by 100 to convert k rate from [m/s] to [cm/s]
                SootGasRates(I,J,1,3) = a_soot*oh*4.8d0*0.36d0*100.0d0*T(I,J)**0.50d0  !multiplied by 100 to convert k rate from [m/s] to [cm/s] !Woolley et. al 2009 rate !Multipied by 4.8 for Adjustment (JUSTIN)
            
            endif
                   
        end subroutine SurfaceRates
        !***********************************************************************************************************************************              

        !*********************************************************************************************************************************************
        SUBROUTINE NuclChem
        !This subroutine calculates the nucleation rates for a given CV, defined by the current I and J global variables, and stores the results
        !in SootNucRates. Nucleation rates are stored as mol/cc/s on a per nucleating species basis in the first half of the array, and in terms of
        !g/cc/s in the second half of the array

            !local variable declarations
            DOUBLE PRECISION :: temprate, BETACOLL, temprate1
            DOUBLE PRECISION, DIMENSION(NumNuc) :: PAHC, PAH_rad, PAHMolCon
            INTEGER :: L, M
            DOUBLE PRECISION :: NucEff, Keq, Kfor, Krev
            DOUBLE PRECISION :: C2H2, C6H6
            DOUBLE PRECISION, PARAMETER :: Cut = 0.d0
            !done local variable declarations

            !set nucleation rates to zero
            SootNucRates(I,J,:) = 0.d0

            if(SootModel.EQ.1) then

                !populate PAHC, PAH_rad, and PAHMolCon
                if(Mechanism.LE.3) then
                    !BAPYR
                    PAHC(1)= 20
                    PAH_rad(1) = RAD_BAPYRS
                    PAHMolCon(1) = RHO(I,J)*SpeciesMF(I,J,IBAPYR)/WT(IBAPYR)

                    !BGHIF
                    PAHC(2)= 18
                    PAH_rad(2) = RAD_BGHIF
                    PAHMolCon(2) = RHO(I,J)*SpeciesMF(I,J,IBGHIF)/WT(IBGHIF)

                    !BAPYRS
                    PAHC(3)= 20
                    PAH_rad(3) = RAD_BAPYR
                    PAHMolCon(3) = RHO(I,J)*SpeciesMF(I,J,IBAPYRS)/WT(IBAPYRS)

                else    !Nucleation based in Pyrene (A4) for non-5-rings PAH chemical kinetics (ABF and MARINOV)
                    PAHC(1) = 16
                    PAH_rad(1) = RAD_A4
                    PAHMolCon(1) = RHO(I,J)*SpeciesMF(I,J,IA4)/WT(IA4)
                end if
                NucEff = 0.0001d0 !* 2.d0    !original: 0.0001d0
                !Keq = 3.0D-15*exp(16148*(1.d0/T(I,J)))*(82.05d0*T(I,J))


                do l=1, NumNuc
                    do m = l ,NumNuc
                        if(Pressure.LT.10d0) then  !if pressure is sufficiently low, use free molecular regime collision theory
                            CALL BETA_FREE(PAH_rad(l), PAHC(l)*C_mass,0,PAH_rad(m), PAHC(m)*C_mass,0,BETACOLL)
                            Kfor = NucEff*BETACOLL*AV**2
                            temprate =  Kfor*PAHMolCon(l)*PAHMolCon(m)

                        else !use transition/continuum regime collision theory

                            CALL BETA_AEROSOL(PAH_rad(l), PAHC(l)*C_mass,0,PAH_rad(m), PAHC(m)*C_mass,0,BETACOLL)
                            Kfor = NucEff*BETACOLL*AV**2
                            temprate =  Kfor*PAHMolCon(l)*PAHMolCon(m)

                        endif

                        SootNucRates(I,J,L) = SootNucRates(I,J,L) + temprate/AV
                        SootNucRates(I,J,M) = SootNucRates(I,J,M) + temprate/AV
                        SootNucRates(I,J,NumNuc+L) = SootNucRates(I,J,NumNuc+L)+ temprate*PAHC(l)*C_mass
                        SootNucRates(I,J,NumNuc+M) = SootNucRates(I,J,NumNuc+M)+ temprate*PAHC(m)*C_mass
                    enddo
                enddo

            else
                if(SootModel.Eq.2) then !if using 2-eq model with only acetylene nucleation
                    C2H2 = RHO(I,J)*SpeciesMF(I,J,IC2H2)/WT(IC2H2)
                    !C2H2-->2C +H2
                    SootNucRates(I,J,1) = 1.35d6*DEXP(-41000.d0/1.987d0/T(I,J))*C2H2  !rate is from Fairweather 1992 model
                
                else !if using the 2-eq model with acetylene and benzene inception
                    C2H2 = RHO(I,J)*SpeciesMF(I,J,IC2H2)/WT(IC2H2)
                    C6H6 = RHO(I,J)*SpeciesMF(I,J,IC6H6)/WT(IC6H6)
                    
                    !calculate rate of soot nucleation based on C2H2 in [mol/cc/sec]. See Woolley et. al 2009 for reference
                    SootNucRates(I,J,1) = 10000.0d0*dexp(-21000.0d0/T(I,J))*C2H2

                    !calculate rate of soot nucleation based on C6H6 in [mol/cc/sec]. See Woolley et. al 2009 for reference
                    SootNucRates(I,J,2) = 75000.0d0*dexp(-21000.0d0/T(I,J))*C6H6
                endif
            endif 
        
        end subroutine NuclChem                       
        !****************************************************************************************************************************************************** 
        
        !***************************************************************************************************************************************************************
        SUBROUTINE GasScrubbing(Scrub)
        !This subroutine takes SootSurfaceRates and SootNucRates, and determines the scrubbing rates for each of the KK species in the mechanism, for a given CV defined
        !by the global I and J variables. Coupling between soot and species can be altered in this subroutine.
        
            !local variable declarations
            DOUBLE PRECISION, DIMENSION(:) :: Scrub
            !done local variable declarations
            
            !initialize scrub
            Scrub = 0.0d0
            
            if(SootModel.Eq.1) then !if using the sectional model
            
                !============================================================================
                !   From Berkeley Soot code, these are the reactions considered for the HACA mechanism
                !   Surface growth is currently described by the following scheme:
                !     1.  Csoot-H + H = Csoot* + H2         (fR1, rR1)
                !     2.  Csoot-H + OH = Csoot* + H2O   (fR2, rR2)
                !     3.  Csoot*  + H -> Csoot              (fR3)
                !     4.  Csoot*  + C2H2 -> Csoot-H + H     (fR4)
                !     5.  Csoot*  + O2 -> products (2CO)    (fR5)
                !     6.  Csoot   + OH -> products (CO)          (fR6)
                !*****************************************************************************
                
                Scrub(IH) = -sum(SootGasRates(I,J,:,1))-sum(SootGasRates(I,J,:,3))+sum(SootGasRates(I,J,:,4))
                Scrub(IH2) = sum(SootGasRates(I,J,:,1))
                Scrub(IOH) = -sum(SootGasRates(I,J,:,2))-sum(SootGasRates(I,J,:,6))
                Scrub(IH2O) =  sum(SootGasRates(I,J,:,2))
                Scrub(IO2) =  -sum(SootGasRates(I,J,:,5))
                Scrub(ICO) =  2.d0*sum(SootGasRates(I,J,:,5)) + sum(SootGasRates(I,J,:,6))            
                Scrub(IC2H2) = -sum(SootGasRates(I,J,:,4))
                
                if(Mechanism.LE.3) then
                    Scrub(IBAPYR) = -sum(SootGasRates(I,J,:,7)) - SootNucRates(I,J,1)
                    Scrub(IBGHIF) = -sum(SootGasRates(I,J,:,8)) - SootNucRates(I,J,2)
                    Scrub(IBAPYRS) = -sum(SootGasRates(I,J,:,9)) - SootNucRates(I,J,3)
                else    !pyrene condensation
                    Scrub(IA4) = -sum(SootGasRates(I,J,:,7)) - SootNucRates(I,J,1)
                end if
            
            else !if using two equation model
            
                Scrub(IC2H2) = -SootGasRates(I,J,1,1) - SootNucRates(I,J,1)      
                Scrub(IH) = SootGasRates(I,J,1,3)
                Scrub(IH2) = SootNucRates(I,J,1)+3.d0*SootNucRates(I,J,2) + SootGasRates(I,J,1,1)
                Scrub(IO2) = -0.5d0*SootGasRates(I,J,1,2)
                Scrub(IOH) = -SootGasRates(I,J,1,3) 
                Scrub(ICO) = SootGasRates(I,J,1,2)+SootGasRates(I,J,1,3) 
                IF(sootmodel.eq.3) Scrub(IC6H6) = -SootNucRates(I,J,2)
            
            endif
        
        end subroutine GasScrubbing
        !**************************************************************************************************************************************************************
        
        !***************************************************************************************************************************************************************                                    
        SUBROUTINE Soot2eqRate(SOURCETERM)                        
        !This subroutine takes the nucleation and surface growth rates and determines the source terms for soot mass fraction and particles for the 2-eq model for a given CV,
        !defined by the current I and J global variables. These values are stored in SOURCETERM

            !local variable declarations
            DOUBLE PRECISION, DIMENSION(:) :: SOURCETERM 
            DOUBLE PRECISION :: d_soot
            !done local variable declarations
            
            !source term for soot mass fraction
            SOURCETERM(1) = 12.d0*(2.d0*SootGasRates(I,J,1,1)-SootGasRates(I,J,1,2)-SootGasRates(I,J,1,3) &
                            +SootNucRates(I,J,1)*2.d0+SootNucRates(I,J,2)*6.d0)
            
            !source term for soot number density
            d_soot =  ((6.0d0*SootSec(I,J,1))/(PI*densityP*SootSec(I,J,2)+smallnum))**(1.0d0/3.0d0)
            SOURCETERM(2) = ((2.0d0*SootNucRates(I,J,1)+6.0d0*SootNucRates(I,J,2))*(AV/C_min_mono) - & 
                             2.0d0*C_aggl*(d_soot)**(1.0d0/2.0d0)*(6.0d0*BOLTZMANN*T(I,J)/DensityP)**(1.0d0/2.0d0) &
                             *(RHO(I,J)*SootSec(I,J,2))**(2.0d0))
        
        end subroutine Soot2eqRate
        !************************************************************************************************************************************************************************
        
        !*********************************************************************************************************************************************************************
        SUBROUTINE SootSecRate(SOURCETERM,DSDYFLAG,CS1,CNS1,CS2,CNS2)                        
        !This subroutine takes the nucleation and surface HACA and condensation growth rates and determines the source terms for soot mass fraction and particles for the 2-eq model
        !for a given CV, defined by the current I and J global variables. These values are stored in SOURCETERM. All rates in #/cc/s

            !local variable declarations
            DOUBLE PRECISION, DIMENSION(:) :: SOURCETERM
            DOUBLE PRECISION, PARAMETER :: CoagEff = 0.8d0 !0.5d0 !MOD iter 57k cond_gamma f. Liu
            DOUBLE PRECISION, DIMENSION(:,:,:) :: CS1, CNS1
            DOUBLE PRECISION, DIMENSION(:,:) :: CS2, CNS2
            DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CollRate
            DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: NucRate, FragRate, SurfRate, SintRate, ObliRate, PNX
            DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: S_frag, PNXOLD, CoagRate
            DOUBLE PRECISION :: TotalNucRate, BetaCoag
            DOUBLE PRECISION :: AINPAH, AINC2H2, AINO2, AINOH, AOUTPAH, AOUTC2H2, AOUTO2, AOUTOH, temp
            DOUBLE PRECISION :: total_rate_for_frag, specific_ox_rate, Area
            DOUBLE PRECISION :: SOOTSECSUM, tao_i_i, tao_i_ip1
            INTEGER :: L, R, S
            DOUBLE PRECISION :: W, WL, WU, UTOT, DUMMY
            DOUBLE PRECISION :: A_frag
            LOGICAL :: DSDYFLAG
            !done local variable declarations
            
            !allocate Matrices
            allocate(NucRate(2*MSection))
            allocate(PNX(MSection))
            allocate(PNXOLD(MSection)) 
            allocate(S_frag(MSection))
            allocate(FragRate(2*MSection))
            allocate(SurfRate(2*MSection))
            allocate(SintRate(2*MSection))
            allocate(CoagRate(2*MSection))
            allocate(ObliRate(2*MSection))
            allocate(CollRate(MSection,MSection))
           
            !initialize all matrices
            NucRate = 0.d0
            SurfRate = 0.d0
            FragRate = 0.d0
            SintRate = 0.d0
            ObliRate = 0.d0
            CoagRate = 0.d0
            Dummy = 10.d0
            
            !determine if any soot is present
            SOOTSECSUM = sum(SootSec(I,J,:))
            
            !determine the number of primary particles per aggregate, and store in PNX
            DO L=1, MSection      
                IF(SootSec(I,J,L).LE.smallnum .OR. SootSec(I,J,L+MSection).LE.smallnum) THEN
                    PNX(L)=1.d0
                ELSE
                    PNX(L)=SootSec(I,J,L+MSection)/(SootSec(I,J,L))
                ENDIF
                IF(PNX(L).LT.1.d0) THEN
                    PNX(L) = 1.d0
                ENDIF  
            ENDDO
            !done determining the number of primary particles per aggregate.    
            
            IF(DSDYFlag) THEN !if calculating DSDY, then calculate PNXOLD
                DO L = 1, MSection
                    IF(SootSecO(I,J,L).LE.smallnum .OR. SootSecO(I,J,L+MSection).LE.smallnum) THEN
                        PNXOLD(L)=1.d0
                    ELSE
                        PNXOLD(L)=SootSecO(I,J,L+MSection)/(SootSecO(I,J,L))
                    ENDIF
                    IF(PNXOLD(L).LT.1.d0) THEN
                        PNXOLD(L) = 1.d0
                    ENDIF 
                ENDDO
            ENDIF                  
            
            !****************************************************************
            !first, determine nucleation rates
            TotalNucRate = 0.d0
            DO L=1,NumNuc
                TotalNucRate = TotalNucRate + SootNucRates(I,J,L+NumNuc)
            ENDDO      
            NucRate(1) = TotalNucRate/DEXP(XS(1))
            NucRate(1+MSection) = TotalNucRate/DEXP(XS(1))
            !done determining nucleation rates
            !*****************************************************************          
            
            !***********************************************************************
            !second, determine coagulation rates, if dsdy is not true and there is soot present
                
            IF(.NOT.DSDYFLAG) THEN!if dsdyflag is false, calculate coagulation pre-factors                                                 
                              !if true, use prefactor values passed into subroutine                                          
                !get collision rates
                DO R=1,MSection
                    DO S=R,MSection
                        CALL BETA_AEROSOL(DUMMY,PNX(R),R,DUMMY,PNX(S),S,CollRate(R,S))
                    ENDDO
                ENDDO           
                     
                !coag. of sections r, s <= L, mass->L
                DO L=1,MSection
                    DO R=1,L
                        DO S=R,L
                            W=DEXP(XS(L))
                            IF (L.EQ.1) THEN
                                WL=W
                            ELSE
                                WL=DEXP(XS(L-1))
                            ENDIF
                            IF (L.EQ.MSection) THEN
                                WU=W*FS
                            ELSE
                                WU=DEXP(XS(L+1))
                            ENDIF
                            
                            UTOT = DEXP(XS(R))+ DEXP(XS(S))

                            IF(UTOT .GT. WL .AND. UTOT .LT. W) THEN
                                CS1(R,S,L)=(WL-UTOT)/(WL-W)*CollRate(R,S)
                                CNS1(R,S,L)=CS1(R,S,L)*W/UTOT
                                IF(R.EQ.S) THEN
                                    CS1(R,S,L)=CS1(R,S,L)/2.d0
                                    CNS1(R,S,L)=CNS1(R,S,L)/2.d0
                                ENDIF
                            ELSE
                                IF(UTOT .GE. W .AND. UTOT .LT. WU) THEN                             
                                    CS1(R,S,L)=(WU-UTOT)/(WU-W)*CollRate(R,S)
                                    CNS1(R,S,L)=CS1(R,S,L)*W/UTOT
                                    IF(R.EQ.S) THEN
                                        CS1(R,S,L)=CS1(R,S,L)/2.d0
                                        CNS1(R,S,L)=CNS1(R,S,L)/2.d0
                                    ENDIF
                                ELSE
                                    CS1(R,S,L)=0.0d0
                                    CNS1(R,S,L)=0.0d0
                                ENDIF
                            ENDIF
                            CS1(S,R,L)=CS1(R,S,L)
                            CNS1(S,R,L)=CNS1(R,S,L)
                        ENDDO
                    ENDDO
                ENDDO 
                
                !destruction of aggregates and primaries due to coagulation
                DO L=1,MSection
                    DO S=L,MSection                  
                        CS2(S,L)=0.d0-CollRate(L,S)
                        CS2(L,S)=CS2(S,L)
                        CNS2(S,L)=CS2(S,L)
                        CNS2(L,S)=CS2(L,S)
                    ENDDO
                ENDDO
            
            ENDIF !done getting rates 
                      
            !now multiply by number of aggregates in the colliding sections
!            IF(.Not.dsdyflag) THEN !if dsdyflag is false, calculate all coagulation rates
                DO L=1,MSection

                    DO R=1,L
                        DO S=R,L
                            CoagRate(L)=CoagRate(L)+CS1(R,S,L)*SootSec(I,J,R)*SootSec(I,J,S)*RHO(I,J)**2.d0
                            CoagRate(MSection+L)=CoagRate(MSection+L)+CNS1(R,S,L)*SootSec(I,J,R)*SootSec(I,J,S)*RHO(I,J)**2.d0 &
                                                 *(PNX(R)+PNX(S))
                        ENDDO
                    ENDDO        

                    DO S=1,MSection
                        CoagRate(L)=CoagRate(L)+CS2(S,L)*SootSec(I,J,S)*SootSec(I,J,L)*RHO(I,J)**2.d0
                        CoagRate(MSection+L)=CoagRate(MSection+L)+CNS2(S,L)*SootSec(I,J,S)*SootSec(I,J,L)*RHO(I,J)**2.d0 &
                                                 *(PNX(L))
                    ENDDO
                    
                ENDDO  
                SootCoagRates(I,J,:) = CoagRate(:)                                                                                                                
            !done determining coagulation rates
            !**************************************************************************************** 
            
            !*************************************************************************************************
            !determine surface growth rates due to the HACA mechanism, and PAH condensation
            
            !   L=1, first section
            AINPAH=0.0d0
            AINC2H2=0.0d0

            !sum up mass added due to PAH condensation
            if(Mechanism.LE.3) then
                temp = SootGasRates(I,J,1,7)*20.d0*C_MW + SootGasRates(I,J,1,8)*18.d0*C_MW + SootGasRates(I,J,1,9)*20.d0*C_MW  !7=BAPYR, 8=BGHIF, 9=BAPYRS
            else
                temp = SootGasRates(I,J,1,7)*16.d0*C_MW  !7=A4
            end if

            AOUTPAH =temp /DEXP(XS(1))/(FS-1) !number of aggregates leaving section 1 to section 2 due to PAH condensation

            AOUTC2H2=SootGasRates(I,J,1,4)*2.d0*C_MW/DEXP(XS(1))/(FS-1) !number of aggregates leaving section 1 to section 2 due to C2H2 addition
            AINO2   =SootGasRates(I,J,2,5)*2.d0*C_MW/DEXP(XS(1))/(FS-1) !number of aggregates coming into section 1 from section 2 due to O2 oxidation
            AINOH   =SootGasRates(I,J,2,6)*1.d0*C_MW /DEXP(XS(1))/(FS-1)!number of aggregates coming into section 1 from section 2 due to OH oxidation
            AOUTO2  =SootGasRates(I,J,1,5)*2.d0*C_MW /DEXP(XS(1))      !number of aggregates leaving section 1 and being destroyed due to O2 oxidation
            AOUTOH  =SootGasRates(I,J,1,6)*1.d0*C_MW /DEXP(XS(1)) !number of aggregates leaving section 1 and being destroyed due to OH oxidation

            SurfRate(1)=AINPAH+AINC2H2+AINO2+AINOH-AOUTPAH-AOUTC2H2-AOUTO2-AOUTOH  !total number of aggregates coming to section 1 (negative number means there is a net leaving)

            !determine the flux of primary particles in section 1 due to flux of aggregates
            SurfRate(MSection+1) =-PNX(1)*AOUTPAH-PNX(1)*AOUTC2H2 +AINO2*PNX(2)-AOUTO2*PNX(1)+AINOH*PNX(2)-AOUTOH*PNX(1)

            !now do section 2 to (MSection-1)            
            DO L=2,MSection-1

                AINPAH=AOUTPAH  !the number of aggregates entering section L due to PAH condensation is the number leaving section L-1 due to PAH condensation

                AINC2H2=AOUTC2H2 !the number of aggregates entering section L due to C2H2 addition is the number leaving section L-1 due to C2H2 addition
                AOUTO2=AINO2    !the number of aggregates leaving section L due to O2 oxidation is the number going to section L-1 due to O2 oxidation                
                AOUTOH=AINOH     !the number of aggregates leaving section L due to OH oxidation is the number going to section L-1 due to OH oxidation
                
                !sum up mass added due to PAH condensation
                if(Mechanism.LE.3) then
                    temp = SootGasRates(I,J,L,7)*20.d0*C_MW + SootGasRates(I,J,L,8)*18.d0*C_MW + SootGasRates(I,J,L,9)*20.d0*C_MW  !7=BAPYR, 8=BGHIF, 9=BAPYRS
                else
                    temp = SootGasRates(I,J,1,7)*16.d0*C_MW  !7=A4
                end if

                AOUTPAH =temp /DEXP(XS(L))/(FS-1) !number of aggregates leaving section L to section L+1 due to PAH condensation

                AOUTC2H2=SootGasRates(I,J,L,4)*2.d0*C_MW/DEXP(XS(L))/(FS-1) !number of aggregates leaving section L to section L+1 due to C2H2 addition
                AINO2   = SootGasRates(I,J,L+1,5)*2.d0*C_MW/DEXP(XS(L))/(FS-1) !number of aggregates coming into section L from section L+1 due to O2 oxidation
                AINOH   =SootGasRates(I,J,L+1,6)*1.d0*C_MW /DEXP(XS(L))/(FS-1)!number of aggregates coming into section L from section L+1 due to OH oxidation

                SurfRate(L)=AINPAH+AINC2H2+AINO2+AINOH-AOUTPAH-AOUTC2H2-AOUTO2-AOUTOH  !total number of aggregates coming to section L (negative number means there is a net leaving)


                !determine the flux of primary particles in section 1 due to flux of aggregates
                SurfRate(MSection+L) = AINPAH*PNX(L-1) -AOUTPAH*PNX(L) + AINC2H2*PNX(L-1)-AOUTC2H2*PNX(L) &
                                       + AINO2*PNX(L+1)  -AOUTO2*PNX(L) + AINOH*PNX(L+1)  -AOUTOH*PNX(L)
            
            ENDDO

            !and now do the final section            
            SurfRate(MSection) = AOUTPAH + AOUTC2H2 - AINO2 - AINOH
            
            !determine the flux of primary particles in section 1 due to flux of aggregates
            SurfRate(2*MSection) = AOUTPAH*PNX(MSection-1)+ AOUTC2H2*PNX(MSection-1)-AINO2*PNX(MSection)-AINOH*PNX(MSection)     
            
            !done getting surface growth rates
            !****************************************************************************************************************************************************
            
            !******************************************************************************************************************
            !get fragmentation rates
            
            !first, get the total surface growth rate due to HACA growth, PAH condensation, and the total oxidation rate
            CALL SootProRate(total_rate_for_frag,specific_ox_rate)
            
            total_rate_for_frag = total_rate_for_frag - specific_ox_rate !this is the total rate due to surface growth and oxidation
            
            !get the total surface area density 
            area = 0.d0
            DO L=1,MSection                                                            
                Area=Area+PNX(L)*4.D0*PI*(3.D0*(DEXP(XS(L))/PNX(L))/4.D0/PI/densityP)**(2.d0/3.d0)*SootSec(I,J,L)*RHO(I,J)   !cm^2/cm^3
            ENDDO
         
            tao_i_i = (FS-2.0d0)/(FS-1.0d0) 
            tao_i_ip1 = FS/(FS-1.0d0)   
       
            if(total_rate_for_frag.LT.0.d0) then !if there is a net amount of oxidation occuring
            
                specific_ox_rate = specific_ox_rate /(area + 1.0d-50)  !(g/cm^3/sec)/(cm^2/cm^3) = g/cm^2/sec
            
                A_frag = 1.0d5*specific_ox_rate ! the overall frag rate A (in S=A*(np)^(1/Df)) is the first order function of the specific oxidation rate 
            
                S_frag(1)=0.0d0
                do L=2,MSection
                    S_frag(L)=A_frag*PNX(L)**(1.0d0/DFRCT)
                    if (PNX(L) .lt. 2.0d0) then
                        S_frag(L)=0.0d0
                    endif
                enddo
                
                FragRate(1)=tao_i_ip1*S_frag(2)*SootSec(I,J,2)*RHO(I,J)
                FragRate(1+MSection)=tao_i_ip1*S_frag(2)*SootSec(I,J,2)*RHO(I,J)*PNX(2)/FS
                do L=2,MSection-1
                    FragRate(L)= (tao_i_i - 1.0d0) * S_frag(L)*SootSec(I,J,L)*RHO(I,J) + tao_i_ip1*S_frag(L+1) &
                                *SootSec(I,J,L+1)*RHO(I,J)
                    FragRate(L+MSection)= (tao_i_i - 1.0d0) * S_frag(L)*SootSec(I,J,L)*RHO(I,J)*PNX(L) + tao_i_ip1 &
                                *S_frag(L+1)*SootSec(I,J,L+1)*RHO(I,J)*PNX(L+1)/FS
                enddo
                FragRate(MSection)= (tao_i_i - 1.0d0) * S_frag(MSection)*SootSec(I,J,MSection)*RHO(I,J)
                FragRate(2*MSection)= (tao_i_i - 1.0d0) * S_frag(MSection)*SootSec(I,J,MSection)*RHO(I,J)*PNX(MSection)

            endif
            !done getting fragmentation rates
            !****************************************************************************************************************************
            
            DO L=1,2*MSection
            
                SOURCETERM(L) = NucRate(L) + CoagEff*CoagRate(L) + SurfRate(L) + FragRate(L) + ObliRate(L) + SintRate(L)    
            
            ENDDO     
            
            deallocate(NucRate)
            deallocate(PNX)
            deallocate(S_frag)
            deallocate(FragRate)
            deallocate(SurfRate)
            deallocate(SintRate)
            deallocate(ObliRate)
            deallocate(CollRate)
            deallocate(CoagRate)
            deallocate(PNXOLD)
            
        end subroutine SootSecRate
        !*********************************************************************************************************************************************************************
        
        !*******************************************************************************************************************************************************************
        SUBROUTINE SootProRate(SootGrowRate,SootOxRate)
        !This subroutine calculates the total rate of soot growth and oxidation production over all sections. The result is returned in g-soot/cm^3/s.
        !A positive oxidation rate means destruction of soot, and a positive growth rate means production of soot. This does not include soot production rates due to nucleation
        
            !local variable declaration
            DOUBLE PRECISION SootGrowRate, SootOxRate
            INTEGER :: L
            !done local variable declaration
            
            !initialize rates
            SootGrowRate = 0.d0
            SootOxRate = 0.d0
        
            DO L=1,MSection
                SootOxRate = SootOxRate + SootGasRates(I,J,L,5)*2.d0*C_MW + SootGasRates(I,J,L,6)*1.d0*C_MW 
                if(Mechanism.LE.3) then
                    SootGrowRate = SootGrowRate + SootGasRates(I,J,L,4)*2.d0*C_MW + SootGasRates(I,J,L,7)*20.d0*C_MW + &
                                                  SootGasRates(I,J,L,8)*18.d0*C_MW + SootGasRates(I,J,L,9)*20.d0*C_MW
                else    !Pyrene condensation
                    SootGrowRate = SootGrowRate + SootGasRates(I,J,L,4)*2.d0*C_MW + SootGasRates(I,J,L,7)*16.d0*C_MW !A4=7
                end if
            ENDDO
        
        end subroutine SootProRate
        !*****************************************************************************************************************************************      
        
        !*********************************************************************************************************************************************
        SUBROUTINE Pdiff(Sec, Pdiffusion)
        !This subroutine calculates the ordinary diffusion coefficient (Pdiffusion) for soot aggregates of the given section, Sec, for the CV defined by
        !the global I and J variables.
        
            !local variable declarations
            DOUBLE PRECISION :: Viscosity, freePath
            DOUBLE PRECISION :: Pdiffusion, PNX, RXP, RFX, RSX, RMFX, RMCX, FDX, CDX, VISCOS, FREEMP
            INTEGER :: Sec
            !done local variable declarations
            
            !get the number of primary particles per aggregate
            IF(SootSec(I,J,Sec).LE.smallnum.OR.SootSec(I,J,Sec+MSection).LE.smallnum) THEN !if number density of aggregates or primaries is low, assume 1 primary per aggregate
                PNX=1.D0
            ELSE
                PNX=SootSec(I,J,MSection+Sec)/SootSec(I,J,Sec) !else, just Np/Na
            ENDIF 
            !done getting primary particles per aggregate 
            
            !get viscosity and mean free path
            VISCOS = viscosity(T(I,J))
            FREEMP = freePath(Pressure,T(I,J))
             
            
            !get primary particle radius
            RXP=(3.d0*(DEXP(XS(Sec))/PNX)/4.d0/PI/densityP)**(1.d0/3.d0)
            
            IF(PNX.GT.1.d0) THEN
                RFX=RXP*(FVOL*PNX)**(1.0d0/DFRCT)
                ! Compute the free-molecule mobility diameter
                ! which is same as the projected area equivalent diameter
                RSX=RXP*PNX**0.43d0
                RMFX=RSX
                ! Now compute the continuum mobility diameter
                RMCX=RFX*(-0.06483d0*DFRCT**2.d0+0.6353d0*DFRCT-0.4898d0)
                IF(RMFX.GT.RMCX) RMCX=RMFX
            ELSE
                RFX=(3.d0*DEXP(XS(Sec))/4.d0/PI/densityP)**(1.d0/3.d0)
                RSX=RFX
                RMFX=RSX
                RMCX=RSX
            END IF
            ! NOW COMPUTE THE DRAG ON THE CLUSTERS
            FDX=VISCOS*RMFX**2.d0/FREEMP
            CDX=1.612d0*VISCOS*RMCX
            !  INTERPOLATE IN THE TRANSITION REGIME BY ADDING MOBILITIES
            Pdiffusion=1.612d0*BOLTZMANN*T(I,J)/(6.d0*PI)*(1.d0/CDX+1.d0/FDX)
        
        end subroutine PDiff
        !******************************************************************************************************************************************
        
        !******************************************************************************************************************************************
        SUBROUTINE BETA_AEROSOL (Rad1, Mass1, flag1, Rad2, Mass2, flag2,Collisions)
        !This subroutine calculates the collision rate between gas species and soot, soot and soot, or two gas species using aerosol collision theory. 
        !The flags define what is to be calculated. Flag = 0 defines that the entity is a gas phase species. Flag \= 0 means it is a aggregate of soot from section Flag. 
        !When Flag = 0, Rad = Radius of colliding gas phase species, and Mass = mass of gas phase species
        !When Flag \= 0, Rad = dummy variable, and Mass = the number of primaries per aggregate in section Flag.
        !The collision rate is returned in the units of Collisions/(cm^3*s)/(#/cm^3)/(#/cm^3), or collision rate per cc on a per number of 
        !each colliding entity.
        
            !local variable declaration
            INTEGER :: flag1, flag2
            DOUBLE PRECISION :: Rad1, Mass1, Rad2, Mass2, Collisions
            DOUBLE PRECISION :: Umass, Vmass, mred
            DOUBLE PRECISION :: Viscosity, freePath
            DOUBLE PRECISION :: VISCOS, FREEMP, R3, C1, C2, C3, C4, ALPHA, FALPHA
            DOUBLE PRECISION :: RCMAX, RCMIN, NMIN, RABS, DIFFN, CMEAN, KND, DAHNEKE
            DOUBLE PRECISION :: RVX, RMCX, RMFX, RXP, DFX, RFX, RSX, RCX, FDX, CDX, DIFFNX
            DOUBLE PRECISION :: RVY, RMCY, RMFY, RYP, DFY, RFY, RSY, RCY, FDY, CDY, DIFFNY
            !done local variable declaration
          
            !get viscosity and mean free path
            VISCOS = viscosity(T(I,J))
            FREEMP = freePath(Pressure,T(I,J))                      
            
            !Get the mass and radius of the entities that are colliding
            IF(Flag1.EQ.0) then
                Umass=Mass1
                RVX = Rad1
                RFX = Rad1
                RCX = Rad1
                RMCX = Rad1
                RMFX = Rad1
            ELSE      
                Umass=DEXP(XS(flag1)) !the mass of the second entity is the mass of an aggregate of section flag1
                RXP = (3.d0*(DEXP(XS(flag1))/Mass1)/4.d0/PI/densityP)**(1.d0/3.d0) !diameter of the primary particles
            ENDIF
            
            IF(Flag2.EQ.0) then
                Vmass=Mass2
                RVY = Rad2
                RFY = Rad2
                RCY = Rad2
                RMCY = Rad2
                RMFY = Rad2
            ELSE      
                Vmass=DEXP(XS(flag2)) !the mass of the second entity is the mass of an aggregate of section flag1
                RYP = (3.d0*(DEXP(XS(flag2))/Mass2)/4.d0/PI/densityP)**(1.d0/3.d0) !diameter of the primary particles
            ENDIF      
            
            IF(Flag1.GT.0) then !if the first entity is a soot aggregate
            
                IF(Mass1.GT.1) THEN !if there is more than one primary per aggregate 
                    DFX = DFRCT
                    RFX=RXP*(FVOL*Mass1)**(1.0d0/DFRCT)
                    ! Compute the free-molecule mobility diameter
                    ! which is same as the projected area equivalent diameter
                    RSX=RXP*Mass1**0.43d0
                    RMFX=RSX
                    ! Now compute the continuum mobility diameter
                    RMCX=RFX*(-0.06483d0*DFRCT**2.d0+0.6353d0*DFRCT-0.4898d0)
                    IF(RMFX.GT.RMCX) RMCX=RMFX
                    RVX=RXP*Mass1**(1.d0/3.d0)
                ELSE
                    DFX = 3.0d0
                    RFX=(3.0d0*Umass/4.0d0/PI/densityP)**(1.d0/3.d0)
                    RSX=RFX
                    RMFX=RSX
                    RMCX=RSX
                    RVX=RFX
                END IF
            ENDIF
            
            IF(Flag2.GT.0) then !if the first entity is a soot aggregate        
!
                IF(Mass2.GT.1) THEN !if there is more than one primary per aggregate
                    DFY = DFRCT
                    RFY=RYP*(FVOL*Mass2)**(1.0d0/DFRCT)
                    RSY=RYP*Mass2**0.43d0
                    RMFY=RSY
                    RMCY=RFY*(-0.06483d0*DFY**2.d0+0.6353d0*DFRCT-0.4898d0)
                    IF(RMFY.GT.RMCY) RMCY=RMFY
                    RVY=RYP*Mass2**(1.d0/3.d0)
                ELSE
                    DFY = 3.0d0
                    RFY=(3.0d0*Vmass/4.0d0/PI/densityP)**(1.d0/3.d0)
                    RSY=RFY
                    RMFY=RSY
                    RMCY=RSY
                    RVY=RFY
                END IF
            ENDIF         
            
            ! Now compute the absorbing sphere radius.
        
            ! First, compute individual collision radius   
            IF(Flag1.GT.0) then
                IF(Mass1.LE.3.d0) THEN
                    R3=RXP*(2.0405d0-0.1492d0*AKF)
                    RCX=RXP+(R3-RXP)*(Mass1-1.d0)/2.d0
                ELSE
                    RCX=RXP*(DSQRT(DFX/(DFX+2.))*(FVOL*Mass1)**(1.d0/DFX)*   &
                        (1.037d0*AKF**0.077d0+((2.0405d0-0.1492d0*AKF)/DSQRT(DFX/(DFX+2.d0))/ &
                        (FVOL*Mass1)**(1.d0/DFX)-1.037d0*AKF**0.077d0)*DEXP(-((Mass1-3.d0)/17.d0)**0.522d0)))
                ENDIF
            ENDIF
        
            IF(Flag2.GT.0) then       
                IF(Mass2.LE.3.d0) THEN
                    R3=RYP*(2.0405d0-0.1492d0*AKF)
                    RCY=RYP+(R3-RYP)*(Mass2-1.d0)/2.d0
                ELSE
                    RCY=RYP*(DSQRT(DFY/(DFY+2.d0))*(FVOL*Mass2)**(1.d0/DFY)* &
                        (1.037d0*AKF**0.077d0+((2.0405d0-0.1492d0*AKF)/DSQRT(DFY/(DFY+2.d0))/ &
                        (FVOL*Mass2)**(1.d0/DFY)-1.037d0*AKF**0.077d0)*DEXP(-((Mass2-3.d0)/17.d0)**0.522d0)))
                ENDIF          
            ENDIF
      
            !if at least one entity is a soot aggregate
            IF(Flag1.GT.0.OR.Flag2.GT.0) then
            
                !compute the penetration factor
                C1=0.296d0
                C2=0.412d0
                C3=4.427d0+3.697d0*(AKF-1.3d0)
                C4=0.897d0+0.013d0*(AKF-1.3d0)
           
                IF(Flag1.GT.0.AND.Flag2.GT.0) then  !if both entities are soot aggregates   
                
                    IF(Umass.GT.Vmass) THEN
                        ALPHA=Mass2/Mass1
                        RCMAX=RCX
                        RCMIN=RCY
                        NMIN=Mass2
                    ELSE
                        ALPHA=Mass1/Mass2
                        RCMAX=RCY
                        RCMIN=RCX
                        NMIN=Mass1
                    END IF
                    
                    !penetrations factor
                    FALPHA=(-ALPHA+DSQRT(4.d0*ALPHA)-1.d0)*DEXP((1.d0-DFX)/C1)+1.d0   &
                            -((-ALPHA+DSQRT(4.d0*ALPHA)-1.d0)*DEXP((1.d0-DFX)/C1)+1.d0 &
                            -C4*ALPHA**(1.d0/DFX-0.5d0))*DEXP(-C3*ALPHA**C2)
                    
                    RABS=RCMAX*FALPHA+RCMIN !absorbing sphere radius with penetration factor 
                    
                    ! Now take into account rotation effect
                    RABS=RABS/(0.085d0*DEXP(-3.42d0*ALPHA**2.d0)+0.9054d0+0.18d0*DEXP(-0.41d0*NMIN**(1.d0/3.d0)-1.d0/ALPHA))
                    
                ELSE
               
                    IF(Flag1.GT.0) then
                        RCMAX = RCX
                        RCMIN = RCY
                        ALPHA = 1.d0/Mass1
                        FALPHA=(-ALPHA+DSQRT(4.d0*ALPHA)-1.d0)*DEXP((1.d0-DFX)/C1)+1.d0   &
                            -((-ALPHA+DSQRT(4.d0*ALPHA)-1.d0)*DEXP((1.d0-DFX)/C1)+1.d0 &
                            -C4*ALPHA**(1.d0/DFX-0.5d0))*DEXP(-C3*ALPHA**C2)
                    ELSE
                        RCMAX = RCY
                        RCMIN = RCX
                        ALPHA = 1.d0/Mass2
                        FALPHA=(-ALPHA+DSQRT(4.d0*ALPHA)-1.d0)*DEXP((1.d0-DFY)/C1)+1.d0   &
                            -((-ALPHA+DSQRT(4.d0*ALPHA)-1.d0)*DEXP((1.d0-DFY)/C1)+1.d0 &
                            -C4*ALPHA**(1.d0/DFY-0.5d0))*DEXP(-C3*ALPHA**C2)
                    ENDIF
                    
                    RABS=RCMAX*FALPHA+RCMIN !absorbing sphere radius with penetration factor 
                                                               
                ENDIF              
                
                ! Outer radius is the upper bound for collision radius
                IF(RABS.GT.(RFX+RFY)) THEN
                    RABS=RFX+RFY
                ! Volume equivalent radius is the lower bound for collision radius
                ELSE IF(RABS.LT.(RVX+RVY)) THEN
                    RABS=RVX+RVY
                ENDIF
                
            ELSE ! if neither entities are soot aggregates
                RABS = RCX + RCY
            ENDIF       

            ! NOW COMPUTE THE DRAG ON THE ENTITIES
            FDX=VISCOS*RMFX**2/FREEMP
            FDY=VISCOS*RMFY**2/FREEMP
            CDX=1.612d0*VISCOS*RMCX
            CDY=1.612d0*VISCOS*RMCY

            !  INTERPOLATE IN THE TRANSITION REGIME BY ADDING MOBILITIES
            DIFFNX=1.612d0*BOLTZMANN*T(I,J)/(6.d0*PI)*(1.0d0/CDX+1.0d0/FDX)
            DIFFNY=1.612d0*BOLTZMANN*T(I,J)/(6.d0*PI)*(1.0d0/CDY+1.0d0/FDY)
            DIFFN=DIFFNX+DIFFNY         !(D1+D2)
            ! COMPUTE THE KNUDSEN NUMBER AND INTERPOLATION FACTOR
            CMEAN=2.2d0*DSQRT(2.d0*BOLTZMANN*T(I,J)/PI*(1.0d0/Umass+1.0d0/Vmass))
            KND=DIFFN/CMEAN/RABS
            DAHNEKE=(1.0d0+KND)/(1.0d0+2.0d0*KND*(1.0d0+KND))   !fD
            ! Now put it all together
            Collisions=4.d0*PI*RABS*DIFFN*DAHNEKE
            
            if(Collisions.LE.0.d0) Collisions = 0.d0      
            
        end subroutine BETA_AEROSOL
        !**************************************************************************************************************************

!******************************************************************************************************************************************
        SUBROUTINE BETA_FREE (Rad1, Mass1, flag1, Rad2, Mass2, flag2,Collisions)
        !This subroutine calculates the collision rate between gas species and soot, soot and soot, or two gas species using free molecular regime theory. 
        !The flags define what is to be calculated. Flag = 0 defines that the entity is a gas phase species. Flag \= 0 means it is a aggregate of soot from section Flag. 
        !When Flag = 0, Rad = Radius of colliding gas phase species, and Mass = mass of gas phase species
        !When Flag \= 0, Rad = dummy variable, and Mass = the number of primaries per aggregate in section Flag.
        !The collision rate is returned in the units of Collisions/(cm^3*s)/(#/cm^3)/(#/cm^3), or collision rate per cc on a per number of 
        !each colliding entity.
        
            !local variable declaration
            INTEGER :: flag1, flag2
            DOUBLE PRECISION :: Rad1, Mass1, Rad2, Mass2, Collisions
            DOUBLE PRECISION :: Umass, Vmass, mred
            DOUBLE PRECISION :: Viscosity, freePath
            DOUBLE PRECISION :: VISCOS, FREEMP, R3, C1, C2, C3, C4, ALPHA, FALPHA
            DOUBLE PRECISION :: RCMAX, RCMIN, NMIN, RABS, DIFFN, CMEAN, KND, DAHNEKE
            DOUBLE PRECISION :: RVX, RMCX, RMFX, RXP, DFX, RFX, RSX, RCX, FDX, CDX, DIFFNX
            DOUBLE PRECISION :: RVY, RMCY, RMFY, RYP, DFY, RFY, RSY, RCY, FDY, CDY, DIFFNY
            !done local variable declaration
          
            !get viscosity and mean free path
            VISCOS = viscosity(T(I,J))
            FREEMP = freePath(Pressure,T(I,J))                      
            
            !Get the mass and radius of the entities that are colliding
            IF(Flag1.EQ.0) then
                Umass=Mass1
                RVX = Rad1
                RFX = Rad1
                RCX = Rad1
                RMCX = Rad1
                RMFX = Rad1
            ELSE      
                Umass=DEXP(XS(flag1)) !the mass of the second entity is the mass of an aggregate of section flag1
                RXP = (3.d0*(DEXP(XS(flag1))/Mass1)/4.d0/PI/densityP)**(1.d0/3.d0) !diameter of the primary particles
            ENDIF
            
            IF(Flag2.EQ.0) then
                Vmass=Mass2
                RVY = Rad2
                RFY = Rad2
                RCY = Rad2
                RMCY = Rad2
                RMFY = Rad2
            ELSE      
                Vmass=DEXP(XS(flag2)) !the mass of the second entity is the mass of an aggregate of section flag1
                RYP = (3.d0*(DEXP(XS(flag2))/Mass2)/4.d0/PI/densityP)**(1.d0/3.d0) !diameter of the primary particles
            ENDIF      
            
            IF(Flag1.GT.0) then !if the first entity is a soot aggregate
            
                IF(Mass1.GT.1) THEN !if there is more than one primary per aggregate 
                    DFX = DFRCT
                    RFX=RXP*(FVOL*Mass1)**(1.0d0/DFRCT)
                    ! Compute the free-molecule mobility diameter
                    ! which is same as the projected area equivalent diameter
                    RSX=RXP*Mass1**0.43d0
                    RMFX=RSX
                    ! Now compute the continuum mobility diameter
                    RMCX=RFX*(-0.06483d0*DFRCT**2.d0+0.6353d0*DFRCT-0.4898d0)
                    IF(RMFX.GT.RMCX) RMCX=RMFX
                    RVX=RXP*Mass1**(1.d0/3.d0)
                ELSE
                    DFX = 3.0d0
                    RFX=(3.0d0*Umass/4.0d0/PI/densityP)**(1.d0/3.d0)
                    RSX=RFX
                    RMFX=RSX
                    RMCX=RSX
                    RVX=RFX
                END IF
            ENDIF
            
            IF(Flag2.GT.0) then !if the first entity is a soot aggregate        
!
                IF(Mass2.GT.1) THEN !if there is more than one primary per aggregate
                    DFY = DFRCT
                    RFY=RYP*(FVOL*Mass2)**(1.0d0/DFRCT)
                    RSY=RYP*Mass2**0.43d0
                    RMFY=RSY
                    RMCY=RFY*(-0.06483*DFY**2.d0+0.6353d0*DFRCT-0.4898d0)
                    IF(RMFY.GT.RMCY) RMCY=RMFY
                    RVY=RYP*Mass2**(1.d0/3.d0)
                ELSE
                    DFY = 3.0d0
                    RFY=(3.0d0*Vmass/4.0d0/PI/densityP)**(1.d0/3.d0)
                    RSY=RFY
                    RMFY=RSY
                    RMCY=RSY
                    RVY=RFY
                END IF
            ENDIF          

            
            ! Now compute the absorbing sphere radius.
        
            ! First, compute individual collision radius   
            IF(Flag1.GT.0) then
                IF(Mass1.LE.3.d0) THEN
                    R3=RXP*(2.0405d0-0.1492d0*AKF)
                    RCX=RXP+(R3-RXP)*(Mass1-1.d0)/2.d0
                ELSE
                    RCX=RXP*(DSQRT(DFX/(DFX+2.))*(FVOL*Mass1)**(1.d0/DFX)*   &
                        (1.037d0*AKF**0.077d0+((2.0405d0-0.1492d0*AKF)/DSQRT(DFX/(DFX+2.d0))/ &
                        (FVOL*Mass1)**(1.d0/DFX)-1.037d0*AKF**0.077d0)*DEXP(-((Mass1-3.d0)/17.d0)**0.522d0)))
                ENDIF
            ENDIF
        
            IF(Flag2.GT.0) then       
                IF(Mass2.LE.3.d0) THEN
                    R3=RYP*(2.0405d0-0.1492d0*AKF)
                    RCY=RYP+(R3-RYP)*(Mass2-1.d0)/2.d0
                ELSE
                    RCY=RYP*(DSQRT(DFY/(DFY+2.d0))*(FVOL*Mass2)**(1.d0/DFY)* &
                        (1.037d0*AKF**0.077d0+((2.0405d0-0.1492d0*AKF)/DSQRT(DFY/(DFY+2.d0))/ &
                        (FVOL*Mass2)**(1.d0/DFY)-1.037d0*AKF**0.077d0)*DEXP(-((Mass2-3.d0)/17.d0)**0.522d0)))
                ENDIF          
            ENDIF
      
            !if at least one entity is a soot aggregate
            IF(Flag1.GT.0.OR.Flag2.GT.0) then
            
                !compute the penetration factor
                C1=0.296d0
                C2=0.412d0
                C3=4.427d0+3.697d0*(AKF-1.3d0)
                C4=0.897d0+0.013d0*(AKF-1.3d0)
           
                IF(Flag1.GT.0.AND.Flag2.GT.0) then  !if both entities are soot aggregates   
                
                    IF(Umass.GT.Vmass) THEN
                        ALPHA=smallnum
                        RCMAX=RCX
                        RCMIN=RCY
                        NMIN=Mass2
                    ELSE
                        ALPHA=smallnum
                        RCMAX=RCY
                        RCMIN=RCX
                        NMIN=Mass1
                    END IF
                    
                    !penetrations factor
                    FALPHA=(-ALPHA+DSQRT(4.d0*ALPHA)-1.d0)*DEXP((1.d0-DFX)/C1)+1.d0   &
                            -((-ALPHA+DSQRT(4.d0*ALPHA)-1.d0)*DEXP((1.d0-DFX)/C1)+1.d0 &
                            -C4*ALPHA**(1.d0/DFX-0.5d0))*DEXP(-C3*ALPHA**C2)
                    
                    RABS=RCMAX*FALPHA+RCMIN !absorbing sphere radius with penetration factor 
                    
                    ! Now take into account rotation effect
                    RABS=RABS/(0.085d0*DEXP(-3.42d0*ALPHA**2d0)+0.9054d0+0.18d0*DEXP(-0.41d0*NMIN**(1.d0/3.d0)-1.d0/ALPHA))
                    
                ELSE
               
                    IF(Flag1.GT.0) then
                        IF(Mass1.GT.1d0) THEN
                            RCMAX = RCX
                            RCMIN = RCY
                            ALPHA = smallnum
                            FALPHA=(-ALPHA+DSQRT(4.d0*ALPHA)-1.d0)*DEXP((1.d0-DFX)/C1)+1.d0   &
                                -((-ALPHA+DSQRT(4.d0*ALPHA)-1.d0)*DEXP((1.d0-DFX)/C1)+1.d0 &
                                -C4*ALPHA**(1.d0/DFX-0.5d0))*DEXP(-C3*ALPHA**C2)
                        ELSE
                            FALPHA = 1.d0
                            RCMAX = RCX
                            RCMIN = RCY
                        ENDIF                
                    ELSE
                        IF(Mass2.GT.1d0) THEN
                            RCMAX = RCY
                            RCMIN = RCX
                            ALPHA = smallnum
                            FALPHA=(-ALPHA+DSQRT(4.d0*ALPHA)-1.d0)*DEXP((1.d0-DFY)/C1)+1.d0   &
                                -((-ALPHA+DSQRT(4.d0*ALPHA)-1.d0)*DEXP((1.d0-DFY)/C1)+1.d0 &
                                -C4*ALPHA**(1.d0/DFY-0.5d0))*DEXP(-C3*ALPHA**C2)
                        ELSE
                            RCMAX = RCY
                            RCMIN = RCX
                            FALPHA = 1.d0
                        ENDIF               
                    ENDIF
                    
                    RABS=RCMAX*FALPHA+RCMIN !absorbing sphere radius with penetration factor 
                                                               
                ENDIF              
                
                ! Outer radius is the upper bound for collision radius
                IF(RABS.GT.(RFX+RFY)) THEN
                    RABS=RFX+RFY
                ! Volume equivalent radius is the lower bound for collision radius
                ELSE IF(RABS.LT.(RVX+RVY)) THEN
                    RABS=RVX+RVY
                ENDIF
                
            ELSE ! if neither entities are soot aggregates
                RABS = RCX + RCY
            ENDIF     
            !get the reduced mass
            mred = Umass*Vmass/(Umass+Vmass)       
            !Now calculate the collision frequency
            Collisions=RABS**2.d0*DSQRT(8.d0*PI*BOLTZMANN*T(I,J)/mred)
            
            if(Collisions.LE.0.d0) Collisions = 0.d0      
            
        end subroutine BETA_FREE
        !**************************************************************************************************************************        
        
        !******************************************************************************************************************************************
        SUBROUTINE DOM(TEMP,XCO,XCO2,XH2O,SFV,QR,XUIN,YVIN,PRESSURE,SOLIDS)       
         
            DOUBLE PRECISION, DIMENSION (:,:) :: TEMP, XCO, XCO2, XH2O, SFV, QR
        
            LOGICAL, DIMENSION (:,:) :: SOLIDS
            
            DOUBLE PRECISION, DIMENSION (:) :: XUIN, YVIN
            
            DOUBLE PRECISION :: PRESSURE       
            
            INTEGER :: LL1, LL2, LL3,imax,jmax, m, l, k, i, j, imu, imu_start, imu_end, im, igas, num_band
            
            DOUBLE PRECISION :: wvnb, dwvnb, wvnbm, tz, tziv, xy1, xy2, wvnb1, wvnb2,factor, term_h2o 
            
            DOUBLE PRECISION :: T1, T2, T3, T4, T5, T6, T7, T8, term_co, term_co2, an, as, ae, aw, volu, vp, up, down 
            
            DOUBLE PRECISION, PARAMETER :: sigma = 5.6696e-8, fs = 0.5, fa = 0.5, epsilon = 1e-3
            
            DOUBLE PRECISION, PARAMETER :: planck=6.626176e-34, c0=2.99792458e8,boltzm=1.380662e-23

            DOUBLE PRECISION, PARAMETER :: c1=planck*c0*c0,c2=planck*c0/boltzm
            
            INTEGER, PARAMETER :: MDIR = 36, MMAX = 12, LMAX = 4, n_band = 10, ngauss = 4
            
            DOUBLE PRECISION, DIMENSION (:,:,:,:), ALLOCATABLE :: xip, xim
            
            DOUBLE PRECISION, DIMENSION (:,:,:), ALLOCATABLE :: xisouth, absor
            
            DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: xi_n, xi_s, xi_e, xi_w, xisum, xib, alpha
            
            DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xinorth, xieast, xiwest, tz0, tz1, tr0, tr1
            
            DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: ez0, ez1, er0, er1, xi0, xmu0, z, yr, rc
            
            INTEGER, DIMENSION (:), ALLOCATABLE :: lvalue
            
            !gauss quadrature point and weight
            DOUBLE PRECISION, DIMENSION (7) :: g, weight
            
            DOUBLE PRECISION, DIMENSION (9,4,9),SAVE :: b_co, b_co2, b_h2o
            DOUBLE PRECISION, DIMENSION (mmax,lmax), SAVE :: xmu, eta, xi, w

            LL1 = NUMP*L3+2
            LL2 = LL1-1
            LL3 = LL1-2
            imax = LL1
            jmax = M1

            ALLOCATE(xip(imax,jmax,mmax,lmax))
            ALLOCATE(xi_n(imax,jmax))
            ALLOCATE(xi_s(imax,jmax))
            ALLOCATE(xi_e(imax,jmax))
            ALLOCATE(xi_w(imax,jmax))
            ALLOCATE(xinorth(imax))
            ALLOCATE(xisouth(imax,mmax,lmax))
            ALLOCATE(xieast(jmax))
            ALLOCATE(xiwest(jmax))
            ALLOCATE(xim(imax,jmax,mmax,lmax+1))
            ALLOCATE(xisum(imax,jmax))
            ALLOCATE(xib(imax-1,jmax-1))
            ALLOCATE(absor(imax-1,jmax-1,ngauss))
            ALLOCATE(tz0(jmax-1))
            ALLOCATE(tz1(jmax-1))
            ALLOCATE(tr0(imax-1))
            ALLOCATE(tr1(imax-1))
            ALLOCATE(ez0(jmax-1))
            ALLOCATE(ez1(jmax-1))
            ALLOCATE(er0(imax-1))
            ALLOCATE(er1(imax-1))
            ALLOCATE(xi0(mmax))
            ALLOCATE(xmu0(mmax))
            ALLOCATE(z(imax))
            ALLOCATE(yr(jmax))
            ALLOCATE(rc(jmax))     
            ALLOCATE(alpha(mmax,lmax+1))
            ALLOCATE(lvalue(mmax))

            !c         r (mu)
            !c         |
            !c         |
            !c         |
            !c         |
            !c         |-------------------------------------------> z (xi)
            !c
            !c
            !c
            !c List of variables;
            !c    inp:   intensity at nodal point along m'th direction (i,j,m), (z,r,dir)
            !c    inm:   intensity at nodal point in the vincinity of the m'th direction (directional term)
            !c           used in the angular derivative term calculation
            !c    isum:  total intensity integrated over the solid angle, used as convergence criterion
            !c    in_n:  surface intensity at north surface of the local control volume
            !
            !c    in_s:  surface intensity at south surface
            !c    in_e:  surface intensity at east surface
            !c    in_w:  surface intensity at west surface
            !c    inorth: intensity at north boundary
            !c    isouth: intensity at south boundary (symmetry line), which is a function of direction
            !c    ieast:  intensity at east boundary
            !c    iwest:  intensity at west boundary
            !c    temp:  temperature field
            !c    ib:    blackbody intensity field
            !c    absor: field of absorption coefficient of gases
            !c    absorp: absorption coefficient of particle (soot)
            !c    xco:   mole fraction of CO
            !c    xco2:  mole fraction of CO2
            !c    xh2o:  mole fraction of H2O
            !c    sfv:   soot volume fraction
            !c    tz0,tz1,tr0,tr1: temperatures of boundaries
            !c    ez0,ez1,er0,er1: emissivities of boundaries
            !c    mu,eta,xi: direction cosines
            !c    w:         weight factor of each direction
            !c    f:         spatial weight factor, 0.5: central, 1, upwind
            !c    z,r,rc: spatial coordinates
            !c    an,as,ae,aw: surface areas of north, south, east, and west surfaces
            !c    v:           volume of the control volume
            !c    alpha: weight associated with the directional derivative term
            !c    press_abm: ambient pressure, in atm
            !c    n_band: number of narrow bands lumped into a wide band
            !c    ngauss: number of Gauss quadrature points used in band integration
            !c
            !c    
            !c  f = 0.5: central differencing scheme
            !c  f = 1.0: upwind differencing scheme
            !c
            !
            !c
            !c read in directional cosines and weights
            !c
            
            !initialize all variables
            xip = 0.d0
            xi_n= 0.d0
            xi_s = 0.d0
            xi_e = 0.d0
            xi_w = 0.d0
            xinorth =0.d0
            xisouth=0.d0
            xieast=0.d0
            xiwest=0.d0
            xim=0.d0
            xisum=0.d0
            xib=0.d0
            absor=0.d0
            tz0=0.d0
            tz1=0.d0
            tr0=0.d0
            tr1=0.d0
            ez0=0.d0
            ez1=0.d0
            er0=0.d0
            er1=0.d0
            xi0=0.d0
            xmu0=0.d0
            z=0.d0
            yr=0.d0
            rc=0.d0
            alpha=0.d0
            
            !c
            !c The following array is valid for T3 only
            !c
            lvalue(1)  = 2
            lvalue(2)  = 2
            lvalue(3)  = 4
            lvalue(4)  = 4
            lvalue(5)  = 2
            lvalue(6)  = 4
            lvalue(7)  = 4
            lvalue(8)  = 2
            lvalue(9)  = 4
            lvalue(10) = 4
            lvalue(11) = 2
            lvalue(12) = 2

            !read in data from t3_2d.dat if this is the first iteration. Since the save flag is used in their declarations, xmu, eta, xi, w will retain their values for subsequent iterations
            if((iter-1)==first) then
                xmu=0.d0
                eta=0.d0
                xi=0.d0
                w=0.d0
                open ( 5, file='t3_2d.dat')
                    do k = 1, mdir
                       read(5,*) m,l, xmu(m,l),eta(m,l),xi(m,l),w(m,l)
                    enddo
                close(5)
                
                !c
                !c modify the weight factor read in from t3_2d.dat
                !c
                do m = 1, mmax
                   do l = 1, lvalue(m)
                      w(m,l) = 2.0*pi*w(m,l)
                   enddo
                enddo
            endif         



            !c
            !c calculate alpha array
            !c  Note: alpha(m,1/2) is denoted by alpha(m,1)
            !c        alpha(m,3/2) is denoted by alpha(m,2)
            !c  in general alpha(m,l-1/2) by alpha(m,l)
            !c             alpha(m,l+1/2) by alpha(m,l+1)
            !c
                    do m = 1, mmax
                       alpha(m,1) = 0.0
                       do l = 2, lvalue(m)+1
                         alpha(m,l) = alpha(m,l-1)+w(m,l-1)*xmu(m,l-1)
                       enddo
                    enddo
            !c
            !c Gauss quadratures
            !c
                    if ( ngauss.eq.2 ) then
                       g(1) = 0.3399810436d0
                       g(2) = 0.8611363116d0
                       weight(1) = 0.6521451549d0
                       weight(2) = 0.3478548451d0
                    else if ( ngauss.eq.3 ) then
                       g(1) = 0.2386191861d0
                       g(2) = 0.6612093865d0
                       g(3) = 0.9324695142d0
                       weight(1) = 0.4679139346d0
                       weight(2) = 0.3607615731d0
                       weight(3) = 0.1713244924d0
                    else if ( ngauss.eq.4 ) then
                       g(1) = 0.18343d0
                       g(2) = 0.52553d0
                       g(3) = 0.79667d0
                       g(4) = 0.96029d0
                       weight(1) = 0.36268d0
                       weight(2) = 0.31371d0
                       weight(3) = 0.22238d0
                       weight(4) = 0.10123d0
                    endif
            !c
            !c grid information
            !c a factor of 100 is used to convert cm to m
            !c
                  do i = 2, LL1
                     z(i) = XUIN(I)/100.0d0
                  enddo

                  do j = 2, M1
                     yr(j) = YVIN(J)/100.0d0
                  enddo

                  do j = 2, M2
                     rc(j) = 0.5d0*(YVIN(J)+YVIN(J+1))/100.0d0
                  enddo
            !c
            !c specifiy boundary properties
            !c
                    do j = 2, M2
                       tz0(j) = 0.1
                       ez0(j) = 1.0
                       tz1(j) = 0.1
                       ez1(j) = 1.0
                    enddo

                    do i = 2, LL2
                       tr0(i) = 0.1
                       er0(i) = 1.0
                       tr1(i) = 0.1
                       er1(i) = 1.0
                    enddo
            !c
            !c
            !c  specify gas properties
            !c

            !c
            !c  wave_number covers 150 - 9300 /1cm
            !c  uniform interval of wave_number at 25 1/cm
            !c
                  dwvnb = 25.0d0
            !c
            !c clear the radiation source term
            !c
                    QR = 0.0d0

                !if on the first iteration, populate b_co, b_co2, b_h2o. Since they have the save flag in their declarations, their values will remain unchanged for remaining iterations
                if((iter-1)==first) then
                
                    !first, initialize all three arrays to zero
                    b_co=0.d0
                    b_co2=0.d0
                    b_h2o=0.d0

                    !c 
                    !c===========================================================
                    !c   coefficients of H2O per mole fraction
                    !c===========================================================
                    !c
                    !c    region 1, point 1
                    !c
                    b_h2o(1,1,1) = -4.6226684975
                    b_h2o(2,1,1) = 0.0606537855
                    b_h2o(3,1,1) = -1.4055711114e-4
                    b_h2o(4,1,1) = 1.670262129e-7
                    b_h2o(5,1,1) = -1.1287316401e-10
                    b_h2o(6,1,1) = 4.3845097989e-14
                    b_h2o(7,1,1) = -9.1365821118e-18
                    b_h2o(8,1,1) = 7.9370555001e-22
                    b_h2o(9,1,1) = -5.8460485685e-28

                    !c     point 2

                    b_h2o(1,2,1) = 23.7540624693
                    b_h2o(2,2,1) = 0.1832848169
                    b_h2o(3,2,1) = -7.7152778965e-4
                    b_h2o(4,2,1) = 1.3788292421e-6
                    b_h2o(5,2,1) = -1.3763632419e-9
                    b_h2o(6,2,1) = 8.2203634721e-13
                    b_h2o(7,2,1) = -2.9234477414e-16
                    b_h2o(8,2,1) = 5.7163028247e-20
                    b_h2o(9,2,1) = -4.7335078595e-24

                    !c      point 3

                    b_h2o(1,3,1) = 179.017279255
                    b_h2o(2,3,1) = 0.6717805387
                    b_h2o(3,3,1) = -3.6195840549e-3
                    b_h2o(4,3,1) = 7.0186792124e-6
                    b_h2o(5,3,1) = -7.381442325e-9
                    b_h2o(6,3,1) = 4.5835457076e-12
                    b_h2o(7,3,1) = -1.6797298393e-15
                    b_h2o(8,3,1) = 3.3618996568e-19
                    b_h2o(9,3,1) = -2.8347921708e-23

                    !c      point 4

                    b_h2o(1,4,1) = 1793.9386779014
                    b_h2o(2,4,1) = 3.0740545097
                    b_h2o(3,4,1) = -0.0220013704
                    b_h2o(4,4,1) = 4.3042466222e-5
                    b_h2o(5,4,1) = -4.4705247608e-8
                    b_h2o(6,4,1) = 2.740887903e-11
                    b_h2o(7,4,1) = -9.9399228138e-15
                    b_h2o(8,4,1) = 1.9725413772e-18
                    b_h2o(9,4,1) = -1.6514099979e-22


                    !      H2O, region 2, point 1

                    b_h2o(1,1,2) = -0.2139832301
                    b_h2o(2,1,2) = 2.0222497415e-3
                    b_h2o(3,1,2) = -7.1614200923e-6
                    b_h2o(4,1,2) = 1.3247619079e-8
                    b_h2o(5,1,2) = -1.417488264e-11
                    b_h2o(6,1,2) = 9.1804318892e-15
                    b_h2o(7,1,2) = -3.4889064957e-18
                    b_h2o(8,1,2) = 7.1469744082e-22
                    b_h2o(9,1,2) = -6.1012215375e-26


                    !c     point 2

                    b_h2o(1,2,2) = -0.3768703374
                    b_h2o(2,2,2) = 3.1705098205e-3
                    b_h2o(3,2,2) = -9.588787036e-6
                    b_h2o(4,2,2) = 1.5084287478e-8
                    b_h2o(5,2,2) = -1.2383768309e-11
                    b_h2o(6,2,2) = 6.4035817327e-15
                    b_h2o(7,2,2) = -2.1030537498e-18
                    b_h2o(8,2,2) = 3.9770603672e-22
                    b_h2o(9,2,2) = -3.2857753346e-26

                    !c      point 3

                    b_h2o(1,3,2) = -0.6360312707
                    b_h2o(2,3,2) = 6.2268798586e-3
                    b_h2o(3,3,2) = -2.3716625641e-5
                    b_h2o(4,3,2) = 5.9960355826e-8
                    b_h2o(5,3,2) = -6.7525447063e-11
                    b_h2o(6,3,2) = 4.3066100424e-14
                    b_h2o(7,3,2) = -1.6182024924e-17
                    b_h2o(8,3,2) = 3.3336988518e-21
                    b_h2o(9,3,2) = -2.9025269295e-25

                    !c     point 4

                    b_h2o(1,4,2) = -4.3928560588
                    b_h2o(2,4,2) = 0.032758876
                    b_h2o(3,4,2) = -1.1402203873e-4
                    b_h2o(4,4,2) = 4.877666771e-7
                    b_h2o(5,4,2) = -6.6795617143e-10
                    b_h2o(6,4,2) = 4.5104314241e-13
                    b_h2o(7,4,2) = -1.6757754798e-16
                    b_h2o(8,4,2) = 3.3061680816e-20
                    b_h2o(9,4,2) = -2.7212493277e-24


                    !c H2O, region 3, point 1

                    b_h2o(1,1,3) = -1.3694166853
                    b_h2o(2,1,3) = 0.0130112743
                    b_h2o(3,1,3) = -3.8316807173e-5
                    b_h2o(4,1,3) = 6.4310398459e-8
                    b_h2o(5,1,3) = -6.4905846683e-11
                    b_h2o(6,1,3) = 4.0108332004e-14
                    b_h2o(7,1,3) = -1.4841069591e-17
                    b_h2o(8,1,3) = 3.0153453839e-21
                    b_h2o(9,1,3) = -2.5843568131e-25

                    !c     point 2

                    b_h2o(1,2,3) = 6.3823015006
                    b_h2o(2,2,3) = 5.8169621735e-3
                    b_h2o(3,2,3) = -6.7526622167e-5
                    b_h2o(4,2,3) = 1.6472312764e-7
                    b_h2o(5,2,3) = -1.9937103983e-10
                    b_h2o(6,2,3) = 1.3692720701e-13
                    b_h2o(7,2,3) = -5.4228698372e-17
                    b_h2o(8,2,3) = 1.1546456474e-20
                    b_h2o(9,2,3) = -1.0238629348e-24

                    !c     point 3

                    b_h2o(1,3,3) = 56.9074867128
                    b_h2o(2,3,3) = -0.1174976228
                    b_h2o(3,3,3) = 2.5105474332e-5
                    b_h2o(4,3,3) = 2.8970156489e-7
                    b_h2o(5,3,3) = -5.0936360229e-10
                    b_h2o(6,3,3) = 4.0524742695e-13
                    b_h2o(7,3,3) = -1.7336510456e-16
                    b_h2o(8,3,3) = 3.8662570041e-20
                    b_h2o(9,3,3) = -3.5333984329e-24

                    !c     point 4

                    b_h2o(1,4,3) = 371.363881695
                    b_h2o(2,4,3) = -0.5630745243
                    b_h2o(3,4,3) = -7.819350497e-4
                    b_h2o(4,4,3) = 3.5684376526e-6
                    b_h2o(5,4,3) = -4.9653182505e-9
                    b_h2o(6,4,3) = 3.6043931172e-12
                    b_h2o(7,4,3) = -1.4651463356e-15
                    b_h2o(8,4,3) = 3.1600169962e-19
                    b_h2o(9,4,3) = -2.8195711554e-23


                    !c H2O, region 4, point 1

                    b_h2o(1,1,4) = -0.032493021
                    b_h2o(2,1,4) = 3.0948642323e-4
                    b_h2o(3,1,4) = -1.1019211949e-6
                    b_h2o(4,1,4) = 2.0715096254e-9
                    b_h2o(5,1,4) = -2.2090700502e-12
                    b_h2o(6,1,4) = 1.40833971e-15
                    b_h2o(7,1,4) = -5.3148749113e-19
                    b_h2o(8,1,4) = 1.0933648743e-22
                    b_h2o(9,1,4) = -9.4377961609e-27

                    !c     point 2

                    b_h2o(1,2,4) = -0.1055071796
                    b_h2o(2,2,4) = 9.0110624648e-4
                    b_h2o(3,2,4) = -3.0888796513e-6
                    b_h2o(4,2,4) = 6.3895785498e-9
                    b_h2o(5,2,4) = -7.4130435868e-12
                    b_h2o(6,2,4) = 5.0171309314e-15
                    b_h2o(7,2,4) = -1.972129223e-18
                    b_h2o(8,2,4) = 4.172318641e-22
                    b_h2o(9,2,4) = -3.6739854121e-26

                    !c     point 3

                    b_h2o(1,3,4) = 0.1117139498
                    b_h2o(2,3,4) = -9.7409898269e-4
                    b_h2o(3,3,4) = 3.7228936396e-6
                    b_h2o(4,3,4) = -3.1097720621e-9
                    b_h2o(5,3,4) = -1.8058802966e-13
                    b_h2o(6,3,4) = 1.7481372382e-15
                    b_h2o(7,3,4) = -1.0946561022e-18
                    b_h2o(8,3,4) = 2.8747352089e-22
                    b_h2o(9,3,4) = -2.8562368635e-26

                    !c     point 4

                    b_h2o(1,4,4) = 1.7136447382
                    b_h2o(2,4,4) = -0.0137557425
                    b_h2o(3,4,4) = 6.6726195471e-5
                    b_h2o(4,4,4) = -1.1798448605e-7
                    b_h2o(5,4,4) = 1.1235098815e-10
                    b_h2o(6,4,4) = -6.3711354759e-14
                    b_h2o(7,4,4) = 2.1614010003e-17
                    b_h2o(8,4,4) = -4.0601843806e-21
                    b_h2o(9,4,4) = 3.2521343556e-25


                    !c H2O, region 5, point 1

                    b_h2o(1,1,5) = -4.0869031122e-4
                    b_h2o(2,1,5) = 5.3760507969e-6
                    b_h2o(3,1,5) = -1.0321069622e-8
                    b_h2o(4,1,5) = 1.0849317876e-11
                    b_h2o(5,1,5) = -1.699681071e-14
                    b_h2o(6,1,5) = 2.4445499249e-17
                    b_h2o(7,1,5) = -1.4759500685e-20
                    b_h2o(8,1,5) = 3.8726795055e-24
                    b_h2o(9,1,5) = -3.7378115213e-28

                    !c     point 2

                    b_h2o(1,2,5) = 0.0220620354
                    b_h2o(2,2,5) = -1.5129712188e-4
                    b_h2o(3,2,5) = 3.6759020939e-7
                    b_h2o(4,2,5) = -1.8134798694e-10
                    b_h2o(5,2,5) = -1.6242534755e-13
                    b_h2o(6,2,5) = 3.1977130067e-16
                    b_h2o(7,2,5) = -1.9760846976e-19
                    b_h2o(8,2,5) = 5.4336715977e-23
                    b_h2o(9,2,5) = -5.6560991735e-27

                    !c      point 3

                    b_h2o(1,3,5) = 0.067013237
                    b_h2o(2,3,5) = 6.9479079513e-4
                    b_h2o(3,3,5) = -3.7148655906e-6
                    b_h2o(4,3,5) = 8.5662030814e-9
                    b_h2o(5,3,5) = -9.6827874048e-12
                    b_h2o(6,3,5) = 6.335176989e-15
                    b_h2o(7,3,5) = -2.4302727265e-18
                    b_h2o(8,3,5) = 5.054916212e-22
                    b_h2o(9,3,5) = -4.3954489834e-26

                    !c     point 4

                    b_h2o(1,4,5) = 1.1265129652
                    b_h2o(2,4,5) = 2.4817789868e-3
                    b_h2o(3,4,5) = -2.1299604694e-5
                    b_h2o(4,4,5) = 5.1003415296e-8
                    b_h2o(5,4,5) = -5.8035545453e-11
                    b_h2o(6,4,5) = 3.7363355372e-14
                    b_h2o(7,4,5) = -1.3943007133e-17
                    b_h2o(8,4,5) = 2.8130937975e-21
                    b_h2o(9,4,5) = -2.3757572315e-25


                    !c H2O, region 6, point 1

                    b_h2o(1,1,6) = -0.4018728038
                    b_h2o(2,1,6) = 3.9720880474e-3
                    b_h2o(3,1,6) = -9.5312920851e-6
                    b_h2o(4,1,6) = 1.498107137e-8
                    b_h2o(5,1,6) = -1.483203175e-11
                    b_h2o(6,1,6) = 9.1924166601e-15
                    b_h2o(7,1,6) = -3.4482973673e-18
                    b_h2o(8,1,6) = 7.1256797197e-22
                    b_h2o(9,1,6) = -6.2066142973e-26

                    !c     point 2

                    b_h2o(1,2,6) = 1.7139016993
                    b_h2o(2,2,6) = 0.0201012362
                    b_h2o(3,2,6) = -9.0617805539e-5
                    b_h2o(4,2,6) = 1.7886230232e-7
                    b_h2o(5,2,6) = -1.9424926819e-10
                    b_h2o(6,2,6) = 1.2446181669e-13
                    b_h2o(7,2,6) = -4.6922516867e-17
                    b_h2o(8,2,6) = 9.6268909749e-21
                    b_h2o(9,2,6) = -8.291538841e-25

                    !c     point 3

                    b_h2o(1,3,6) = 29.4022178839
                    b_h2o(2,3,6) = -0.0286385869
                    b_h2o(3,3,6) = -1.2126350698e-4
                    b_h2o(4,3,6) = 4.0074426767e-7
                    b_h2o(5,3,6) = -5.2203650788e-10
                    b_h2o(6,3,6) = 3.6724976602e-13
                    b_h2o(7,3,6) = -1.4642731477e-16
                    b_h2o(8,3,6) = 3.1164611651e-20
                    b_h2o(9,3,6) = -2.7538176881e-24

                    !c     point 4

                    b_h2o(1,4,6) = 247.5812567309
                    b_h2o(2,4,6) = -0.4360339762
                    b_h2o(3,4,6) = -2.8477986693e-4
                    b_h2o(4,4,6) = 1.9661978824e-6
                    b_h2o(5,4,6) = -2.8983000371e-9
                    b_h2o(6,4,6) = 2.1536623626e-12
                    b_h2o(7,4,6) = -8.8631165063e-16
                    b_h2o(8,4,6) = 1.9256311931e-19
                    b_h2o(9,4,6) = -1.7261174866e-23

                    !c H2O, region 7, point 1

                    b_h2o(1,1,7) = -2.9005901156e-4
                    b_h2o(2,1,7) = 2.7107223474e-5
                    b_h2o(3,1,7) = -1.2006657102e-7
                    b_h2o(4,1,7) = 2.3668105063e-10
                    b_h2o(5,1,7) = -2.5868402036e-13
                    b_h2o(6,1,7) = 1.6765489461e-16
                    b_h2o(7,1,7) = -6.3864313278e-20
                    b_h2o(8,1,7) = 1.3214437873e-23
                    b_h2o(9,1,7) = -1.146561603e-27

                    !c     point 2

                    b_h2o(1,2,7) = -0.030543143
                    b_h2o(2,2,7) = 3.6523551698e-4
                    b_h2o(3,2,7) = -1.3182718088e-6
                    b_h2o(4,2,7) = 2.5816534441e-9
                    b_h2o(5,2,7) = -2.9223399572e-12
                    b_h2o(6,2,7) = 1.9620159847e-15
                    b_h2o(7,2,7) = -7.6685442383e-19
                    b_h2o(8,2,7) = 1.6107220101e-22
                    b_h2o(9,2,7) = -1.4058459839e-26

                    !c     point 3

                    b_h2o(1,3,7) = -1.4883428457
                    b_h2o(2,3,7) = 0.0136517399
                    b_h2o(3,3,7) = -3.085290117e-5
                    b_h2o(4,3,7) = 3.7187253968e-8
                    b_h2o(5,3,7) = -2.6520534521e-11
                    b_h2o(6,3,7) = 1.1602019298e-14
                    b_h2o(7,3,7) = -3.0621304511e-18
                    b_h2o(8,3,7) = 4.4645887491e-22
                    b_h2o(9,3,7) = -2.7454693884e-26

                    !c     point 4

                    b_h2o(1,4,7) = 86.1946344823
                    b_h2o(2,4,7) = -0.2429700133
                    b_h2o(3,4,7) = 2.839801217e-4
                    b_h2o(4,4,7) = -3.7029549114e-8
                    b_h2o(5,4,7) = -2.4508267028e-10
                    b_h2o(6,4,7) = 2.6811105923e-13
                    b_h2o(7,4,7) = -1.2927012361e-16
                    b_h2o(8,4,7) = 3.0672681203e-20
                    b_h2o(9,4,7) = -2.9102961297e-24

                    !c H2O, region 8, point 1

                    b_h2o(1,1,8) = 0.0105506187
                    b_h2o(2,1,8) = -7.4985344619e-5
                    b_h2o(3,1,8) = 2.2120071526e-7
                    b_h2o(4,1,8) = -2.8806416557e-10
                    b_h2o(5,1,8) = 1.2131376962e-13
                    b_h2o(6,1,8) = 8.209696771e-17
                    b_h2o(7,1,8) = -9.3757235236e-20
                    b_h2o(8,1,8) = 3.1020843167e-23
                    b_h2o(9,1,8) = -3.5408267454e-27

                    !c     point 2

                    b_h2o(1,2,8) = 0.0235643653
                    b_h2o(2,2,8) = -1.0973346008e-4
                    b_h2o(3,2,8) = 8.7022911858e-8
                    b_h2o(4,2,8) = 6.5793834719e-10
                    b_h2o(5,2,8) = -1.3394394396e-12
                    b_h2o(6,2,8) = 1.1854831364e-15
                    b_h2o(7,2,8) = -5.4934647632e-19
                    b_h2o(8,2,8) = 1.2966423584e-22
                    b_h2o(9,2,8) = -1.2321362809e-26

                    !c     point 3

                    b_h2o(1,3,8) = -0.3787032087
                    b_h2o(2,3,8) = 5.2769046554e-3
                    b_h2o(3,3,8) = -1.6384060506e-5
                    b_h2o(4,3,8) = 2.7088844129e-8
                    b_h2o(5,3,8) = -2.6196406412e-11
                    b_h2o(6,3,8) = 1.543422985e-14
                    b_h2o(7,3,8) = -5.4611288123e-18
                    b_h2o(8,3,8) = 1.0664977305e-21
                    b_h2o(9,3,8) = -8.832853432e-26

                    !c     point 4

                    b_h2o(1,4,8) = 7.3435971122
                    b_h2o(2,4,8) = -0.0166185875
                    b_h2o(3,4,8) = 2.7948182365e-5
                    b_h2o(4,4,8) = -3.1159565604e-8
                    b_h2o(5,4,8) = 2.3081316536e-11
                    b_h2o(6,4,8) = -1.130808688e-14
                    b_h2o(7,4,8) = 3.5275198127e-18
                    b_h2o(8,4,8) = -6.3388553938e-22
                    b_h2o(9,4,8) = 4.980644459e-26

                    !c H2O, region 9, point 1

                    b_h2o(1,1,9) = -7.6177742872e-4
                    b_h2o(2,1,9) = 1.0291489468e-5
                    b_h2o(3,1,9) = -3.7992793345e-8
                    b_h2o(4,1,9) = 7.1063323458e-11
                    b_h2o(5,1,9) = -7.6934537724e-14
                    b_h2o(6,1,9) = 5.0786491704e-17
                    b_h2o(7,1,9) = -1.9782403024e-20
                    b_h2o(8,1,9) = 4.1622146992e-24
                    b_h2o(9,1,9) = -3.6475203331e-28

                    !c     point 2

                    b_h2o(1,2,9) = -1.3723603491e-3
                    b_h2o(2,2,9) = 7.9641250619e-5
                    b_h2o(3,2,9) = -2.548471554e-7
                    b_h2o(4,2,9) = 4.4041930514e-10
                    b_h2o(5,2,9) = -4.4288093334e-13
                    b_h2o(6,2,9) = 2.7556075812e-16
                    b_h2o(7,2,9) = -1.0334716682e-19
                    b_h2o(8,2,9) = 2.1280645608e-23
                    b_h2o(9,2,9) = -1.8434142083e-27

                    !c     point 3

                    b_h2o(1,3,9) = 0.1833223942
                    b_h2o(2,3,9) = 1.0096514972e-3
                    b_h2o(3,3,9) = -5.2524009192e-6
                    b_h2o(4,3,9) = 1.0892479176e-8
                    b_h2o(5,3,9) = -1.2172569421e-11
                    b_h2o(6,3,9) = 7.9611584861e-15
                    b_h2o(7,3,9) = -3.0469941324e-18
                    b_h2o(8,3,9) = 6.3206644733e-22
                    b_h2o(9,3,9) = -5.4881622176e-26

                    !c     point 4

                    b_h2o(1,4,9) = 7.5419183797
                    b_h2o(2,4,9) = -0.0110206448
                    b_h2o(3,4,9) = -2.0889422455e-5
                    b_h2o(4,4,9) = 8.7633317845e-8
                    b_h2o(5,4,9) = -1.2107102788e-10
                    b_h2o(6,4,9) = 8.7778932146e-14
                    b_h2o(7,4,9) = -3.5671890686e-17
                    b_h2o(8,4,9) = 7.692720215e-21
                    b_h2o(9,4,9) = -6.8628057212e-25
                    !c
                    !c======================================================
                    !c CO2, region 2, point 1
                    !c  First clear b_co and b_co2 arrays
                    !c======================================================
                    !c

                    b_co2(1,1,2) = -1.1417440756e-3
                    b_co2(2,1,2) = 2.5427475395e-5
                    b_co2(3,1,2) = -1.7574870474e-7
                    b_co2(4,1,2) = 5.0317251839e-10
                    b_co2(5,1,2) = -6.1774564504e-13
                    b_co2(6,1,2) = 3.9978592136e-16
                    b_co2(7,1,2) = -1.4518436803e-19
                    b_co2(8,1,2) = 2.8220967974e-23
                    b_co2(9,1,2) = -2.2987649877e-27

                    !c     point 2

                    b_co2(1,2,2) = 0.0869306415
                    b_co2(2,2,2) = -5.3122138969e-4
                    b_co2(3,2,2) = -3.2738066305e-7
                    b_co2(4,2,2) = 7.2840045627e-9
                    b_co2(5,2,2) = -1.3486894459e-11
                    b_co2(6,2,2) = 1.1164132409e-14
                    b_co2(7,2,2) = -4.865310589e-18
                    b_co2(8,2,2) = 1.0912391303e-21
                    b_co2(9,2,2) = -9.9590226566e-26

                    !c     point 3

                    b_co2(1,3,2) = 2.059085323
                    b_co2(2,3,2) = -0.0243824223
                    b_co2(3,3,2) = 1.0102794982e-4
                    b_co2(4,3,2) = -1.2983188665e-7
                    b_co2(5,3,2) = 7.0862170439e-11
                    b_co2(6,3,2) = -1.0788373478e-14
                    b_co2(7,3,2) = -5.5074485928e-18
                    b_co2(8,3,2) = 2.5597438762e-21
                    b_co2(9,3,2) = -3.1086026704e-25

                    !c     point 4

                    b_co2(1,4,2) = 194.8751431764
                    b_co2(2,4,2) = -0.299131129
                    b_co2(3,4,2) = 1.2507567027e-5
                    b_co2(4,4,2) = 7.802841776e-7
                    b_co2(5,4,2) = -1.3862058992e-9
                    b_co2(6,4,2) = 1.1494688181e-12
                    b_co2(7,4,2) = -5.1146271007e-16
                    b_co2(8,4,2) = 1.1783422818e-19
                    b_co2(9,4,2) = -1.1053381452e-23

                    !c CO2, region 4, point 1

                    b_co2(1,1,4) = -0.0791812871
                    b_co2(2,1,4) = 7.0169424086e-4
                    b_co2(3,1,4) = -2.498794983e-6
                    b_co2(4,1,4) = 4.851937988e-9
                    b_co2(5,1,4) = -5.4740272761e-12
                    b_co2(6,1,4) = 3.6712502598e-15
                    b_co2(7,1,4) = -1.443376703e-18
                    b_co2(8,1,4) = 3.0854227768e-22
                    b_co2(9,1,4) = -2.7694545639e-26

                    !c     point 2

                    b_co2(1,2,4) = -15.704628314
                    b_co2(2,2,4) = 0.1550552055
                    b_co2(3,2,4) = -6.1590360598e-4
                    b_co2(4,2,4) = 1.2949975966e-6
                    b_co2(5,2,4) = -1.5778934756e-9
                    b_co2(6,2,4) = 1.1407685232e-12
                    b_co2(7,2,4) = -4.7653532286e-16
                    b_co2(8,2,4) = 1.058098031e-19
                    b_co2(9,2,4) = -9.6544319952e-24

                    !c     point 3

                    b_co2(1,3,4) = 298.8728956676
                    b_co2(2,3,4) = -2.7078089662
                    b_co2(3,3,4) = 0.01018576
                    b_co2(4,3,4) = -1.7343136089e-5
                    b_co2(5,3,4) = 1.6665382858e-8
                    b_co2(6,3,4) = -9.6816932218e-12
                    b_co2(7,3,4) = 3.3781765575e-15
                    b_co2(8,3,4) = -6.521160649e-19
                    b_co2(9,3,4) = 5.3526552942e-23

                    !c     point 4

                    b_co2(1,4,4) = 4465.1602778195
                    b_co2(2,4,4) = 1.9992429275
                    b_co2(3,4,4) = -0.0462310369
                    b_co2(4,4,4) = 1.1342303425e-4
                    b_co2(5,4,4) = -1.3624822029e-7
                    b_co2(6,4,4) = 9.2968821849e-11
                    b_co2(7,4,4) = -3.6608481183e-14
                    b_co2(8,4,4) = 7.7524228749e-18
                    b_co2(9,4,4) = -6.8383555707e-22


                    !c CO2, region 6, point 1

                    b_co2(1,1,6) = 0.0218195207
                    b_co2(2,1,6) = -2.4604394691e-5
                    b_co2(3,1,6) = -8.2235434821e-7
                    b_co2(4,1,6) = 3.9540060661e-9
                    b_co2(5,1,6) = -7.6142401035e-12
                    b_co2(6,1,6) = 7.4184074934e-15
                    b_co2(7,1,6) = -3.7640890577e-18
                    b_co2(8,1,6) = 9.5515184916e-22
                    b_co2(9,1,6) = -9.616187378e-26

                    !c     point 2

                    b_co2(1,2,6) = 2.0642210602
                    b_co2(2,2,6) = -0.018545259
                    b_co2(3,2,6) = 6.1941965699e-5
                    b_co2(4,2,6) = -8.7184431636e-8
                    b_co2(5,2,6) = 7.014483915e-11
                    b_co2(6,2,6) = -3.4914803762e-14
                    b_co2(7,2,6) = 1.0687422434e-17
                    b_co2(8,2,6) = -1.849294783e-21
                    b_co2(9,2,6) = 1.3851364232e-25

                    !c      point 3

                    b_co2(1,3,6) = 6.4848181434
                    b_co2(2,3,6) = 0.0452302182
                    b_co2(3,3,6) = -1.9144964213e-4
                    b_co2(4,3,6) = 3.8357330514e-7
                    b_co2(5,3,6) = -4.3842429076e-10
                    b_co2(6,3,6) = 2.9688860615e-13
                    b_co2(7,3,6) = -1.1764475143e-16
                    b_co2(8,3,6) = 2.5196729562e-20
                    b_co2(9,3,6) = -2.2520578581e-24

                    !c     point 4

                    b_co2(1,4,6) = 187.1519936448
                    b_co2(2,4,6) = -0.4648190631
                    b_co2(3,4,6) = 2.5735851933e-4
                    b_co2(4,4,6) = 6.9980414378e-7
                    b_co2(5,4,6) = -1.4539411058e-9
                    b_co2(6,4,6) = 1.2209512984e-12
                    b_co2(7,4,6) = -5.3885262716e-16
                    b_co2(8,4,6) = 1.2284842442e-19
                    b_co2(9,4,6) = -1.1420594233e-23


                    !c CO2, region 8, point 1

                    b_co2(1,1,8) = 4.0395041622e-3
                    b_co2(2,1,8) = -3.3724414209e-5
                    b_co2(3,1,8) = 7.3378764261e-8
                    b_co2(4,1,8) = 3.720278744e-11
                    b_co2(5,1,8) = -2.0249982863e-13
                    b_co2(6,1,8) = 2.0617692868e-16
                    b_co2(7,1,8) = -9.8187212252e-20
                    b_co2(8,1,8) = 2.3104967224e-23
                    b_co2(9,1,8) = -2.1711884818e-27

                    !c     point 2

                    b_co2(1,2,8) = -0.0335907794
                    b_co2(2,2,8) = 3.3247080584e-4
                    b_co2(3,2,8) = -5.9961914783e-7
                    b_co2(4,2,8) = 7.7943841922e-10
                    b_co2(5,2,8) = -8.4086298153e-13
                    b_co2(6,2,8) = 6.1444292714e-16
                    b_co2(7,2,8) = -2.6575769361e-19
                    b_co2(8,2,8) = 6.107962336e-23
                    b_co2(9,2,8) = -5.7516843534e-27

                    !c     point 3

                    b_co2(1,3,8) = 0.2993174259
                    b_co2(2,3,8) = -8.6267844768e-5
                    b_co2(3,3,8) = -1.4207228799e-6
                    b_co2(4,3,8) = 4.2669891884e-9
                    b_co2(5,3,8) = -5.8710574935e-12
                    b_co2(6,3,8) = 4.412205421e-15
                    b_co2(7,3,8) = -1.8616611021e-18
                    b_co2(8,3,8) = 4.1461384952e-22
                    b_co2(9,3,8) = -3.79779695e-26

                    !c     point 4


                    b_co2(1,4,8) = 2.7739806496
                    b_co2(2,4,8) = -3.740507148e-3
                    b_co2(3,4,8) = -1.0452141509e-5
                    b_co2(4,4,8) = 3.8920338391e-8
                    b_co2(5,4,8) = -5.2956487048e-11
                    b_co2(6,4,8) = 3.8330441347e-14
                    b_co2(7,4,8) = -1.5600127459e-17
                    b_co2(8,4,8) = 3.3714221405e-21
                    b_co2(9,4,8) = -3.013982123e-25
                    !c
                    !c============================================================
                    !c CO, region 3
                    !c
                    !c   Index order: polynomail order index, 
                    !c                quadrature point index, and 
                    !c                spectral band index
                    !c============================================================
                    !c     point 3

                    b_co(1,3,3) = 0.0287136627
                    b_co(2,3,3) = -2.5869700823e-4
                    b_co(3,3,3) = 9.193456307e-7
                    b_co(4,3,3) = -1.6983259108e-9
                    b_co(5,3,3) = 1.8034686423e-12
                    b_co(6,3,3) = -1.1397078395e-15
                    b_co(7,3,3) = 4.2240030695e-19
                    b_co(8,3,3) = -8.4376299604e-23
                    b_co(9,3,3) = 6.98428958e-27

                    !c      point 4


                    b_co(1,4,3) = 0.7090123163
                    b_co(2,4,3) = -6.4309780009e-3
                    b_co(3,4,3) = 2.2996119549e-5
                    b_co(4,4,3) = -4.2617246869e-8
                    b_co(5,4,3) = 4.5117666965e-11
                    b_co(6,4,3) = -2.819592878e-14
                    b_co(7,4,3) = 1.0299017613e-17
                    b_co(8,4,3) = -2.028578694e-21
                    b_co(9,4,3) = 1.6595276939e-25

                    !c CO, region 4, point 1

                    !c     point 2

                    b_co(1,2,4) = 1.1355743124
                    b_co(2,2,4) = -0.010786347
                    b_co(3,2,4) = 3.7343454033e-5
                    b_co(4,2,4) = -6.021765595e-8
                    b_co(5,2,4) = 5.4611082233e-11
                    b_co(6,2,4) = -2.9562741932e-14
                    b_co(7,2,4) = 9.4839617225e-18
                    b_co(8,2,4) = -1.6638293793e-21
                    b_co(9,2,4) = 1.2295222733e-25

                    !c     point 3

                    b_co(1,3,4) = 15.8056189052
                    b_co(2,3,4) = -0.031949975
                    b_co(3,3,4) = 2.4866643413e-5
                    b_co(4,3,4) = 1.243610046e-8
                    b_co(5,3,4) = -3.6183397134e-11
                    b_co(6,3,4) = 2.7453977922e-14
                    b_co(7,3,4) = -1.033187206e-17
                    b_co(8,3,4) = 1.9749226988e-21
                    b_co(9,3,4) = -1.5316125754e-25

                    !c     point 4

                    b_co(1,4,4) = 438.4111286318
                    b_co(2,4,4) = -1.3532216857
                    b_co(3,4,4) = 2.0671699692e-3
                    b_co(4,4,4) = -1.6609328455e-6
                    b_co(5,4,4) = 6.516728412e-10
                    b_co(6,4,4) = -4.5252345995e-14
                    b_co(7,4,4) = -5.7247481593e-17
                    b_co(8,4,4) = 1.9724587719e-20
                    b_co(9,4,4) = -2.0281788188e-24
                    !c
                    !c CO region 6, CO only has two weak narrow bands at 3775 and 3800 cm-1
                    !c    The contribution from CO can be neglected.
                    !c
                    !c CO , region 7, point 1
                    !c

                    !c     point 2

                    b_co(1,2,7) = 1.4668878362e-3
                    b_co(2,2,7) = -1.322103e-5
                    b_co(3,2,7) = 4.6687806789e-8
                    b_co(4,2,7) = -8.4560817156e-11
                    b_co(5,2,7) = 8.639137025e-14
                    b_co(6,2,7) = -5.1522051397e-17
                    b_co(7,2,7) = 1.7884245925e-20
                    b_co(8,2,7) = -3.3551838374e-24
                    b_co(9,2,7) = 2.629704909e-28

                    !c     point 3

                    b_co(1,3,7) = -0.0149755996
                    b_co(2,3,7) = 1.1938235925e-4
                    b_co(3,3,7) = -2.763202796e-7
                    b_co(4,3,7) = 3.4123943192e-10
                    b_co(5,3,7) = -2.3624374726e-13
                    b_co(6,3,7) = 9.2525116665e-17
                    b_co(7,3,7) = -1.8965230056e-20
                    b_co(8,3,7) = 1.5173058552e-24
                    b_co(9,3,7) = 1.6035385572e-29

                    !c      point 4

                    b_co(1,4,7) = 0.8142685015
                    b_co(2,4,7) = -2.381048839e-3
                    b_co(3,4,7) = 3.7029531371e-6
                    b_co(4,4,7) = -3.043757962e-9
                    b_co(5,4,7) = 1.4027182897e-12
                    b_co(6,4,7) = -3.4205748313e-16
                    b_co(7,4,7) = 3.2053171516e-20
                    b_co(8,4,7) = 1.9723739733e-24
                    b_co(9,4,7) = -4.2458226186e-28
                    
                endif      

            !c
            !c Loop over the entire wave number domain
            !c
                  do 101 imu = 1, 9
                     
                    if ( imu.eq.1 ) then
                         imu_start = 1
                         imu_end   = 12
                    else if ( imu.eq.2 ) then
                         imu_start = 13
                         imu_end   = 43
                    else if ( imu.eq.3 ) then
                         imu_start = 44
                         imu_end   = 72
                    else if ( imu.eq.4 ) then
                         imu_start = 73
                         imu_end   = 93
                    else if ( imu.eq.5 ) then
                         imu_start = 94
                         imu_end   = 126
                    else if ( imu.eq.6 ) then
                         imu_start = 127
                         imu_end   = 147
                    else if ( imu.eq.7 ) then
                         imu_start = 148
                         imu_end   = 182
                    else if ( imu.eq.8 ) then
                         imu_start = 183
                         imu_end   = 205
                    else if ( imu.eq.9 ) then
                         imu_start = 206
                         imu_end   = 367
                    endif
            !c
                    num_band = imu_end-imu_start + 1
            !c
            !c  num_band is the number of narrow-bands inside each spectral band
            !csssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
            !c   Calculate the mean value of Planck blackbody intensity
            !c   over each spectral band
            !csssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

                       do i = 2, imax-1
                          do j = 2, jmax-1
                             xib(i,j) = 0.0d0
                          enddo
                       enddo
                       
                       do im = 1, imu_end-imu_start+1
                          wvnb  = 150.0d0+(im-1)*dwvnb+(imu_start-1)*dwvnb
                          wvnbm = 100.0d0*wvnb

                           do i = 2, imax-1
                              do j = 2, jmax-1
                                 tz=temp(i,j)
                                 tziv=1.0d0/tz
                                 xy1 = 2.0d0*c1*wvnbm*wvnbm*wvnbm
                                 xy2 = c2*wvnbm*tziv
                                 xib(i,j) = xib(i,j)+xy1/(dexp(xy2)-1.0d0)/num_band
                              enddo
                           enddo
                       enddo
            !c
            !c----------------------------------------------------------------
            !c   wavenumber at the centre of each spectral band, required 
            !c   for soot absorption coefficient and wall emission
            !c----------------------------------------------------------------
            !c
                    wvnb1 = 150.0d0+(imu_start-1)*dwvnb
                    wvnb2 = 150.0d0+(imu_end  -1)*dwvnb
                    wvnbm  = 0.5d0*(wvnb1+wvnb2)*100.0d0

                    do 102 igas = 1, ngauss
                         factor = dwvnb*100.0d0*weight(igas)*num_band
            !c
            !c... calculate the mixture absorption coefficients contributed
            !c    by CO, CO2, H2O, and soot particles
            !c
                      do i = 2, imax-1
                          do j = 2, jmax-1
                             T1 = temp(i,j)
                             T2 = temp(i,j)*T1
                             T3 = temp(i,j)*T2
                             T4 = temp(i,j)*T3
                             T5 = temp(i,j)*T4
                             T6 = temp(i,j)*T5
                             T7 = temp(i,j)*T6
                             T8 = temp(i,j)*T7
                             term_co = b_co(1,igas,imu) &
                                    + b_co(2,igas,imu)*T1  &
                                    + b_co(3,igas,imu)*T2  &
                                    + b_co(4,igas,imu)*T3  &
                                    + b_co(5,igas,imu)*T4  &
                                    + b_co(6,igas,imu)*T5  &
                                    + b_co(7,igas,imu)*T6  &
                                    + b_co(8,igas,imu)*T7  &
                                    + b_co(9,igas,imu)*T8

                             term_co = max(term_co, 0.0d0)

                             term_co2 = b_co2(1,igas,imu) &
                                     + b_co2(2,igas,imu)*T1 &
                                     + b_co2(3,igas,imu)*T2 &
                                     + b_co2(4,igas,imu)*T3 &
                                     + b_co2(5,igas,imu)*T4 &
                                     + b_co2(6,igas,imu)*T5 &
                                     + b_co2(7,igas,imu)*T6 &
                                     + b_co2(8,igas,imu)*T7 &
                                     + b_co2(9,igas,imu)*T8 

                             term_co2 = max(term_co2, 0.0d0)

                             term_h2o = b_h2o(1,igas,imu)  &
                                     + b_h2o(2,igas,imu)*T1 &
                                     + b_h2o(3,igas,imu)*T2 &
                                     + b_h2o(4,igas,imu)*T3 &
                                     + b_h2o(5,igas,imu)*T4&
                                     + b_h2o(6,igas,imu)*T5&
                                     + b_h2o(7,igas,imu)*T6&
                                     + b_h2o(8,igas,imu)*T7&
                                     + b_h2o(9,igas,imu)*T8

                             term_h2o = max(term_h2o, 0.0d0)

                             absor(i,j,igas) = (term_co*xco(i,j) &
                                            + term_co2*xco2(i,j)  &
                                            + term_h2o*xh2o(i,j))*PRESSURE &
                                            + 5.5*sfv(i,j)*wvnbm

                 
                         if(SOLIDS(I,J)) then 
                            absor(i,j,igas)=1.0d+10
                        endif

                         enddo
                      enddo
            !c
            !c   Boundary condition at the north surface of the domain
            !c    intensities leaving the north surface
            !c
                       do m = 1, mmax
                          do l = 1, lvalue(m)
                                if ( xmu(m,l).lt.0.0 ) then
                                   do i = 2, imax-1
                                     tz=tr1(i)
                                     tziv=1.0d0/tz
                                     xy1=2.0d0*c1*wvnbm*wvnbm*wvnbm
                                     xy2=c2*wvnbm*tziv
                                     xi_n(i,jmax-1) = xy1/(dexp(xy2)-1.0d0)
                                     xinorth(i)     = xy1/(dexp(xy2)-1.0d0)
                                   enddo
                                endif
                          enddo
                     enddo
            !c
            !c Boundary conditions at the west and east surfaces of the domain
            !c
                     do j = 2, jmax-1
                          do m = 1, mmax
                             do l = 1, lvalue(m)
                                if ( xi(m,l).gt.0.0 ) then
                                   tz=tz0(j)
                                     tziv=1.0d0/tz
                                     xy1=2.0d0*c1*wvnbm*wvnbm*wvnbm
                                     xy2=c2*wvnbm*tziv
                                     xi_w(2,j) =  xy1/(dexp(xy2)-1.0d0)
                                     xiwest(j)  =  xy1/(dexp(xy2)-1.0d0)
                                else
                                     tz=tz1(j)
                                     tziv=1.0d0/tz
                                     xy1=2.0d0*c1*wvnbm*wvnbm*wvnbm
                                     xy2=c2*wvnbm*tziv
                                     xi_e(imax-1,j) =  xy1/(dexp(xy2)-1.0d0)
                                     xieast(j)     =  xy1/(dexp(xy2)-1.0d0)
                                endif
                             enddo
                          enddo
                       enddo

            !c  calculate intensities I_(m,1/2) (notation used by Truelove(1978) for m = 1, 2, 3, 4

                       do m = 1, mmax
                          xi0(m) = xi(m,1)
                          xmu0(m) = -dsqrt(1.0d0-xi0(m)**2)
                       enddo

            !c Calulate the intensities I_(m,1/2), denoted in this program by inm(i,j,m)

                     do m = 1, mmax
                          if ( xi0(m).lt.0.0 ) then
                           do i = imax-1, 2, -1
                                do j = jmax-1,2,-1
                                   an = 2.0d0*pi*(z(i+1)-z(i  ))*yr(j+1)
                                   as = 2.0d0*pi*(z(i+1)-z(i  ))*yr(j)
                                   ae = 2.0d0*pi*rc(j)*(yr(j+1)-yr(j  ))
                                   aw = ae
                                   volu= 2.0d0*pi*rc(j)*(yr(j+1)-yr(j))*(z(i+1)-z(i))
                                   vp = 2.0d0*pi*(yr(j+1)-yr(j))*(z(i+1)-z(i))


                                   up = absor(i,j,igas)*volu*xib(i,j)  &
                                   - xmu0(m)*(an+as*(1.0d0-fs)/fs)*xi_n(i,j) &
                                   - xi0(m)*(ae+aw*(1.0d0-fs)/fs)*xi_e(i,j)
                               down = absor(i,j,igas)*volu  &
                                   - xmu0(m)*(as/fs+vp)   &
                                   - xi0(m)*aw/fs
                               xim(i,j,m,1) = up/down

                               xi_n(i,j-1) = xim(i,j,m,1)/fs - (1.0d0-fs)/fs*xi_n(i,j)
                               xi_e(i-1,j) = xim(i,j,m,1)/fs - (1.0d0-fs)/fs*xi_e(i,j)
                               
                                enddo
                             enddo
                          else
                             do i = 2, imax-1
                                do j = jmax-1,2,-1
                                   an = 2.0d0*pi*(z(i+1)-z(i))*yr(j+1)
                                   as = 2.0d0*pi*(z(i+1)-z(i))*yr(j)
                                   ae = 2.0d0*pi*rc(j)*(yr(j+1)-yr(j))
                                   aw = ae
                                   volu= 2.0d0*pi*rc(j)*(yr(j+1)-yr(j))*(z(i+1)-z(i))
                                   vp = 2.0d0*pi*(yr(j+1)-yr(j))*(z(i+1)-z(i))

                                   up = absor(i,j,igas)*volu*xib(i,j) &
                                   - xmu0(m)*(an+as*(1.0d0-fs)/fs)*xi_n(i,j) &
                                  + xi0(m)*(aw+ae*(1.0d0-fs)/fs)*xi_w(i,j)
                               down = absor(i,j,igas)*volu  &
                                   - xmu0(m)*(as/fs+vp)    &
                                  + xi0(m)*ae/fs
                               xim(i,j,m,1) = up/down

                               xi_n(i,j-1) = xim(i,j,m,1)/fs - (1.0d0-fs)/fs*xi_n(i,j)
                               xi_w(i+1,j) = xim(i,j,m,1)/fs - (1.0d0-fs)/fs*xi_w(i,j)

                               enddo
                            enddo
                        endif

                       enddo

                      
            !c
            !c Now calculate the radiation intensitiy field
            !c
            !c   Since radiation intensities at the central line (r = 0) are unavailable
            !c   the calculations should start from r = r1 to r = 0
            !c
            !c Case1: xi > 0, mu < 0
            !c

                  do m = 1, mmax
                       do l = 1, lvalue(m)
                          if ( xi(m,l).gt.0.0 .and. xmu(m,l).lt.0.0 ) then
                             do i = 2, imax-1
                                do j = jmax-1, 2, -1
                                   if ( j.eq.(jmax-1) ) xi_n(i,j) = xinorth(i)
                                   if ( i.eq.2    )     xi_w(i,j) = xiwest(j)
                                   an = 2.0d0*pi*(z(i+1)-z(i))*yr(j+1)
                                   as = 2.0d0*pi*(z(i+1)-z(i))*yr(j)
                                   ae = 2.0d0*pi*rc(j)*(yr(j+1)-yr(j))
                                   aw = ae
                                   volu= 2.0d0*pi*rc(j)*(yr(j+1)-yr(j))*(z(i+1)-z(i))
                                   up = absor(i,j,igas)*volu*xib(i,j)  &
                                   - xmu(m,l)*(an+as*(1.0d0-fs)/fs)*xi_n(i,j) &
                                   +  xi(m,l)*(aw+ae*(1.0d0-fs)/fs)*xi_w(i,j) &
                                   - (an-as)/w(m,l) &
                                   *(alpha(m,l+1)*(1.0d0-fa)/fa+alpha(m,l)) &
                                   *xim(i,j,m,l)

                                 down = absor(i,j,igas)*volu  &
                                   - xmu(m,l)*as/fs  &
                                   + xi(m,l)*ae/fs  &
                                   - (an-as)/w(m,l)/fa*alpha(m,l+1)


                             xip(i,j,m,l)  = up/down
                             xi_n(i,j-1)   = xip(i,j,m,l)/fs-(1.0d0-fs)/fs*xi_n(i,j)
                             xi_w(i+1,j)   = xip(i,j,m,l)/fs-(1.0d0-fs)/fs*xi_w(i,j)
                         xim(i,j,m,l+1) = xip(i,j,m,l)/fa-(1.0d0-fa)/fa*xim(i,j,m,l)

                                if ( j.eq.2 ) then
                                   xisouth(i,m,lvalue(m)+1-l) = xi_n(i,j-1)
                                endif

                                 enddo
                              enddo
                           endif
                       enddo
                    enddo


            !c Obtain boundary conditions at the south boundary for xi > 0 and mu > 0


            !c Case 2: xi < 0, mu < 0

                  do m = 1, mmax
                       do l = 1, lvalue(m)
                          if ( xi(m,l).lt.0.0 .and. xmu(m,l).lt.0.0 ) then
                             do i = imax-1, 2, -1
                                do j = jmax-1, 2, -1
                                   if ( j.eq.(jmax-1) ) xi_n(i,j) = xinorth(i)
                                   if ( i.eq.(imax-1) ) xi_e(i,j) = xieast(j)
                                   an = 2.0d0*pi*(z(i+1)-z(i))*yr(j+1)
                                   as = 2.0d0*pi*(z(i+1)-z(i))*yr(j)
                                   ae = 2.0d0*pi*rc(j)*(yr(j+1)-yr(j))
                                   aw = ae
                                   volu= 2.0d0*pi*rc(j)*(yr(j+1)-yr(j))*(z(i+1)-z(i))
                                   up = absor(i,j,igas)*volu*xib(i,j)  &
                                   -xmu(m,l)*(an+as*(1.0d0-fs)/fs)*xi_n(i,j) &
                                   - xi(m,l)*(ae+aw*(1.0d0-fs)/fs)*xi_e(i,j) &
                                   - (an-as)/w(m,l) &
                                   *(alpha(m,l+1)*(1.0d0-fa)/fa+alpha(m,l)) &
                                   *xim(i,j,m,l) 
                                 down = absor(i,j,igas)*volu &
                                   -xmu(m,l)*as/fs &
                                   - xi(m,l)*aw/fs &
                                   - (an-as)/w(m,l)/fa*alpha(m,l+1)
                          xip(i,j,m,l)   = up/down
                          xi_n(i,j-1)    = xip(i,j,m,l)/fs-(1.0d0-fs)/fs*xi_n(i,j)
                          xi_e(i-1,j)    = xip(i,j,m,l)/fs-(1.0d0-fs)/fs*xi_e(i,j)
                          xim(i,j,m,l+1) = xip(i,j,m,l)/fa-(1.0d0-fa)/fa*xim(i,j,m,l)

                                if ( j.eq.2 ) then
                                   xisouth(i,m,lvalue(m)+1-l) = xi_n(i,j-1)
                                endif

                                 enddo
                              enddo
                           endif
                       enddo
                    enddo


            !c Case 3: xi > 0 and mu > 0

                  do m = 1, mmax
                       do l = 1, lvalue(m)
                          if ( xi(m,l).gt.0.0 .and. xmu(m,l).gt.0.0 ) then
                             do i = 2, imax-1
                                do j = 2, jmax-1
                                   if ( i.eq.2 ) xi_w(i,j) = xiwest(j)
                                   if ( j.eq.2 ) xi_s(i,j) = xisouth(i,m,l)
                                   an = 2.0d0*pi*(z(i+1)-z(i))*yr(j+1)
                                   as = 2.0d0*pi*(z(i+1)-z(i))*yr(j)
                                   ae = 2.0d0*pi*rc(j)*(yr(j+1)-yr(j))
                                   aw = ae
                                   volu= 2.0d0*pi*rc(j)*(yr(j+1)-yr(j))*(z(i+1)-z(i))
                                   up = absor(i,j,igas)*volu*xib(i,j) &
                                   +xmu(m,l)*(as+an*(1.0d0-fs)/fs)*xi_s(i,j) &
                                   + xi(m,l)*(aw+ae*(1.0d0-fs)/fs)*xi_w(i,j) &
                                   - (an-as)/w(m,l) &
                                   *(alpha(m,l+1)*(1.0d0-fa)/fa+alpha(m,l)) &
                                   *xim(i,j,m,l)
                                 down = absor(i,j,igas)*volu &
                                   +xmu(m,l)*an/fs &
                                   + xi(m,l)*ae/fs &
                                   - (an-as)/w(m,l)/fa*alpha(m,l+1)
                          xip(i,j,m,l)   = up/down
                          xi_w(i+1,j)    = xip(i,j,m,l)/fs-(1.0d0-fs)/fs*xi_w(i,j)
                          xi_s(i,j+1)    = xip(i,j,m,l)/fs-(1.0d0-fs)/fs*xi_s(i,j)
                          xim(i,j,m,l+1) = xip(i,j,m,l)/fa-(1.0d0-fa)/fa*xim(i,j,m,l)
                                 enddo
                              enddo
                           endif
                       enddo
                    enddo

            !c Case 4: xi < 0 , mu > 0

                  do m = 1, mmax
                       do l = 1, lvalue(m)
                          if ( xi(m,l).lt.0.0 .and. xmu(m,l).gt.0.0 ) then
                             do i = imax-1, 2, -1
                                do j = 2, jmax-1
                                   if ( i.eq.(imax-1) ) xi_e(i,j) = xieast(j)
                                   if ( j.eq.2    )     xi_s(i,j) = xisouth(i,m,l)
                                   an = 2.0d0*pi*(z(i+1)-z(i))*yr(j+1)
                                   as = 2.0d0*pi*(z(i+1)-z(i))*yr(j)
                                   ae = 2.0d0*pi*rc(j)*(yr(j+1)-yr(j))
                                   aw = ae
                                   volu= 2.0d0*pi*rc(j)*(yr(j+1)-yr(j))*(z(i+1)-z(i))
                                   up = absor(i,j,igas)*volu*xib(i,j) &
                                   +xmu(m,l)*(as+an*(1.0d0-fs)/fs)*xi_s(i,j) &
                                   - xi(m,l)*(ae+aw*(1.0d0-fs)/fs)*xi_e(i,j) &
                                   - (an-as)/w(m,l) &
                                   *(alpha(m,l+1)*(1.0d0-fa)/fa+alpha(m,l)) &
                                   *xim(i,j,m,l)
                                 down = absor(i,j,igas)*volu &
                                   +xmu(m,l)*an/fs &
                                   - xi(m,l)*aw/fs &
                                   - (an-as)/w(m,l)/fa*alpha(m,l+1)
                          xip(i,j,m,l)   = up/down
                          xi_e(i-1,j)    = xip(i,j,m,l)/fs-(1.0d0-fs)/fs*xi_e(i,j)
                          xi_s(i,j+1)    = xip(i,j,m,l)/fs-(1.0d0-fs)/fs*xi_s(i,j)
                          xim(i,j,m,l+1) = xip(i,j,m,l)/fa-(1.0d0-fa)/fa*xim(i,j,m,l)
                                 enddo
                              enddo
                           endif
                       enddo
                    enddo
            !c
            !c Now construct the radiation sink term, qr = absor*(isum-4*pi*ib)
            !c
            !c First calculate isum
            !c

                  do i = 2, imax-1
                       do j = 2, jmax-1
                          xisum(i,j) = 0.0d0
                       enddo
                    enddo

                  
                          do i = 2, imax-1
                             do j = 2, jmax-1
                             do m = 1, mmax
                       do l = 1, lvalue(m)
                                xisum(i,j) = xisum(i,j) + xip(i,j,m,l)*w(m,l)
                                
                             enddo
                          enddo
                       enddo
                    enddo 

            !c calculate qr

                   do i = 2, imax-1
                       do j = 2, jmax-1
                          qr(i,j) = qr(i,j) &
                               + absor(i,j,igas) &
                               *(xisum(i,j)-4.0d0*pi*xib(i,j))*factor
                       enddo
                    enddo
                    
                    

             102  continue
             101  continue

            DEALLOCATE(xip)
            DEALLOCATE(xi_n)
            DEALLOCATE(xi_s)
            DEALLOCATE(xi_e)
            DEALLOCATE(xi_w)
            DEALLOCATE(xinorth)
            DEALLOCATE(xisouth)
            DEALLOCATE(xieast)
            DEALLOCATE(xiwest)
            DEALLOCATE(xim)
            DEALLOCATE(xisum)
            DEALLOCATE(xib)
            DEALLOCATE(absor)
            DEALLOCATE(tz0)
            DEALLOCATE(tz1)
            DEALLOCATE(tr0)
            DEALLOCATE(tr1)
            DEALLOCATE(ez0)
            DEALLOCATE(ez1)
            DEALLOCATE(er0)
            DEALLOCATE(er1)
            DEALLOCATE(xi0)
            DEALLOCATE(xmu0)
            DEALLOCATE(z)
            DEALLOCATE(yr)
            DEALLOCATE(rc)     
            DEALLOCATE(alpha)
            DEALLOCATE(lvalue)
            
        end subroutine DOM
        !*************************************************************************************************************
        
        !*****************************************************************************************************************************************
        SUBROUTINE SootMassFrac(SootAggDens, GramPerAgg, SootMassFraction)
        !This subroutine takes in soot aggregate number density per section, and the mass of each aggregate per section, and returns the total soot
        !mass fraction

            !variable declarations
            Double precision, dimension(:), intent(in) :: SootAggDens, GramPerAgg !soot aggregate number density (#/g-air), and grams per aggregate (g-s/#)
            Double precision, intent(out) :: SootMassFraction !soot mass fraction (g-s/g-air)
            integer :: p !loop parameter
            !done variable declarations
            
            SootMassFraction = 0.0d0 !initialize soot mass fraction to zero
            
            do p=1,MSection !loop over every section
                SootMassFraction = SootMassFraction + SootAggDens(p)*dexp(GramPerAgg(p)) !add up the cumulative mass fraction (#/g-air)*(g-s/#) = g-s/g-air
            enddo
            
        end subroutine SootMassFrac
        !***************************************************************************************************************************************************
                                    
! --------------------
END PROGRAM CoFlame


!*****************************************************************************************************************************************
!DECLARATION OF EXTERNAL SUBROUTINES & FUNCTIONS
! ****************************************************************************************************************************************************
SUBROUTINE DIFLOW(Flow,Diff,Acof)
!Calculates coefficient matrix values using the Power Law Scheme,
!as outlined in Patankar
    !variable declarations
    double precision Flow, Diff, Acof, GLOBdiflow
    !done variable declarations

    ACOF=DIFF
    IF(FLOW .EQ.0.0d0)RETURN
    GLOBdiflow=DIFF-DABS(FLOW)*0.1d0
    ACOF=0.0d0
    IF(GLOBdiflow .LE. 0.0d0 ) RETURN
    GLOBdiflow=GLOBdiflow/DIFF
    ACOF=DIFF*GLOBdiflow**5

end subroutine diflow
!*************************************************************************************************************************************************

! ************************************************************************************************************************************************
SUBROUTINE SOOTHCP(T,CPSOOT,HSOOT)
!This subroutine calculates the heat capacity and enthalpy of soot on a per carbon atom basis for a given temperature
    
    !local variable declarations
    implicit none
    double precision T,CPSOOT,HSOOT
    integer I
    double precision CC(8),CH(5)
    !done local variable declarations

    CC(1) = -0.2353267028d0
    CC(2) = 0.0159031622d0
    CC(3) = 8.5789194415d-5
    CC(4) = -1.7801434827d-7
    CC(5) = 1.5035398957d-10
    CC(6) = -6.4985704805d-14
    CC(7) = 1.4163210608d-17
    CC(8) = -1.2340325792d-21

    CH(1) = -0.1699516916d0
    CH(2) = -7.9088080344d-4
    CH(3) = 1.9343756673d-5
    CH(4) = -6.5338725677d-9
    CH(5) = 8.3384715475d-13

    CPSOOT = 0.0d0
    DO I = 1, 8
       CPSOOT = CPSOOT + CC(I)*T**(I-1.0d0)
    ENDDO

    HSOOT = 0.0d0
    DO I = 1, 5
       HSOOT = HSOOT + CH(I)*T**(I-1.0d0)
    ENDDO

    CPSOOT = CPSOOT*1.0E7/12.0d0
    HSOOT =  HSOOT*1.0E10

end subroutine SOOTHCP
!********************************************************************************************************************************************
SUBROUTINE Alpha_function(T)
!Alpha_surf funcion form by test A. Jerez (2016) for NDF of A. Fuentes C&F 160 (2013)
    !variable declarations
    double precision T, alpha_surf
    !done variable declarations

    alpha_surf = 0.078D0 * (0.0040d0*dexp(7000.0d0/T))   !Test 11 paper coke
    IF(alpha_surf.LT.0.) alpha_surf=0.0D0
    IF(alpha_surf.GT.0.078D0) alpha_surf= 0.078D0

end subroutine Alpha_function

!**********************************************************************************************************************************************
double precision function freePath(P, T)
!This function calculates the mean free path for air in cm (from AK, myf 1997)
    
    !local variable declaration
    implicit none
    double precision P,T
    !done local variable declarations

    freePath = 2.3701D-03 * T / P / 1.103D+05

    return
end
!*********************************************************************************************************************************************

!************************************************************************************************************************************************
double precision function viscosity(T)
!This function calculates the viscosity of air as a function of temperature (K) (from AK, myf 1997)

    !local variable declarations
    implicit none
    double precision T

    !done local variable declarations

    viscosity = 14.58D-06 * (T**1.5D0) / (T + 100.4D0)
    return
end
!***********************************************************************************************************************************
