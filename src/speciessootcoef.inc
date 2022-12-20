!This file contains the code for calculating the matrix coefficients for the species and soot equations   
   
    DO J=2,M2
        
        !calculate Asouth for the first CV on each CPU
        !if(.NOT.SOLIDCVP(1,J)) then  !if the CV below is not in a solid region
            if(myid.eq.0)then
                DIFF=ARX(J)*GAMS(1,J,K)/XDIF(2)
            else      
                DIFF=ARX(J)*2.d0*GAMS(1,J,K)*GAMS(2,J,K)/((XCV(1)*GAMS(2,J,K)+ XCV(2)*GAMS(1,J,K)+smallnum))     
            endif
            CALL DIFLOW(FLOWN(1,J),DIFF,ACOF)
            AIMS(2,J,K)=ACOF+DMAX1(0.0D0,FLOWN(1,J))          

        DO I=2,L2
            !calculate Asouth for the CV above and Anorth for the current CV for all other CVs
            IF(myid .eq. (NUMP-1) .and. I .EQ. L2) then !special case if at the last CV in the domain
                DIFF=ARX(J)*GAMS(L1,J,K)/XDIF(L1)
            ELSE                       	               
                DIFF=ARX(J)*2.d0*GAMS(I,J,K)*GAMS(I+1,J,K)/((XCV(I)*GAMS(I+1,J,K)+ XCV(I+1)*GAMS(I,J,K)+smallnum))  
            ENDIF
            CALL DIFLOW(FLOWN(I,J),DIFF,ACOF)
            AIMS(I+1,J,K)=ACOF+DMAX1(0.0D0,FLOWN(I,J))
            AIPS(I,J,K)=AIMS(I+1,J,K)-FLOWN(I,J)

            !calculate Awest for the CV to the right, and Aeast for the current CV
            AREA=RMN(J+1)*XCV(I)
            if(J .EQ. M2) then !special case when at the end of the domain
                DIFF=AREA*GAMS(I,M1,K)/YDIF(M1)
            else                             	                        
                DIFF=AREA*2.d0*GAMS(I,J,K)*GAMS(I,J+1,K)/(YCV(J)*GAMS(I,J+1,K)+ YCV(J+1)*GAMS(I,J,K)+smallnum)     
            endif
            CALL DIFLOW(FLOWE(I,J),DIFF,ACOF)
            AJMS(I,J+1,K)=ACOF+DMAX1(0.0D0,FLOWE(I,J))
            AJPS(I,J,K)=AJMS(I,J+1,K)-FLOWE(I,J)  
            
            !Determinne Ap and source term for each species' equation
            VOL=YCVR(J)*XCV(I)  !volume of the control volume
            APT=RHO(I,J)/DT
            APS(I,J,K)=APS(I,J,K)-APT !APS partially calculated in gammaANDsource