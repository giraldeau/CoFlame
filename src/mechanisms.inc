!This file contains the code to assign species indicies based on the selected mechanism

if(Mechanism.Eq.1) then ! DLR 2009 - N.A Slavinskaya and P. Frank C&F 156 (2009).
    KK=94
    IC2H2 = 3
    IBAPYR = 89
    IBAPYRS = 87
    IBGHIF = 88
    IO2 = 14
    IH2O = 16
    IH = 24
    IH2 = 1
    IC2H6 = 5
    IN2 = 46
    ICH4 = 2
    IO = 37
    IOH = 38
    ICO = 18
    ICO2 = 19
    IC2H4 = 4
    IHE = 94
    IAR = 47
    IC6H6 = 48
	IA4 = 73
    
    species_header='VARIABLES ="r (cm)","z (cm)","U","V","T","Density","Pressure","PC", &
        "H2","CH4","C2H2","C2H4","C2H6","C3H4","C3H6","C4H2",&
		"C4H4","C4H6","C4H7","C4H8","C3H8","O2","C2","H2O",&
		"H2O2","CO","CO2","CH2O","CH2CO","CH3CHO","C","H",& 
		"CH","CH2","CH2(S)","CH3","C2H","C2H3","C2H5","C3H2",&
		"H2CCCH","H2CCCCH","C3H5","i-C4H5","O","OH","HO2","HCO",&
		"CH3O","CH2OH","HCCO","CH2HCO","CH3CO","N2","AR","A1",&
		"A1-","C6H5OH","C6H5O","C5H5","C5H6","A1C2H","C7H8","A1C2H3",&
		"A1C2H3*","n-C8H7","C7H7","A1C2H-","A2","A2-","INDENE","INDENYL",&  
		"P2","P2-","A2C2H*","A2R5","A2C2H","A2R5-","A3","A3-",&
		"A4","A4-","C18H12","C18H11","A3C2H","A3C2H*","C4H","C4",&
		"C6H2","C8H2","C6H","C8H","A4C2H","A4C2H*","BAPYR*S","BGHIF",&
		"BAPYR","A2CH3","A2CH2","A3CH3","A3CH2","HE"'    		
endif    

if(Mechanism.Eq.8) then ! DLR 2011 - S. Dworkin et al. C&F 158 (2011)
    KK=94
	IC2H2 = 3
    IBAPYR = 89
    IBAPYRS = 87
    IBGHIF = 88
	IO2 = 14
	IH2O = 16
	IH = 24
	IH2 = 1
	IC2H6 = 5
	IN2 = 46
	ICH4 = 2
	IO = 37
	IOH = 38
	ICO = 18
	ICO2 = 19
	IC2H4 = 4
	IHE = 94
	IAR = 47
	IC6H6 = 48
	IA4 = 73

    species_header='VARIABLES ="r (cm)","z (cm)","U","V","T","Density","Pressure","PC", &
        "H2","CH4","C2H2","C2H4","C2H6","C3H4","C3H6","C4H2",&
		"C4H4","C4H6","C4H7","C4H8","C3H8","O2","C2","H2O",&
		"H2O2","CO","CO2","CH2O","CH2CO","CH3CHO","C","H",&
		"CH","CH2","CH2(S)","CH3","C2H","C2H3","C2H5","C3H2",&
		"H2CCCH","H2CCCCH","C3H5","i-C4H5","O","OH","HO2","HCO",&
		"CH3O","CH2OH","HCCO","CH2HCO","CH3CO","N2","AR","A1",&
		"A1-","C6H5OH","C6H5O","C5H5","C5H6","A1C2H","C7H8","A1C2H3",&
		"A1C2H3*","n-C8H7","C7H7","A1C2H-","A2","A2-","INDENE","INDENYL",&
		"P2","P2-","A2C2H*","A2R5","A2C2H","A2R5-","A3","A3-",&
		"A4","A4-","C18H12","C18H11","A3C2H","A3C2H*","C4H","C4",&
		"C6H2","C8H2","C6H","C8H","A4C2H","A4C2H*","BAPYR*S","BGHIF",&
		"BAPYR","A2CH3","A2CH2","A3CH3","A3CH2","HE"'
endif

if(Mechanism.Eq.3) then !DLR 2014 - V. Chernov et al. C&F 161 (2014).
!    KK=102
    KK=94
    IC2H2 = 3
    IBAPYR = 89
    IBAPYRS = 87
    IBGHIF = 88
    IO2 = 14
    IH2O = 16
    IH = 24
    IH2 = 1
    IC2H6 = 5
    IN2 = 46
    ICH4 = 2
    IO = 37
    IOH = 38
    ICO = 18
    ICO2 = 19
    IC2H4 = 4
    IHE = 94
    IAR = 47
    IC6H6 = 48
    IA4 = 73

    IC3H6 = 7
    IC4H6 = 10
    IC4H8 = 12
    IC3H8 = 13

    species_header='VARIABLES ="r (cm)","z (cm)","U","V","T","Density","Pressure","PC", &
    		"H2","CH4","C2H2","C2H4","C2H6","C3H4","C3H6","C4H2","C4H4","C4H6", &
                "C4H7","C4H8","C3H8", &
    		"O2","C2","H2O","H2O2","CO","CO2","CH2O","CH2CO","CH3CHO","C","H", &
    		"CH","CH2","CH2(S)","CH3","C2H","C2H3","C2H5","C3H2","H2CCCH", &
    		"H2CCCCH","C3H5","i-C4H5","O","OH","HO2","HCO","CH3O","CH2OH", &
    		"HCCO","CH2HCO","CH3CO","N2","AR","A1",  &
    		"A1-","C6H5OH","C6H5O","C5H5","C5H6","A1C2H","C7H8","A1C2H3",  &
    		"A1C2H3*","n-C8H7","C7H7","A1C2H-","A2","A2-","INDENE","INDENYL", &
    		"P2","P2-","A2C2H*","A2R5","A2C2H","A2R5-","A3","A3-","A4","A4-", &
    		"C18H12","C18H11","A3C2H","A3C2H*","C4H","C4","C6H2","C8H2","C6H", &
    		"C8H","A4C2H","A4C2H*","BAPYR*S","BGHIF","BAPYR","A2CH3","A2CH2", &
    		"A3CH3","A3CH2","HE"'
endif

if(Mechanism.Eq.4) then !ABF - J. Appel, H. Bockhorn and M. Frenklach C&F 121 (2000).
    KK=101
	IC2H2 = 23
    !IBAPYR = 89
    !IBAPYRS = 87
    !IBGHIF = 88
    IO2 = 4
    IH2O = 6
    IH = 2
    IH2 = 1
    IC2H6 = 27
    IN2 = 33
    ICH4 = 14
    IO = 3
    IOH = 5
    ICO = 15
    ICO2 = 16
    IC2H4 = 25
    !IHE = 94
    IAR = 31
	IC6H6 = 63
	IA4 = 95

    species_header='VARIABLES ="r (cm)","z (cm)","U","V","T","Density", &
    	       "H2","H","O","O2","OH","H2O","HO2","H2O2", &
    	       "C","CH","CH2","CH2(S)","CH3","CH4","CO","CO2", &
    	       "HCO","CH2O","CH2OH","CH3O","CH3OH","C2H","C2H2","C2H3", &
    	       "C2H4","C2H5","C2H6","HCCO","CH2CO","HCCOH", &
    	       "AR","C2H3O","N2","C2O","C3H2","C3H3","PC3H4","AC3H4", &
    	       "C4H","C4H2","H2C4O","n-C4H3","i-C4H3","C4H4","n-C4H5","i-C4H5", &
    	       "C4H6","C4H612","C5H2","C5H3","C6H","C6H2","C6H3","l-C6H4", &
    	       "c-C6H4","n-C6H5","i-C6H5","l-C6H6","n-C6H7","i-C6H7","c-C6H7", &
    	       "C6H8","A1","A1-","C6H5O","C6H5OH","C5H6","C5H5","C5H5O", &
    	       "C5H4OH","C5H4O","n-A1C2H2","i-A1C2H2","A1C2H3","A1C2H3*", &
    	       "A1C2H","A1C2H*","A1C2H-","A1C2H)2","A2","A2-1", &
    	       "A2-2","naphthyne","A2C2HA","A2C2HB","A2C2HA*","A2C2HB*", &
    	       "A2C2H2","A2C2H)2","A3-1","A3-4","A3","A3C2H","A3C2H2", &
    	       "A4","P2","P2-","P2-H","A4-","A2R5","A2R5-"'
endif

if(Mechanism.Eq.5) then !Marinov (1998) Marinov, Combustion and Flame 114 (1998).
    KK=155
	IC2H2 = 19
    !IBAPYR = 89
    !IBAPYRS = 87
    !IBGHIF = 88
    IO2 = 11
    IH2O = 16
    IH = 2
    IH2 = 1
    IC2H6 = 23
    IN2 = 155
    ICH4 = 3
    IO = 12
    IOH = 13
    ICO = 10
    ICO2 = 9
    IC2H4 = 21
    !IHE = 94
    IAR = 154
	IC6H6 = 43
	IA4 = 130

	species_header='VARIABLES ="r (cm)","z (cm)","U","V","T","Density","Pressure","PC", &
			"h2","h", &
			"ch4","ch3","ch2","ch","ch2o","hco","co2","co", &
			"o2","o","oh","ho2","h2o2","h2o",&
			"c2h","hcco","c2h2","c2h3","c2h4","c2h5","c2h6", &
			"ch2oh","ch3o","hccoh","h2ccch", &
			"c3h2","ch2(s)","ch2co","c4h2","c5h2","c5h3","c6h2","c2","c2o","c","hcoh", &
			"ch2chcch","hcchcch","h2cccch","c6h5","c6h6","c6h5o", &
			"ch2chchch","ch2chcch2", &
			"h2c4o","ch2chchch2","ch3oh","ch2hco", &
			"c-c5h5","c-c5h5o","c-c5h6","c6h5oh","c10h8","c-c5h4oh","c-c5h4o", &
			"c4h10","pc4h9","sc4h9","c4h8-1","c4h8-2","c4h7", &
			"c3h6","ac3h5","pc3h5", &
			"sc3h5","ch2chcho","pc3h4","ac3h4","ch3co","c6h5ch2","c6h5ch3", &
			"ch2chco","ch3chco","ch3hco","chocho", &
			"c3h8","ic3h7","nc3h7","ch3chcch2","ch3ccch2", &
			"c6h5c2h5","c10h10","c10h9","c10h7", &
			"c6h4c2h3","c6h4c2h","c6h5c2h3","c6h5c2h","c10h7o", &
			"c6h5cho","c6h5co","hoc6h4ch3","c6h5ch2oh","oc6h4o", &
			"oc6h4ch3","c10h7ch2","c10h7ch3", &
			"ch3c6h4ch3","indene","indenyl","acenphthln","phnthrn", &
			"ch3chcch","ch3ch2cch","ch3c6h4ch2", &
			"phnthryl-1","phnthryl-9","flrnthn", &
			"kkk","kk0", &
			"c-2*4c6h6o","c-c5h7","l-c5h8","l-c5h7", &
			"c10h7c2h5","c10h7c2h3","c10h7oh","c6h5cco", &
			"c6h5chch","c6h5cch2","c10h7cch2","c10h7cch","c10h6cch", &
			"anthracn","ch3indene","ch3indenyl","kk2","pyrene", &
			"fluoryl","fluorene","h2cccch2","chchcho","hcccho","hccco", &
			"kk3","kk4","phnthroxy-1","kk6", &
			"bz(a)ndene","bz(a)ndnyl","phnthrol-1","kk22", &
			"biphenyl","kk1", &
			"bz(ghi)fln","fulvene","kk7","kk8", &
			"kk9","cyc6h7", &
			"fulvenyl","ar","n2" '
endif

if(Mechanism.Eq.2) then !Hybrid (2017), a combined from DLR for PAH and DME (2017).
    KK=143
    IC2H2 = 3
    IBAPYR = 89
    IBAPYRS = 87
    IBGHIF = 88
    IO2 = 14
    IH2O = 16
    IH = 24
    IH2 = 1
    IC2H6 = 5
    IN2 = 46
    ICH4 = 2
    IO = 37
    IOH = 38
    ICO = 18
    ICO2 = 19
    IC2H4 = 4
    IHE = 94
    IAR = 47
!
    IC6H6 = 48
    IA4 = 73
    IDME = 125

        species_header='VARIABLES ="r (cm)","z (cm)","U","V","T","Density","Pressure","PC", &
                        "H2","CH4","C2H2","C2H4","C2H6", &
                        "C3H4","C3H6","C4H2","C4H4","C4H6","C4H7","C4H8","C3H8", &      
                        "O2","C2","H2O","H2O2","CO","CO2",&
                        "CH2O","CH2CO","CH3CHO","C","H","CH","CH2", &          
                        "CH2S","CH3","C2H","C2H3", &
                        "C2H5","C3H2","H2CCCH","H2CCCCH","C3H5","i-C4H5","O","OH","HO2",&
                        "HCO","CH3O",&
                        "CH2OH","HCCO","CH2HCO","CH3CO","N2","AR","A1","A1-",&
                        "C6H5OH","C6H5O","C5H5","C5H6",&
                        "A1C2H","C7H8","A1C2H3","A1C2H3*","n-C8H7","C7H7","A1C2H-" , &
                        "A2","A2-", &
                        "INDENE","INDENYL","P2","P2-","A2C2H*","A2R5","A2C2H","A2R5-",     &
                        "A3","A3-", &
                        "A4","A4-","C18H12", &
                        "C18H11","A3C2H","A3C2H*","C4H","C4","C6H2","C8H2","C6H", &
                        "C8H", "A4C2H", &
                        "A4C2H*","BAPYRS","BGHIF","BAPYR","A2CH3","A2CH2","A3CH3",  &
                        "A3CH2","HE", &
                        "CH3OH","C2H5OH","PC2H4OH","SC2H4OH","CH2CHO","CH3COCH3",                 &
                        "CH3COCH2","C2H5CHO","C2H5CO","C2H5O","CH3O2",                     &
                        "C2H5O2","CH3O2H","C2H5O2H","C2H3O1-2", &
                        "CH3CO2","C2H4O1-2","C2H4O2H","O2C2H4OH","CH3CO3",        &
                        "CH3CO3H","C2H3CO","C2H3CHO","CH3COCH2O2","CH3COCH2O2H","CH3COCH2O",  &
                        "HCO3H","HCO3","HCO2","O2C2H4O2H","CH3OCH3",     &
                        "CH3OCH2","CH3OCH2O2","CH2OCH2O2H","CH3OCH2O2H","CH3OCH2O",                 &
                        "O2CH2OCH2CHO","HO2CH2OCHO","OCH2OCHO","HOCH2OCO","HOCH2O","HCO2H",                     &
                        "CH3OCHO","CH3OCO","CH2OCHO","CH3OCH2OH","HOCH2O2H","OCH2O2H",          &
                        "HOCH2O2" '
 endif

if(Mechanism.Eq.7) then !Hybrid (2017), a combined from DLR for PAH and Kaust for main oxidation Combustion and Flame Consalvi et al. (2017).
    KK=175
    IC2H2 = 25
    IBAPYR = 101
    IBAPYRS = 169
    IBGHIF = 170
    IO2 = 8
    IH2O = 7
    IH = 3
    IH2 = 4
    IC2H6 = 29
    IN2 = 1
    ICH4 = 16
    IO = 5
    IOH = 6
    ICO = 17
    ICO2 = 23
    IC2H4 = 27
    IHE = 175
    IAR = 2
!
    IC6H6 = 74
    IA4 = 112
    IC7H16 = 118
    IC7H8 = 83        ! C6H5CH3

        species_header='VARIABLES ="r (cm)","z (cm)","U","V","T","Density","Pressure","PC", &
                        "N2","AR","H","H2","O", &
                        "OH","H2O","O2","HO2","H2O2","C","CH","CH2", &
                        "CH2S","CH3","CH4","CO","HCO","CH2O",&
                        "CH2OH","CH3O","CH3OH","CO2","C2H","C2H2","C2H3", &
                        "C2H4","C2H5","C2H6","HCCO", &
                        "CH2CO","HCCOH","C2O","C2H3O","CH3CHO","C3H2","C3H3","AC3H4","PC3H4",&
                        "C3H5-A","C3H5-S",&
                        "C3H5-T","C3H6","n-C3H7","i-C3H7","C4H","C4H2","n-C4H3",& 
                        "i-C4H3","C4H4","n-C4H5",&
                        "i-C4H5","H2C4O","C4H6-13","C4H6-12","i-C4H8","i-C4H9","t-C4H9" , &
                        "C5H2","C5H3", &
                        "C5H5","C5H6","C5H4O","C5H4OH","C5H5O","C6H","C6H2","C6H3",     &
                        "l-C6H4","c-C6H4", &
                        "n-C6H5","i-C6H5","l-C6H6", & 
                        "A1","A1-","C6H5O","C6H5OH","n-C6H7","i-C6H7","c-C6H7","C6H8", &
                        "C6H5CH2","C6H5CH3","A1C2H","A1C2H-","A1C2H3","A1C2H3*","C2",  &
                        "FH2","H2CCCCH", &
                        "C3H4","A2","A2-","INDENE","INDENYL","C4H7",                 &
                        "A2C2H*","A2C2H","A2CH3","A2CH2","BAPYR",                     &
                        "A3","A3-","A3C2H","A3C2H*", &
                        "FCO2","P2","P2-","n-C8H7","A2R5",        &
                        "A2R5-","A4","A4-","C4","C8H","CH2HCO","CH3CO",                        &
                        "C7H16","C7H15","C7H14","C5H11-1","C5H10-1",                           &
                        "C5H9","C5H8","C4H9-1","C4H9-2","C4H8-1","C4H8-2",                     &
                        "n-C4H7","C4H7-B","i-C8H18","M224Pyl1","M224Pyl3","M244Pyl2",          &
                        "M244Pyl1","M24Pen1","M24Pen2","M44Pen2","M44Pen1",                    &
                        "M2Pdn14","M4Pdn13","M2Pdn13","neoC5H11","M2Bten1","M2Bten2",          &
                        "a-C5H9","b-C5H9","C5H8-13","C5H8-14","M2Bdn13","C4H7-S",              &
                        "i-C4H7","C10H22","C10H21","C10H20","C6H13-1","C6H13_2","C6H12-1",     &
                        "C3H8","C10H10","C10H9","c-C5H5","C5H5-2","c-C5H6","A3CH3","A3CH2",    &
                        "A4C2H","A4C2H*","BAPYR*S","BGHIF","C7H7","C18H11","C18H12","C8H2",    &
                        "HE" '
 endif  
