%list of parameters used in Lui_Model_Solver
%from Holder-Pearson et al 2018 and Van Cauter et al 1992

function val=initiate_val_struct()
    val.d1=-log(0.5)/20;
    
    val.k2=0.0104; %rate constant defining insulin diffusion from the subcutaneous space into local interstitium
    val.k3=0.060; %rate constant defining insulin being absorbed into the plasma
    val.kdi=0.006; %rate constant defining the breakdown of insulin in the local interstitium

    val.pg=0.004; %non insulin mediated uptake - LOTZ
    val.EGP=0.96; %endogenous glucose production
    val.CNS=0.3; %glucose consumption attributed to the central nervous system
    val.alpha_G=0.0154; %insulin binding saturation constant
    val.alpha_I=0.0017; %hepatic insulin saturation constant
    val.Vp=4; %distribution volume plasma
    val.umin=16.7; %min endogenous insulin secretion
    val.umax=267; %max endogenous insulin secretion
    val.ksec=14.9;
    val.koffset=-50; 

    a=log(2)/4.95;
    b=log(2)/32.4;
    F=0.76;
    val.ck2=F*(b-a)+a;
    val.ck3=(a*b)/val.ck2;
    val.ck1=a+b-val.ck2-val.ck3;
    val.Vq=val.Vp*(val.ck1/val.ck2);
    val.nK=val.ck3;
    val.nI=val.Vq*val.ck2;
    gamma=0.5;
    val.nC=(val.nI/val.Vq)*((1/gamma)-1);
    val.Vg=(val.Vp+val.Vq)*1.2;
    
    val.S=[]; %C-pep secretion rate (pmol/L/min)
    val.Y_Sec=[]; %peripheral cpep calculated with S (pmol/L)
    val.Uen=[]; %C-pep secretion rate (mU/min)
    
    val.T_Cpep=[]; %Time vector for c-pep (min)
    val.C=[]; %Central c-pep concentration (pmol/L)
    val.Y_for=[]; %peripheral c-pep concentration forward calc (pmol/L)
    
    val.nL=[]; %hepatic insulin clearance rate (1/min)
    val.xL=[]; %first pass constant as endogenous secretion
    
    val.SI=[]; %Insulin sensitivity (L/mU*min)
    val.da=[];
    val.d2=[];
    
    val.P1=[];
    val.P2=[];
    val.Isc=[]; %subcutaneous space into which the insulin is injected (mU)
    val.Qlocal=[]; %local interstitium (mU)
    val.G=[]; %plasma glucose (mmol/L)
    val.I=[]; %plasma insulin (mU/L)
    val.Q=[]; %interstitial insulin (mU/L)
    val.T_GIQ=[]; %Time vector for GIQ data (min)
    
    val.SSRg=[];
    val.SSRi=[];
    
    val.errg=[];
    val.erri=[];
end 