function [EE,R,Ptot] = EE_R_Ptot(Kmax,K,Mc,Md,PLO,PLI,p)
% The size of Md should be 1x18 
if Mc>K
    % Rc Calculation
    B = 20e6;
    Tc = 1800;
    Bsigma2_dBm = -96; % dBm
    Bsigma2 = 10^((Bsigma2_dBm-30)/10);
    Md_temp = repmat(Md,size(PLI,1),1);
    
    meanI = Bsigma2*(1./PLO) + sum((p*PLI.*Md_temp),2)./PLO;
    alpha1= p*Mc/K./meanI;
    % Avreaging user rate over all test point
    R = mean( B*(1-Kmax/Tc).*log2(1+alpha1*(Mc-K)));
    
    % BaseBand Power Parameters
    P_COD = 0.1e-9;
    P_DEC = 0.8e-9;
    P_BT = 0.25e-9; %Power required for backhaul traffic (W/(bit/s))
    P_SYN = 2;
    L_BS = 12.8e9;
    P_BS= 1;
    
%     A = P_COD + P_DEC;
    A = P_COD + P_DEC + P_BT; %Parameter A in the power consumption model
    
    C00=P_SYN;
    C01=0;
    C02=0;
    C03=B/(3*Tc*L_BS);
    
    C10=P_BS;
    C11=(B/L_BS)*(2+1/Tc);
    % C12=(3*B)/L_BS;
    C12=(B*(3-1))/(Tc*L_BS);
    
    
    C0BB = A*K*R + (C00 + C01*K^1 + C02*K^2 + C03*K^3);
    C1BB = (C10 + C11*K^1 + C12*K^2);
    
    Poth = 18;
    Pmax_PA = p*10^0.8;
    eta = 0.80;
    PWR_PA =(1/eta)*sqrt(p*Pmax_PA);
    
    C0 = C0BB + Poth;
    C1 = C1BB + PWR_PA;
    
    Ptot = C0 + C1*Mc;
    
    EE = K*R/Ptot;
    
else
    R=0; Ptot=0; EE=0;
end

end