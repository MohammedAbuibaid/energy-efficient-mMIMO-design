function[Mopt,Mavg,EEopt_avg,Rc,Ropt_avg,EE] = EE_Optimizer(PLO,PLI,MgOpt,KgOpt,pgOpt,Md,loading,lambdaS)

EE = zeros(MgOpt,KgOpt);
Rc = zeros(MgOpt,KgOpt);

for K=1:KgOpt
    for M=1:MgOpt
        Mc = M;
        [EE(M,K),Rc(M,K)] = EE_R_Ptot_PA(PLO,PLI,KgOpt,K,Mc,Md*ones(1,18),pgOpt,[]);
    end
end

[EEopt, Mopt]=max(EE);

Ropt = zeros(1,KgOpt);
for K=1:KgOpt
    Ropt(K) =  Rc(Mopt(K),K);
end

Pi = zeros(1,KgOpt);
for K=1:KgOpt
    [Pi(K),~] = MGmm_SD_Queue(K,KgOpt,(loading/100)*lambdaS,1,Ropt);
end

Mavg = ceil(sum(Mopt.*Pi));
EEopt_avg = sum(EEopt.*Pi);
Ropt_avg = sum(Ropt.*Pi);

end

