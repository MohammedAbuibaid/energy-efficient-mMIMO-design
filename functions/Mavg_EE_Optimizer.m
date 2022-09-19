function [Mopt_EEopt,Mavg_EEopt,EEopt,Ravg_EEopt]= Mavg_EE_Optimizer(PLO,PLI,KgOpt,MgOpt,pgOpt,loading,lambdaS)

EE_Converged=false;
M_iterate=MgOpt;

while(~EE_Converged)
    
    Md = M_iterate;
    [Mopt_EEopt,Mavg_EEopt,EEopt,~,Ravg_EEopt] = EE_Optimizer(PLO,PLI,MgOpt,KgOpt,pgOpt,Md,loading,lambdaS);
    disp('$$$ EE Optimization CALL $$$')
    disp(['Mavg = ' num2str(Mavg_EEopt)])
    disp(['Ropt_avg = ' num2str(Ravg_EEopt/1e6) ' Mbps' ])
    
    if (Mavg_EEopt == M_iterate)
        EE_Converged = true;
        disp('%% EE Optimization COMPLETED %%')
    else 
        M_iterate = Mavg_EEopt;
        
    end
    
end

% Outage = Pi(end);

end

