function [Mopt,Mavg,EEopt_avg,Ropt_avg]= Mavg_JOptimizer(PLO,PLI,KgOpt,MgOpt,p,loading,lambdaS,GDR)
%For convergence, one of below condition should be satisfied
%   1\ EE Optimizer reaches the Nash Equilibrium while Ropt>GDR
%   2\ Ropt becomes lower than GDR
%   For the 2nd condition, the Rc Recovery process will be called to
%   conpansate the the loss in in Ropt.

EE_Converged=false;
RR_Converged=false;

M_iterate=MgOpt;

while(~EE_Converged && ~RR_Converged)
    
    Md = M_iterate;
    [Mopt,Mavg,EEopt_avg,Rc,Ropt_avg,EE] = EE_Optimizer(PLO,PLI,MgOpt,KgOpt,p,Md,loading,lambdaS);
    disp('$$$ EE Optimization CALL $$$')
    disp(['Mavg = ' num2str(Mavg)])
    disp(['Ropt_avg = ' num2str(Ropt_avg/1e6) ' Mbps' ])
    
    
    if (Ropt_avg < GDR)
        RR_Converged = true;
        disp('** EE Optimization STOPPED & CALLING RR Recovery **')
        [Mopt,Mavg,Ropt_avg,EEopt_avg] = RR_Optimizer(PLO,PLI,MgOpt,KgOpt,p,loading,lambdaS,GDR,Mopt,Rc,Mavg,EE);
                                              
        disp(['Mavg = ' num2str(Mavg)])
        disp(['Ropt_avg = ' num2str(Ropt_avg/1e6) ' Mbps' ])
        
    
    elseif (Ropt_avg >= GDR && Mavg == M_iterate)
        EE_Converged = true;
        disp('%% EE Optimization COMPLETED %%')
        
    
    else 
        M_iterate = Mavg;
        
    end
    
end

% Outage = Pi(end);

end

