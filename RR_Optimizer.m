function [Mopt,Mavg,Ropt_avg,EEopt_avg] = RR_Optimizer(PLO,PLI,MgOpt,KgOpt,pgOpt,loading,lambdaS,GDR,Mopt,Rc,Mavg,EE)
converged = false;

a = 0;
b = MgOpt - Mopt(1);
c=0;

c_m1=b;


Mx = Mopt;

while (~converged)
    
    Ropt = zeros(1,KgOpt);
    EEopt= zeros(1,KgOpt);
    for K=1:KgOpt
        Ropt(K) =  Rc(Mx(K),K);
        EEopt(K)=  EE(Mx(K),K);
    end
    
    Pi = zeros(1,KgOpt);
    for K=1:KgOpt
        [Pi(K),~] = MGmm_SD_Queue(K,KgOpt,(loading/100)*lambdaS,1,Ropt);
    end
    
    Ropt_avg = sum(Ropt.*Pi);
    EEopt_avg = sum(EEopt.*Pi);
    TOL = Ropt_avg - GDR;
    if c>0
        %         disp('-------------');
        %         disp(['Step C = ' num2str(c) ])
        %         disp(['Mavg = ' num2str(Mavg) ])
        %         disp(['AvgURperBS = ' num2str(AvgURperBS/1e6) ' Mbps']);
        %         disp(['TOL = ' num2str(TOL/1e6) ' Mbps'])
    end
    
    
    
    if (c == c_m1 || (TOL>=0 && TOL<=0.3*1e6)|| range(Mx) == 0)
        
        converged = true;
        disp('%% RR Recovery COMPLETED %%')
        disp(['C Value  = ' num2str(c) ])
        %         disp(['Mavg = ' num2str(Mavg) ])
        %         disp(['AvgURperBS = ' num2str(AvgURperBS/1e6) ' Mbps']);
        
        
        
    else
        
        if TOL<0
            a = c;
        else
            b = c;
        end
        
        c_m1 =c;
        c = ceil((b+a)/2);
        
        Mx = Mopt+c;
        Mx(Mx > MgOpt)=MgOpt;
        
        Ropt = zeros(1,KgOpt);
        for K=1:KgOpt
            Ropt(K) =  Rc(Mx(K),K);
        end
        Pi = zeros(1,KgOpt);
        for K=1:KgOpt
            [Pi(K),~] = MGmm_SD_Queue(K,KgOpt,(loading/100)*lambdaS,1,Ropt);
        end
        
        Mavg = floor(sum(Mx.*Pi));
        
        
        Md_temp= Mavg*ones(1,18);
        Rc = zeros(MgOpt,KgOpt);
        EE = zeros(MgOpt,KgOpt);
        for K=1:KgOpt
            for M=1:MgOpt
                Mc = M; Md = Md_temp;
                [EE(M,K),Rc(M,K),~] = EE_R_Ptot_PA(PLO,PLI,KgOpt,K,Mc,Md,pgOpt,[]);
            end
        end
        
        
        
    end
    
    
    
end

% disp(['AvgURperBS = ' num2str(AvgURperBS/1e6) ' Mbps']);
% disp(['Mavg = ' num2str(Mavg) ])




end

