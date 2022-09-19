function lambdaS = Searching_lambdaS(Rc,Kmax,GoS,plott)
% Rc should be in Mbps in order to speed up finding lambdaS  

WindowSize = 100e6;  % 100 Mbps
lam_str=0;
lam_end =WindowSize;
    
while(WindowSize>100)
    [Pi_Kmax,~] = MGmm_SD_Queue(Kmax,Kmax,lam_end,1,Rc);
    while (Pi_Kmax<GoS)
        lam_end = lam_end + WindowSize;
        lam_str = lam_end - WindowSize;
        % disp(['Current lambdaS interval = [ ' num2str(lam_str/1e9) ' Gbps, ' num2str(lam_end/1e9) ' Gbps ]']);
        [Pi_Kmax,~] = MGmm_SD_Queue(Kmax,Kmax,lam_end,1,Rc);
        % disp(['Pi_Kmax = ' num2str(Pi_Kmax)]);

    end
    WindowSize = WindowSize/100;
    lam_end = lam_str + WindowSize;
end

%disp(['Finding lambdaS in the interval = [' num2str(lam_str) ', ' num2str(lam_end) ']' ' bps']);

% Sampling the window into 500 steps and determining value of lambdaS
n_steps = 1000;
Pi_Kmax = zeros(1,n_steps);
lam_value = zeros(1,n_steps);
for lam_Ind = 1:1:n_steps
    lam_value(1,lam_Ind) = lam_str+((lam_end-lam_str)/n_steps)*lam_Ind;
    [Pi_Kmax(1,lam_Ind),~] = MGmm_SD_Queue(Kmax,Kmax,lam_value(lam_Ind),1,Rc);
end

[Pi_error,lambdaS_Ind ] = min(abs(Pi_Kmax - GoS));
lambdaS = lam_value(lambdaS_Ind);
Pi_lambdaS = Pi_Kmax(lambdaS_Ind);
% disp(['lambdaS value = ' num2str(lambdaS/1e9) ' Gbps' '; Pi error = ' num2str(Pi_error)]);






if plott == true
    figure
    plot(lam_value,Pi_Kmax);
	grid on; hold on;
    plot(lambdaS,Pi_lambdaS,'ro','linewidth',1)
    title('LambdaS Approximation')
    xlabel('LambdaS, Mbps')
    ylabel('GoS')
end

end