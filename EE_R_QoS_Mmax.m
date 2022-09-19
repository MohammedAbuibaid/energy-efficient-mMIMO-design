function [EE_avg,Rc_avg] = EE_R_QoS_Mmax(PLO,PLI,KgOpt,MgOpt,p,loading,lambdaS,GDR)

%EE_R_QoS_Mmax Data Rate and QoS @ Mmax
%   Detailed explanation goes here

Mc = MgOpt; Md = MgOpt*ones(1,18);

EE = zeros(1,KgOpt);
Rc = zeros(1,KgOpt);
QoS = zeros(1,KgOpt);

for K=1:KgOpt
    [EE(1,K),Rc(1,K),~,QoS(1,K)] = EE_R_Ptot_PA(PLO,PLI,KgOpt,K,Mc,Md,p,GDR);
end

Pi = zeros(1,KgOpt);
for K=1:KgOpt
    [Pi(1,K),~] = MGmm_SD_Queue(K,KgOpt,(loading/100)*lambdaS,1,Rc);
end

EE_avg = floor(sum(EE.*Pi));
Rc_avg = floor(sum(Rc.*Pi));


disp('------- Avg @ Mmax ---------')
disp(['R_Mmax_avg = ' num2str(Rc_avg/1e6) ' Mbps' ])
% QoS_avg = floor(sum(1e2*QoS.*Pi))/1e2;
% disp(['QoS_Mmax_avg = ' num2str(QoS_avg)])
disp('----------------------------')

end

