KgOpt = 100;

load('UserPDF.mat')
%% Traffic model: M/G/m/m state-dependent queue
loading = [10 20 30 40 50 60 70 80 90 100];
GoS = 0.02;
% Searching for lambdaS (FileSize)
Rc = zeros(1,KgOpt);
Mc = MgOpt; Md = MgOpt*ones(1,18);

for K= 1:KgOpt
    [~,Rc(K),~] = EE_R_Ptot_PA(PLO,PLI,KgOpt,K,Mc,Md,pgOpt,[]);
end
lambdaS = Searching_lambdaS(Rc,KgOpt,GoS,false);
disp(['lambdaS = ' num2str(lambdaS/1e9) ' GB']);



%% User Probability Dist at 10,50,100 % loading
Pi_10 = zeros(1,KgOpt);
Pi_40 = zeros(1,KgOpt);
Pi_70 = zeros(1,KgOpt);
Pi_100 = zeros(1,KgOpt);
for K=1:KgOpt
    Pi_10(K) = MGmm_SD_Queue(K,KgOpt,0.10*lambdaS,1,Rc);
    Pi_40(K) = MGmm_SD_Queue(K,KgOpt,0.40*lambdaS,1,Rc);
    Pi_70(K) = MGmm_SD_Queue(K,KgOpt,0.70*lambdaS,1,Rc);
    Pi_100(K) = MGmm_SD_Queue(K,KgOpt,lambdaS,1,Rc);
end

%% User Probability Dist at 10,50,100 % loading

figure('Name','User Distribution at different cell loadings')
hold on; box on;
plot(1:KgOpt,Pi_10,'.-m',1:KgOpt,Pi_40,'.-b',1:KgOpt,Pi_70,'.-g',1:KgOpt,Pi_100,'.-r','linewidth',1.3);
ylim([0 0.225])
set(gca,'TickLength',[0 0])
xlabel('Number of Users ({\itK})')
ylabel('Probability')
legend('10% load','40% load','70% load','100% load', 'location', 'best')