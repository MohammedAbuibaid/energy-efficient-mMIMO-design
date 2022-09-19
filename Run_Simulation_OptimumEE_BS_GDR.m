% Copyright: Mohammed Abuibaid
% Carleton University
% Sep 19, 2022

close all; clear all; clc;
addpath(genpath('functions'))


%% Initializers: TestPoints, BS, ...

% 19-Cell Hexagonal Cluster
Rmax = 500; %Cell radius (distance to a vertex of hexagonal cell)
ISD = Rmax*sqrt(3);% Inter-Site-Distance
Rmin = 35; %Users inside this circle will not considered in simulations

%Coordinates of all BSs
u = [0 1 0 -1 -1  0  1  2 2 1 0 -1 -2 -2 -2 -1  0  1  2]; % 30-Degree axis
v = [0 0 1  1  0 -1 -1 -1 0 1 2  2  2  1  0 -1 -2 -2 -2]; % Vertical axis

%Placing the BSs around the cell in the origin.
BSLocations = sqrt(3).*(ISD/2+1i*Rmax/2).*u + (0+1i*ISD).*v;

% Generating TestPoints distributed uniformaly over the cell coerage area
TestPoints = 10000;
Ri = [Rmin Rmax];
ULD_Uni = 1 ;
UELocations = UE_insertion_MonteCarlo_HexCell(TestPoints,ULD_Uni,Ri,false);
[PLO,PLI] = Wrap_Around_PLO_PLI(BSLocations,UELocations,1,Rmax,false);

% % Searching for the Global Optimum triplet (MgOpt, KgOpt, pgOpt)
Mmax = 300; Kmax=150; p_vec = 0.06:0.0125:0.1475;
% % [EEglobalOpt,MgOpt,KgOpt,pgOpt] = GlobalOptimum_EE(p_vec,Mmax,Kmax,PLO,PLI,true);
% load Golbal Optimum mMIMO design variables (can be found by uncommeting above line GlobalOptimum_EE..m)
pgOpt = 0.0725 ;  KgOpt = 100 ;  MgOpt = 198; 
disp(['pgOpt = ' num2str(pgOpt), ' ,  KgOpt = ' num2str(KgOpt) ' ,  MgOpt = ' num2str(MgOpt)])

% % Traffic model: User Profile
eMBB = [" 40%", " 60%", " 80%"];
BS_GDR = [69 95 120]*1e6;

% % Traffic model: M/G/m/m state-dependent queue
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
%save('UserPDF.mat','KgOpt','Rc','lambdaS');




%% ^^^^^^^^ Section IV-B: Optimum Number of Antennas for Guaranteed BS Average User Rate ^^^^^^^^
Mavg_GDR = zeros(length(BS_GDR),length(loading));
EEavg_GDR = zeros(length(BS_GDR),length(loading));
Ravg_GDR = zeros(length(BS_GDR),length(loading));
Mopt_GDR= cell(length(BS_GDR),length(loading));

EE_Mmax_avg= zeros(1,length(loading));
Rc_Mmax_avg= zeros(1,length(loading));

disp('XXXXXXXXXXXXXX EE & Rc Joint Optimization Process XXXXXXXXXXXXXX');
for L = 1:length(loading)
    %Avg User Rate @ Mmax
    fprintf('\n')
    [EE_Mmax_avg(1,L),Rc_Mmax_avg(1,L)] = EE_R_QoS_Mmax(PLO,PLI,KgOpt,MgOpt,pgOpt,loading(L),lambdaS,[]);
    for r = 1:length(BS_GDR)
        fprintf('\n')
        disp(['++++++++  Loading = ' num2str(loading(L)) ' % ++++++++ GDR = ' num2str(BS_GDR(r)/1e6) ' Mbps'])
        if(Rc_Mmax_avg(1,L) > BS_GDR(r))
            %Calling Mavg_JOptimizer to tailor number of antennas to given GDR
            [Mopt_GDR{r,L},Mavg_GDR(r,L),EEavg_GDR(r,L),Ravg_GDR(r,L)]= Mavg_JOptimizer(PLO,PLI,KgOpt,MgOpt,pgOpt,loading(L ),lambdaS,BS_GDR(r));
        else
            Mopt_GDR{r,L} = MgOpt*ones(1,KgOpt);
            Mavg_GDR(r,L) = MgOpt;
            EEavg_GDR(r,L) = EE_Mmax_avg(1,L);
            Ravg_GDR(r,L) = Rc_Mmax_avg(1,L);
        end
    end    
end
% save('Section_IV_B.mat','Mopt_GDR','EEavg_GDR','Ravg_GDR','Mavg_GDR','EE_Mmax_avg','Rc_Mmax_avg','loading','BS_GDR','KgOpt','MgOpt','pgOpt','lambdaS','PLO','PLI');




%% EE Optimization Process
fprintf('\n')
fprintf('\n')
disp('$$$$$$$$$$$$$$ EE Optimization Process $$$$$$$$$$$$$$');
Mavg_EEopt = zeros(1,length(loading));
EEopt = zeros(1,length(loading));
Ravg_EEopt = zeros(1,length(loading));
Mopt_EEopt= cell(1,length(loading));
for L = 1:length(loading)
    disp(['IIIIIIII Loading = ' num2str(loading(L)) '% IIIIIIII'])
    [Mopt_EEopt{1,L},Mavg_EEopt(1,L),EEopt(1,L),Ravg_EEopt(1,L)]= Mavg_EE_Optimizer(PLO,PLI,KgOpt,MgOpt,pgOpt,loading(L),lambdaS);    
end    

% % save('Mavg_EE_Optimizer.mat','Mopt_EEopt','EEopt','Ravg_EEopt','Mavg_EEopt','loading','KgOpt','MgOpt','pgOpt','lambdaS','PLO','PLI');



