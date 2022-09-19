close all; clear all; clc;

%% 19-Cell Hexagonal Cluster
Rmax = 500; %Cell radius (distance to a vertex of hexagonal cell)
ISD = Rmax*sqrt(3);% Inter-Site-Distance
Rmin = 35; %Users inside this circle will not considered in simulations
%Coordinates of all BSs
u = [0 1 0 -1 -1  0  1  2 2 1 0 -1 -2 -2 -2 -1  0  1  2]; % 30-Degree axis
v = [0 0 1  1  0 -1 -1 -1 0 1 2  2  2  1  0 -1 -2 -2 -2]; % Vertical axis
%Placing the BSs around the cell in the origin.
BSLocations = sqrt(3).*(ISD/2+1i*Rmax/2).*u + (0+1i*ISD).*v;


%% Generating TestPoints distributed uniformaly over the cell coerage area
TestPoints = 10000;
Ri = [Rmin Rmax];
ULD_Uni = 1 ;
UELocations = UE_insertion_MonteCarlo_HexCell(TestPoints,ULD_Uni,Ri,false);
[PLO,PLI] = Wrap_Around_PLO_PLI(BSLocations,UELocations,1,Rmax,false);


%% Searching for the Global Optimum triplet (MgOpt, KgOpt, pgOpt)
%Mmax = 300; Kmax=150; p_vec = 0.06:0.0125:0.1475;
%[EEglobalOpt,MgOpt,KgOpt,pgOpt] = GlobalOptimum_EE(p_vec,Mmax,Kmax,PLO,PLI,true);
%disp(['pgOpt = ' num2str(pgOpt), ' ,  KgOpt = ' num2str(KgOpt) ' ,  MgOpt = ' num2str(MgOpt)])
pgOpt = 0.0725 ;  KgOpt = 100 ;  MgOpt = 198;


%% ^^^^^^^^ Section IV-A: Energy Efficiency vs. User Rate Trade-off ^^^^^^^^
EE = zeros(MgOpt,KgOpt);
Rc = zeros(MgOpt,KgOpt);
Ptot = zeros(MgOpt,KgOpt);

for K=1:KgOpt
	for M=1:MgOpt
		Mc = M; Md = M*ones(1,18);
		[EE(M,K),Rc(M,K),Ptot(M,K)] = EE_R_Ptot_PA(PLO,PLI,KgOpt,K,Mc,Md,pgOpt,[]);
	end
end
 
%%
%% Plot EE 3D 
figure('Name','EE 3D Plot');
hold on; box on; grid on;
gridDensity = 20;
surface(1:KgOpt,1:MgOpt,EE/1e6,'EdgeColor','none'); %Plot the 3d surface
ColorMap = [0.9769 0.9839 0.0805;0.97212299187524 0.969621566540854 0.0963678477832249;0.967849009059724 0.955336600285485 0.110275856699262;0.964361468914577 0.940951409544263 0.122439517634256;0.961918415614537 0.926405523191179 0.133006049292428;0.960768487262726 0.911671936089652 0.142037723677499;0.961257726206098 0.896719881475539 0.149788585027059;0.963667834445408 0.881536394114721 0.156789064180283;0.967944593593052 0.866167883409773 0.163653529793876;0.973601256311507 0.850716461245219 0.170774647540138;0.979914978799366 0.835277871249971 0.178171603139242;0.98610236930393 0.819926811881292 0.185746539215814;0.991146055584862 0.804785179392497 0.193808883739653;0.994360581102825 0.789981013643174 0.202912818345569;0.995893880208333 0.775548779832749 0.212985804675874;0.995620091074989 0.761673116307032 0.223041777674357;0.992014391290574 0.74912330758231 0.231234020850772;0.983385379922958 0.738942957142421 0.235182304772877;0.969235727083115 0.732158587619236 0.232018959990002;0.950619342163631 0.729374483013153 0.219904550824846;0.929448589888073 0.730133107921055 0.201234700738816;0.907126160780589 0.732949683012281 0.182928962035406;0.883277173242115 0.736991673142569 0.169535965438116;0.858138956283388 0.741888382066999 0.160773694619678;0.831826359013149 0.747383108084542 0.156779965427944;0.804258465099335 0.753271756335667 0.158801735510145;0.775389459364755 0.759348915263585 0.167246756989615;0.74532498028619 0.765434276199341 0.180884940583365;0.714233823485602 0.771385402497791 0.198022207605271;0.682249372491382 0.777092750599271 0.21761453517278;0.649465922110421 0.782459189305987 0.2393410981587;0.61598338400977 0.787373647367387 0.263128171970731;0.581916234915597 0.791714353579567 0.288668990798224;0.547272939722878 0.795321201265426 0.315260912545522;0.512107476479667 0.798279137120928 0.341872973932539;0.476595661245074 0.800559236839839 0.368408408260345;0.441526878247942 0.801984348969232 0.395111759058634;0.407632272693089 0.802498882034847 0.42210967672893;0.375337330354963 0.802058169582912 0.44906389062064;0.344250199753898 0.800743684237344 0.475282361330305;0.313883032444545 0.798738102613177 0.50031673954555;0.285016038894653 0.79608519852502 0.524372908455985;0.258792765017918 0.792765104947771 0.547659205109732;0.236247613375527 0.788755297170367 0.570130869334085;0.217907648195539 0.78411284626552 0.591640941456386;0.202725832880111 0.77904852545602 0.61210716280256;0.188755620783851 0.773765713391985 0.631569319261823;0.173843109707605 0.768429277065822 0.65020089931488;0.156179502868652 0.763131499808175 0.668200321306501;0.13475427688417 0.757891194579715 0.685748499089196;0.109269114848546 0.752658963285174 0.702976742281233;0.0805105970337277 0.74732697666259 0.719948253472646;0.0515038315727597 0.741759554817563 0.736671324829828;0.0272504555656797 0.73585723069963 0.753115878600166;0.0126249648684547 0.72955817145393 0.769231759325663;0.00889499538966589 0.722874210126059 0.785012271595001;0.0154244401840937 0.715801874687558 0.800437088049026;0.0308687896501451 0.708305377301716 0.815432985192253;0.0525218651635306 0.70029336479732 0.829809122494289;0.0731053780828203 0.691657078674861 0.843234235422952;0.0888660669054304 0.682422207042149 0.855457711846488;0.1011896430833 0.672696663434165 0.866363348702022;0.11121659937359 0.662519426046099 0.875854506610689;0.119169865635463 0.653057289777483 0.883191512625558];
colormap(ColorMap);
c = colorbar('FontName','Times New Roman','FontSize',14);
c.Label.String = 'Energy Efficiency [Mbit/Joule]';
%Plot lines on top of the 3d surface, to make it easier to see the shape
for m = [1 gridDensity:gridDensity:MgOpt MgOpt]
	plot3(1:KgOpt,m*ones(1,KgOpt),EE(m,:)/1e6,'k-');
end
for k = [1 gridDensity:gridDensity:KgOpt]
	plot3(k*ones(1,MgOpt),1:MgOpt,EE(:,k)/1e6,'k-');
end
view([-25 50]);
axis([0 KgOpt 0 MgOpt 0 15]);
xlabel('Number of Users ({\itK})','Rotation',14,'FontName','Times New Roman','FontSize',14,'Position',[47,-25,-1.5]);
ylabel('Number of Antennas ({\itM})','Rotation',-50,'FontName','Times New Roman','FontSize',14,'Position',[-15,101.5,-0.5]);
zlabel({'Energy Efficiency';'[Mbit/Joule]'},'FontName','Times New Roman','FontSize',14,'Position',[-5,210,8]);
set(gca,'FontName','Times New Roman','FontSize',14);

% %% Save Fig as pdf
% h = gcf;
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'Fig.01.A.pdf','-dpdf')



%% Plot Rc 3D
figure('Name','Rc 3D Plot');
hold on; box on; grid on;
gridDensity = 20;
surface(1:KgOpt,1:MgOpt,Rc/1e6,'EdgeColor','none'); %Plot the 3d surface
ColorMap = [0.9769 0.9839 0.0805;0.97212299187524 0.969621566540854 0.0963678477832249;0.967849009059724 0.955336600285485 0.110275856699262;0.964361468914577 0.940951409544263 0.122439517634256;0.961918415614537 0.926405523191179 0.133006049292428;0.960768487262726 0.911671936089652 0.142037723677499;0.961257726206098 0.896719881475539 0.149788585027059;0.963667834445408 0.881536394114721 0.156789064180283;0.967944593593052 0.866167883409773 0.163653529793876;0.973601256311507 0.850716461245219 0.170774647540138;0.979914978799366 0.835277871249971 0.178171603139242;0.98610236930393 0.819926811881292 0.185746539215814;0.991146055584862 0.804785179392497 0.193808883739653;0.994360581102825 0.789981013643174 0.202912818345569;0.995893880208333 0.775548779832749 0.212985804675874;0.995620091074989 0.761673116307032 0.223041777674357;0.992014391290574 0.74912330758231 0.231234020850772;0.983385379922958 0.738942957142421 0.235182304772877;0.969235727083115 0.732158587619236 0.232018959990002;0.950619342163631 0.729374483013153 0.219904550824846;0.929448589888073 0.730133107921055 0.201234700738816;0.907126160780589 0.732949683012281 0.182928962035406;0.883277173242115 0.736991673142569 0.169535965438116;0.858138956283388 0.741888382066999 0.160773694619678;0.831826359013149 0.747383108084542 0.156779965427944;0.804258465099335 0.753271756335667 0.158801735510145;0.775389459364755 0.759348915263585 0.167246756989615;0.74532498028619 0.765434276199341 0.180884940583365;0.714233823485602 0.771385402497791 0.198022207605271;0.682249372491382 0.777092750599271 0.21761453517278;0.649465922110421 0.782459189305987 0.2393410981587;0.61598338400977 0.787373647367387 0.263128171970731;0.581916234915597 0.791714353579567 0.288668990798224;0.547272939722878 0.795321201265426 0.315260912545522;0.512107476479667 0.798279137120928 0.341872973932539;0.476595661245074 0.800559236839839 0.368408408260345;0.441526878247942 0.801984348969232 0.395111759058634;0.407632272693089 0.802498882034847 0.42210967672893;0.375337330354963 0.802058169582912 0.44906389062064;0.344250199753898 0.800743684237344 0.475282361330305;0.313883032444545 0.798738102613177 0.50031673954555;0.285016038894653 0.79608519852502 0.524372908455985;0.258792765017918 0.792765104947771 0.547659205109732;0.236247613375527 0.788755297170367 0.570130869334085;0.217907648195539 0.78411284626552 0.591640941456386;0.202725832880111 0.77904852545602 0.61210716280256;0.188755620783851 0.773765713391985 0.631569319261823;0.173843109707605 0.768429277065822 0.65020089931488;0.156179502868652 0.763131499808175 0.668200321306501;0.13475427688417 0.757891194579715 0.685748499089196;0.109269114848546 0.752658963285174 0.702976742281233;0.0805105970337277 0.74732697666259 0.719948253472646;0.0515038315727597 0.741759554817563 0.736671324829828;0.0272504555656797 0.73585723069963 0.753115878600166;0.0126249648684547 0.72955817145393 0.769231759325663;0.00889499538966589 0.722874210126059 0.785012271595001;0.0154244401840937 0.715801874687558 0.800437088049026;0.0308687896501451 0.708305377301716 0.815432985192253;0.0525218651635306 0.70029336479732 0.829809122494289;0.0731053780828203 0.691657078674861 0.843234235422952;0.0888660669054304 0.682422207042149 0.855457711846488;0.1011896430833 0.672696663434165 0.866363348702022;0.11121659937359 0.662519426046099 0.875854506610689;0.119169865635463 0.653057289777483 0.883191512625558];
colormap(ColorMap);
c = colorbar('FontName','Times New Roman','FontSize',14);
c.Label.String = 'User Data Rate [Mbps]';
%Plot lines on top of the 3d surface, to make it easier to see the shape
for m = [1 gridDensity:gridDensity:MgOpt MgOpt]
	plot3(1:KgOpt,m*ones(1,KgOpt),Rc(m,:)/1e6,'k-');
	end
for k = [1 gridDensity:gridDensity:KgOpt]
	plot3(k*ones(1,MgOpt),1:MgOpt,Rc(:,k)/1e6,'k-');
end
view([-25 50]);
axis([0 KgOpt 0 MgOpt 0 200]);
xlabel('Number of Users ({\itK})','Rotation',14,'FontName','Times New Roman','FontSize',14,'Position',[47,-25,-1.5]);
ylabel('Number of Antennas ({\itM})','Rotation',-43,'FontName','Times New Roman','FontSize',14,'Position',[-15,101.5,-0.5]);
zlabel('User Data Rate [Mbps]','FontName','Times New Roman','FontSize',14,'Position',[-5,210,8]);
set(gca,'FontName','Times New Roman','FontSize',14);

% %% Save Fig as pdf
% h = gcf;
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'Fig.01.B.pdf','-dpdf')


%% Plot Ptot 3D
% figure('Name','Ptot 3D Plot');
% hold on; box on; grid on;
% gridDensity = 20;
% surface(1:KgOpt,1:MgOpt,Ptot,'EdgeColor','none'); %Plot the 3d surface
% ColorMap = [0.9769 0.9839 0.0805;0.97212299187524 0.969621566540854 0.0963678477832249;0.967849009059724 0.955336600285485 0.110275856699262;0.964361468914577 0.940951409544263 0.122439517634256;0.961918415614537 0.926405523191179 0.133006049292428;0.960768487262726 0.911671936089652 0.142037723677499;0.961257726206098 0.896719881475539 0.149788585027059;0.963667834445408 0.881536394114721 0.156789064180283;0.967944593593052 0.866167883409773 0.163653529793876;0.973601256311507 0.850716461245219 0.170774647540138;0.979914978799366 0.835277871249971 0.178171603139242;0.98610236930393 0.819926811881292 0.185746539215814;0.991146055584862 0.804785179392497 0.193808883739653;0.994360581102825 0.789981013643174 0.202912818345569;0.995893880208333 0.775548779832749 0.212985804675874;0.995620091074989 0.761673116307032 0.223041777674357;0.992014391290574 0.74912330758231 0.231234020850772;0.983385379922958 0.738942957142421 0.235182304772877;0.969235727083115 0.732158587619236 0.232018959990002;0.950619342163631 0.729374483013153 0.219904550824846;0.929448589888073 0.730133107921055 0.201234700738816;0.907126160780589 0.732949683012281 0.182928962035406;0.883277173242115 0.736991673142569 0.169535965438116;0.858138956283388 0.741888382066999 0.160773694619678;0.831826359013149 0.747383108084542 0.156779965427944;0.804258465099335 0.753271756335667 0.158801735510145;0.775389459364755 0.759348915263585 0.167246756989615;0.74532498028619 0.765434276199341 0.180884940583365;0.714233823485602 0.771385402497791 0.198022207605271;0.682249372491382 0.777092750599271 0.21761453517278;0.649465922110421 0.782459189305987 0.2393410981587;0.61598338400977 0.787373647367387 0.263128171970731;0.581916234915597 0.791714353579567 0.288668990798224;0.547272939722878 0.795321201265426 0.315260912545522;0.512107476479667 0.798279137120928 0.341872973932539;0.476595661245074 0.800559236839839 0.368408408260345;0.441526878247942 0.801984348969232 0.395111759058634;0.407632272693089 0.802498882034847 0.42210967672893;0.375337330354963 0.802058169582912 0.44906389062064;0.344250199753898 0.800743684237344 0.475282361330305;0.313883032444545 0.798738102613177 0.50031673954555;0.285016038894653 0.79608519852502 0.524372908455985;0.258792765017918 0.792765104947771 0.547659205109732;0.236247613375527 0.788755297170367 0.570130869334085;0.217907648195539 0.78411284626552 0.591640941456386;0.202725832880111 0.77904852545602 0.61210716280256;0.188755620783851 0.773765713391985 0.631569319261823;0.173843109707605 0.768429277065822 0.65020089931488;0.156179502868652 0.763131499808175 0.668200321306501;0.13475427688417 0.757891194579715 0.685748499089196;0.109269114848546 0.752658963285174 0.702976742281233;0.0805105970337277 0.74732697666259 0.719948253472646;0.0515038315727597 0.741759554817563 0.736671324829828;0.0272504555656797 0.73585723069963 0.753115878600166;0.0126249648684547 0.72955817145393 0.769231759325663;0.00889499538966589 0.722874210126059 0.785012271595001;0.0154244401840937 0.715801874687558 0.800437088049026;0.0308687896501451 0.708305377301716 0.815432985192253;0.0525218651635306 0.70029336479732 0.829809122494289;0.0731053780828203 0.691657078674861 0.843234235422952;0.0888660669054304 0.682422207042149 0.855457711846488;0.1011896430833 0.672696663434165 0.866363348702022;0.11121659937359 0.662519426046099 0.875854506610689;0.119169865635463 0.653057289777483 0.883191512625558];
% colormap(ColorMap);
% %Plot lines on top of the 3d surface, to make it easier to see the shape
% for m = [1 gridDensity:gridDensity:MgOpt MgOpt]
% 	plot3(1:KgOpt,m*ones(1,KgOpt),Ptot(m,:),'k-');
% 	end
% for k = [1 gridDensity:gridDensity:KgOpt]
% 	plot3(k*ones(1,MgOpt),1:MgOpt,Ptot(:,k),'k-');
% end
% view([-25 50]);
% axis([0 KgOpt 0 MgOpt 0 350]);
% xlabel('Number of Users ({\itK})','Rotation',10);
% ylabel('Number of Antennas ({\itM})','Rotation',-50);
% zlabel('Total Power [W]');

%% Save Fig as pdf
% h = gcf;
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'Ptot_3D','-dpdf')


%% ^^^^^^^^^^^^^^^^^^^^^^^^ End of Section I ^^^^^^^^^^^^^^^^^^^^^^^^


%% Traffic model: User Profile
eMBB = [" 40%", " 60%", " 80%"];
BS_GDR = [69 95 120]*1e6;


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

%save('UserPDF.mat','KgOpt','Rc','lambdaS');


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

%save('UserDist.mat','KgOpt','Pi_10','Pi_40','Pi_70','Pi_100');

%% Probability Dist Plot
load('UserDist.mat')
figure('Name','User Distribution at different cell loadings')
hold on; box on; grid on
plot(1:KgOpt,Pi_10,'.-m',1:KgOpt,Pi_40,'.-b',1:KgOpt,Pi_70,'.-g',1:KgOpt,Pi_100,'.-r','linewidth',1.3);
ylim([0 0.225])
set(gca,'TickLength',[0 0])
xlabel('Number of Users ({\itK})','FontName','Times New Roman','FontSize',14)
ylabel('Probability','FontName','Times New Roman','FontSize',14)
legend('10% loading','40% loading','70% loading','100% loading', 'location', 'best','FontName','Times New Roman','FontSize',14)
set(gca,'FontName','Times New Roman','FontSize',14);
% %% Save Fig as pdf
% h = gcf;
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'UserDist','-dpdf')

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
% fprintf('\n')
% fprintf('\n')
% disp('$$$$$$$$$$$$$$ EE Optimization Process $$$$$$$$$$$$$$');
% Mavg_EEopt = zeros(1,length(loading));
% EEopt = zeros(1,length(loading));
% Ravg_EEopt = zeros(1,length(loading));
% Mopt_EEopt= cell(1,length(loading));
% for L = 1:length(loading)
%     disp(['IIIIIIII Loading = ' num2str(loading(L)) '% IIIIIIII'])
%     [Mopt_EEopt{1,L},Mavg_EEopt(1,L),EEopt(1,L),Ravg_EEopt(1,L)]= Mavg_EE_Optimizer(PLO,PLI,KgOpt,MgOpt,pgOpt,loading(L),lambdaS);    
% end    
% 
% save('Mavg_EE_Optimizer.mat','Mopt_EEopt','EEopt','Ravg_EEopt','Mavg_EEopt','loading','KgOpt','MgOpt','pgOpt','lambdaS','PLO','PLI');





%% Plot Mavg for GDR and for Opt EE
load('Section_IV_B.mat')
figure('Name','Optimum Number of Antennas for GDR and Opt EE');
hold on;  box on, grid on
plot(loading,MgOpt*ones(1,length(loading)),'-.k','linewidth',1.3)
plot(loading,Mavg_GDR(3,:),'-.+m','linewidth',1.3)
plot(loading,Mavg_GDR(2,:),'-.*b','linewidth',1.3)
plot(loading,Mavg_GDR(1,:),'-.og','linewidth',1.3)
load('Mavg_EE_Optimizer.mat');
plot(loading,Mavg_EEopt,'-.xr','linewidth',1.3)
xlim([0 100])
ylim([0 220])
ylabel('Average Number of Antennas','FontName','Times New Roman','FontSize',14);
xlabel('Network load [%]','FontName','Times New Roman','FontSize',14);
legend('Fixed Antenna System', strcat('\rho_1= 80 %',' , GDR = 120 Mbps'), strcat('\rho_1= 60 %',' , GDR = 95 Mbps'),strcat('\rho_1= 40 %', ', GDR = 69 Mbps'),'Optimum EE Algorithm [7]','location','southeast')
set(gca,'FontName','Times New Roman','FontSize',14);

% %% Save Fig as pdf
% h = gcf;
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'Fig02','-dpdf')




%% Plot EE for GDR and for Opt EE
%load('Section_IV_B.mat')
figure('Name','EE for GDR and OptEE');
hold on;  box on, grid on
plot(loading,EE_Mmax_avg/1e6,'-.k','linewidth',1.3)
plot(loading,EEavg_GDR(3,:)/1e6,'-.+m','linewidth',1.3)
plot(loading,EEavg_GDR(2,:)/1e6,'-.*b','linewidth',1.3)
plot(loading,EEavg_GDR(1,:)/1e6,'-.og','linewidth',1.3)
%load('Mavg_EE_Optimizer.mat');
plot(loading,EEopt/1e6,'-.xr','linewidth',1.3)
xlim([0 100])
ylim([0 15])
ylabel('Energy Efficiency [Mbit/Joule]','FontName','Times New Roman','FontSize',14);
xlabel('Network load [%]','FontName','Times New Roman','FontSize',14);
legend('Fixed Antenna System', strcat('\rho_1= 80 %',' , GDR = 120 Mbps'), strcat('\rho_1= 60 %',' , GDR = 95 Mbps'),strcat('\rho_1= 40 %', ', GDR = 69 Mbps'),'Optimum EE Algorithm [7]','location','southeast')
set(gca,'FontName','Times New Roman','FontSize',14);


% %% Save Fig as pdf
% h = gcf;
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'Fig03','-dpdf')

%% Plot Ropt for GDR and for Opt EE
%load('Section_IV_B.mat')
figure('Name','EE for GDR and OptEE');
hold on;  box on, grid on
plot(loading,Rc_Mmax_avg/1e6,'-.k','linewidth',1.3)
plot(loading,Ravg_GDR(3,:)/1e6,'-.+m','linewidth',1.3)
plot(loading,Ravg_GDR(2,:)/1e6,'-.*b','linewidth',1.3)
plot(loading,Ravg_GDR(1,:)/1e6,'-.og','linewidth',1.3)
%load('Mavg_EE_Optimizer.mat');
plot(loading,Ravg_EEopt/1e6,'-.xr','linewidth',1.3)
xlim([0 100])
% ylim([0 15])
ylabel('User Data Rate [Mbps]','FontName','Times New Roman','FontSize',14);
xlabel('Network load [%]','FontName','Times New Roman','FontSize',14);
legend('Fixed Antenna System', strcat('\rho_1= 80 %',' , GDR = 120 Mbps'), strcat('\rho_1= 60 %',' , GDR = 95 Mbps'),strcat('\rho_1= 40 %', ', GDR = 69 Mbps'),'Optimum EE Algorithm [7]','location','northeast')
set(gca,'FontName','Times New Roman','FontSize',14);

% %% Save Fig as pdf
% h = gcf;
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'Fig04','-dpdf')
