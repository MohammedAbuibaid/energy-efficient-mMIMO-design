function [ EEglobalOpt, MglobalOpt, KglobalOpt, p ] = GlobalOptimum_EE(p_vec,Mmax,Kmax,PLO,PLI,plott)


EE = zeros(Mmax,Kmax);
Kopt_vec = zeros(1,length(p_vec));
Mopt_vec = zeros(1,length(p_vec));
EE_max_vec = zeros(1,length(p_vec));

for pInd = 1:length(p_vec)
    disp(['Current p = ' num2str(p_vec(pInd)) ',    p_max = ' num2str(p_vec(end))]);
    for M=1:Mmax
        for K=1:Kmax
            [EE(M,K),~,~] =  EE_R_Ptot_PA(PLO,PLI,Kmax,K,M,M*ones(1,18),p_vec(pInd),[]);
        end
    end
    [EEvalues,MInd]  = max(EE);
    [EE_max, Kopt] = max(EEvalues);
    Mopt = MInd(Kopt);
    Kopt_vec(pInd)=Kopt;
    Mopt_vec(pInd)=Mopt;
    EE_max_vec(pInd)=EE_max;
end

[EEglobalOpt, Indg] = max(EE_max_vec);
MglobalOpt = Mopt_vec(Indg);
KglobalOpt = Kopt_vec(Indg);
p = p_vec(Indg);








if plott ==true
    figure
    subplot(3,1,1);
    plot(p_vec,EE_max_vec)
    title('Max EE')
    hold on;grid on
    plot(p,EEglobalOpt, 'ro')
    subplot(3,1,2);
    plot(p_vec,Mopt_vec)
    title('Optimal M for max EE')
    hold on;grid on
    plot(p,MglobalOpt, 'ro')
    subplot(3,1,3);
    plot(p_vec,Kopt_vec)
    title('Optimal K for max EE')
    hold on;grid on
    plot(p,KglobalOpt, 'ro')
    xlabel('Average power per Antenna, p [W]')
    
    %%
%     clear all;close all;clc;
%     load('Ref_0_Environment.mat');
%     Mmax=300; Kmax=150; p = 0.0975;
%     [PLO,PLI] = Wrap_Around_PLO_PLI(BSLocations,UELocations,1,Scale,false);
%     

 
    % Illustrating Global Optimum point at EE, R, and Ptot.
    EE = zeros(Mmax,Kmax);
    R  = zeros(Mmax,Kmax);
    Ptot=zeros(Mmax,Kmax);
    for M=1:Mmax
%         disp(['Current M = ' num2str(M) '   Mmax = ' num2str(Mmax)]);
        for K=1:Kmax
            [EE(M,K),R(M,K),Ptot(M,K)] = EE_R_Ptot(Kmax,K,M,M*ones(1,18),PLO,PLI,p);
        end
    end
    
    [EEvalues,indM] = max(EE(1:200,1:100));
    [GlobalEE,indK] = max(EEvalues);
    GlobalM = indM(indK);
    GlobalK = indK;
    GlobalR = R(GlobalM,GlobalK);
    GlobalPtot = Ptot(GlobalM,GlobalK);
    
    
    %% Energy Efficiency [Mbit/Joule]
    figure
    hold on; box off; grid on;
    title('Energy Efficiency [Mbit/Joule]')
    gridDensity = 25;
    surface(1:Kmax,1:Mmax,EE/1e6,'EdgeColor','none'); %Plot the 3d surface
    colormap(autumn);
    %Plot lines on top of the 3d surface, to make it easier to see the shape
    for m = [1 gridDensity:gridDensity:Mmax]
        plot3(1:Kmax,m*ones(1,Kmax),EE(m,:)/1e6,'k-');
    end
    for k = [1 gridDensity:gridDensity:Kmax]
        plot3(k*ones(1,Mmax),1:Mmax,EE(:,k)/1e6,'k-');
    end
    %Compute and plot the optimal point
    
    plot3(GlobalK,GlobalM,GlobalEE/1e6,'b*','MarkerSize',13);
    txt = {'Global Optimum:',['M = ' num2str(GlobalM) ', K = ' num2str(GlobalK)],['EE = ' num2str(GlobalEE/1e6) ' Mbit/J']};
    text(150,200,GlobalEE/1e6,txt,'HorizontalAlignment','left','BackgroundColor', [1 1 1],'EdgeColor','b')
    view([-46 24]);
    axis([0 Kmax 0 Mmax]);
    xlabel('Number of Users (K)');
    ylabel('Number of Antennas (M)');
    zlabel('Energy Efficiency [Mbit/Joule]');
    
    % Average User Rate [Mbps]
    %%
    figure
    hold on; box off; grid on;
    title('Average User Rate [Mbps]')
    gridDensity = 25;
    surface(1:Kmax,1:Mmax,R/1e6,'EdgeColor','none'); %Plot the 3d surface
    colormap(autumn);
    %Plot lines on top of the 3d surface, to make it easier to see the shape
    for m = [1 gridDensity:gridDensity:Mmax]
        plot3(1:Kmax,m*ones(1,Kmax),R(m,:)/1e6,'k-');
    end
    for k = [1 gridDensity:gridDensity:Kmax]
        plot3(k*ones(1,Mmax),1:Mmax,R(:,k)/1e6,'k-');
    end
    %plot the optimal point
    plot3(GlobalK,GlobalM,GlobalR/1e6,'b*','MarkerSize',13);
    txtR = {'Global Optimum:',['M = ' num2str(GlobalM) ', K = ' num2str(GlobalK)],['R = ' num2str(GlobalR/1e6) ' Mbps']};
    text(100,150,GlobalR/1e6+30,txtR,'HorizontalAlignment','left','BackgroundColor', [1 1 1],'EdgeColor','b')
    view([-46 24]);
    axis([0 Kmax 0 Mmax]);
    ylabel('Number of Antennas (M)');
    xlabel('Number of Users (K)');
    zlabel('Average User Rate [Mbps]');
   
    % Total Power Consumption [Joule]
    %%
    figure
    hold on; box off; grid on;
    title('Total Power Consumption [Joule/s]')
    gridDensity = 25;
    surface(1:Kmax,1:Mmax,Ptot,'EdgeColor','none'); %Plot the 3d surface
    colormap(autumn);
    %Plot lines on top of the 3d surface, to make it easier to see the shape
    for m = [1 gridDensity:gridDensity:Mmax]
        plot3(1:Kmax,m*ones(1,Kmax),Ptot(m,:),'k-');
    end
    for k = [1 gridDensity:gridDensity:Kmax]
        plot3(k*ones(1,Mmax),1:Mmax,Ptot(:,k),'k-');
    end
    %plot the optimal point
    plot3(GlobalK,GlobalM,GlobalPtot,'b*','MarkerSize',13);
    txtP = {'Global Optimum:',['M = ' num2str(GlobalM) ', K = ' num2str(GlobalK)],['Ptot = ' num2str(GlobalPtot) ' Joule/s']};
    text(100,140,GlobalPtot+50,txtP,'HorizontalAlignment','left','BackgroundColor', [1 1 1],'EdgeColor','b')
    view([-46 24]);
    axis([0 Kmax 0 Mmax]);
    ylabel('Number of Antennas (M)');
    xlabel('Number of Users (K)');
    zlabel('Total Power Consumption [Joule/s]');
    
end
%%
end