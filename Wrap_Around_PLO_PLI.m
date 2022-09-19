function [PLO,PLI] = Wrap_Around_PLO_PLI(BSlocations,UElocations,BSj,Rmax,plott)

% % INPUTS:
% UElocations should be column vector of testpoints at the origin
% BSLoc : Location of the BS under study


%%
% %////////////////////////////////////////////////////////////////////////
% clear all;close all;clc;
% Rmin = 35; %Users inside this circle will not considered in simulations
% TestPoints = 15000;  % 15000, 3740, 71060
% %Test Points Generation for 19 Hexagonal mutli-cell senario
% Rmax = 500; %Cell radius (distance to a vertex of hexagonal cell)
% UElocations = HexCell_SquareGrid_TestPoints(Rmax,Rmin,TestPoints,false);
% Locations of 19 BSs in a Hexagonal Cellular Network 
% ISD = Rmax*sqrt(3);
% % % Coordinates of all BSs in the Hexagonal Network
% u = [0 1 0 -1 -1  0  1 2 1 0 -1 -2 -2 -2 -1  0  1  2  2]; % 30-Degree axis
% v = [0 0 1  1  0 -1 -1 0 1 2  2  2  1  0 -1 -2 -2 -2 -1]; % Vertical axis
% % Setting the BSs in their locations as seen from the origin.
% BSlocations = sqrt(3).*(ISD/2+1i*Rmax/2).*u + (0+1i*ISD).*v;
% plott = true;
% %////////////////////////////////////////////////////////////////////////


% Light Colors for plotting section
Sky_Blue=[0.5294 0.8078 0.9804];
Light_Salmon=[1.0000 0.6250 0.4766];
LightGreen = [0.5625 0.9297 0.5625];
COLORS = [Sky_Blue; Light_Salmon;LightGreen];

%Path-loss Model Parameters
kappa = 3.76;
dbar = 10^(-3.53);

% PLO Calculation
PLO = dbar./(abs(UElocations).^kappa);

% Wrap-Around Translation Points
Po1 = 3.0+1i*8*sqrt(3)/2;
Po2 = 4.5-1i*7*sqrt(3)/2;
Po3 = 7.5 + 1i*sqrt(3)/2;
Replica = Rmax*[0; Po1; -Po1; Po2; -Po2; Po3; -Po3];

%Calculate the interference-distance to all BS in the Network (except serving BS) using Wrap-Around Model
BSLoc = BSlocations (BSj);
minD = NaN*ones(length(UElocations),19);
InteBSLoc = NaN*ones(length(UElocations),19);
for Ind = 1:19%length(BSlocations)
    if BSlocations(Ind)~=BSLoc
        D = zeros(length(UElocations),7);
        for pt = 1:7
            D(:,pt)=abs((UElocations+BSLoc)-(BSlocations(Ind)+Replica(pt)));
        end
        [minD(:,Ind), minD_Pt] = min(D,[],2);
        %         for tp = 1:length(UElocations)
        %         InteBSLoc(tp,Ind) = BSlocations(Ind) + Replica(minD_Pt(tp));
        %         end
        InteBSLoc(:,Ind) = BSlocations(Ind) + Replica(minD_Pt);
        
    end
end

% PLI Calculation
minD(:,isnan(minD(1,:)))=[];
PLI = dbar./(minD).^kappa;


if plott == true
    for InteBSInd = 1:19;
        if ~isnan(InteBSLoc(1,InteBSInd))
            figure
            InteBSLoc_temp = unique(InteBSLoc(:,InteBSInd),'stable');
            for ind2 = 1:length(InteBSLoc_temp)
                if length(InteBSLoc_temp)<=3
                    %Interfering BSs
                    Vplot=InteBSLoc_temp(ind2);
                    plot(real(Vplot),imag(Vplot),'.','MarkerSize',50,'color',COLORS(ind2,:));
                    hold on
                    % Cell under study
                    Vplot=UElocations(InteBSLoc(:,InteBSInd)==InteBSLoc_temp(ind2))+BSLoc;
                    plot(real(Vplot),imag(Vplot),'.','color',COLORS(ind2,:));
                    hold on
                else
                    error('Error: Need more Colors');
                end
            end
            
            for rep_ind =1:7
                Vplot = BSlocations +(Replica(rep_ind));
                plot(real(Vplot),imag(Vplot),'ko','MarkerSize', 15);
                hold on
                for Ind = 1:length(BSlocations)
                    x1 = real(BSlocations(Ind)+(Replica(rep_ind)));
                    y1 = imag(BSlocations(Ind)+(Replica(rep_ind)));
                    text(x1,y1,num2str(Ind),'HorizontalAlignment','Center');
                    hold on
                end
            end
        end
    end
    
end
%%
end % function


% if plott == true
%     
%     figure % Plotting 19 Hexagonal Cells with different colors
%     %Generate Distinguishable Colors
%     colors = distinguishable_colors(19*2);
%     colorIndex =1;
%     for j = 1:1:19
%         plot(UElocations+BSlocations(j),'o','color',colors(colorIndex,:));
%         hold on;
%         colorIndex = colorIndex +2;
%     end
%     grid on;
%     title('Multi-cell simulation scenario')
%     xlabel('Distance m') % x-axis label
%     ylabel('Distance m') % y-axis label
% end

% 
% 
%     figure % Plotting Serving Cell
%     plot(UElocations,'.','linewidth',1);
%     title('Uniform UE distribution withn a cell')
%     xlabel('Distance m') % x-axis label
%     ylabel('Distance m') % y-axis label
%     
%     figure % Plotting 19 Hexagonal Cells with different colors
%     %Generate Distinguishable Colors
%     colors = distinguishable_colors(25*2);
%     colorIndex =1;
%     for j = 1:1:25
%         plot(UElocations+BSlocations(j),'.','color',colors(colorIndex,:));
%         hold on
%         colorIndex = colorIndex +2;
%     end
%     title('Multi-cell simulation scenario')
%     xlabel('Distance m') % x-axis label
%     ylabel('Distance m') % y-axis label
% end