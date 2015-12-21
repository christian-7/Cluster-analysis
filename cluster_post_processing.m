%% Load data from txt files

% 1 number of cluster
% 2 area from convex hull
% 3 area from delaunay triangle
% 4 density as mol/area
% 5 mean diameter
% 6 stdev diameter
% 7 radave from HK
% 8 aspratio from HK
% 9 clustarea of diamter

MDCK_EGFR=struct('Unchanged', [], 'Starved' , [], 'StarvedGM', []);
MDCK_EGFR.Unchanged=dlmread('EGFR_unchanged_MDCK.txt');
MDCK_EGFR.Starved=dlmread('EGFR_starved_MDCK.txt');
MDCK_EGFR.StarvedGM=dlmread('EGFR_starved_GM_MDCK.txt');

A549_EGFR=struct('Unchanged', [], 'Starved' , [], 'StarvedGM', []);
A549_EGFR.Unchanged=dlmread('EGFR_unchanged.txt');
A549_EGFR.Starved=dlmread('EGFR_starved.txt');
A549_EGFR.StarvedGM=dlmread('EGFR_ST_GM.txt');

%% Compare Area/Density between conditions

% MDCK EGFR Area

bins=0:0.001:0.02;

f=transpose(hist(nonzeros(MDCK_EGFR.Unchanged(:,2)),bins)); 
f2=transpose(hist(nonzeros(MDCK_EGFR.Starved(:,2)),bins)); 
f3=transpose(hist(nonzeros(MDCK_EGFR.StarvedGM(:,2)),bins)); 
f4=[f/sum(f) f2/sum(f2) f3/sum(f3)];

figure('Position',[200 300 800 600])
set(gcf,'numbertitle','off','name','Compare Area/Density EGFR','PaperOrientation','landscape') % Title of the figure

subplot(2,2,1)
bar(bins,f4,1,'grouped');
title('Hist of Size MDCK EGFR');
xlabel('Area (?m^2)');
ylabel('norm counts');
axis([-0.001 0.02 0 0.3]);
leg1=legend('Unchanged','Starved','StarvedGM');
set(leg1,'FontSize',12);

clear x x2 x3 f f2 f3 f4 bins;

% MDCK EGFR Density

bins=0:5e3:1e5;

f=transpose(hist(nonzeros(MDCK_EGFR.Unchanged(:,4)),bins)); 
f2=transpose(hist(nonzeros(MDCK_EGFR.Starved(:,4)),bins)); 
f3=transpose(hist(nonzeros(MDCK_EGFR.StarvedGM(:,4)),bins)); 
f4=[f/sum(f) f2/sum(f2) f3/sum(f3)];

subplot(2,2,2)
bar(bins,f4,1,'grouped');
title('Hist of Density MDCK EGFR');
xlabel('Density (mol/?m^2)');
ylabel('norm counts');
axis([0 1e5 0 0.25])
leg2=legend('Unchanged','Starved','StarvedGM');

set(leg2,'FontSize',12);


% MDCK A549 Area

bins=0:0.001:0.02;

f=transpose(hist(nonzeros(A549_EGFR.Unchanged(:,2)),bins)); 
f2=transpose(hist(nonzeros(A549_EGFR.Starved(:,2)),bins)); 
f3=transpose(hist(nonzeros(A549_EGFR.StarvedGM(:,2)),bins)); 
f4=[f/sum(f) f2/sum(f2) f3/sum(f3)];

subplot(2,2,3)
bar(bins,f4,1,'grouped');
title('Hist of Size A549 EGFR');
xlabel('Area (?m^2)');
ylabel('norm counts');
axis([-0.001 0.02 0 0.3])
leg1=legend('Unchanged','Starved','StarvedGM');
set(leg1,'FontSize',12);

clear x x2 x3 f f2 f3 f4 bins;

% A549 EGFR Density

bins=0:5e3:1e5;

f=transpose(hist(nonzeros(A549_EGFR.Unchanged(:,4)),bins)); 
f2=transpose(hist(nonzeros(A549_EGFR.Starved(:,4)),bins)); 
f3=transpose(hist(nonzeros(A549_EGFR.StarvedGM(:,4)),bins)); 
f4=[f/sum(f) f2/sum(f2) f3/sum(f3)];

subplot(2,2,4)
bar(bins,f4,1,'grouped');

title('Hist of Density A549 EGFR');
xlabel('Density (mol/?m^2)');
ylabel('norm counts');
axis([0 1e5 0 0.25]);
leg2=legend('Unchanged','Starved','StarvedGM');
set(leg2,'FontSize',12);

%% 

% A549 EGFR Molecules per cluster

bins=0:50:500;

f=transpose(hist((A549_EGFR.Unchanged(:,4)).*A549_EGFR.Unchanged(:,2),bins)); 
f2=transpose(hist((A549_EGFR.Starved(:,4)).*A549_EGFR.Starved(:,2),bins)); 
f3=transpose(hist((A549_EGFR.StarvedGM(:,4)).*A549_EGFR.StarvedGM(:,2),bins)); 
f4=[f/sum(f) f2/sum(f2) f3/sum(f3)];

% subplot(2,2,4)
bar(bins./15,f4,1,'grouped');

title('Molecules per Cluster');
xlabel('Molecules per cluster');
ylabel('norm counts');
axis([0 50 0 0.5]);
leg2=legend('Unchanged','Starved','StarvedGM');
set(leg2,'FontSize',12);

%% Compare Area/Density between cells

% Area MDCK vs. A549 

bins=0:0.003:0.1;

f=transpose(hist(nonzeros(MDCK_EGFR.Unchanged(:,5)),bins)); 
f2=transpose(hist(nonzeros(A549_EGFR.Unchanged(:,5)),bins)); 
f4=[f/sum(f) f2/sum(f2)];

figure('Position',[200 300 800 300])
set(gcf,'numbertitle','off','name','Compare A549 MDCK ','PaperOrientation','landscape') % Title of the figure


subplot(1,2,1)
bar(bins,f4,1,'grouped');
title('Hist of Size MDCK A549');
xlabel('Diameter, \mum','FontSize',14);
ylabel('norm counts','FontSize',14);
axis([0 0.1 0 0.12]);
leg1=legend('MDCK','A549');
set(leg1,'FontSize',14);
set(gca,'FontSize',14);
box on

clear x x2 x3 f f2 f3 f4 bins;

% Density MDCK vs. A549 

bins=0:3e3:1e5;

f=transpose(hist(nonzeros(MDCK_EGFR.Unchanged(:,4)),bins)); 
f2=transpose(hist(nonzeros(A549_EGFR.Unchanged(:,4)),bins)); 
f4=[f/sum(f) f2/sum(f2)];

subplot(1,2,2)
bar(bins,f4,1,'grouped');
title('Hist of Density MDCK A549');
xlabel('Cluster density (mol/\mum^2)','FontSize',14);
ylabel('norm counts','FontSize',14);
axis([0 1e5 0 0.2])
leg2=legend('MDCK','A549');
set(leg2,'FontSize',12);
set(gca,'FontSize',14);
box on

%% 

scatter(MDCK_EGFR.Unchanged(:,2),MDCK_EGFR.Unchanged(:,4),1);hold on;
scatter(A549_EGFR.Unchanged(:,2),A549_EGFR.Unchanged(:,4),1,'red');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Load single clusters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  structure3=[0 0 0];
 
subset(:,3)=subset(:,3)+max(structure3(:,3));
structure3=cat(1,structure3, subset);        

clear index subset

%% 

%% NN search for individual clusters

NoN={};
Gradient={};

tic


 for    index=1:1:max(structure3(:,3));

        vx=find(structure3(:,3) == index);
        cluster=structure3(vx,1);
        cluster(:,2)=structure3(vx,2);
        
                if length(vx)>50; %& length(vx)< 100000 ;                              % only cluster with more than n points
                    
                idx = rangesearch(cluster,cluster,0.03);
        
                        for i=1:1:length(idx);
                        NoN{index,1}(i,1)=length(idx{i,1});     % count the total number of neighbors for each point in dataset
                        end
         
                scatter(cluster(:,1),cluster(:,2),2,NoN{index,1});hold on;
                axis([0 12 0 12]);
                colorbar
                colormap(jet)
%                 axis([4.2 5.8 4.3 5.9]);
                
                Gradient{index,3}=max(NoN{index,1})/min(NoN{index,1});  % fold increase
                Gradient{index,2}=max(NoN{index,1});                    % max NN
                Gradient{index,1}=min(NoN{index,1});                    % min NN
                
                else
            
                end
        
        clear vx idx cluster;
        
end
toc

%% Plot specific cluster

clusNbr=1759;

        vx=find(structure3(:,3) == clusNbr);
        cluster=structure3(vx,1);
        cluster(:,2)=structure3(vx,2);
        
        idx = rangesearch(cluster,cluster,0.03);
        
for i=1:1:length(idx);
NoNind(i,1)=length(idx{i,1});     % count the total number of neighbors for each point in dataset
end

figure
%  scatter(cluster(:,1),cluster(:,2),2,NoNind);hold on;
 scatter3(cluster(:,1),cluster(:,2),NoNind,3,NoNind);hold on;
% axis([min(cluster(:,1)) max(cluster(:,1)) min(cluster(:,2)) max(cluster(:,2))]);
% axis([0 12 0 10]);
colorbar
colormap(jet)
                
   
% clear NoNind

%% Plot histogram of ratio between min/max --> fold increase of density

binCenters = 0:5:200;

% x=transpose(hist(cell2mat(Gradient_unchanged(:,end)),binCenters)); 
% x2=transpose(hist(cell2mat(Gradient_starved(:,end)),binCenters)); 
% x3=transpose(hist(cell2mat(Gradient_starved_GM(:,end)),binCenters)); 

x=transpose(hist(cell2mat(Gradient_unchanged_A549(:,end)),binCenters)); 
x2=transpose(hist(cell2mat(Gradient_starved_A549(:,end)),binCenters)); 
x3=transpose(hist(cell2mat(Gradient_starved_GM_A549(:,end)),binCenters));

x4=[x/sum(x) x2/sum(x2) x3/sum(x3)];

% subplot(1,2,2)
bar(binCenters,x4);hold on;
% bar(x2,f2/sum(f2), 0.3);
axis([0 100 0 0.3]);
title('Internal density gradient EGFR A549');
xlabel('density ratio (min/max) ');
ylabel('norm counts');
leg2=legend('unchanged','starved', 'starved GM');
set(leg2,'FontSize',12);

MedianUnchanged=median(x)
MeanUnchanged=mean(x)
MedianStarved=median(x2)
MeanStarved=mean(x2)
MedianStarvedGM=median(x3)
MeanStarvedGM=mean(x3)

%% 

bins=0:400:10000;

f=transpose(hist(MDCK_EGFR.Unchanged(MDCK_EGFR.Unchanged(:,2)<0.01,2)*10^6,bins));
f2=transpose(hist(A549_EGFR.Unchanged(A549_EGFR.Unchanged(:,2)<0.01,2)*10^6,bins)); 
f4=[f/sum(f) f2/sum(f2)];

figure('Position',[200 300 800 400])
set(gcf,'numbertitle','off','name','Compare A549 MDCK ','PaperOrientation','landscape') % Title of the figure

subplot(1,2,1)
bar(bins,f4,1,'grouped');
% title('Hist of Size MDCK A549');
xlabel('Area, nm^2','FontSize',14);
ylabel('norm counts','FontSize',14);
axis([0 10000 0 0.2]);
leg1=legend('MDCK','A549');
set(leg1,'FontSize',14);
set(gca,'FontSize',14);
box on

bins=0:5:100;

f=transpose(hist(MDCK_EGFR.Unchanged(MDCK_EGFR.Unchanged(:,5)<0.5,5)*1e3,bins));
f2=transpose(hist(A549_EGFR.Unchanged(A549_EGFR.Unchanged(:,5)<0.5,5)*1e3,bins)); 
f4=[f/sum(f) f2/sum(f2)];

% 
% figure('Position',[200 300 400 400])
% set(gcf,'numbertitle','off','name','Compare A549 MDCK ','PaperOrientation','landscape') % Title of the figure

subplot(1,2,2)

bar(bins,f4,1,'grouped');
% title('Hist of Size MDCK A549');
xlabel('Radius, nm','FontSize',14);
ylabel('norm counts','FontSize',14);
axis([0 100 0 0.25]);
leg1=legend('MDCK','A549');
set(leg1,'FontSize',14);
set(gca,'FontSize',14);
box on



