clear, close all, clc

%% Load Data

% For individual cluster analysis

structure3=[];

a=load('DBSCAN_clusters_A549_EGF_A647_2000mW_10ms__2_MMStack_locResults_DC.mat');
structure3=cat(1,structure3, a.subset);

b=load('DBSCAN_clusters_A549_EGF_A647_2000mW_10ms__3_MMStack_locResults_DC.mat');
b.subset(:,3)=b.subset(:,3)+max(a.subset(:,3));
structure3=cat(1,structure3, b.subset);

c=load('DBSCAN_clusters_A549_EGF_A647_2000mW_10ms__4_MMStack_locResults_DC.mat'); 
c.subset(:,3)=c.subset(:,3)+max(b.subset(:,3));
structure3=cat(1,structure3, c.subset);

d=load('DBSCAN_clusters_A549_EGF_A647_2000mW_10ms__1_MMStack_locResults_DC.mat');
d.subset(:,3)=d.subset(:,3)+max(c.subset(:,3));
structure3=cat(1,structure3, d.subset);

e=load('DBSCAN_clusters_A549_EGF_A647_2000mW_10ms__6_MMStack_locResults_DC.mat');
e.subset(:,3)=e.subset(:,3)+max(d.subset(:,3));
structure3=cat(1,structure3, e.subset);

f=load('DBSCAN_clusters_A549_EGF_A647_2000mW_10ms__9_MMStack_locResults_DC.mat');
f.subset(:,3)=f.subset(:,3)+max(e.subset(:,3));
structure3=cat(1,structure3, f.subset);

EGFR_ind=[a;b;c;d;e;f];
       
clear index subset

clear a b c d e f 

% For population cluster analysis

% 1 number of cluster
% 2 area from convex hull
% 3 area from delaunay triangle
% 4 density as mol/area
% 5 mean diameter
% 6 stdev diameter
% 7 Mean Mol

a=dlmread('DBSCAN_delaunay_HK_A549_EGF_A647_2000mW_10ms__2_MMStack_locResults_DC.txt');
b=dlmread('DBSCAN_delaunay_HK_A549_EGF_A647_2000mW_10ms__3_MMStack_locResults_DC.txt');
c=dlmread('DBSCAN_delaunay_HK_A549_EGF_A647_2000mW_10ms__4_MMStack_locResults_DC.txt');
d=dlmread('DBSCAN_delaunay_HK_A549_EGF_A647_2000mW_10ms__1_MMStack_locResults_DC.txt');
e=dlmread('DBSCAN_delaunay_HK_A549_EGF_A647_2000mW_10ms__6_MMStack_locResults_DC.txt');
% f=dlmread('DBSCAN_delaunay_HK_A549_EGF_A647_2000mW_10ms__10_MMStack_locResults_DC.txt');

EGFR=[a;b;c;d;e];

clear a b c d e 

%% Population Statistics

figure('Position',[500 400 1000 600])
h=gcf;
set(h,'PaperOrientation','landscape');

binCenters = 0:5:200;
x=transpose(hist(EGFR(:,5),binCenters)); 

subplot(2,2,1)
bar(binCenters,x/sum(x));hold on;
axis([0 100 0 0.15]);
title('Cluster Diameter');
xlabel('cluster diameter [nm] ');
ylabel('norm counts');

binCenters = 0:0.01:0.2;
x=transpose(hist(EGFR(:,4),binCenters)); 

subplot(2,2,2)
bar(binCenters,x/sum(x));hold on;
axis([0 0.2 0 0.3]);
title('Molecule Density');
xlabel('molecule density [mol/nm^2] ');
ylabel('norm counts');


binCenters = 0:10:400;

x=transpose(hist(EGFR(:,7),binCenters)); 

subplot(2,2,3)
bar(binCenters,x/sum(x));hold on;
axis([0 400 0 0.15]);
title('Number of Localizations');
xlabel('# of Localizations');
ylabel('norm counts');


binCenters = 0:500:8000;
x=transpose(hist(EGFR(:,2),binCenters)); 

subplot(2,2,4)
bar(binCenters,x/sum(x));hold on;
axis([0 8000 0 0.3]);
title('Cluster Area');
xlabel('cluster area [nm^2] ');
ylabel('norm counts');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Single clusters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NN search for individual clusters

NoN={};
Gradient={};

tic


 for    index=1:100 % max(structure3(:,3));

        vx=find(structure3(:,3) == index);
        cluster=structure3(vx,1);
        cluster(:,2)=structure3(vx,2);
        
                if length(vx)>50; %& length(vx)< 100000 ;                              % only cluster with more than n points
                    
                idx = rangesearch(cluster,cluster,30);
        
                        for i=1:1:length(idx);
                        NoN{index,1}(i,1)=length(idx{i,1});     % count the total number of neighbors for each point in dataset
                        end
         
                scatter(cluster(:,1),cluster(:,2),2,NoN{index,1});hold on;
%                 axis([0 12 0 12]);
                colorbar
                colormap(jet)
                
                Gradient{index,3}=max(NoN{index,1})/min(NoN{index,1});  % fold increase
                Gradient{index,2}=max(NoN{index,1});                    % max NN
                Gradient{index,1}=min(NoN{index,1});                    % min NN
                
                else
            
                end
        
        clear vx idx cluster;
        
end
toc

%% Plot density histogram of one specific cluster

clusNbr=10;
cluster=[];
NoNind=[];

        vx=find(structure3(:,3) == clusNbr);
        cluster=structure3(vx,1);
        cluster(:,2)=structure3(vx,2);
        
        idx = rangesearch(cluster,cluster,30);
        
for i=1:length(idx);
NoNind(i,1)=length(idx{i,1});     % count the total number of neighbors for each point in dataset
end

figure
scatter3(cluster(:,1),cluster(:,2),NoNind,3,NoNind);hold on;
colorbar
colormap(jet)


%% Plot histogram of ratio between min/max --> fold increase of density

binCenters = 0:3:50;
x=transpose(hist(cell2mat(Gradient(:,end)),binCenters)); 

figure('Position',[200 300 800 400])
set(gcf,'numbertitle','off','name','Internal Density Gradient','PaperOrientation','landscape') % Title of the figure

bar(binCenters,x/sum(x));hold on;
axis([-10 50 0 0.6]);
title('Internal density gradient EGFR');
xlabel('density ratio (min/max) ');
ylabel('norm counts');

Median=median(x)
Mean=mean(x)

%% Plot clusters with color according to frame

figure

for    index=1:20 % max(structure3(:,3));

        vx=find(structure3(:,3) == index);
        cluster=structure3(vx,1);
        cluster(:,2)=structure3(vx,2);
        cluster(:,3)=structure3(vx,4);
        
        scatter(cluster(:,1),cluster(:,2),1,(cluster(:,3))); hold on;
        title('identified Clusters (color -> frame)')
        colorbar
        axis on
%         axis([0 max(cluster(:,1)) 0 max(cluster(:,2))])  
        
        clear cluster
                    
end
