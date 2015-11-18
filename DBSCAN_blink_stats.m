%% Performs DBSCAN and analysis (Delauney and HK) on full or subset input image

% eg for well clustered data, easy to find clusters --> EGF

%%%%%%%%%%INPUT%%%%%%%%%%

% filename_peaks='cell_1_GPF.prm'             % filename of PS output file
% k=10;                                       % minimum size cluster
% Eps=0.035;                                  % minimum distance between points   


%%%%%%%%%%OUTPUT%%%%%%%%%%

%   newCon with 8 columns


% 1 number of cluster
% 2 area from convex hull
% 3 area from delaunay triangle
% 4 density as mol/area
% 5 mean diameter
% 6 stdev diameter
% 6 radave from HK
% 7 aspratio from HK
% 8 clustarea of diamter

% OUTPUT NAME: DBSCAN_delaunay_HK_"filename_peaks"

%% Clear All
clear all, close all, clc, clear all
%%  %%%%%%%%%%INPUT%%%%%%%%%%

filename_peaks='locResults_FOV1_30ms_200mW_1';                 % filename of PS output file
k=7;                                         % minimum size cluster
Eps=0.020;                                   % minimum distance between points   

%% Load and Plot data

filename2=[filename_peaks '.txt'];
peaks=dlmread(filename2,',',1,0);

pix=0.1;
sdx=pix.*nonzeros(peaks(:,2));% 3,20
sdy=pix.*nonzeros(peaks(:,3));% 4,21
all(:,1)=sdx;
all(:,2)=sdy;

all=unique(all,'rows');

% figure
% scatter(sdx,sdy,1);

fprintf('\n -- Data loaded --\n')

%% %% Cluster DBSCAN

% input image --> all
tic
[class,type]=DBSCAN(all,k,Eps);     % uses parameters specified at input
class2=transpose(class);            % class - vector specifying assignment of the i-th object to certain cluster (m,1)
type2=transpose(type);              % (core: 1, border: 0, outlier: -1)
% toc
% 
% fprintf('\n -- DBSCAN computed --\n')

fprintf(' -- DBSCAN computed in %f sec -- \n',toc)

%% %% Find core points --> type = 1 or 0 and define them as subset

coreBorder=find(type2 >= 0);

subset(:,1)=all(coreBorder,1);
subset(:,2)=all(coreBorder,2);
subset(:,3)=class2(coreBorder);

figure('Position',[700 600 600 250])
subplot(1,2,1)
scatter(sdx,sdy,1);
title('Raw Data')
axis on
% axis([0 13 0 13])

subplot(1,2,2)
scatter(subset(:,1),subset(:,2),2,mod(subset(:,3),10))
title('Found Clusters')
axis on
% axis([0 13 0 13])                        

inCluster = length(subset)/length(all) % percent of molecules in cluster

%% Save clusters from 1st and 2nd DBSCAN identification

filenamec1=['Step_DBSCAN_clusters_1st_' filename_peaks];

save(filenamec1,'subset');

fprintf('\n -- Clusters saved --\n')

%% Calculate Area (delaunay and convhull), density and diameter for all clusters 


numbr=subset(:,3);
maxk=max(numbr);
index=1;

for numbr=1:maxk;
    
    vi=find(subset(:,3)==numbr);
    
    x=subset(vi,1);
    y=subset(vi,2);
    
    clusterxCenter=sum(x)/length(x);
    clusteryCenter=sum(y)/length(y);
    
    [conh,Aconh]=convhull(x,y);
    Tri=delaunay(x,y);
    areaTri =polyarea(x(Tri'),y(Tri'));
    
    
    outerPoints(:,1)=x(conh);
    outerPoints(:,2)=y(conh);

    Dx=pdist(outerPoints); 
    Dx=transpose(Dx);
    
    newCon(index,1)=numbr;                  % number of cluster
    newCon(index,2)=Aconh;                  % area from convex hull
    newCon(index,3)=sum(areaTri);           % area from delaunay triangle
    newCon(index,4)=length(vi)/Aconh;       % density as mol/area
    newCon(index,5)=mean(Dx);               % mean diameter
    newCon(index,6)=std(Dx);                % stdev of diamter
    newCon(index,7)=length(vi);             % number of molecules
    
    Centers(index,1)= clusterxCenter;
    Centers(index,2)= clusteryCenter;
    
   
    
    clear outerPoints;
    
    index=index+1;

end 

[KCell,AreaCell]=convhull(Centers(:,1),Centers(:,2));

%% For Blinking analysis

binCenters = 0:5:100;
[N,X] = hist(newCon(:,7),binCenters);

figure

N = N./length(newCon(:,7));
bar(X(2:end-1),N(2:end-1));
title('Number of Molecules')


MeanMol=mean(newCon(:,7))


%% Plot results from newCon

figure('Position',[900 50 900 600])

subplot(2,3,1)
hist(newCon(:,2),30)
title('Area from Conv Hull')

subplot(2,3,2)
hist(newCon(:,3),30)
title('Area from Delaunay')

subplot(2,3,3)
hist(newCon(:,4),30)
title('density (local/area)')

subplot(2,3,4)
hist(newCon(:,5),30)
title('Mean Diameter')

subplot(2,3,5)
hist(newCon(:,6),30)
title('Stdev Diameter')

subplot(2,3,6)
plot(Centers(KCell,1),Centers(KCell,2),'r-',Centers(:,1),Centers(:,2),'b+')
title('Cluster Centers and ConvHull')

%% Calculate Mean distance between clusters
   % Calculate distance distribution of cluster centers and fit to normal distribution

tic

D=pdist(Centers,'euclidean');

[f,xi] = ksdensity(D);

figure('Position',[900 50 900 200])
subplot(1,3,1)
hist(D,50);hold on;
title('pdist of centers in pxl')
subplot(1,3,2)
plot(xi,f);
title('pdf of dist centers in pxl')
subplot(1,3,3)
scatter(Centers(:,1),Centers(:,2),3)
title('Cluster Centers')

toc

MeanClusDis=fitdist(D(:),'normal')

ClusterDensity=length(Centers)/AreaCell  % Clusters per µm^2



%% Save Output for HK clustering file

all_out(:,1)=all(coreBorder);
all_out(:,2)=all(coreBorder,2);
all_out(:,3)=class2(coreBorder);

%subset(:,3)=class2;
%save('clusterDBSCAN_cell1.txt', 'all_out','-ascii','-tabs' );
%setColTitles(all_out, {'Xpos', 'Ypos','cluster number'});

filename = ['For_HK_OUT_' filename_peaks '.txt'] ;
dlmwrite(filename, all_out)


%% Post processing with Hoshen-Kopelman

%filename_clusters='clusterDBSCAN_cell3.txt'

  maxr_weight_loc(j)=0;
  meanr_weight(j)=0;
% take out line 104

clusterfile = filename; 
Cluster_analysis_simulation  % in K:\Christian\matlab\new clustering

%% Export newCon

% 1 number of cluster
% 2 area from convex hull
% 3 area from delaunay triangle
% 4 density as mol/area
% 5 mean diameter
% 6 stdev diameter
% 6 radave from HK
% 7 aspratio from HK
% 8 clustarea from HK

 newCon(:,7)=radave;                    % radius from HK
 newCon(:,8)=aspratio;                  % aspratio from HK
 newCon(:,9)=clustarea;                 % clustarea from HK

filename = ['DBSCAN_delaunay_HK_' filename_peaks '.txt'];
dlmwrite(filename, newCon)



