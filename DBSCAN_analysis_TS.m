%%  DBSCAN and basic analysis from HTP data (ThunderSTORM output)

%%%%%%%%%%INPUT%%%%%%%%%%

% filename_peaks='cell_1_GPF'                 % filename of TS output file
% k=10;                                       % minimum size cluster
% Eps=0.035;                                  % minimum distance between points   
% thresh=200;                                 % loglikelihood threshold

%%%%%%%%%%OUTPUT%%%%%%%%%%

%   newCon with 8 columns

% 1 number of cluster
% 2 area from convex hull
% 3 area from delaunay triangle
% 4 density as mol/area
% 5 mean diameter
% 6 stdev diameter
% 7 Mean Mol


% OUTPUT NAME: DBSCAN_delaunay_HK_"filename_peaks"

% clear, close all, clc, clear

function [newCon,inCluster,MeanMol,ClusterDensity,MeanClusDis]=DBSCAN_analysis_TS(filename_peaks,figures) 
%%  %%%%%%%%%%INPUT%%%%%%%%%%

% filename_peaks%='locResults';                        % filename of TS output file
filename_peaks2=[filename_peaks '.txt'];
peaks=dlmread(filename_peaks2,',',1,0);

fprintf('\n -- Data loaded --\n')
%% Select Parameters
                                       % loglikelihood threshold
k=3;                                               % minimum number of neighbors within Eps
Eps=0.020;                                             % minimum distance between points, nm

fprintf('\n -- Parameters selected --\n')
%% Filter for loglikelihood and Plot data

% Plot  data

dataDBS(:,1)=peaks(:,2)/1000; % x in mum
dataDBS(:,2)=peaks(:,3)/1000; % y in mum

scatter(dataDBS(:,1),dataDBS(:,2),1);

fprintf('\n -- Data plotted --\n')


%% %% Cluster DBSCAN

% input image --> all
tic
[class,type]=DBSCAN(dataDBS,k,Eps);     % uses parameters specified at input
class2=transpose(class);            % class - vector specifying assignment of the i-th object to certain cluster (m,1)
type2=transpose(type);              % (core: 1, border: 0, outlier: -1)

fprintf(' -- DBSCAN computed in %f sec -- \n',toc)

%% Find core points --> type = 1 or 0 and define them as subset

coreBorder=find(type2 >= 0);

subset(:,1)=dataDBS(coreBorder,1);
subset(:,2)=dataDBS(coreBorder,2);
subset(:,3)=class2(coreBorder);

if figures==1;

figure('Position',[700 600 900 400])
subplot(1,2,1)
scatter(dataDBS(:,1),dataDBS(:,2),1);
title('Raw Data')
axis on
axis([0 max(dataDBS(:,1)) 0 max(dataDBS(:,2))])

subplot(1,2,2)
scatter(subset(:,1),subset(:,2),1,mod(subset(:,3),10))
title('identified Clusters')
axis on
axis([0 max(dataDBS(:,1)) 0 max(dataDBS(:,2))])  

else 
end

inCluster = length(subset)/length(dataDBS) % percent of molecules in cluster

%% Save clusters from 1st and 2nd DBSCAN identification

% filenamec1=['Step_DBSCAN_clusters_1st_' filename_peaks];
% 
% save(filenamec1,'subset');
% 
% fprintf('\n -- Clusters saved --\n')

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

if figures==1;

figure

N = N./length(newCon(:,7));
bar(X(2:end-1),N(2:end-1));
title('Number of Molecules per cluster')
else
end

MeanMol=mean(newCon(:,7))


%% Plot results from newCon

if figures==1;

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

else
end

%% Calculate Mean distance between clusters
   % Calculate distance distribution of cluster centers and fit to normal distribution

tic

D=pdist(Centers,'euclidean');

[f,xi] = ksdensity(D);

if figures==1;

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

else 
end    

toc

MeanClusDis=fitdist(D(:),'normal')

ClusterDensity=length(Centers)/AreaCell  % Clusters per µm^2


%% Export newCon

% 1 number of cluster
% 2 area from convex hull
% 3 area from delaunay triangle
% 4 density as mol/area
% 5 mean diameter
% 6 stdev diameter

filename = ['DBSCAN_delaunay_HK_' filename_peaks '.txt'];
dlmwrite(filename, newCon);

end