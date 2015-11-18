%%  DBSCAN and basic analysis from HTP data (ThunderSTORM output)

%%%%%%%%%%INPUT%%%%%%%%%%

% filename_peaks='cell_1_GPF'                 % filename of TS output file
% k1+2      =10;                              % minimum size cluster
% Eps1+2    =0.035;                           % minimum distance between points   
% thresh    =500;                             % loglikelihood threshold

%%%%%%%%%%OUTPUT%%%%%%%%%%

%   newCon with 8 columns

% 1 number of cluster
% 2 area from convex hull
% 3 area from delaunay triangle
% 4 density as mol/area
% 5 mean diameter
% 6 stdev diameter
% 7 Mean Mol
%
% OUTPUT NAME: 'Step_DBSCAN_' filename_peaks  --> Structure that contains 3
%
% Structures            1. newCon_1st
%                       2. newCon_2nd
%                       3. both

function [structure, structure2, structure3,inCluster,MeanMol,ClusterDensity,MeanClusDis]=Stepwise_DBSCAN_TS(fig,filename) 

%% Select input parameters

%%%%%%%%%%% Select DBSCAN Parameters %%%%%%%%%%%%%%%%%%%%%%
%%%%% DBSCAN parameter example 15/0.015 and 15/0.020%%%%%%%

k1=3;
Eps1=0.020;
k2=3;
Eps2=0.020;

fprintf('\n -- Parameters selected --\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load file, apply LL threshold, plot data

filename_peaks=filename;
filename2=[filename_peaks '.txt'];
peaks=dlmread(filename2,',',1,0);

fprintf('\n -- Data loaded --\n')

% Plot  data

dataDBS(:,1)=peaks(:,2)/1000; % x in mum
dataDBS(:,2)=peaks(:,3)/1000; % y in mum

scatter(dataDBS(:,1),dataDBS(:,2),1);

fprintf('\n -- Data plotted --\n')

%% Select region of the image 
%%%%%%%%%%% select region of interest %%%%%%%%%%%

upperx=max(dataDBS(:,1));
lowerx=0;

uppery=max(dataDBS(:,2));
lowery=0;

%%%%%%%%%%% crop region %%%%%%%%%%%

vx=find(dataDBS(:,1) < upperx & dataDBS(:,1) > lowerx);
subset=dataDBS(vx);
subset(:,2)=dataDBS(vx,2);

vy=find(subset(:,2) < uppery & subset(:,2) > lowery);
subset2=subset(vy);
subset2(:,2)=subset(vy,2);

if fig==1;

figure
scatter(subset2(:,1), subset2(:,2),1,'black'), hold on

else end

fprintf('\n -- ROI selected --\n')

%% Cluster DBSCAN

% input image --> subset2
%%%%%%%%%%%%%%%%%%%%%%% 1st DBSCAN %%%%%%%%%%%%%%%%%%%%%%
tic
[class,type]=DBSCAN(subset2,k1,Eps1); % 
class2_1st=transpose(class);
type2=transpose(type);

target_1st=find(type2 >= 0);
result_1st=subset2(target_1st,:);

%%%%%%%%%%%%%%%%%%%%%%% 2nd DBSCAN %%%%%%%%%%%%%%%%%%%%%%

target_2nd=find(type2 < 0);
subset_2nd=subset2(target_2nd,:);

[class,type]=DBSCAN(subset_2nd,k2,Eps2);
class2_2nd=transpose(class);
type2=transpose(type);

target_2nd=find(type2 >= 0);
result_2nd=subset_2nd(target_2nd,:);

fprintf(' -- DBSCAN computed in %f min -- \n',toc/60)

%%%%%%%%%%%%%%%%%%%%%%%% Plot found cluster %%%%%%%%%%%%%%%%%%%%%%%

if fig==1;

 figure
 scatter(result_1st(:, 1), result_1st(:, 2),1,mod(class2_1st(target_1st),10));hold on
 scatter(result_2nd(:, 1), result_2nd(:, 2),1,mod(class2_2nd(target_2nd),10)); hold on
 legend('1st round','2nd round');
 
 
 figure
 scatter(result_1st(:, 1), result_1st(:, 2),1,'red');hold on
 scatter(result_2nd(:, 1), result_2nd(:, 2),1,'blue'); hold on
 legend('1st round','2nd round');
 
else end


%% Save clusters from 1st and 2nd DBSCAN identification
% For single cluster analysis
%     density, gradient etc.

filenamec1=['Step_DBSCAN_clusters_1st_DBSCAN_' filename_peaks];
filenamec2=['Step_DBSCAN_clusters_2nd_DBSCAN' filename_peaks];


save(filenamec1,'result_1st');
save(filenamec2,'result_2nd');

fprintf('\n -- Clusters Saved --\n')

%% Calculate Histogram and save as image
 
% c=hist3([result_1st(:, 1), result_1st(:, 2)],[600 600]); % heigth x width
% d=hist3([result_2nd(:, 1), result_2nd(:, 2)],[600 600]);
% 
% imwrite(c,'DBSCAN_1.tiff');
% imwrite(d,'DBSCAN_2.tiff');
 
%% Assign cluster numbers for 1st and 2nd results

result_1st(:,3)=class2_1st(target_1st);  % x,y of the cluster and cluster number from class
result_2nd(:,3)=class2_2nd(target_2nd);

inCluster=(length(result_1st)+length(result_2nd))/length(subset2)

%%  Analyse First Results and put them into newCon_1st

numbr=result_1st(:,3);
maxk=max(numbr);
index=1;

% Plot the cropped ROI and overlay with found clusters

if fig==1;
    
figure
scatter(subset2(:,1), subset2(:,2),1,'black'), hold on;

for numbr=1:maxk;
    
    vi=find(result_1st(:,3)==numbr);
    
    x=result_1st(vi,1);
    y=result_1st(vi,2);
    
    clusterxCenter=sum(x)/length(x);
    clusteryCenter=sum(y)/length(y);
    
    [conh,Aconh]=convhull(x,y);
    Tri=delaunay(x,y);
    areaTri =polyarea(x(Tri'),y(Tri'));
    
    plot(x(conh),y(conh),'red'); hold on % Overlay plot over cropped ROI  
    
    outerPoints(:,1)=x(conh);
    outerPoints(:,2)=y(conh);

    Dx=pdist(outerPoints); 
    Dx=transpose(Dx);
    
    newCon_1st(index,1)=numbr;                  % number of cluster
    newCon_1st(index,2)=Aconh;                  % area from convex hull
    newCon_1st(index,3)=sum(areaTri);           % area from delaunay triangle
    newCon_1st(index,4)=length(vi)/Aconh;       % density as mol/area
    newCon_1st(index,5)=mean(Dx);               % mean diameter
    newCon_1st(index,6)=std(Dx);                % stdev of diamter
    newCon_1st(index,7)=length(vi);             % number of molecules

    Centers_1st(index,1)= clusterxCenter;
    Centers_1st(index,2)= clusteryCenter;
    
    clear outerPoints;
    
    index=index+1;

end 

else
  

for numbr=1:maxk;
    
    vi=find(result_1st(:,3)==numbr);
    
    x=result_1st(vi,1);
    y=result_1st(vi,2);
    
    clusterxCenter=sum(x)/length(x);
    clusteryCenter=sum(y)/length(y);
    
    [conh,Aconh]=convhull(x,y);
    Tri=delaunay(x,y);
    areaTri =polyarea(x(Tri'),y(Tri'));
    
    
    outerPoints(:,1)=x(conh);
    outerPoints(:,2)=y(conh);

    Dx=pdist(outerPoints); 
    Dx=transpose(Dx);
    
    newCon_1st(index,1)=numbr;                  % number of cluster
    newCon_1st(index,2)=Aconh;                  % area from convex hull
    newCon_1st(index,3)=sum(areaTri);           % area from delaunay triangle
    newCon_1st(index,4)=length(vi)/Aconh;       % density as mol/area
    newCon_1st(index,5)=mean(Dx);               % mean diameter
    newCon_1st(index,6)=std(Dx);                % stdev of diamter
    newCon_1st(index,7)=length(vi);             % number of molecules

    Centers_1st(index,1)= clusterxCenter;
    Centers_1st(index,2)= clusteryCenter;
       
    
    clear outerPoints;
    
    index=index+1;

end     

end



%% Analyse second Results and put them into newCon_2nd

numbr=result_2nd(:,3);
maxk=max(numbr);
index=1;

if fig==1;
for numbr=1:maxk;
    
    vi=find(result_2nd(:,3)==numbr);
    
    x=result_2nd(vi,1);
    y=result_2nd(vi,2);
    
    clusterxCenter=sum(x)/length(x);
    clusteryCenter=sum(y)/length(y);
    
    [conh,Aconh]=convhull(x,y);
    Tri=delaunay(x,y);
    areaTri =polyarea(x(Tri'),y(Tri'));
    
    plot(x(conh),y(conh),'blue'); hold on % Overlay plot over cropped ROI    
    
    outerPoints(:,1)=x(conh);
    outerPoints(:,2)=y(conh);

    Dx=pdist(outerPoints); 
    Dx=transpose(Dx);
    
    newCon_2nd(index,1)=numbr;                  % number of cluster
    newCon_2nd(index,2)=Aconh;                  % area from convex hull
    newCon_2nd(index,3)=sum(areaTri);           % area from delaunay triangle
    newCon_2nd(index,4)=length(vi)/Aconh;       % density as mol/area
    newCon_2nd(index,5)=mean(Dx);               % mean diameter
    newCon_2nd(index,6)=std(Dx);                % stdev of diamter
    newCon_2nd(index,7)=length(vi);             % number of molecules
    clear outerPoints;
    
    Centers_2nd(index,1)= clusterxCenter;
    Centers_2nd(index,2)= clusteryCenter;
    
    index=index+1;

end

else

for numbr=1:maxk;
    
    vi=find(result_2nd(:,3)==numbr);
    
    x=result_2nd(vi,1);
    y=result_2nd(vi,2);
    
    clusterxCenter=sum(x)/length(x);
    clusteryCenter=sum(y)/length(y);
    
    [conh,Aconh]=convhull(x,y);
    Tri=delaunay(x,y);
    areaTri =polyarea(x(Tri'),y(Tri'));
     
    
    outerPoints(:,1)=x(conh);
    outerPoints(:,2)=y(conh);

    Dx=pdist(outerPoints); 
    Dx=transpose(Dx);
    
    newCon_2nd(index,1)=numbr;                  % number of cluster
    newCon_2nd(index,2)=Aconh;                  % area from convex hull
    newCon_2nd(index,3)=sum(areaTri);           % area from delaunay triangle
    newCon_2nd(index,4)=length(vi)/Aconh;       % density as mol/area
    newCon_2nd(index,5)=mean(Dx);               % mean diameter
    newCon_2nd(index,6)=std(Dx);                % stdev of diamter
    newCon_2nd(index,7)=length(vi);             % number of molecules
    clear outerPoints;
    
    Centers_2nd(index,1)= clusterxCenter;
    Centers_2nd(index,2)= clusteryCenter;
    
    index=index+1;

end

end


%% Mean Number of Molecules, Cluster density, mean distance

MeanMol=mean(cat(1,newCon_1st(:,7),newCon_2nd(:,7)))

Centers=[Centers_1st; Centers_2nd];

[KCell,AreaCell]=convhull(Centers(:,1),Centers(:,2));

ClusterDensity=length(Centers)/AreaCell

D=pdist(Centers,'euclidean');
[f,xi] = ksdensity(D);

if fig==1;

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


MeanClusDis=fitdist(D(:),'normal')

%% Plot results from newCon_1st and newCon_2nd

if fig==1;

figure('Position',[900 50 900 600])

subplot(2,3,1)
hist(newCon_2nd(:,2),30);hold on;
hist(newCon_1st(:,2),30);hold on;
h = findobj(gca,'Type','patch');
set(h(1),'Facecolor',[1 0 0],'EdgeColor','k');
set(h(2),'Facecolor',[0 0 1],'EdgeColor','k'); 
title('Area from Conv Hull')

subplot(2,3,2)
hist(newCon_2nd(:,3),30);hold on;
hist(newCon_1st(:,3),30);
h = findobj(gca,'Type','patch');
set(h(1),'Facecolor',[1 0 0],'EdgeColor','k');
set(h(2),'Facecolor',[0 0 1],'EdgeColor','k'); 
title('Area from Delaunay')

subplot(2,3,3)
hist(newCon_1st(:,4),30);hold on;
hist(newCon_2nd(:,4),30);
h = findobj(gca,'Type','patch');
set(h(1),'Facecolor',[1 0 0],'EdgeColor','k');
set(h(2),'Facecolor',[0 0 1],'EdgeColor','k'); 
title('density (local/area)')

subplot(2,3,4)
hist(newCon_2nd(:,5),30);hold on;
hist(newCon_1st(:,5),30);
h = findobj(gca,'Type','patch');
set(h(1),'Facecolor',[1 0 0],'EdgeColor','k');
set(h(2),'Facecolor',[0 0 1],'EdgeColor','k'); 
title('Mean Diameter')

subplot(2,3,5)
hist(newCon_2nd(:,6),30);hold on;
hist(newCon_1st(:,6),30);
h = findobj(gca,'Type','patch');
set(h(1),'Facecolor',[1 0 0],'EdgeColor','k');
set(h(2),'Facecolor',[0 0 1],'EdgeColor','k'); 
title('Stdev Diameter')

else end

%% Generate structure from result

%%%    generates mycell using newCon_1st and newCon_2nd 
mycell={};
mycell{1,1}=newCon_1st;
mycell{1,2}=newCon_2nd;
index=1;

%%%    generates 3 structures from mycell 
%
% 1. newCon_1st
% 2. newCon_2nd
% 3. both
%

structure=struct('ConvHull', [], 'Delaunay', [], 'Density', [], 'MeanDiam', [], 'StDevDiam', [],'MeanMol', []);  % newCon_1st
structure2=struct('ConvHull', [], 'Delaunay', [], 'Density', [], 'MeanDiam', [], 'StDevDiam', [],'MeanMol', []); % newCon_2nd
structure3=struct('ConvHull', [], 'Delaunay', [], 'Density', [], 'MeanDiam', [], 'StDevDiam', [],'MeanMol', []); % both

for     index=1:length(mycell(:,1));
       
        structure.ConvHull=cat(1,structure.ConvHull, mycell{index,1}(:,2));
        structure.Delaunay=cat(1,structure.Delaunay, mycell{index,1}(:,3));
        structure.Density=cat(1,structure.Density, mycell{index,1}(:,4));
        structure.MeanDiam=cat(1,structure.MeanDiam, mycell{index,1}(:,5));
        structure.StDevDiam=cat(1,structure.StDevDiam, mycell{index,1}(:,5));
        structure.MeanMol=cat(1,structure.MeanMol, mycell{index,1}(:,6));
        
        structure2.ConvHull=cat(1,structure2.ConvHull, mycell{index,2}(:,2));
        structure2.Delaunay=cat(1,structure2.Delaunay, mycell{index,2}(:,3));
        structure2.Density=cat(1,structure2.Density, mycell{index,2}(:,4));
        structure2.MeanDiam=cat(1,structure2.MeanDiam, mycell{index,2}(:,5));
        structure2.StDevDiam=cat(1,structure2.StDevDiam, mycell{index,2}(:,6));
        structure2.MeanMol=cat(1,structure2.MeanMol, mycell{index,2}(:,6));
        
        structure3.ConvHull=cat(1,structure3.ConvHull, mycell{index,1}(:,2));
        structure3.Delaunay=cat(1,structure3.Delaunay, mycell{index,1}(:,3));
        structure3.Density=cat(1,structure3.Density, mycell{index,1}(:,4));
        structure3.MeanDiam=cat(1,structure3.MeanDiam, mycell{index,1}(:,5));
        structure3.StDevDiam=cat(1,structure3.StDevDiam, mycell{index,1}(:,6));
        structure3.MeanMol=cat(1,structure3.MeanMol, mycell{index,1}(:,6));
        
        structure3.ConvHull=cat(1,structure3.ConvHull, mycell{index,2}(:,2));
        structure3.Delaunay=cat(1,structure3.Delaunay, mycell{index,2}(:,3));
        structure3.Density=cat(1,structure3.Density, mycell{index,2}(:,4));
        structure3.MeanDiam=cat(1,structure3.MeanDiam, mycell{index,2}(:,5));
        structure3.StDevDiam=cat(1,structure3.StDevDiam, mycell{index,2}(:,6));
        structure3.MeanMol=cat(1,structure3.MeanMol, mycell{index,2}(:,6));
       
end


% 1 number of cluster
% 2 area from convex hull
% 3 area from delaunay triangle
% 4 density as mol/area
% 5 mean diameter
% 6 stdev diameter
% 7 Mean number of molecules


%% Save Structures Analysis from ConvHull and delaunay

filename=['Step_DBSCAN_' filename_peaks]

save(filename,'structure')
save(filename,'structure2','-append')
save(filename,'structure3','-append')

end
