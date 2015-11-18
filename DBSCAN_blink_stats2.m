%% Performs DBSCAN cluster and analyzes the blinking statistics of individual clusters  

%%%%%%%%%%INPUT%%%%%%%%%%

% filename_peaks='cell_1_GPF.prm'             % filename of PS output file
% k=10;                                       % minimum number of localizations in cluster
% Eps=0.035;                                  % minimum distance between points   


%%%%%%%%%%OUTPUT%%%%%%%%%%

%           Mean and Median of:
% 
%                   Number of Gaps
%                   Number of Blinks
%                   Length of Blink 
% 
%                   Gaussian fit of scatter along X and Y dimension 
%                   Localization precision (sigma from gaussian)
% 
%                   Histograms of Gaps and Blinks
%                   XY Scatter and histogram along each axis

clear, close all, clc, clear all
%% 
%%%%%%%%%%INPUT%%%%%%%%%%

filename_peaks='locResults_FOV1_30ms_200mW_2_new';% filename of PS output file
filename_peaks2=[filename_peaks '.txt'];

k=5;                                       % minimum size cluster
Eps=20;                                     % minimum distance between points   

%% Load and Plot data

peaks=dlmread(filename_peaks2,',',1,0);

pix=1; % from Thunderstorm --> data will be in nm

sdx=pix.*(peaks(:,1));% 3,20
sdy=pix.*(peaks(:,2));% 4,21
frame=peaks(:,4);% 4,21
all(:,1)=sdx-min(sdx);
all(:,2)=sdy-min(sdy);
all(:,3)=frame;

all=unique(all,'rows');

all2(:,1)=sdx;
all2(:,2)=sdy;
all2(:,3)=frame;

figure
scatter(sdx,sdy,1);

%% %% Cluster DBSCAN

% input image --> all2

tic

[class,type]=DBSCAN(all2,k,Eps);     % uses parameters specified at input
class2=transpose(class);            %class - vector specifying assignment of the i-th object to certain cluster (m,1)
type2=transpose(type);              %(core: 1, border: 0, outlier: -1)

toc

%% %% Find core points --> type = 1 or 0 and define them as subset

subset=[];
subset2=[];
n=1;
for m=1:length(type2);
    
    if type2(m,1)==0 | type2(m,1)==1
        
        subset(n,1)=all2(m,1); 
        subset(n,2)=all2(m,2);
        subset(n,3)=class2(m);
        subset(n,4)=all2(m,3);      
    else
        subset(n,1)=NaN;
        subset(n,2)=NaN; 
        subset(n,3)=NaN; 
        subset(n,4)=NaN; 
    end
    n=n+1;
end

subset2(:,1)=subset(~isnan(subset(:,1)),1);
subset2(:,2)=subset(~isnan(subset(:,2)),2);
subset2(:,3)=subset(~isnan(subset(:,3)),3);
subset2(:,4)=subset(~isnan(subset(:,4)),4);

figure('Position',[700 600 1200 500])
subplot(1,2,1)
scatter(sdx,sdy,1);
axis([0 max(sdx) 0 max(sdy)])
axis on


subplot(1,2,2)
scatter(subset2(:,1),subset2(:,2),2,subset2(:,3))
axis([0 max(sdx) 0 max(sdy)])
axis on
                        
inCluster = length(subset2)/length(all)

%% Generate track variable

        % generate tracks variable from subset
        % generate inTracks structure
        % isolate and normalize blinkframes

tracks(:,1)=subset2(:,3); % ID
tracks(:,2)=subset2(:,4); % frame
tracks(:,3)=subset2(:,1); % X
tracks(:,4)=subset2(:,2); % Y


indTracks=struct('tracks',[],'blinkframes',[],'gaplength',[],'nbrofblinks',[],'blinklength',[]);

for number=1:max(tracks(:,1));

    target=find(tracks(:,1)==number);
    
    indTracks.tracks{number,1}=tracks(target);          % Track ID, i.e. number of trajectory
    indTracks.tracks{number,1}(:,2)=tracks(target,2);   % Frame
    indTracks.tracks{number,1}(:,3)=tracks(target,3);   % x position
    indTracks.tracks{number,1}(:,4)=tracks(target,4);   % y position
    
    indTracks.blinkframes{number,1}=sortrows(tracks(target,2))-min(tracks(target,2));
        
end

%% Find gaps and calculate length

gaps=[];

for track=1:length(indTracks.blinkframes);
    
    for index=1:length(indTracks.blinkframes{track,1})-1
    
    indTracks.gaplength{track,1}(index,1)=indTracks.blinkframes{track,1}(index+1,1)-indTracks.blinkframes{track,1}(index,1)-1;
    
    end
    
    gaps=vertcat(gaps,nonzeros(indTracks.gaplength{track,1}));
    
end

%% Number of blinks and blink length

blinklengthHist=[];

for index1=1:length(indTracks.blinkframes);

count=2;   
    
blinks=indTracks.blinkframes{index1,1}+1;
    
for index=1:length(blinks); 
   
 if index==length(blinks); 
     
     if     blinks(index-1)+1 == blinks(index);
            
            blinks(index,count)=1;
                                   
     else 
     blinks(index,count)=0;       
                    

     end
     
 else          
        if  blinks(index+1) == blinks(index)+1;
            
            blinks(index,count)=1;
                      
             
        else 
            
            blinks(index,count)=1;
            numberBl(count-1,1)=sum(blinks(:,count));
            count=count+1;
            

        end
end
end

numberBl(count-1,1)=sum(blinks(:,count));

indTracks.nbrofblinks{index1,1}=length(numberBl(:,1));
indTracks.blinklength{index1,1}=numberBl(:,1);

blinklengthHist=cat(1,blinklengthHist,numberBl(:,1));

clear numberBl
end


%% Plot results from Blink and Gap Calculation

figure('Position',[50 300 800 200],'name','Blink and dark time Stats')

binrange=0:50;%max(gaps);
[bincounts] = histc(gaps,binrange);

subplot(1,3,1)
bar(binrange,bincounts,'histc')
title('dark time (frames)')
xlabel('dark time (frames)');
ylabel('counts');
hold on;

subplot(1,3,2)
histogram(cell2mat(indTracks.nbrofblinks))
title('Number of blinks (frames)');
xlabel('count');
ylabel('number of blinks');

subplot(1,3,3)
histogram(blinklengthHist)
title('Blink length(frames)');
xlabel('count');
ylabel('blink length');

MeanGap=mean(gaps)
MedianGap=median(gaps)

MeanNoB=mean(cell2mat(indTracks.nbrofblinks))
MedianNoB=median(cell2mat(indTracks.nbrofblinks))

MeanBlinkLength=mean(blinklengthHist)
MedianBlinkLength=median(blinklengthHist)
 
%% %% Normalize each cluster to its center of mass, plot 2D Histogram

figure('Position',[1000 300 700 600],'name','XY spread of individual clusters')

m=floor(sqrt(length(indTracks.tracks)))+1;
allclustersCx=[];
allclustersCy=[];

clusterx=[];
clustery=[];

for track=1:length(indTracks.tracks);

clusterx=indTracks.tracks{track,1}(:,3);
clustery=indTracks.tracks{track,1}(:,4);

clusterxC=sum(clusterx)/length(clusterx);
clusteryC=sum(clustery)/length(clustery);
clusterx=clusterx-clusterxC;
clustery=clustery-clusteryC;

allclustersCx=vertcat(allclustersCx,clusterx);
allclustersCy=vertcat(allclustersCy,clustery);

clear clusterx clustery
             
end

% c=hist3([allclustersCx*100, allclustersCy*100],[20 20]);

subplot(2,2,1)
scatter(allclustersCx, allclustersCy)
axis([-100 100 -100 100])
title('XY scatter')
xlabel('x dimension [nm]');
ylabel('y dimension [nm]');
box on;

subplot(2,2,2)
hist3([allclustersCx, allclustersCy],[20 20])
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
title('2D Hist')

pdx=fitdist(allclustersCx(:),'normal')
pdy=fitdist(allclustersCy(:),'normal')

subplot(2,2,3)
histfit(allclustersCx, 20)
axis([-50 50 0 2000])
title('x dimension')
xlabel('x dimension [nm]');
ylabel('counts');

subplot(2,2,4)
histfit(allclustersCy, 20)
axis([-50 50 0 2000])
title('y dimension')
xlabel('y dimension [nm]');
ylabel('counts');



