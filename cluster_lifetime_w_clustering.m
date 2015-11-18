%% Clear Workspace

clear, close all, clc

%% Select file from TS output

filenameC1='cell1_A549_gain500_30ms_lowerPM_1_TS_filtered';                      % -->  transformed far red channel, i.e. from Trans_2D_after_TS.m
filename_peaks1=[filenameC1 '.txt'];
locs=dlmread(filename_peaks1,',',1,0);

scatter(locs(:,2), locs(:,3), 1)

fprintf('\n -- Data Loaded --\n')

%% Select ROI for DBSCAN

%%%%%%%%%%%%%%%%%%% Select ROI %%%%%%%%%%%%%%%%%%%%%%%%% 

upperx= 3.4e4  %max(locs(:,2));
lowerx= 2e4 %0;

uppery= 2.5e4 %max(locs(:,3));
lowery=1e4 %0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vx=find(locs(:,2) < upperx & locs(:,2) > lowerx);
subset=locs(vx,1:8);

vy=find(subset(:,3) < uppery & subset(:,3) > lowery);
subset2=subset(vy,1:8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
scatter(subset2(:,2), subset2(:,3),1,'red'); hold on;

forDBS=[];
forDBS(:,1)=subset2(:,2);
forDBS(:,2)=subset2(:,3);

fprintf('\n -- ROI selected --\n')

%% Cluster DBSCAN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% input image --> subset2
tic

[class,type]=DBSCAN(forDBS,10,40);          % 
class2_1st=transpose(class);                % cluster number
type2=transpose(type);                      % cluster, border, outlyer

target_1st=find(type2 >= 0);
result_1st=subset2(target_1st,1:8);

figure
scatter(result_1st(:, 2), result_1st(:, 3),1,mod(class2_1st(target_1st),10));hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target_2nd=find(type2 < 0);
subset_2nd=subset2(target_2nd,:);

forDBS_2nd=[];
forDBS_2nd(:,1)=subset_2nd(:,2);
forDBS_2nd(:,2)=subset_2nd(:,3);

[class,type]=DBSCAN(forDBS_2nd,10,35);
class2_2nd=transpose(class);
type2=transpose(type);

target_2nd=find(type2 >= 0);
result_2nd=subset_2nd(target_2nd,1:8);

figure
scatter(result_2nd(:, 2), result_2nd(:, 3),1,class2_2nd(target_2nd),'black'); hold on

toc

 figure
 scatter(result_1st(:, 2), result_1st(:, 3),1,mod(class2_1st(target_1st),10));hold on
 scatter(result_2nd(:, 2), result_2nd(:, 3),1,mod(class2_2nd(target_2nd),10)); hold on
 
fprintf(' -- DBSCAN computed in %f sec -- \n',toc)

 %% Generate locs variable
 
 locs1=[];
 locs2=[];
 
 
 locs1(:,1)=result_1st(:,2);             % X
  locs1(:,2)=result_1st(:,3);            % Y
   locs1(:,3)=class2_1st(target_1st);    % cluster number
    locs1(:,4)=result_1st(:,1);          % frame number
 
 locs2(:,1)=result_2nd(:,2);             % X
  locs2(:,2)=result_2nd(:,3);            % Y
   locs2(:,3)=class2_2nd(target_2nd);    % cluster number
    locs2(:,4)=result_2nd(:,1);          % frame number
 
%% Plot localization counter for the selected cluster

figure('Position',[500 600 300 300])
ind_cluster={};
lifetime=[];


for m=1:max(locs1(:,3)); % for all clusters
    
% generate count variabe and put 0 if there is no localization in the cluster 

count=zeros(length(subset2),1); 
% count(int16(subset2(find(locs(:,3)==n),4)))=1;
target=find(locs1(:,3)==m);
cluster=locs1(target,1:4);
count(int16(cluster(:,4)))=1;

% generate life variabe. If frame is positive add 1 to the localization
% counter

life=zeros(length(subset2),1);
life(1,1)=1;                                    % localization conter
life(1,2)=min(cluster(:,4));                    % frame number

for n=2:length(count);
    
    if  count(n)==1;
        life(n,1)=life(n-1,1)+1;
        life(n,2)=n;
    else
        life(n)=life(n-1);
        life(n,2)=n;
    end
    
end
for n=2:length(count);
    
    if  count(n)==1;
        life(n,1)=life(n-1,1)+1;
        life(n,2)=n;
    else
        life(n)=life(n-1);
        life(n,2)=n;
    end
    
end


plot(life(:,2)*0.04,life(:,1)); hold on;

ind_cluster{m}=life;
lifetime=cat(1,lifetime,length(cluster));

clear target cluster

end


for p=1:max(locs2(:,3)); % for all clusters
    
% generate count variabe and put 0 if there is no localization in the cluster 

count=zeros(length(subset2),1); 
% count(int16(subset2(find(locs(:,3)==n),4)))=1;
target=find(locs2(:,3)==p);
cluster=locs1(target,1:4);
count(int16(cluster(:,4)))=1;

% generate life variabe. If frame is positive add 1 to the localization
% counter

life=zeros(length(subset2),1);
life(1,1)=1;                                    % localization conter
life(1,2)=min(cluster(:,4));                    % frame number

for n=2:length(count);
    
    if  count(n)==1;
        life(n,1)=life(n-1,1)+1;
        life(n,2)=n;
    else
        life(n)=life(n-1);
        life(n,2)=n;
    end
    
end
for n=2:length(count);
    
    if  count(n)==1;
        life(n,1)=life(n-1,1)+1;
        life(n,2)=n;
    else
        life(n)=life(n-1);
        life(n,2)=n;
    end
    
end

ind_cluster{m+p}=life;
lifetime=cat(1,lifetime,length(cluster));

plot(life(:,2)*0.04,life(:,1)); hold on;
title(['Cluster Time Evolution n=', num2str(m+p)])
xlabel('time (s)');
ylabel('cumulative count');
clear target cluster
end

figure('Position',[1000 600 300 300])
hist(lifetime);
title(['Lifetime Histogramm n=', num2str(m+p)])
xlabel('Cluster lifetime (sec)');
ylabel('count');
