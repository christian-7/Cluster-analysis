%% Sort by length

filename_peaks2='A647_EGF_10ms_1500mW_COT_Au__1_MMStack_locResults_DC.dat';
peaks=dlmread(filename_peaks2,',',1,0);
subsetLL=peaks;


%% Plot frame vs. ID after filtering by length

% x position
% y position
% cluster ID
% frame


figure
hold on;

for index=1:max(subset(:,3));
    
    cluster=find(subset(:,3)==index);
    subset(cluster,5)=length(cluster);
    
    if length(cluster)>142;
        
        scatter(subset(cluster,4)-min(subset(cluster,4)),subset(cluster,3),2);hold on;
        
    else end
    
    xlabel('frame');
    ylabel('cluster ID');
    
    clear cluster
    
end

%% Plot frame vs. ID

figure
hold on;

for index=1:max(subset(:,3));
    
    cluster=find(subset(:,3)==index);
  
    scatter(subset(cluster,4),subset(cluster,3),1);hold on;
        
    xlabel('frame');
    ylabel('cluster ID');
    
    clear cluster
    
end


%% Plot frame vs. Y coordinate

figure
hold on;

for index=1:max(subset(:,3));
    
    cluster=find(subset(:,3)==index);
  
    scatter(subset(cluster,4),subset(cluster,2),1);hold on;
        
    xlabel('frame');
    ylabel('y [nm]');
    
    clear cluster
    
end

%%  Plot cumulative number of locs

% figure
hold on;
frame=[];

for index=min(subset(:,4)):max(subset(:,4));
    
    target=find(subset(:,4)==index);
  
    if index==min(subset(:,4));
    
    frame(index,1)=length(target); % numbr of locs per frame
    frame(index,2)=index;
    
    else 
        
    frame(index,1)=length(target)+frame(index-1,1); % numbr of locs per frame
    frame(index,2)=index;    
    
    end
end
    
    
    scatter(frame(:,2),frame(:,1),1);hold on;
        
    xlabel('frame');
    ylabel('localizations');
    
    clear cluster
    



%% Plot Clusters 

figure

scatter(subset(:,1),subset(:,2),1,mod(subset(:,3),10))
title('all identified Clusters')
axis on
axis([min(subset(:,1)) max(subset(:,1)) min(subset(:,2)) max(subset(:,2))])  

figure

scatter(subset(:,1),subset(:,2),1,(subset(:,4)));
title('identified Clusters (color -> frame)')
colorbar
axis on
axis([min(subset(:,1)) max(subset(:,1)) min(subset(:,2)) max(subset(:,2))])  

%% 

figure
hold on;

for index=1:max(subset(:,3));
    
    cluster=find(subset(:,3)==index);
    subset(cluster,5)=length(cluster);
    
    if length(cluster)>140;
        
        scatter(subset(cluster,1),subset(cluster,2),1,mod(subset(cluster,3),10))
        title('identified Clusters after filter')
        axis on
        axis([min(subset(:,1)) max(subset(:,1)) min(subset(:,2)) max(subset(:,2))])  ;hold on;
        
    else end
    
    xlabel('frame');
    ylabel('cluster ID');
    
    clear cluster
    
end


figure
hold on;

for index=1:max(subset(:,3));
    
    cluster=find(subset(:,3)==index);
    subset(cluster,5)=length(cluster);
    
    if length(cluster)>140;
        
        scatter(subset(cluster,1),subset(cluster,2),1,subset(cluster,4))
        title('identified Clusters after filter (color -> frame)')
        axis on
        axis([min(subset(:,1)) max(subset(:,1)) min(subset(:,2)) max(subset(:,2))])  ;hold on;
        
    else end
    
    xlabel('frame');
    ylabel('cluster ID');
    
    clear cluster
    
end


    