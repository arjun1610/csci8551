function edgedetectACO
%
% lekhak: varsh007@umn.edu
%
% Input:
% gray image with a square size
%
global eta visibility lambda pheromone alpha beta rho phi steps_M iterations L
global ants_K ant_movement T
close all;
clear;
clc;

%parameter setting, 
alpha = 1;      
beta = 0.1; 
rho = 0.1;  
phi = 0.01; 
iterations = 10;
L = 100;
% image loading
filename = 'lenac';
%for reading the RGB images 
img=rgb2gray(imread([filename '.png']));
%for reading bw images
% img=imread([filename '.png']);
img = double(img)./255;
[rows, cols] = size(img);
%initialization
% we have K ants in total, equal to the the width of the image, due to an assumption that each
% column will have one ant.
ants_K = cols;
ant_movement = zeros(ants_K, L);
ants_position_coordinates = zeros(ants_K, 2); 
% pheromone matrix initialization
pheromone = 0.0001 .* ones(size(img));
% initializing it here due to the fact of second update rule
pheromone_init=pheromone;
%visiblity function initialization
eta = zeros(size(img));
for X =1:rows
    for Y=1:cols        
        vis1 = [X-2 Y-1; X-2 Y+1; X-1 Y-2; X-1 Y-1; X-1 Y; X-1 Y+1; X-1 Y+2; X Y-1];
        vis2 = [X+2 Y+1; X+2 Y-1; X+1 Y+2; X+1 Y+1; X+1 Y; X+1 Y-1; X+1 Y-2; X Y+1];
        vis = find(vis1(:,1)>=1 & vis1(:,1)<=rows & vis1(:,2)>=1 & vis1(:,2)<=cols & vis2(:,1)>=1 & vis2(:,1)<=rows & vis2(:,2)>=1 & vis2(:,2)<=cols);        
        vis1_in_range = vis1(vis, :);
        vis2_in_range = vis2(vis, :);
        visibility = zeros(size(vis1_in_range,1),1);
        for xy = 1:size(vis1_in_range,1)
            % modulus of visibility function values, 8 values for 8
            % connectivity neighborhood
            visibility(xy) = abs(img(vis1_in_range(xy,1), vis1_in_range(xy,2))-img(vis2_in_range(xy,1), vis2_in_range(xy,2)));
        end
        
        if size(vis1_in_range,1) == 0
            eta(X, Y) = 0; 
        else
            lambda = 1;
            visibility = sin(pi .* visibility./2./lambda);
            % l2 norm of the function f(x) = sin(pi*(x)/2lambda)
            % l1 norm of the function f(x) 
            % find the maximum of the I differences and then multiply it by variance of the neighborhood
            % do most manipulations here to find the best heuristic 
            eta(X, Y) = (std(visibility)^2)*max(visibility);
        end
    end
end

eta = eta./max(max(eta));  %do normalization dividing by max value of the overall matrix

%initialize the position of ants at maximum value of visibility
%function
[~,index]=max(eta);
%initialize the position of ants at minimum value of visibility
% [~,index]=min(eta);
%initialize the position of ants at mean value of visibility function 
% meanValue=mean(eta,1);
% [~,index]=min(abs(eta-repmat(meanValue,rows,1)));
ants_position_coordinates(1:cols,2)= 1:cols;
% for twice the sqrt(M1xM2) ants
%     in_range_values=0:cols-1;
%     vis1=ones(1,cols)*cols;
%     ants_position_coordinates(cols+1:end,1)=vis1-in_range_values;
ants_position_coordinates(1:rows,1)=index;
%     ants_position_coordinates(rows+1:end,2)=index;

for iterator = 1: iterations
    tic;
    delta_p = zeros(rows, cols);  
    for step_iterator = 1: L
        delta_p_current = zeros(rows, cols);
        for ant_iterator = 1:ants_K
            Kth_ant_row = ants_position_coordinates(ant_iterator,1);
            Kth_ant_col = ants_position_coordinates(ant_iterator,2);
            % how to move from I,J to K,L - Neighborhood move
            X = Kth_ant_row;
            Y = Kth_ant_col;
            neighborhood = [X-1 Y-1; X-1 Y; X-1 Y+1; X Y-1; X Y+1; X+1 Y-1; X+1 Y; X+1 Y+1];
            ant_vis = neighborhood(neighborhood(:,1)>=1 & neighborhood(:,1)<=rows & neighborhood(:,2)>=1 & neighborhood(:,2)<=cols,:);
            %calculate the transit probability function as defined in the
            %paper
            ant_visibility_trp = zeros(size(ant_vis,1),1);
            ant_pheromone_trp = zeros(size(ant_vis,1),1);
            for xy = 1:size(ant_vis,1)
                value = (ant_vis(xy,1)-1)*rows + (ant_vis(xy,2));                
                if isempty(ant_movement(ant_movement(ant_iterator,:)==value))
                    ant_visibility_trp(xy) = eta(ant_vis(xy,1), ant_vis(xy,2));
                    ant_pheromone_trp(xy) = pheromone(ant_vis(xy,1), ant_vis(xy,2));
                end
            end
            if (sum(sum(ant_visibility_trp))==0) || (sum(sum(ant_pheromone_trp))==0)
                for xy = 1:size(ant_vis,1)
                    ant_visibility_trp(xy) = eta(ant_vis(xy,1), ant_vis(xy,2));
                    ant_pheromone_trp(xy) = pheromone(ant_vis(xy,1), ant_vis(xy,2));
                end
            end
            trp = (ant_visibility_trp.^alpha) .* (ant_pheromone_trp.^beta) ./ (sum(sum((ant_visibility_trp.^alpha) .* (ant_pheromone_trp.^beta))));
            
            value = find(cumsum(trp)>=rand(1),1);
            ant_next_row_idx = ant_vis(value,1);
            ant_next_col_idx = ant_vis(value,2);
            if isempty(length(ant_next_row_idx))
                ant_next_row_idx = Kth_ant_row;
                ant_next_col_idx = Kth_ant_col;
            end
            ants_position_coordinates(ant_iterator,1) = ant_next_row_idx;
            ants_position_coordinates(ant_iterator,2) = ant_next_col_idx;
            %record the delta_p_current
            delta_p_current(ants_position_coordinates(ant_iterator,1), ants_position_coordinates(ant_iterator,2)) = 1;
            % record the new position into ant's memory
            ant_movement(ant_iterator,step_iterator) = (ants_position_coordinates(ant_iterator,1)-1)*cols + ants_position_coordinates(ant_iterator,2);
        % local update     
        % first part if i,j belongs to current ant 
        % second part otherwise 
        pheromone = ((1-rho).*pheromone + rho.*eta.*delta_p_current).*delta_p_current + pheromone.*(abs(1-delta_p_current));
        end % end of ant_iterator
        % update the pheromone function globally       
        delta_p = (delta_p + (delta_p_current>0))>0;
        pheromone = ((1-phi).*pheromone.*delta_p + phi.*pheromone_init);
    end % end of step_iterator
    % create the image at every iteration to see the difference
    fprintf('pheromone sum of iteration  %d: %2f\n', iterator, sum(sum(pheromone)));
    fprintf('deltaP sum of iteration %d: %2f\n', iterator, sum(sum(delta_p)));
    T = classifier(pheromone); 
    imwrite(uint8(abs((pheromone>=T).*255-255)), gray(256), [filename '_edge_aco_' num2str(iterator) '.bmp'], 'bmp');
    timespent=toc;
    fprintf('timespent for %d iteration: %2f\n', iterator, timespent);
end % end of iterator
