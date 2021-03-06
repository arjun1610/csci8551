function edgedetectACO
%
% lekhak: varsh007@umn.edu
%
%
global eta alpha beta rho phi N L init ants_K
global Threshold pheromone img rows cols
 
alpha = 1;      
beta = 0.1; 
rho = 0.1;  
phi = 0.01; 
N = 10;
L = 400;
init = 0.0001;

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
pheromone = init.* ones(size(img));
% initializing it here due to the fact of second update rule
pheromone_init=pheromone;
eta = visibility();
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
for iterator = 1: N
    tic;
    for step_iterator = 1: L
        last_deposited_pheromone = zeros(rows, cols);
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
                value = (ant_vis(xy,1)-1)*cols + (ant_vis(xy,2));                
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
	    total_sum=(sum(sum((ant_visibility_trp.^alpha) .* (ant_pheromone_trp.^beta))));
            trp = (ant_visibility_trp.^alpha) .* (ant_pheromone_trp.^beta) ./ total_sum;
            value = find(cumsum(trp)>=rand(1),1);
            ant_next_row_idx = ant_vis(value,1);
            ant_next_col_idx = ant_vis(value,2);
            if isempty(length(ant_next_row_idx))
                ant_next_row_idx = Kth_ant_row;
                ant_next_col_idx = Kth_ant_col;
            end
            ants_position_coordinates(ant_iterator,1) = ant_next_row_idx;
            ants_position_coordinates(ant_iterator,2) = ant_next_col_idx;
            %the last step coordinate taken by the ant
            last_deposited_pheromone(ants_position_coordinates(ant_iterator,1), ants_position_coordinates(ant_iterator,2)) = eta(ants_position_coordinates(ant_iterator,1), ants_position_coordinates(ant_iterator,2));
            % record the new position into ant's movement coordinates
            % store the flattened coordinate
            ant_movement(ant_iterator,step_iterator) = (ants_position_coordinates(ant_iterator,1)-1)*cols + ants_position_coordinates(ant_iterator,2);
        % local update     
        % first part if i,j belongs to current ant 
        % second part otherwise 
        pheromone = ((1-rho).*pheromone + rho.*eta.*last_deposited_pheromone).*last_deposited_pheromone + pheromone.*(abs(1-last_deposited_pheromone));
        end % end of ant_iterator
        % global update
        pheromone = ((1-phi).*pheromone + phi.*pheromone_init);
    end % end of step_iterator
    % create the image at every iteration to see the difference
    fprintf('pheromone sum of iteration  %d: %2f\n', iterator, sum(sum(pheromone)));    
    Threshold = classifier(pheromone); 
    imwrite(uint8(abs((pheromone>=Threshold).*255-255)), gray(256), [filename '_edge_aco_' num2str(iterator) '.bmp'], 'bmp');
    timespent=toc;
    fprintf('timespent for %d iteration: %2f\n', iterator, timespent);
end % end of iterator
