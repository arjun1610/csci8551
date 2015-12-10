function edgedetectACO
%
% lekhak: varsh007@umn.edu
%
% Input:
% gray image with a square size
%
close all;
clear;
clc;
% image loading
filename = 'lenac';
%for reading the RGB images 
img=rgb2gray(imread([filename '.png']));
%for reading bw images
% img=imread([filename '.png']);
img = double(img)./255;
[rows, cols] = size(img);

%visiblity function initialization
v = zeros(size(img));
for i =1:rows
    for j=1:cols
        %definition of clique
        temp1 = [i-2 j-1; i-2 j+1; i-1 j-2; i-1 j-1; i-1 j; i-1 j+1; i-1 j+2; i j-1];
        temp2 = [i+2 j+1; i+2 j-1; i+1 j+2; i+1 j+1; i+1 j; i+1 j-1; i+1 j-2; i j+1];
        temp0 = find(temp1(:,1)>=1 & temp1(:,1)<=rows & temp1(:,2)>=1 & temp1(:,2)<=cols & temp2(:,1)>=1 & temp2(:,1)<=rows & temp2(:,2)>=1 & temp2(:,2)<=cols);
        
        temp11 = temp1(temp0, :);
        temp22 = temp2(temp0, :);
        
        temp00 = zeros(size(temp11,1),1);
        for kk = 1:size(temp11,1)
            temp00(kk) = abs(img(temp11(kk,1), temp11(kk,2))-img(temp22(kk,1), temp22(kk,2)));
        end
        
        if size(temp11,1) == 0
            v(i, j) = 0; 
        else
            lambda = 1;
            temp00 = sin(pi .* temp00./2./lambda);
            %l2 norm of the function f(x) = sin(pi*(x)/2lambda)
	    %l1 norm of the function f(x) 
	    %find the maximum of the I differences and then multiply it by variance of the neighborhood
            v(i, j) = (std(temp00)^2)*max(temp00);
        end
    end
end

v = v./max(max(v));  %do normalization dividing by max value of the overall matrix
% pheromone matrix initialization
p = 0.0001 .* ones(size(img));
%initializing it here due to the fact of second update rule
p_init=p;

%parameter setting, 
alpha = 1;      
beta = 0.1; 
rho = 0.1;  
phi = 0.01; 


ant_total_num = round(sqrt(rows*cols));

ant_pos_idx = zeros(ant_total_num, 2); % record the location of ant

% initialize the positions of ants

%     rand('state', sum(clock));
%     temp = rand(ant_total_num, 2);
%     ant_pos_idx(:,1) = round(1 + (rows-1) * temp(:,1)); %row index
%     ant_pos_idx(:,2) = round(1 + (cols-1) * temp(:,2)); %column index

%initialize the position of ants at maximum value of visibility
%function
%    [~,index]=max(v);
%initialize the position of ants at minimum value of visibility
% [~,index]=min(v);
%initialize the position of ants at mean value of visibility function 
meanValue=mean(v,1);
[~,index]=min(abs(v-repmat(meanValue,rows,1)));
    ant_pos_idx(1:cols,2)= 1:cols;
%
%       temp=0:cols-1;
%     temp1=ones(1,cols)*cols;
%     ant_pos_idx(cols+1:end,1)=temp1-temp;
    ant_pos_idx(1:rows,1)=index;
%     ant_pos_idx(rows+1:end,2)=index;
search_clique_mode = '8';   %Figure 1

% define the memory length, the positions in ant's memory are
% non-admissible positions for the next movement
memory_length = 20;
% record the positions in ant's memory, convert 2D position-index (row, col) into
% 1D position-index
ant_memory = zeros(ant_total_num, memory_length);

% System setup
total_step_num = 900;
total_iteration_num = 7;

for iteration_idx = 1: total_iteration_num
    tic;
    %record the positions where ant have reached in the last 'memory_length' iterations
    delta_p = zeros(rows, cols);  
    for step_idx = 1: total_step_num
        delta_p_current = zeros(rows, cols);
        for ant_idx = 1:ant_total_num
            ant_current_row_idx = ant_pos_idx(ant_idx,1);
            ant_current_col_idx = ant_pos_idx(ant_idx,2);
            
            % find the neighborhood of current position
            if search_clique_mode == '4'
                i = ant_current_row_idx;
                j = ant_current_col_idx;
                ant_search_range_temp = [i-1 j; i j+1; i+1 j; i j-1];
            elseif search_clique_mode == '8'
                i = ant_current_row_idx;
                j = ant_current_col_idx;
                ant_search_range_temp = [i-1 j-1; i-1 j; i-1 j+1; i j-1; i j+1; i+1 j-1; i+1 j; i+1 j+1];
            end
            
            %remove the positions out of the image's range
            temp = find(ant_search_range_temp(:,1)>=1 & ant_search_range_temp(:,1)<=rows & ant_search_range_temp(:,2)>=1 & ant_search_range_temp(:,2)<=cols);
            ant_search_range = ant_search_range_temp(temp, :);
            
            %calculate the transit prob. to the neighborhood of current
            %position
            ant_transit_prob_v = zeros(size(ant_search_range,1),1);
            ant_transit_prob_p = zeros(size(ant_search_range,1),1);
            
            for kk = 1:size(ant_search_range,1)
                
                temp = (ant_search_range(kk,1)-1)*cols + ant_search_range(kk,2);
                
                if isempty(ant_memory(ant_memory(ant_idx,:)==temp)) %not in ant's memory
                    ant_transit_prob_v(kk) = v(ant_search_range(kk,1), ant_search_range(kk,2));
                    ant_transit_prob_p(kk) = p(ant_search_range(kk,1), ant_search_range(kk,2));
                else %in ant's memory
                    ant_transit_prob_v(kk) = 0;
                    ant_transit_prob_p(kk) = 0;
                end
            end
            
            % if all neighborhood are in memory, then the permissible search range is RE-calculated.
            if (sum(sum(ant_transit_prob_v))==0) || (sum(sum(ant_transit_prob_p))==0)
                for kk = 1:size(ant_search_range,1)
                    %  temp = (ant_search_range(kk,1)-1)*cols + ant_search_range(kk,2);
                    ant_transit_prob_v(kk) = v(ant_search_range(kk,1), ant_search_range(kk,2));
                    ant_transit_prob_p(kk) = p(ant_search_range(kk,1), ant_search_range(kk,2));
                end
            end
            
            ant_transit_prob = (ant_transit_prob_v.^alpha) .* (ant_transit_prob_p.^beta) ./ (sum(sum((ant_transit_prob_v.^alpha) .* (ant_transit_prob_p.^beta))));
            
            % generate a random number to determine the next position.
            rand('state', sum(100*clock));
            temp = find(cumsum(ant_transit_prob)>=rand(1), 1);
            
            ant_next_row_idx = ant_search_range(temp,1);
            ant_next_col_idx = ant_search_range(temp,2);
            
            if isempty(length(ant_next_row_idx))
                ant_next_row_idx = ant_current_row_idx;
                ant_next_col_idx = ant_current_col_idx;
            end
            
            ant_pos_idx(ant_idx,1) = ant_next_row_idx;
            ant_pos_idx(ant_idx,2) = ant_next_col_idx;
            
            %record the delta_p_current
            delta_p_current(ant_pos_idx(ant_idx,1), ant_pos_idx(ant_idx,2)) = 1;
            
            % record the new position into ant's memory
            if step_idx <= memory_length
                
                ant_memory(ant_idx,step_idx) = (ant_pos_idx(ant_idx,1)-1)*cols + ant_pos_idx(ant_idx,2);
                
            elseif step_idx > memory_length
                ant_memory(ant_idx,:) = circshift(ant_memory(ant_idx,:),[0 -1]);
                ant_memory(ant_idx,end) = (ant_pos_idx(ant_idx,1)-1)*cols + ant_pos_idx(ant_idx,2);
                
            end
             p = ((1-rho).*p + rho.*v.*delta_p_current).*delta_p_current + p.*(abs(1-delta_p_current));
            
        end % end of ant_idx
        
        % update the pheromone function        
        delta_p = (delta_p + (delta_p_current>0))>0;
        p = ((1-phi).*p.*delta_p + phi.*p_init);

                
    end % end of step_idx
    % create the image at every iteration to see the difference
    fprintf('pheromone sum of iteration  %d: %2f\n', iteration_idx, sum(sum(p)));
    fprintf('deltaP sum of iteration %d: %2f\n', iteration_idx, sum(sum(delta_p)));
    T = func_seperate_two_class(p); %eq. (13)-(21), Calculate the threshold to seperate the edge map into two class
    imwrite(uint8(abs((p>=T).*255-255)), gray(256), [filename '_edge_aco_' num2str(iteration_idx) '.bmp'], 'bmp');
    
    timespent=toc;
    fprintf('timespent for %d iteration: %2f\n', iteration_idx, timespent);
    
end % end of iteration_idx

% generate edge map matrix
% It uses pheromone function to determine edge?
fprintf('Done!\n');

function level = func_seperate_two_class(I)

I = I(:);

% STEP 1: Compute mean intensity of image from histogram, set T=mean(I)
[counts, N]=hist(I,256);
i=1;
mu=cumsum(counts);
T(i)=(sum(N.*counts))/mu(end);

% STEP 2: compute Mean above T (MAT) and Mean below T (MBT) using T from
% step 1
mu2=cumsum(counts(N<=T(i)));
MBT=sum(N(N<=T(i)).*counts(N<=T(i)))/mu2(end);

mu3=cumsum(counts(N>T(i)));
MAT=sum(N(N>T(i)).*counts(N>T(i)))/mu3(end);
i=i+1;
T(i)=(MAT+MBT)/2;

% STEP 3 to n: repeat step 2 if T(i)~=T(i-1)
Threshold=T(i);
while abs(T(i)-T(i-1))>=1
    mu2=cumsum(counts(N<=T(i)));
    MBT=sum(N(N<=T(i)).*counts(N<=T(i)))/mu2(end);
    
    mu3=cumsum(counts(N>T(i)));
    MAT=sum(N(N>T(i)).*counts(N>T(i)))/mu3(end);
    
    i=i+1;
    T(i)=(MAT+MBT)/2;
    Threshold=T(i);
end

% Normalize the threshold to the range [i, 1].
level = Threshold;
