function main
 
close all; clear all; clc;
 
img= double(imread('cm.tif'))./255;
%img =imnoise(i,'gaussian');
[nrow, ncol] = size(img);
 
v = zeros(size(img));
v_norm = 0;
for rr =1:nrow
    for cc=1:ncol
        %definition of clique
        temp1 = [rr-2 cc-1; rr-2 cc+1; rr-1 cc-2; rr-1 cc-1; rr-1 cc; rr-1 cc+1; rr-1 cc+2; rr cc-1];
        temp2 = [rr+2 cc+1; rr+2 cc-1; rr+1 cc+2; rr+1 cc+1; rr+1 cc; rr+1 cc-1; rr+1 cc-2; rr cc+1];
 
        temp0 = find(temp1(:,1)>=1 & temp1(:,1)<=nrow & temp1(:,2)>=1 & temp1(:,2)<=ncol & temp2(:,1)>=1 & temp2(:,1)<=nrow & temp2(:,2)>=1 & temp2(:,2)<=ncol);
 
        temp11 = temp1(temp0, :);
        temp22 = temp2(temp0, :);
 
        temp00 = zeros(size(temp11,1));
        for kk = 1:size(temp11,1)
            temp00(kk) = abs(img(temp11(kk,1), temp11(kk,2))-img(temp22(kk,1), temp22(kk,2)));
        end
 
        if size(temp11,1) == 0
            v(rr, cc) = 0;
            v_norm = v_norm + v(rr, cc);
        else
            lambda = 1;
            temp00 = sin(pi .* temp00./2./lambda);
            v(rr, cc) = sum(sum(temp00));
            v_norm = v_norm + v(rr, cc);
        end
    end
end
v = v./v_norm; 
%v = v.*100;
% pheromone function initialization
p = 0.0001 .* ones(size(img));    
 
%paramete setting
alpha = 10;     
beta = 0.01;    
rho = 0.1;     
phi = 0.05;    
 
ant_total_num = round(sqrt(nrow*ncol));
ant_pos_idx = zeros(ant_total_num, 2); % record the location of ant
 
% initialize the positions of ants
temp = rand(ant_total_num, 2);
ant_pos_idx(:,1) = round(1 + (nrow-1) * temp(:,1)); %row index
ant_pos_idx(:,2) = round(1 + (ncol-1) * temp(:,2)); %column index
 
search_clique_mode = '8';
memory_length = 40;
 
% record the positions in ant's memory
ant_memory = zeros(ant_total_num, memory_length);
total_step_num = 300;
total_iteration_num = 4;
 
for iteration_idx = 1: total_iteration_num
 
    %record the positions where ant have reached in the last 'memory_length' iterations   
    delta_p = zeros(nrow, ncol);
 
    for step_idx = 1: total_step_num
 
        delta_p_current = zeros(nrow, ncol);
 
        for ant_idx = 1:ant_total_num
 
            ant_current_row_idx = ant_pos_idx(ant_idx,1);
            ant_current_col_idx = ant_pos_idx(ant_idx,2);
 
            % find the neighborhood of current position
            if search_clique_mode == '4'
                rr = ant_current_row_idx;
                cc = ant_current_col_idx;
                ant_search_range_temp = [rr-1 cc; rr cc+1; rr+1 cc; rr cc-1];
            elseif search_clique_mode == '8'
                rr = ant_current_row_idx;
                cc = ant_current_col_idx;
                ant_search_range_temp = [rr-1 cc-1; rr-1 cc; rr-1 cc+1; rr cc-1; rr cc+1; rr+1 cc-1; rr+1 cc; rr+1 cc+1];
            end
 
            %remove the positions out of the image's range
            temp = find(ant_search_range_temp(:,1)>=1 & ant_search_range_temp(:,1)<=nrow & ant_search_range_temp(:,2)>=1 & ant_search_range_temp(:,2)<=ncol);
            ant_search_range = ant_search_range_temp(temp, :);
 
            %calculate the transit prob. to the neighborhood of current
            %position
            ant_transit_prob_v = zeros(size(ant_search_range,1),1);
            ant_transit_prob_p = zeros(size(ant_search_range,1),1);
 
            for kk = 1:size(ant_search_range,1)
 
                temp = (ant_search_range(kk,1)-1)*ncol + ant_search_range(kk,2);
 
                if length(find(ant_memory(ant_idx,:)==temp))==0 %not in ant's memory
                    ant_transit_prob_v(kk) = v(ant_search_range(kk,1), ant_search_range(kk,2));
                    ant_transit_prob_p(kk) = p(ant_search_range(kk,1), ant_search_range(kk,2));
                else %in ant's memory  
                    ant_transit_prob_v(kk) = 0;
                    ant_transit_prob_p(kk) = 0;                   
                end
            end
 
            % if all neighborhood are in memory, then the permissible search range is RE-calculated.
            if (sum(sum(ant_transit_prob_v))==0) && (sum(sum(ant_transit_prob_p))==0)               
                for kk = 1:size(ant_search_range,1)
                    temp = (ant_search_range(kk,1)-1)*ncol + ant_search_range(kk,2);
                    ant_transit_prob_v(kk) = v(ant_search_range(kk,1), ant_search_range(kk,2));
                    ant_transit_prob_p(kk) = p(ant_search_range(kk,1), ant_search_range(kk,2));
                end
            end                       
 
            ant_transit_prob = (ant_transit_prob_v.^alpha) .* (ant_transit_prob_p.^beta) ./ (sum(sum((ant_transit_prob_v.^alpha) .* (ant_transit_prob_p.^beta))));      
 
            % generate a random number to determine the next position. 
            temp = find(cumsum(ant_transit_prob)>=rand(1), 1);
 
            ant_next_row_idx = ant_search_range(temp,1);
            ant_next_col_idx = ant_search_range(temp,2);
 
            if length(ant_next_row_idx) == 0
                ant_next_row_idx = ant_current_row_idx;
                ant_next_col_idx = ant_current_col_idx;
            end
 
            ant_pos_idx(ant_idx,1) = ant_next_row_idx;
            ant_pos_idx(ant_idx,2) = ant_next_col_idx;
 
            %record the delta_p_current
            delta_p_current(ant_pos_idx(ant_idx,1), ant_pos_idx(ant_idx,2)) = 1;
 
            % record the new position into ant's memory
            if step_idx <= memory_length
                ant_memory(ant_idx,step_idx) = (ant_pos_idx(ant_idx,1)-1)*ncol + ant_pos_idx(ant_idx,2);
            elseif step_idx > memory_length
                ant_memory(ant_idx,:) = circshift(ant_memory(ant_idx,:),[0 -1]);
                ant_memory(ant_idx,end) = (ant_pos_idx(ant_idx,1)-1)*ncol + ant_pos_idx(ant_idx,2);
 
            end
 
            %update the pheromone function
            p = ((1-rho).*p + rho.*delta_p_current.*v).*delta_p_current + p.*(abs(1-delta_p_current));
 
        end % end of ant_idx
 
        % update the pheromone function
               p = (1-phi).*p;
 
    end % end of step_idx
 
end % end of iteration_idx
 
T = func_seperate_two_class(p);
imshow(uint8(abs((p>=T).*255-255)));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Inner Function  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
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
