function eta = visibility()
global lambda img
lambda = 1;
%visiblity function initialization
[rows, cols]=size(img);
eta = zeros(rows,cols);
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
            visibility = sin(pi .* visibility./2./lambda);
            % l2 norm of the function f(x) = sin(pi*(x)/2lambda)
            % l1 norm of the function f(x) 
            % find the maximum of the I differences and then multiply it by variance of the neighborhood
            % do most manipulations here to find the best heuristic 
            eta(X, Y) = max(visibility);
        end
    end
end
eta = eta.* repmat(var(eta),rows,1);
eta = eta./max(max(eta));  %do normalization dividing by max value of the overall matrix
end
