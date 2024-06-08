function [ map, map_gt, Indian ] = pre_processing(map, map_gt )
% Preserve marked pixels
    Indian = map;
    Indian_gt = map_gt;
    for i =1: size(Indian,1)
        for j =1: size(Indian,2)
            if Indian_gt(i,j) == 0
                Indian(i,j,:) =0;
            end
        end
    end    
    [n1,n2,n3] = size(map);
    map = reshape(map, n1*n2, n3)';
    map_gt = reshape(map_gt, 1, n1*n2);
    background_index = find(map_gt == 0);
    map(:,background_index) = [];
    map_gt(background_index) = [];
end

