function [D_coordinates, Y_coordinates] = create_coordinates(Indian_gt,D_index,Y,test_sample_index)
    Gt_map = Indian_gt;
    no_zero_nums = size( find(Gt_map~=0), 1);
    Gt_map(Gt_map~=0) = 1:no_zero_nums;
    D_coordinates = zeros(2, size(D_index,2));
    Y_coordinates = zeros(2, size(Y,2));
    [Y_coordinates(1,:), Y_coordinates(2,:)] = find(Gt_map~=0);
    for i=1:size(D_index,2)
        [temp_row, temp_column] = find(Gt_map == D_index(i));
        D_coordinates(1,i) = temp_row;
        D_coordinates(2,i) = temp_column;
    end
    Y_coordinates = Y_coordinates(:,test_sample_index);
end

