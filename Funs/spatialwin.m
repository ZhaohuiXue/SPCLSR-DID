function [b_map, b_0_map, neighbord_pixel2,map] = spatialwin( map,gt_map, win_size, D_index, sim_T)
    [n1,n2,n3] = size(map);
    sample_index = find(gt_map~= 0);
    sample_index(D_index) = [];
    b_map = zeros(n1,n2);
    b_map(sample_index) = 1:size(sample_index,1);
    b_0_map = zeros(n1+win_size-1,n2+win_size-1);
    b_0_map( (win_size+1)/2:(win_size+1)/2+n1-1, (win_size+1)/2:(win_size+1)/2+n2-1 ) = b_map;
    for i = 1:n1   
        for j = 1:n2 
            if b_map(i,j) == 0 
                neighbord_pixel(i,j).id = 0;
                continue;
            end
            temp = b_0_map( i:i+win_size-1 , j: j+win_size-1);
            temp(temp == 0) = [];
            temp = reshape(temp, 1, size(temp,1)*size(temp,2));
            neighbord_pixel(i,j).id = temp;
        end
    end
    jsp = 1;
    for i = sample_index'
        neighbord_pixel2(jsp).id =  neighbord_pixel(i).id;
        jsp = jsp+1;
    end
    map = reshape(map,n1*n2, n3)';
    map = map./repmat(sqrt(sum(map.^2,1)),size(map,1),1 ) ;
    map = map(:,sample_index); 
    for i = 1:size(sample_index,1)
        joint_pixels_index =  neighbord_pixel2(i).id;
        joint_pixels = map(:,joint_pixels_index);
        center_pixel = map(:,i);
        temp_cos = center_pixel'* joint_pixels;
        temp_ed_1 = repmat(center_pixel, 1,size(joint_pixels,2));
        temp_ed_2 = joint_pixels - temp_ed_1;
        temp_ed_3 = sqrt(  sum(temp_ed_2 .* temp_ed_2)  );
        temp_ed = exp(-temp_ed_3);
        temp_ed_cos = temp_cos;
        temp_ed_cos(temp_ed_cos < sim_T) = 0;
        for j = 1:size(temp_ed_cos,2)    
            if temp_ed_cos(j) == 0
                neighbord_pixel2(i).id(j) =0;
            end        
        end
        neighbord_pixel2(i).id(neighbord_pixel2(i).id == 0) = [];
    end
end

