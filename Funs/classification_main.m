function [ predict_label, res_value, res_struct_mat] = classification_main(Re_Y,Re_by,neighbord_pixel, num_class, D_label,struct_mat)
    Re_Y =Re_Y./repmat(sqrt(sum(Re_Y.^2,1)),size(Re_Y,1),1 ) ; 
    for i = 1:num_class
        Re_i_Y = Re_by(i).class;
        Re_i_Y =Re_i_Y./repmat(sqrt(sum(Re_i_Y.^2,1)),size(Re_i_Y,1),1 );
        res_value(i,:) = sum( abs(Re_i_Y - Re_Y) );      
    end
    res_struct_mat = ones(num_class, size(Re_Y,2));
    for i = 1: num_class
        temp_index = D_label==i;
        struct_mat_i = struct_mat(temp_index,:) ;
        if size(struct_mat_i,1) == 1
            res_struct_mat(i,:) = struct_mat_i;
        else
            res_struct_mat(i,:) = min(struct_mat_i);
        end
    end
    res_value = res_value .* res_struct_mat;  
    for i = 1:size(Re_Y,2)
        neighbord_index = neighbord_pixel(i).id;
        neighbord_res = res_value(:,neighbord_index);
        min_res = min( min(neighbord_res) );
        [label_row, ~] = find(neighbord_res == min_res);
        predict_label(i) = label_row(1);
    end  
end


