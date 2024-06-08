function [ Re_Y,D_i,Re_by ] = Re_HSI(D,D_label,Z,num_class)
    for i =1:num_class
        i_class_index = find(D_label == i);
        D_temp = zeros(size(D));
        D_temp(:,i_class_index) = D(:,i_class_index);
        D_i(i).class = D_temp;
        Re_by(i).class = D_temp*Z;
    end
    Re_Y = D*Z;
end
