function [D,D_label, D_class,D_index,test_sample_index,num_class] = sampling( Y,Y_label, rate)
    num_class = max(Y_label);
    D_index = [];
    for i = 1:num_class
        i_class_index = find( Y_label == i );
        size_i_class_index = size(i_class_index,2);
        if rate < 1
            D_i_class_num = ceil( size_i_class_index * rate );
        elseif  1 <= rate  && rate <= 0.5*size_i_class_index
            D_i_class_num = rate;
        else 
            D_i_class_num = ceil(0.5 * size_i_class_index);
        end
        relect_element =  randperm(size_i_class_index);
        relect_element = relect_element(1:D_i_class_num);
        i_class_select = i_class_index(relect_element);
        D_class(i).index = i_class_select;
        D_index = [D_index, i_class_select];
    end
    D = Y(:,D_index);
    D_label = Y_label(D_index);
    test_sample_index = 1:size(Y_label,2);
    test_sample_index(D_index) = [];
end

