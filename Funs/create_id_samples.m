function [ID,ID_index,  ID_truth_labels, ID_predict_labels, ID_class,OA_mid_coff, OA_low_coff, OA_high_coff] = create_id_samples(incremental_samples_index ,D_class2 , ...
                                                    num_class,Y,Y_test, Y_test_label,predict_label, coff_min_max, Y_label, test_sample_index,ID_rate)
    
    ID_index = [];
    ID_label = [];
    ID_class = struct('index',[]);
    ID_class(num_class).index = [];
    ID_class_low = struct('index',[]);
    ID_class_low(num_class).index = [];
    ID_class_high = struct('index',[]);
    ID_class_high(num_class).index = [];
    for i = 1:size(incremental_samples_index,1)
        if incremental_samples_index(i,1) == 1
            temp_label = incremental_samples_index(i,2);
            ID_class(temp_label).index(end+1) = i;
        end
    end
    for i = 1:num_class
        i_D_index = D_class2(i).index;
        i_ID_index = ID_class(i).index;
        if sum(i_ID_index) > 0
            i_D = Y(:,i_D_index);
            i_ID = Y_test(:,i_ID_index);
            [~, i_coffident_index,i_coffident_index_low, i_coffident_index_high ] = cal_corrcoff(i_D,i_ID, coff_min_max);
            if sum(i_coffident_index) ~=0
                temp = ID_class(i).index(i_coffident_index);
                temp_low = ID_class(i).index(i_coffident_index_low);
                temp_high = ID_class(i).index(i_coffident_index_high);
                ID_class(i).index = temp;
                ID_class_low(i).index = temp_low;
                ID_class_high(i).index = temp_high;
            end
        end
    end
    [OA_mid_coff, OA_low_coff, OA_high_coff] = cal_high_low_OA(ID_class,ID_class_low, ID_class_high,  num_class,Y_test_label,predict_label);
    for i =1:num_class
            size_i_class_index = size(ID_class(i).index,2);
            ID_i_class_num = ceil( size_i_class_index * ID_rate );
            relect_element =  randperm(size_i_class_index);
            relect_element =  relect_element(1:ID_i_class_num);
            i_class_select = ID_class(i).index(relect_element);  
            ID_class(i).index = i_class_select;
            ID_class(i).index = [D_class2(i).index, test_sample_index( ID_class(i).index)];  
            ID_index = [ID_index, ID_class(i).index];   
    end
    ID = Y(:,ID_index);
    ID_truth_labels = Y_label(ID_index);
    Y_label(test_sample_index)=predict_label;
    ID_predict_labels = Y_label(ID_index);
end

%%
function [mean_coff, coffident_index,coffident_index_low, coffident_index_high ] = cal_corrcoff(D, ED, coff_min_max)
    D_del_mean = D - mean(D);
    ED_del_mean = ED - mean(ED);
    up_cof = ED_del_mean'*D_del_mean;
    D_sqrt_var = sum( D_del_mean.*D_del_mean ) ;
    ED_sqrt_var = sum(ED_del_mean.*ED_del_mean);
    down_cof = sqrt(ED_sqrt_var'*D_sqrt_var);
    coff = abs( up_cof ./ down_cof );
    if size(coff,2) >1 
        mean_coff = mean(transpose(coff));
    else
        mean_coff = transpose(coff);
    end
    coffident_index =[];
    coffident_index_low =[];
    coffident_index_low = [];
    coffident_index_high =[];
    coffident_index = [];
    if size(mean_coff,2) > 5
        num_mencoff_col = size(mean_coff, 2);
        [~, sort_index] = sort(mean_coff);
        start_num =  max(ceil(num_mencoff_col*coff_min_max(1)), 1);
        mid_num = floor(num_mencoff_col*coff_min_max(2));
        max_num = num_mencoff_col;
        coffident_index_low = sort_index(1:start_num);
        coffident_index = sort_index(start_num:mid_num);
        coffident_index_high = sort_index(mid_num:max_num);
    end
end

function [OA_mid_coff, OA_low_coff,OA_high_coff] = cal_high_low_OA(ID_class,ID_class_low, ID_class_high, num_class,Y_test_label,predict_label)
    low_index = [];
    mid_index =[];
    high_index = [];
    for i = 1:num_class
        mid_index = [mid_index, ID_class(i).index];
        low_index = [low_index, ID_class_low(i).index];
        high_index = [high_index, ID_class_high(i).index];
    end
    low_pred_label = predict_label(low_index);
    low_truth_label = Y_test_label(low_index);
    mid_pred_label = predict_label(mid_index);
    mid_truth_label = Y_test_label(mid_index);
    high_pred_label = predict_label(high_index);
    high_truth_label = Y_test_label(high_index);
    % disp('%%%%%%%%%%%%%%%%%%%%%%%%%%-%%%%%%%%%%%%%%%%%%%%%%%%%%-%%%%%%%%%%%%%%%%%%%%%%%%%%')
    % disp(['low_avg_nums: ', num2str(numel(low_truth_label))])
    % disp(['med_avg_nums: ', num2str(numel(mid_truth_label))])
    % disp(['hig_avg_nums: ', num2str(numel(high_truth_label))])
    OA_low_coff = sum(low_pred_label == low_truth_label)/numel(low_pred_label);
    OA_mid_coff = sum(mid_pred_label == mid_truth_label)/numel(mid_pred_label);
    OA_high_coff = sum(high_pred_label == high_truth_label)/numel(high_pred_label);

end