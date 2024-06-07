clear
clc
close all
addpath("Data\",'Funs\')
load Indian.mat
%% 
%
% "Structure-Prior-Constrained Low-Rank and Sparse Representation With Discriminative Incremental Dictionary for Hyperspectral Image Classification," in IEEE Transactions on Geoscience and Remote Sensing, vol. 62, pp. 1-19, 2024
%
%% Settings
sampling_rate = 5/100;         
incremental_rate = 20/100;                 
coff_min_max = [0.1, 0.9];   
lamda1 = 1;
lamda2 = 0.02;
win_size = 1;             
win_size_change = 3;      
win_size2 = max([win_size_change-4,3]);  
similarity_threshold  = 0.96;
iters = 100;   % 50-200
%% Main Program
% Preserve marked pixels
[Y,Y_label,Indian] = pre_processing(Indian, Indian_gt);
% Build dictionary set
[D,D_label,D_class,D_index,test_sample_index,num_class] = sampling(Y, Y_label, sampling_rate);
% Build test set
Y_test = Y(:,test_sample_index);  Y_test_label = Y_label(test_sample_index);
% Keep initial dictionary and test sample information
D2 = D;   D_label2 = D_label; D_class2 = D_class; D_index2 = D_index; test_sample_index2 = test_sample_index;
% Pre-classification and re-classification
for sq = 1:2
    if sq == 1
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%-Pre-classification in progress-%%%%%%%%%%%%%%%%%%%%%%%%%%')
    else
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%-Re-classification in progress -%%%%%%%%%%%%%%%%%%%%%%%%%%')
    end
    % Build structure matrix
    [D_coordinates,Y_test_coordinates] = create_coordinates(Indian_gt, D_index,Y, test_sample_index);
    struct_mat = create_structure_mat(Y_test,D,D_class,D_coordinates,Y_test_coordinates);
    % Calculation coefficient
    Z= lrsr(Y_test,D, lamda1,lamda2,iters,struct_mat);
    % Subspace data restoration
    [ Re_Y,D_i,Re_by] = re_hsi(D,D_label,Z,num_class);
    % Spatial sample
    [b_map, b_0_map, neighbord_pixel,map] = spatialwin( Indian, Indian_gt, win_size, D_index,0);
    win_size = win_size2;
    % Classification
    [predict_label, res_value,res_struct_mat] = classification_main(Re_Y,Re_by,neighbord_pixel, num_class, D_label , struct_mat);
    if sq==1
        predict_label2 = predict_label;
    end
    % Classification accuracy analysis
    final_predict_label = Y_label;
    final_predict_label(test_sample_index2) = predict_label2;
    final_predict_label(test_sample_index) = predict_label;
    final_predict_label(D_index2) =[];

    eva_OA = sum(final_predict_label == Y_test_label) / length(Y_test_label); 
    if sq == 1
        disp(['Pre-classification ACC = ',num2str(eva_OA)]);
    else
         disp(['Re-classification ACC = ',num2str(eva_OA)]);
    end
    % Pre-classification map
    temp_label = Y_label;
    temp_label(test_sample_index2)=final_predict_label;
    predict_map= Indian_gt;
    predict_map(predict_map~=0)= temp_label;
    figure(sq)
    imshow(color_picture(predict_map));
    % Termination condition
    if sq == 2
        break
    end
    % Discriminative incremental dictionary learning(DID)
    [incremental_samples_index] = did(Indian, Indian_gt, win_size_change, similarity_threshold, predict_label,D_index2);
    [ID,ID_index, ID_truth_labels, ID_preditc_labels,ID_class, OA_mid_coff, OA_low_coff, OA_high_coff] = create_id_samples(incremental_samples_index,D_class2 , ...
                                                            num_class,Y, Y_test,Y_test_label, predict_label, coff_min_max,Y_label, test_sample_index, incremental_rate);
    % dictionary update
    D = ID; D_index = ID_index; D_label = ID_preditc_labels;D_class = ID_class;   
    test_sample_index = 1:size(Y,2); test_sample_index(D_index) = []; Y_test = Y(:,test_sample_index);
end
