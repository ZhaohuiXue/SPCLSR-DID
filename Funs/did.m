function [incremental_samples] = did(Indian, Indian_gt, win_size_change, similarity_threshold, predict_label,D_index)
    [~,~,neighbord_pixel,~] = spatialwin( Indian, Indian_gt, win_size_change, D_index, similarity_threshold);
    incremental_samples= zeros(size(predict_label,2), 2);
    for i = 1:size(predict_label,2)
        first_label = predict_label(i);
        neighbor_labels = predict_label(neighbord_pixel(i).id);
        if sum(  neighbor_labels == first_label ) == size(neighbor_labels, 2)
            incremental_samples(i,1) = 1;
            incremental_samples(i,2) = first_label;
        end    
    end
end

