function [ structure_mat] = create_structure_mat( Y, D, ~,D_coordinates,Y_coordinates)
    n1 = size(D,2);
    n2 = size(Y,2);
    structure_mat =zeros(n1,n2);
    structure_coods_mat =zeros(n1,n2);
    for i = 1:n1
        for j = 1:n2        
            structure_mat(i,j) = norm( D(:,i) - Y(:,j),2);  
            structure_coods_mat(i,j) = norm( D_coordinates(:,i) - Y_coordinates(:,j),2 );
        end
    end
    structure_mat = 1 - ( 1 - structure_mat / max(max(structure_mat)) ).^2;
    structure_coods_mat = structure_coods_mat / max(max(structure_coods_mat));
    structure_mat = structure_coods_mat.*structure_mat;
end

