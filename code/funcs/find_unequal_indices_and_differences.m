function [unequal_indices, differences] = find_unequal_indices_and_differences(array1, array2)
    % Check if the sizes of the arrays are the same
    if ~isequal(size(array1), size(array2))
        error('Arrays must be of the same size.');
    end
    
    % Find logical arrays where elements are not equal and not NaN
    unequal_logical = ~(isequaln(array1, array2) | (isnan(array1) & isnan(array2)));
    
    % Find the indices where the arrays are not equal
    unequal_indices = find(unequal_logical);
    
    % Calculate the differences at the unequal indices
    differences = array1(unequal_logical) - array2(unequal_logical);
end
