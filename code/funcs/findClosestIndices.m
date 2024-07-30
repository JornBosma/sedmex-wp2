function indices = findClosestIndices(A, B)
    % Function to find the indices of the values in A closest to the values in B
    
    % Initialize an array to store the indices
    indices = zeros(size(B));
    
    % Loop through each element in B
    for i = 1:length(B)
        % Calculate the absolute differences between the current B element and all elements in A
        [~, idx] = min(abs(A - B(i)));
        
        % Store the index of the closest element in A
        indices(i) = idx;
    end
end
