function numbers = extractSieveNumbers(inputTable)
    % Extract the column headers from the table
    columnHeaders = inputTable.Properties.VariableNames;
    
    % Initialize an array to hold the numbers
    numbers = [];
    
    % Define the regular expression pattern for matching 'Sieve_###mu_g'
    pattern = '^Sieve_(\d+)mu_g$';
    
    % Loop through the column headers
    for i = 1:length(columnHeaders)
        header = columnHeaders{i};
        
        % Check if the header matches the pattern
        tokens = regexp(header, pattern, 'tokens');
        
        if ~isempty(tokens)
            % Extract the number part and convert to double
            numStr = tokens{1}{1};
            number = str2double(numStr);
            
            % Append the number to the array
            numbers(end + 1) = number;
        elseif strcmp(header, 'Sieve_Pan_g')
            % If the header is 'Sieve_Pan_g', add 0 to the array
            numbers(end + 1) = 0;
        end
    end
end
