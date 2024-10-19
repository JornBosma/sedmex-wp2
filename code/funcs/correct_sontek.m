function [t2, corrected_par2] = correct_sontek(t2, par2, par_name)
    % t2: time vector for instrument 2 (input)
    % par2: wave height data for instrument 2 (input)
    % t1, par1: time and wave height data for instrument 1 (used for correction)
    % t3, par3: time and wave height data for instrument 3 (used for correction)
    % Output: t2 (unchanged), corrected_par2 (NaN-filled data for instrument 2)
  
    dataPath = [filesep 'Volumes' filesep 'T7 Shield' filesep 'DataDescriptor' filesep 'hydrodynamics' filesep];
    instruPath1 = [dataPath 'ADV' filesep 'L2C4VEC' filesep 'tailored_L2C4VEC.nc'];
    instruPath3 = [dataPath 'pressuresensors' filesep 'L2C6OSSI' filesep 'tailored_L2C6OSSI.nc'];

    t1 = ncread(instruPath1, 't');
    t3 = ncread(instruPath3, 't');

    par1 = ncread(instruPath1, par_name);
    par3 = ncread(instruPath3, par_name);
    
    % Reference start date for time data
    t0 = datetime('2021-09-01 00:00:00','InputFormat','yyyy-MM-dd HH:mm:ss');
    
    % Convert relative time to datetime format
    time2 = t0 + seconds(t2);
    time1 = t0 + seconds(t1);
    time3 = t0 + seconds(t3);
    
    % Find the common time range across the three instruments
    common_time = intersect(time1, time2);
    common_time = intersect(common_time, time3);
    
    % Preallocate corrected_par2 with NaNs
    corrected_par2 = nan(size(par2));
    
    % Find indices where the time matches common_time
    [~, idx2] = ismember(common_time, time2);
    [~, idx1] = ismember(common_time, time1);
    [~, idx3] = ismember(common_time, time3);

    % Filter the data for common times
    filtered_par2 = par2(idx2);
    filtered_par1 = par1(idx1);
    filtered_par3 = par3(idx3);
    
    % Correct Instrument 2 data based on interpolation from Instrument 1 and 3
    for j = 1:length(common_time)
        val1 = filtered_par1(j);  % Instrument 1
        val2 = filtered_par2(j);  % Instrument 2
        val3 = filtered_par3(j);  % Instrument 3

        % Handle the logic based on availability of values
        if isnan(val1) && isnan(val3)
            % corrected_par2(j) = val2;  % If both are NaN, keep original value
            corrected_par2(idx2(j)) = NaN;  % If both are NaN, keep NaN
        elseif isnan(val1)
            corrected_par2(idx2(j)) = val3;  % If Instrument 1 is NaN, use Instrument 3
        elseif isnan(val3)
            corrected_par2(idx2(j)) = val1;  % If Instrument 3 is NaN, use Instrument 1
        else
            corrected_par2(idx2(j)) = (val1 + val3) / 2;  % Interpolate between Instrument 1 and 3
        end
    end
end