function convert_time_to_seconds(nc_file)
    % convert_time_to_seconds - Converts the time variable 't' from 'minutes since 2021-09-10 18:50:00'
    % to 'seconds since 2021-09-01' in a NetCDF file.
    %
    % Syntax: convert_time_to_seconds(nc_file)
    %
    % Input:
    %   nc_file - The path to the NetCDF file (string).
    %
    % Example:
    %   convert_time_to_seconds('/path/to/your/file.nc');

    % Step 1: Open the NetCDF file in write mode
    ncid = netcdf.open(nc_file, 'NC_WRITE');
    
    try
        % Step 2: Read the 't' variable and its units
        varid_t = netcdf.inqVarID(ncid, 't');
        t = netcdf.getVar(ncid, varid_t);

        % Step 3: Get the units of the 't' variable
        units_t = netcdf.getAtt(ncid, varid_t, 'units');
        ref_date_str = '2021-09-10 18:50:00';
        ref_date = datetime(ref_date_str, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

        % Step 4: Convert 't' from minutes to datetime
        t_datetime = ref_date + minutes(t);

        % Step 5: Define new reference date '2021-09-01 00:00:00'
        new_ref_date = datetime('2021-09-01 00:00:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

        % Step 6: Convert 't' to seconds since '2021-09-01'
        t_seconds = seconds(t_datetime - new_ref_date);

        % Step 7: Update the 't' variable with new values in seconds
        netcdf.putVar(ncid, varid_t, t_seconds);

        % Step 8: Update the units of 't' to 'seconds since 2021-09-01 00:00:00'
        new_units = 'seconds since 2021-09-01 00:00:00';
        netcdf.putAtt(ncid, varid_t, 'units', new_units);

        disp('Conversion completed successfully.');
    catch ME
        % If an error occurs, display the error message
        disp('Error occurred while processing the NetCDF file:');
        disp(ME.message);
    end

    % Step 9: Close the NetCDF file
    netcdf.close(ncid);
end