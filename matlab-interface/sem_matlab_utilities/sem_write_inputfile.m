function sem_write_inputfile( input_file_name, data )
% sem_write_inputfile( input_file_name, data )
%
% Writes the contents of data to an input file for the SEM code.
%
% Takes 2 arguments:
%
%   input_file_name - String indicating the field file to read from disk.
%   data            - Struct of variable name/values to write (c.f. as
%                     returned by sem_read_inputfile()).
%
% Returns nothing.
%
% 11 Nov 2014
% Sumedh Joshi

    % Specify the order of fields to write.  Any fields specified in data will
    % be written with these first, and in the order below, with the rest in a
    % random order.
    ordered_fields = { 'fname_runname', ...
                       'fname_init', ...
                       'fname_mesh', ...
                       'fname_setup', ...
                       'n', ...
                       'nsubx', ...
                       'nsubz', ...
                       'dt', ...
                       'tend', ...
                       'timesteps_between_writes', ...
                       'nu', ...
                       'nu_d', ...
                       'rho0', ...
                       'facrobin', ...
                       'facrobin_ppe', ...
                       'filter_order', ...
                       'gmres_maxit_poisson', ...
                       'gmres_maxit_viscous', ...
                       'gmres_tol_poisson', ...
                       'gmres_tol_viscous', ...
                       'gmres_restart_poisson', ...
                       'gmres_restart_viscous', ...
                       'bc_diffusion', ...
                       'bc_viscous_x', ...
                       'bc_viscous_z', ...
                       'read_bcs_from_initfile', ...
                       'check_null_error', ...
                       'check_numerical_error', ...
                       'force_direct_solve', ...
                       'do_interfacial_averaging', ...
                       'use_capacitance_preconditioner', ...
                       'use_deflation', ...
                       'use_parallel_gmres', ...
                       'read_from_setupfile', ...
                       'write_to_setupfile', ...
                       'setup_and_stop' };

    % Open a file stream.
    fid = fopen( input_file_name, 'w+' );

    % Write a time-stamp at the top of the file.
    fprintf( fid, '#####################\n' );
    fprintf( fid, '# Input file for SEM code.\n# Generated: %s\n', datestr( now ) );
    fprintf( fid, '#####################\n' );

    % Get a list of fields of the data struct and determine the longest so we
    % can left justify the output correctly.
    F          = fieldnames( data );
    pad_length = max( cellfun( @length, F ) );

    % Get a cell array of the values of the data struct.
    V = struct2cell( data );

    % Compute the order in which we walk through our input structure's fields
    % so that the output follows a logical ordering.  Any additional fields
    % that aren't in ordered_fields will be handled separately afterwards.
    ordered_indices = cellfun( @(x) find( strcmp( x, F ) ), ordered_fields, ...
                               'UniformOutput', false );
    empty_indices   = cellfun( @isempty, ordered_indices );
    ordered_indices = [ordered_indices{~empty_indices}];

    % Loop over the data/value pairs, writing to the input file.
    for ii = ordered_indices
        write_configuration( fid, F{ii}, V{ii}, pad_length );
    end

    % Remove all of the data/value pairs we just wrote so we can handle
    % any that remain.
    F(ordered_indices) = [];
    V(ordered_indices) = [];

    % Handle any configuration parameters that aren't in the ordered_fields
    % list above.
    for ii = 1:numel( F )
        write_configuration( fid, F{ii}, V{ii}, pad_length );
    end

    % Close the file.
    fclose( fid );

return

function write_configuration( fid, field, value, pad_length )
% write_configuration( fid, field, value, pad_length )
%
% Writes a single configuration field to the file descriptor specified.  The
% field name will be padded with spaces (that is, left justified) before the
% equals sign and value are written.  Care is taken to identify the value's
% type and write out a value appropriate for reading within the solver.
%
% Takes 4 values:
%
%   fid        - File descriptor to write the configuration to.
%   field      - String specifying the name of the configuration.
%   value      - Value to be written.
%   pad_length - Width of the field to write.  field will be left justified,
%                padded with spaces, with pad_length.
%
% Returns nothing.

    % NOTE: We check for logical values first since Octave treats logicals
    %       as numeric 1's and 0's.  If this case wasn't considered before
    %       a pure numeric, logicals would cause the SEM solver to bark
    %       when it parses the input file.
    if islogical( value )
        % This is a logical - print true or false.
        if value
            fprintf( fid, '%-*s = .true.\n', pad_length, field );
        else
            fprintf( fid, '%-*s = .false.\n', pad_length, field );
        end
    elseif isnumeric( value )
        % If this is just one number, print and move on.
        if length( value ) == 1
            fprintf( fid, '%-*s = %g\n', pad_length, field, value );
        elseif length( value ) == 2
            % This is an array of length two, which is probably
            % domain boundaries.
            fprintf( fid, '%-*s = %g,%g\n', ...
                     pad_length, field, value );
        else
            % This is an array, print it in Fortran format.
            %
            % NOTE: Currently this assumes that we're printing out the 4
            %       external boundaries.
            fprintf( fid, '%-*s = %g,%g,%g,%g\n', ...
                     pad_length, field, value );
        end
    elseif ischar( value )
        % If it's a string, just print it.
        fprintf( fid, '%-*s = %s\n', pad_length, field, value );
    end

return
