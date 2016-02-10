function data = sem_read_inputfile( input_file_name )
% data = sem_read_fieldfile( field_file_name )
%
% Reads the contents of the input file used by the SEM solver code.%
%
% Takes 1 argument:
%
%   input_file_name - String indicating the field file to read from disk.
%
% Returns 1 value:
%
%   data - struct which contains all variables set by the input file.  All
%          of the field names will be in lower-case.

    % Open a file stream.
    fid = fopen( input_file_name );

    % Figure out how many lines there are.
    nlines = 0;
    while(1)
        junk = fgetl( fid );
        if junk == -1
            break;
        end

        nlines = nlines + 1;
    end

    % Loop over the lines, parsing data.
    field_names  = {};
    field_values = {};
    count        = 0;
    frewind( fid );
    for ii = 1:nlines

       % Get this line, and do some tidying of the string.
       iiline = fgetl( fid );
       iiline = lower( strtrim( iiline ) );

       % See if this line is empty.  If it is, skip it.
       if length( iiline ) == 0
           continue;
       end

       % See if this line is a comment.  If it is, skip it.
       if iiline(1) == '#'
           continue;
       end

       % Look for an equals sign, and split the string between the variable
       % name and its value.
       eqndx    = strfind( iiline, '=' );
       variable = strtrim( iiline(1:eqndx-1) );
       value    = strtrim( iiline(eqndx+1:end) );

       % See if there is a trailing comment on this line.  If so, strip it.
       hashndx = strfind( value, '#' );
       if ~isempty( hashndx )
           value = strtrim( value(1:hashndx-1) );
       end

       % Process the data to cast into the appropriate variable type.

       % See if the value is boolean.
       if strfind( value, '.true.' )
           value = true;
       end

       if strfind( value, '.false.' )
           value = false;
       end

       % If the value is not a filename, then see if it has numbers.
       if isempty( strfind( value, '_' ) )

           % Check to see if it has three commas.
           if length( strfind( value, ',' ) ) == 3

               % This is a vector of boundary condition flags.
               value2   = cell2mat( regexp( value, '\d+', 'match' ) );
               value    = zeros( 4, 1 );
               value(1) = str2num( value2(1) );
               value(2) = str2num( value2(2) );
               value(3) = str2num( value2(3) );
               value(4) = str2num( value2(4) );

           else
               % What we're left with now is just single digit data.
               if isstr(value)
                   value = str2num(value);
               end
           end
       end

       % Append the data.
       field_names{count+1}  = variable;
       field_values{count+1} = value;
       count                 = count + 1;
    end

    % Make the output struct.
    data = cell2struct( field_values, field_names, 2 );

    % Close the filestream.
    fclose( fid );

end
