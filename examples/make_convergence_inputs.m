%
% Make input files for the convergence tests.
%
% Sumedh Joshi
% 13 June 2015.

% Assumes that sigma, the Gaussian width, is set to 25.0.

% Add some paths.
addpath( genpath( '/home/smj96/Dropbox/repos/spectral-element-method-acoustics' ) );

% Set some static inputs.
input.Lx = [-400, 400];
input.Ly = [-400, 400];
input.Lz = [-400, 400];
input.mu = 0.0;
input.t_final = 1.0;
input.dt = 1.0e3;                        % We're going to let the SEMA code pick the time-step.
input.report_every_n_steps = 1;
input.timesteps_between_writes = 100000; % Basically forget about this; the SEMA code will write the last time-step.

% Build the search space.
n = [ 4 6 8 10 12 14 16 ];
m = [ 5 10 ];

% Open a file stream for writing an execution script.
fid    = fopen( 'run_convergence_test.sh', 'w+' );
fprintf( fid, 'OMP_NUM_THREADS=1\n' );

% Loop over the search space building input files.
for ii = 1:length(n)
   for jj = 1:length(m)

      % Build the run name.
      iiname = [ 'convergence_N' num2str(n(ii)) '_M' num2str(m(jj)) ];

      % Set some filenames.
      input.fname_runname = iiname;
      infile_name = [iiname '_in'];

      % Set some grid values.
      input.n  = n(ii);
      input.mx = m(jj);
      input.my = m(jj);
      input.mz = m(jj);

      % Write the input file.
      sem_write_inputfile( infile_name, input );

      % Write the line to execute this case into the run script.
      fprintf( fid, [ '/usr/bin/mpiexec -np ' num2str( m(jj) ) ' ./sem_acoustics ' infile_name ' 2>' iiname '_err  | tee ' iiname '_log\n'] );

   end
end

% Close the execution script.
fclose(fid);
