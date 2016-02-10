%
% Read output files for the convergence tests.
% Compare the answers with theory.
%
% Sumedh Joshi
% 13 June 2015.

% Assumes that sigma, the Gaussian width, is set to 25.0.

% Add some paths.
addpath( genpath( '/home/smj96/Dropbox/repos/spectral-element-method-acoustics' ) );

% Set some constants related to the theoretical solution.
sigma = 25.0;
t = 1.0;
c = 343.0;

% Build the search space.
n = [ 4 6 8 10 12 14 16 ]; %12 16 20 ];
m = [ 5 10 ];

% Loop over the search space building input files.
errL2 = zeros( length(n), length(m) );
errLi = zeros( length(n), length(m) );
L2    = zeros( length(n), length(m) );
Li    = zeros( length(n), length(m) );
for ii = 1:length(n)
   for jj = 1:length(m)

      % Build the run name.
      iiname = [ 'convergence_N' num2str(n(ii)) '_M' num2str(m(jj)) ];

      % If the data is available, read and process it.
      try

         % Load the data.
         data = sem_read_fieldfile( [ iiname '.h5' ], -1 );

         % Build the radial vector.
         x = data.grid.x; y = data.grid.y; z = data.grid.z;
         r = sqrt( x.^2 + y.^2 + z.^2 );

         % Build the exact solution.
         se = (r - c*t).*exp( -(r - c*t).^2 / sigma^2 ) ./ ( 2 * r );

         % Get the simulated data.
         s = data.field.s;

         % Get the norm of this data.
         L2(ii,jj) = norm( s(:) );
         Li(ii,jj) = max(abs(s(:)));

         % Compute the error in both L2 and Linf norm.
         errL2(ii,jj) = norm( s( ~isinf( se ) ) - se( ~isinf( se ) ) );
         errLi(ii,jj) = max( abs( s( ~isinf( se ) ) - se( ~isinf( se ) ) ) );

      end

   end
end

clear r s se x y z;
save('convergence_data');
