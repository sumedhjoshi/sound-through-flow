%%
%
% Adjoint simulation for the sound speed in a 1D wave equation model with 
% a uniform flow. -
% 9/11/14.

% Set some grid parameters.
n    = 1024;
xlim = [-1000 1000];
xlim = [-1000 2000];
xlim = [-2000 2000];
n    = 1024 + 1024;

% Set some constants.
nadjoint = 48;
x0       = 0;

% Build the grid.
x = linspace(xlim(1),xlim(2),n)';
h = x(2) - x(1);

% Set some time stepping parameters.
dt   = h/1500/10;
tmax = 0.40;
t    = linspace(0,tmax,round(tmax/dt));

% Set the sound speed with a Gaussian bump.
c = 1500;
V = 500;
rho0 = 1000;
beta = rho0 * c.^2;

% Set some constants.
sigma = (xlim(2) - xlim(1))/25;
sigma = 40;

% Set some initial conditions.

   % Stationary pulse at zero.
   u0 = zeros(n,1);
   s0 = 1e-3*2*exp( -(x - x0).^2 / sigma^2);

    % Right moving pulse.
    %s00 = 1e-3;
    %s0 = s00*exp( -(x - x0).^2 / sigma^2 );
    %u0 = c.*s0;

    %
    % The hyperbolic conservation laws are:
    %
    %       1. dsdt = - div( u )
    %       2. dudt = - beta * grad(s) / rho0
    %
    % Translated to finite difference approximations, they are:
    %
    %       1. s(k+1) = s(k) - dt * div( u(k) )
    %       2. u(k+1) = u(k) - dt * beta * grad( s(k) ) / rho0
    %
    % Where beta is the fluid bulk modulus.

% Build the differentiation matrices.

    % First order central difference.
    D1 = diag(0*ones(n,1)) + diag(1*ones(n-1,1),1) - diag(1*ones(n-1,1),-1);
    D1 = (1/(2*h))*D1;
    sD1  = sparse(D1);

    % Build Dirichlet-Dirichlet and Neumann-Neumann versions of these
    % operators.
    sD1d = sD1;
    sD1n = sD1;
    sD1n(1,1) = 1/h;
    sD1n(1,2) = -1/h;
    sD1n(end,end) = 1/h;
    sD1n(end,end-1) = -1/h;

    %
    % Second order central difference.
    D2 = diag(-2*ones(n,1)) + diag(1*ones(n-1,1),1) + diag(1*ones(n-1,1),-1);
    D2 = (1/h^2)*D2;
    sD2 = sparse(D2);

% Assemble the data.
Ud = zeros(n,length(t));
Sd = zeros(n,length(t));
Sd = zeros(n,length(t));
Ud(:,1) = u0;
Sd(:,1) = s0;
A = sparse( [ speye(n,n) - dt * V * sD1d, -dt*sD1n; -dt*c.^2*sD1d, speye(n,n) - dt * V .* sD1n ] );
for ii = 2:length(t)
    iiphi = A*[Sd(:,ii-1); Ud(:,ii-1)];
    Sd(:,ii) = iiphi(1:n);
    Ud(:,ii) = iiphi(n+1:end);
end

% Set the initial guess for the adjoint simulation.
vk    = 100;
Vk(1) = vk;

% Begin the adjoint simulation.
G    = zeros(nadjoint,1);
Uadj = zeros(n,nadjoint);
Sadj = zeros(n,nadjoint);
Vk   = zeros(1,nadjoint);
Vk(1) = vk;
for jj = 1:nadjoint

    fprintf('Beginning adjoint simulation number %3.0i\n',jj);

    % Build the time-propagator matrix.
    A = sparse( [ speye(n,n) - dt * vk * sD1d, -dt*sD1n; -dt*diag(c.^2)*sD1d, speye(n,n)  - dt * diag(vk) * sD1n ]);

    % Solve the forward problem with the current guess.
    U = zeros(n,length(t));
    S = zeros(n,length(t));
    U(:,1) = u0;
    S(:,1) = s0;
    for ii = 2:length(t)
        iiphi = A*[ S(:,ii-1); U(:,ii-1) ];
        S(:,ii) = iiphi(1:n);
        U(:,ii) = iiphi(n+1:end);
    end

    % Compute the data-model mismatch functional G.
    phi  = [S(:,end); U(:,end)];
    phid = [Sd(:,end); Ud(:,end)];
    G(jj) = sum( (phi - phid).^2 );
    Uadj(:,jj) = U(:,end);
    Sadj(:,jj) = S(:,end);

    % Print the value of the data-model mismatch.
    fprintf('\t\t Data-model mismatch functional value: %8.4f \n',G(jj));

    % Construct the "initial" condition for the adjoint simulation.
    lambda        = zeros(2*n,length(t));
    lambdaN       = 2*[phi(:,end) - phid(:,end)];

    %
    % Filter with a butterworth filter.
    Wn = 0.08;
    Wn = 0.8;
    Fo = 16;
    [Bbutter Abutter] = butter(Fo,Wn);
%     lambdaN(1:n)      = filter(Bbutter,Abutter,lambdaN(1:n));
%     lambdaN(n+1:end)  = filter(Bbutter,Abutter,lambdaN(n+1:end)c);

    % Set the initial condition.
    lambda(:,end) = lambdaN;

    % Solve the adjoint problem.
    for ii = length(t)-1:-1:1
        lambda(:,ii) = A'*lambda(:,ii+1);
    end

    % Build the estimate of the gradient.
    gradG = 0;
    kernel = zeros(1,length(t));
    for kk = 1:length(t) % These bounds are made up; how do I find these?  What should they be?
        iiAc         = A*[S(:,kk); U(:,kk)];
        iiAc_phik    = [ -dt * sD1d * S(:,kk); -dt * sD1n * U(:,kk) ];
        gradG        = gradG + lambda(:,kk)'*iiAc_phik;
        kernel(:,kk) = lambda(:,kk)'*iiAc_phik;
    end
    keyboard;

    % Filter the gradient.
    gradGu = gradG;

    % Update the guess for background flow with the gradient.
    %step_size = 1500; % In meters/second.
    step_size = 100*G(jj);
    %step_size = 1000;
    vk = vk + step_size.*gradG';
    Vk(jj) = vk;
    if G(jj) < 0.1
        keyboard;
    end

end
