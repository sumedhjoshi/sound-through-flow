%%
%
% 2D finite difference time domain simulation of the inhomogenous wave
% equation computed from the hyperbolic conservation laws derived
% previously.  - 3/18/2013 

%
% Set some grid parameters. 
n    = 150; 
xlim = [-250 250];
ylim = [-250 250];

%
% Set the bulk modulus. 
c    = 343; 
rho0 = 1.02; 
beta = rho0 * c^2; 

%
% Build the grid. 
x = linspace(xlim(1),xlim(2),n)'; 
y = linspace(ylim(1),ylim(2),n)';
[Y X] = meshgrid(y,x);
X = X(:); 
Y = Y(:);
hx = x(2) - x(1); 
hy = y(2) - y(1); 
Lx = hx * n; 
Ly = hy * n;

%
% Set some time stepping parameters. 
dt   = hx/c/15; 
tmin = 0.0;
% tmax = (Lx / 2) / c; 
tmax = (Lx) / c; 
t    = linspace(tmin,tmax,round((tmax-tmin)/dt)); 

%
% Set some constants.
sigma = (xlim(2) - xlim(1))/50; 

%
% Set some initial conditions. 
u0x = zeros(n*n,1); 
u0y = zeros(n*n,1); 
x0 = -0; 
y0 = 0;
s0 = exp( - ((X - x0).^2 + (Y - y0).^2 ) / sigma^2); 


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

%
% Build the differentiation matrices. 

    %
    % First order central difference. 
    D1  = diag(0*ones(n,1)) + diag(1*ones(n-1,1),1) - diag(1*ones(n-1,1),-1);
    D1x = sparse(kron(sparse(eye(n)),D1)*(1/(2*hx)));
    D1y = sparse(kron(D1,sparse(eye(n)))*(1/(2*hy)));
    
%     %
%     % Second order central difference.
%     D2  = diag(-2*ones(n,1)) + diag(1*ones(n-1,1),1) + diag(1*ones(n-1,1),-1);
%     D2x = kron(eye(n),D2)*(1/hx^2); 
%     D2y = kron(D2,eye(n))*(1/hy^2);

%
% Set the background field. 
Vx = zeros(n^2,length(t)); 
Vy = zeros(n^2,length(t));

%     %
%     % Build a complex potential flow. 
%     Z = X + 1i*Y; 
%     A = 4*1e-5; 
%     nexp = 3; 
%     W = A*Z.^nexp; 
%     Wz = nexp*A*Z.^(nexp-1); 

%
% Set some constants for the Rankine vortex.
Gamma = 500; 
Rv    = max(X(:))/20;
nu    = 10e-6;
x0 = 175;
theta = atan(Y./(X-x0));
R     = sqrt((X-x0).^2 + Y.^2);
U1    = Gamma.*R/(2*pi*Rv^2);
U2    = Gamma./(2*pi*R);
U(R<Rv) = U1(R<Rv); 
U(R>=Rv) = U2(R>=Rv); 
 
for ii = 1:length(t)
%    Vx(:,ii) = (Y); % Shear layer.     
%    Vx(:,ii) = sign(Y-500)*500; % Shifted shear layer. 
%    Vx(:,ii) = 900*sin(10*pi*Y/Ly); % Alternating shear layers.
%    Vy(:,ii) = 900*sin(10*pi*X/Lx); % Alternating shear layers. (totally not physical).
%    Vx(:,ii) =  900*sin(4*pi*(X+0)/Lx).*cos(4*pi*(Y+0)/Ly); % Steady state Taylor vortex.
%    Vy(:,ii) = -900*cos(4*pi*(X+0)/Lx).*sin(4*pi*(Y+0)/Ly); % Steady state Taylor vortex.
%     Vx(:,ii) = real(Wz); % Complex potential flow. 
%     Vy(:,ii) = imag(Wz); % Complex potential flow. 
%     Vx(:,ii) = -cos(theta).*(rho0./R).*(Gamma^2./(4*pi^2*R.^2).*exp(-2*R.^2/(4*nu*t(ii)))); 
%     Vy(:,ii) = sin(theta).*(rho0./R).*(Gamma^2./(4*pi^2*R.^2).*exp(-2*R.^2/(4*nu*t(ii))));
    Vx(:,ii) =  sin(theta).*U'.*sign(X-x0); % Rankine vortex.
    Vy(:,ii) =  cos(theta).*U'.*sign(X-x0);  % Rankine: vortex.  
end

% %
% % Display the velocity. 
% figure; 
% imagesc(x,y,reshape(Vx(:,1),n,n)'); 
% set(gca,'xtick',[]); 
% set(gca,'ytick',[]);
% title('u_{0x}'); 
% % print_graphics(gcf,'background_velocity_x',1,0,0,0);
% 
% figure; 
% imagesc(x,y,reshape(Vy(:,1),n,n)'); 
% set(gca,'xtick',[]); 
% set(gca,'ytick',[]);
% title('u_{0y}');
% % print_graphics(gcf,'background_velocity_y',1,0,0,0);

%
% Start time stepping. 
Ux = zeros(n^2,length(t)); 
Uy = zeros(n^2,length(t)); 
S  = zeros(n^2,length(t)); 
Ux(:,1) = u0x; 
Uy(:,1) = u0y; 
S(:,1)  = s0; 
for ii = 2:length(t)
    p(ii,length(t),1);
    
    %
    % Build the condensation.
    S(:,ii)  =  S(:,ii-1) - dt * ( D1x * Ux(:,ii-1) + D1y * Uy(:,ii-1) );% ...
                          %- dt * S(:,ii-1) * ( );
    
    %
    % Build the velocity in x. 
    Ux(:,ii) = Ux(:,ii-1) - dt * beta * D1x * S(:,ii-1) / rho0  ...
                          - dt * Vx(:,ii) .* (D1x * Ux(:,ii-1)) ... 
                          - dt * Vy(:,ii) .* (D1y * Ux(:,ii-1)) ...
                          - dt * ( (D1x * S(:,ii-1)).*Vx(:,ii) + (D1x * S(:,ii-1)).*Vx(:,ii) ).*Vx(:,ii) ...
                          - dt * S(:,ii-1) .* (Vx(:,ii) - Vx(:,ii-1)) / dt; 

    %
    % Build the velocity in y. 
    Uy(:,ii) = Uy(:,ii-1) - dt * beta * D1y * S(:,ii-1) / rho0   ... 
                          - dt * Vx(:,ii) .* ( D1x * Uy(:,ii-1)) ... 
                          - dt * Vy(:,ii) .* ( D1y * Uy(:,ii-1)) ... 
                          - dt * ( (D1x * S(:,ii-1)).*Vx(:,ii) + (D1y * S(:,ii-1)).*Vy(:,ii) ).*Vy(:,ii) ...                          
                          - dt * S(:,ii-1) .* (Vy(:,ii) - Vy(:,ii-1)) / dt; 
end

%
% Build the velocity norm. 
U = Ux.^2 + Uy.^2;
U = (U).^(1/2);

%
% Get the center velocity and condensation.
uxcenter = Ux(logical((X==min(abs(X))).*(Y==min(abs(Y)))),:);
uycenter = Uy(logical((X==min(abs(X))).*(Y==min(abs(Y)))),:);
scenter  = S(logical((X==min(abs(X))).*(Y==min(abs(Y)))),:);

%
% Set the minimum and maximum values. 
smin = min(S(:)); 
smax = max(S(:));
umin = min(U(:,1));
umax = max(U(:,1));

keyboard; 

%
% Animate the velocity. 
figure;
filename = 'acoustic.gif';
for k = 1:4:length(t)
      imagesc(x,y,reshape(U(:,k),n,n)'); 
      set(gca,'ydir','normal');
      set(gcf,'color','w');
%       hline(500,'w--');
      xlabel('meters');
      ylabel('meters');
      caxis([umin umax]);            
      colorbar; 
      drawnow
      frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if k == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0);
      end
end

% %
% % Make a plot at four snapshot times. 
% figure; 
% nsnaps = 4; 
% for ii = 1:nsnaps
%     tndx = round(length(t) * (ii/nsnaps));
%     iiU = U(:,tndx); 
%     iiU = reshape(iiU,n,n)'; 
%     subplot(1,nsnaps,ii);
%     imagesc(x,y,iiU); 
%     set(gca,'ydir','normal');
%     caxis([umin,umax]); 
%     hline(500,'w--');
%     set(gca,'ytick',[]); 
%     set(gca,'xtick',[]); 
%     box on; 
% end
% set(gcf,'position',[ 23         402        1120         215]);
% set(gcf,'color','w');
% fixfig_pres;
% print_graphics(gcf,'shearflow_snapshots',1,0,0,0); 

% %
% % Make a 4x4 plot at 16 snapshot times. 
% figure; 
% nsnaps = 16; 
% for ii = 1:nsnaps
%     tndx = round(length(t) * (ii/nsnaps));
%     iiU = U(:,tndx); 
%     iiU = reshape(iiU,n,n)'; 
%     subplot(sqrt(nsnaps),sqrt(nsnaps),ii);
%     imagesc(x,y,iiU); 
%     set(gca,'ydir','normal');
%     caxis([umin,umax]); 
% %     hline(500,'w--');
%     set(gca,'ytick',[]); 
%     set(gca,'xtick',[]); 
%     box on; 
% end
% set(gcf,'position',[2374           2         985         950]); 
% set(gcf,'color','w');
% fixfig_pres;
% print_graphics(gcf,'gridplot',1,0,0,0); 






