% Lattice Boltzmann Implementation for NERS 570 Final Project using the
% D2Q9 Lattice for a variety of flow physics problems. Currently problems
% supported are the Lid-Driven Cavity (stationary sides, moving top),
% Couette flow (moving top, static bottom, periodic sides), and channel flow
% (stationary tops, periodic sides)

clear, clc, close all

% Define knobs
flow_type = "Channel"; % Options: Lid_Driven_Cavity, Couette, Channel
Nx = 100; % Nodes in x
Ny = 100; % Nodes in y
Ma = 0.1; % Desired Mach number
Re = 100; % Desired Reynolds number
max_it = 50000; % Maximum number of iterations
tol = 1e-7; % Steady-state tolerance

% Initial condition
rho_init = 1.0;
ux_init = 0;
uy_init = 0;

% Derive characteristics of the flow physics
L = Ny - 1; % Length of the domain in lattice units
cs = 1 / sqrt(3); % Speed of sound
cs2 = 1 / 3; % Speed of sound squared
U_lid = Ma*cs; % Lid velocity
nu = (U_lid*L) / Re; % Kinematic viscosity
tau = ( nu/(cs2) ) + 0.5; % Relaxation parameter
omega = 1 / tau; 

% Add body force if flow type is channel
if flow_type == "Channel"
    Fx = 1e-4; % Model pressure gradient with a body force
else
    Fx = 0.0; % No body force otherwise
end

% Define Lattice - Start at rest, go east, and move counterclockwise
ex = [0, 1, 1, 0, -1, -1, -1, 0, 1];
ey = [0, 0, 1, 1, 1, 0, -1, -1, -1];
w = [4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36]; % Weights for Maxwellian Distro

% Allocate Macroscopic fields
rho = ones(Ny, Nx) * rho_init;
ux = ones(Ny, Nx) * ux_init;
uy = ones(Ny, Nx) * uy_init;

% Allocate distribution functions
f = zeros(Ny, Nx, 9);
feq = zeros(Ny, Nx, 9);
fnew = zeros(Ny, Nx, 9);

% Initialize to equilibrium
u_sq = ux.^2 + uy.^2;
for i = 1:9
   e_dot_u = ex(i)*ux + ey(i)*uy; % e_i dot u
   f(:, :, i) = w(i) * rho .* (1 + 3*e_dot_u + 4.5*e_dot_u.^2 - 1.5*u_sq);


end


% Boundary indices
ixL = 1;
ixR = Nx;
iyB = 1;
iyT = Ny;


% Start Update Loop
prevUx = ux;
it = 0; % Current iteration

while(it < max_it)

    % Update iteration
    it = it + 1;

    % Compute macroscopic quantities from distribution
    rho = sum(f, 3); % Compute Density
    ux  = sum(f .* reshape(ex, 1, 1, 9), 3) ./ rho; % Compute ux 
    uy  = sum(f .* reshape(ey, 1, 1, 9), 3) ./ rho; % Compute uy

    % Apply body force in x direction (if channel flow)
    ux = ux + (Fx * tau ./ rho); 

    % Collision
    u_sq = ux.^2 + uy.^2;
    for i = 1:9
        e_dot_u = ex(i)*ux + ey(i)*uy; % e_i dot u
        feq(:, :, i) = w(i) * rho .* (1 + 3*e_dot_u + 4.5*e_dot_u.^2 - 1.5*u_sq);
    end
    fcoll = f - omega * (f - feq); % BGK collision model



    % Advect
    f_stream = zeros(size(f));
    for i = 1:9 
        f_stream(:,:,i) = circshift(fcoll(:,:,i), [ey(i), ex(i)]); 
    end
    f = f_stream;


    if flow_type == "Lid_Driven_Cavity"
        % Left Wall Bounce Back BCs
        f(:, ixL, [9 2 3]) = f(:,ixL, [5 6 7]); % Left Boundary (new right = previous left)
    
        % Right Wall Bounce Back BCs
        f(:, ixR, [5 6 7]) = f(:,ixR, [9 2 3]); % Right Boundary (new left = previous right)
    end
   

    % Bottom Wall Bounce Back BCs
    f(iyB, :, [3 4 5]) = f(iyB,:, [7 8 9]); % Bottom Boundary (new upward = previous downward)

    % Top Wall Bounce Back BCs
    f(iyT, :, [7 8 9]) = f(iyT,:, [3 4 5]);


    if flow_type == "Lid_Driven_Cavity" || flow_type == "Couette" 
        % Enforce Moving Wall BC with velocity correction
        rho_top = sum(f(iyT,:, :), 3); 
        f(iyT,:,7) = f(iyT,:,7)  - 2*w(7).*rho_top.*(U_lid)/cs2;
        f(iyT,:,9) = f(iyT,:,9) + 2*w(9).*rho_top.*(U_lid)/cs2;
    end


    if mod(it,100)==0
        % after bounce-back + moving-wall correction
        rho = sum(f,3);
        ux  = sum(f .* reshape(ex,1,1,9),3) ./ rho;

        err = max(abs(ux(:)-prevUx(:))); 
        fprintf('it=%d, err=%.3e\n',it,err);
        if err<tol, break; end
        prevUx=ux;
    end



end








%% --- Contour and Velocity Vector Plot ---

% Compute velocity magnitude
u_mag = sqrt(ux.^2 + uy.^2);

% Create grid coordinates
[X, Y] = meshgrid(1:Nx, 1:Ny);

% Contour of velocity magnitude
figure('Name','Velocity Magnitude Contour');
contourf(X, Y, u_mag, 30, 'LineColor','none');
colorbar;
xlabel('x'); ylabel('y');
title('Velocity Magnitude Contour');
set(gca,'YDir','normal');


box on;
ax = gca;
ax.Position = [0.1, 0.15, 0.72, 0.75];   % leave space for colorbar
set(gca, 'LooseInset', max(get(gca,'TightInset'), [0.05 0.05 0.05 0.05]));

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [7, 4.5]);         % wider canvas for colorbar
set(gcf, 'PaperPosition', [0, 0, 7, 4.5]);

print('-dpdf', 'LDC_cont.pdf');





delta = 2;
x_vec = 1:delta:Nx;
y_vec = 1:delta:Ny;
figure;
quiver(x_vec, y_vec, ux(1:delta:end,1:delta:end), uy(1:delta:end,1:delta:end), 'AutoScale','on');
set(gca,'YDir','normal');
xlabel('x'); ylabel('y');
title('Velocity Vectors');


box on;
ax = gca;
ax.Position = [0.13, 0.15, 0.78, 0.8]; 
set(gca, 'LooseInset', max(get(gca,'TightInset'), [0, 0, 0, 0]));
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6, 4]);
set(gcf, 'PaperPosition', [0, 0, 6, 4]);
print('-dpdf', 'LDC_vec.pdf');


