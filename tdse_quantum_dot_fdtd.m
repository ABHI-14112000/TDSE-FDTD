% TDSE for a 1D quantum dot using an FDTD-style finite-difference grid in MATLAB
% Method: Crank-Nicolson time stepping (unconditionally stable)
%
% You can run this file directly in MATLAB.

clear; clc; close all;

%% Physical constants (SI)
hbar = 1.054571817e-34;      % J*s
m0   = 9.1093837015e-31;     % kg
q    = 1.602176634e-19;      % C (J/eV conversion)

%% Simulation setup
m_eff = 0.067 * m0;          % Effective mass (GaAs-like)
L     = 80e-9;               % Domain length (80 nm)
Nx    = 800;                 % Number of spatial points
x     = linspace(-L/2, L/2, Nx).';
dx    = x(2)-x(1);

dt      = 1e-17;             % Time step (s)
Nt      = 5000;              % Number of time steps
nPlot   = 50;                % Plot every nPlot steps

%% Quantum dot potential (finite well)
% Dot width and barrier height
wellWidth = 20e-9;           % 20 nm
V0_eV     = 0.30;            % 0.30 eV barrier
V         = V0_eV*q * ones(Nx,1);
V(abs(x) <= wellWidth/2) = 0;

%% Initial wave packet (inside dot)
x0     = -5e-9;              % Initial center (m)
sigma  = 2e-9;               % Width (m)
k0     = 2.0e9;              % Central wave number (1/m)
psi    = exp(-(x-x0).^2/(2*sigma^2)) .* exp(1i*k0*x);

% Normalize
psi = psi / sqrt(trapz(x, abs(psi).^2));

%% Build Crank-Nicolson matrices: (I + i dt/(2hbar) H) psi^{n+1} = (I - i dt/(2hbar) H) psi^n
% Hamiltonian: H = -(hbar^2/2m) d2/dx2 + V
alpha = hbar^2/(2*m_eff*dx^2);

diagH = 2*alpha + V;
offH  = -alpha * ones(Nx-1,1);

% Dirichlet boundaries psi = 0 at ends
I = speye(Nx);
H = spdiags([offH diagH offH], -1:1, Nx, Nx);

A = I + 1i*dt/(2*hbar) * H;
B = I - 1i*dt/(2*hbar) * H;

% Enforce hard-wall boundaries directly in matrices
A(1,:) = 0; A(1,1) = 1;
A(end,:) = 0; A(end,end) = 1;
B(1,:) = 0; B(1,1) = 1;
B(end,:) = 0; B(end,end) = 1;
psi(1) = 0; psi(end) = 0;

% Pre-factorization for speed
[Lfac, Ufac, Pfac, Qfac, Rfac] = lu(A);

%% Figure setup
figure('Color','w');
subplot(2,1,1);
hProb = plot(x*1e9, abs(psi).^2, 'b', 'LineWidth', 1.5); hold on;
hPot  = plot(x*1e9, V/(q), 'r--', 'LineWidth', 1.2);
ylabel('|\psi|^2 (arb.) / V (eV)');
legend('|\psi|^2', 'V(x)', 'Location', 'best');
grid on; xlim([min(x), max(x)]*1e9);
title('TDSE in a quantum dot (1D finite well)');

subplot(2,1,2);
hRe = plot(x*1e9, real(psi), 'k', 'LineWidth', 1.2); hold on;
hIm = plot(x*1e9, imag(psi), 'm', 'LineWidth', 1.2);
xlabel('x (nm)'); ylabel('\psi components');
legend('Re(\psi)', 'Im(\psi)', 'Location', 'best');
grid on; xlim([min(x), max(x)]*1e9);

%% Time evolution
for n = 1:Nt
    rhs = B * psi;
    psi = Qfac * (Ufac \ (Lfac \ (Pfac * (Rfac \ rhs))));

    % Re-enforce boundaries and renormalize
    psi(1) = 0; psi(end) = 0;
    psi = psi / sqrt(trapz(x, abs(psi).^2));

    if mod(n, nPlot) == 0
        set(hProb, 'YData', abs(psi).^2);
        set(hRe,   'YData', real(psi));
        set(hIm,   'YData', imag(psi));
        subplot(2,1,1);
        title(sprintf('TDSE in a quantum dot (t = %.2f fs)', n*dt*1e15));
        drawnow;
    end
end

%% Diagnostics
probability = trapz(x, abs(psi).^2);
E_kin = real(trapz(x, conj(psi) .* (-(hbar^2/(2*m_eff)) * secondDerivative(psi, dx))));
E_pot = real(trapz(x, conj(psi) .* (V .* psi)));
E_tot = E_kin + E_pot;

fprintf('Final normalization integral: %.8f\n', probability);
fprintf('Final total energy: %.6e J (%.6f eV)\n', E_tot, E_tot/q);

%% Local function
function d2psi = secondDerivative(psi, dx)
    d2psi = zeros(size(psi));
    d2psi(2:end-1) = (psi(3:end) - 2*psi(2:end-1) + psi(1:end-2)) / dx^2;
    d2psi(1) = 0;
    d2psi(end) = 0;
end
