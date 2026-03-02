% TDSE for a 1D quantum dot using an explicit FDTD leapfrog scheme in MATLAB
% -----------------------------------------------------------------------
% Equation: i*hbar*dpsi/dt = -(hbar^2/2m)*d2psi/dx2 + V(x)*psi
%
% FDTD form (real/imag split):
%   psi = psiR + i*psiI
%   dpsiR/dt = -(1/hbar) * [-(hbar^2/2m)*d2(psiI)/dx2 + V*psiI]
%   dpsiI/dt = +(1/hbar) * [-(hbar^2/2m)*d2(psiR)/dx2 + V*psiR]
%
% Time staggering:
%   psiR at n*dt, psiI at (n+1/2)*dt (leapfrog)
%
% This script is intentionally self-contained and easy to modify.

clear; clc; close all;

%% Physical constants (SI)
hbar = 1.054571817e-34;      % J*s
m0   = 9.1093837015e-31;     % kg
q    = 1.602176634e-19;      % J/eV

%% Material and domain
mEff = 0.067*m0;             % GaAs-like effective mass
L    = 120e-9;               % Domain length (m)
Nx   = 1200;                 % Number of grid points
x    = linspace(-L/2, L/2, Nx).';
dx   = x(2) - x(1);

%% Time step (explicit stability guidance)
% Practical choice: keep dt comfortably below the kinetic scale hbar/(2*alpha)
% where alpha = hbar^2/(2m*dx^2). Use a safety factor < 1.
alpha       = hbar^2/(2*mEff*dx^2);
dtStability = hbar/(2*alpha);
safety      = 0.20;
dt          = safety*dtStability;
Nt          = 12000;
plotEvery   = 120;

%% Quantum-dot potential (finite square well in a barrier)
wellWidth = 20e-9;           % m
V0eV      = 0.30;            % eV
V         = V0eV*q*ones(Nx,1);
V(abs(x) <= wellWidth/2) = 0;

%% Optional absorbing layer near boundaries (reduces reflections)
absWidth = 15e-9;            % absorbing layer width at each edge
etaMax   = 8e13;             % 1/s damping strength
eta      = absorbingProfile(x, L, absWidth, etaMax);

%% Initial wave packet (inside/near dot)
x0    = -6e-9;
sigma = 2.2e-9;
k0    = 1.7e9;
psi0  = exp(-((x-x0).^2)/(2*sigma^2)) .* exp(1i*k0*x);
psi0  = psi0 / sqrt(trapz(x, abs(psi0).^2));

% Leapfrog state variables
psiR = real(psi0);                                % at t = n*dt
psiI = imag(psi0);                                % approx at t = n*dt initially

% Kick imag part by half-step to enforce staggering
lapR = secondDerivative(psiR, dx);
psiI = psiI + 0.5*dt*( (alpha/hbar)*lapR - (V/hbar).*psiR );

% Hard-wall endpoints + absorber on top
psiR([1 end]) = 0;
psiI([1 end]) = 0;

%% Plot setup
figure('Color','w');
subplot(2,1,1);
hProb = plot(x*1e9, psiR.^2 + psiI.^2, 'b', 'LineWidth', 1.4); hold on;
hPot  = plot(x*1e9, V/q, 'r--', 'LineWidth', 1.1);
ylabel('|\psi|^2 (arb.) / V (eV)');
legend('|\psi|^2', 'V(x)', 'Location', 'best');
grid on;

subplot(2,1,2);
hRe = plot(x*1e9, psiR, 'k', 'LineWidth', 1.1); hold on;
hIm = plot(x*1e9, psiI, 'm', 'LineWidth', 1.1);
xlabel('x (nm)');
ylabel('\psi components');
legend('Re(\psi)', 'Im(\psi)', 'Location', 'best');
grid on;

%% Time stepping (explicit leapfrog FDTD)
for n = 1:Nt
    % 1) update real part using imag part
    lapI = secondDerivative(psiI, dx);
    psiR = psiR - dt*( (alpha/hbar)*lapI - (V/hbar).*psiI );

    % 2) absorber damping for real part
    psiR = psiR .* exp(-eta*dt);

    % 3) update imag part using new real part
    lapR = secondDerivative(psiR, dx);
    psiI = psiI + dt*( (alpha/hbar)*lapR - (V/hbar).*psiR );

    % 4) absorber damping for imag part
    psiI = psiI .* exp(-eta*dt);

    % 5) boundary conditions
    psiR([1 end]) = 0;
    psiI([1 end]) = 0;

    if mod(n, plotEvery) == 0
        prob = psiR.^2 + psiI.^2;
        set(hProb, 'YData', prob);
        set(hRe,   'YData', psiR);
        set(hIm,   'YData', psiI);
        subplot(2,1,1);
        title(sprintf('1D TDSE FDTD (t = %.2f fs, dt = %.3g fs)', n*dt*1e15, dt*1e15));
        drawnow;
    end
end

%% Diagnostics
psi = psiR + 1i*psiI;
normInt = trapz(x, abs(psi).^2);
Ekin = real(trapz(x, conj(psi).*(-(hbar^2/(2*mEff))*secondDerivative(psi, dx))));
Epot = real(trapz(x, conj(psi).*(V.*psi)));
Etot = Ekin + Epot;

fprintf('dx = %.3e m, dt = %.3e s (stability ref %.3e s)\n', dx, dt, dtStability);
fprintf('Final normalization integral = %.8f\n', normInt);
fprintf('Final total energy = %.6e J (%.6f eV)\n', Etot, Etot/q);

%% Local helpers
function d2f = secondDerivative(f, dx)
    d2f = zeros(size(f));
    d2f(2:end-1) = (f(3:end) - 2*f(2:end-1) + f(1:end-2))/dx^2;
    d2f(1) = 0;
    d2f(end) = 0;
end

function eta = absorbingProfile(x, L, absWidth, etaMax)
    eta = zeros(size(x));
    xL = -L/2;
    xR =  L/2;

    leftIdx  = x < (xL + absWidth);
    rightIdx = x > (xR - absWidth);

    sL = ((x(leftIdx) - (xL + absWidth))/absWidth);
    sR = ((x(rightIdx) - (xR - absWidth))/absWidth);

    % Polynomial CAP profile
    eta(leftIdx)  = etaMax*(sL.^2);
    eta(rightIdx) = etaMax*(sR.^2);
end
