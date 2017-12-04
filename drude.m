function epsr = drude(f, N_d, T)
% DRUDE Calculate the relative permittivity of GaAs at a given frequencies
% 
% epsr = DRUDE_GAAS(freq)
%
% PARAMETERS:
%   f   - Frequencies in THz (can be an array)
%   N_d - Doping density [cm^{-3}]
%   T   - Temperature [K]
%
% RETURNS:
%   epsr - Relative complex permittivity at each frequency
%
% See Chapter 13 "Optical Waveguides" in P. Harrison & A. Valavanis, 
%    "Quantum Wells, Wires & Dots", 4th Ed., Wiley (2016).
%
% (c) Alex Valavanis <a.valavanis@leeds.ac.uk>
%     University of Leeds, 2017

%% Specify fundamental constants

% Free-space permittivity [F/m]
eps0 = 8.854187817e-12;

% Reduced Planck constant [J/s]
hbar = 1.054571800e-34;

% Electron charge [C]
e0 = 1.6021766209e-19;

% Free-electron rest mass [kg]
m0 = 9.10938356e-31;

%% Specify material constants

% TODO: These should not be hardcoded GaAs values!

% Effective mass of GaAs [kg]
mass = 0.067 * m0;

% Relative permittivity of GaAs
epsr_bulk = 10.89;

% Phonon frequencies [rad/s]
w_LO = 36.22 * 0.001 * e0 / hbar;
w_TO = 33.32 * 0.001 * e0 / hbar;

% Phonon damping factor [rad/s]
w_phonon_damp = 0.30 * 0.001 * e0 / hbar;

% Caughey-Thomas mobility parameters for GaAs
mu_L_300   = 8500; % cm^2 / (Vs)
mu_min_300 = 800;  % cm^2 / (Vs)
N_ref_300  = 1e17; % cm^{-3}

% Caughey-Thomas scaling constants
gamma_CT   = [-2.2 -0.9 6.2];

%% Calculate useful input parameters
% Radiation frequency [rad/s]
w = f * 1e12 * 2*pi;

%% Calculate the phonon contribution to the relative permittivity
epsr_phonon = epsr_bulk * (w_LO^2 - w_TO^2) ./ (w_TO^2 - w.^2 - 1i*w*w_phonon_damp);

%% Find the plasma contribution to the relative permittivity

% Find the plasma frequency [rad/s]
% Note that N_d is converted to m^{-3} here
w_p = sqrt(N_d * 100^3 * e0^2 / (eps0 * epsr_bulk * mass));

% Find mobility using Caughey-Thomas model
% Note that the N_ref value can stay in cm^{-3}
mu_L   = mu_L_300   * (T/300)^gamma_CT(1);
mu_min = mu_min_300 * (T/300)^gamma_CT(2);
N_ref  = N_ref_300  * (T/300)^gamma_CT(3);

mobility = mu_min + (mu_L - mu_min) / (1 + sqrt(N_d / N_ref));

% Find the plama damping factor [rad/s]
w_plasma_damp = e0 / (mass * mobility);

% Find the plasma contribution to the relative permittivity
%w_p/(2*pi)
epsr_plasma = epsr_bulk * (1 - w_p^2 ./(w.^2 + 1i * w * w_plasma_damp));

% Find the total permittivity
epsr = epsr_phonon + epsr_plasma;