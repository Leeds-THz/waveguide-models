function epsr = drude(material, f, varargin)
% DRUDE Calculate relative material permittivity at a given frequencies
% 
% epsr = DRUDE(material, freq)
%
% REQUIRED PARAMETERS:
%   material - The name of the material ('Au', 'GaAs' etc)
%   f        - Frequencies in THz (can be an array)
%
% OPTIONAL PARAMETERS:
%   'doping'      - Doping density [cm^{-3}]
%   'temperature' - Temperature [K]
%
% RETURNS:
%   epsr - Relative complex permittivity at each frequency
%
% EXAMPLES:
%   epsr = DRUDE('Au', 3.5);
%     Find the permittivity of gold at 3.5 THz.
%
%   epsr = DRUDE('GaAs', 3, 'doping', 1e18, 'temperature', 200);
%     Find the permittivity of GaAs with 1e18 cm^{-3} doping at 3 THz at
%     a temperature of 200 K.
%
% See Chapter 13 "Optical Waveguides" in P. Harrison & A. Valavanis, 
%    "Quantum Wells, Wires & Dots", 4th Ed., Wiley (2016).
%
% (c) Alex Valavanis <a.valavanis@leeds.ac.uk>
%     University of Leeds, 2017

%% Check for suitable MATLAB version
if verLessThan('matlab', '8.2')
    error('MATLAB 2013b or higher is required');
end

%% Handle input arguments
parser = inputParser;

% Define REQUIRED function arguments
addRequired( parser, 'material', @ischar);
addRequired( parser, 'f',        @isnumeric);

% Define OPTIONAL function arguments 
default_doping      = 0;   % [cm^{-3}] (Assume an undoped material)
default_temperature = 300; % [K] (Assume room temperature)

addParameter(parser, 'doping',      default_doping,      @isnumeric);
addParameter(parser, 'temperature', default_temperature, @isnumeric);

% Parse all function arguments
parse(parser, material, f, varargin{:});
N_d = parser.Results.doping;
T   = parser.Results.temperature;

%% Specify fundamental constants

% Free-space permittivity [F/m]
eps0 = 8.854187817e-12;

% Reduced Planck constant [J/s]
hbar = 1.054571800e-34;

% Electron charge [C]
e0 = 1.6021766209e-19;

% Free-electron rest mass [kg]
m0 = 9.10938356e-31;

%% Create lookup tables for material constants

% List of known materials
materials = {'GaAs',...
             'Si',...
             'Au'};

% Effective mass [kg]
mass_table = containers.Map(materials, m0 * [0.067,... %GaAs
                                             0.156,... %Si, assume electron mass at Gamma
                                             0]);

% Relative permittivity
epsr_bulk_table = containers.Map(materials, [10.89,...   %GaAs
                                             11.6964,... %Si
                                             1]);

% Phonon frequencies [rad/s]
w_LO_table = containers.Map(materials,...
                            0.001 * e0 / hbar * [36.22,... %GaAs
                                                 0,...     %Si
                                                 0]);
w_TO_table = containers.Map(materials,...
                            0.001 * e0 / hbar * [33.32,... %GaAs
                                                 0,...     %Si
                                                 0]);

% Phonon damping factor [rad/s]
w_phonon_damp_table = containers.Map(materials,...
                                     0.001 * e0 / hbar * [0.30,... %GaAs
                                                          0,...    %Si
                                                          0]);

% Plasma frequencies [rad/s]
% For semiconductors, enter '0' to calculate these using the doping
w_p_table = containers.Map(materials,...
                           e0 / hbar * [0,... % GaAs
                                        0,... % Si
                                        9.02]);

% Plasma damping frequencies (if known) [rad/s]
% For semiconductors, enter '0' to calculate these using the
% Caughey-Thomas mobility model
w_plasma_damp_table = containers.Map(materials,...
                           0.001 * e0 / hbar * [0,... % GaAs
                                                0,... % Si
                                                26.67]);
                                            
% Caughey-Thomas mobility parameters for GaAs
% cm^2 / (Vs)
mu_L_300_table = containers.Map(materials,...
                                [8500,...  % GaAs
                                 460,...   % Si
                                 0]);

mu_min_300_table = containers.Map(materials,...
                                  [800,... % GaAs
                                   45,...  % Si
                                   0]);

% cm^{-3}
N_ref_300_table = containers.Map(materials,...
                                 [1e17,...    %GaAs
                                  2.23e17,... %Si
                                  0]);

% Caughey-Thomas scaling constants
gamma_0_table = containers.Map(materials,...
                                [-2.2,...  % GaAs
                                 -2.18,... % Si
                                 0]);

gamma_1_table = containers.Map(materials,...
                                [-0.9,... % GaAs
                                 3.2,...  % Si
                                 0]);

gamma_2_table = containers.Map(materials,...
                                [6.2,...  % GaAs
                                 0.72,... % Si
                                 0]);

%% Read material values from lookup tables
mass = mass_table(material);
epsr_bulk = epsr_bulk_table(material);
w_LO = w_LO_table(material);
w_TO = w_TO_table(material);
w_phonon_damp = w_phonon_damp_table(material);
mu_L_300 = mu_L_300_table(material);
mu_min_300 = mu_min_300_table(material);
N_ref_300 = N_ref_300_table(material);
gamma_0 = gamma_0_table(material);
gamma_1 = gamma_1_table(material);
gamma_2 = gamma_2_table(material);
w_p = w_p_table(material);
w_plasma_damp = w_plasma_damp_table(material);
                             
%% Calculate useful input parameters
% Radiation frequency [rad/s]
w = f * 1e12 * 2*pi;

%% Calculate the phonon contribution to the relative permittivity
epsr_phonon = epsr_bulk * (w_LO^2 - w_TO^2) ./ (w_TO^2 - w.^2 - 1i*w*w_phonon_damp);

%% Find the plasma contribution to the relative permittivity

if (w_p == 0)
    % Find the plasma frequency [rad/s]
    % Note that N_d is converted to m^{-3} here
    w_p = sqrt(N_d * 100^3 * e0^2 / (eps0 * epsr_bulk * mass));
    
    % Find mobility using Caughey-Thomas model
    % Note that the N_ref value can stay in cm^{-3}
    mu_L   = mu_L_300   * (T/300)^gamma_0;
    mu_min = mu_min_300 * (T/300)^gamma_1;
    N_ref  = N_ref_300  * (T/300)^gamma_2;
    
    mobility = mu_min + (mu_L - mu_min) / (1 + sqrt(N_d / N_ref));
    
    % Find the plama damping factor [rad/s]
    w_plasma_damp = e0 / (mass * mobility);
end

% Find the plasma contribution to the relative permittivity
epsr_plasma = epsr_bulk * (1 - w_p^2 ./(w.^2 + 1i * w * w_plasma_damp));
    
% Find the total permittivity
epsr = epsr_phonon + epsr_plasma;