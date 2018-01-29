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
%     University of Leeds, 2017-18

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

% Electron charge [C]
e0 = 1.6021766209e-19;

%% Create lookup tables for material constants
material_list = materiallibrary();

%% Read material values from lookup tables
mass          = material_list(material).eff_mass;
epsr_bulk     = material_list(material).eps_r;
w_LO          = material_list(material).w_LO;
w_TO          = material_list(material).w_TO;
w_phonon_damp = material_list(material).w_phonon_damp;
                             
%% Calculate useful input parameters
% Radiation frequency [rad/s]
w = f * 1e12 * 2*pi;

%% Calculate the phonon contribution to the relative permittivity
epsr_phonon = epsr_bulk * (w_LO^2 - w_TO^2) ...
              ./ (w_TO^2 - w.^2 - 1i*w*w_phonon_damp);

%% Find the plasma contribution to the relative permittivity

% First, try to grab the plasma frequency and damping factor from the
% material library
w_p           = material_list(material).w_p;
w_plasma_damp = material_list(material).w_plasma_damp;

% If we couldn't find the plasma frequency, calculate it
if (w_p == 0)
    % Find the plasma frequency [rad/s]
    % Note that N_d is converted to m^{-3} here
    w_p = sqrt(N_d * 100^3 * e0^2 / (eps0 * epsr_bulk * mass));
end

% If we couldn't find the plasma frequency, calculate it from mobility
if (w_plasma_damp == 0)
    % Try to read the bulk, room-temperature, undoped value
    mobility = material_list(material).mobility;

    % If we don't have a fixed mobility value, calculate it from the
    % Caughey-Thomas model
    if (mobility == 0)
        % Read Caughey-Thomas parameters
        mu_L_300      = material_list(material).mu_L_300;
        mu_min_300    = material_list(material).mu_min_300;
        N_ref_300     = material_list(material).N_ref_300;
        gamma_0       = material_list(material).gamma_0;
        gamma_1       = material_list(material).gamma_1;
        gamma_2       = material_list(material).gamma_2;

        % Find mobility using Caughey-Thomas model
        % Note that the N_ref value can stay in cm^{-3}
        mu_L   = mu_L_300   * (T/300)^gamma_0;
        mu_min = mu_min_300 * (T/300)^gamma_1;
        N_ref  = N_ref_300  * (T/300)^gamma_2;
    
        mobility = mu_min + (mu_L - mu_min) / (1 + sqrt(N_d / N_ref))
    end

    % Find the plama damping factor [rad/s]
    w_plasma_damp = e0 / (mass * mobility);
end

% Find the plasma contribution to the relative permittivity
epsr_plasma = epsr_bulk * (1 - w_p^2 ./(w.^2 + 1i * w * w_plasma_damp));
    
% Find the total permittivity
epsr = epsr_phonon + epsr_plasma;