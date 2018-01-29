classdef Material < handle
    properties (Access = private)
        % Free-space permittivity [F/m]
        eps0 = 8.854187817e-12;

        % Reduced Planck constant [J/s]
        hbar = 1.054571800e-34;

        % Electron charge [C]
        e0 = 1.6021766209e-19;

        % Free-electron rest mass [kg]
        m0 = 9.10938356e-31;
    end
    
    properties
        eff_mass  = 9.10938356e-31;  % Effective mass [kg]
        eps_r     = 1;  % Relative permittivity
        
        % Phonon parameters
        w_LO          = 0; % LO phonon frequency [rad/s]
        w_TO          = 0; % TO phonon frequency [rad/s]
        w_phonon_damp = 0; % Phonon damping frequency [rad/s]
        
        % Plasma parameters
        w_p           = 0; % Plasma frequency [rad/s]
        w_plasma_damp = 0; % Plasma damping frequency [rad/s]

        % Electron mobility for pure, undoped room-temperature material
        mobility   = 0; % cm^2 / (Vs)        
        
        % Caughey-Thomas mobility parameters
        mu_L_300   = 0; % cm^2 / (Vs)
        mu_min_300 = 0; % cm^2 / (Vs)
        N_ref_300  = 0; % cm^{-3}
        gamma_0    = 0;
        gamma_1    = 0;
        gamma_2    = 0;
    end
    
    methods
        function obj = set_eff_mass(obj, eff_mass)
            obj.eff_mass = obj.m0 * eff_mass;
        end

        %% Set phonon parameters (only needed for semiconductors)
        % Note that the parameters are specified as energies (meV)
        function obj = set_phonon_params(obj, E_LO, E_TO, E_phonon_damp)
            obj.w_LO          = E_LO          * 0.001 * obj.e0 / obj.hbar;
            obj.w_TO          = E_TO          * 0.001 * obj.e0 / obj.hbar;
            obj.w_phonon_damp = E_phonon_damp * 0.001 * obj.e0 / obj.hbar;
        end
        
        %% set plasma parameters (only if you don't have Caughey-Thomas params)
        function obj = set_plasma_params(obj, E_p, E_plasma_damp)
            obj.w_p           = E_p           * obj.e0 / obj.hbar;
            obj.w_plasma_damp = E_plasma_damp * 0.001 * obj.e0 / obj.hbar;
        end

        %% set Caughey-Thomas parameters
        function obj = set_caughey_thomas_params(obj, mu_L_300,   ...
                                                      mu_min_300, ...
                                                      N_ref_300,  ...
                                                      gamma_0,    ...
                                                      gamma_1,    ...
                                                      gamma_2)
            obj.mu_L_300   = mu_L_300;
            obj.mu_min_300 = mu_min_300;
            obj.N_ref_300  = N_ref_300;
            obj.gamma_0    = gamma_0;
            obj.gamma_1    = gamma_1;
            obj.gamma_2    = gamma_2;
        end              
    end
end