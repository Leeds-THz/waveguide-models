function material_list = materiallibrary()

% Create an empty table to contain all the parameters
material_list = containers.Map;

%% Set properties for GaAs
GaAs = Material;
GaAs.set_eff_mass(0.067);
GaAs.eps_r    = 10.89;
GaAs.set_phonon_params(36.22, 33.32, 0.30);
GaAs.set_caughey_thomas_params(8500, 800, 1e17, -2.2, -0.9, 6.2);
material_list('GaAs') = GaAs;

%% Set properties for Si
Si = Material;
Si.set_eff_mass(0.156); %Assume electron mass at Gamma
Si.eps_r = 11.6964;
Si.set_caughey_thomas_params(460, 45, 2.23e17, -2.18, 3.2, 0.72);
material_list('Si') = Si;

%% Set properties for Au
Au = Material;
Au.set_plasma_params(9.02, 26.67);
material_list('Au') = Au;