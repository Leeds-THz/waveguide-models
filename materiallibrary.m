function material_list = materiallibrary()

% Create an empty table to contain all the parameters
material_list = containers.Map;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SEMICONDUCTORS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
AlAs = Material;
AlAs.set_eff_mass(0.150);
AlAs.eps_r    = 8.48;
AlAs.set_phonon_params(49.78, 44.86, 0.99);
AlAs.set_caughey_thomas_params(410, 10, 1e17, -2.1, 0.0, 0.0);
material_list('AlAs') = AlAs;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% METALS %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set properties for Ag
Ag = Material;
Ag.set_plasma_params(9.01, 18.0);
material_list('Ag') = Ag;

%% Set properties for Al
Al = Material;
Al.set_plasma_params(14.8, 81.8);
material_list('Al') = Al;

%% Set properties for Au
Au = Material;
Au.set_plasma_params(9.02, 26.67);
material_list('Au') = Au;

%% Set properties for Co
Co = Material;
Co.set_plasma_params(3.97, 36.6);
material_list('Co') = Co;

%% Set properties for Cu
Cu = Material;
Cu.set_plasma_params(7.39, 9.08);
material_list('Cu') = Cu;

%% Set properties for Fe
Fe = Material;
Fe.set_plasma_params(4.09, 18.22);
material_list('Fe') = Fe;

%% Set properties for Mo
Mo = Material;
Mo.set_plasma_params(7.46, 51.08);
material_list('Mo') = Mo;

%% Set properties for Ni
Ni = Material;
Ni.set_plasma_params(4.88, 43.64);
material_list('Ni') = Ni;

%% Set properties for Pb
Pb = Material;
Pb.set_plasma_params(7.36, 202.1);
material_list('Pb') = Pb;

%% Set properties for Pd
Pd = Material;
Pd.set_plasma_params(5.46, 15.37);
material_list('Pd') = Pd;

%% Set properties for Pt
Pt = Material;
Pt.set_plasma_params(5.15, 69.2);
material_list('Pt') = Pt;

%% Set properties for Ti
Ti = Material;
Ti.set_plasma_params(2.52, 47.4);
material_list('Ti') = Ti;

%% Set properties for V
V = Material;
V.set_plasma_params(5.16, 60.6);
material_list('V') = V;

%% Set properties for W
W = Material;
W.set_plasma_params(6.41, 60.4);
material_list('W') = W;