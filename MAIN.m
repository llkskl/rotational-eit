%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code package produces the simulations and reconstructions in our
% article "Rotational electrical impedance tomography using electrodes with
% limited boundary coverage provides window for multimodal sensing".
%
% This MAIN file first creates a mesh for a phantom with two set of four
% electrodes symmetrically on both sides of the sample. Then inclusions are
% inserted into the phantom and forward problem is simulated. Inverse
% problem is solved in finer mesh from noisy data.
%
% Licenced GPL v. 3
% Olli Koskela, olli.koskela@tut.fi, 2017
% 
% This package is developed using EIDORS v 3.8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Before running this file you should have started EIDORS with startup.m
% provided in your EIDORS distribution folder.

%% Create measurement geometry and phantoms
% electrode geometry and stimulation and measurement patterns
electrode_positions = [-11.25-22.5 -11.25 11.25 11.25+22.5 180-11.25-22.5 180-11.25 180+11.25 180+11.25+22.5];
stim = EightElectrodesAllCombinations();
name = 'Symmetrical 8 ele';

% diamater of the, electrode size and aqueous layer width
diam = 2;
elect_size = 0.14;
max_mesh_elem_size = 0.03;
max_mesh_elem_size_rim = 0.03;
max_mesh_elem_size_elec = 0.03;
outerlayer_width_percent = 0.05; % the aqueous solution

% create mesh for forward simulation
emptyPhantom = CreateMesh( diam, elect_size, electrode_positions, ...
    max_mesh_elem_size, max_mesh_elem_size_rim, max_mesh_elem_size_elec, ...
    outerlayer_width_percent );

% homogenous phantom
emptyPhantom.fwd_model.stimulation = stim;
emptyPhantom.fwd_model.get_all_meas = 1;
emptyPhantom.name = 'empty';
emptyPhantom.elem_data(:) = 1;

% in homogenous phantom: small object and large object
TwoObjectPhantom = emptyPhantom;
select_fun = @(x,y,z) (x-0.4).^2+(y-0.4).^2<0.1^2;
memb_ship = elem_select(TwoObjectPhantom.fwd_model, select_fun);
for i = 1:numel(memb_ship)
    if memb_ship(i) > 0
        TwoObjectPhantom.elem_data(i) = 0.1;
    end
end
select_fun = @(x,y,z) (x- (-0.5)).^2+(y-0.5).^2<0.2^2;
memb_ship = elem_select(TwoObjectPhantom.fwd_model, select_fun);
for i = 1:numel(memb_ship)
    if memb_ship(i) > 0
        TwoObjectPhantom.elem_data(i) = 0.1;
    end
end
TwoObjectPhantom.name = 'TOP';

% plot both phantoms
figure('Name','Homogenous phantom');
show_fem(emptyPhantom);
figure('Name','Phantom');
show_fem(TwoObjectPhantom);

%% FORWARD SIMULATION

% rotational measurement positions
thetas = linspace(0,360,17);

% Homogenous data
datahmg_NF = RotationalEITmeasurement(emptyPhantom ,thetas);
% Data w/ inclusions
data_NF = RotationalEITmeasurement(TwoObjectPhantom ,thetas);

% Returned data is matrix of NxM where N is the amount measurements in
% total for all stimulations and M = numel(thetas).

%% INVERSION

% regularization parameter
regularization_parameter = 0.08;

% create inversion mesh
max_mesh_elem_size = 0.04;
max_mesh_elem_size_rim = 0.04;
max_mesh_elem_size_elec = 0.04;
outerlayer_width_percent = 0.05;
inverseMesh = CreateMesh( diam, elect_size, electrode_positions, ...
    max_mesh_elem_size, max_mesh_elem_size_rim, max_mesh_elem_size_elec, ...
    outerlayer_width_percent );

% inversion parameters
inverseMesh.solve = 'rotational_inv_solve_diff_GN_one_step';
inverseMesh.reconst_type = 'difference';
inverseMesh.jacobian_bkgnd.value = 1;
inverseMesh.hyperparameter.value = regularization_parameter; 
inverseMesh.type = 'inv_model';
inverseMesh.parameters.max_iterations = 10;
inverseMesh.RtR_prior = @prior_laplace; priorname = 'laplace';

% create mapping matrix
inverseMesh.thetas = thetas;
M = CreateRotationMatrixCW(inverseMesh,inverseMesh.thetas);
inverseMesh.rotationalmapping = M;

% cross correlation matrix (set to identiety matrix, no cross correlation)
nmeas = numel(data_NF);
inverseMesh.meas_icov = speye(nmeas);
inverseMesh.fwd_model.stimulation = stim;

%% add noise to data
        
% common error level
common_err = 1e-4*max(data_NF(:));
% componentwise error level
component_err = 1e-3;
% noisy datas
data = AddTwoComponentNoise(data_NF,common_err,component_err);
datahmg = AddTwoComponentNoise(datahmg_NF,common_err,component_err);

% Solve noise free case
%rec_NF = inv_solve(inverseMesh,datahmg_NF,data_NF);
% Solve noisy case
rec = inv_solve(inverseMesh,datahmg,data);

% show reconstruction
figure('Name','Reconstruction');
show_fem(rec);
