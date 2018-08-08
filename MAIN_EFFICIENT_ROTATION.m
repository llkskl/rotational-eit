%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is part of the code for our articles "Finite Element Mapping for
% Efficient Image Reconstruction in Rotational Electrical Impedance 
% Tomography" and "Rotational electrical impedance tomography using 
% electrodes with limited boundary coverage provides window for multimodal
% sensing"
%
% An example use case of the codes related to above articles is presented
% in this file.
%
% Licenced GPL v. 3
% Olli Koskela, olli.koskela@tut.fi, 2017
% 
% This package is developed using EIDORS v 3.8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Before running this file you should have started EIDORS with startup.m
% provided in your EIDORS distribution folder.

%% Create rotational empty rotational mesh

n_rotpos = 32;
electrode_size = 0.14;
electrode_centres = [-33.75 -11.25 11.25 33.75  146.25  168.75  191.25  213.75];
surface_r = 1;
object_r = 0.9*surface_r;
max_vertex_length = 0.03;

rotationalMesh = CreateRotationalMesh(n_rotpos,electrode_size,electrode_centres,...
    surface_r, object_r, max_vertex_length);

figure('Name','Mesh'),
show_fem(rotationalMesh)

%% Create phantoms

homogPhantom = rotationalMesh;
homogPhantom.elem_data(:) = 1;
homogPhantom.fwd_model.stimulation = EightElectrodesAllCombinations;

% in homogenous phantom: small object and large object
TwoObjectPhantom = homogPhantom;
select_fun = @(x,y,z) (x-0.4).^2+(y-0.25).^2<0.1^2;
memb_ship = elem_select(TwoObjectPhantom.fwd_model, select_fun);
for i = 1:numel(memb_ship)
    if memb_ship(i) > 0
        TwoObjectPhantom.elem_data(i) = 0.1;
    end
end
select_fun = @(x,y,z) (x- (-0.25)).^2+(y-0.25).^2<0.2^2;
memb_ship = elem_select(TwoObjectPhantom.fwd_model, select_fun);
for i = 1:numel(memb_ship)
    if memb_ship(i) > 0
        TwoObjectPhantom.elem_data(i) = 0.1;
    end
end
TwoObjectPhantom.name = 'TOP';

% plot both phantoms
figure('Name','Homogenous phantom');
show_fem(homogPhantom);
figure('Name','Phantom');
show_fem(TwoObjectPhantom);

%% Weighted mapping data, as in MeasSciTech paper

% rotational measurement positions
thetas = linspace(0,360,n_rotpos+1); thetas(end) = [];

tic
% Homogenous data
datahmg_NF = RotationalEITmeasurement(homogPhantom ,thetas);
% Data w/ inclusions
data_NF = RotationalEITmeasurement(TwoObjectPhantom ,thetas);
t1 = toc;

eidors_cache( 'clear_all' )

%% Clicking mapping data, as in IFMBE proceedings

tic
% Homogenous
data_homog_click = ClickingRotationalMeasurement(TwoObjectPhantom,n_rotpos);
% Data w/ inclusions
data_click = ClickingRotationalMeasurement(homogPhantom,n_rotpos);
t2 = toc;

eidors_cache( 'clear_all' )

%% Create inverse mesh

electrode_size = 0.14;
electrode_centres = [-33.75 -11.25 11.25 33.75  146.25  168.75  191.25  213.75];
% tarkista tämä vielä elektrodilla paikassa 1 aste
surface_r = 1;
object_r = 0.9*surface_r;
max_vertex_length = 0.06;

inverseMesh = CreateRotationalMesh(n_rotpos,electrode_size,electrode_centres,...
    surface_r, object_r, max_vertex_length);

figure,
show_fem(inverseMesh)

%
regularization_parameter = 0.03;
inverseMesh.solve = 'rotational_inv_solve_diff_GN_one_step';
inverseMesh.reconst_type = 'difference';
inverseMesh.jacobian_bkgnd.value = 1;
inverseMesh.hyperparameter.value = regularization_parameter; 
inverseMesh.type = 'inv_model';
inverseMesh.parameters.max_iterations = 10;
inverseMesh.RtR_prior = @prior_laplace; priorname = 'laplace';

% cross correlation matrix (set to identiety matrix, no cross correlation)
nmeas = numel(data_homog_click);
inverseMesh.meas_icov = speye(nmeas);
inverseMesh.fwd_model.stimulation = EightElectrodesAllCombinations;

%% Solve using weighted mapping, as in MeasSciTech paper

tic
% create mapping matrix
inverseMesh.thetas = thetas;
M = CreateRotationMatrixCW(inverseMesh,inverseMesh.thetas);
inverseMesh.rotationalmapping = M;
tM = toc;

eidors_cache( 'clear_all' )

% common error level
common_err = 1e-4*max(data_NF(:));
% componentwise error level
component_err = 1e-3;
% noisy datas
data = AddTwoComponentNoise(data_NF,common_err,component_err);
datahmg = AddTwoComponentNoise(datahmg_NF,common_err,component_err);

tic
rec = inv_solve(inverseMesh,datahmg,data);
tRecOld = toc;

eidors_cache( 'clear_all' )

% show reconstruction
figure('Name','Reconstruction');
show_fem(rec);

%% Solve using weighted mapping, as in IFMBE proceedings

inverseMesh2 = inverseMesh;
inverseMesh2.solve = 'rotational2_inv_solve_diff_GN_one_step';


% common error level
common_err = 1e-4*max(data_NF(:));
% componentwise error level
component_err = 1e-3;
% noisy datas
data2 = AddTwoComponentNoise(data_homog_click,common_err,component_err);
datahmg2 = AddTwoComponentNoise(data_click,common_err,component_err);

tic
rec = inv_solve(inverseMesh2,datahmg2,data2);
tRecNew = toc;

eidors_cache( 'clear_all' )
fprintf('Using weighted mapping, 2X forward problem takes %i seconds.\n',t1);
fprintf('Using clicking mapping, 2X forward problem takes %i seconds.\n',t2);
fprintf('Creation of M takes %i seconds.\n',tM);
fprintf('Weighted mapping rec takes %i seconds.\n',tRecOld);
fprintf('Clicking mapping rec takes %i seconds.\n',tRecNew);

% show reconstruction
figure('Name','Reconstruction new');
show_fem(rec);


