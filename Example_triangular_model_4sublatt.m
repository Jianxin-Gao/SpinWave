clear; clc;
close all;

% Model define
cell_size = 4;
bond_dir{1} = [1; 0]; % direction of interaction bond
bond_dir{2} = [-1/2; -sqrt(3)/2];
bond_dir{3} = [-1/2; sqrt(3)/2];

Jxy = 0.1;
Jz = 1;
Jpd = 0.1;
JGamma = 0.3;
h = [0; 0; 2.7];

params = Model(cell_size);
params.add_two_site_intr(1, 3, Jmat(Jxy, Jz, Jpd, JGamma, 0), bond_dir{1});
params.add_two_site_intr(1, 2, Jmat(Jxy, Jz, Jpd, JGamma, 2*pi/3), bond_dir{3});
params.add_two_site_intr(1, 4, Jmat(Jxy, Jz, Jpd, JGamma, -2*pi/3), bond_dir{2});

params.add_two_site_intr(2, 3, Jmat(Jxy, Jz, Jpd, JGamma, -2*pi/3), bond_dir{2});
params.add_two_site_intr(2, 4, Jmat(Jxy, Jz, Jpd, JGamma, 0), bond_dir{1});
params.add_two_site_intr(2, 1, Jmat(Jxy, Jz, Jpd, JGamma, 2*pi/3), bond_dir{3});

params.add_two_site_intr(3, 4, Jmat(Jxy, Jz, Jpd, JGamma, 2*pi/3), bond_dir{3});
params.add_two_site_intr(3, 2, Jmat(Jxy, Jz, Jpd, JGamma, -2*pi/3), bond_dir{2});
params.add_two_site_intr(3, 1, Jmat(Jxy, Jz, Jpd, JGamma, 0), bond_dir{1});

params.add_two_site_intr(4, 1, Jmat(Jxy, Jz, Jpd, JGamma, -2*pi/3), bond_dir{2});
params.add_two_site_intr(4, 2, Jmat(Jxy, Jz, Jpd, JGamma, 0), bond_dir{1});
params.add_two_site_intr(4, 3, Jmat(Jxy, Jz, Jpd, JGamma, 2*pi/3), bond_dir{3});


params.add_magnetic_field(h);

% Optimize classical energy
[E, conf] = params.opt_energy([rand(1,cell_size)*pi, rand(1,cell_size)*2*pi]);
params.conf = conf;
params.plot_spin_configuration();

% Spin wave spec
nk = 100;
Gamma = [0, 0];
M_point = [pi, pi/sqrt(3)];
K_point = [4*pi/3, 0];
path = [
    [linspace(Gamma(1), K_point(1), nk)', linspace(Gamma(2), K_point(2), nk)'];[linspace(K_point(1), M_point(1), nk)', linspace(K_point(2), M_point(2), nk)'];[linspace(M_point(1), Gamma(1), nk)', linspace(M_point(2), Gamma(2), nk)']
    ];
path = unique(path, 'rows', 'stable');
kx = path(:,1)'; ky = path(:,2)';

omega = params.spin_wave_spec(path);

figure;
plot(1:size(path, 1), omega', 'LineWidth', 3);
hold on;

xticks([1, nk, 2*nk-1, size(path, 1)]);
xticklabels({'\Gamma', 'K', 'M', '\Gamma'});
ylabel('\omega');
xlim([1, size(path, 1)]);
ylim([0, max(omega(1,:))*1.1])
set(gca, 'LineWidth', 2, 'FontSize', 20);



function J = Jmat(Jxy, Jz, Jpd, JGamma, phi)
J =[
    Jxy + 2*Jpd*cos(phi), -2*Jpd*sin(phi), -JGamma*sin(phi);
    -2*Jpd*sin(phi), Jxy - 2*Jpd*cos(phi), JGamma*cos(phi);
    -JGamma*sin(phi), JGamma*cos(phi), Jz
    ];
end

