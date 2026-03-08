clear; clc;
close all;

% Model define
cell_size = 3;
bond_dir{1} = [1; 0]; % direction of interaction bond
bond_dir{2} = [-1/2; -sqrt(3)/2];
bond_dir{3} = [-1/2; sqrt(3)/2];

Jxy = 1;
Jz = 1;
h = [0; 0; 0];

Jmat = diag([Jxy, Jxy, Jz]);

params = Model(cell_size);
params.add_two_site_intr(1, 2, Jmat, bond_dir{1});
params.add_two_site_intr(1, 2, Jmat, bond_dir{2});
params.add_two_site_intr(1, 2, Jmat, bond_dir{3});

params.add_two_site_intr(2, 3, Jmat, bond_dir{1});
params.add_two_site_intr(2, 3, Jmat, bond_dir{2});
params.add_two_site_intr(2, 3, Jmat, bond_dir{3});

params.add_two_site_intr(3, 1, Jmat, bond_dir{1});
params.add_two_site_intr(3, 1, Jmat, bond_dir{2});
params.add_two_site_intr(3, 1, Jmat, bond_dir{3});

params.add_magnetic_field(h);



% Optimize classical energy
alpha = Jxy/Jz;
ctheta = 1/(1 + alpha);
stheta = sqrt(1-ctheta^2);
theta = acos(ctheta);

% conf{1} = [stheta; 0; ctheta];
% conf{2} = [-stheta; 0; ctheta];
% conf{3} = [0; 0; -1];

% [E, conf] = params.opt_energy([theta, -theta, pi, 0, 0, 0]);
[E, conf] = params.opt_energy([pi/2, pi/2, pi/2, 2*pi/3, -2*pi/3, 0]);

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



Hk = params.get_Hk(K_point);
omega = params.spin_wave_spec(path);


figure;
plot(1:size(path, 1), omega', 'LineWidth', 3);
hold on;

xticks([1, nk, 2*nk-1, size(path, 1)]);
xticklabels({'\Gamma', 'K', 'M', '\Gamma'});
ylabel('\omega');
xlim([1, size(path, 1)]);
set(gca, 'LineWidth', 2, 'FontSize', 20);

