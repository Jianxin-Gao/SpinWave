clear; clc;
close all;

% Model define
cell_size = 3;
bond_dir{1} = [1; 0; 0]; % direction of interaction bond
bond_dir{2} = [-1/2; -sqrt(3)/2; 0];
bond_dir{3} = [-1/2; sqrt(3)/2; 0];

Jxy = 1;
Jz = 1;
h = [0; 0; 0];

Jmat = diag([Jxy, Jxy, Jz]);

sys = Model(cell_size);

% sys.pos{1} = [-1/2; 0; 0];
% sys.pos{2} = [1/2; 0; 0];
% sys.pos{3} = [0; sqrt(3)/2; 0];


sys.add_two_site_intr(1, 2, Jmat, bond_dir{1});
sys.add_two_site_intr(1, 2, Jmat, bond_dir{2});
sys.add_two_site_intr(1, 2, Jmat, bond_dir{3});

sys.add_two_site_intr(2, 3, Jmat, bond_dir{1});
sys.add_two_site_intr(2, 3, Jmat, bond_dir{2});
sys.add_two_site_intr(2, 3, Jmat, bond_dir{3});

sys.add_two_site_intr(3, 1, Jmat, bond_dir{1});
sys.add_two_site_intr(3, 1, Jmat, bond_dir{2});
sys.add_two_site_intr(3, 1, Jmat, bond_dir{3});

sys.add_magnetic_field(h);



% Optimize classical energy
% % % Y-state ground state
% alpha = Jxy/Jz;
% ctheta = 1/(1 + alpha);
% stheta = sqrt(1-ctheta^2);
% theta = acos(ctheta);
% conf{1} = [stheta; 0; ctheta];
% conf{2} = [-stheta; 0; ctheta];
% conf{3} = [0; 0; -1];
% 120-order ground state
theta = 2*pi/3;
conf{1} = [sin(theta); -cos(theta); 0];
conf{2} = [-sin(theta); -cos(theta); 0];
conf{3} = [0; -1; 0];

% You can choose to optimize configuration as follows:
% [E, conf] = params.opt_energy([theta, -theta, pi, 0, 0, 0]);

sys.conf = conf;
sys.plot_spin_configuration();


% Spin wave spec
nk = 100;
Gamma = [0, 0, 0];
M_point = [pi, pi/sqrt(3), 0];
K_point = [4*pi/3, 0, 0];

path = [
    [linspace(Gamma(1), K_point(1), nk)', linspace(Gamma(2), K_point(2), nk)', linspace(Gamma(3), K_point(3), nk)'];
    [linspace(K_point(1), M_point(1), nk)', linspace(K_point(2), M_point(2), nk)', linspace(K_point(3), M_point(3), nk)'];
    [linspace(M_point(1), Gamma(1), nk)', linspace(M_point(2), Gamma(2), nk)', linspace(M_point(3), Gamma(3), nk)']
    ];
path = unique(path, 'rows', 'stable');
kx = path(:,1)'; ky = path(:,2)';


sys.spin_wave_spec(path, 'dynamical', true, 'dynamical_broad', true, 'sigma', 0.025, 'n_omega', 3000);

figure;
sys.plot_spin_wave_spec();

xticks([1, nk, 2*nk-1, size(path, 1)]);
xticklabels({'\Gamma', 'K', 'M', '\Gamma'});


Sk_tot = sys.struc_fac_utils('INS');

figure;
sys.plot_dynamics(Sk_tot);
xticks([1, nk, 2*nk-1, size(path, 1)]);
xticklabels({'\Gamma', 'K', 'M', '\Gamma'});
% caxis([0,1])
% figure; hold on;
% params.plot_spin_wave_spec('sigma', 0.03);
% caxis([0, 2.5]);
% cbar = colorbar;
% cbar.Label.String = 'S(q, \omega)';
% cbar.Label.FontSize = 20;
% xticks([1, nk, 2*nk-1, size(path, 1)]);
% xticklabels({'\Gamma', 'K', 'M', '\Gamma'});
% ylabel('\omega');
% xlim([1, size(path, 1)]);
% set(gca, 'LineWidth', 2, 'FontSize', 20);

