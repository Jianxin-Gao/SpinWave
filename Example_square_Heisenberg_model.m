clear; clc;
close all;

% Model define
cell_size = 2;
bond_dir{1} = [1; 0; 0]; % direction of interaction bond
bond_dir{2} = [0; 1; 0];

J = 1;
h = [0; 0; 0];

Jmat = diag([J, J, J]);

sys = Model(cell_size);

% sys.pos{1} = [0; 0; 0];
% sys.pos{2} = [1; 0; 0];

sys.add_two_site_intr(1, 2, Jmat, bond_dir{1});
sys.add_two_site_intr(1, 2, Jmat, bond_dir{2});
sys.add_two_site_intr(1, 2, Jmat, -bond_dir{1});
sys.add_two_site_intr(1, 2, Jmat, -bond_dir{2});

sys.add_magnetic_field(h);



% Optimize classical energy
sys.conf{1} = [0; 0; 1];
sys.conf{2} = [0; 0; -1];

sys.plot_spin_configuration();




% Spin wave spec
nk = 100;
Gamma = [0, 0, 0];
X_point = [pi, 0, 0];
M_point = [pi, pi, 0];

path = [
    [linspace(Gamma(1), X_point(1), nk)', linspace(Gamma(2), X_point(2), nk)', linspace(Gamma(3), X_point(3), nk)'];
    [linspace(X_point(1), M_point(1), nk)', linspace(X_point(2), M_point(2), nk)', linspace(Gamma(3), X_point(3), nk)'];
    [linspace(M_point(1), Gamma(1), nk)', linspace(M_point(2), Gamma(2), nk)', linspace(Gamma(3), X_point(3), nk)']
    ];
path = unique(path, 'rows', 'stable');
kx = path(:,1)'; ky = path(:,2)';



omega_exact = 4 * J * sys.S * sqrt(1-1/4*(cos(kx)+cos(ky)).^2); % exact solution

% Plot
sys.spin_wave_spec(path, 'dynamical', true);
sys.plot_spin_wave_spec(); hold on;
plot(1:size(path, 1), omega_exact', 'o', 'LineWidth', 2);
xticks([1, nk, 2*nk-1, size(path, 1)]);
xticklabels({'\Gamma', 'X', 'M', '\Gamma'});


Sk_tot = sys.struc_fac_utils('total');
sys.plot_dynamics(Sk_tot);
xticks([1, nk, 2*nk-1, size(path, 1)]);
xticklabels({'\Gamma', 'X', 'M', '\Gamma'});


