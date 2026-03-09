clear; clc;
close all;

% Model define
cell_size = 2;
bond_dir{1} = [1; 0; 0]; % direction of interaction bond
bond_dir{2} = [0; 1; 0];

J = 1;
h = [0; 0; 0];

Jmat = diag([J, J, J]);

params = Model(cell_size);

params.pos{1} = [0; 0; 0];
params.pos{2} = [1; 0; 0];

params.add_two_site_intr(1, 2, Jmat, bond_dir{1});
params.add_two_site_intr(1, 2, Jmat, bond_dir{2});
params.add_two_site_intr(1, 2, Jmat, -bond_dir{1});
params.add_two_site_intr(1, 2, Jmat, -bond_dir{2});

params.add_magnetic_field(h);



% Optimize classical energy
params.conf{1} = [0; 0; 1];
params.conf{2} = [0; 0; -1];

params.plot_spin_configuration();




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

[omega, intensity] = params.spin_wave_spec(path);


omega_exact = 4 * J * params.S * sqrt(1-1/4*(cos(kx)+cos(ky)).^2); % exact solution

% Plot
figure;
params.plot_spin_wave_spec('sigma', 0.05); hold on;
plot(1:size(path, 1), omega_exact', 'r--', 'LineWidth', 2);
xticks([1, nk, 2*nk-1, size(path, 1)]);
xticklabels({'\Gamma', 'X', 'M', '\Gamma'});
cbar = colorbar;
cbar.Label.String = 'S(q, \omega)';
cbar.Label.FontSize = 20;
xticks([1, nk, 2*nk-1, size(path, 1)]);
xticklabels({'\Gamma', 'X', 'M', '\Gamma'});
ylabel('\omega');
xlim([1, size(path, 1)]);
% ylim([0, max(omega(1, :))*1.3])
set(gca, 'LineWidth', 2, 'FontSize', 20);

