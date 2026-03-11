clear; clc;
close all;

% Model define
cell_size = 2;
bond_dir{1} = [1; 0; 0]; % direction of interaction bond
bond_dir{2} = [-1; 0; 0];

J = 1;
h = [0; 0; 0];

Jmat = J*eye(3);

sys = Model(cell_size);

% sys.pos{1} = [-1/2; 0; 0];
% sys.pos{2} = [1/2; 0; 0];
% sys.pos{3} = [0; sqrt(3)/2; 0];


sys.add_two_site_intr(1, 2, Jmat, bond_dir{1});
sys.add_two_site_intr(1, 2, Jmat, bond_dir{2});

sys.add_two_site_intr(2, 1, Jmat, bond_dir{1});
sys.add_two_site_intr(2, 1, Jmat, bond_dir{2});

sys.add_magnetic_field(h);

% AFM ground state
conf{1} = [0; 0; 1];
conf{2} = [0; 0; -1];

sys.conf = conf;


% Spin wave spec
nk = 101;


path = [linspace(-pi, pi, nk)', zeros(nk, 2)];

kx = path(:,1)';
spec_exact = 4 * J * sys.S * abs(sin(kx));

figure;
sys.spin_wave_spec(path);
sys.plot_spin_wave_spec(); hold on;
plot(1:nk, spec_exact, 'o')
xticks([1, (nk+1)/2, nk]);
xticklabels({'-\pi', '0', '\pi'});