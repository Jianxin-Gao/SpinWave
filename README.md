# SpinWave: Linear Spin Wave Theory Solver in MATLAB 🌊🧭

![Language](https://img.shields.io/badge/Language-MATLAB-blue.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)

A lightweight, object-oriented MATLAB toolkit for calculating the spin-wave excitation spectrum (magnon dispersion) of magnetic materials based on **Linear Spin Wave Theory (LSWT)**. 

## ✨ Features

- **Object-Oriented Design**: Built around the `@Model` class for clean, intuitive, and modular model construction.
- **Arbitrary Magnetic Structures**: Naturally supports complex multi-sublattice ($N$-sublattice) and non-collinear magnetic ground states. Local reference frames are calculated automatically.
- **Universal Interactions**: Define generic $3 \times 3$ interaction matrices to easily implement Heisenberg exchange, XXZ anisotropy, Dzyaloshinskii-Moriya Interactions (DMI), and external magnetic fields.
- **Classical Energy Optimization**: Built-in numerical optimizer (`opt_energy`) to find the exact classical magnetic ground state configuration.
- **Built-in Visualization**: Quickly plot classical spin configurations to verify your lattice and magnetic unit cell (`plot_spin_configuration`).
- **Bosonic Diagonalization**: Accurately diagonalizes the bosonic Hamiltonian $\mathcal{H}(\mathbf{k})$ using the Colpa/Metric matrix approach to strictly preserve bosonic commutation relations.

## 📂 Project Structure

```text
.
├── @Model/                                  # MATLAB Class folder for the SpinWave Model
│   ├── Model.m                              # Core class definition, properties, and basic methods
│   ├── get_Hk.m                             # Constructs the 2N x 2N Hamiltonian matrix H(k)
│   └── spin_wave_spec.m                     # Diagonalizes H(k) and returns the magnon dispersion
├── Example_triangular_XXZ_model.m           # Example: Y-state/Umbrella state on a triangular XXZ model
├── Example_triangular_model_4sublatt.m      # Example: Complex 4-sublattice order on a triangular lattice
├── LICENSE                                  # MIT License
└── README.md                                # Project documentation
```

## 🚀 Quick Start: Triangular Lattice XXZ Model

Below is a complete example of calculating the spin wave dispersion of a 3-sublattice triangular antiferromagnet with XXZ anisotropy.

```matlab
% 1. Initialize the Model (3 sublattices)
cell_size = 3;
params = Model(cell_size);

% 2. Define Interaction Matrix and Bond Directions
Jxy = 0.5; Jz = 1;
Jmat = diag([Jxy, Jxy, Jz]);

% Triangular lattice nearest-neighbor vectors
bond_dir = {[1; 0], [-1/2; -sqrt(3)/2],[-1/2; sqrt(3)/2]};

% 3. Add the 9 unique physical bonds (3 directions * 3 sublattice pairs)
% The class automatically handles the Hermitian conjugate bonds (-d, Jmat.')
for i = 1:3
    params.add_two_site_intr(1, 2, Jmat, bond_dir{i});
    params.add_two_site_intr(2, 3, Jmat, bond_dir{i});
    params.add_two_site_intr(3, 1, Jmat, bond_dir{i});
end

% Add external magnetic field
params.add_magnetic_field([0; 0; 0]);

% 4. Define Classical Ground State (Analytical Umbrella State)
alpha = Jxy/Jz;
ctheta = 1/(1 + alpha);
stheta = sqrt(1-ctheta^2);

conf{1} =[stheta; 0; ctheta];
conf{2} = [-stheta; 0; ctheta];
conf{3} = [0; 0; -1];

params.conf = conf;
params.plot_spin_configuration(); % Verify configuration visually

% Note: You can also use the built-in optimizer:
% [E, conf] = params.opt_energy([theta, -theta, pi, 0, 0, 0]);

% 5. Define k-path (Gamma -> K -> M -> Gamma)
nk = 100;
Gamma = [0, 0];
M_point = [pi, pi/sqrt(3)];
K_point = [4*pi/3, 0];

path =[
    linspace(Gamma(1), K_point(1), nk)', linspace(Gamma(2), K_point(2), nk)';
    linspace(K_point(1), M_point(1), nk)', linspace(K_point(2), M_point(2), nk)';
    linspace(M_point(1), Gamma(1), nk)', linspace(M_point(2), Gamma(2), nk)'
];
path = unique(path, 'rows', 'stable'); % Remove duplicate connection points

% 6. Calculate Spectrum and Plot
omega = params.spin_wave_spec(path);

figure;
plot(1:size(path, 1), omega', 'LineWidth', 3);
xticks([1, nk, 2*nk-1, size(path, 1)]);
xticklabels({'\Gamma', 'K', 'M', '\Gamma'});
ylabel('Energy \omega');
xlim([1, size(path, 1)]);
set(gca, 'LineWidth', 2, 'FontSize', 20);
```

## 🗺️ Roadmap (Future Works)

- [ ] Support calculation of the **Dynamic Structure Factor $S(\mathbf{q}, \omega)$** to simulate Inelastic Neutron Scattering (INS) cross-sections.
- [ ] Add an automated neighbor-finding algorithm based on atomic spatial distances (`add_interactions_by_distance`).
- [ ] Extend to 3D lattices (e.g., Pyrochlore, Kagome stacks).

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

