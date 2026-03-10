function spin_wave_spec(obj, path, varargin)
% Local change obj.omega and obj.dynamic_struc
% omega (cell_size \times n_k): linear spin wave spectrum
% dynamic_struc (3 \times 3 \times n_k \times n_omega): S^{\alpha \beta}(k, omega)

p = inputParser;
addParameter(p, 'sigma', 0.01); % Gauss sigma factor
addParameter(p, 'n_omega', 200); % Number of omega in dynamical calculation
addParameter(p, 'dynamical', false); % if calculate dynamical structure factor or not
parse(p, varargin{:});

N = obj.cell_size;
obj.omega = zeros(N, size(path,1));
obj.path = path;

g_mat = diag([ones(1, N), -ones(1, N)]);
V_sorted = cell(1, size(path, 1));

for itk = 1:size(path, 1)
    k = path(itk, :);
    Hk = obj.get_Hk(k);
    [V, D] = eig(g_mat * Hk);
    eig_vals = diag(D);
    evals = real(eig_vals);
    [evals_sorted, idx] = sort(evals, 'descend');
    
    obj.omega(:, itk) = evals_sorted(1:N);
    V_sorted{itk} = V(:,idx);
    if max(abs(imag(eig_vals))) > 1e-10
        warning('Imagary part of eig(g*H_k) is larger than 1e-10!')
    end
    
    
end


if p.Results.dynamical
    
    Rmat = obj.get_rotation_matrix_();
    
    P_mat = zeros(3, 2 * N);
    omega_max = max(max(obj.omega)) * 1.1;
    num_omega = p.Results.n_omega;
    omega_grid = linspace(0, omega_max, num_omega); 
    sigma = p.Results.sigma;
    obj.dynamic_struc = zeros(3, 3, size(path, 1), num_omega);
    obj.omega_list = omega_grid;
    for j = 1:N
        R = Rmat{j};
        u_j = sqrt(obj.S / 2) * (R(:, 1) - 1i * R(:, 2));
        v_j = sqrt(obj.S / 2) * (R(:, 1) + 1i * R(:, 2));
        P_mat(:, j) = u_j;
        P_mat(:, j + N) = v_j;
    end
    
    for itk = 1:size(path, 1)
        omegas = obj.omega(:, itk);
        
        
        Tq = zeros(2*N,2*N);
        for i = 1:2*N
            g_norm = V_sorted{itk}(:,i)' * g_mat * V_sorted{itk}(:,i);
            Tq(:,i) = V_sorted{itk}(:,i)/sqrt(abs(g_norm));
        end
        W = P_mat*Tq;
        
        I_tensor = zeros(3, 3, N);
        for nu = 1: N
            I_tensor(:,:,nu) = W(:,nu) * W(:,nu)';
        end
        S_ab = zeros(3, 3, num_omega);
        for nu = 1:N
            dw = omega_grid - omegas(nu);
            broadening_profile = (1 / (sigma * sqrt(2*pi))) * exp(-0.5 * (dw / sigma).^2);
            for alpha = 1:3
                for beta = 1:3
                    intensity = I_tensor(alpha, beta, nu) * broadening_profile;
                    S_ab(alpha, beta, :) = S_ab(alpha, beta, :) + reshape(intensity, 1, 1, num_omega);
                end
            end
            
        end
        
        obj.dynamic_struc(:, :, itk, :) = S_ab;
    end
end




end