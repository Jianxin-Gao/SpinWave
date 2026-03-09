function [omega, intensity] = spin_wave_spec(obj, path)
N = obj.cell_size;
omega = zeros(N, size(path,1));
intensity = zeros(N, size(path,1));

g_mat = diag([ones(1, N), -ones(1, N)]);
Rmat = obj.get_rotation_matrix_();

P_mat = zeros(3, 2 * N);
for j = 1:N
    R = Rmat{j};
    u_j = sqrt(obj.S / 2) * (R(:, 1) - 1i * R(:, 2));
    v_j = sqrt(obj.S / 2) * (R(:, 1) + 1i * R(:, 2));
    P_mat(:, j) = u_j;
    P_mat(:, j + N) = v_j;
end

for itk = 1:size(path, 1)
    k = path(itk, :);
    q = k(:);
    
    Hk = obj.get_Hk(k);
    
    [V, D] = eig(g_mat * Hk);
    eig_vals = diag(D);
    if max(abs(imag(eig_vals))) > 1e-10
        warning('Imagary part is larger than 1e-10')
    end
    phys_modes = zeros(1, N);
    T_phys = zeros(2 * N, N);
    mode_idx = 1;
    
    for m = 1:(2*N)
        norm_val = real(V(:, m)' * g_mat * V(:, m));
        if norm_val > 1e-6 && mode_idx <= N
            phys_modes(mode_idx) = real(eig_vals(m));
            T_phys(:, mode_idx) = V(:, m) / sqrt(norm_val);
            mode_idx = mode_idx + 1;
        end
    end
    
    [omega(:, itk), sort_idx] = sort(phys_modes);
    T_phys = T_phys(:, sort_idx);
    
    Lambda = eye(2*N);
    for j = 1:N
        if ~isempty(obj.pos{j})
            phase = exp(-1i * dot(q, obj.pos{j}));
            Lambda(j, j) = phase;
            Lambda(j + N, j + N) = phase;
        end
    end
    
    V_mat = P_mat * Lambda * T_phys;
    
    q_norm = norm(q);
    if q_norm < 1e-8
        Pol = eye(3);
    else
        q_dir = q / q_norm;
        Pol = eye(3) - q_dir * q_dir'; % delta_ab - q_a*q_b
    end
    
    for n = 1:N
        intensity(n, itk) = real(V_mat(:, n)' * Pol * V_mat(:, n));
    end
end

obj.omega = omega;
obj.intensity = intensity;
obj.path = path;
end