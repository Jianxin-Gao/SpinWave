function omega = spin_wave_spec(obj, path)
omega = zeros(obj.cell_size, size(path,1));
g_mat = diag([ones(1, obj.cell_size), -ones(1, obj.cell_size)]);

for itk = 1:size(path, 1)
    k = path(itk, :);
    Hk = obj.get_Hk(k);
    spec = eig(g_mat * Hk);
    real_spec = sort(real(spec));
    omega(:, itk) = real_spec(obj.cell_size+1:end);
end
end

