function dynamic_struc_fac = struc_fac_utils(obj, type)
path = obj.path;
nk = size(path, 1); nomega = length(obj.omega_list);

switch type
    case 'total'
        % S(q, omega) = trace( S^{\alpha, \beta}(q, omega) )
        dynamic_struc_fac = obj.dynamic_struc(1,1,:,:) + obj.dynamic_struc(2,2,:,:) + obj.dynamic_struc(3,3,:,:);
        dynamic_struc_fac = permute(dynamic_struc_fac, [3, 4, 1, 2]);
    case 'INS'
        % S(q, omega) = \sum_{\alpha, \beta} (\delta_{\alpha\beta} - q_\alpha q_\beta/|q|^2 S^{\alpha, \beta}(q, omega))
        dynamic_struc_fac = zeros(nk, nomega);
        for itk = 1:nk
            k_vec = obj.path(itk, :).';
            if norm(k_vec) > 1e-10
                k_vec = k_vec/norm(k_vec);
            else
                k_vec = [0;0;0];
            end
            D = eye(3) - k_vec * k_vec';
            Sab_k = squeeze(obj.dynamic_struc(:, :, itk, :));
            dynamic_struc_fac(itk, :) = squeeze(sum(sum(D .* Sab_k, 1), 2));
            
        end
    case 'para-z'
        % S(q, omega) = S^{zz}(q, omega)
        dynamic_struc_fac = squeeze(obj.dynamic_struc(3, 3, :, :));
    case 'perp-z'
        % S(q, omega) = S^{xx}(q, omega) + S^{yy}(q, omega)
        dynamic_struc_fac = squeeze(obj.dynamic_struc(1, 1, :, :) + obj.dynamic_struc(2, 2, :, :));
        
end