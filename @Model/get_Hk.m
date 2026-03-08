function [Hk] = get_Hk(obj, k)




% A*(a a) + B*(ad a) + C*(ad ad) + D*(a ad)
Ak = zeros(obj.cell_size);
Bk = zeros(obj.cell_size);
Ck = zeros(obj.cell_size);
Dk = zeros(obj.cell_size);
Rmat = get_rotation_matrix_(obj);
h = obj.intr_mat.magnetic_field;

% local term: cons * a_{im}^\dagger * a_{im}
for site0 = 1:obj.cell_size
    cons = h(1) * Rmat{site0}(1, 3) + h(2) * Rmat{site0}(2, 3) + h(3) * Rmat{site0}(3, 3); % Zeeman term
    for it = 1:length(obj.intr_mat.two_site)
        if obj.intr_mat.two_site(it).site1 == site0 % J_{mnij} S_{im}^z S_{jn}^z term
            site2 = obj.intr_mat.two_site(it).site2;
            J = Rmat{site0}' * obj.intr_mat.two_site(it).Jmat * Rmat{site2};
            cons = cons - obj.S * J(3, 3);
        end
    end
    Bk(site0, site0) = Bk(site0, site0) + cons;
    Dk(site0, site0) = Dk(site0, site0) + cons;
end

for it = 1:length(obj.intr_mat.two_site)
    % two-site interaction
    % H_{ijmn} = A_{ijmn}a_{im}a_{jn}
    %          + B_{ijmn}a_{im}^\dagger a_{jn}
    %          + C_{ijmn}a_{im}^\dagger a_{jn}^\dagger
    %          + D_{ijmn}a_{im} a_{jn}^\dagger
    
    
    site1 = obj.intr_mat.two_site(it).site1;
    site2 = obj.intr_mat.two_site(it).site2;
    J = Rmat{site1}' * obj.intr_mat.two_site(it).Jmat * Rmat{site2};
    A = obj.S/2 * (J(1, 1) - J(2, 2) - 1i*J(1, 2) - 1i*J(2, 1));
    B = obj.S/2 * (J(1, 1) + J(2, 2) - 1i*J(1, 2) + 1i*J(2, 1));
    C = obj.S/2 * (J(1, 1) - J(2, 2) + 1i*J(1, 2) + 1i*J(2, 1));
    D = obj.S/2 * (J(1, 1) + J(2, 2) + 1i*J(1, 2) - 1i*J(2, 1));
    
    d = obj.intr_mat.two_site(it).d;
    
    Ak(site2, site1) = Ak(site2, site1) + A * exp(-1i * dot(k, d));
    Bk(site1, site2) = Bk(site1, site2) + B * exp(1i * dot(k, d));
    Ck(site1, site2) = Ck(site1, site2) + C * exp(1i * dot(k, d));
    Dk(site1, site2) = Dk(site1, site2) + D * exp(1i * dot(k, d));
    
    
end
Hk = [Bk, Ck; Ak, Dk];

end