classdef Model < handle
    properties
        cell_size
        intr_mat
        conf
        S = 1/2
    end
    
    methods
        function obj = Model(cell_size)
            obj.cell_size = cell_size;
            obj.intr_mat = struct('two_site', struct('site1', [], 'site2', [], 'Jmat', [], 'd', []), ...
                'magnetic_field', [0;0;0]);
        end
        
        function add_two_site_intr(obj, site1, site2, Jmat, d)
            % two-site spin-spin interaction
            cnt = length([obj.intr_mat.two_site.site1]);
            obj.intr_mat.two_site(cnt + 1).site1 = site1;
            obj.intr_mat.two_site(cnt + 1).site2 = site2;
            obj.intr_mat.two_site(cnt + 1).Jmat = Jmat;
            obj.intr_mat.two_site(cnt + 1).d = d;
            
            obj.intr_mat.two_site(cnt + 2).site1 = site2;
            obj.intr_mat.two_site(cnt + 2).site2 = site1;
            obj.intr_mat.two_site(cnt + 2).Jmat = Jmat.';
            obj.intr_mat.two_site(cnt + 2).d = -d;
        end
        
        function add_magnetic_field(obj, h)
            obj.intr_mat.magnetic_field = h;
        end
        
        function E = energy_func(obj, params)
            % E = 1/2 * \sum_{ijmn}{S_{im}.' * J_{ijmn} * S_{jn}} - \sum_{im}hS{im}
            lstheta = params(1:obj.cell_size); lsphi = params(obj.cell_size+1:end);
            h = obj.intr_mat.magnetic_field;
            
            E = 0;
            conf_input = cell(1, obj.cell_size);
            
            % Zeeman term
            for it = 1:obj.cell_size
                conf_input{it} = [sin(lstheta(it))*cos(lsphi(it)); sin(lstheta(it))*sin(lsphi(it)); cos(lstheta(it))];
                E = E - dot(h, conf_input{it});
            end
            
            % two-site term
            for it = 1:length(obj.intr_mat.two_site)
                site1 = obj.intr_mat.two_site(it).site1;
                site2 = obj.intr_mat.two_site(it).site2;
                J = obj.intr_mat.two_site(it).Jmat;
                E = E + 1/2 * conf_input{site1}.' * J * conf_input{site2};
            end
            
        end
        function [E, conf] = opt_energy(obj, params0)
            conf = cell(1, obj.cell_size);
            options = optimset('Display','iter', ...
                'MaxFunEvals', 60000, ...
                'MaxIter', 20000, ...
                'PlotFcns', @optimplotfval, ...
                'TolFun', 1e-32, ...
                'TolX', 1e-32);
            [params, E] = fminsearch(@(params) energy_func(obj, params), params0, options);
            lstheta = params(1:obj.cell_size);
            lsphi = params(obj.cell_size+1:end);
            
            for it = 1:obj.cell_size
                conf{it} = [sin(lstheta(it))*cos(lsphi(it)); sin(lstheta(it))*sin(lsphi(it)); cos(lstheta(it))];
            end
        end
        
        function plot_spin_configuration(obj, varargin)
            if ~isempty(varargin)
                conf_input = varargin{1};
            else
                conf_input = obj.conf;
            end
            
            figure; hold on;
            for it = 1:obj.cell_size
                conf_input{it} = conf_input{it}/norm(conf_input{it});
                quiver3(0,0,0, conf_input{it}(1), conf_input{it}(2), conf_input{it}(3), ...
                    'Color', [0,0,0], ...
                    'LineWidth', 2, ...
                    'MaxHeadSize', 0.2);
            end
            xlabel('x')
            ylabel('y')
            zlabel('z')
            view(-37.5, 30)
            axis equal;
        end
        
    end
    
    
    methods (Access = private)
        function Rmat = get_rotation_matrix_(obj)
            Rmat = cell(1, length(obj.conf));
            for it = 1:length(Rmat)
                loc_S = obj.conf{it} / norm(obj.conf{it});
                theta = acos(loc_S(3));
                phi = atan2(loc_S(2), loc_S(1));
                Rmat{it} = [
                    cos(theta)*cos(phi), -sin(phi), sin(theta)*cos(phi);
                    cos(theta)*sin(phi),  cos(phi), sin(theta)*sin(phi);
                    -sin(theta),           0,        cos(theta)
                    ];
            end
        end
    end
end
