classdef getDestructuredResults
    methods(Static)
        function [vir, drug, drugdp, tar, inf, chem, vir_b, drug_b, drugdp_b, tar_b, inf_b] = getResults_LES_Drug(sol, Mesh, dir)
            % no chemokine, 4 compartments: lumen, epi, stroma, blood
            % dir gives the direction of the sol input. 
            % dir = 1: r x t
            % dir = 2: t x r
            if dir == 1
                vir = sol(1 : Mesh.numX, :);
                drug = sol(Mesh.numX + 1 : Mesh.numX + Mesh.numX, :);
                drugdp = sol(2 * Mesh.numX + 1 : 2 * Mesh.numX + Mesh.numT, :);
                tar = sol(2 * Mesh.numX + Mesh.numT + 1 : 2 * Mesh.numX + Mesh.numT + Mesh.numS, :); 
                inf = sol(2 * Mesh.numX + Mesh.numT + Mesh.numS + 1 : 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS, :); 
                chem = [];
                
                vir_b = sol(2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 1, :);
                drug_b = sol(2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 2, :);
                drugdp_b = sol(2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 3, :);
                tar_b = sol(2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 4, :);
                inf_b = sol(2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 5, :);
            elseif dir == 2
                vir = sol(:, 1 : Mesh.numX);
                drug = sol(:, Mesh.numX + 1 : Mesh.numX + Mesh.numX);
                drugdp = sol(:, 2 * Mesh.numX + 1 : 2 * Mesh.numX + Mesh.numT);
                tar = sol(:, 2 * Mesh.numX + Mesh.numT + 1 : 2 * Mesh.numX + Mesh.numT + Mesh.numS); 
                inf = sol(:, 2 * Mesh.numX + Mesh.numT + Mesh.numS + 1 : 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS); 
                chem = [];
                
                vir_b = sol(:, 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 1);
                drug_b = sol(:, 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 2);
                drugdp_b = sol(:, 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 3);
                tar_b = sol(:, 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 4);
                inf_b = sol(:, 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 5);
            end
        end 
        function [vir, drug, drugdp, tar, inf, chem, vir_b, drug_b, drugdp_b, tar_b, inf_b] = getResults_LES(sol, Mesh, dir)
            % no chemokine, 4 compartments: lumen, epi, stroma, blood
            % dir gives the direction of the sol input. 
            % dir = 1: r x t
            % dir = 2: t x r
            if (dir == 1)
                vir = sol(1 : Mesh.numX, :);
                tar = sol(Mesh.numX + 1 : Mesh.numX + Mesh.numS, :); 
                inf = sol(Mesh.numX + Mesh.numS + 1 : Mesh.numX + (2 * Mesh.numS), :); 
                vir_b = sol(Mesh.numX + (2 * Mesh.numS) + 1, :);
                tar_b = sol(Mesh.numX + (2 * Mesh.numS) + 2, :);
                inf_b = sol(Mesh.numX + (2 * Mesh.numS) + 3, :);
            elseif (dir == 2)
                vir = sol(:, 1 : Mesh.numX);
                tar = sol(:, Mesh.numX + 1 : Mesh.numX + Mesh.numS); 
                inf = sol(:, Mesh.numX + Mesh.numS + 1 : Mesh.numX + (2 * Mesh.numS)); 
                vir_b = sol(:, Mesh.numX + (2 * Mesh.numS) + 1);
                tar_b = sol(:, Mesh.numX + (2 * Mesh.numS) + 2);
                inf_b = sol(:, Mesh.numX + (2 * Mesh.numS) + 3);
            end
            chem = [];
            drug = [];
            drugdp = [];
            drug_b = [];
            drugdp_b = [];
            

        end
        function [vir, drug, drugdp, tar, inf, chem, vir_b, drug_b, drugdp_b, tar_b, inf_b] = getResults_LT(sol, Mesh, dir)
            % no chemokine, 3 compartments: lumen, tissue, blood
            % dir gives the direction of the sol input. 
            % dir = 1: r x t
            % dir = 2: t x r
            if dir == 1
                vir = sol(1 : Mesh.numX, :);
                tar = sol(Mesh.numX + 1 : Mesh.numX + Mesh.numT, :); 
                inf = sol(Mesh.numX + Mesh.numT + 1 : Mesh.numX + (2 * Mesh.numT), :); 
                vir_b = sol(Mesh.numX + (2 * Mesh.numT) + 1, :);
                tar_b = sol(Mesh.numX + (2 * Mesh.numT) + 2, :);
                inf_b = sol(Mesh.numX + (2 * Mesh.numT) + 3, :);
            elseif dir == 2
                vir = sol(:, 1 : Mesh.numX);
                tar = sol(:, Mesh.numX + 1 : Mesh.numX + Mesh.numT); 
                inf = sol(:, Mesh.numX + Mesh.numT + 1 : Mesh.numX + (2 * Mesh.numT)); 
                vir_b = sol(:, Mesh.numX + (2 * Mesh.numT) + 1);
                tar_b = sol(:, Mesh.numX + (2 * Mesh.numT) + 2);
                inf_b = sol(:, Mesh.numX + (2 * Mesh.numT) + 3);
            end
            chem = [];
            drug = [];
            drugdp = [];
            drug_b = [];
            drugdp_b = [];
           
        end
    end
end

