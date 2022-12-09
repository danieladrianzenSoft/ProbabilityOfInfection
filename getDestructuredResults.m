classdef getDestructuredResults
    methods(Static)
        function [vir, tar, inf, chem, vir_b, tar_b, inf_b] = getResults_LES(sol, Mesh)
            % no chemokine, 4 compartments: lumen, epi, stroma, blood
            vir = sol(1 : Mesh.numX, :);
            tar = sol(Mesh.numX + 1 : Mesh.numX + Mesh.numS, :); 
            inf = sol(Mesh.numX + Mesh.numS + 1 : Mesh.numX + (2 * Mesh.numS), :); 
            chem = 0;
            
            vir_b = sol(Mesh.numX + (2 * Mesh.numS) + 1, :);
            tar_b = sol(Mesh.numX + (2 * Mesh.numS) + 2, :);
            inf_b = sol(Mesh.numX + (2 * Mesh.numS) + 3, :);
        end
        function [vir, tar, inf, chem, vir_b, tar_b, inf_b] = getResults_LT(sol, Mesh)
            % no chemokine, 3 compartments: lumen, tissue, blood
            vir = sol(1 : Mesh.numX, :);
            tar = sol(Mesh.numX + 1 : Mesh.numX + Mesh.numT, :); 
            inf = sol(Mesh.numX + Mesh.numT + 1 : Mesh.numX + (2 * Mesh.numT), :); 
            chem = 0;
            
            vir_b = sol(Mesh.numX + (2 * Mesh.numT) + 1, :);
            tar_b = sol(Mesh.numX + (2 * Mesh.numT) + 2, :);
            inf_b = sol(Mesh.numX + (2 * Mesh.numT) + 3, :);
        end
    end
end

