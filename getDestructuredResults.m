classdef getDestructuredResults
    methods(Static)
        function [vir, drug, drugdp, tar, inf, chem, vir_b, drug_b, drugdp_b, tar_b, inf_b, vt_t, it_t, vt_b, it_b] = getResults_LES_IVR(sol, Mesh, dir, includeOutputs)
            % no chemokine, 4 compartments: lumen, epi, stroma, blood
            % dir gives the direction of the sol input. 
            % dir = 1: r x t
            % dir = 2: t x r
            arguments
               sol,
               Mesh,
               dir,
               includeOutputs {mustBeNumeric} = 0
            end

            if dir == 1
                vir = sol(1 : Mesh.numX, :);
                drug = sol(Mesh.numX + 1 : Mesh.numR + Mesh.numX + Mesh.numX, :);
                drugdp = sol(Mesh.numR + 2 * Mesh.numX + 1 : Mesh.numR + 2 * Mesh.numX + Mesh.numT, :);
                tar = sol(Mesh.numR + 2 * Mesh.numX + Mesh.numT + 1 : Mesh.numR + 2 * Mesh.numX + Mesh.numT + Mesh.numS, :); 
                inf = sol(Mesh.numR + 2 * Mesh.numX + Mesh.numT + Mesh.numS + 1 : Mesh.numR + 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS, :); 
                chem = [];
                
                vir_b = sol(Mesh.numR + 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 1, :);
                drug_b = sol(Mesh.numR + 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 2, :);
                drugdp_b = sol(Mesh.numR + 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 3, :);
                tar_b = sol(Mesh.numR + 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 4, :);
                inf_b = sol(Mesh.numR + 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 5, :);

                if includeOutputs == 1
                    vt_t = sol(Mesh.numR + 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 6,:); % #virus-cell infections in tissue
                    it_t = sol(Mesh.numR + 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 7,:); % #cell-cell infections in tissue
                    vt_b = sol(Mesh.numR + 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 8,:); % #virus-cell infections in blood
                    it_b = sol(Mesh.numR + 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 9,:); % #cell-cell infections in blood 
                else
                    vt_t = [];
                    it_t = [];
                    vt_b = [];
                    it_b = [];
                end
            elseif dir == 2
                vir = sol(:, 1 : Mesh.numX);
                drug = sol(:, Mesh.numX + 1 : Mesh.numR + Mesh.numX + Mesh.numX);
                drugdp = sol(:, Mesh.numR + 2 * Mesh.numX + 1 : Mesh.numR + 2 * Mesh.numX + Mesh.numT);
                tar = sol(:, Mesh.numR + 2 * Mesh.numX + Mesh.numT + 1 : Mesh.numR + 2 * Mesh.numX + Mesh.numT + Mesh.numS); 
                inf = sol(:, Mesh.numR + 2 * Mesh.numX + Mesh.numT + Mesh.numS + 1 : Mesh.numR + 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS); 
                chem = [];
                
                vir_b = sol(:, Mesh.numR + 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 1);
                drug_b = sol(:, Mesh.numR + 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 2);
                drugdp_b = sol(:, Mesh.numR + 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 3);
                tar_b = sol(:, Mesh.numR + 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 4);
                inf_b = sol(:, Mesh.numR + 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 5);

                if includeOutputs == 1
                    vt_t = sol(:, Mesh.numR + 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 6); % #virus-cell infections in tissue
                    it_t = sol(:, Mesh.numR + 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 7); % #cell-cell infections in tissue
                    vt_b = sol(:, Mesh.numR + 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 8); % #virus-cell infections in blood
                    it_b = sol(:,Mesh.numR + 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 9); % #cell-cell infections in blood     
                else
                    vt_t = [];
                    it_t = [];
                    vt_b = [];
                    it_b = [];
                end
            end
        end 
        function [vir, drug, drugdp, tar, inf, chem, vir_b, drug_b, drugdp_b, tar_b, inf_b, vt_t, it_t, vt_b, it_b] = getResults_LES_Gel(sol, Mesh, dir, includeOutputs)
            % no chemokine, 4 compartments: lumen, epi, stroma, blood
            % dir gives the direction of the sol input. 
            % dir = 1: r x t
            % dir = 2: t x r
            arguments
               sol,
               Mesh,
               dir,
               includeOutputs {mustBeNumeric} = 0
            end
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

                if includeOutputs == 1
                    vt_t = sol(2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 6,:); % #virus-cell infections in tissue
                    it_t = sol(2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 7,:); % #cell-cell infections in tissue
                    vt_b = sol(2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 8,:); % #virus-cell infections in blood
                    it_b = sol(2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 9,:); % #cell-cell infections in blood 
                else
                    vt_t = [];
                    it_t = [];
                    vt_b = [];
                    it_b = [];
                end
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

                if includeOutputs == 1
                    vt_t = sol(:,2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 6); % #virus-cell infections in tissue
                    it_t = sol(:,2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 7); % #cell-cell infections in tissue
                    vt_b = sol(:,2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 8); % #virus-cell infections in blood
                    it_b = sol(:,2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 9); % #cell-cell infections in blood     
                else
                    vt_t = [];
                    it_t = [];
                    vt_b = [];
                    it_b = [];
                end
            end
        end 
        function [vir, drug, drugdp, tar, inf, chem, vir_b, drug_b, drugdp_b, tar_b, inf_b, vt_t, it_t, vt_b, it_b] = getResults_LES(sol, Mesh, dir, includeOutputs)
            % no chemokine, 4 compartments: lumen, epi, stroma, blood
            % dir gives the direction of the sol input. 
            % dir = 1: r x t
            % dir = 2: t x r
            arguments
               sol,
               Mesh,
               dir,
               includeOutputs {mustBeNumeric} = 0
            end
            if (dir == 1)
                vir = sol(1 : Mesh.numX, :);
                tar = sol(Mesh.numX + 1 : Mesh.numX + Mesh.numS, :); 
                inf = sol(Mesh.numX + Mesh.numS + 1 : Mesh.numX + (2 * Mesh.numS), :); 
                vir_b = sol(Mesh.numX + (2 * Mesh.numS) + 1, :);
                tar_b = sol(Mesh.numX + (2 * Mesh.numS) + 2, :);
                inf_b = sol(Mesh.numX + (2 * Mesh.numS) + 3, :);
                
                if includeOutputs == 1
                    vt_t = sol(Mesh.numX + 2 * Mesh.numS + 4,:); % #virus-cell infections in tissue
                    it_t = sol(Mesh.numX + 2 * Mesh.numS + 5,:); % #cell-cell infections in tissue
                    vt_b = sol(Mesh.numX + 2 * Mesh.numS + 6,:); % #virus-cell infections in blood
                    it_b = sol(Mesh.numX + 2 * Mesh.numS + 7,:); % #cell-cell infections in blood   
                else
                    vt_t = [];
                    it_t = [];
                    vt_b = [];
                    it_b = [];
                end
            elseif (dir == 2)
                vir = sol(:, 1 : Mesh.numX);
                tar = sol(:, Mesh.numX + 1 : Mesh.numX + Mesh.numS); 
                inf = sol(:, Mesh.numX + Mesh.numS + 1 : Mesh.numX + (2 * Mesh.numS)); 
                vir_b = sol(:, Mesh.numX + (2 * Mesh.numS) + 1);
                tar_b = sol(:, Mesh.numX + (2 * Mesh.numS) + 2);
                inf_b = sol(:, Mesh.numX + (2 * Mesh.numS) + 3);

                if includeOutputs == 1
                    vt_t = sol(:,Mesh.numX + 2 * Mesh.numS + 4); % #virus-cell infections in tissue
                    it_t = sol(:,Mesh.numX + 2 * Mesh.numS + 5); % #cell-cell infections in tissue
                    vt_b = sol(:,Mesh.numX + 2 * Mesh.numS + 6); % #virus-cell infections in blood
                    it_b = sol(:,Mesh.numX + 2 * Mesh.numS + 7); % #cell-cell infections in blood  
                else
                    vt_t = [];
                    it_t = [];
                    vt_b = [];
                    it_b = [];
                end

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

