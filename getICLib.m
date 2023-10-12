classdef getICLib
    methods(Static)
        function IC = getIC_EpiStroma_Drug_TVD(Mesh, V, T, D, T_VD, sol, getResults)
            if T_VD == 0
                fprintf('/nT_VD = 0, second IC unnecessary\n')
            elseif T_VD < 0
                [vir, drug, drugdp, tar, inf, chem, vir_b, drug_b, drugdp_b, tar_b, inf_b] = getResults(sol, Mesh, 2);
                IC = [(V.V0_L/2).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numEpithelium), 0.*ones(1,Mesh.numStroma)... %virions
                      drug(end, 1:Mesh.numL)/2, drug(end, Mesh.numL+1:Mesh.numX)... %drug TFV
                      drugdp(end, 1:Mesh.numE+Mesh.numS)...% drug TFV-DP
                      tar(end,:),... %target cells
                      inf(end,:),... %cells (target, infected)
                      0, drug_b(end), drugdp_b(end), tar_b(end), inf_b(end)]; %V, D(TFV), D(TFV-DP), T, I in blood
            elseif T_VD > 0
                [vir, drug, drugdp, tar, inf, chem, vir_b, drug_b, drugdp_b, tar_b, inf_b] = getResults(sol, Mesh, 2);
                IC = [vir(end,1:Mesh.numL)/2, vir(end,Mesh.numL+1:Mesh.numX)... %virions
                      (D.C0_L/2).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %drug TFV
                      0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)...% drug TFV-DP
                      tar(end,:),... %target cells
                      inf(end,:),... %cells (target, infected)
                      vir_b(end), 0, 0, tar_b(end), inf_b(end)]; %V, D(TFV), D(TFV-DP), T, I in blood
            else
                fprintf('/nError trying to create IC for T_VD != 0\n')
            end

        end
        function IC = getIC_EpiStroma_Drug(Mesh, V, T, D, T_VD)
            if T_VD == 0
                % SAME TIME, V AND D ARE DILUTED 1:1
                IC = [(V.V0_L/2).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %virions
                      (D.C0_L/2).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %drug TFV
                      0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)...% drug TFV-DP
                      T.T0_T.*ones(1,Mesh.numS),... %target cells
                      0.*ones(1,Mesh.numS)... %cells (target, infected)
                      0, 0, 0, T.T0_B, 0]; %V, D(TFV), D(TFV-DP), T, I in blood
            elseif T_VD < 0
                % D APPLIED BEFORE V, V=0 Initially
                IC = [(0).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %virions
                  (D.C0_L).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %drug TFV
                  0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)...% drug TFV-DP
                  T.T0_T.*ones(1,Mesh.numS),... %target cells
                  0.*ones(1,Mesh.numS)... %cells (target, infected)
                  0, 0, 0, T.T0_B, 0]; %V, D(TFV), D(TFV-DP), T, I in blood
            elseif T_VD > 0
                % D APPLIED AFTER V, D=0 Initially
                IC = [(V.V0_L).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %virions
                  (0).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %drug TFV
                  0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)...% drug TFV-DP
                  T.T0_T.*ones(1,Mesh.numS),... %target cells
                  0.*ones(1,Mesh.numS)... %cells (target, infected)
                  0, 0, 0, T.T0_B, 0]; %V, D(TFV), D(TFV-DP), T, I in blood
            else
                fprintf('/nError trying to create IC\n')
            end
        end
        function IC = getIC_EpiStroma(Mesh, V, T, D, T_VD)
            IC = [V.V0_L.*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %virions
                  T.T0_T.*ones(1,Mesh.numS),... %target cells
                  0.*ones(1,Mesh.numS)... %cells (target, infected)
                  0, T.T0_B, 0]; %V, T, I in blood
        end
        function IC = getIC_tissue(Mesh, V, T, D, T_VD)
            IC = [V.V0_L.*ones(1,Mesh.numL), 0.*ones(1,Mesh.numT),... %virions
                  T.T0_T.*ones(1,Mesh.numT), 0.*ones(1,Mesh.numT),... %cells (target, infected)
                  0, T.T0_B, 0]; %V, T, I in blood
        end
    end
end