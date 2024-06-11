classdef getICLib
    methods(Static)
        function IC = getIC_EpiStroma_IVR_TVD(Mesh, V, T, D, T_VD, sol, getResults, includeOutputs)
            arguments
               Mesh,
               V,
               T
               D,
               T_VD,
               sol,
               getResults,
               includeOutputs {mustBeInteger} = 0
            end

            if T_VD == 0 % no need to use this function
                fprintf('/nT_VD = 0, second IC unnecessary\n')
                return
            end

            % get destructured results, whether we want to include
            % additional outputs (i.e. # virus-cell infections in tissue, etc) or not
            if (includeOutputs == 1) 
                [vir, drug, drugdp, tar, inf, chem, vir_b, drug_b, drugdp_b, tar_b, inf_b, vt_t, it_t, vt_b, it_b] = getResults(sol, Mesh, 2, includeOutputs);
            else
                [vir, drug, drugdp, tar, inf, chem, vir_b, drug_b, drugdp_b, tar_b, inf_b, ~] = getResults(sol, Mesh, 2, includeOutputs);
            end

            if T_VD < 0
                IC = [(V.V0_L.*V.dF_L).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numEpithelium), 0.*ones(1,Mesh.numStroma)... %virions
                      drug(end, 1:Mesh.numR), drug(end, Mesh.numR + 1 : Mesh.numR + Mesh.numL).*D.dF_L, drug(end, Mesh.numR + Mesh.numL + 1 : Mesh.numR + Mesh.numX)... %drug TFV
                      drugdp(end, 1:Mesh.numE+Mesh.numS)...% drug TFV-DP
                      tar(end,:),... %target cells
                      inf(end,:),... %cells (target, infected)
                      0, drug_b(end), drugdp_b(end), tar_b(end), inf_b(end)]; %V, D(TFV), D(TFV-DP), T, I in blood
            elseif T_VD > 0
                IC = [vir(end,1:Mesh.numL).*V.dF_L, vir(end,Mesh.numL+1:Mesh.numX)... %virions
                      (D.C0).*ones(1,Mesh.numR), (0).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %drug TFV
                      0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)...% drug TFV-DP
                      tar(end,:),... %target cells
                      inf(end,:),... %cells (target, infected)
                      vir_b(end), 0, 0, tar_b(end), inf_b(end)]; %V, D(TFV), D(TFV-DP), T, I in blood
            else
                fprintf('/nError trying to create IC for T_VD != 0\n')
            end
            
            % append IC for additional outputs if required
            if (includeOutputs == 1)
                IC = [IC, vt_t(end), it_t(end), vt_b(end), it_b(end)];
            end
        end
        function IC = getIC_EpiStroma_IVR(Mesh, V, T, D, T_VD, includeOutputs)
            arguments
               Mesh,
               V,
               T
               D,
               T_VD,
               includeOutputs {mustBeInteger} = 0
            end
            if T_VD == 0
                % SAME TIME, D IN R, NOT IN L
                IC = [(V.V0_L.*V.dF_L).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %virions
                      (D.C0).*ones(1,Mesh.numR), 0.*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %drug TFV
                      0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)...% drug TFV-DP
                      T.T0_T.*ones(1,Mesh.numS),... %target cells
                      0.*ones(1,Mesh.numS)... %cells (target, infected)
                      0, 0, 0, T.T0_B, 0]; %V, D(TFV), D(TFV-DP), T, I in blood
            elseif T_VD < 0
                % D APPLIED BEFORE V, V=0 Initially
                IC = [(0).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %virions
                  (D.C0).*ones(1,Mesh.numR), (0).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %drug TFV
                  0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)...% drug TFV-DP
                  T.T0_T.*ones(1,Mesh.numS),... %target cells
                  0.*ones(1,Mesh.numS)... %cells (target, infected)
                  0, 0, 0, T.T0_B, 0]; %V, D(TFV), D(TFV-DP), T, I in blood
            elseif T_VD > 0
                % D APPLIED AFTER V, D=0 Initially
                IC = [(V.V0_L.*V.dF_L).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %virions
                  (0).*ones(1,Mesh.numR), (0).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %drug TFV
                  0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)...% drug TFV-DP
                  T.T0_T.*ones(1,Mesh.numS),... %target cells
                  0.*ones(1,Mesh.numS)... %cells (target, infected)
                  0, 0, 0, T.T0_B, 0]; %V, D(TFV), D(TFV-DP), T, I in blood
            else
                fprintf('/nError trying to create IC\n')
            end
            % append IC for additional outputs if required
            if (includeOutputs == 1)
                IC = [IC, 0, 0, 0, 0];
            end
        end
        function IC = getIC_EpiStroma_Gel_TVD(Mesh, V, T, D, T_VD, sol, getResults, includeOutputs)
            arguments
               Mesh,
               V,
               T
               D,
               T_VD,
               sol,
               getResults,
               includeOutputs {mustBeInteger} = 0
            end

            if T_VD == 0 % no need to use this function
                fprintf('/nT_VD = 0, second IC unnecessary\n')
                return
            end

            % get destructured results, whether we want to include
            % additional outputs (i.e. # virus-cell infections in tissue, etc) or not
            if (includeOutputs == 1) 
                [vir, drug, drugdp, tar, inf, chem, vir_b, drug_b, drugdp_b, tar_b, inf_b, vt_t, it_t, vt_b, it_b] = getResults(sol, Mesh, 2, includeOutputs);
            else
                [vir, drug, drugdp, tar, inf, chem, vir_b, drug_b, drugdp_b, tar_b, inf_b, ~] = getResults(sol, Mesh, 2, includeOutputs);
            end

            if T_VD < 0
                IC = [(V.V0_L.*V.dF_L).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numEpithelium), 0.*ones(1,Mesh.numStroma)... %virions
                      drug(end, 1:Mesh.numL).*D.dF_L, drug(end, Mesh.numL+1:Mesh.numX)... %drug TFV
                      drugdp(end, 1:Mesh.numE+Mesh.numS)...% drug TFV-DP
                      tar(end,:),... %target cells
                      inf(end,:),... %cells (target, infected)
                      0, drug_b(end), drugdp_b(end), tar_b(end), inf_b(end)]; %V, D(TFV), D(TFV-DP), T, I in blood
            elseif T_VD > 0
                IC = [vir(end,1:Mesh.numL).*V.dF_L, vir(end,Mesh.numL+1:Mesh.numX)... %virions
                      (D.C0.*D.dF_L).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %drug TFV
                      0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)...% drug TFV-DP
                      tar(end,:),... %target cells
                      inf(end,:),... %cells (target, infected)
                      vir_b(end), 0, 0, tar_b(end), inf_b(end)]; %V, D(TFV), D(TFV-DP), T, I in blood
            else
                fprintf('/nError trying to create IC for T_VD != 0\n')
            end
            
            % append IC for additional outputs if required
            if (includeOutputs == 1)
                IC = [IC, vt_t(end), it_t(end), vt_b(end), it_b(end)];
            end
        end
        function IC = getIC_EpiStroma_Gel(Mesh, V, T, D, T_VD, includeOutputs)
            arguments
               Mesh,
               V,
               T
               D,
               T_VD,
               includeOutputs {mustBeInteger} = 0
            end
            if T_VD == 0
                % SAME TIME, V AND D ARE DILUTED 1:1
                IC = [(V.V0_L.*V.dF_L).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %virions
                      (D.C0.*D.dF_L).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %drug TFV
                      0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)...% drug TFV-DP
                      T.T0_T.*ones(1,Mesh.numS),... %target cells
                      0.*ones(1,Mesh.numS)... %cells (target, infected)
                      0, 0, 0, T.T0_B, 0]; %V, D(TFV), D(TFV-DP), T, I in blood
            elseif T_VD < 0
                % D APPLIED BEFORE V, V=0 Initially
                IC = [(0).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %virions
                  (D.C0.*D.dF_L).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %drug TFV
                  0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)...% drug TFV-DP
                  T.T0_T.*ones(1,Mesh.numS),... %target cells
                  0.*ones(1,Mesh.numS)... %cells (target, infected)
                  0, 0, 0, T.T0_B, 0]; %V, D(TFV), D(TFV-DP), T, I in blood
            elseif T_VD > 0
                % D APPLIED AFTER V, D=0 Initially
                IC = [(V.V0_L.*V.dF_L).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %virions
                  (0).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %drug TFV
                  0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)...% drug TFV-DP
                  T.T0_T.*ones(1,Mesh.numS),... %target cells
                  0.*ones(1,Mesh.numS)... %cells (target, infected)
                  0, 0, 0, T.T0_B, 0]; %V, D(TFV), D(TFV-DP), T, I in blood
            else
                fprintf('/nError trying to create IC\n')
            end
            % append IC for additional outputs if required
            if (includeOutputs == 1)
                IC = [IC, 0, 0, 0, 0];
            end
        end
        function IC = getIC_EpiStroma(Mesh, V, T, D, T_VD, includeOutputs)
            arguments
               Mesh,
               V,
               T
               D,
               T_VD,
               includeOutputs {mustBeInteger} = 0
            end
            IC = [(V.V0_L.*V.dF_L).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %virions
                  T.T0_T.*ones(1,Mesh.numS),... %target cells
                  0.*ones(1,Mesh.numS)... %cells (target, infected)
                  0, T.T0_B, 0]; %V, T, I in blood
            % append IC for additional outputs if required
            if (includeOutputs == 1)
                IC = [IC, 0, 0, 0, 0];
            end
        end
        function IC = getIC_tissue(Mesh, V, T, D, T_VD)
            IC = [(V.V0_L.*V.dF_L).*ones(1,Mesh.numL), 0.*ones(1,Mesh.numT),... %virions
                  T.T0_T.*ones(1,Mesh.numT), 0.*ones(1,Mesh.numT),... %cells (target, infected)
                  0, T.T0_B, 0]; %V, T, I in blood
        end
    end
end