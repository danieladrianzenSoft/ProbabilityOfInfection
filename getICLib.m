classdef getICLib
    methods(Static)
        function IC = getIC_EpiStroma(Mesh, V, T)
            IC = [V.V0_L.*ones(1,Mesh.numL), 0.*ones(1,Mesh.numE), 0*ones(1,Mesh.numS)... %virions
                  T.T0_T.*ones(1,Mesh.numS),... %target cells
                  0.*ones(1,Mesh.numS)... %cells (target, infected)
                  0, T.T0_B, 0]; %V, T, I in blood
        end
        function IC = getIC_tissue(Mesh, V, T)
            IC = [V.V0_L.*ones(1,Mesh.numL), 0.*ones(1,Mesh.numT),... %virions
                  T.T0_T.*ones(1,Mesh.numT), 0.*ones(1,Mesh.numT),... %cells (target, infected)
                  0, T.T0_B, 0]; %V, T, I in blood
        end
    end
end