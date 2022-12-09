classdef getSparsityLib
    methods(Static)
        function S = noDrugCellMigration(Mesh, totalSize)
            e = ones(totalSize,1);
            A = spdiags([e e e],-1:1,totalSize,totalSize);
            S = full(A);

            V_SIndices = Mesh.numL + Mesh.numE + 1 : Mesh.numX;
            T_SIndices = Mesh.numX + 1 : Mesh.numX + Mesh.numS;
            I_SIndices = Mesh.numX + Mesh.numS + 1 : Mesh.numX + 2 * Mesh.numS;
            V_BIndices = Mesh.numX + 2 * Mesh.numS + 1;
            %T_BIndices = Mesh.numX + 2 * Mesh.numS + 2;
            I_BIndices = Mesh.numX + 2 * Mesh.numS + 3;

            % V_SIndices
            S(V_SIndices,I_SIndices) = 1;  

            % T_SIndices
            S(T_SIndices,V_SIndices) = 1;
            S(T_SIndices,I_SIndices) = 1;

            % I_SIndices
            S(I_SIndices,V_SIndices) = 1;
            S(I_SIndices,T_SIndices) = 1;

            % V_BIndices
            S(V_BIndices,I_BIndices) = 1;
            %S(V_BIndices,V_SIndices) = 1;
            
            % I_BIndices
            S(I_BIndices,V_BIndices) = 1;

        end
        function S = defaultSparsity(mesh, totalSize)
            % B(numx+1:2*numx+indE+4*indS,1:indG+indE) = 0; %a
            % B(1:indG+indE,numx+1:2*numx+indE+4*indS) = 0; %a
            % B(numx+1:numx+indG,2*numx+1:2*numx+indE+4*indS) = 0; %b
            % B(2*numx+1:2*numx+indE+4*indS,numx+1:numx+indG) = 0; %b
            % B(2*numx+indE+1:2*numx+indE+4*indS,numx+indG+1:numx+indG+indE) = 0; %c
            % B(numx+indG+1:numx+indG+indE,2*numx+indE+1:2*numx+indE+4*indS) = 0; %c
            % B(2*numx+indE+indS+1:2*numx+indE+4*indS,numx+indG+1:numx+indG+indE)=0; %d
            % B(numx+indG+1:numx+indG+indE,2*numx+indE+indS+1:2*numx+indE+4*indS)=0; %d
            % B(2*numx+indE+1:2*numx+indE+4*indS,2*numx+1:2*numx+indE)=0; %e
            % B(2*numx+1:2*numx+indE,2*numx+indE+1:2*numx+indE+4*indS)=0; %e
            % B(2*numx+1:2*numx+indE,numx+indG+indE+1:2*numx)=0; %f
            % B(numx+indG+indE+1:2*numx,2*numx+1:2*numx+indE)=0; %f
            S = ones(totalSize, totalSize);
        end
    end
end