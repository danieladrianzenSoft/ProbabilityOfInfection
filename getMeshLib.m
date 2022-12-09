classdef getMeshLib
    methods(Static)
        function [Mesh, totalSize] = createMesh_LES(geometry, meshSize)
            Mesh = Discretization1D;
            Mesh.numElements = meshSize;
            Mesh.numLumen = ceil((geometry.h_L/(geometry.h_L+geometry.h_T))*meshSize);
            Mesh.numEpithelium = ceil(((geometry.h_L + geometry.h_E) / (geometry.h_L + geometry.h_T)) * meshSize) - ...
               (ceil((geometry.h_L / (geometry.h_L + geometry.h_T) * meshSize)));
            Mesh.numStroma = meshSize - (Mesh.numLumen + Mesh.numEpithelium);
            Mesh.numTissue = Mesh.numEpithelium + Mesh.numStroma;
            Mesh.x = [linspace(0, geometry.h_L, Mesh.numLumen), ...
                 linspace(geometry.h_L + (geometry.h_E / Mesh.numEpithelium), geometry.h_L + geometry.h_E, Mesh.numEpithelium) ...
                 linspace(geometry.h_L + geometry.h_E + (geometry.h_S / Mesh.numStroma), geometry.h_L + geometry.h_T, Mesh.numStroma)];
            Mesh = Mesh.setSymbols();
            totalSize = Mesh.numX + Mesh.numS + Mesh.numS + 3; 
        end
        function [Mesh, totalSize] = createMesh_LT(geometry, meshSize)
            Mesh = Discretization1D;
            Mesh.numElements = meshSize;
            Mesh.numLumen = ceil((geometry.h_L/(geometry.h_L+geometry.h_T))*meshSize);
            Mesh.numEpithelium = ceil(((geometry.h_L + geometry.h_E) / (geometry.h_L + geometry.h_T)) * meshSize) - ...
               (ceil((geometry.h_L / (geometry.h_L + geometry.h_T) * meshSize)));
            Mesh.numStroma = meshSize - (Mesh.numLumen + Mesh.numEpithelium);
            Mesh.numTissue = Mesh.numEpithelium + Mesh.numStroma;
            Mesh.x = [linspace(0, geometry.h_L, Mesh.numLumen), ...
                 linspace(geometry.h_L + (geometry.h_E / Mesh.numEpithelium), geometry.h_L + geometry.h_E, Mesh.numEpithelium) ...
                 linspace(geometry.h_L + geometry.h_E + (geometry.h_S / Mesh.numStroma), geometry.h_L + geometry.h_T, Mesh.numStroma)];
            Mesh = Mesh.setSymbols();
            totalSize = Mesh.numX + Mesh.numT + Mesh.numT + 3; 
        end
    end
end