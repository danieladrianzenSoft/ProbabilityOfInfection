classdef getMeshLib
    methods(Static)
        function [Mesh, totalSize] = createMesh_LES_IVR(geometry,meshSize,includeOutputs)
            arguments
               geometry,
               meshSize,
               includeOutputs {mustBeNumeric} = 0
            end

            Mesh = Discretization1D;
 
            %Mesh.numEpithelium = ceil(((geometry.h_L + geometry.h_E) / (geometry.h_L + geometry.h_T)) * meshSize) - ...
            %   (ceil((geometry.h_L / (geometry.h_L + geometry.h_T) * meshSize)));
            Mesh.numEpithelium = ceil((geometry.h_E)/(geometry.h_L + geometry.h_T) * meshSize);
            if Mesh.numEpithelium < 40
                Mesh.numEpithelium = 40;
            end
            Mesh.numLumen = ceil((geometry.h_L/(geometry.h_L+geometry.h_T))*meshSize);
            if Mesh.numLumen < 40
                Mesh.numLumen = 40;
            end
            Mesh.numResin = 200;

            Mesh.numStroma = meshSize - (Mesh.numResin + Mesh.numLumen + Mesh.numEpithelium);
            Mesh.numTissue = Mesh.numEpithelium + Mesh.numStroma;
            Mesh.numElements = meshSize - Mesh.numResin;

            %%% CHECK IF Mesh.x NEEDS TO BE EDITED FOR PLOTTING PURPOSES
            Mesh.x = [linspace(0, geometry.h_L, Mesh.numLumen), ...
                 linspace(geometry.h_L + (geometry.h_E / Mesh.numEpithelium), geometry.h_L + geometry.h_E, Mesh.numEpithelium) ...
                 linspace(geometry.h_L + geometry.h_E + (geometry.h_S / Mesh.numStroma), geometry.h_L + geometry.h_T, Mesh.numStroma)];
            Mesh = Mesh.setSymbols();

            if includeOutputs == 0
                totalSize = Mesh.numR + 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 5; 
            else
                totalSize = Mesh.numR + 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 9; 
            end

        end
        function [Mesh, totalSize] = createMesh_LES_Gel(geometry,meshSize,includeOutputs)
            arguments
               geometry,
               meshSize,
               includeOutputs {mustBeNumeric} = 0
            end

            Mesh = Discretization1D;
            Mesh.numElements = meshSize;
            Mesh.numEpithelium = ceil(((geometry.h_L + geometry.h_E) / (geometry.h_L + geometry.h_T)) * meshSize) - ...
               (ceil((geometry.h_L / (geometry.h_L + geometry.h_T) * meshSize)));
            if Mesh.numEpithelium < 40
                Mesh.numEpithelium = 40;
            end
            Mesh.numLumen = ceil((geometry.h_L/(geometry.h_L+geometry.h_T))*meshSize);
            if Mesh.numLumen < 40
                Mesh.numLumen = 40;
            end

            Mesh.numStroma = meshSize - (Mesh.numLumen + Mesh.numEpithelium);
            Mesh.numTissue = Mesh.numEpithelium + Mesh.numStroma;
            Mesh.x = [linspace(0, geometry.h_L, Mesh.numLumen), ...
                 linspace(geometry.h_L + (geometry.h_E / Mesh.numEpithelium), geometry.h_L + geometry.h_E, Mesh.numEpithelium) ...
                 linspace(geometry.h_L + geometry.h_E + (geometry.h_S / Mesh.numStroma), geometry.h_L + geometry.h_T, Mesh.numStroma)];
            Mesh = Mesh.setSymbols();

            if includeOutputs == 0
                totalSize = 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 5; 
            else
                totalSize = 2 * Mesh.numX + Mesh.numT + 2 * Mesh.numS + 9; 
            end

        end
        function [Mesh, totalSize] = createMesh_LES(geometry, meshSize, includeOutputs)
            arguments
               geometry,
               meshSize,
               includeOutputs {mustBeNumeric} = 0
            end

            Mesh = Discretization1D;
            Mesh.numElements = meshSize;
            Mesh.numEpithelium = ceil(((geometry.h_L + geometry.h_E) / (geometry.h_L + geometry.h_T)) * meshSize) - ...
               (ceil((geometry.h_L / (geometry.h_L + geometry.h_T) * meshSize)));
            if Mesh.numEpithelium < 40
                Mesh.numEpithelium = 40;
            end
            Mesh.numLumen = ceil((geometry.h_L/(geometry.h_L+geometry.h_T))*meshSize);
            if Mesh.numLumen < 40
                Mesh.numLumen = 40;
            end

            Mesh.numStroma = meshSize - (Mesh.numLumen + Mesh.numEpithelium);
            Mesh.numTissue = Mesh.numEpithelium + Mesh.numStroma;
            Mesh.x = [linspace(0, geometry.h_L, Mesh.numLumen), ...
                 linspace(geometry.h_L + (geometry.h_E / Mesh.numEpithelium), geometry.h_L + geometry.h_E, Mesh.numEpithelium) ...
                 linspace(geometry.h_L + geometry.h_E + (geometry.h_S / Mesh.numStroma), geometry.h_L + geometry.h_T, Mesh.numStroma)];
            Mesh = Mesh.setSymbols();

            if includeOutputs == 0
                totalSize = Mesh.numX + Mesh.numS + Mesh.numS + 3; 
            else
                totalSize = Mesh.numX + Mesh.numS + Mesh.numS + 7; 
            end
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