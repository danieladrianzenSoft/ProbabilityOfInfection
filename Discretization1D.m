classdef Discretization1D
    properties
        numElements {mustBeNumeric}
        numResin {mustBeNumeric}
        numLumen {mustBeNumeric}
        numEpithelium {mustBeNumeric}
        numStroma {mustBeNumeric}
        numTissue {mustBeNumeric}
        x {mustBeNumeric}
        numX {mustBeNumeric}
        numR {mustBeNumeric}
        numL {mustBeNumeric}
        numE {mustBeNumeric}
        numS {mustBeNumeric}
        numT {mustBeNumeric}
    end
    methods
        function obj = setSymbols(obj)
            obj.numX = obj.numElements;
            obj.numR = obj.numResin;
            obj.numL = obj.numLumen;
            obj.numE = obj.numEpithelium;
            obj.numS = obj.numStroma;
            obj.numT = obj.numTissue;
        end
    end
end

