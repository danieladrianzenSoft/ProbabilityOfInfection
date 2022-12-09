classdef Vagina
    properties
        thicknessLumen {mustBeNumeric}
        thicknessEpithelium {mustBeNumeric}
        thicknessStroma {mustBeNumeric}
        thicknessTissue {mustBeNumeric}
        width {mustBeNumeric}
        length {mustBeNumeric}
        h_L {mustBeNumeric}
        h_E {mustBeNumeric}
        h_S {mustBeNumeric}
        h_T {mustBeNumeric}
        W {mustBeNumeric}
        L {mustBeNumeric}
    end
    methods
        function obj = setSymbols(obj)
            obj.h_L = obj.thicknessLumen;
            obj.h_E = obj.thicknessEpithelium;
            obj.h_S = obj.thicknessStroma;
            obj.h_T = obj.thicknessTissue;
            obj.W = obj.width;
            obj.L = obj.length;
        end
    end
end