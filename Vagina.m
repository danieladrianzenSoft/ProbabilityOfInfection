classdef Vagina
    properties
        thicknessLumen {mustBeNumeric}
        thicknessEpithelium {mustBeNumeric}
        thicknessStroma {mustBeNumeric}
        thicknessTissue {mustBeNumeric}
        width {mustBeNumeric}
        length {mustBeNumeric}
        volumeFractionCellsEpithelium {mustBeNumeric};
        volumeFractionCellsStroma {mustBeNumeric};
        h_L {mustBeNumeric}
        h_E {mustBeNumeric}
        h_S {mustBeNumeric}
        h_T {mustBeNumeric}
        W {mustBeNumeric}
        L {mustBeNumeric}
        phi_E {mustBeNumeric}
        phi_S {mustBeNumeric}
    end
    methods
        function obj = setSymbols(obj)
            obj.h_L = obj.thicknessLumen;
            obj.h_E = obj.thicknessEpithelium;
            obj.h_S = obj.thicknessStroma;
            obj.h_T = obj.thicknessTissue;
            obj.W = obj.width;
            obj.L = obj.length;
            obj.phi_E = obj.volumeFractionCellsEpithelium;
            obj.phi_S = obj.volumeFractionCellsStroma;
        end
    end
end