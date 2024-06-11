classdef Vagina
    properties
        thicknessLumen {mustBeNumeric}
        thicknessEpithelium {mustBeNumeric}
        thicknessStroma {mustBeNumeric}
        thicknessTissue {mustBeNumeric}
        surfaceArea {mustBeNumeric}
        width {mustBeNumeric}
        length {mustBeNumeric}
        volumeFractionCellsEpithelium {mustBeNumeric}
        volumeFractionCellsStroma {mustBeNumeric}
        volumeVaginalFluid {mustBeNumeric}
        h_L {mustBeNumeric}
        h_E {mustBeNumeric}
        h_S {mustBeNumeric}
        h_T {mustBeNumeric}
        W {mustBeNumeric}
        L {mustBeNumeric}
        SA {mustBeNumeric}
        Vf {mustBeNumeric}
        phi_E {mustBeNumeric}
        phi_S {mustBeNumeric}
    end
    methods
        function obj = setSymbols(obj)
            obj.h_L = obj.thicknessLumen;
            obj.h_E = obj.thicknessEpithelium;
            obj.h_S = obj.thicknessStroma;
            obj.h_T = obj.thicknessTissue;
            obj.Vf = obj.volumeVaginalFluid;
            obj.W = obj.width;
            obj.L = obj.length;
            obj.SA = obj.surfaceArea;
            obj.phi_E = obj.volumeFractionCellsEpithelium;
            obj.phi_S = obj.volumeFractionCellsStroma;
        end
    end
end