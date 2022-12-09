classdef ICell
    properties
        productionRateVirus {mustBeNumeric}
        deathRateTissue {mustBeNumeric}
        deathRateBlood {mustBeNumeric}
        diffCoeffTissue {mustBeNumeric}
        lossTissueToBlood {mustBeNumeric}
        cellToCellInfection {mustBeNumeric}
        rho {mustBeNumeric}
        dI_T {mustBeNumeric}
        dI_B {mustBeNumeric}
        DI_T {mustBeNumeric}
        kI_TB {mustBeNumeric}
        w {mustBeNumeric}
    end
    methods
        function obj = setSymbols(obj)
            obj.rho = obj.productionRateVirus;
            obj.dI_T = obj.deathRateTissue;
            obj.dI_B = obj.deathRateBlood;
            obj.DI_T = obj.diffCoeffTissue;
            obj.kI_TB = obj.lossTissueToBlood;
            obj.w = obj.cellToCellInfection;
        end
    end
end