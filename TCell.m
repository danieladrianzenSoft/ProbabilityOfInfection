classdef TCell
    properties
        initialTissueConcentration {mustBeNumeric}
        initialBloodConcentration {mustBeNumeric}
        volumeFractionTissue {mustBeNumeric}
        deathRateTissue {mustBeNumeric}
        deathRateBlood {mustBeNumeric}
        productionRateTissue {mustBeNumeric}
        productionRateBlood {mustBeNumeric}
        radius {mustBeNumeric}
        diffCoeffTissue {mustBeNumeric}
        uptakeRateFromBlood {mustBeNumeric}
        maxChemotaxSpeed {mustBeNumeric}
        sensitivityToChemokine {mustBeNumeric}
        T0_T {mustBeNumeric}
        T0_B {mustBeNumeric}
        phiT {mustBeNumeric}
        dT_T {mustBeNumeric}
        dT_B {mustBeNumeric}
        lambda_T {mustBeNumeric}
        lambda_B {mustBeNumeric}
        r {mustBeNumeric}
        DT_T {mustBeNumeric}
        Ku {mustBeNumeric}
        Stax {mustBeNumeric}
        Sigma {mustBeNumeric}
    end
    methods
        function obj = setSymbols(obj)
            obj.T0_T = obj.initialTissueConcentration;
            obj.T0_B = obj.initialBloodConcentration;
            obj.phiT = obj.volumeFractionTissue;
            obj.dT_T = obj.deathRateTissue;
            obj.dT_B = obj.deathRateBlood;
            obj.lambda_T = obj.productionRateTissue;
            obj.lambda_B = obj.productionRateBlood;
            obj.r = obj.radius;
            obj.DT_T = obj.diffCoeffTissue;
            obj.Ku = obj.uptakeRateFromBlood;
            obj.Stax = obj.maxChemotaxSpeed;
            obj.Sigma = obj.sensitivityToChemokine;
        end
        function volume = getVolume(obj)
            volume = 4/3 * pi * (obj.radius)^3;
        end
    end
end

