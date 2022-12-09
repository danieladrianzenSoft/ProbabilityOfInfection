classdef Chemokine
    properties
        lossRateBinding {mustBeNumeric}
        lossRateDegradation {mustBeNumeric}
        secretionRate {mustBeNumeric}
        diffCoeffTissue {mustBeNumeric}
        kT {mustBeNumeric}
        kU {mustBeNumeric}
        kS {mustBeNumeric}
        DT {mustBeNumeric}
    end
    methods
        function obj = setSymbols(obj)
            obj.kT = obj.lossRateBinding;
            obj.kU = obj.lossRateDegradation;
            obj.kS = obj.secretionRate;
            obj.DT = obj.diffCoeffTissue;
        end
    end
end