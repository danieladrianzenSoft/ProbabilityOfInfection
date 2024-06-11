classdef Virus
    properties
        diffCoeffLumen {mustBeNumeric}
        diffCoeffEpithelium {mustBeNumeric}
        diffCoeffStroma {mustBeNumeric}
        diffCoeffTissue {mustBeNumeric}
        dilutionFractionLumen {mustBeNumeric}
        lossDilutionInLumen {mustBeNumeric}
        lossTissueToBlood {mustBeNumeric}
        lossClearanceInBlood {mustBeNumeric}
        partitionLumenTissue {mustBeNumeric}
        volumeOfDistribution {mustBeNumeric}
        volumeSemen {mustBeNumeric}
        initialLumenConcentration {mustBeNumeric}
        radius {mustBeNumeric}
        ratioCollisionsCausingInfection {mustBeNumeric}
        infectivity {mustBeNumeric}
        DV_L {mustBeNumeric}
        DV_E {mustBeNumeric}
        DV_S {mustBeNumeric}
        DV_T {mustBeNumeric}
        dF_L {mustBeNumeric}
        kD {mustBeNumeric}
        kB {mustBeNumeric}
        kL {mustBeNumeric}
        phi_LT {mustBeNumeric}
        Vb {mustBeNumeric}
        Vs {mustBeNumeric}
        V0_L {mustBeNumeric}
        r {mustBeNumeric}
        s {mustBeNumeric}
        beta {mustBeNumeric}
    end
    methods
        function obj = setSymbols(obj)
            obj.DV_L = obj.diffCoeffLumen;
            obj.DV_T = obj.diffCoeffTissue;
            obj.DV_E = obj.diffCoeffEpithelium;
            obj.DV_S = obj.diffCoeffStroma;
            obj.dF_L = obj.dilutionFractionLumen;
            obj.kD = obj.lossDilutionInLumen;
            obj.kB = obj.lossTissueToBlood;
            obj.kL = obj.lossClearanceInBlood;
            obj.phi_LT = obj.partitionLumenTissue;
            obj.Vb = obj.volumeOfDistribution;
            obj.Vs = obj.volumeSemen;
            obj.V0_L = obj.initialLumenConcentration;
            obj.r = obj.radius;
            obj.s = obj.ratioCollisionsCausingInfection;
            obj.beta = obj.infectivity;
        end
    end
end

