classdef Drug
    properties
        diffCoeffLumen {mustBeNumeric}
        diffCoeffEpithelium {mustBeNumeric}
        diffCoeffStroma {mustBeNumeric}
        lossDilutionInLumen {mustBeNumeric}
        lossStromaToBlood {mustBeNumeric}
        lossClearanceInBlood {mustBeNumeric}
        partitionLumenEpthelium {mustBeNumeric}
        partitionEpitheliumStroma {mustBeNumeric}
        volumeOfDistribution {mustBeNumeric}
        initialLumenConcentration {mustBeNumeric}
        rateActivation {mustBeNumeric}
        rateDeactivation {mustBeNumeric}
        ratioActivationFromDeactivated {mustBeNumeric}
        Dd_L {mustBeNumeric}
        Dd_E {mustBeNumeric}
        Dd_S {mustBeNumeric}
        kd_D {mustBeNumeric}
        kd_B {mustBeNumeric}
        kd_L {mustBeNumeric}
        phi_LE {mustBeNumeric}
        phi_ES {mustBeNumeric}
        Vb {mustBeNumeric}
        C0_L {mustBeNumeric}
        Kon {mustBeNumeric}
        Koff {mustBeNumeric}
        r
    end
%phid_GE,phid_ES,kd_D,kd_B,kd_L,Vb

    methods
        function obj = setSymbols(obj)
            obj.Dd_L = obj.diffCoeffLumen;
            obj.Dd_E = obj.diffCoeffEpithelium;
            obj.Dd_S = obj.diffCoeffStroma;
            obj.kd_D = obj.lossDilutionInLumen;
            obj.kd_B = obj.lossStromaToBlood;
            obj.kd_L = obj.lossClearanceInBlood;
            obj.phi_LE = obj.partitionLumenEpthelium;
            obj.phi_ES = obj.partitionEpitheliumStroma;
            obj.Vb = obj.volumeOfDistribution;
            obj.C0_L = obj.initialLumenConcentration;
            obj.Kon = obj.rateActivation;
            obj.Koff = obj.rateDeactivation;
            obj.r = obj.ratioActivationFromDeactivated;
        end
    end
end