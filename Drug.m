classdef Drug
    properties
        thicknessResin {mustBeNumeric}
        effectiveSurfaceArea {mustBeNumeric}
        diffCoeffResin {mustBeNumeric}
        diffCoeffLumen {mustBeNumeric}
        diffCoeffEpithelium {mustBeNumeric}
        diffCoeffStroma {mustBeNumeric}
        dilutionFractionLumen {mustBeNumeric}
        lossDilutionInLumen {mustBeNumeric}
        lossStromaToBlood {mustBeNumeric}
        lossClearanceInBlood {mustBeNumeric}
        partitionResinLumen {mustBeNumeric}
        partitionLumenEpthelium {mustBeNumeric}
        partitionEpitheliumStroma {mustBeNumeric}
        volumeGel {mustBeNumeric}
        volumeOfDistribution {mustBeNumeric}
        initialConcentration {mustBeNumeric}
        rateActivation {mustBeNumeric}
        rateDeactivation {mustBeNumeric}
        ratioActivationFromDeactivated {mustBeNumeric}
        h_R {mustBeNumeric}
        SA {mustBeNumeric}
        Dd_R {mustBeNumeric}
        Dd_L {mustBeNumeric}
        Dd_E {mustBeNumeric}
        Dd_S {mustBeNumeric}
        dF_L {mustBeNumeric}
        kd_D {mustBeNumeric}
        kd_B {mustBeNumeric}
        kd_L {mustBeNumeric}
        phi_RL {mustBeNumeric}
        phi_LE {mustBeNumeric}
        phi_ES {mustBeNumeric}
        Vb {mustBeNumeric}
        Vg {mustBeNumeric}
        C0 {mustBeNumeric}
        Kon {mustBeNumeric}
        Koff {mustBeNumeric}
        r
    end
%phid_GE,phid_ES,kd_D,kd_B,kd_L,Vb

    methods
        function obj = setSymbols(obj)
            obj.h_R = obj.thicknessResin;
            obj.SA = obj.effectiveSurfaceArea;
            obj.Dd_R = obj.diffCoeffResin;
            obj.Dd_L = obj.diffCoeffLumen;
            obj.Dd_E = obj.diffCoeffEpithelium;
            obj.Dd_S = obj.diffCoeffStroma;
            obj.dF_L = obj.dilutionFractionLumen;
            obj.kd_D = obj.lossDilutionInLumen;
            obj.kd_B = obj.lossStromaToBlood;
            obj.kd_L = obj.lossClearanceInBlood;
            obj.phi_RL = obj.partitionResinLumen;
            obj.phi_LE = obj.partitionLumenEpthelium;
            obj.phi_ES = obj.partitionEpitheliumStroma;
            obj.Vb = obj.volumeOfDistribution;
            obj.Vg = obj.volumeGel;
            obj.C0 = obj.initialConcentration;
            obj.Kon = obj.rateActivation;
            obj.Koff = obj.rateDeactivation;
            obj.r = obj.ratioActivationFromDeactivated;
        end
    end
end