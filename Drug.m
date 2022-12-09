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
    end
end