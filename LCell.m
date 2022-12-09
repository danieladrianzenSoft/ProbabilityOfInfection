classdef LCell
    properties
        deathRateStroma {mustBeNumeric}
        deathRateBlood {mustBeNumeric}
        ratioInfectionsResultingInLatency {mustBeNumeric}
        transitionRateLtoI {mustBeNumeric}
        proliferationRate {mustBeNumeric}
    end
end

