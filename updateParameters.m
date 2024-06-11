function [V, T, I, U, D, Geometry] = updateParameters(tSpan,V,T,I,U,D,Geometry,params)
    if params.C0 == 0 || params.drugType == 0
        V.dilutionFractionLumen = params.Vs / (params.Vf + params.Vs);
        D.dilutionFractionLumen = 0;
        Geometry.thicknessLumen = 1/2 * ((params.Vf + params.Vs) / Geometry.surfaceArea);
    elseif params.drugType == 1
        if params.T_VD < 0
            if (tSpan(1) == 0)
                V.dilutionFractionLumen = 0;
                D.dilutionFractionLumen = params.Vg / (params.Vf + params.Vg);
                Geometry.thicknessLumen = 1/2 * ((params.Vf + params.Vg) / Geometry.surfaceArea);
            else
                V.dilutionFractionLumen = params.Vs / (params.Vf + params.Vs + params.Vg);
                D.dilutionFractionLumen = params.Vg / (params.Vf + params.Vs + params.Vg);
                Geometry.thicknessLumen = 1/2 * ((params.Vf + params.Vs + params.Vg) / Geometry.surfaceArea);
            end
        elseif params.T_VD > 0
            if (tSpan(1) == 0)
                V.dilutionFractionLumen = params.Vs / (params.Vf + params.Vs);
                D.dilutionFractionLumen = 0;
                Geometry.thicknessLumen = 1/2 * ((params.Vf + params.Vs) / Geometry.surfaceArea);
            else
                V.dilutionFractionLumen = params.Vs / (params.Vf + params.Vs + params.Vg);
                D.dilutionFractionLumen = params.Vg / (params.Vf + params.Vs + params.Vg);
                Geometry.thicknessLumen = 1/2 * ((params.Vf + params.Vs + params.Vg) / Geometry.surfaceArea);
            end
        else
            V.dilutionFractionLumen = params.Vs / (params.Vf + params.Vs + params.Vg);
            D.dilutionFractionLumen = params.Vg / (params.Vf + params.Vs + params.Vg);
            Geometry.thicknessLumen = 1/2 * ((params.Vf + params.Vs + params.Vg) / Geometry.surfaceArea);
        end
    elseif params.drugType == 2
        if params.T_VD < 0
            if (tSpan(1) == 0)
                V.dilutionFractionLumen = 0;
                D.dilutionFractionLumen = 0;
                Geometry.thicknessLumen = 1/2 * ((params.Vf) / Geometry.surfaceArea);
            else
                V.dilutionFractionLumen = params.Vs / (params.Vf + params.Vs);
                D.dilutionFractionLumen = 0;
                Geometry.thicknessLumen = 1/2 * ((params.Vf + params.Vs) / Geometry.surfaceArea);
            end
        elseif params.T_VD > 0
            % if (tSpan(1) == 0)
            %     V.dilutionFractionLumen = 0;
            %     D.dilutionFractionLumen = params.Vs / (params.Vf + params.Vs);
            %     Geometry.thicknessLumen = 1/2 * ((params.Vf + params.Vs) / Geometry.surfaceArea);
            % else
                V.dilutionFractionLumen = 0;
                D.dilutionFractionLumen = params.Vs / (params.Vf + params.Vs);
                Geometry.thicknessLumen = 1/2 * ((params.Vf + params.Vs) / Geometry.surfaceArea);
            %end
        else
            V.dilutionFractionLumen = params.Vs / (params.Vf + params.Vs);
            D.dilutionFractionLumen = 0;
            Geometry.thicknessLumen = 1/2 * ((params.Vf + params.Vs) / Geometry.surfaceArea);
        end

    end
    V = V.setSymbols();
    D = D.setSymbols();
    Geometry = Geometry.setSymbols();
end