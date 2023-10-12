function [cm] = getCustomColormap(numColors,options)
    arguments
        numColors (1,1) {mustBeNumeric}
        options.colormap = ""
    end

    if options.colormap == "teal"
        cm_vec = [161,241,209;108,218,174;84,173,152;75,157,145;50,101,121;48,95,118]/255;
    elseif options.colormap == "purple"
        cm_vec = [80,64,77;72,61,139;102,51,153;106,90,205;150,111,214;177,156,217]/255;
    elseif options.colormap == "purple-teal"
        cm_vec = [80,64,77;106,90,205;177,156,217;161,241,209;84,173,152;48,95,118]/255;
    elseif options.colormap == "purple-teal-pink"
        %cm_vec = [80,64,77;106,90,205;177,156,217;48,95,118;84,173,152;161,241,209]/255;
        cm_vec = [80,64,77;106,90,205;177,156,217;161,241,209;84,173,152;48,95,118;]/255;
        %cm_vec = [80,64,77;106,90,205;177,156,217;227,151,127]/255;
    else
        cm_vec = [80,64,77;227,151,127;108,160,220;175,203,169;196,174,173;150,111,214]/255;
    end
    vec = linspace(100,0,length(cm_vec));
    %cm = interp1(vec,cm_vec,linspace(max(vec),min(vec),numColors));
    if numColors > 6
        cm = interp1(vec,cm_vec,linspace(max(vec),min(vec),numColors));
    else
        cm = cm_vec(1:numColors,:);
    end
end 