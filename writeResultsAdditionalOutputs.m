function writeResultsAdditionalOutputs(filenameOutputs,overwrite,outputs)
    arguments
       filenameOutputs {mustBeText},
       overwrite {mustBeNumeric},
       outputs.sim_id (:,:) {mustBeNumeric} = [],
       outputs.t (:,:) {mustBeNumeric} = [],
       outputs.vir_tar_t (:,:) {mustBeNumeric} = [],
       outputs.inf_tar_t (:,:) {mustBeNumeric} = [],
       outputs.vir_tar_b (:,:) {mustBeNumeric} = [],
       outputs.inf_tar_b (:,:) {mustBeNumeric} = []
    end

    %variableNames = [];
    %outputTableValues = [];
    outputTable = table();
    if ~isempty(outputs.sim_id)
        %outputTableValues = [outputTableValues,outputs.t]; 
        %variableNames = [variableNames,"t"];
        outputTable = [outputTable,table(outputs.sim_id,VariableNames={'sim_id'})];
    end
    if ~isempty(outputs.t)
        %outputTableValues = [outputTableValues,outputs.t]; 
        %variableNames = [variableNames,"t"];
        outputTable = [outputTable,table(outputs.t,VariableNames={'t'})];
    end
    if ~isempty(outputs.vir_tar_t)
        %outputTableValues = [outputTableValues,outputs.vir_tar_t]; 
        %variableNames = [variableNames,"vir_tar_t"];
        outputTable = [outputTable,table(outputs.vir_tar_t,VariableNames={'vir_tar_t'})];
    end
    if ~isempty(outputs.inf_tar_t)
        %outputTableValues = [outputTableValues,outputs.inf_tar_t]; 
        %variableNames = [variableNames,"inf_tar_t"];
        outputTable = [outputTable,table(outputs.inf_tar_t,VariableNames={'inf_tar_t'})];
    end
    if ~isempty(outputs.vir_tar_b)
        %outputTableValues = [outputTableValues,outputs.vir_tar_b]; 
        %variableNames = [variableNames,"vir_tar_b"];
        outputTable = [outputTable,table(outputs.vir_tar_b,VariableNames={'vir_tar_b'})];
    end
    if ~isempty(outputs.inf_tar_b)
        %outputTableValues = [outputTableValues,outputs.inf_tar_b]; 
        %variableNames = [variableNames,"inf_tar_b"];
        outputTable = [outputTable,table(outputs.inf_tar_b,VariableNames={'inf_tar_b'})];
    end

    %outputTable = table(outputTableValues,'VariableNames',variableNames);


if overwrite == 1
    writetable(outputTable,filenameOutputs,'WriteMode','overwritesheet',...
    'WriteRowNames',true)
else
    writetable(outputTable,filenameOutputs,'WriteMode','Append',...
    'WriteVariableNames',false,'WriteRowNames',true) 
end

end
