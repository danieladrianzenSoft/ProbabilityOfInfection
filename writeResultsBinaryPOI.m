%function writeResultsBinaryPOI(params,infected,t_inf,filenameResults,overwrite)
function writeResultsBinaryPOI(params,isInfected,filenameResults,overwrite)

results = [params, table(isInfected,'VariableNames',{'isInfected'})];

if overwrite == 1
    writetable(results,filenameResults,'WriteMode','overwritesheet',...
    'WriteRowNames',true)
else
    writetable(results,filenameResults,'WriteMode','Append',...
    'WriteVariableNames',false,'WriteRowNames',true) 
end

end
