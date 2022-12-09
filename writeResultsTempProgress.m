function writeResultsTempProgress(infected,filenameResults)

results = table(infected,'VariableNames',{'infected'});

writetable(results,filenameResults,'WriteMode','Append',...
    'WriteVariableNames',false,'WriteRowNames',true) 

end
