% Repeat execution of dbscan_sv_synaptosome_20230809_nomalD.m 100 times
for i = 1:100
    % Generate the output file name (MaxD_norm1, MaxD_norm2, ..., MaxD_norm100)
    outputFileName = ['MaxD_norm', num2str(i), '.xlsx'];

    % Execute dbscan_sv_synaptosome_20211008_nomalD.m with the desired output file name
    dbscan_sv_synaptosome_20230809_nomalD;

    % Save the output file with the generated name
    movefile('MaxD_norm.xlsx', outputFileName);

    % Optionally, you can clear variables to avoid conflicts in the next iteration
    clearvars -except i
end

