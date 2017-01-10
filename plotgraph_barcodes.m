function [] = plotgraph_barcodes(tfs, con, filename)
nnodes = length(tfs);
fid = fopen(filename, 'w'); %should be .sif file
for i=1:nnodes
    %for j=(i+1):nnodes
    for j=1:nnodes
        if con(i,j)== -1
            fprintf(fid, '%s  %s  %s\n', tfs{i}, 'neg', tfs{j});
        elseif con(i,j) == 1
            fprintf(fid, '%s  %s  %s\n', tfs{i}, 'pos', tfs{j});
        elseif con(i,j) == 2
            fprintf(fid, '%s  %s  %s\n', tfs{i}, 'sym', tfs{j});
        elseif con(i,j) == 4
            fprintf(fid, '%s  %s  %s\n', tfs{i}, 'sym', tfs{j});
        elseif con(i,j) == 3
            fprintf(fid, '%s  %s  %s\n', tfs{i}, 'sn', tfs{j});
        end
    end
end
fclose(fid);

end