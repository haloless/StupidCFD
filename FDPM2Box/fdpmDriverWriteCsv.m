
% write header
filename = 'hoge00.csv';
fid = fopen(filename,'w');
% fprintf(fid, 'x,y,z,s11,s22,s12,s33\n');
fprintf(fid, 'x,y,z,u,v,s11,s22,s12,s33\n');
fclose(fid);

clear fid;


dlmwrite(filename, [nodeCoord',zeros(numNodes,1), nodeDisp', nodeSigma'], ...
'-append', 'precision','%0.5e', 'delimiter',',');


disp(['Saved to ', filename]);

