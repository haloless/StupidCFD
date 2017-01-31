function [] = SaveObjFile(filename, node,elem)

fileid = fopen(filename,'w');

nnode = size(node,2);
nelem = size(elem,1);

fprintf(fileid,'# sphere grid\n');

for i = 1:nnode
	fprintf(fileid, 'v %f %f %f\n', node(1,i),node(2,i),node(3,i));
end

for i = 1:nelem
	fprintf(fileid, 'f %d %d %d\n', elem(i,1),elem(i,2),elem(i,3));
end

fclose(fileid);

return
end


