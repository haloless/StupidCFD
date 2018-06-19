function [groups,nodes] = gmshLoad(filename)

fp = fopen(filename);

% MeshFormat
begin_section(fp, 'MeshFormat');
fgetl(fp);
end_section(fp, 'MeshFormat');

%
begin_section(fp, 'PhysicalNames');
groups = readPhysicalNames(fp);
end_section(fp, 'PhysicalNames');

%
begin_section(fp, 'Nodes');
nodes = readNodes(fp);
end_section(fp, 'Nodes');

% 
begin_section(fp, 'Elements');
numElems = fscanf(fp, '%d\n', [1,1]);
for i = 1:numElems
    tline = fgetl(fp);
    items = sscanf(tline, '%d'); items = items.';
    elm_num = items(1);
    elm_type = items(2);
    elm_ntag = items(3);
    elm_phys = items(4);
    elm_geom = items(5);
    elm_conn = items(4+elm_ntag:end);
    groups(elm_phys).conn(end+1,:) = elm_conn;
end
end_section(fp, 'Elements');

%
if 1
    numNodes = size(nodes, 1);
    for i = 1:length(groups)
        tmp = zeros(numNodes,1);
        conn = groups(i).conn;
        for j = 1:size(conn,2)
            tmp(conn(:,j)) = 1;
        end
        groups(i).nodeIds = find(tmp == 1);
    end
end



%
fclose(fp);

return
end


function [] = begin_section(fp, name)
    tline = fgetl(fp);
    disp(tline);
end

function [] = end_section(fp, name)
    tline = fgetl(fp);
    disp(tline);
end

function [groups] = readPhysicalNames(fp)
    tline = fgetl(fp);
    numPhys = sscanf(tline, '%d');
    
    groups = struct;
    for i = 1:numPhys
        tline = fgetl(fp);
        items = sscanf(tline, '%d %d %s\n').';
        ndim = items(1);
        iphys = items(2);
        name = char(items(3:end)); name = name(2:end-1); % strip double quotes
        groups(i).ndim = ndim;
        groups(i).name = name;
        groups(i).conn = [];
    end  
end

function [nodes] = readNodes(fp)
    tline = fgetl(fp); 
    numNodes = sscanf(tline, '%d');
    
    nodes = fscanf(fp, '%d %f %f %f\n', [4, numNodes]);
    nodes = nodes.';
    nodes = nodes(:,2:end); % strip index
end


