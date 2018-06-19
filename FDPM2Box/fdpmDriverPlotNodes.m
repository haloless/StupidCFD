%fdpmDriverPlotNodes: Plot fdpm node mesh
% It can also plot boundary nodes where BC is applied.

clear args_data args_legend;

args_data = {nodeX,nodeY,'.k'};
args_legend = {'node'};

if exist('fixed_nodes_x','var') && ~isempty(fixed_nodes_x)
    args_data{end+1} = nodeX(fixed_nodes_x);
    args_data{end+1} = nodeY(fixed_nodes_x);
    args_data{end+1} = '>';
    args_legend{end+1} = 'dispx';
end
if exist('fixed_nodes_y','var') && ~isempty(fixed_nodes_y)
    args_data{end+1} = nodeX(fixed_nodes_y);
    args_data{end+1} = nodeY(fixed_nodes_y);
    args_data{end+1} = '^';
    args_legend{end+1} = 'dispy';
end
if exist('loaded_nodes','var') && ~isempty(loaded_nodes)
    args_data{end+1} = nodeX(loaded_nodes);
    args_data{end+1} = nodeY(loaded_nodes);
    args_data{end+1} = '+';
    args_legend{end+1} = 'trac';
end

figure;
plot(args_data{:});
legend(args_legend{:});
axis equal;

clear args_data args_legend;



