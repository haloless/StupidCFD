function [g] = gmshGetGroupByName(groups, name)

g = [];
for i = 1:length(groups)
    if strcmp(groups(i).name, name)
        g = groups(i);
        break;
    end
end

return
end

