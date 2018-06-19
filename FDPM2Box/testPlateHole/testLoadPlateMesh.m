
clear;

[groups,nodes] = gmshLoad('plate.msh');


figure;
plot(nodes(:,1),nodes(:,2),'.k');
hold on;
group = gmshGetGroupByName(groups,'hole');
conn = group.conn;
if size(conn,2) == 2
    for i = 1:size(conn, 1)
        plot(nodes(conn(i,:),1), nodes(conn(i,:),2));
    end
else
    triplot(conn, nodes(:,1), nodes(:,2));
end
hold off;
axis equal;




