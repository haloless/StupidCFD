%

plot(nodeX,nodeY,'o', ...
nodeX(fixed_nodes),nodeY(fixed_nodes),'x', ...
nodeX(loaded_nodes),nodeY(loaded_nodes),'s');
legend('node','disp','trac');
axis('equal');
