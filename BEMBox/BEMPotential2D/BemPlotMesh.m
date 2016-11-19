
figure;
hold on;
for i = 1:n
    j1 = node(1,i);
    j2 = node(2,i);
    plot([y(1,j1),y(1,j2)],[y(2,j1),y(2,j2)],'.-b');
    text(x(1,i),x(2,i),num2str(bc(:,i)));
end
plot(x(1,:),x(2,:),'.r');
% quiver(x(1,:),x(2,:),dnorm(1,:),dnorm(2,:));
plot(xfield(1,:),xfield(2,:),'.g');
hold off;
axis equal;
