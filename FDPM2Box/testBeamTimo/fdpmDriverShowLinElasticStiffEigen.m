
if 1
    if ~exist('plot_disp_scale','var')
        plot_disp_scale = 50.0;
    end
    
    figure;
    neig = 10;
    [v,d] = eigs(Ktan, neig, 'sm');
    for ieig = 1:neig
        vv = v(:,ieig);
        vv = reshape(vv, 2, []).';
        triplot(tri, nodeX+plot_disp_scale.*vv(:,1), nodeY+plot_disp_scale.*vv(:,2));
        title(['eig(',int2str(ieig),')=',num2str(d(ieig,ieig))]);
        axis equal;
        pause;
    end
end


