% Stress and strain
nodes_strain = zeros(n_nodes, 7);
nodes_stress = zeros(n_nodes, 7);
for ele=1:n_elements
    % Assign average stress and strain of all connected elements to each node
    for node=1:4
        nodes_strain(connectivity(ele,node),1:6) = elements_strain(ele,:);
        nodes_strain(connectivity(ele,node),7) = nodes_strain(connectivity(ele,1),7) + 1;
        nodes_stress(connectivity(ele,node),1:6) = elements_stress(ele,:);
        nodes_stress(connectivity(ele,node),7) = nodes_stress(connectivity(ele,1),7) + 1;
    end
end
nodes_strain = nodes_strain(:,1:6)./nodes_strain(:,7);
nodes_stress = nodes_stress(:,1:6)./nodes_stress(:,7);
nodes_von_mises = zeros(n_nodes, 1);
for node=1:n_nodes
    xx = nodes_stress(node,1);
    yy = nodes_stress(node,2);
    zz = nodes_stress(node,3);
    xy = nodes_stress(node,4);
    yz = nodes_stress(node,5);
    xz = nodes_stress(node,6);
    nodes_von_mises(node,:) = sqrt(0.5*((xx-yy)^2.0+(yy-zz)^2+(zz-xx)^2+6.0*((xy^2)+(yz^2)+(xz^2))));
end

% Stress XX
figure(7)
colormap jet
trisurf(ELEMENT_TRIANGULATION.freeBoundary, 1000.0*nodes(:,1), 1000.0*nodes(:,2), 1000.0*nodes(:,3), 1.0e-6*nodes_stress(:,1), 'EdgeAlpha', '0.1', 'FaceColor', 'interp');
title("X Normal Stress",'FontSize',16)
c = colorbar('eastoutside');
c.Label.String = 'σxx (MPa)';
c.Label.FontSize = 12;
c.FontSize = 12;
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])
h=gca;
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];
h.ZAxis.TickLength = [0 0];
axis equal;
view(30,30);
set(gcf,'Position',[200 55 750 1000])
saveas(gcf,'stress_xx.png',"png")

% Stress YY
figure(8)
colormap jet
trisurf(ELEMENT_TRIANGULATION.freeBoundary, 1000.0*nodes(:,1), 1000.0*nodes(:,2), 1000.0*nodes(:,3), 1.0e-6*nodes_stress(:,2), 'EdgeAlpha', '0.1', 'FaceColor', 'interp');
title("Y Normal Stress",'FontSize',16)
c = colorbar('eastoutside');
c.Label.String = 'σyy (MPa)';
c.Label.FontSize = 12;
c.FontSize = 12;
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])
h=gca;
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];
h.ZAxis.TickLength = [0 0];
axis equal;
view(30,30);
set(gcf,'Position',[1098 55 750 1000])
saveas(gcf,'stress_yy.png',"png")

% Stress ZZ
figure(9)
colormap jet
trisurf(ELEMENT_TRIANGULATION.freeBoundary, 1000.0*nodes(:,1), 1000.0*nodes(:,2), 1000.0*nodes(:,3), 1.0e-6*nodes_stress(:,3), 'EdgeAlpha', '0.1', 'FaceColor', 'interp');
title("Z Normal Stress",'FontSize',16)
c = colorbar('eastoutside');
c.Label.String = 'σzz (MPa)';
c.Label.FontSize = 12;
c.FontSize = 12;
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])
h=gca;
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];
h.ZAxis.TickLength = [0 0];
axis equal;
view(30,30);
set(gcf,'Position',[200 55 750 1000])
saveas(gcf,'stress_zz.png',"png")

% Stress XY
figure(10)
colormap jet
trisurf(ELEMENT_TRIANGULATION.freeBoundary, 1000.0*nodes(:,1), 1000.0*nodes(:,2), 1000.0*nodes(:,3), 1.0e-6*nodes_stress(:,4), 'EdgeAlpha', '0.1', 'FaceColor', 'interp');
title("XY Shear Stress",'FontSize',16)
c = colorbar('eastoutside');
c.Label.String = 'σxy (MPa)';
c.Label.FontSize = 12;
c.FontSize = 12;
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])
h=gca;
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];
h.ZAxis.TickLength = [0 0];
axis equal;
view(30,30);
set(gcf,'Position',[1098 55 750 1000])
saveas(gcf,'stress_xy.png',"png")

% Stress YZ
figure(11)
colormap jet
trisurf(ELEMENT_TRIANGULATION.freeBoundary, 1000.0*nodes(:,1), 1000.0*nodes(:,2), 1000.0*nodes(:,3), 1.0e-6*nodes_stress(:,5), 'EdgeAlpha', '0.1', 'FaceColor', 'interp');
title("YZ Shear Stress",'FontSize',16)
c = colorbar('eastoutside');
c.Label.String = 'σyz (MPa)';
c.Label.FontSize = 12;
c.FontSize = 12;
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])
h=gca;
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];
h.ZAxis.TickLength = [0 0];
axis equal;
view(30,30);
set(gcf,'Position',[200 55 750 1000])
saveas(gcf,'stress_yz.png',"png")

% Stress XZ
figure(12)
colormap jet
trisurf(ELEMENT_TRIANGULATION.freeBoundary, 1000.0*nodes(:,1), 1000.0*nodes(:,2), 1000.0*nodes(:,3), 1.0e-6*nodes_stress(:,6), 'EdgeAlpha', '0.1', 'FaceColor', 'interp');
title("XZ Shear Stress",'FontSize',16)
c = colorbar('eastoutside');
c.Label.String = 'σxz (MPa)';
c.Label.FontSize = 12;
c.FontSize = 12;
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])
h=gca;
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];
h.ZAxis.TickLength = [0 0];
axis equal;
view(30,30);
set(gcf,'Position',[1098 55 750 1000])
saveas(gcf,'stress_xz.png',"png")

% Von Mises Stress
figure(13)
colormap jet
trisurf(ELEMENT_TRIANGULATION.freeBoundary, 1000.0*nodes(:,1), 1000.0*nodes(:,2), 1000.0*nodes(:,3), 1.0e-6*nodes_von_mises(:,1), 'EdgeAlpha', '0.1', 'FaceColor', 'interp');
title("Von Mises Stress",'FontSize',16)
c = colorbar('eastoutside');
c.Label.String = 'σv (MPa)';
c.Label.FontSize = 12;
c.FontSize = 12;
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])
h=gca;
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];
h.ZAxis.TickLength = [0 0];
axis equal;
view(30,30);
set(gcf,'Position',[200 55 750 1000])
saveas(gcf,'von_mises.png',"png")