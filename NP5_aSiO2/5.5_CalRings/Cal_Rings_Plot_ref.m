function [positions_Si,positions_O] = Cal_Rings_Plot_ref(Atom_type_Ring,Atom_pos_Ring,AtomTtypeNum,Psize,az,el)
n=length(Atom_type_Ring);

Atom_pos_Ring(n+1,:)=Atom_pos_Ring(1,:);
Atom_type_Ring(n+1)=Atom_type_Ring(1);


ID_O=find(Atom_type_Ring==AtomTtypeNum(1));
ID_Si=find(Atom_type_Ring==AtomTtypeNum(2));
positions_Si=Atom_pos_Ring(ID_Si,:)*1;
positions_O=Atom_pos_Ring(ID_O,:)*1;
%%
radii = [0.5,0.35]*1.5; % Radius for each atom
colors = 	[0.9, 0.3, 0.3];%[251 67 62]/255; % RGB color for each atom
colors2 = [1, 1, 0.6];%[1, 0.75, 0];%[1, 1, 0.6];%[250 250 0]/255; % RGB color for each atom
%colors3 = [0 102 102]/255; % RGB color for each atom
%colors4 = [220 111 0]/255; % RGB color for each atom

%figure('color','white');
%set(gcf,'Position',[100,100,210*2,210*2]); 
hold on;
% Draw each atom as a sphere
for i = 1:size(positions_Si, 1)
    [X, Y, Z] = sphere(150); % Sphere with finer mesh
    surf(positions_Si(i,1) + radii(1)*X, ...
         positions_Si(i,2) + radii(1)*Y, ...
         positions_Si(i,3) + radii(1)*Z, ...
         'FaceColor', colors2, 'EdgeColor', 'none','FaceAlpha',1);
end

for i = 1:size(positions_O, 1)
    [X, Y, Z] = sphere(150); % Sphere with finer mesh
    surf(positions_O(i,1) + radii(2)*X, ...
         positions_O(i,2) + radii(2)*Y, ...
         positions_O(i,3) + radii(2)*Z, ...
         'FaceColor', colors, 'EdgeColor', 'none','FaceAlpha',1.0);
end

% for i=1:n
% X=Atom_pos_Ring(i,1);
% Y=Atom_pos_Ring(i,2);
% Z=Atom_pos_Ring(i,3);
% U=Atom_pos_Ring(i+1,1)-Atom_pos_Ring(i,1);
% V=Atom_pos_Ring(i+1,2)-Atom_pos_Ring(i,2);
% W=Atom_pos_Ring(i+1,3)-Atom_pos_Ring(i,3);
% quiver3(X,Y,Z,U,V,W,0,'k','LineWidth',2);
% %quiver3(Rings(1,1),Rings(1,2),Rings(1,3),Rings(2,1)-Rings(1,1),Rings(2,2)-Rings(1,2),Rings(2,3)-Rings(1,3),0)
% end

hold on

camlight('headlight'); % Light comes from camera direction
%camlight('right'); % Add extra light from right side
lighting phong; % Use Phong lighting for smooth shiny surfaces
material shiny; % Make the surface shinier
%light('Position',[0 0 -1],'Style','infinite');  % 放在右上
%material dull      
%material metal
camva(Psize);
%camzoom(1.5);        % 放大 1.5 倍
%camdolly(0.2, 0, 0);
axis equal;
axis off;
axis tight;
%view(3);
grid off;
%[caz,cel] = view
view(az,el)

end

    