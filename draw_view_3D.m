% This function draw the elements in 3D 
function [] = draw_view_3D(O,S,D,P_plane,inter_plane,I)
	O = [0;0;0];
	% Get a square that contains only points of the plane that are not too far from the optical center
	[h,l,~] = size(I);
	% Just to remember the size of P
	%[x,y] = meshgrid(1:l,1:h)
	%P = [x(:);y(:);d(:)] % Older version
	drawable_pts = reshape(inter_plane.*(abs(P_plane(1,:))+abs(P_plane(2,:))+abs(P_plane(3,:))<20),h,l);
	ind_valid = find(drawable_pts);
	% The idea of this function is to find the biggest square area that contains points drawable (not too far and intersection with plane) in P_plane
	% Extract the smallest column index so that for each line (row), a first "1" already appears (first column '1' appearance for each row smaller than result)
	% For now we remove this type of visualization
		%i_min = max(arrayfun(@(i) min(find(drawable_pts(i,:))),1:size(drawable_pts,1)));
		%i_max = min(arrayfun(@(i) max(find(drawable_pts(i,:))),1:size(drawable_pts,1)));
		%j_min = max(arrayfun(@(j) min(find(transpose(drawable_pts(:,j)))),i_min:i_max));
		%j_max = min(arrayfun(@(j) max(find(transpose(drawable_pts(:,j)))),i_min:i_max));
	X_1 = reshape(P_plane(1,:),h,l);
	Y_1 = reshape(P_plane(2,:),h,l);
	Z_1 = reshape(P_plane(3,:),h,l);
	% Use Delaunay triangulation to get a mesh from the drawable points
	T = delaunay(X_1(ind_valid),Y_1(ind_valid));
	figure(1);%...'Name','3D visualization of plane and source');
	hold off;
	%mesh(	X_1(j_min:j_max,i_min:i_max),...
	%	Y_1(j_min:j_max,i_min:i_max),...
	%	Z_1(j_min:j_max,i_min:i_max),...
	%	I(j_min:j_max,i_min:i_max,:));
	trimesh(T,X_1(ind_valid),Y_1(ind_valid),Z_1(ind_valid),I(ind_valid));
	hold on;
	% Plot camera center and source
	%plot3(O(1),O(2),O(3),'O');
	R = [1 0 0; 0 1 0;0 0 1];
	t = transpose(O);
	pose = rigid3d(R,t);
	cam = plotCamera('AbsolutePose',pose,'Opacity',0.4,'Color','b','Size',0.2);
	plot3(S(1),S(2),S(3),'+r');
	% Show SpotLight direction if any
	if norm(D)>0.5
		% Calculate a possible rotation matrix for the spotlight	
		V = null(D*transpose(D));
		if dot(cross(V(:,1),V(:,2)),D)>0
			R_D = [V(:,1),V(:,2),D];
		else
			R_D = [V(:,2),V(:,1),D];
		end	
		pose_D = rigid3d(transpose(R_D),transpose(S));
		spot = plotCamera('AbsolutePose',pose_D,'Opacity',0.4,'Color','r','Size',0.02);
		%plot3(S(1)+[0,D(1)],S(2)+[0,D(2)],S(3)+[0,D(3)],'-r');
	end
	axis equal;
end
