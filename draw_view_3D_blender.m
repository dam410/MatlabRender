% This function draw the elements in 3D 
function [] = draw_view_3D_blender(O,S,D,mesh_0,T_cam,P,I)
	O = [0;0;0];
	% Get a square that contains only points of the plane that are not too far from the optical center
	[h,l,~] = size(I);
	figure('Name','3D visualization of plane and source');
	hold on;
	inv_T_cam = inv(T_cam);
	for i_mesh=1:length(mesh_0)
		mesh_ = mesh_0(i_mesh);
		T_mesh = mesh_.pose;
		for i = 1:length(mesh_.data.polygons)
			P_i = inv_T_cam*T_mesh*...
				[transpose(mesh_.data.vertices(mesh_.data.polygons(i).vertices_index+1,:));...
				ones(1,length(mesh_.data.polygons(i).vertices_index))];
			P_i = P_i./P_i(4,:);
			T_i = delaunay(P_i(1,:),P_i(2,:));
			trimesh(T_i,P_i(1,:),P_i(2,:),P_i(3,:));
		end
	plot3(P(1,:),P(2,:),P(3,:),'+b');
	% Plot camera center and source
	%plot3(O(1),O(2),O(3),'O');
	R = [1 0 0; 0 1 0;0 0 1];
	t = transpose(O);
	pose = rigid3d(R,t);
	cam = plotCamera('AbsolutePose',pose,'Opacity',0.4,'Color','b','Size',0.02);
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
	end
	axis equal;
end
