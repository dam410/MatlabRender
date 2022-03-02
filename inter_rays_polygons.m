% Calculate the intersection of rays with a mesh
function [P,N_p,inter,kd,ks] = inter_rays_polygons(P_camera,O,T_cam,mesh_array)
	% Extract polygons from the mesh
	current_zbuffer = Inf*ones(1,size(P_camera,2));
	current_P = zeros(3,size(P_camera,2));
	current_N_p = zeros(3,size(P_camera,2));
	current_inter = zeros(1,size(P_camera,2));
	current_kd = zeros(1,size(P_camera,2));
	current_ks = zeros(1,size(P_camera,2));
	inv_T_cam = inv(T_cam);
	% Calculate for each polygon its intersection with the rays
	for i_mesh =1:length(mesh_array)	
		polys = mesh2polygons(mesh_array(i_mesh),T_cam);
		T_mesh = mesh_array(i_mesh).pose;
		for i=1:size(polys,1)
			N = polys{i,1}(:,3);
			d = polys{i,2};
			R_poly = polys{i,1};
			t = (d+transpose(N)*O)./(transpose(N)*O-transpose(N)*P_camera);
			% Check if the point of intersection is in front of the camera
			% We remove this condition right now
			%inter = (dot(P_camera-O,N*ones(1,size(P_camera,2)))<0);
			inter = ones(1,size(P_camera,2));
			% Calculate the point of intersection
			P = O + t.*P_camera-t.*O;
			% Check if the point of intersection is inside the polygon
			coord_pt = transpose(R_poly)*P+[0;0;d];
			%figure('Name','Debuggin');
			%plot(coord_pt(1,:),coord_pt(2,:),'+r');
			%hold on;
			%plot(polys{i,3},polys{i,4},'+g');
			inside = transpose(inpolygon(transpose(coord_pt(1,:)),transpose(coord_pt(2,:)),polys{i,3},polys{i,4}));
			% Old function to check our plane projection
			% figure('Name','Show the polygon inside');
			% plot(coord_pt(1,:),coord_pt(2,:),'+r',polys{i,3},polys{i,4},'-b')
			zbuffer = sum((P-O).^2,1);
			index_valid = find(zbuffer<current_zbuffer & inter & inside);
			current_zbuffer(index_valid) = zbuffer(index_valid);
			current_P(:,index_valid) = P(:,index_valid);
			current_N_p(:,index_valid) = N*ones(1,length(index_valid));
			current_inter(index_valid) = ones(1,length(index_valid));
			material = mesh_array(i_mesh).data.materials(mesh_array(i_mesh).data.polygons(i).material_index(1)+1);
			% Right now only deal with gray level
			current_kd(index_valid) = material.diffuse_intensity*material.diffuse_color(1);
			current_ks(index_valid) = material.specular_intensity*material.specular_color(1);
		end
	end
	%figure('Name','Show the intersection');
	%imshow(reshape(inter,3072,4608));
	P = current_P;
	N_p = current_N_p;
	inter = current_inter;
	kd = current_kd;
	ks = current_ks;
end
