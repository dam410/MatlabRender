% This function calculate for each polygon in the mesh an intrinsic frame on it, to test the intersection with our rays
function [polys] = mesh2polygons(mesh_,T_cam)
	polys = [];
	for i=1:length(mesh_)
		mesh_i = mesh_(i);
		% Maybe we should vectorize but for now let's just make it quick and readable
		T_mesh = mesh_i.pose;
		inv_T_cam = inv(T_cam);
		for i = 1:length(mesh_i.data.polygons)
			% Calculate the two first points and normal in world coordinates
			P1 = inv_T_cam*T_mesh*[transpose(mesh_i.data.vertices(mesh_i.data.polygons(i).vertices_index(1)+1,:));1];
			P1 = P1/P1(4);
			P2 = inv_T_cam*T_mesh*[transpose(mesh_i.data.vertices(mesh_i.data.polygons(i).vertices_index(2)+1,:));1];
			P2 = P2/P2(4);
			N = inv_T_cam*T_mesh*[mesh_i.data.polygons(i).normal;0];
			N = N/norm(N);
			% Calculate the rotation matrix to get into the plane embedding the polygon
			V1 = (P2(1:3)-P1(1:3))/norm((P2(1:3)-P1(1:3)));
			R_poly = [V1,cross(N(1:3),V1),N(1:3)];
			d_poly = -dot(P1,N);
			% Calculate the 2D coordinates of all the vertex on this plane (We suppose that they are on the plane)
			% Put all the points into camera frame
			P_polygon = inv_T_cam*T_mesh*[transpose(mesh_i.data.vertices(mesh_i.data.polygons(i).vertices_index+1,:));ones(1,length(mesh_i.data.polygons(i).vertices_index))];
			P_polygon = P_polygon(1:3,:)./P_polygon(4,:);
			%axis equal;
			coord_pt = transpose(R_poly)*P_polygon(1:3,:)+[0;0;d_poly];
			x_pt = [transpose(coord_pt(1,:));coord_pt(1,1)];
			y_pt = [transpose(coord_pt(2,:));coord_pt(2,1)];
			polys = [polys;{R_poly,d_poly,x_pt,y_pt}];

			% Old function to verify our plane projection with R_poly and d_poly
			%figure('Name','Polygon alone in space in camera frame');
			%plot3(P_polygon(1,:),P_polygon(2,:),P_polygon(3,:),'ob');
		end
	end
end
