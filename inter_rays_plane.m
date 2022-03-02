function [P,N_p,inter] = inter_rays_plane(P_camera,O,N,d)
        % Return the points intersection of the ray with the plane
        t = (d+transpose(N)*O)./(transpose(N)*O-transpose(N)*P_camera);
        inter = (dot(P_camera-O,N*ones(1,size(P_camera,2)))<0);
        P = O + t.*P_camera-t.*O;
	N_p = repmat(N,1,length(P));
end

