function [K,O,T_cam,h,l] = get_camera_matrices(cam_struct,renderer_struct)
	focal = cam_struct.data.focal;
	sx = cam_struct.data.sx;
	sy = cam_struct.data.sy;
	x0 = cam_struct.data.x0;
	y0 = cam_struct.data.y0;
	if nargin<2
		h = cam_struct.data.height;
		l = cam_struct.data.width;
		scale = cam_struct.data.scale*l/h;
		pixel_aspect_ratio = cam_struct.data.pixel_aspect_ratio;
	else
		h = renderer_struct.data.height;
		l = renderer_struct.data.width;
		scale = renderer_struct.data.scale*l/h;
		pixel_aspect_ratio = renderer_struct.data.pixel_aspect_ratio;
	end
	if strcmp(cam_struct.data.sensor_fit,'VERTICAL');
		su = 2*h/l*focal*scale/sx/pixel_aspect_ratio;
		sv = 2*h/l*focal*scale/sy;
	else
		su = 2*h/l*focal*scale/sx;
		sv = 2*h/l*focal*scale/sy;
	end
	K = [-sv, 0, y0; 0, su, x0; 0, 0, 1];
	% Even if the camera is not centered in the scene, this is our reference frame
	% the other scene elements will be transformed
	O = zeros(3,1);
	% The way we model the camera (facing negative z axis), we have to invert it
	% from the blender file
	R_cam = cam_struct.pose(1:3,1:3);
	T_cam = cam_struct.pose;
	T_cam(1:3,1:3) = R_cam*diag([1,1,1]);
	T_cam = T_cam*[1,0,0,0;0,1,0,0;0,0,-1,0;0,0,0,1];
	H_image = H_image_fcn(h,l);
end
