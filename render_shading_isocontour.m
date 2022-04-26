function [I_test,I_raw,D_source,Lj,P_camera] = render_shading_isocontour(h,l,options)
	arguments
		h (1,1) double
		l (1,1) double
		options.Surface (1,1) string = 'Plane'
		options.LightType (1,1) string = 'PLS'
		options.Scattering (1,1) string = 'Phong'
		options.SurfaceParameters (:,:) double = [0;1;0;2]
		options.Meshes struct = struct()
		options.Camera struct = struct()
		options.Light struct = struct()
		options.Renderer struct = struct()
		options.LightParameters (:,:) double = [0;0.5;9;0.6]
		options.ScatteringParameters (:,:) double = [0,1.0,1]
		options.CameraIntrinsic (3,3) double = eye(3)
		options.CameraResponseFcn (1,1) function_handle = @(x) camera_response(x)
		options.InverseCameraResponseFcn (1,1) function_handle = @(x) inv_camera_response(x)
		options.UseImageTransformation (1,1) double = 1;
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%% RENDERER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%
	default_linear = @(x) double(uint8(x.*255.075));
	if length(fieldnames(options.Renderer))
		%disp('Renderer has been given');
		h = options.Renderer.data.height;
		l = options.Renderer.data.width;
		scale = options.Renderer.data.scale*l/h;
		pixel_aspect_ratio = options.Renderer.data.pixel_aspect_ratio;
		% Here we will check that every unknown stuff in Blender renderer
		% is set to raw or default
		renderer = options.Renderer.data;
		check_display = strcmp(renderer.display_device,'None');
		check_look = strcmp(renderer.look,'None');
		check_sequencer_color = strcmp(renderer.sequencer_color,'Raw');
		check_view = strcmp(renderer.view_transform,'Default');
		if check_display*check_look*check_sequencer_color*check_view
			%disp('No special blender function');
			% See if blender use a curve
			if ~renderer.use_curve
				%disp('Curve not used');
				camResponseFcn = default_linear;
			else
				%disp('Curve in used');
				% Count the point in the Bezier curve
				% If all curve have only two points
				% and those points are [0,0] and [1,1]
				% the response is linear
				mat_linear = [0;1;0;1];
				check_linear_B = (renderer.curve_B.nb_pt == 2) && all(renderer.curve_B.pts_location(:) == mat_linear);
				check_linear_C = (renderer.curve_C.nb_pt == 2) && all(renderer.curve_C.pts_location(:) == mat_linear);
				check_linear_R = (renderer.curve_R.nb_pt == 2) && all(renderer.curve_R.pts_location(:) == mat_linear);
				check_linear_G = (renderer.curve_G.nb_pt == 2) && all(renderer.curve_G.pts_location(:) == mat_linear);
				if check_linear_B*check_linear_G*check_linear_R*check_linear_C
					%disp('All curves are linears');
					camResponseFcn = default_linear;
				else
					%disp('Not all curves are linears');
					camResponseFcn = options.CameraResponseFcn;
				end
			end
		else
			% In this case we do not really know how the
			% camera response function work so we use our own
			%disp('Some unknown blender function are applied');
			camResponseFcn = options.CameraResponseFcn;
		end
	else
		%disp('No renderer data given');
		camResponseFcn = options.CameraResponseFcn;
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMERA PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% If camera is not given, assume a generic camera with given intrinsic parameters
	if ~length(fieldnames(options.Camera))
		K = options.CameraIntrinsic;
		O = zeros(3,1);
		T_cam = eye(4);
	else
	% Otherwise get the parameters of the camera using the information given by
	%	the struct
		if length(fieldnames(options.Renderer))
			[K,O,T_cam] = get_camera_matrices(options.Camera,options.Renderer);
		else
			[K,O,T_cam] = get_camera_matrices(options.Camera);
		end
	end
	

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SURFACE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Calculate the projected points (Inversion cause of Meshgrid ... maybe)
	[u_image,v_image] = meshgrid(1:l,1:h);
	if options.UseImageTransformation
		P_camera = image_to_camera([transpose(v_image(:));transpose(u_image(:))],K,h,l);
	else
		P_camera = image_to_camera_direct([transpose(v_image(:));transpose(u_image(:))],K,h,l);
	end
	%  Intersection with the surface (Might later add new surface type)
	% P are the 3D points of intersection in camera frame
	% N_p are the normals to the surface in thoses points
	% inter are booleans indicating visibility of the points
	coef = 1;
	if strcmp(options.Surface,'Plane')
		N = options.SurfaceParameters(1:3);
		d = options.SurfaceParameters(4);
		[P,N_p,inter] = inter_rays_plane(P_camera,O,N,d);
		if strcmp(options.Scattering,'Phong')
			kd = options.ScatteringParameters(1)*ones(1,h*l);
			ks = options.ScatteringParameters(2)*ones(1,h*l);
			if length(options.ScatteringParameters)>2
				coef = options.ScatteringParameters(3);
			else
				coef = 1;
			end
		end
	elseif strcmp(options.Surface, 'Meshes')
		if strcmp(options.Scattering,'Phong')
			[P,N_p,inter,kd,ks] = inter_rays_polygons(P_camera,O,T_cam,options.Meshes);
		end
	end


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LIGHT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if ~length(fieldnames(options.Light))
		% Calculating Light intensity
		if strcmp(options.LightType,'PLS')
			S = options.LightParameters(1:3);
			Lj = options.LightParameters(4)*ones(1,size(P,2));
			D = zeros(3,1);
		elseif strcmp(options.LightType,'SLS')
			S = options.LightParameters(1:3);
			D = options.LightParameters(4:6);
			mu = options.LightParameters(8);
			S_X = (S-P)./sqrt(dot(S-P,S-P));
			Lj = options.LightParameters(7)*exp(-mu*(1+dot(S_X,D*ones(1,h*l))));
		end
	else
		% Calculating Light parameters and intensity using given struct
		if strcmp(options.Light.data.type,'point')
			inv_T_cam = inv(T_cam);
			S_4 = inv(T_cam)*options.Light.pose*[0;0;0;1];
			S = S_4(1:3)/S_4(4);
			if strcmp(options.Light.data.data.falloff_type,'CONSTANT')
				Lj = options.Light.data.energy*ones(1,size(P,2));
			elseif strcmp(options.Light.data.data.falloff_type,'INVERSE_SQUARE')
				r2_X = sum((S-P).^2);
				D_2 = options.Light.data.distance^2;
				Lj = inter.*options.Light.data.energy.*(D_2./(D_2+options.Light.data.data.quadratic_attenuation*r2_X));
			end
			%Lj = 0.5*ones(1,size(P,2));
			D = zeros(3,1);
		end
	end

	%% Calculate intensities with the shading equation
	if strcmp(options.Scattering,'Phong')
		% Vector d_p corresponding to the distance of each plane to camera
		d_p = -sum(P.*N_p,1);
		R = S-2*(d_p+dot(N_p,S*ones(1,h*l))).*N_p;
		S_X = (S-P)./sqrt(dot(S-P,S-P));
		R_X = (R-P)./sqrt(dot(R-P,R-P));
		O_X = (O-P)./sqrt(dot(O-P,O-P));
		L_P = inter.*(kd.*Lj.*abs(dot(S_X,N_p))+ks.*Lj.*power(max(0,-dot(R_X,O_X)),coef));
		% Create the image
		I_P = reshape(camResponseFcn(L_P),h,l);
		I = zeros(h,l,3);
		I(:,:,1) = I_P;
		I(:,:,2) = I_P;
		I(:,:,3) = I_P;
		%% Check that the formulation with S instead of R is equivalent
		%L_P_2 = inter.*(kd*Lj*dot(S_X,N*ones(1,h*l))+...
		%        ks*Lj*max(0,(-(dot(S_X,O_X) - 2*(dot(S_X,N*ones(1,h*l))).*dot(O_X,N*ones(1,h*l))))));
		%I_P_2 = reshape(camera_response(L_P_2),h,l);
		%diff_formulation = max(max(abs(I_P_2-I_P)))

	end
	I_test = uint8(I);
	I_raw = reshape(L_P,h,l);
	S_P = S-P;
	D_source = reshape(sqrt(sum((S_P).^2)),h,l);
	Lj = reshape(Lj,h,l);
	if strcmp(options.Surface,'Plane')
		%draw_view_3D(O,S,D,P,inter,I_test);
	elseif strcmp(options.Surface, 'Meshes')
		draw_view_3D_blender(O,S,D,options.Meshes,T_cam,P,I_test);
	end
end
