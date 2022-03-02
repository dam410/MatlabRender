% This function extract from a json file the scene elements information
function [meshes,cameras,lights,renders] = decode_json_blender(file)
	text = fileread(file);
	data = jsondecode(text);
	% Initialize the arrays of struct
	meshes = [];
	cameras = [];
	lights = [];
	renders = [];
	% Treat the scene elements separately
	fields_name = fieldnames(data);
	for i=1:length(fields_name)
		data_obj = getfield(data,fields_name{i});
		if strcmp(data_obj.type,'camera')
			cameras = [cameras,data_obj];
		elseif strcmp(data_obj.type,'lamp')
			lights = [lights,data_obj];
		elseif strcmp(data_obj.type,'mesh')
			meshes = [meshes,data_obj];
		elseif strcmp(data_obj.type,'render')
			renders = [renders,data_obj];
		end
	end
end

