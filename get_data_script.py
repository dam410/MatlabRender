from mathutils import *; 
from math import *; 
from numpy import * 
import json 
import itertools 
import math
import bpy

# This function gets the pose of any object and return it as a 4x4 matrix
def get_pose(obj):
	# Get the euler angle and convert to matrix
	obj.rotation_mode = 'XYZ';
	theta_x = obj.rotation_euler[0];
	theta_y = obj.rotation_euler[1];
	theta_z = obj.rotation_euler[2];
	R_x = matrix([[1, 0, 0 ], [0, cos(theta_x), -sin(theta_x)],[ 0, sin(theta_x), cos(theta_x)]]);
	R_y = matrix([[cos(theta_y), 0, sin(theta_y) ],[ 0,1,0 ],[ -sin(theta_y), 0, cos(theta_y) ]]);
	R_z = matrix([[cos(theta_z), -sin(theta_z), 0],[ sin(theta_z), cos(theta_z),0],[ 0, 0, 1]]);
	R = R_z*R_y*R_x;
	# Get the position vector
	t = matrix(obj.location).transpose();
	# Check the scale
	scale = 1.0;
	if obj.scale[0]!=obj.scale[1] or obj.scale[1]!=obj.scale[2]:
		print('Error : Scale of the object is not uniform');
	else:
		scale = obj.scale[0];
	# Create the transformation matrix
	T = concatenate((concatenate((R,t),axis=1),matrix([0,0,0,scale])),axis=0);
	return T;

def rotation_matrix(theta_x,theta_y,theta_z):
	R_x = matrix([[1, 0, 0 ], [0, cos(theta_x), -sin(theta_x)],[ 0, sin(theta_x), cos(theta_x)]]);
	R_y = matrix([[cos(theta_y), 0, sin(theta_y) ],[ 0,1,0 ],[ -sin(theta_y), 0, cos(theta_y) ]]);
	R_z = matrix([[cos(theta_z), -sin(theta_z), 0],[ sin(theta_z), cos(theta_z),0],[ 0, 0, 1]]);
	return R_z*R_y*R_x;

# Get specific properties of the light source in a dict
def get_lamp_dict(lamp):
	lamp_dict = {};
	lamp_dict['energy'] = lamp.data.energy;
	lamp_dict['distance'] = lamp.data.distance;
	if bpy.app.version < (2,80,0):
		lamp_dict['use_specular'] = lamp.data.use_specular;
		lamp_dict['use_diffuse'] = lamp.data.use_diffuse;
	else:
		lamp_dict['use_specular'] = lamp.data.specular_factor > 0.0;
		lamp_dict['use_diffuse'] = lamp.data.diffuse_factor > 0.0;
	# Get its intensity and its type
	if lamp.data.type == 'POINT':
		lamp_dict['type'] = 'point';
		lamp_point_dict = {};
		lamp_point_dict['falloff_type'] = lamp.data.falloff_type;
		lamp_point_dict['quadratic_attenuation'] = lamp.data.quadratic_attenuation;
		lamp_dict['data'] = lamp_point_dict;
	elif lamp.data.type == 'SPOT':
		lamp_dict['type'] = 'spot';
		lamp_spot_dict = {};
		lamp_spot_dict['falloff_type'] = lamp.data.falloff_type;
		lamp_spot_dict['quadratic_attenuation'] = lamp.data.quadratic_attenuation;
		lamp_spot_dict['spot_size'] = lamp.data.spot_size;
		lamp_spot_dict['spot_blend'] = lamp.data.spot_blend;
		lamp_dict['data'] = lamp_spot_dict;
	else:
		print('Error : Type of lamp not recognized');
	return lamp_dict;

# Get the specific properties of a camera in a dict
def get_cam_dict(cam):
	cam_dict = {};
	cam_dict['focal'] = cam.data.lens;
	cam_dict['type'] = cam.data.type;
	cam_dict['sensor_fit'] = cam.data.sensor_fit;
	cam_dict['sx'] = cam.data.sensor_height;
	cam_dict['sy'] = cam.data.sensor_width;
	cam_dict['x0'] = cam.data.shift_x;
	cam_dict['y0'] = cam.data.shift_y;
	#scene = bpy.context.scene
	#cam_dict['height'] = scene.render.resolution_y;
	#cam_dict['width'] = scene.render.resolution_x;
	#cam_dict['scale'] = scene.render.resolution_percentage / 100.0;
	#cam_dict['pixel_aspect_ratio'] = scene.render.pixel_aspect_x/scene.render.pixel_aspect_y;
	return cam_dict;

# Get the specific properties of a mesh in a dict
def get_mesh_dict(mesh):
	mesh_dict = {};
	# Accumulate the vertices and their normals in a list
	vertices = [];
	normals = [];
	for i in range(0,size(mesh.data.vertices)):
		vertices.append(mesh.data.vertices[i].co);
		normals.append(mesh.data.vertices[i].normal);
	vertices_mat = matrix(vertices);
	normals_mat = matrix(normals);
	mesh_dict['vertices'] = vertices_mat.tolist();
	mesh_dict['normals'] = normals_mat.tolist();
	# Adding the polygons (a dict containing: list of vertex + a material + a normal)
	#	Add a list of vertex indices and a list of materials
	polygons = [];
	for i in range(0,size(mesh.data.polygons)):
		polygon_dict = {};
		polygons_vertices_index = [];
		polygon_dict['material_index'] = mesh.data.polygons[i].material_index;
		polygon_dict['normal'] = list(mesh.data.polygons[i].normal.to_tuple());
		for j in mesh.data.polygons[i].loop_indices:
			polygons_vertices_index.append(mesh.data.loops[j].vertex_index);
		polygon_dict['vertices_index'] = polygons_vertices_index;
		polygons.append(polygon_dict);
	mesh_dict['polygons'] = polygons;
	# Adding the materials (a dict containing properties)
	materials = [];
	for i in range(0,size(mesh.data.materials)):
		material_dict = {};
		material_dict['diffuse_color'] = list(mesh.data.materials[i].diffuse_color);
		material_dict['specular_color'] = list(mesh.data.materials[i].specular_color);
		if bpy.app.version > (2,80,0):
			d_col = mesh.data.materials[i].diffuse_color;
			material_dict['diffuse_intensity'] = (d_col[0]+d_col[1]+d_col[2])/3.0;
			material_dict['specular_shader'] = 'COOKTORR';
			material_dict['diffuse_shader'] = 'LAMBERT'
			material_dict['ambient'] = 0.0
		else:
			material_dict['diffuse_intensity'] = mesh.data.materials[i].diffuse_intensity;
			material_dict['specular_shader'] = mesh.data.materials[i].specular_shader;
			material_dict['diffuse_shader'] = mesh.data.materials[i].diffuse_shader;
			material_dict['ambient'] = mesh.data.materials[i].ambient;
		material_dict['specular_intensity'] = mesh.data.materials[i].specular_intensity;
		materials.append(material_dict);
	mesh_dict['materials'] = materials;
	return mesh_dict;

# This function transform a scene element into a dict python object
def get_obj_dict(obj):
	obj_dict = {};
	# Get its pose and add it to the dict
	obj_dict['pose'] = get_pose(obj).tolist();
	# Check the type
	if obj.type == 'LAMP' or obj.type == 'LIGHT':
		obj_dict['type'] = 'lamp';
		obj_dict['data'] = get_lamp_dict(obj);
	elif obj.type == 'CAMERA':
		obj_dict['type'] = 'camera';
		obj_dict['data'] = get_cam_dict(obj);
	elif obj.type == 'MESH':
		obj_dict['type'] = 'mesh';
		obj_dict['data'] = get_mesh_dict(obj);
	else:
		print('Error : One object is not recognized, it will be ignored')
	return obj_dict;

# New object to represent the general property of the scene,
# Image resolution and color management
def get_render_property(scene):
	upper_obj_dict = {};	
	upper_obj_dict['type'] = 'render';
	obj_dict = {};
	obj_dict['height'] = scene.render.resolution_y;
	obj_dict['width'] = scene.render.resolution_x;
	obj_dict['scale'] = scene.render.resolution_percentage / 100.0;
	obj_dict['filepath'] = scene.render.filepath;
	obj_dict['pixel_aspect_ratio'] = scene.render.pixel_aspect_x/scene.render.pixel_aspect_y;
	obj_dict['display_device'] = scene.display_settings.display_device;
	obj_dict['view_transform'] = scene.view_settings.view_transform;
	obj_dict['look'] = scene.view_settings.look;
	if hasattr(scene.node_tree,'nodes') and scene.node_tree.nodes.find("NoiseLevel"):
		obj_dict['noise_level'] = scene.node_tree.nodes["NoiseLevel"].outputs[0].default_value;
	else:
		obj_dict['noise_level'] = 0;
	# TO DO: Get the sequencer color when it will work, right now just ensure that
	# it is set on raw manually
	obj_dict['sequencer_color'] = scene.sequencer_colorspace_settings.name;
	obj_dict['exposure'] = scene.view_settings.exposure;
	obj_dict['gamma'] = scene.view_settings.gamma;
	obj_dict['use_curve'] = scene.view_settings.use_curve_mapping;
	if obj_dict['use_curve']:
		obj_dict['curve_R'] = get_curvemap(scene.view_settings.curve_mapping.curves[0]);
		obj_dict['curve_G'] = get_curvemap(scene.view_settings.curve_mapping.curves[1]);
		obj_dict['curve_B'] = get_curvemap(scene.view_settings.curve_mapping.curves[2]);
		obj_dict['curve_C'] = get_curvemap(scene.view_settings.curve_mapping.curves[3]);
	upper_obj_dict['data'] = obj_dict;
	return upper_obj_dict;

def get_curvemap(curve):
	curve_dict = {};
	pts_location = [];
	pts_type = [];
	for i in range(0,size(curve.points)):
		pts_location.append(list(curve.points[i].location));
		pts_type.append(curve.points[i].handle_type);
	curve_dict['pts_location'] = pts_location;
	curve_dict['pts_type'] = pts_type;
	curve_dict['nb_pt'] = size(curve.points);
	return curve_dict;

# Generate a JSON file to export the scene elements in Matlab
def gen_file_scene(bpy,json_output,png_output,folder):
	D = bpy.data;
	D.scenes[0].render.filepath = folder+'/'+png_output;
	f = open(folder+'/'+json_output,'w');
	scene_dict = {};
	nb_obj_added = 0;
	for i in range(0,size(D.objects)):
		if not D.objects[i].hide_render:
			obj_dict = get_obj_dict(D.objects[i]);
			scene_dict[i] = obj_dict;
			nb_obj_added = nb_obj_added + 1;
	scene_dict[nb_obj_added] = get_render_property(D.scenes[0]);
	json_string = json.dumps(scene_dict, sort_keys=True, indent=4);
	f.write(json_string);
	f.close();
	bpy.ops.render.render(write_still=True);	

def generate_data_one_plane_auto(bpy,folder):
	D = bpy.data;
	plane = D.objects.get('Plane');
	cam = D.objects.get('Camera');
	source = D.objects.get('Lamp');
	cam.location = Vector([0,0,4.0]);
	cam.rotation_euler = Euler([0.0,0.0,0.0],'XYZ');
	# Range of variation choosed for our first synthetic experiment
	list_val_a = [0,15,30];
	list_val_b = [0,15,30];
	list_val_h = [0.5,1.0,1.5];
	list_val_x = [-1,0,1];
	list_val_y = [-0.5,0,0.5];
	for v_vectors in itertools.product(list_val_a,list_val_b,list_val_h,list_val_x,list_val_y):
		val_a = v_vectors[0];
		val_b = v_vectors[1];
		val_h = v_vectors[2];
		val_x = v_vectors[3];
		val_y = v_vectors[4];
		plane.location = Vector([val_x,val_y,0.0]);
		plane.rotation_euler = Euler([val_a*math.pi/180.0,val_b*math.pi/180.0,0],'XYZ');
		N = Vector(get_pose(plane)[range(0,3),:][:,2]);
		source.location = val_h*N+plane.location;
		string_a = 'a' + "{:0>+2d}".format(int(val_a));
		string_b = 'b' + "{:0>+2d}".format(int(val_b));
		string_h = 'h' + "{:0>+2d}".format(int(val_h*10));
		string_x = 'x' + "{:0>+2d}".format(int(val_x*10));
		string_y = 'y' + "{:0>+2d}".format(int(val_y*10));
		generated_filename = 'scene_'+string_a+string_b+string_h+string_x+string_y;
		gen_file_scene(bpy,generated_filename+'.json',generated_filename+'.png',folder);
	
def generate_data_two_planes_auto(bpy,folder):
	D = bpy.data;
	plane1 = D.objects.get('Plane1');
	plane2 = D.objects.get('Plane2');
	cam = D.objects.get('Camera');
	source = D.objects.get('Lamp');
	cam.location = Vector([0,0,4.0]);
	cam.rotation_euler = Euler([0.0,0.0,0.0],'XYZ');
	source.location = Vector([0,0,1]);
	# Range of variation choosed for our first synthetic experiment
	list_val_h1 = [1.0,1.5];
	#list_val_h1 = [0.5,1.0,1.5];
	list_val_h2 = [1.0,1.5];
	#list_val_h2 = [0.5,1.0,1.5];
	list_val_a1 = [-20,0,20];
	#list_val_a1 = [10,30,40];
	list_val_a2 = [-20,0,20]
	#list_val_a2 = [-10,-30,-40]
	list_val_b1 = [-35,-55];
	#list_val_b1 = [-30,-15,0,15,30];
	list_val_b2 = [35,55];
	#list_val_b2 = [-30,-15,0,15,30];
	for v_vectors in itertools.product(list_val_h1,list_val_h2,list_val_a1,list_val_a2,list_val_b1,list_val_b2):
		val_h1 = v_vectors[0];
		val_h2 = v_vectors[1];
		val_a1 = v_vectors[2];
		val_a2 = v_vectors[3];
		val_b1 = v_vectors[4];
		val_b2 = v_vectors[5];
		# Rotate the plane then translate it to be at the specified distance from the source
		plane1.rotation_euler = Euler([val_a1*math.pi/180.0,val_b1*math.pi/180.0,0],'XYZ');
		N1 = Vector(get_pose(plane1)[range(0,3),:][:,2]);
		plane1.location = source.location-val_h1*N1;
		# Rotate the plane then translate it to be at the specified distance from the source
		plane2.rotation_euler = Euler([val_a2*math.pi/180.0,val_b2*math.pi/180.0,0],'XYZ');
		N2 = Vector(get_pose(plane2)[range(0,3),:][:,2]);
		plane2.location = source.location-val_h2*N2;
		# Naming the file
		string_a1 = '_a1_' + "{:0>+2d}".format(int(val_a1));
		string_a2 = '_a2_' + "{:0>+2d}".format(int(val_a2));
		string_b1 = '_b1_' + "{:0>+2d}".format(int(val_b1));
		string_b2 = '_b2_' + "{:0>+2d}".format(int(val_b2));
		string_h1 = '_h1_' + "{:0>+2d}".format(int(val_h1*10));
		string_h2 = '_h2_' + "{:0>+2d}".format(int(val_h2*10));
		generated_filename = 'scene'+string_a1+string_a2+string_b1+string_b2+string_h1+string_h2;
		gen_file_scene(bpy,generated_filename+'.json',generated_filename+'.png',folder);


def rotation_euler(R):
	theta_x = math.atan2(R[2,1],R[2,2]);
	theta_y = math.atan2(-R[2,0],math.sqrt(R[2,1]*R[2,1]+R[2,2]*R[2,2]));
	theta_z = math.atan2(R[1,0],R[0,0]);
	return Vector([theta_x,theta_y,theta_z]);

# Place the four planes all passing through the point pt_vec
# 	alpha is the angle between planes
#	beta is the angle between planes axis and principal axis
def place_plane(bpy,alpha,beta,source,distance_source):
	D = bpy.data;
	plane1 = D.objects.get('Plane1');
	plane2 = D.objects.get('Plane2');
	plane3 = D.objects.get('Plane3');
	plane4 = D.objects.get('Plane4');
	#pt_vec = Vector([0,0,-1]);
	pt_vec = Vector([random.rand()*0.2-0.1,random.rand()*0.2-0.1,-1]);
	# Adding noise and placing the central point
	# dx = cos(math.pi/4)*sin(alpha/2)*math.sqrt(2)*0.5;
	dX = sin(alpha/2)*math.sqrt(2)*0.5;
	dz = math.cos(alpha/2)*math.sqrt(2)*0.5;
	# Place all the centers of the planes
	plane1.location = Vector([0,dX,dz]);
	plane2.location = Vector([dX,0,dz]);
	plane3.location = Vector([-dX,0,dz]);
	plane4.location = Vector([0,-dX,dz]);
	#plane1.location = Vector([-dx,dx,dz]);
	#plane2.location = Vector([dx,dx,dz]);
	#plane3.location = Vector([-dx,-dx,dz]);
	#plane4.location = Vector([dx,-dx,dz]);
	# Rotate them arround pt_vec (local origin here)
	R_b = rotation_matrix(0,beta,0);
	plane1.location = Vector(R_b*transpose(matrix(plane1.location)));
	plane2.location = Vector(R_b*transpose(matrix(plane2.location)));
	plane3.location = Vector(R_b*transpose(matrix(plane3.location)));
	plane4.location = Vector(R_b*transpose(matrix(plane4.location)));
	# Translate them at pt_vec
	plane1.location = pt_vec + plane1.location;
	plane2.location = pt_vec + plane2.location;
	plane3.location = pt_vec + plane3.location;
	plane4.location = pt_vec + plane4.location;
	R_0 = rotation_matrix(0,0,math.pi/4);
	R_alpha = rotation_matrix(pi/2-alpha/2,0,0);
	R_1 = rotation_matrix(0,0,math.pi/2);
	# R_1 = rotation_matrix(0,0,math.pi/4);
	# R_2 = rotation_matrix(0,0,-math.pi/4);
	# Rotate the planes a bit less
	#plane1.rotation_euler = rotation_euler(R_b*R_1*R_alpha*R_0);
	#plane2.rotation_euler = rotation_euler(R_b*R_2*R_alpha*R_0);
	#plane3.rotation_euler = rotation_euler(R_b*R_1*R_1*R_1*R_alpha*R_0);
	#plane4.rotation_euler = rotation_euler(R_b*R_1*R_1*R_1*R_1*R_1*R_alpha*R_0);
	plane1.rotation_euler = rotation_euler(R_b*R_alpha*R_0);
	plane2.rotation_euler = rotation_euler(R_b*R_1*R_1*R_1*R_alpha*R_0);
	plane3.rotation_euler = rotation_euler(R_b*R_1*R_alpha*R_0);
	plane4.rotation_euler = rotation_euler(R_b*R_1*R_1*R_alpha*R_0);

	# Source is placed at distance source
	source.location = pt_vec + Vector(R_b*transpose(matrix([0,0,distance_source])));


def angle_variation(bpy,nb_angle,nb_samples,folder):
	# Fecthing scene elements
	D = bpy.data;
	cam = D.objects.get('Camera');
	source = D.objects.get('Lamp');
	plane1 = D.objects.get('Plane1');
	plane2 = D.objects.get('Plane2');
	plane3 = D.objects.get('Plane3');
	plane4 = D.objects.get('Plane4');
	# Default parameters
	cam.location = Vector([0,0,4.0]);
	beta_def = math.pi/180.0*30;
	distance_source = 1.0;
	cam.rotation_euler = Euler([0.0,0.0,0.0],'XYZ');
	source.location = Vector([0,0,1]);
	plane1.hide_render = False;
	plane2.hide_render = True;
	plane3.hide_render = True;
	plane4.hide_render = False;
	f = open(folder+"/filenames_params_angle.txt", "a")
	name_gen = "scene_angles";
	for i in range(0,nb_angle):
		angle = round(20+145*random.rand(),2);
		for j in range(0,nb_samples):
			alpha = angle*math.pi/180.0;
			place_plane(bpy,alpha,beta_def,source,distance_source);
			i_str = "_"+"{:0>2d}".format(int(i));
			j_str = "_"+"{:0>2d}".format(int(j));
			generated_filename = name_gen + i_str + j_str
			gen_file_scene(bpy,generated_filename+'.json',generated_filename+'.png',folder);
			f.write(generated_filename+".json"+" "+"{:.2f}".format(angle)+"\n");
	f.close()


def nb_plane_variation(bpy,nb_samples,folder):
	# Fetching scene elements
	D = bpy.data;
	cam = D.objects.get('Camera');
	source = D.objects.get('Lamp');
	plane1 = D.objects.get('Plane1');
	plane2 = D.objects.get('Plane2');
	plane3 = D.objects.get('Plane3');
	plane4 = D.objects.get('Plane4');
	# Default parameters
	beta_def = math.pi/180.0*30;
	distance_source = 1.0;
	angle = 90;
	cam.location = Vector([0,0,4.0]);
	cam.rotation_euler = Euler([0.0,0.0,0.0],'XYZ');
	plane1.hide_render = False;
	plane2.hide_render = True;
	plane3.hide_render = True;
	plane4.hide_render = False;
	f = open(folder+"/filenames_params_nb_planes.txt", "a")
	name_gen = "scene_exp_1_n_plane";
	for i in range(0,3):
		if i==1:
			plane2.hide_render = False;
		elif i==2:
			plane3.hide_render = False;
		for j in range(0,nb_samples):
			place_plane(bpy,angle*math.pi/180.0,beta_def,source,distance_source)
			i_str = "_"+"{:0>2d}".format(int(i));
			j_str = "_"+"{:0>2d}".format(int(j));
			generated_filename = name_gen + i_str + j_str
			gen_file_scene(bpy,generated_filename+'.json',generated_filename+'.png',folder);
			f.write(generated_filename+".json"+" "+"{:.2f}".format(i+2)+"\n");
	f.close()

def distance_source_variation(bpy,nb_source,nb_samples,folder):
	# Fetching scene elements
	D = bpy.data;
	cam = D.objects.get('Camera');
	source = D.objects.get('Lamp');
	plane1 = D.objects.get('Plane1');
	plane2 = D.objects.get('Plane2');
	plane3 = D.objects.get('Plane3');
	plane4 = D.objects.get('Plane4');
	# Default parameters
	beta_def = math.pi/180.0*30;
	angle = 90;
	cam.location = Vector([0,0,4.0]);
	cam.rotation_euler = Euler([0.0,0.0,0.0],'XYZ');
	plane1.hide_render = False;
	plane2.hide_render = True;
	plane3.hide_render = True;
	plane4.hide_render = False;
	f = open(folder+"/filenames_params_distance_source.txt", "a")
	name_gen = "scene_distance_source";
	for i in range(0,nb_source):
		distance_source = random.rand()+0.5;
		for j in range(0,nb_samples):
			place_plane(bpy,angle*math.pi/180.0,beta_def,source,distance_source)
			i_str = "_"+"{:0>2d}".format(int(i));
			j_str = "_"+"{:0>2d}".format(int(j));
			generated_filename = name_gen + i_str + j_str
			gen_file_scene(bpy,generated_filename+'.json',generated_filename+'.png',folder);
			f.write(generated_filename+".json"+" "+"{:.2f}".format(distance_source)+"\n");
	f.close()

def noise_level_variation(bpy,level_max,nb_samples,folder):
	# Fetching scene elements
	D = bpy.data;
	cam = D.objects.get('Camera');
	source = D.objects.get('Lamp');
	cam.location = Vector([0,0,4.0]);
	cam.rotation_euler = Euler([0.0,0.0,0.0],'XYZ');
	source.location = Vector([0,0,1]);
	# Default parameters
	distance_source = 1.0;
	beta_def = math.pi/180.0*30;
	angle = 90;
	plane1 = D.objects.get('Plane1');
	plane2 = D.objects.get('Plane2');
	plane3 = D.objects.get('Plane3');
	plane4 = D.objects.get('Plane4');
	plane1.hide_render = False;
	plane2.hide_render = True;
	plane3.hide_render = True;
	plane4.hide_render = False;
	f = open(folder+"/filenames_params_noise_level.txt", "a")
	name_gen = "scene_noise_level";
	for i in range(0,level_max):
		for j in range(0,nb_samples):
			place_plane(bpy,angle*math.pi/180.0,beta_def,source,distance_source)
			i_str = "_"+"{:0>2d}".format(int(i));
			j_str = "_"+"{:0>2d}".format(int(j));
			generated_filename = name_gen + i_str + j_str
			gen_file_scene(bpy,generated_filename+'.json',generated_filename+'.png',folder);
			f.write(generated_filename+".json"+" "+"{:.2f}".format(i)+"\n");
	f.close()

def create_realdata_problem_replicate(bpy,folder):
	# Fetching scene elements
	D = bpy.data;
	cam = D.objects.get('Camera');
	source = D.objects.get('Lamp');
	cam.location = Vector([0,0,4.0]);
	cam.rotation_euler = Euler([0.0,0.0,0.0],'XYZ');
	source.location = Vector([0,0,1]);
	# Default parameters
	distance_source = 1.0;
	beta_def = math.pi/180.0*30;
	angle = 90;
	plane1 = D.objects.get('Plane1');
	plane2 = D.objects.get('Plane2');
	plane3 = D.objects.get('Plane3');
	plane4 = D.objects.get('Plane4');
	plane1.hide_render = False;
	plane2.hide_render = False;
	plane3.hide_render = False;
	plane4.hide_render = False;
	name_gen = "scene_exp_3_realdata_replicate";
	place_plane(bpy,angle*math.pi/180.0,beta_def,source,distance_source)
	original_source_location = source.location;
	for i_source_dx in range(0,10):
		for i_noise_level in range(0,10):
			i_str = "_source_dx_"+"{:0>2d}".format(int(i_source_dx));
			j_str = "_noise_level_"+"{:0>2d}".format(int(i_noise_level));
			generated_filename = name_gen + i_str + j_str
			source.location = original_source_location + Vector([i_source_dx/100.0,0.0,0.0]);
			bpy.context.scene.node_tree.nodes["NoiseLevel"].outputs[0].default_value = i_noise_level/20.0;
			gen_file_scene(bpy,generated_filename+'.json',generated_filename+'.png',folder);


# To use execute those lines in blender
# sys.path.append('/home/dam/Documents/PostDoc_Damien/LightModel/blender/')
# import get_data_script
# get_data_script.gen_file_scene(D,'/home/dam/Documents/PostDoc_Damien/LightModel/matlab/test_X.json');
