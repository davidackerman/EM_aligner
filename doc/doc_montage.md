# Montaging a section
We are assuming that the Renderer and point-match services (and associated database) are set up and accessible at http://tem-services.int.janelia.org.
Also, you are using Matlab 2015a and above with toolboxes: Computer Vision Systems (or Video and Blockset), ImageProcessing, Statistics and (optional) Parallel computing . The EM_aligner directory and subdirectories are on your Matlab path.

A full montage of a specific section (a set of tiles sharing the same z value) defined by "source_collection" and "section_number" will calculate point-matches using SURF features between tile pairs, persist those point-matches in a point-match database collection "target_point_match_collection", solve the registration problem using "solver_options", and persist the resulting transformations into the Renderer collection "target_collection". 

at the Matlab prompt:

[1] fn = 'full_path_to_json_configuration_file.json';  %% specify the input json file

[2] montage_section_SL_prll(fn); %% perform montage

If CATMAID dynamic rendering is set up, you can view your registered montage using a URL similar to this:
http://tem-services.int.janelia.org:8080/render-ws/view/stacks.html?owner=flyTEM&project=test&dynamicRenderHost=renderer:8080&catmaidHost=renderer-catmaid:8000



An example json input is provided below.


'''json

{
    "__comments": {
        "comment_0": "IMPORTANT: Before running check: source_collection, target_collection ...",
        "comment_1": "... target_point_match_collection, scratch and section_number",
        "comment_2": "Parameters to tweak I : solver_options.scale, solver_options.lambda",
        "comment_3": "Parameters to tweak II: solver_options.min_points, solver_options.max_points",
        "comment_4": " common values: solver_options.scale 0.3-0.5, solver_options.lambda = 0.0001 - 100",
        "comment_5": " common values: solver_options.min_points = 3, solver_options.max_points = 100"
    },
	"solver_options": {
		"min_tiles": 3,
		"degree": 1,
		"outlier_lambda": 1000,
		"solver": "backslash",
		"min_points": 10,
        "max_points": 50,
		"nbrs": 0,
		"xs_weight": 0.05,
		"stvec_flag": 0,
		"conn_comp": 1,
		"distributed": 0,
		"lambda": 0.01,
		"edge_lambda": 0.01,
		"small_region_lambda": 10,
		"small_region": 5,
		"calc_confidence": 1,
		"translation_fac": 1,
        "dthresh_factor": 1.2,
        "scale": 0.5,
        "use_peg": 1,
        "peg_weight": 0.0001,
        "peg_npoints": 5
	},
    "SURF_options": {
		"SURF_NumOctaves": 2,
		"SURF_NumScaleLevels": 3,
        "SURF_MetricThreshold": 1000,
        "SURF_MaxFeatures": 100
	},
	"source_collection": {
		"stack": "v1_acquire",
		"owner": "flyTEM",
		"project": "spc",
		"service_host": "tem-services.int.janelia.org:8080",
		"baseURL": "http://tem-services.int.janelia.org:8080/render-ws/v1",
		"verbose": 1
	},
	"target_collection": {
		"stack": "v1_montage_SURF_kk04_P2",
		"owner": "flyTEM",
		"project": "spc",
		"service_host": "tem-services.int.janelia.org:8080",
		"baseURL": "http://tem-services.int.janelia.org:8080/render-ws/v1",
		"verbose": 1,
        "versionNotes": "experiments to find optimal parameters for Allen dataset"
	},
	"target_point_match_collection": {
		"server": "http://tem-services.int.janelia.org:8080/render-ws/v1",
		"owner": "flyTEM",
		"match_collection": "spc_test"
	},
	"section_number": 3357,
    "image_filter": "true",
    "scratch": "/gpfs1/scratch/spc/matlab_work/montage/scratch",
    "renderer_client": "/groups/flyTEM/flyTEM/render/bin/render.sh",
	"verbose": 2,
    "ncpus": 8,
    "complete": 0

}

'''

