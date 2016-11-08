##Prerequisites
We are assuming that the Renderer and point-match services (and associated database) are set up and accessible for example at http://tem-services.int.janelia.org.
Also, you are using Matlab 2015a and above with toolboxes: Computer Vision Systems (or Video and Blockset), ImageProcessing, Statistics, (optional) Matlab compiler and (optional) Parallel computing. The EM_aligner directory and subdirectories are on your Matlab path.

Additionally:
A Renderer collection of montaged (contiguous) sections exists. Separation in z between these sections should be small enough that features can be detected across sections.

## Option 1: Rough alignment using SIFT point-matches runs within a Matlab session, but calls external Java code for point-matching between montage scapes

Rough allignment across montaged sections for a specified range of z values. Calls external Java code to determine point-matches between montage sections, then solves the registration problem for the montage scapes using these point-matches, applies the resulting transformation to all tiles in the montage collection within the specified z range, and persists the resulting transformations into a "rough-aligned" Renderer collection. 

at the Matlab prompt:

[1] fn = 'full_path_to_json_configuration_file.json';  %% specify the input json file

[2] montage_section_SL_prll(fn); %% perform montage

If CATMAID dynamic rendering is set up, you can view your rough alignment using a URL for example similar to this:
http://tem-services.int.janelia.org:8080/render-ws/view/stacks.html?owner=flyTEM&project=test&dynamicRenderHost=renderer:8080&catmaidHost=renderer-catmaid:8000

Rough alignment result should be very close to the expected final quality.

An example json input file is provided below.


```json

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
		"lambda": 0.1,
		"edge_lambda": 0.1,
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

```
## Montage using deployed application: deployed point-match generation + solve
Under directory matlab_compiled, open the file "compile_montage_section" for editing.
Edit lines 4 and 5 and run the script.
Now the "compiled" output file can be called from a the shell (or used in a qsub system) by passing it the path to the json input file without need for additional Matlab licenses.
An example for deploying multiple montage jobs (using a Matlab script) is provided [here] (/template_production_scripts/SPC_scripts/script_to_montage_multiple_sections_using_deployed_montager.m).




