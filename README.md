# EM_aligner
A set of Matlab tools for aligning EM images into a coherent image volume in two and three dimensions. This library works in conjunction with the "Renderer" ecosystem of tools. 

Status: In production use at Janelia. This is a nascent set of tools that is undergoing large changes. We consider the library suitable for use by our collaborators as well as other research groups. Due to limited staffing, we do not guarantee support for outside groups.

Prerequisites 
- 	Renderer service and its dependencies: (available freely and documented here: https://github.com/saalfeldlab/render)
-	(Optional) Bash scripts dependent on a system call that launches a java process
	-	Client-side rendering: This is relevant to all steps requiring generation of point-matches; full section montage and cross-layer point-match generation.
	-	Rough section alignment: This script
	-	Point-match generation (two scripts)

Steps relevant to all datasets:
1-	Install Renderer service (http…..) and its dependencies
2-	Ingest image metadata (http….)
a.	Rules and assumptions for tile-spec ingestion
3-	Run montage of every section (same z-coordinate value)
4-	Run rough alignment of montage collection
5-	Run point-match generation across layers (point out locations in the code for custom point-match generation).
6-	Run grand solve
7-	Post-stitching steps (Render images, intensity correction and CATMAID staging) will not be described here.

Steps relevant to large datasets:
1-	slab definition
2-	collection fusion
3-	solver limitations

