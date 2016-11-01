# EM_aligner
A set of Matlab tools for aligning EM images into a coherent image volume in two and three dimensions. This library works in conjunction with the "Renderer" ecosystem of tools. 

Status: In production use at Janelia. This is a nascent set of tools that is undergoing large changes. We consider the library suitable for use by our collaborators as well as other research groups. Due to limited staffing, we do not guarantee support for outside groups.

Prerequisites 
- 	Renderer and point-match services and dependencies: (available freely and documented here: https://github.com/saalfeldlab/render)
-	(Optional) Bash scripts dependent on a system call that launches a java process
	-	Client-side rendering: This is relevant to all steps requiring generation of point-matches; full section montage and cross-layer point-match generation.
	-	Rough section alignment: This script
	-	Point-match generation (two scripts)

Main steps for stitching small-to-moderate datasets:
- 	Install Renderer and point-match services and dependencies as indicated above
-	Ingest image metadata:
	-	See rules and assumptions for tile-spec ingestion (coming soon)
-	Run montage of every section (same z-coordinate value)
-	Run rough alignment of montage collection
-	Run point-match generation across layers (point out locations in the code for custom point-match generation).
-	Run full volume solve
-	Post-stitching steps (Render images, intensity correction and CATMAID staging) will not be described here.

Steps relevant to large datasets:
- 	Install Renderer service and its dependencies as indicated above
-	Ingest image metadata:
	-	See rules and assumptions for tile-spec ingestion (coming soon)
-	Run montage of every section (same z-coordinate value)
- 	Slab definition for rough alignment
-	Run rough alignment of montage collection
-	Run point-match generation across layers
-	(Optional) fuse rough slabs
-	(Optional) redefine slabs for fine alignment
-	Run individual slab solve
-	Fuse fine slabs
-	Post-stitching steps (Render images, intensity correction and CATMAID staging) will not be described here.

