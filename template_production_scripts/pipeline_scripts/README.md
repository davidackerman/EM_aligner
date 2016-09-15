# Building the binaries

Before building the binaries make sure that your matlab environment variables is set properly.

`
export matlabroot="/usr/local/matlab-2016a"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:\
$matlabroot/bin/glnxa64:\
$matlabroot/runtime/glnxa64:\
$matlabroot/sys/os/glnxa64:\
$matlabroot/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:\
$matlabroot/sys/java/jre/glnxa64/jre/lib/amd64/server:\
$matlabroot/sys/java/jre/glnxa64/jre/lib/amd64"
export XAPPLRESDIR="$matlabroot/X11/app-defaults"
export MCR_INHIBIT_CTF_LOCK=0
`

Then simply run:
     `make`

If you want to build only a specific tool you can run `make target` where target_name
is one or more from the list:

* generate\_diagnostic\_stats
* Register\_montage
* Register\_fine
* Register\_rough
* Collection\_fusion
* generate\_slab\_definitions

For example to build Register\_fine and Register\_rough you can run:

    `make Register_rough Register_fine`

# Running the compiled binaries

If you want to run the binaries directly you need to have the matlab toolbox libraries in your LD\_LIBRARY\_PATH so you need to set your environment just as if you compile these tools. The main advantage of compiling the tools is that the matlab licenses are only needed at compile time and after that you are no longer constraint by the available licenses


## Running Register_montage

### Command line:

`
bin/Register_montage start_section end_section json_options_file
`

### Parameters:
* start\_section - first section to register.
* end\_section - last section to register.
* json\_options\_file - JSON config file containing source, target, point match and solver options

e.g.
Run montage for sections 1 to 50 using the options defined in stitching_options.json

`
bin/Register_montage 1 50 stitching_options.json
`

## Running Register_rough

### Command line:

`
bin/Register_rough slabs json_slab_definitions rough_slabs_dir store_dir json_options_file
`

### Parameters:
* slabs - array of slabs to be rough aligned
* json\_slab\_definitions - JSON file that contains the slab definitions
* rough\_slabs\_dir - directory to store serialized rough slab matrices
* store\_dir - directory used for intermediate results (defaults to '/nobackup/flyTEM/khairy/FAFB00v13/montage\_scape\_pms')
* json\_options\_file - JSON config file

e.g.
Run rough alignment for slabs 1 to 10 using the stitching_options.json and the slabs defined in rough_slab_defs.json

`
bin/Register_rough "1:10" rough_slab_defs.json "" stitching_options.json
`

## Running Register_fine

### Command line:

`
bin/Register_fine slabs json_slab_definitions slabs_dir json_options_file
`

### Parameters:
* slabs - array of slabs to be rough aligned
* json\_slab\_definitions - JSON file that contains the slab definitions
* slabs\_dir - directory used for storing solved slabs (defaults to '/nobackup/flyTEM/khairy/FAFB00v13/matlab_slabs')
* json\_options\_file - JSON config file

e.g.
Run fine alignment for slabs 1 to 10 using fine_slabs_def.json and stitching_options.json

`
bin/Register_fine "1:10" fine_slabs_def.json "" stitching_options.json
`

## Running Collection_fusion

### Command line:

`
bin/Collection_fusion start_slab end_slab json_slab_definitions json_options_file
`

### Parameters:
* start_slab - first slab to fuse.
* end_slab - last slab to fuse.
* json\_slab\_definitions - JSON file that contains the slab definitions
* json\_options\_file - JSON config file

e.g.
`
bin/Collection_fusion 2 78 fine_slabs_def.json stitching_options.json
`


## Running generate\_diagnostic\_stats

### Command line:

`
bin/generate_diagnostic_stats zstart zend todo json_options_file
`

### Parameters:
* zstart - first section to be checked - defaults to 1
* zend - last section to be checked - defaults to 7062
* todo - specifies what stats to generate - defaults to 15
* json\_options\_file - JSON config file

_todo_ is a bitwise value that specifes what data to generate. The meaning of the bits is the following:

	- 0x1 - generate the area and perimeter median plots
	- 0x2 - save diagnostic data table in a CSV file
	- 0x4 - generate the residual diagrams
	- 0x8 - generate the deformation diagrams

The default value for _todo_ is 15, i.e., generate all stats that it can generate

e.g.
Generate the area and perimeter plots for the first 20 sections from 'FULL_FAFB_FUSED' and also
save the diagnostic data in a CSV file.

`
bin/generate_diagnostic_stats 1 20 3 stitching_options.json
`
