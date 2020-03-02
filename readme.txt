Summary of analysis scripts for Optical Trapping data

Data handling:
	Imstack = load_imstack( dir, treat, speed, trap, spots, rep)
		- Wrapper for bioformats file opener
		- Looks for directory 'dir' within Documents/data/OpTrap (can be changed)
		- Other input arguments are parts of filename system I use to identify the dataset
		- Returns a cell array Imstack{m}{n_frames,2} where m has always been 1
	data = load_infos(treat, fieldList)
		- Loops over files in Documents/data/OpTrap/infos to find ones starting with 'treat'
		- Loads the mat file, extracts specified fields from the data struct 
		- Returns a cell array with headed columns
	view_stack(Imstack, [crop], [pause], [repeats])
		- Creates a slideshow of the input Imstack
		- Other input arguments optional and buggy
		- Crop = [y1; y2; x1; x2] will draw a box to show a crop area on each frame. Can be a vector or a 4 x n_frames matrix
		- Pause: number of seconds to wait between frames
		- Repeats: number of times to loop over whole stack
	PostProcessDataCurrent
		- Legacy code for plotting results of PostProcessCellDeform (v1)
		- I opened it for the first time when writing this document
	
Complete analysis functions:
	PostProcessCellDeform
		- Original analysis pipeline
		- Designed to work on .avi data, adapted to work on stack loaded from .tif
		- Not very accurate
	PostProcessCellDeform_v2
		- Updated analysis pipeline
		- Different options for which modules to use
		- Returns info struct with all processed data from the image series
		- Can take input options to be passed to modules, and to choose module versions
	PostProcessCellDeformImage
		- Another version of original analysis pipeline
	PostProcessCellDeformBatch
		- Batch processing version of original analysis pipeline
	PostProcessCellDeformBatch_v2(dirIn)
		- Batch version of new pipeline
		- Loads all tif datasets in rootdir/dirIn (rootdir is specified in function)
		- Runs processing with the same processing options
		- Saves to infos folder
		- Parameters for Analysis must be specified by editing the file

Analysis modules:
	find_cell
		- Find cell using canny edge detection, gaussian filter then Hough transform
		- Not very accurate
		- Returns a struct with data for each frame
	find_cell_opt
		- Parameter optimisation version of find_cell
		- Returns number of frames where at least 1 good circle was found
	find_cell_v2
		- Find cell using gaussian filter, Laplacian edge enhancement, then Hough transform
		- Returns struct
	find_cell_v2_opt
		- Parameter optimisation version of find_cell_v2
		- Returns number of frames where at least 1 good circle was found
	iterate_find_cell
		- Iterative solver for use with find_cell
		- NB: I've not actually had any success with the iterative solvers, the result is too nonsmooth for iterative methods to work easily
	segment_cell
		- Segments cell from image using edge detection
		- Returns binary mask
		- Doesn't work very well for cells
		- Legacy code
	segment_cell_v2
		- Segments cell from image using edge detection and some further filtering
		- Has several tuneable parameters
		- Returns binary mask
	segment_cell_v2_opt
		- Parameter optimisation version of segment_cell_v2
		- Returns a metric based on regionprops
	iterate_segment_cell_v3
		- Iterative solver for use with segment_cell_v3
		- NB: I've not actually had any success with the iterative solvers, the result is too nonsmooth for iterative methods to work easily
	cell_radii
		- Find cell radii using bright halo around cell
		- Assumes cell is orientated along x or y axis
		- Not very good
	ellipseDetection
		- Fits an ellipse by trying every pair of points to be major axes, and then using Hough to find an approoriate minor axis
		- I've never used it
		- Claims to be memory intensive
	ellipse_fit_fun_final
		- Fits ellipses to binary masks
		- Does some kinda fancy stuff
		- I've never used it
	ellipse_mask
		- Creates an ellipsoidal mask when given the dimensions 

	
