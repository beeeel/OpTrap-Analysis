# OpTrap-Analysis

## Folder hierachy:

The cell defomation analysis pipeline is made of modules. Some of these require inputs from previous modules (unwrap_cell needs the cell's location from find_cell or LineMaxima, for example). The modules are within their own folder. These are run by the runners, in their own folder. Some of the modules and runners use functions from helpers, depending on the file format you are loading.

## Summary of analysis scripts for Optical Trapping data


Runners:

    PPCD_runner
        - Holds some option arguments for modules
        - Handles bulk loading, processing, and saving of analysed data
        - Also produces useful summary plots

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

Modules:

	find_cell
		- Find cell using canny edge detection, gaussian filter then Hough transform
		- Not very accurate
		- Returns a struct with data for each frame

	find_cell_v2
		- Find cell using gaussian filter, Laplacian edge enhancement, then Hough transform
		- Returns struct

    find_cell_v3
        - Contrast enhance, Hough transform. Attempt to speed up v2 (unsuccessful)

	segment_cell
		- Segments cell from image using edge detection
		- Returns binary mask
		- Doesn't work very well for cells
		- Legacy code

	segment_cell_v2
		- Segments cell from image using edge detection and some further filtering
		- Has several tuneable parameters
		- Returns binary mask

    segment_cell_v3
        - Contrast enhance, threshold, dilate to complete circle and remove small objects

    segment_cell_v4
        - Simple Otsu threshold

    segment_cell_v5
        - Contrast enhancement and ActiveContour

    segment_cell_v6
        - Not very good

    unwrap_cell_v1
        - Interpolation to create polar co-ordinates image of cell, fit equation measure ellipse.
        - Sensitive to poor centering

    unwrap_cell_v2
        - As v1 but with centering improvement option

    LineMaxima_v1
        - Find centre of cell from centre of bright ring around cell.

	cell_radii
		- Find cell radii using bright halo around cell
		- Assumes cell is orientated along x or y axis
		- Not very good

	ellipseDetection
		- Fits an ellipse by trying every pair of points to be major axes, and then using Hough to find an approoriate minor axis
		- Memory intensive - pass it a mask of the cell outline, not the whole cell

Helpers:

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
