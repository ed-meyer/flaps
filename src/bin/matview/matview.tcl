#
# $Id: matview.tcl,v 1.93 1999/03/11 22:08:17 kohl Exp $
#

#
#                       MatView Version 1.0:
#                  Simple Scalable Matrix Viewer
#           Oak Ridge National Laboratory, Oak Ridge TN.
#                    Author:  James Arthur Kohl
#                   (C) 1998 All Rights Reserved
#
#                               NOTICE
#
#  Permission to use, copy, modify, and distribute this software and
#  its documentation for any purpose and without fee is hereby granted
#  provided that the above copyright notice appear in all copies and
#  that both the copyright notice and this permission notice appear
#  in supporting documentation.
#
#  Neither the Institution, Oak Ridge National Laboratory, nor the
#  Authors make any representations about the suitability of this
#  software for any purpose.  This software is provided ``as is''
#  without express or implied warranty.
#
#  MatView was funded by the U.S. Department of Energy.
#

#
# Main TCL Source File
#

#
# Hide Main Window
#

wm withdraw .

#
# Setup Main TCL Directories
#

set ckenv [ info exists env(MATVIEW_ROOT) ]

if { $ckenv == 1 } \
	{ set matview_directory "$env(MATVIEW_ROOT)" } \
\
else \
{
	set matview_directory [ get_matview_default_dir ]

#	puts "Warning: the MATVIEW_ROOT environment variable is not set!"
#	puts "Using \"$matview_directory\"."
#	puts "(Consult the README file for help.)"
}

#
# Source Main TCL Procs
#

set ckdir_procs [ file exists procs.tcl ]

if { $ckdir_procs == 1 } \
	{ source procs.tcl } \
\
else \
{
	set ckdir_procs [ file exists $matview_directory/procs.tcl ]

	if { $ckdir_procs == 1 } \
		{ source $matview_directory/procs.tcl } \
\
	else \
	{
		puts "Error: Matview 'procs.tcl' File Not Found!"
		puts "MATVIEW_ROOT not set correctly?"
		puts "Bailing Out..."

		exit
	}
}

#
# Source TCL Utility Procs
#

set ckdir_util [ file exists util.tcl ]

if { $ckdir_util == 1 } \
	{ source util.tcl } \
\
else \
{
	set ckdir_util [ file exists $matview_directory/util.tcl ]

	if { $ckdir_util == 1 } \
		{ source $matview_directory/util.tcl } \
\
	else \
	{
		puts "Error: Matview 'util.tcl' File Not Found!"
		puts "MATVIEW_ROOT not set correctly?"
		puts "Bailing Out..."

		exit
	}
}

#
# Startup Message
#

# puts ""

# puts -nonewline "Initializing MatView [ title_info ]..."
flush stdout

#
# Proc Debug Flag
#

set proc_debug FALSE
#set proc_debug TRUE

#
# Set Main Global Names
#

set MV .matview

set MV_M $MV.map

set MV_C $MV.canvas

set MV_I $MV.image

set MV_N $MV.navigate

set MV_G $MV.graph

set MV_P $MV.progress

#
# Thumb Globals
#

set MV_TC $MV.thumb_canvas

set MV_TI $MV.thumb_image

#
# Create Toplevel Window
#

set main_height 750
set main_width 650

toplevel $MV

wm withdraw $MV

set_geometry -1 -1 125 100

set min_main_height 380
set min_main_width 650

wm minsize $MV $min_main_width $min_main_height

set depth [ winfo depth $MV ]

set save_main_height -1

#
# Set up Colors
#

# eem puts -nonewline "."
# eem flush stdout

if { $depth == 1 } \
{
	set active_fg_color white

	set selector_color black

	set fg_color black

	set bg_color white
} \
\
else \
{
	set active_fg_color purple

	set selector_color purple

	set fg_color blue

	set bg_color #d9d9d9
}

#
# Reset Main Window Background
#

$MV configure -background $bg_color

#
# Set Default Data Font (Fixed Width)
#

set data_font "-*-courier-bold-r-normal--14-*-*-*-*-*-*-*"

#
# Create Dummy Frame for Text Sizing (see get_real_text_size{})
#

# eem puts -nonewline "."
# eem flush stdout

frame $MV.dummy

get_base_text_size

#
# Set Up Other Spacing Constants
#

set border_space 10

set frame_border 2

#
# Create Title & Version Number
#

# eem puts -nonewline "."
# eem flush stdout

set title_str "M a t V i e w   [ title_info ]"

label $MV.title -text $title_str \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.title "Label"

place $MV.title -relx 0.5 -y $border_space -anchor n

update

#
# Get Main Font for Entry Widgets
#

# eem puts -nonewline "."
# eem flush stdout

set main_font [ $MV.title cget -font ]

#
# Create Main Panel
#

# eem puts -nonewline "."
# eem flush stdout

label $MV.message -text {Status:   Welcome to MatView} \
	-foreground $fg_color -background $bg_color \
	-font $main_font -relief sunken -anchor nw

restrict_bindings $MV.message "Label"

label $MV.help_msg -text {Help:} \
	-foreground $fg_color -background $bg_color \
	-font $main_font -relief sunken -anchor nw

restrict_bindings $MV.help_msg "Label"

set msg_list ""

#
# Main Menu
#

# eem puts -nonewline "."
# eem flush stdout

set file_menu_list \
	[ list \
		[ list "Load W/O Values"	statebutton	noload_values	OFF \
			noload_mat_values \
	"Toggle Loading of Matrix Data Values When Reading Matrix File." ] \
		[ list "Open Matrix File..."	command		open \
			"Load Matrix File Into MatView." ] \
		[ list "Reload Matrix File"	command		reload \
			"Reload Same Matrix File Into MatView." ] \
		[ list "Print Image..."		command		print \
			"Open Print Dialog to Print or Save Matrix Views." ] \
		[ list "About MatView..."	command		about \
			"Display MatView Overview & Author Contact Information." ] \
		[ list "Quit MatView"		command		quit \
			"Quit MatView." ] \
	]

makeMenu $MV.file_menu fileHandle $file_menu_list lower {} none

set file_lower_subs \
	[ list $MV.mode_menu $MV.reduce_menu $MV.options_menu \
		$MV.complex_menu ]

set cmd [ list raiseMenu $MV.file_menu $MV.file $MV $file_lower_subs ]

button $MV.file -text "File..." -padx 1 -pady 1 \
	-font $main_font \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color \
	-command $cmd

restrict_bindings $MV.file ""

#
# Visualization Modes Menu
#

# eem puts -nonewline "."
# eem flush stdout

set mode_menu_list \
	[ list \
		[ list "Data Modes:"	label	datamodes ] \
		[ list "Matrix Mode"	radiobutton		matrix		data_mode \
			"Set Data Mode to Matrix Semantics." ] \
		[ list "Field Mode"		radiobutton		field		data_mode \
			"Set Data Mode to Field Semantics." ] \
		[ list "Visualization Modes:"	label	vismodes ] \
		[ list "Normal Mode"	radiobutton		normal		vis_mode \
			"Set Visualization Mode to Normal." ] \
		[ list "Differential"	radiobutton		diff		vis_mode \
			"Set Visualization Mode to Differential." ] \
		[ list "Differential Mode Ops:"	label	diffops ] \
		[ list "Save Current View"			command		set_current \
			"Save Current View Field for Differential Comparison." ] \
		[ list "Free Saved View"			command		free_current \
			"Free Currently Saved View Field." ] \
		[ list "Write Saved View to File"	command		save_current \
			"Write Currently Saved Field to a File." ] \
		[ list "Load View from File"		command		load_current \
			"Load Field from a File for Differential Comparison." ] \
	]

set data_mode [ get_data_mode ]

set vis_mode [ get_vis_mode ]

makeMenu $MV.mode_menu modeHandle $mode_menu_list lower {} none

set mode_lower_subs \
	[ list $MV.file_menu $MV.reduce_menu $MV.options_menu \
		$MV.complex_menu ]

set cmd [ list raiseMenu $MV.mode_menu $MV.modes $MV $mode_lower_subs ]

button $MV.modes -text "Modes..." -padx 1 -pady 1 \
	-font $main_font \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color \
	-command $cmd

restrict_bindings $MV.modes ""

#
# Data Reduction Function Menu
#

# eem puts -nonewline "."
# eem flush stdout

set reduce_menu_list \
	[ list \
		[ list "Avg"	radiobutton		avg		reduction_function \
			"Set Data Reduction Function to Average Value." ] \
		[ list "Min"	radiobutton		min		reduction_function \
			"Set Data Reduction Function to Minimum Value." ] \
		[ list "Max"	radiobutton		max		reduction_function \
			"Set Data Reduction Function to Maximum Value." ] \
		[ list "Sum"	radiobutton		sum		reduction_function \
			"Set Data Reduction Function to Sum of Values." ] \
		[ list "Midpoint"	radiobutton	midpoint	reduction_function \
			"Set Data Reduction Function to Midpoint Value." ] \
		[ list "Count"	radiobutton		count	reduction_function \
		"Set Data Reduction Function to Number of Non-Zero Values." ] \
		[ list "Std"	radiobutton		std		reduction_function \
	"Set Data Reduction Function to Standard Deviation of Values." ] \
		[ list "Var"	radiobutton		var		reduction_function \
			"Set Data Reduction Function to Variance of Values." ] \
	]

set reduction_function [ get_reduction ]

makeMenu $MV.reduce_menu reduceHandle $reduce_menu_list lower {} none

set reduce_lower_subs \
	[ list $MV.file_menu $MV.mode_menu $MV.options_menu \
		$MV.complex_menu ]

set cmd [ list raiseMenu $MV.reduce_menu $MV.reduce $MV \
	$reduce_lower_subs ]

button $MV.reduce -text "Reduction..." -padx 1 -pady 1 \
	-font $main_font \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color \
	-command $cmd

restrict_bindings $MV.reduce ""

#
# Options Menu
#

# eem puts -nonewline "."
# eem flush stdout

set options_menu_list \
	[ list \
		[ list "Auto Color Map"			statebutton	auto_color	OFF \
			state_auto_color \
		"Toggle Auto Color Map Ranging (By Data Value Min & Max)." ] \
		[ list "Greyscale"				statebutton	greyscale	OFF \
			state_greyscale \
			"Toggle Greyscale Versus Color View Drawing." ] \
		[ list "Filter Values"			statebutton	\
			filter_values	OFF \
			state_filter_values \
			"Toggle Filtering of Values Outside of Color Range." ] \
		[ list "View Absolute Values"	statebutton	absvals		OFF \
			view_abs_vals \
			"Toggle Viewing of Absolute Values of Matrix Values." ] \
		[ list "View Pattern Only"		statebutton	pattern		OFF \
			pattern_only \
		"Toggle Viewing of Matrix Pattern (Ignore Element Values)." ] \
		[ list "Transpose Matrix"		statebutton	transpose	OFF \
			transpose_matrix \
			"Toggle Transposing of Matrix." ] \
		[ list "Separate View Window"	statebutton	sep_view	OFF \
			use_sep_view \
			"Toggle Use of Separate Window for Viewing Matrix." ] \
		[ list "Show Min & Max"			statebutton	min_max		OFF \
			state_min_max \
			"Toggle Display of Data Value Min & Max." ] \
		[ list "Show Cumulative"		statebutton	cumulative	OFF \
			state_cumulative \
			"Toggle Display of Cumulative Min & Max Stats." ] \
		[ list "Reset Cumulative"		command		reset_cumulative \
			"Reset Cumulative Min & Max Stats." ] \
	]

makeMenu $MV.options_menu optionHandle $options_menu_list lower {} none

set options_lower_subs \
	[ list $MV.file_menu $MV.mode_menu $MV.reduce_menu \
		$MV.complex_menu ]

set cmd [ list raiseMenu $MV.options_menu $MV.options $MV \
	$options_lower_subs ]

button $MV.options -text "Options..." -padx 1 -pady 1 \
	-font $main_font \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color \
	-command $cmd

restrict_bindings $MV.options ""

#
# Complex View Types
#

# eem puts -nonewline "."
# eem flush stdout

set complex_menu_list \
	[ list \
		[ list "Re"		radiobutton		re		complex_view \
			"Set Complex View Type to Real Component Value." ] \
		[ list "Im"		radiobutton		im		complex_view \
			"Set Complex View Type to Imaginary Component Value." ] \
		[ list "Magn"	radiobutton		magn	complex_view \
			"Set Complex View Type to Magnitude of Complex Vector." ] \
		[ list "Phase"	radiobutton		phase	complex_view \
			"Set Complex View Type to Phase of Complex Vector." ] \
	]

set complex_view [ get_complex_view ]

set is_complex OFF

makeMenu $MV.complex_menu complexHandle $complex_menu_list lower {} none

set complex_lower_subs \
	[ list $MV.file_menu $MV.mode_menu $MV.reduce_menu \
		$MV.options_menu ]

set cmd [ list raiseMenu $MV.complex_menu $MV.complex $MV \
	$complex_lower_subs ]

button $MV.complex -text "Complex Views..." -padx 1 -pady 1 \
	-font $main_font \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color \
	-command $cmd

restrict_bindings $MV.complex ""

#
# Reduction Status Indicator
#

# eem puts -nonewline "."
# eem flush stdout

label $MV.rstat_label -text {Reduction:  none} \
	-foreground $fg_color -background $bg_color \
	-font $main_font -relief flat -anchor nw

restrict_bindings $MV.rstat_label "Label"

set rstat_width 125

#
# Matrix Grid Interface
#

# eem puts -nonewline "."
# eem flush stdout

checkbutton $MV.grid -text {Grid} \
	-padx 1 -pady 1 \
	-command gridHandle \
	-onvalue ON -offvalue OFF -variable grid_status \
	-bd $frame_border -relief raised \
	-selectcolor $selector_color -foreground $fg_color \
	-activeforeground $active_fg_color

restrict_bindings $MV.grid "Checkbutton"

label $MV.grid_row -text {Rows:} \
	-foreground $fg_color -background $bg_color \
	-font $main_font -relief flat -anchor nw

restrict_bindings $MV.grid_row "Label"

entry $MV.grid_row_entry -font $main_font \
	-bd $frame_border -relief sunken \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.grid_row_entry "Entry"

entry_setup $MV.grid_row_entry "gridHandle"

$MV.grid_row_entry insert 0 100

label $MV.grid_roff -text {Off:} \
	-foreground $fg_color -background $bg_color \
	-font $main_font -relief flat -anchor nw

restrict_bindings $MV.grid_roff "Label"

entry $MV.grid_roff_entry -font $main_font \
	-bd $frame_border -relief sunken \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.grid_roff_entry "Entry"

entry_setup $MV.grid_roff_entry "gridHandle"

$MV.grid_roff_entry insert 0 0

label $MV.grid_col -text {Cols:} \
	-foreground $fg_color -background $bg_color \
	-font $main_font -relief flat -anchor nw

restrict_bindings $MV.grid_col "Label"

entry $MV.grid_col_entry -font $main_font \
	-bd $frame_border -relief sunken \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.grid_col_entry "Entry"

entry_setup $MV.grid_col_entry "gridHandle"

$MV.grid_col_entry insert 0 100

label $MV.grid_coff -text {Off:} \
	-foreground $fg_color -background $bg_color \
	-font $main_font -relief flat -anchor nw

restrict_bindings $MV.grid_coff "Label"

entry $MV.grid_coff_entry -font $main_font \
	-bd $frame_border -relief sunken \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.grid_coff_entry "Entry"

entry_setup $MV.grid_coff_entry "gridHandle"

$MV.grid_coff_entry insert 0 0

set grid_status OFF

#
# Default Zero Color Setting
#

# eem puts -nonewline "."
# eem flush stdout

checkbutton $MV.zero_color -text {Zero Color:} \
	-padx 1 -pady 1 \
	-command handleZeroColor \
	-onvalue ON -offvalue OFF -variable fixed_zero_color \
	-bd $frame_border -relief flat \
	-selectcolor $selector_color -foreground $fg_color \
	-activeforeground $active_fg_color

restrict_bindings $MV.zero_color "Checkbutton"

entry $MV.zero_entry -font $main_font \
	-bd $frame_border -relief sunken \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.zero_entry "Entry"

entry_setup $MV.zero_entry "handleZeroColor"

$MV.zero_entry insert 0 white

set fixed_zero_color OFF

#
# Matrix Selection
#

# eem puts -nonewline "."
# eem flush stdout

button $MV.mat -text "Matrix File:" \
	-padx 1 -pady 1 -font $main_font \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color \
	-command "matHandle browse"

restrict_bindings $MV.mat ""

entry $MV.mat_entry -font $main_font \
	-bd $frame_border -relief sunken \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.mat_entry "Entry"

entry_setup $MV.mat_entry "matHandle set"

$MV.mat_entry insert 0 [ get_mat_name ]

#
# Plane Selection
#

# eem puts -nonewline "."
# eem flush stdout

label $MV.plane_label -text {Plane Selection:} \
	-foreground $fg_color -background $bg_color \
	-font $main_font -relief flat -anchor nw

restrict_bindings $MV.plane_label "Label"

radiobutton $MV.xy_plane -text { X-Y Plane } \
	-padx 1 -pady 1 \
	-command "handlePlane" \
	-bd $frame_border -relief sunken -value "XY" \
	-variable plane_choice \
	-selectcolor $selector_color \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color

restrict_bindings $MV.xy_plane ""

radiobutton $MV.yz_plane -text { Y-Z Plane } \
	-padx 1 -pady 1 \
	-command "handlePlane" \
	-bd $frame_border -relief sunken -value "YZ" \
	-variable plane_choice \
	-selectcolor $selector_color \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color

restrict_bindings $MV.yz_plane ""

radiobutton $MV.xz_plane -text { X-Z Plane } \
	-padx 1 -pady 1 \
	-command "handlePlane" \
	-bd $frame_border -relief sunken -value "XZ" \
	-variable plane_choice \
	-selectcolor $selector_color \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color

restrict_bindings $MV.xz_plane ""

set plane_choice "XY"

#
# Vis Region Bounds & Cell
#

# eem puts -nonewline "."
# eem flush stdout

# "X" Axis

label $MV.xrange -text "X Axis Range:" \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.xrange "Label"

label $MV.xrange_lb -text "(" \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.xrange_lb "Label"

label $MV.xrange_cm -text "," \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.xrange_cm "Label"

label $MV.xrange_rb -text ")" \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.xrange_rb "Label"

entry $MV.xmin_entry -font $main_font \
	-bd $frame_border -relief sunken \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.xmin_entry "Entry"

entry_setup $MV.xmin_entry \
	"handleVisRegion $MV.xmin_entry $MV.xmax_entry"

$MV.xmin_entry insert 0 0

entry $MV.xmax_entry -font $main_font \
	-bd $frame_border -relief sunken \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.xmax_entry "Entry"

entry_setup $MV.xmax_entry \
	"handleVisRegion $MV.xmax_entry $MV.xcell_entry"

$MV.xmax_entry insert 0 0

label $MV.xcell -text "Cell Size:" \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.xcell "Label"

entry $MV.xcell_entry -font $main_font \
	-bd $frame_border -relief sunken \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.xcell_entry "Entry"

entry_setup $MV.xcell_entry \
	"handleVisRegion $MV.xcell_entry $MV.ymin_entry"

$MV.xcell_entry insert 0 1

# "Y" Axis

label $MV.yrange -text "Y Axis Range:" \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.yrange "Label"

label $MV.yrange_lb -text "(" \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.yrange_lb "Label"

label $MV.yrange_cm -text "," \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.yrange_cm "Label"

label $MV.yrange_rb -text ")" \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.yrange_rb "Label"

entry $MV.ymin_entry -font $main_font \
	-bd $frame_border -relief sunken \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.ymin_entry "Entry"

entry_setup $MV.ymin_entry \
	"handleVisRegion $MV.ymin_entry $MV.ymax_entry"

$MV.ymin_entry insert 0 0

entry $MV.ymax_entry -font $main_font \
	-bd $frame_border -relief sunken \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.ymax_entry "Entry"

entry_setup $MV.ymax_entry \
	"handleVisRegion $MV.ymax_entry $MV.ycell_entry"

$MV.ymax_entry insert 0 0

label $MV.ycell -text "Cell Size:" \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.ycell "Label"

entry $MV.ycell_entry -font $main_font \
	-bd $frame_border -relief sunken \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.ycell_entry "Entry"

entry_setup $MV.ycell_entry \
	"handleVisRegion $MV.ycell_entry $MV.zval_entry"

$MV.ycell_entry insert 0 1

# "Z" Axis

label $MV.zval -text "Z Axis Value:" \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.zval "Label"

entry $MV.zval_entry -font $main_font \
	-bd $frame_border -relief sunken \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.zval_entry "Entry"

entry_setup $MV.zval_entry \
	"handleVisRegion $MV.zval_entry $MV.xmin_entry"

$MV.zval_entry insert 0 0

#
# Min & Max Stats
#

# eem puts -nonewline "."
# eem flush stdout

label $MV.minval -text "Min:" \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.minval "Label"

label $MV.maxval -text "Max:" \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.maxval "Label"

#
# Navigation Canvas
#

# eem puts -nonewline "."
# eem flush stdout

canvas $MV_N -bd $frame_border -relief sunken -confine 1 \
	-background $bg_color

restrict_bindings $MV_N "Canvas"

set nav_border [ expr 2 * $frame_border ]

set navigate_size [ expr 32 + ( 2 * $nav_border ) ]

set ckdir_compass [ file exists $matview_directory/compass.xbm ]

if { $ckdir_compass == 1 } \
{
	load_bitmap_file $matview_directory/compass.xbm

	$MV_N create bitmap $nav_border $nav_border \
		-anchor nw -foreground purple \
		-bitmap "@$matview_directory/compass.xbm"
} \
\
else \
{
	set ckdir_compass [ file exists compass.xbm ]

	if { $ckdir_compass == 1 } \
	{
		load_bitmap_file compass.xbm

		$MV_N create bitmap $nav_border $nav_border \
			-anchor nw -foreground purple \
			-bitmap "@compass.xbm"
	} \
\
	else \
	{
		puts "Warning: Matview 'compass.xbm' File Not Found!"
		puts "MATVIEW_ROOT not set correctly?"
	}
}

#
# Color Map Canvas
#

# eem puts -nonewline "."
# eem flush stdout

canvas $MV_M -bd $frame_border -relief sunken -confine 1 \
	-background $bg_color

restrict_bindings $MV_M "Canvas"

entry $MV.color_min -font $main_font \
	-bd $frame_border -relief sunken \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.color_min "Entry"

entry_setup $MV.color_min "handleColor set"

$MV.color_min insert 0 0

entry $MV.color_max -font $main_font \
	-bd $frame_border -relief sunken \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.color_max "Entry"

entry_setup $MV.color_max "handleColor set"

$MV.color_max insert 0 0

scrollbar $MV.color_scroll -orient horiz \
	-bd $frame_border -relief sunken \
	-command "handleColor scroll"

restrict_bindings $MV.color_scroll "Scrollbar"

set color_scroll_value 0.495

#
# Color Function Graph Canvas
#

# eem puts -nonewline "."
# eem flush stdout

canvas $MV_G -bd $frame_border -relief sunken -confine 1 \
	-background $bg_color

restrict_bindings $MV_G "Canvas"

#
# Main Canvas
#

# eem puts -nonewline "."
# eem flush stdout

canvas $MV_C -bd $frame_border -relief sunken -confine 1 \
	-background "#ffffff"

restrict_bindings $MV_C "Canvas"

image create photo $MV_I

set image_height [ winfo screenheight $MV ]
set image_width [ winfo screenwidth $MV ]

# eem puts -nonewline "."
# eem flush stdout

set rect_id [ $MV_C create rectangle \
	0 0 $image_width $image_height -fill white ]

# eem puts -nonewline "."
# eem flush stdout

set image_id [ $MV_C create image 0 0 -image $MV_I -anchor nw ]

# eem puts -nonewline "."
# eem flush stdout

#
# Thumb Canvas
#

# eem puts -nonewline "."
# eem flush stdout

canvas $MV_TC -bd $frame_border -relief sunken -confine 1 \
	-background "#ffffff"

restrict_bindings $MV_TC "Canvas"

image create photo $MV_TI

set thumb_height 64
set thumb_width 64

set thumb_rect_id [ $MV_TC create rectangle \
	0 0 $thumb_width $thumb_height -fill white ]

set thumb_image_id [ $MV_TC create image 0 0 -image $MV_TI -anchor nw ]

#
# View Region Size Display
#

# eem puts -nonewline "."
# eem flush stdout

label $MV.region_size -text {Current Region Size:} \
	-foreground $fg_color -background $bg_color

restrict_bindings $MV.region_size "Label"

#
# Progress Bar
#

# eem puts -nonewline "."
# eem flush stdout

canvas $MV_P -bd $frame_border -relief sunken -confine 1 \
	-background $bg_color

restrict_bindings $MV_P "Canvas"

set progress_height [ expr $border_space + $frame_border ]

set progress_width 150

set progress_bar \
	[ $MV_P create rectangle 0 0 0 $progress_height -fill purple ]

#
# Initialize Matrix Resize Lock */
#

set resize_mat_lock 0

#
# Layout Main Panel
#

# eem puts -nonewline "."
# eem flush stdout

layout_main_panel

#
# Pass In Default Color Map Settings
#

handleColor init $color_scroll_value

#
# Create Main Query Pop-Up
#

# eem puts -nonewline "."
# eem flush stdout

frame $MV.query -bd 4 -relief raised

restrict_bindings $MV.query "Frame"

label $MV.query.label -text "XXX" \
	-foreground $fg_color

restrict_bindings $MV.query.label "Label"

place $MV.query.label -relx 0.50 -y [ expr $border_space / 2 ] -anchor n

update

place forget $MV.query

#
# Separate Toplevel Viewing Window
#

# eem puts -nonewline "."
# eem flush stdout

set VF .view_field

set VF_C $VF.canvas

set VF_P $VF.progress

set view_height 600
set view_width 600

toplevel $VF

wm withdraw $VF

wm geometry $VF "${view_width}x${view_height}+200+150"

wm minsize $VF 25 25

$VF configure -background $bg_color

canvas $VF_C -bd $frame_border -relief sunken -confine 1 \
	-background $bg_color

restrict_bindings $VF_C "Canvas"

set vf_height [ winfo screenheight $VF ]
set vf_width [ winfo screenwidth $VF ]

# eem puts -nonewline "."
# eem flush stdout

set vf_rect_id [ $VF_C create rectangle \
	0 0 $vf_width $vf_height -fill white ]

# eem puts -nonewline "."
# eem flush stdout

set vf_image_id [ $VF_C create image 0 0 -image $MV_I -anchor nw ]

# eem puts -nonewline "."
# eem flush stdout

#
# Separate View Region Size Display
#

label $VF.region_size -text {Current Region Size:} \
	-foreground $fg_color -background $bg_color

restrict_bindings $VF.region_size "Label"

#
# Separate Progress Bar
#

canvas $VF_P -bd $frame_border -relief sunken -confine 1 \
	-background $bg_color

restrict_bindings $VF_P "Canvas"

set sep_progress_bar \
	[ $VF_P create rectangle 0 0 0 $progress_height -fill purple ]

layout_vf_panel

#
# Create Separate Viewing Window Query Pop-Up
#

# eem puts -nonewline "."
# eem flush stdout

frame $VF.query -bd 4 -relief raised

restrict_bindings $VF.query "Frame"

label $VF.query.label -text "XXX" \
	-foreground $fg_color

restrict_bindings $VF.query.label "Label"

place $VF.query.label -relx 0.50 -y [ expr $border_space / 2 ] -anchor n

update

place forget $VF.query

#
# Pop-Up Print Dialog
#

# eem puts -nonewline "."
# eem flush stdout

set PD .print_image

set PD_P $PD.progress

set pd_height 400
set pd_width 400

toplevel $PD

wm withdraw $PD

wm geometry $PD "${pd_width}x${pd_height}+250+325"

wm minsize $PD 100 100

$PD configure -background $bg_color

label $PD.dest_label -text "Print To:" \
	-foreground $fg_color

restrict_bindings $PD.dest_label "Label"

radiobutton $PD.dest_printer -text { Printer } \
	-padx 1 -pady 1 \
	-command "printHandle dest" \
	-bd $frame_border -relief sunken -value "printer" \
	-variable print_dest \
	-selectcolor $selector_color \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color

restrict_bindings $PD.dest_printer ""

radiobutton $PD.dest_file -text { File } \
	-padx 1 -pady 1 \
	-command "printHandle dest" \
	-bd $frame_border -relief sunken -value "file" \
	-variable print_dest \
	-selectcolor $selector_color \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color

restrict_bindings $PD.dest_file ""

set print_dest "printer"

label $PD.pcmd_label -text "Printer Command:" \
	-foreground $fg_color

restrict_bindings $PD.pcmd_label "Label"

entry $PD.pcmd_entry -font $main_font \
	-bd $frame_border -relief sunken \
	-foreground $fg_color -background $bg_color

restrict_bindings $PD.pcmd_entry "Entry"

entry_setup $PD.pcmd_entry "printHandle pcmd"

$PD.pcmd_entry insert 0 "lpr -Plw"

label $PD.file_label -text "File Name:" \
	-foreground $fg_color

restrict_bindings $PD.file_label "Label"

entry $PD.file_entry -font $main_font \
	-bd $frame_border -relief sunken \
	-foreground $fg_color -background $bg_color

restrict_bindings $PD.file_entry "Entry"

entry_setup $PD.file_entry "printHandle filename"

$PD.file_entry insert 0 "image.ps"

button $PD.file_browse -text "Browse..." -padx 1 -pady 1 \
	-font $main_font \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color \
	-command "printHandle file_browse"

restrict_bindings $PD.file_browse ""

label $PD.format_label -text "File Format:" \
	-foreground $fg_color

restrict_bindings $PD.format_label "Label"

radiobutton $PD.format_ps -text { PostScript } \
	-padx 1 -pady 1 \
	-command "printHandle format" \
	-bd $frame_border -relief sunken -value "ps" \
	-variable print_format \
	-selectcolor $selector_color \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color

restrict_bindings $PD.format_ps ""

radiobutton $PD.format_ppm -text { PPM } \
	-padx 1 -pady 1 \
	-command "printHandle format" \
	-bd $frame_border -relief sunken -value "ppm" \
	-variable print_format \
	-selectcolor $selector_color \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color

restrict_bindings $PD.format_ppm ""

set print_format "ps"

label $PD.orient_label -text "Orientation:" \
	-foreground $fg_color

restrict_bindings $PD.orient_label "Label"

radiobutton $PD.orient_portrait -text { Portrait } \
	-padx 1 -pady 1 \
	-command "printHandle orient" \
	-bd $frame_border -relief sunken -value "portrait" \
	-variable print_orient \
	-selectcolor $selector_color \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color

restrict_bindings $PD.orient_portrait ""

radiobutton $PD.orient_landscape -text { Landscape } \
	-padx 1 -pady 1 \
	-command "printHandle orient" \
	-bd $frame_border -relief sunken -value "landscape" \
	-variable print_orient \
	-selectcolor $selector_color \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color

restrict_bindings $PD.orient_landscape ""

set print_orient "portrait"

label $PD.media_label -text "Media:" \
	-foreground $fg_color

restrict_bindings $PD.media_label "Label"

radiobutton $PD.media_color -text { Color } \
	-padx 1 -pady 1 \
	-command "printHandle media" \
	-bd $frame_border -relief sunken -value "color" \
	-variable print_media \
	-selectcolor $selector_color \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color

restrict_bindings $PD.media_color ""

radiobutton $PD.media_greyscale -text { Greyscale } \
	-padx 1 -pady 1 \
	-command "printHandle media" \
	-bd $frame_border -relief sunken -value "greyscale" \
	-variable print_media \
	-selectcolor $selector_color \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color

restrict_bindings $PD.media_greyscale ""

set print_media "color"

button $PD.preview -text "Preview" -padx 1 -pady 1 \
	-font $main_font \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color \
	-command "printHandle preview"

restrict_bindings $PD.preview ""

button $PD.print -text "Print" -padx 1 -pady 1 \
	-font $main_font \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color \
	-command "printHandle print"

restrict_bindings $PD.print ""

button $PD.cancel -text "Cancel" -padx 1 -pady 1 \
	-font $main_font \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color \
	-command "printHandle cancel"

restrict_bindings $PD.cancel ""

canvas $PD_P -bd $frame_border -relief sunken -confine 1 \
	-background $bg_color

restrict_bindings $PD_P "Canvas"

set pd_progress_bar \
	[ $PD_P create rectangle 0 0 0 $progress_height -fill purple ]

layout_print_dialog

#
# Printing View Toplevel
#

# eem puts -nonewline "."
# eem flush stdout

set PV .print_view

set PV_C $PV.canvas

set PV_I $PV.image

set pview_field_height 600
set pview_field_width 600

toplevel $PV

wm withdraw $PV

wm geometry $PV "${pview_field_width}x${pview_field_height}+200+250"

wm resizable $PV 0 0

$PV configure -background $bg_color

canvas $PV_C -bd $frame_border -relief sunken -confine 1 \
	-background $bg_color

restrict_bindings $PV_C "Canvas"

image create photo $PV_I

set pv_height [ winfo screenheight $PV ]
set pv_width [ winfo screenwidth $PV ]

# eem puts -nonewline "."
# eem flush stdout

set pv_image_id -1

button $PV.print -text "Print" -padx 1 -pady 1 \
	-font $main_font \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color \
	-command "printHandle preview_print"

restrict_bindings $PV.print ""

button $PV.refresh -text "Refresh" -padx 1 -pady 1 \
	-font $main_font \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color \
	-command "printHandle preview"

restrict_bindings $PV.refresh ""

button $PV.props -text "Properties" -padx 1 -pady 1 \
	-font $main_font \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color \
	-command "print"

restrict_bindings $PV.props ""

button $PV.cancel -text "Cancel" -padx 1 -pady 1 \
	-font $main_font \
	-foreground $fg_color -background $bg_color \
	-activeforeground $active_fg_color \
	-command "printHandle cancel"

restrict_bindings $PV.cancel ""

layout_pv_panel

#
# Main Event Bindings
#

bind $MV.mat_entry		<Enter>		"focus $MV.mat_entry"
bind $MV.mat_entry		<Leave>		"matHandle set"

bind $MV.xmin_entry		<Enter>		"focus $MV.xmin_entry"
bind $MV.xmin_entry		<Leave> \
	"handleVisRegion $MV.xmin_entry $MV.xmax_entry"
bind $MV.xmin_entry		<Tab> \
	"handleVisRegion $MV.xmin_entry $MV.xmax_entry"

bind $MV.xmax_entry		<Enter>		"focus $MV.xmax_entry"
bind $MV.xmax_entry		<Leave> \
	"handleVisRegion $MV.xmax_entry $MV.xcell_entry"
bind $MV.xmax_entry		<Tab> \
	"handleVisRegion $MV.xmax_entry $MV.xcell_entry"

bind $MV.xcell_entry	<Enter>		"focus $MV.xcell_entry"
bind $MV.xcell_entry	<Leave> \
	"handleVisRegion $MV.xcell_entry $MV.ymin_entry"
bind $MV.xcell_entry	<Tab> \
	"handleVisRegion $MV.xcell_entry $MV.ymin_entry"

bind $MV.ymin_entry		<Enter>		"focus $MV.ymin_entry"
bind $MV.ymin_entry		<Leave> \
	"handleVisRegion $MV.ymin_entry $MV.ymax_entry"
bind $MV.ymin_entry		<Tab> \
	"handleVisRegion $MV.ymin_entry $MV.ymax_entry"

bind $MV.ymax_entry		<Enter>		"focus $MV.ymax_entry"
bind $MV.ymax_entry		<Leave> \
	"handleVisRegion $MV.ymax_entry $MV.ycell_entry"
bind $MV.ymax_entry		<Tab> \
	"handleVisRegion $MV.ymax_entry $MV.ycell_entry"

bind $MV.ycell_entry	<Enter>		"focus $MV.ycell_entry"
bind $MV.ycell_entry	<Leave> \
	"handleVisRegion $MV.ycell_entry $MV.xmin_entry"
bind $MV.ycell_entry	<Tab> \
	"handleVisRegion $MV.ycell_entry $MV.xmin_entry"

bind $MV.zval_entry		<Enter>		"focus $MV.zval_entry"
bind $MV.zval_entry		<Leave> \
	"handleVisRegion $MV.zval_entry $MV.xmin_entry"
bind $MV.zval_entry		<Tab> \
	"handleVisRegion $MV.zval_entry $MV.xmin_entry"

bind $MV.color_min		<Enter>		"focus $MV.color_min"
bind $MV.color_min		<Leave>		"handleColor set"

bind $MV.color_max		<Enter>		"focus $MV.color_max"
bind $MV.color_max		<Leave>		"handleColor set"

bind $MV_N <ButtonPress>			"canvas_handle nav_press %x %y"

set nav_help \
	{Click Mouse Button to Navigate Laterally Around Matrix.}

set nav_help_cmd [ list setTmpMsg Help $nav_help ]

bind $MV_N <Enter>					"focus $MV_N; $nav_help_cmd"

bind $MV_N <Leave>					"popMsgs"

bind $MV_C <ButtonPress-1>			"canvas_handle query_press %x %y"
bind $MV_C <Button1-Motion>			"canvas_handle query_slide %x %y"
bind $MV_C <ButtonRelease-1>		"canvas_handle query_release %x %y"

bind $MV_C <ButtonPress-2>			"canvas_handle zoom_press %x %y"
bind $MV_C <Button2-Motion>			"canvas_handle zoom_slide %x %y"
bind $MV_C <ButtonRelease-2>		"canvas_handle zoom_release %x %y"

bind $MV_C <3>						"canvas_handle pop %x %y"

set canvas_help \
	{Mouse Buttons: Left=Query, Middle_Sweep=Zoom_Area, Right=Unzoom.}

set canvas_help_cmd [ list setTmpMsg Help $canvas_help ]

bind $MV_C <Enter>					"focus $MV_C; \
										canvas_handle enter %x %y; \
										$canvas_help_cmd"

bind $MV_C <Motion>					"focus $MV_C; \
										canvas_handle move %x %y"

bind $MV_C <Leave>					"canvas_handle leave %x %y; popMsgs"

bind $MV_C <KeyPress>				"canvas_handle %K %x %y"

bind $MV_TC <ButtonPress-1>			"thumb_nav press %x %y"
bind $MV_TC <Button1-Motion>		"thumb_nav slide %x %y"
bind $MV_TC <ButtonRelease-1>		"thumb_nav release %x %y"

bind $MV_TC <ButtonPress-2>			"thumb_nav region_press %x %y"
bind $MV_TC <Button2-Motion>		"thumb_nav region_slide %x %y"
bind $MV_TC <ButtonRelease-2>		"thumb_nav region_release %x %y"

bind $MV_TC <3>						"canvas_handle pop %x %y"

set thumb_canvas_help \
	{Click Mouse Button to Move Current View Window.}

set thumb_canvas_help_cmd [ list setTmpMsg Help $thumb_canvas_help ]

bind $MV_TC <Enter>					"focus $MV_TC; $canvas_help_cmd"

bind $MV_TC <Motion>				"focus $MV_TC"

bind $MV_TC <Leave>					"popMsgs"

bind $MV <Configure>				"resize_main_panel"

bind $MV <Destroy>					"quit"

#
# View Field Event Bindings
#

bind $VF_C <ButtonPress-1>			"canvas_handle query_press %x %y"
bind $VF_C <Button1-Motion>			"canvas_handle query_slide %x %y"
bind $VF_C <ButtonRelease-1>		"canvas_handle query_release %x %y"

bind $VF_C <ButtonPress-2>			"canvas_handle zoom_press %x %y"
bind $VF_C <Button2-Motion>			"canvas_handle zoom_slide %x %y"
bind $VF_C <ButtonRelease-2>		"canvas_handle zoom_release %x %y"

bind $VF_C <3>						"canvas_handle pop %x %y"

bind $VF_C <Enter>					"focus $VF_C; \
										canvas_handle enter %x %y; \
										$canvas_help_cmd"

bind $VF_C <Motion>					"focus $VF_C; \
										canvas_handle move %x %y"

bind $VF_C <Leave>					"canvas_handle leave %x %y; popMsgs"

bind $VF_C <KeyPress>				"canvas_handle %K %x %y"

bind $VF <Configure>				"resize_vf_panel"

#
# Button Help Bindings
#

# eem puts -nonewline "."
# eem flush stdout

# main panel buttons

butt_help $MV.message label "MatView Status Messages..."

butt_help $MV.help_msg label \
	"This is where Help Messages are Displayed..."

butt_help $MV.file button \
"Raise / Lower Main File Menu - to Exit MatView or Set Saved Fields."

butt_help $MV.modes button \
	"Raise / Lower Modes Menu - to Set Normal or Differential Viewing."

butt_help $MV.options button \
	"Raise / Lower Options Menu - to Configure Matview Options."

butt_help $MV.reduce button \
"Raise / Lower Reduction Function Menu - to Choose How Data is Reduced."

butt_help $MV.complex button \
"Raise / Lower Complex Views Menu - to Control Views of Complex Values."

# grid stuff

butt_help $MV.grid label \
	"Toggle Use of Matrix Grid."

butt_help $MV.grid_row label \
	"Number of Matrix Rows Per Grid Block."

butt_help $MV.grid_row_entry entry \
	[ list "gridHandle" \
		"Enter the Number of Matrix Rows Per Grid Block." ]

butt_help $MV.grid_roff label \
	"Offset of Matrix Rows for First Grid Block."

butt_help $MV.grid_roff_entry entry \
	[ list "gridHandle" \
		"Enter the Offset of Matrix Rows for the First Grid Block." ]

butt_help $MV.grid_col label \
	"Number of Matrix Columns Per Grid Block."

butt_help $MV.grid_col_entry entry \
	[ list "gridHandle" \
		"Enter the Number of Matrix Columns Per Grid Block." ]

butt_help $MV.grid_coff label \
	"Offset of Matrix Columns for First Grid Block."

butt_help $MV.grid_coff_entry entry \
	[ list "gridHandle" \
		"Enter the Offset of Matrix Columns for the First Grid Block." ]

# zero color stuff

butt_help $MV.zero_color label \
	"Toggle Setting of Fixed Zero Color."

butt_help $MV.zero_entry entry \
	[ list "handleZeroColor" \
		"Enter the Fixed Zero Value Color (Text Name or #rrggbb)." ]

# region size display

butt_help $MV.region_size label \
	"Shows the Size of the Currently Displayed Portion of the Matrix."

# progress bar & reduction status

butt_help $MV_P label \
"Progress Indicator - Shows Completed/Remaining Internal Processing."

butt_help $MV.rstat_label label \
	"Reduction Status - Displays Current Data Reduction Function."

# matrix file

butt_help $MV.mat button \
	"Browse to Set the Matrix File to be Displayed by MatView."

butt_help $MV.mat_entry entry \
	[ list "matHandle set" \
	"Enter the Name of the Matrix File to be Displayed by MatView." ]

# slice plane selection

butt_help $MV.plane_label label \
	"Selection of the 2-D Plane to be Sliced Along for Viewing."

butt_help $MV.xy_plane button \
	"Select the X-Y Plane as the Slice Plane for Viewing."

butt_help $MV.yz_plane button \
	"Select the Y-Z Plane as the Slice Plane for Viewing."

butt_help $MV.xz_plane button \
	"Select the X-Z Plane as the Slice Plane for Viewing."

# primary axis range

butt_help $MV.xrange label \
	"Range of Element Indices for the Primary Axis."

butt_help $MV.xmin_entry entry \
	[ list "handleVisRegion $MV.xmin_entry $MV.xmax_entry" \
	"Enter the Minimum Element Index for the Primary Axis." ]

butt_help $MV.xmax_entry entry \
	[ list "handleVisRegion $MV.xmax_entry $MV.xcell_entry" \
	"Enter the Maximum Element Index for the Primary Axis." ]

butt_help $MV.xcell label \
	"Cell Size (Resolution) of Data Samples for Primary Axis."

butt_help $MV.xcell_entry entry \
	[ list "handleVisRegion $MV.xcell_entry $MV.ymin_entry" \
	"Enter the Data Sample Cell Size for the Primary Axis." ]

# secondary axis range

butt_help $MV.yrange label \
	"Range of Element Indices for the Secondary Axis."

butt_help $MV.ymin_entry entry \
	[ list "handleVisRegion $MV.ymin_entry $MV.ymax_entry" \
	"Enter the Minimum Element Index for the Secondary Axis." ]

butt_help $MV.ymax_entry entry \
	[ list "handleVisRegion $MV.ymax_entry $MV.ycell_entry" \
	"Enter the Maximum Element Index for the Secondary Axis." ]

butt_help $MV.ycell label \
	"Cell Size (Resolution) of Data Samples for Secondary Axis."

butt_help $MV.ycell_entry entry \
	[ list "handleVisRegion $MV.ycell_entry $MV.zval_entry" \
	"Enter the Data Sample Cell Size for the Secondary Axis." ]

# slice plane index

butt_help $MV.zval label \
	"Element Index for the Location of the Slice Plane."

butt_help $MV.zval_entry entry \
	[ list "handleVisRegion $MV.zval_entry $MV.xmin_entry" \
	"Enter the Element Index for the Location of the Slice Plane." ]

# color map

butt_help $MV.color_min entry \
	[ list "handleColor set" \
	"Enter the Minimum Data Value for the Color Mapping Range." ]

butt_help $MV.color_max entry \
	[ list "handleColor set" \
	"Enter the Maximum Data Value for the Color Mapping Range." ]

butt_help $MV_M label \
	"Color Map for Data Value Range."

butt_help $MV.color_scroll label \
	"Slide Scrollbar to Adjust Color Map Progression Function."

butt_help $MV_G label \
	"Color Map Distribution Function for Mapping Data Values to Colors."

# separate region size display

butt_help $VF.region_size label \
	"Shows the Size of the Currently Displayed Portion of the Matrix."

# separate progress bar

butt_help $VF_P label \
"Progress Indicator - Shows Completed/Remaining Internal Processing."

# print dialog

butt_help $PD.dest_label label \
	"Select Destination for Image - Print or Save to File."

butt_help $PD.dest_printer button \
	"Send Image to Printer."

butt_help $PD.dest_file button \
	"Save Image to File."

butt_help $PD.pcmd_label label \
	"Command for Sending Image to Printer."

butt_help $PD.pcmd_entry entry \
	[ list "printHandle pcmd" \
		"Enter Command for Sending Image to Printer." ]

butt_help $PD.file_label label \
	"Name of File for Saving Image."

butt_help $PD.file_entry entry \
	[ list "printHandle filename" \
		"Enter Name of File for Saving Image." ]

butt_help $PD.file_browse button \
	"Browse to Set File for Saving Image."

butt_help $PD.format_label label \
	"Select File Format for Saving Image."

butt_help $PD.format_ps button \
	"Save Image to PostScript File."

butt_help $PD.format_ppm button \
	"Save Image to PPM File."

butt_help $PD.orient_label label \
	"Select Orientation for Printing / Saving Image."

butt_help $PD.orient_portrait button \
	"Print / Save Image in Portrait Orientation."

butt_help $PD.orient_landscape button \
	"Print / Save Image in Landscape Orientation."

butt_help $PD.media_label label \
	"Select Media for Printing / Saving Image (Color vs. Greyscale)."

butt_help $PD.media_color button \
	"Print / Save Image for Color Media."

butt_help $PD.media_greyscale button \
	"Print / Save Image in Greyscale for Black & White Media."

butt_help $PD.print button \
	"Actually Print / Save Image Now."

butt_help $PD.preview button \
	"Preview Image to be Printed / Saved."

butt_help $PD.cancel button \
	"Cancel Print / Save Image - Lower Print Dialog."

butt_help $PD_P label \
"Progress Indicator - Shows Completed/Remaining Internal Processing."

# print view

butt_help $PV.print button \
	"Actually Print / Save Image Now."

butt_help $PV.refresh button \
	"Refresh Preview to Correspond to Current View Settings."

butt_help $PV.props button \
	"Return to Print Dialog to Set Print / Save Properties."

butt_help $PV.cancel button \
	"Cancel Print / Save Image - Lower Print View Window."

#
# Cursor Stuff...
#

set toplist [ list $MV $VF $PD $PV ]

set save_cursor "none"

#
# Defaults Globals
#

set nodefault_zero_color FALSE
set nodefault_color_range FALSE
set nodefault_grid FALSE
set nodefault_data_mode FALSE
set nodefault_complex_view FALSE
set nodefault_print_command FALSE
set nodefault_progress_ticks TRUE

#
# Actually Bring Up Main MatView Window
#

wm deiconify $MV

#
# Reset Welcome Message
#

# eem puts "done."
# eem flush stdout

setMsg "Welcome to MatView"

update

