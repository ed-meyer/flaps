#
# $Id: procs.tcl,v 1.77 1999/03/11 22:10:57 kohl Exp $
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


# eem puts -nonewline "\[procs.tcl\]"
# eem flush stdout


#
# MatView TCL Procs
#


proc fileHandle { index label menuvar foo } \
{
	#debug_proc fileHandle entry

	global MV

	if { $index == "noload_values" } \
	{
		global $menuvar

		set_opt $index [ set $menuvar ]
	} \
\
	elseif { $index == "open" } \
	{
		lowerMenu $MV.file_menu
		
		matHandle browse
	} \
\
	elseif { $index == "reload" } \
		{ reload_mat } \
\
	elseif { $index == "print" } \
		{ print } \
\
	elseif { $index == "about" } \
		{ about } \
\
	elseif { $index == "quit" } \
		{ quit } \
\
	else \
		{ puts "fileHandle{}:  Unknown Index $index..." }

	#debug_proc fileHandle exit
}


proc about { } \
{
	#debug_proc about entry

	global AM

	global active_fg_color
	global border_space
	global frame_border
	global main_font
	global fg_color
	global bg_color

	set AM .about_matview

	set height 400
	set width 450

	set ck [ winfo exists $AM ]

	if { $ck == 1 } \
	{
		raise $AM
		return
	}

	toplevel $AM
	
	wm geometry $AM ${width}x${height}+100+100

	wm resizable $AM 0 0

	$AM configure -background $bg_color

	button $AM.close -text {Close} -padx 1 -pady 1 \
		-font $main_font \
		-foreground $fg_color -background $bg_color \
		-activeforeground $active_fg_color \
		-command {destroy $AM}

	restrict_bindings $AM.close "Button"

	place $AM.close \
		-x [ expr $width - $border_space ] \
		-y [ expr $height - $border_space ] \
		-anchor se

	update

	set cht [ expr [ above $AM.close ] - ( 4 * $border_space ) ]

	set cwt [ expr $width - ( 4 * $border_space ) ]

	scrollbar $AM.vscroll -orient vert \
		-bd $frame_border -relief sunken \
		-command "$AM.canvas yview"

	restrict_bindings $AM.vscroll "Scrollbar"

	place $AM.vscroll \
		-x $border_space -y $border_space -anchor nw \
		-height $cht -width [ expr 2 * $border_space ]

	scrollbar $AM.hscroll -orient horiz \
		-bd $frame_border -relief sunken \
		-command "$AM.canvas xview"

	restrict_bindings $AM.hscroll "Scrollbar"

	place $AM.hscroll \
		-x [ expr 3 * $border_space ] \
		-y [ expr [ above $AM.close ] - ( 3 * $border_space ) ] \
		-height [ expr 2 * $border_space ] -width $cwt \
		-anchor nw

	canvas $AM.canvas -bd $frame_border -relief sunken -confine 1 \
		-background $bg_color \
		-xscrollcommand "$AM.hscroll set" \
		-yscrollcommand "$AM.vscroll set"

	restrict_bindings $AM.canvas "Canvas"

	place $AM.canvas \
		-x [ expr 3 * $border_space ] -y $border_space \
		-height $cht -width $cwt -anchor nw

	update

	set dims [ fill_about_text $AM.canvas $fg_color $main_font ]

	set tht [ lindex $dims 0 ]
	set twt [ lindex $dims 1 ]

	$AM.canvas configure -scrollregion [ list 0 0 $twt $tht ]

	#debug_proc about exit
}


proc print { } \
{
	#debug_proc print entry

	global MV

	global PD
	global PV

	set ckpd [ winfo ismapped $PD ]

	if { $ckpd == 1 } \
		{ wm withdraw $PD } \
\
	else \
		{ wm deiconify $PD }

	wm withdraw $PV

	lowerMenu $MV.file_menu

	#debug_proc print exit
}


proc printHandle { cmd } \
{
	#debug_proc printHandle entry

	global PD
	global PV

	global print_format
	global print_orient
	global print_media
	global print_dest

	global pv_image_id

	if { $cmd == "dest" } \
	{
		layout_print_dialog
	} \
\
	elseif { $cmd == "file_browse" } \
	{
		set ifile [ $PD.file_entry get ]

		set file [ tk_getSaveFile -initialfile $ifile \
			-parent $PD -title {Enter File for Printing Current View} ]
		
		if { $file != "" } \
		{
			set wd [ pwd ]
			set d [ file dirname $file ]

			if { $d == $wd } \
				{ set file [ file tail $file ] }

			set_print_file $file
		}
	} \
\
	elseif { $cmd == "format" } \
	{
		set file [ $PD.file_entry get ]

		set ext [ file extension $file ]

		if { $ext == ".ps" || $ext == ".ppm" } \
		{
			set root [ file rootname $file ]

			set file "$root.$print_format"

			set_print_file $file
		}
	} \
\
	elseif { $cmd == "print" \
		|| $cmd == "preview" \
		|| $cmd == "preview_print" } \
	{
		if { $pv_image_id == -1 } \
		{
			global PV_C
			global PV_I

			set_cursor busy

			set pv_image_id \
				[ $PV_C create image 0 0 -image $PV_I -anchor nw ]
			
			set_cursor reset
		}

		if { $print_dest == "printer" } \
		{
			set pcmd [ $PD.pcmd_entry get ]

			do_print $cmd printer $print_orient $print_media $pcmd
		} \
\
		elseif { $print_dest == "file" } \
		{
			set file [ $PD.file_entry get ]

			do_print $cmd file $print_orient $print_media \
				$file $print_format
		}

		wm withdraw $PD

		if { $cmd == "preview" } \
		{
			layout_pv_panel

			wm deiconify $PV
		} \
\
		else \
			{ wm withdraw $PV }
	} \
\
	elseif { $cmd == "cancel" } \
	{
		do_print $cmd $print_dest $print_orient $print_media xxx

		wm withdraw $PD
		wm withdraw $PV
	}

	#debug_proc printHandle exit
}


proc set_print_file { file } \
{
	#debug_proc set_print_file entry

	global PD

	$PD.file_entry delete 0 end
	$PD.file_entry insert 0 $file

	set xv [ $PD.file_entry xview ]

	set first [ lindex $xv 0 ]
	set last [ lindex $xv 1 ]

	if { $last != 1 } \
	{
		set diff [ expr 1.0 - $last ]

		set fix [ expr $first + $diff ]

		$PD.file_entry xview moveto $fix
	}

	#debug_proc set_print_file exit
}


proc modeHandle { index label menuvar foo } \
{
	#debug_proc modeHandle entry

	global MV

	if { $index == "set_current" } \
	{
		save_field set

		lowerMenu $MV.mode_menu
	} \
\
	elseif { $index == "free_current" } \
	{
		save_field free

		lowerMenu $MV.mode_menu
	} \
\
	elseif { $index == "save_current" } \
	{
		set file [ tk_getSaveFile -initialfile {savefield} \
			-parent $MV -title {Enter File to Save Current Field} ]
		
		if { $file != "" } \
		{
			save_field save $file
		}

		lowerMenu $MV.mode_menu
	} \
\
	elseif { $index == "load_current" } \
	{
		set file [ tk_getOpenFile \
			-parent $MV -title {Enter File to Load as Current Field} ]
		
		if { $file != "" } \
		{
			save_field load $file
		}

		lowerMenu $MV.mode_menu
	} \
\
	else \
	{
		global data_mode
		global vis_mode

		set_data_mode $data_mode

		set_vis_mode $vis_mode
	}

	#debug_proc modeHandle exit
}


proc reduceHandle { index label menuvar foo } \
{
	#debug_proc reduceHandle entry

	global reduction_function

	set_reduction $reduction_function

	#debug_proc reduceHandle exit
}


proc reduceStatus { status } \
{
	#debug_proc reduceStatus entry

	global MV

	global fg_color
	global bg_color

	if { $status == "none" } \
	{
		set fg $fg_color
		set bg $bg_color
	} \
\
	else \
	{
		set fg $bg_color
		set bg $fg_color
	}

	$MV.rstat_label configure -text "Reduction:  $status" \
		-foreground $fg -background $bg

	#debug_proc reduceStatus exit
}


proc optionHandle { index label menuvar foo } \
{
	#debug_proc optionHandle entry

	if { $index == "auto_color" \
		|| $index == "filter_values" \
		|| $index == "min_max" \
		|| $index == "cumulative" \
		|| $index == "absvals" \
		|| $index == "pattern" } \
	{
		global $menuvar

		set_opt $index [ set $menuvar ]

		layout_main_panel
	} \
\
	elseif { $index == "greyscale" } \
	{
		global print_media
		global $menuvar

		if { [ set $menuvar ] == "ON" } \
			{ set print_media "greyscale" } \
		else \
			{ set print_media "color" }

		set_opt $index [ set $menuvar ]

		layout_main_panel
	} \
\
	elseif { $index == "reset_cumulative" } \
	{
		set_opt $index none
	} \
\
	elseif { $index == "sep_view" } \
	{
		global VF

		global $menuvar

		global resize_mat_lock

		global save_main_height
		global min_main_height
		global matview_height
		global main_height

		set state [ set $menuvar ]

		if { $state == "ON" } \
		{
			global MV
			global VF_C
			global vf_rect_id

			wm deiconify $VF

			# sneak in a fake hide_image...
			$VF_C raise $vf_rect_id

			layout_vf_panel

			set ht [ expr $main_height - $matview_height ]

			if { $ht < $min_main_height } \
				{ set ht $min_main_height }

			set save_main_height $main_height

			set resize_mat_lock 1
			set_geometry -1 $ht -1 -1
			set resize_mat_lock 0

			set matview_height [ expr $matview_height \
				+ $ht - [ winfo height $MV ] ]
		} \
\
		else \
		{
			global MV_C
			global rect_id

			wm withdraw $VF

			# sneak in a fake hide_image...
			$MV_C raise $rect_id

			if { $save_main_height > 0 } \
				{ set ht $save_main_height } \
			else \
				{ set ht [ expr $main_height + $matview_height ] }

			set resize_mat_lock 1
			set_geometry -1 $ht -1 -1
			set resize_mat_lock 0
		}

		update

		set_sep_view $state
	} \
\
	elseif { $index == "transpose" } \
	{
		global $menuvar

		set_opt $index [ set $menuvar ]
	} \
\
	else \
		{ puts "optionHandle{}:  Unknown Option \"$index\"." }

	#debug_proc optionHandle exit
}


proc complexHandle { index label menuvar foo } \
{
	#debug_proc complexHandle entry

	global MV

	global complex_view

	lowerMenu $MV.complex_menu

	$MV.complex configure -text "Complex View: $label"

	set_complex_view $complex_view

	#debug_proc complexHandle exit
}


proc set_complex { onoff } \
{
	#debug_proc complexHandle entry

	global is_complex

	if { $is_complex != $onoff } \
	{
		set is_complex $onoff

		layout_main_panel
	}

	#debug_proc complexHandle exit
}


proc matHandle { op } \
{
	#debug_proc matHandle entry

	global MV

	if { $op == "set" } \
	{
		set mat [ $MV.mat_entry get ]

		set_mat_name $mat
	} \
\
	elseif { $op == "browse" } \
	{
		set file [ tk_getOpenFile -initialfile {matrixfile} \
			-parent $MV -title {Enter Matrix File to Display} ]
		
		if { $file != "" } \
		{
			set wd [ pwd ]
			set d [ file dirname $file ]

			if { $d == $wd } \
				{ set file [ file tail $file ] }

			$MV.mat_entry delete 0 end
			$MV.mat_entry insert 0 $file

			set_mat_name $file
		}
	} \
\
	else \
		{ puts "matHandle{}:  Error Unknown Op = $op" }

	#debug_proc matHandle exit
}


proc handlePlane { } \
{
	#debug_proc handlePlane entry

	global plane_choice

	layout_main_panel

	set_plane_choice $plane_choice

	handleVisRegion none none

	#debug_proc handlePlane exit
}


proc handleVisRegion { current next } \
{
	#debug_proc handleVisRegion entry

	global MV

	set xmin	[ $MV.xmin_entry get ]
	set xmax	[ $MV.xmax_entry get ]
	set xcell	[ $MV.xcell_entry get ]

	set ymin	[ $MV.ymin_entry get ]
	set ymax	[ $MV.ymax_entry get ]
	set ycell	[ $MV.ycell_entry get ]

	set zval	[ $MV.zval_entry get ]

	set_vis_region $xmin $xmax $xcell $ymin $ymax $ycell $zval

	set fw [ focus -displayof $MV ]

	if { $fw == $current } \
		{ focus $next }

	#debug_proc handleVisRegion exit
}


proc handleColor { op { cmd none } { value -1 } { units units } } \
{
	#debug_proc handleColor entry

	global MV

	global color_scroll_value

	if { $op == "set" } \
	{
		set cmin	[ $MV.color_min get ]
		set cmax	[ $MV.color_max get ]

		set changed [ set_color_range $cmin $cmax ]

		if { $changed == 1 } \
			{ draw_color_map }
	} \
\
	elseif { $op == "scroll" } \
	{
		if { $cmd == "moveto" } \
		{
			set color_scroll_value $value
		} \
\
		elseif { $cmd == "scroll" } \
		{
			if { $units == "units" } \
				{ set factor 0.01 } \
			elseif { $units == "pages" } \
				{ set factor 0.10 } \
			else \
				{ set factor 0.01 }

			set color_scroll_value \
				[ expr $color_scroll_value + ( $value * $factor ) ]
		} \
\
		else \
			{ puts "handleColor{}:  unknown scroll cmd=$cmd" }

		if { $color_scroll_value < 0.0 } \
			{ set color_scroll_value 0.0 } \
		elseif { $color_scroll_value > 1.0 } \
			{ set color_scroll_value 1.0 }

		$MV.color_scroll set $color_scroll_value $color_scroll_value

		set changed [ set_color_function $color_scroll_value ]

		if { $changed == 1 } \
			{ draw_color_map }
	} \
\
	elseif { $op == "init" } \
	{
		set cmin	[ $MV.color_min get ]
		set cmax	[ $MV.color_max get ]

		set chclr [ set_color_range $cmin $cmax ]

		set val $cmd

		$MV.color_scroll set $val $val

		set chfn [ set_color_function $val ]

		if { $chclr == 1 || $chfn == 1 } \
			{ draw_color_map }
	} \
\
	else \
		{ puts "handleColor{}:  unknown op=$op" }

	#debug_proc handleColor exit
}


proc gridHandle { } \
{
	#debug_proc gridHandle entry

	global MV

	global grid_status

	if { $grid_status == "ON" } \
	{
		set rows [ $MV.grid_row_entry get ]
		set roff [ $MV.grid_roff_entry get ]
		set cols [ $MV.grid_col_entry get ]
		set coff [ $MV.grid_coff_entry get ]

		set_grid on $rows $roff $cols $coff
	} \
\
	else \
		{ set_grid off -1 -1 -1 -1 }

	#debug_proc gridHandle exit
}


proc handleZeroColor { } \
{
	#debug_proc handleZeroColor entry

	global MV

	global fixed_zero_color

	if { $fixed_zero_color == "ON" } \
	{
		set color [ $MV.zero_entry get ]

		if { $color == "" } \
			{ set_zero_color none } \
		else \
			{ set_zero_color $color }
	} \
\
	else \
		{ set_zero_color none }

	#debug_proc handleZeroColor exit
}


proc set_view_size { top nx ny px py } \
{
	#debug_proc set_view_size entry

	if { $nx < 0 && $ny < 0 } \
		{ $top.region_size configure -text "Current Region Size:" } \
\
	else \
	{
		$top.region_size configure \
			-text "Current Region Size:  $nx x $ny  ($px x $py pixels)"
	}

	#debug_proc set_view_size exit
}


proc clear_image { image wt ht } \
{
	#debug_proc clear_image entry

	$image put white -to 0 0 $wt $ht

	#debug_proc clear_image exit
}


proc query_handle { top canvas x y str } \
{
	#debug_proc query_handle entry

	global border_space
	global row_height
	global col_width

	if { $x != -1 && $y != -1 } \
	{
		# Set Query String

		$top.query.label configure -text $str

		# Calc Query Pop-Up Width

		set len [ string length $str ]

		set ht [ expr $row_height + ( 2 * $border_space ) ]
		set wt [ expr ( $len * $col_width ) + ( 2 * $border_space ) ]

		# Calc Pop-Up Coords

		set cx [ expr $x + [ winfo x $canvas ] + $border_space ]
		set cy [ expr $y + [ winfo y $canvas ] + $border_space ]

		# Manual Check-In-Main...

		set ckx [ expr $cx + $wt ]
		set twt [ winfo width $top ]

		if { $ckx > $twt } \
		{
			set cx [ expr $twt - $wt ]

			if { $cx < 0 } { set cx 0 }
		}

		set cky [ expr $cy + $ht ]
		set tht [ winfo height $top ]

		if { $cky > $tht } \
		{
			set cy [ expr $tht - $ht ]

			if { $cy < 0 } { set cy 0 }
		}

		place $top.query -x $cx -y $cy -height $ht -width $wt \
			-anchor nw
	} \
\
	else \
	{
		place forget $top.query
	}

	update

	#debug_proc query_handle exit
}


#
# User Default Procs
#


proc default_complex_view { view } \
{
	#debug_proc default_complex_view entry

	global MV

	global nodefault_complex_view
	global complex_view

	if { $nodefault_complex_view == "TRUE" } \
		{ return }

	set complex_view $view

	if { $view == "re" } \
		{ set label "Re" } \
	elseif { $view == "im" } \
		{ set label "Im" } \
	elseif { $view == "magn" } \
		{ set label "Magn" } \
	elseif { $view == "phase" } \
		{ set label "Phase" } \

	$MV.complex configure -text "Complex View: $label"

	set_complex_view $complex_view

	#debug_proc default_complex_view exit
}


proc default_zero_color { onoff color } \
{
	#debug_proc default_zero_color entry

	global MV

	global nodefault_zero_color
	global fixed_zero_color

	if { $nodefault_zero_color == "TRUE" } \
		{ return }

	$MV.zero_entry delete 0 end
	$MV.zero_entry insert 0 $color

	set fixed_zero_color $onoff

	handleZeroColor

	#debug_proc default_zero_color exit
}


proc default_color_range { min max } \
{
	#debug_proc default_color_range entry

	global MV

	global nodefault_color_range

	if { $nodefault_color_range == "TRUE" } \
		{ return }

	$MV.color_min delete 0 end
	$MV.color_min insert 0 $min

	$MV.color_max delete 0 end
	$MV.color_max insert 0 $max

	handleColor set

	#debug_proc default_color_range exit
}


proc default_grid { onoff nrow roff ncol coff } \
{
	#debug_proc default_grid entry

	global MV

	global nodefault_grid
	global grid_status

	if { $nodefault_grid == "TRUE" } \
		{ return }

	$MV.grid_row_entry delete 0 end
	$MV.grid_row_entry insert 0 $nrow

	$MV.grid_roff_entry delete 0 end
	$MV.grid_roff_entry insert 0 $roff

	$MV.grid_col_entry delete 0 end
	$MV.grid_col_entry insert 0 $ncol

	$MV.grid_coff_entry delete 0 end
	$MV.grid_coff_entry insert 0 $coff

	set grid_status $onoff

	gridHandle

	#debug_proc default_grid exit
}


proc default_data_mode { mode } \
{
	#debug_proc default_data_mode entry

	global nodefault_data_mode
	global data_mode

	if { $nodefault_data_mode == "TRUE" } \
		{ return }

	set data_mode $mode

	set_data_mode $data_mode

	#debug_proc default_data_mode exit
}


proc default_print_command { pcmd } \
{
	#debug_proc default_print_command entry

	global PD

	global nodefault_print_command

	if { $nodefault_print_command == "TRUE" } \
		{ return }

	$PD.pcmd_entry delete 0 end
	$PD.pcmd_entry insert 0 $pcmd

	printHandle pcmd

	#debug_proc default_print_command exit
}


proc network_optimize { } \
{
	#debug_proc network_optimize entry

	global network_display

	set network_display TRUE

	default_interface_param progress_ticks 10

	#debug_proc network_optimize exit
}


proc default_interface_param { what value } \
{
	#debug_proc default_interface_param entry

	global nodefault_progress_ticks

	if { $what == "progress_ticks" \
			&& $nodefault_progress_ticks == "TRUE" } \
		{ return }

	interfaceParamHandle $what $value

	#debug_proc default_interface_param exit
}


#
# Resize & Layout Procs
#


proc resize_vf_panel { } \
{
	#debug_proc resize_vf_panel entry

	global VF

	global resize_mat_lock

	global view_height
	global view_width

	global use_sep_view

	set ht [ winfo height $VF ]
	set wt [ winfo width $VF ]

	if { $view_width != $wt || $view_height != $ht } \
	{
		set view_height $ht
		set view_width $wt

		layout_vf_panel

		if { $resize_mat_lock == 0 && $use_sep_view == "ON" } \
			{ resize_mat }
	}

	#debug_proc resize_vf_panel exit
}


proc resize_print_dialog { } \
{
	#debug_proc resize_print_dialog entry

	global PD

	global pd_height
	global pd_width

	set ht [ winfo height $PD ]
	set wt [ winfo width $PD ]

	if { $pd_width != $wt || $pd_height != $ht } \
	{
		set pd_height $ht
		set pd_width $wt

		layout_print_dialog
	}

	#debug_proc resize_print_dialog exit
}


proc resize_main_panel { } \
{
	#debug_proc resize_main_panel entry

	global MV

	global resize_mat_lock

	global main_height
	global main_width

	global use_sep_view

	set ht [ winfo height $MV ]
	set wt [ winfo width $MV ]

	if { $main_width != $wt || $main_height != $ht } \
	{
		set main_height $ht
		set main_width $wt

		layout_main_panel

		if { $resize_mat_lock == 0 && $use_sep_view == "OFF" } \
			{ resize_mat }
	}

	#debug_proc resize_main_panel exit
}


proc layout_vf_panel { } \
{
	#debug_proc layout_vf_panel entry

	global view_field_height
	global view_field_width
	global view_height
	global view_width

	global progress_height
	global progress_width

	global border_space
	global frame_border
	global row_height
	global col_width

	global use_sep_view

	global VF
	global VF_C
	global VF_P

	global MV_I

	#
	# View Region Size Display
	#

	set y [ expr $view_height - $frame_border ]

	place $VF.region_size \
		-x $border_space -y $y -anchor sw \
		-height $progress_height

	#
	# Progress Bar
	#

	place $VF_P \
		-x [ expr $view_width - $border_space ] -y $y -anchor se \
		-width $progress_width -height $progress_height

	update

	#
	# Canvas
	#

	set view_field_height [ expr $view_height - $border_space \
		- $progress_height - ( 2 * $frame_border ) ]
	set view_field_width [ expr $view_width - ( 2 * $border_space ) ]

	place $VF_C -x $border_space -y $border_space \
		-width $view_field_width -height $view_field_height

	if { $use_sep_view == "ON" } \
	{
		$MV_I configure \
			-width $view_field_width -height $view_field_height
	}

	update

	#debug_proc layout_vf_panel exit
}


proc layout_print_dialog { } \
{
	#debug_proc layout_print_dialog entry

	global pd_height
	global pd_width

	global progress_height
	global progress_width

	global border_space
	global frame_border
	global row_height
	global col_width

	global print_format
	global print_orient
	global print_media
	global print_dest

	global PD
	global PD_P

	#
	# Destination
	#

	place $PD.dest_label \
		-x $border_space -y $border_space -anchor nw
	
	update

	place $PD.dest_printer \
		-x [ expr 2 * $border_space ] \
		-y [ below $PD.dest_label ] \
		-anchor nw

	update

	place $PD.dest_file \
		-x [ expr 2 * $border_space ] \
		-y [ below $PD.dest_printer ] \
		-anchor nw

	update

	#
	# Orientation
	#

	place $PD.orient_label \
		-x [ expr $pd_width / 3 ] -y $border_space -anchor nw

	update

	set x [ expr [ left $PD.orient_label ] + $border_space ]

	place $PD.orient_portrait \
		-x $x -y [ below $PD.orient_label ] -anchor nw

	update

	place $PD.orient_landscape \
		-x $x -y [ below $PD.orient_portrait ] -anchor nw

	update

	#
	# Media
	#

	place $PD.media_label \
		-x [ expr 2 * $pd_width / 3 ] -y $border_space -anchor nw

	update

	set x [ expr [ left $PD.media_label ] + $border_space ]

	place $PD.media_color \
		-x $x -y [ below $PD.media_label ] -anchor nw

	update

	place $PD.media_greyscale \
		-x $x -y [ below $PD.media_color ] -anchor nw

	update

	#
	# Printer Command / File Interface
	#

	set y [ expr [ below $PD.dest_file ] + $border_space \
		+ ( $row_height / 2 ) ]

	if { $print_dest == "printer" } \
	{
		place $PD.pcmd_label -x $border_space -y $y -anchor w

		update

		set x [ right $PD.pcmd_label ]

		place $PD.pcmd_entry -x $x -y $y -anchor w \
			-width [ expr $pd_width - $x - $border_space ]

		place forget $PD.file_label
		place forget $PD.file_entry
		place forget $PD.file_browse

		place forget $PD.format_label
		place forget $PD.format_ps
		place forget $PD.format_ppm

		update

		set y [ expr [ below $PD.pcmd_entry ] + $border_space ]
	} \
\
	else \
	{
		place $PD.file_label -x $border_space -y $y -anchor w

		place $PD.file_browse \
			-x [ expr $pd_width - $border_space ] -y $y -anchor e

		update

		set x1 [ right $PD.file_label ]
		set x2 [ expr [ left $PD.file_browse ] - $frame_border ]

		place $PD.file_entry -x $x1 -y $y -anchor w \
			-width [ expr $x2 - $x1 ]

		update

		set y [ expr [ below $PD.file_entry ] + $border_space \
			+ ( $row_height / 2 ) ]

		place $PD.format_label -x $border_space -y $y -anchor w

		update

		place $PD.format_ps \
			-x [ right $PD.format_label ] -y $y -anchor w

		update

		place $PD.format_ppm\
			-x [ right $PD.format_ps ] -y $y -anchor w

		place forget $PD.pcmd_label
		place forget $PD.pcmd_entry

		update

		set y [ expr [ below $PD.format_ppm ] + $border_space ]
	}

	#
	# Print / Cancel Buttons
	#

	place $PD.cancel \
		-x [ expr $pd_width - $border_space ] -y $y -anchor ne

	update

	place $PD.preview \
		-x [ expr [ left $PD.cancel ] - $frame_border ] -y $y -anchor ne

	update

	place $PD.print \
		-x [ expr [ left $PD.preview ] - $frame_border ] \
		-y $y -anchor ne

	update

	#
	# Progress Bar
	#

	place $PD_P \
		-x $border_space -y [ below $PD.print ] -anchor sw \
		-width $progress_width -height $progress_height

	update

	#
	# Resize Window
	#

	set pd_height [ expr [ below $PD.print ] + $border_space ]

	wm geometry $PD "${pd_width}x${pd_height}"

	update

	#debug_proc layout_print_dialog exit
}


proc layout_pv_panel { } \
{
	#debug_proc layout_pv_panel entry

	global pview_field_height
	global pview_field_width
	global pview_height
	global pview_width

	global border_space
	global frame_border
	global row_height
	global col_width

	global PV
	global PV_C
	global PV_I

	#
	# Determine Preview Window Size
	#

	set butt_ht [ expr $row_height + ( 2 * $frame_border ) ]

	set pview_height [ expr $pview_field_height \
		+ $border_space + ( 2 * $frame_border ) + $butt_ht ]
	
	set pview_width [ expr $pview_field_width + ( 2 * $border_space ) ]

	# minimum width for buttons... :-)

	if { $pview_width < 200 } \
		{ set pview_width 200 }

	wm geometry $PV "${pview_width}x${pview_height}"

	#
	# Bottom Buttons
	#

	set y [ expr $pview_height - $frame_border ]

	place $PV.cancel \
		-x [ expr $pview_width - $border_space ] -y $y -anchor se \
		-height $butt_ht 

	update

	place $PV.props -x [ left $PV.cancel ] -y $y -anchor se \
		-height $butt_ht 

	update

	place $PV.refresh -x [ left $PV.props ] -y $y -anchor se \
		-height $butt_ht 

	update

	place $PV.print -x [ left $PV.refresh ] -y $y -anchor se \
		-height $butt_ht 

	update

	#
	# Canvas
	#

	place $PV_C -x $border_space -y $border_space \
		-width $pview_field_width -height $pview_field_height

	$PV_I configure \
		-width $pview_field_width -height $pview_field_height

	update

	#debug_proc layout_pv_panel exit
}


proc layout_main_panel { } \
{
	#debug_proc layout_main_panel entry

	global matview_height
	global matview_width
	global thumb_height
	global thumb_width
	global main_height
	global main_width
	global color_size

	global border_space
	global frame_border
	global row_height
	global col_width

	global plane_choice

	global state_cumulative
	global state_min_max

	global use_sep_view

	global is_complex

	global data_mode

	global MV
	global MV_M
	global MV_C
	global MV_I
	global MV_N
	global MV_G
	global MV_P

	global MV_TC
	global MV_TI

	global progress_height
	global progress_width

	global navigate_size

	global rstat_width

	#
	# Place Messages
	#

	set msg_y [ expr [ below $MV.title ] + $border_space ]

	set msg_wt [ expr $main_width - (2 * $border_space) ]

	place $MV.message -x $border_space -y $msg_y -width $msg_wt \
		-anchor nw

	update

	place $MV.help_msg -x $border_space -y [ below $MV.message] \
		-width $msg_wt -anchor nw

	update

	#
	# Main Menu Buttons
	#

	set y [ expr [ below $MV.help_msg ] + $border_space ]

	place $MV.file -x $border_space -y $y -anchor nw

	update

	place $MV.modes -x [ right $MV.file ] -y $y -anchor nw

	update

	place $MV.reduce -x [ right $MV.modes ] -y $y -anchor nw

	update

	place $MV.options -x [ right $MV.reduce ] -y $y -anchor nw

	#
	# Reduction Status
	#

	place $MV.rstat_label \
		-x [ expr $main_width - $border_space ] -y $y -anchor ne \
		-width $rstat_width

	update

	#
	# Complex View Types
	#

	if { $is_complex == "ON" } \
	{
		set cplx_x [ expr [ right $MV.options ] \
			+ ( ( [ left $MV.rstat_label ] - [ right $MV.options ] ) \
				/ 2 ) ]

		place $MV.complex -x $cplx_x -y $y -anchor n
	} \
\
	else \
		{ place forget $MV.complex }

	#
	# Layout Main Panel
	#

	set y [ expr [ below $MV.file ] + $border_space ]

	place $MV.mat -x $border_space -y $y -anchor nw

	update

	place $MV.mat_entry -x [ right $MV.mat ] -y $y -anchor nw \
		-width [ expr $main_width - [ winfo width $MV.mat ] \
			- ( 2 * $border_space ) ] \
		-height [ winfo height $MV.mat ]

	update

	set y [ expr [ below $MV.mat_entry ] + $border_space ]

	if { $data_mode == "field" } \
	{
		place $MV.plane_label -x $border_space -y $y -anchor nw

		update

		set y [ expr $y + ( [ winfo height $MV.plane_label ] / 2 ) ]

		place $MV.xy_plane -x [ right $MV.plane_label ] -y $y -anchor w

		update

		place $MV.yz_plane \
			-x [ expr $border_space + [ right $MV.xy_plane ] ] -y $y \
			-anchor w

		update

		place $MV.xz_plane \
			-x [ expr $border_space + [ right $MV.yz_plane ] ] -y $y \
			-anchor w

		update

		set y [ expr [ below $MV.xz_plane ] + ( $border_space / 2 ) ]
	} \
\
	else \
	{
		place forget $MV.plane_label

		place forget $MV.xy_plane
		place forget $MV.yz_plane
		place forget $MV.xz_plane
	}

	set wt [ expr ( $main_width - ( 2 * $border_space ) ) / 2 ]

	if { $data_mode == "field" } \
	{
		if { $plane_choice == "XY" } \
		{
			$MV.xrange configure -text "X Axis Range:"
			$MV.yrange configure -text "Y Axis Range:"

			$MV.zval configure -text "Z Axis Value:"
		} \
\
		elseif { $plane_choice == "YZ" } \
		{
			$MV.xrange configure -text "Y Axis Range:"
			$MV.yrange configure -text "Z Axis Range:"

			$MV.zval configure -text "X Axis Value:"
		} \
\
		elseif { $plane_choice == "XZ" } \
		{
			$MV.xrange configure -text "X Axis Range:"
			$MV.yrange configure -text "Z Axis Range:"

			$MV.zval configure -text "Y Axis Value:"
		} \
\
		else \
		{
			puts "Error:  Unknown Plane Selection...  Assuming X-Y."

			$MV.xrange configure -text "X Axis Range:"
			$MV.yrange configure -text "Y Axis Range:"

			$MV.zval configure -text "Z Axis Value:"
		}
	} \
\
	else \
	{
		$MV.xrange configure -text "Col Range:"
		$MV.yrange configure -text "Row Range:"
	}

	#
	# Axis Labels
	#

	if { $data_mode == "matrix" } \
	{
		set yy $y
		place $MV.yrange -x $border_space -y $yy -anchor nw
		update

		set xy [ expr [ below $MV.yrange ] + $frame_border ]
		place $MV.xrange -x $border_space -y $xy -anchor nw
		update

		set y [ expr [ below $MV.xrange ] + $border_space ]
	} \
\
	else \
	{
		set xy $y
		place $MV.xrange -x $border_space -y $xy -anchor nw
		update

		set yy [ expr [ below $MV.xrange ] + $frame_border ]
		place $MV.yrange -x $border_space -y $yy -anchor nw
		update

		set y [ expr [ below $MV.yrange ] + $border_space ]
	}

	set xx [ max_len [ right $MV.xrange ] [ right $MV.yrange ] ]

	place $MV.xrange_lb -x $xx -y $xy -anchor nw
	place $MV.yrange_lb -x $xx -y $yy -anchor nw
	update

	set pwt [ winfo width $MV.xrange_lb ]

	set ewt [ expr ( $wt - $xx - ( 3 * $pwt ) ) / 2 ]

	#
	# "X" Entries
	#

	place $MV.xmin_entry -x [ right $MV.xrange_lb ] -y $xy -anchor nw \
		-width $ewt
	update

	place $MV.xrange_cm -x [ right $MV.xmin_entry ] -y $xy -anchor nw
	update

	place $MV.xmax_entry -x [ right $MV.xrange_cm ] -y $xy -anchor nw \
		-width $ewt
	update

	place $MV.xrange_rb -x [ right $MV.xmax_entry ] -y $xy -anchor nw
	update

	place $MV.xcell \
		-x [ expr [ right $MV.xrange_rb ] + $frame_border ] \
		-y $xy -anchor nw
	update

	place $MV.xcell_entry -x [ right $MV.xcell ] -y $xy -anchor nw \
		-width [ expr $ewt / 2 ]
	update

	#
	# "Y" Entries
	#

	place $MV.ymin_entry -x [ right $MV.yrange_lb ] -y $yy -anchor nw \
		-width $ewt
	update

	place $MV.yrange_cm -x [ right $MV.ymin_entry ] -y $yy -anchor nw
	update

	place $MV.ymax_entry -x [ right $MV.yrange_cm ] -y $yy -anchor nw \
		-width $ewt
	update

	place $MV.yrange_rb -x [ right $MV.ymax_entry ] -y $yy -anchor nw
	update

	place $MV.ycell \
		-x [ expr [ right $MV.yrange_rb ] + $frame_border ] \
		-y $yy -anchor nw
	update

	place $MV.ycell_entry -x [ right $MV.ycell ] -y $yy -anchor nw \
		-width [ expr $ewt / 2 ]
	update

	#
	# "Z" Value
	#

	if { $data_mode == "field" } \
	{
		set zy [ winfo y $MV.xz_plane ]

		place $MV.zval \
			-x [ expr [ right $MV.xz_plane ] + ( 2 * $border_space ) ] \
			-y $zy -anchor nw
		update

		place $MV.zval_entry -x [ right $MV.zval ] -y $zy -anchor nw \
			-width $ewt
		update
	} \
\
	else \
	{
		place forget $MV.zval
		place forget $MV.zval_entry
	}

	#
	# Min & Max Stats
	#

	if { $state_min_max == "ON" || $state_cumulative == "ON" } \
	{
		set mmx [ expr [ right $MV.xcell_entry ] + $border_space ]

		if { $data_mode == "matrix" } \
		{
			set miny $yy
			set maxy $xy
		} \
		else \
		{
			set miny $xy
			set maxy $yy
		}

		place $MV.minval -x $mmx -y $miny
		place $MV.maxval -x $mmx -y $maxy
	} \
\
	else \
	{
		place forget $MV.minval
		place forget $MV.maxval
	}

	#
	# Thumb Canvas
	#

	place $MV_TC -x $border_space -y $y \
		-width $thumb_width -height $thumb_height

	$MV_TI configure -width $thumb_width -height $thumb_height

	update

	#
	# Grid Interface
	#

	set xx [ expr [ right $MV_TC ] + $frame_border ]

	place $MV.grid -x $xx -y $y -anchor nw

	update

	place $MV.grid_row \
		-x [ right $MV.grid ] -y $y -anchor nw

	update

	place $MV.grid_row_entry \
		-x [ right $MV.grid_row ] -y $y -anchor nw \
		-width 50

	update

	place $MV.grid_roff \
		-x [ right $MV.grid_row_entry ] -y $y -anchor nw

	update

	place $MV.grid_roff_entry \
		-x [ right $MV.grid_roff ] -y $y -anchor nw \
		-width 50

	update

	place $MV.grid_col \
		-x [ right $MV.grid_roff_entry ] -y $y -anchor nw

	update

	place $MV.grid_col_entry \
		-x [ right $MV.grid_col ] -y $y -anchor nw \
		-width 50

	update

	place $MV.grid_coff \
		-x [ right $MV.grid_col_entry ] -y $y -anchor nw

	update

	place $MV.grid_coff_entry \
		-x [ right $MV.grid_coff ] -y $y -anchor nw \
		-width 50

	#
	# Fixed Zero Color
	#

	place $MV.zero_entry \
		-x [ expr $main_width - $border_space ] -y $y -anchor ne \
		-width 75

	update

	place $MV.zero_color -x [ left $MV.zero_entry ] -y $y -anchor ne

	update

	#
	# Navigate Canvas
	#

	set y [ expr [ below $MV.zero_entry ] + $frame_border ]

	place $MV_N \
		-x $xx -y $y -anchor nw \
		-width $navigate_size -height $navigate_size

	update

	#
	# Color Function Graph Canvas
	#

	place $MV_G \
		-x [ expr $main_width - $border_space ] -y $y -anchor ne \
		-width $navigate_size -height $navigate_size

	update

	#
	# Color Map Canvas
	#

	set yy [ expr $y + [ winfo height $MV_N ] ]

	place $MV.color_min -x [ right $MV_N ] -y $yy -anchor sw -width 100

	place $MV.color_max -x [ left $MV_G ] -y $yy -anchor se -width 100

	update

	set mwt [ expr [ left $MV.color_max ] - [ right $MV.color_min ] ]

	place $MV_M -x [ right $MV.color_min ] -y $yy -anchor sw \
		-width $mwt -height [ winfo height $MV.color_min ]

	update

	place $MV.color_scroll -x [ right $MV_N ] -y $y \
		-width [ expr $mwt + 200 ] \
		-height [ expr [ above $MV_M ] - $y ]

	set color_size [ winfo width $MV_M ]

	draw_color_map

	#
	# View Region Size Display
	#

	set y [ expr $main_height - $frame_border ]

	place $MV.region_size \
		-x $border_space -y $y -anchor sw \
		-height $progress_height

	#
	# Progress Bar
	#

	place $MV_P \
		-x [ expr $main_width - $border_space ] -y $y -anchor se \
		-width $progress_width -height $progress_height

	update

	#
	# Main Canvas
	#

	if { $use_sep_view == "OFF" } \
	{
		set y [ expr [ below $MV_M ] + $border_space ]

		set matview_height [ expr $main_height - $y \
			- $progress_height - ( 2 * $frame_border ) ]
		set matview_width [ expr $main_width - (2 * $border_space) ]

		place $MV_C -x $border_space -y $y \
			-width $matview_width -height $matview_height

		$MV_I configure -width $matview_width -height $matview_height
	} \
\
	else \
		{ place forget $MV_C }

	update

	#debug_proc layout_main_panel exit
}

