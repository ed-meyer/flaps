#
# $Id: util.tcl,v 1.7 1998/08/28 17:42:14 kohl Exp $
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

# eem puts -nonewline "\[util.tcl\]"
# eem flush stdout

#
# TCL Utility Procs
#

proc makeMenu { name cmd info_list close leavesubs prefix } \
{
	#debug_proc makeMenu entry

	global active_fg_color
	global selector_color
	global border_space
	global frame_border
	global main_font
	global fg_color

	frame $name -bd $frame_border -relief raised

	restrict_bindings $name "Frame"

	set root [ rootName $name ]

	set mil "${root}_list"

	global $mil

	set $mil ""

	set LAST none

	set width 0

	foreach i $info_list \
	{
		set label	[ lindex $i 0 ]
		set type	[ lindex $i 1 ]
		set index	[ lindex $i 2 ]

		#
		# Create Button / Label
		#

		set fix_label [ strip_label $label ]

		if { $type == "checkbutton" } \
		{
			#
			# Get State From Info List
			#

			set state	[ lindex $i 3 ]
			set help	[ lindex $i 4 ]

			#
			# Create checkbutton state var
			#

			if { $prefix == "none" } \
				{ set prefix $root }

			set menuvar "${prefix}_state_$fix_label"

			global $menuvar

			if { $state == "ON" } \
				{ set $menuvar ON } \
\
			else \
				{ set $menuvar OFF }

			register_menuvar $name $index $menuvar

			#
			# Set checkbutton command / create
			#

			set command [ list $cmd $index $label $menuvar FALSE ]

			set NAME "$name.ckbutt_$fix_label"

			checkbutton $NAME -text $label -command $command \
				-font $main_font -padx 1 -pady 1 \
				-bd 1 -relief raised -anchor nw \
				-onvalue ON -offvalue OFF -variable $menuvar \
				-selectcolor $selector_color -foreground $fg_color \
				-activeforeground $active_fg_color

			restrict_bindings $NAME ""
		} \
\
		elseif { $type == "statebutton" } \
		{
			#
			# Get State From Info List
			#

			set state	[ lindex $i 3 ]
			set menuvar	[ lindex $i 4 ]
			set help	[ lindex $i 5 ]

			global $menuvar

			if { $state == "ON" } \
				{ set $menuvar ON } \
\
			else \
				{ set $menuvar OFF }

			register_menuvar $name $index $menuvar

			#
			# Set checkbutton command / create
			#

			set command [ list $cmd $index $label $menuvar FALSE ]

			set NAME "$name.ckbutt_$fix_label"

			checkbutton $NAME -text $label -command $command \
				-font $main_font -padx 1 -pady 1 \
				-bd 1 -relief raised -anchor nw \
				-onvalue ON -offvalue OFF -variable $menuvar \
				-selectcolor $selector_color -foreground $fg_color \
				-activeforeground $active_fg_color

			restrict_bindings $NAME ""
		} \
\
		elseif { $type == "radiobutton" } \
		{
			set var		[ lindex $i 3 ]
			set help	[ lindex $i 4 ]

			set command [ list $cmd $index $label none FALSE ]

			set NAME "$name.rdbutt_$fix_label"

			radiobutton $NAME -text $label \
				-font $main_font -padx 1 -pady 1 \
				-command $command \
				-value $index -variable $var \
				-bd $frame_border -relief raised -anchor nw \
				-selectcolor $selector_color -foreground $fg_color \
				-activeforeground $active_fg_color

			restrict_bindings $NAME ""
		} \
\
		elseif { $type == "command" } \
		{
			set help	[ lindex $i 3 ]

			set command [ list $cmd $index $label none FALSE ]

			set NAME "$name.butt_$fix_label"

			button $NAME -text $label -command $command \
				-font $main_font -padx 1 -pady 1 \
				-bd 1 -relief raised \
				-foreground $fg_color -activeforeground $active_fg_color

			restrict_bindings $NAME ""
		} \
\
		elseif { $type == "particle" } \
		{
			global particle_radius
			global particle_color

			set help	[ lindex $i 3 ]

			set NAME "$name.part_$fix_label"

			frame $NAME -bd 1 -relief raised

			restrict_bindings $NAME Frame

			radiobutton $NAME.radius -text "R" \
				-font $main_font -padx 1 -pady 1 \
				-value $index \
				-variable particle_radius \
				-bd 2 -relief sunken -anchor nw \
				-selectcolor $selector_color -foreground $fg_color \
				-activeforeground $active_fg_color

			restrict_bindings $NAME.radius ""

			butt_help $NAME.radius button $help

			radiobutton $NAME.color -text "C" \
				-font $main_font -padx 1 -pady 1 \
				-value $index \
				-variable particle_color \
				-bd 2 -relief sunken -anchor nw \
				-selectcolor $selector_color -foreground $fg_color \
				-activeforeground $active_fg_color

			restrict_bindings $NAME.color ""

			butt_help $NAME.color button $help

			label $NAME.label -text $label \
				-foreground $fg_color -font $main_font

			restrict_bindings $NAME.label ""

			place $NAME.radius -x 0 -y 0
			update

			place $NAME.color \
				-x [ expr [ right $NAME.radius ] + $frame_border ] -y 0
			update

			place $NAME.label \
				-x [ expr [ right $NAME.color ] + $frame_border ] -y 0
			update

			place $NAME -x -100 -y -100 \
				-width [ expr [ right $NAME.label ] + $border_space ] \
				-height [ expr [ below $NAME.radius ] + $frame_border ]
		} \
\
		elseif { $type == "submenu" } \
		{
			set focus_entry [ lindex $i 3 ]
			set subs		[ lindex $i 4 ]
			set help		[ lindex $i 5 ]

			set NAME "$name.butt_$fix_label"

			set command [ list raiseSubMenu $index $name $NAME \
				$focus_entry $subs ]

			button $NAME -text $label -command $command \
				-font $main_font -padx 1 -pady 1 \
				-bd 1 -relief raised -foreground $fg_color \
				-activeforeground $active_fg_color

			restrict_bindings $NAME ""
		} \
\
		elseif { $type == "exchange" } \
		{
			set button	[ lindex $i 3 ]
			set parent	[ lindex $i 4 ]
			set subs	[ lindex $i 5 ]
			set help	[ lindex $i 6 ]

			set command [ list exchangeMenu $index $name \
				$button $parent $subs ]

			set NAME "$name.butt_$fix_label"

			button $NAME -text $label -command $command \
				-font $main_font -padx 1 -pady 1 \
				-bd 1 -relief raised -foreground $fg_color \
				-activeforeground $active_fg_color

			restrict_bindings $NAME ""
		} \
\
		elseif { $type == "exchangeIndirect" } \
		{
			set bname	[ lindex $i 3 ]
			set parent	[ lindex $i 4 ]
			set subs	[ lindex $i 5 ]
			set help	[ lindex $i 6 ]

			set command [ list exchangeMenuIndirect $index $name \
				$bname $parent $subs ]

			set NAME "$name.butt_$fix_label"

			button $NAME -text $label -command $command \
				-font $main_font -padx 1 -pady 1 \
				-bd 1 -relief raised -foreground $fg_color \
				-activeforeground $active_fg_color

			restrict_bindings $NAME ""
		} \
\
		else \
		{
			set NAME "$name.label_$fix_label"

			label $NAME -text $label \
				-foreground $fg_color -font $main_font

			restrict_bindings $NAME ""
		}

		if { $type != "label" && $type != "particle" } \
			{ butt_help $NAME button $help }

		#
		# Place Button / Label on Menu Frame
		#

		if { $LAST != "none" } \
			{ set mby [ below $LAST ] } \
\
		else \
			{ set mby $frame_border }

		place $NAME -x $frame_border -y $mby

		update

		#
		# Process Button Width
		#

		set width [ max_len [ winfo width $NAME ] $width ]

		#
		# Save Button Name
		#

		set LAST $NAME

		lappend $mil $NAME
	}

	#
	# Create "Done" Button
	#

	set NAME $name.butt_Done

	if { "$close" == "lower" } \
		{ set cmd [ list lowerMenu $name ] } \
\
	else \
		{ set cmd "$close" }

	button $NAME -text "Done" \
		-font $main_font -padx 1 -pady 1 \
		-bd 1 -relief raised -command $cmd \
		-foreground $fg_color -activeforeground $active_fg_color

	restrict_bindings $NAME ""

	bind $name <Leave> [ list leaveMenu $name %x %y $leavesubs $cmd ]

	if { $LAST != "none" } \
		{ set dy [ below $LAST ] } \
\
	else \
		{ set dy $frame_border }

	place $NAME -x $frame_border -y $dy

	butt_help $NAME close "Lower Menu"

	lappend $mil $NAME

	update

	#
	# Process Final Button Width
	#

	set width [ max_len [ winfo width $NAME ] $width ]

	#
	# Set Main Menu Frame Size
	#

	set menu_height [ expr [ below $NAME ] + 4 ]

	set menu_width [ expr $width + (4 * $frame_border) ]

	$name configure -width $menu_width -height $menu_height

	place forget $name

	#
	# Adjust Menu Button Widths
	#

	foreach i [ set $mil ] \
		{ place $i -width $width }

	#debug_proc makeMenu exit
}

proc raiseMenu { name button parent subs } \
{
	#debug_proc raiseMenu entry

	foreach i $subs \
		{ place forget $i }

	set ckm [ winfo ismapped $name ]

	if { $ckm == 1 } \
		{ place forget $name } \
\
	else \
		{ placeMenu $name $button $parent }

	#debug_proc raiseMenu exit
}

proc raiseMenuIndirect { name bname button parent subs } \
{
	#debug_proc raiseMenuIndirect entry

	global $bname

	foreach i $subs \
		{ place forget $i }

	set ckm [ winfo ismapped $name ]

	if { $ckm == 1 } \
		{ place forget $name } \
\
	else \
	{
		set $bname $button

		placeMenu $name $button $parent
	}

	#debug_proc raiseMenuIndirect exit
}

proc placeMenu { name button parent } \
{
	#debug_proc placeMenu entry

	global main_height
	global main_width

	set px [ winfo rootx $parent ]
	set py [ winfo rooty $parent ]

	set x [ expr [ winfo rootx $button ] - $px ]

	set y [ expr [ winfo rooty $button ] \
		+ [ winfo height $button ] - $py ]

	place $name -x $x -y $y

	raise $name

	update

	check_in_main $name $x $y $main_width $main_height

	update

	#debug_proc placeMenu exit
}

proc replaceMenu { name button parent } \
{
	#debug_proc replaceMenu entry

	set ckup [ winfo ismapped $name ]

	if { $ckup == 1 } \
		{ placeMenu $name $button $parent }

	#debug_proc replaceMenu exit
}

proc raiseSubMenu { submenu menu button focus_entry subs } \
{
	#debug_proc raiseSubMenu entry

	foreach i $subs \
		{ place forget $i }

	set ckup [ winfo ismapped $submenu ]

	if { $ckup == 1 } \
		{ place forget $submenu } \
\
	else \
		{ placeSubMenu $submenu $menu $button $focus_entry }

	#debug_proc raiseSubMenu exit
}

proc placeSubMenu { submenu menu button focus_entry } \
{
	#debug_proc placeSubMenu entry

	global main_height
	global main_width

	set x [ winfo x $menu ]
	set y [ winfo y $menu ]

	set x [ expr $x + [ left $button ] ]
	set y [ expr $y + [ below $button ] ]

	place $submenu -x $x -y $y

	raise $submenu

	update

	check_in_main $submenu $x $y $main_width $main_height

	if { $focus_entry != "none" } \
		{ focus $focus_entry }

	update

	#debug_proc placeSubMenu exit
}

proc replaceSubMenu { submenu menu button focus_entry } \
{
	#debug_proc replaceSubMenu entry

	set ckup [ winfo ismapped $submenu ]

	if { $ckup == 1 } \
		{ placeSubMenu $submenu $menu $button $focus_entry }

	#debug_proc replaceSubMenu exit
}

proc exchangeMenu { submenu menu button parent subs } \
{
	#debug_proc exchangeMenu entry

	set cksm [ winfo ismapped $submenu ]

	if { $cksm == 1 } \
	{
		place forget $submenu

		raiseMenu $menu $button $parent $subs
	} \
\
	else \
	{
		place forget $menu

		raiseMenu $submenu $button $parent $subs
	}

	#debug_proc exchangeMenu exit
}

proc exchangeMenuIndirect { submenu menu bname parent subs } \
{
	#debug_proc exchangeMenuIndirect entry

	global $bname

	set cksm [ winfo ismapped $submenu ]

	if { $cksm == 1 } \
	{
		place forget $submenu

		raiseMenu $menu [ set $bname ] $parent $subs
	} \
\
	else \
	{
		place forget $menu

		raiseMenu $submenu [ set $bname ] $parent $subs
	}

	#debug_proc exchangeMenuIndirect exit
}

proc check_in_main { name x y width height } \
{
	#debug_proc check_in_main entry

	set wt [ winfo width $name ]

	set ckx [ expr $x + $wt ]

	if { $ckx > $width } \
	{
		set x [ expr $width - $wt ]

		if { $x < 0 } { set x 0 }

		place $name -x $x
	}

	set ht [ winfo height $name ]

	set cky [ expr $y + $ht ]

	if { $cky > $height } \
	{
		set y [ expr $height - $ht ]

		if { $y < 0 } { set y 0 }

		place $name -y $y
	}

	#debug_proc check_in_main exit
}

proc leaveMenu { name x y subs cmd } \
{
	#debug_proc leaveMenu entry

	set ht [ winfo height $name ]
	set wt [ winfo width $name ]

	if { $x < 0 || $y < 0 || $x >= $wt || $y >= $ht } \
	{
		foreach w $subs \
		{
			set ckw [ winfo ismapped $w ]

			if { $ckw == 1 } \
				{ return }
		}

		eval $cmd
	}

	#debug_proc leaveMenu exit
}

proc lowerMenu { name } \
{
	#debug_proc lowerMenu entry

	place forget $name

	#debug_proc lowerMenu exit
}

proc setMsg { text } \
{
	#debug_proc setMsg entry

	global MV

	$MV.message configure -text "Status:   $text"

	update

	#debug_proc setMsg exit
}

proc setTmpMsg { hdr text } \
{
	#debug_proc setMsg entry

	global MV

	global msg_list

	$MV.help_msg configure -text "$hdr:   $text"

	set msg_list [ linsert $msg_list 0 [ list 1 $hdr $text ] ]

	update

	#debug_proc setMsg exit
}

proc popMsgs { } \
{
	#debug_proc popMsgs entry

	global MV

	global msg_list

	set cnt 0

	foreach m $msg_list \
	{
		set type [ lindex $m 0 ]

		if { $type == 1 } \
		{
			set msg_list [ lreplace $msg_list $cnt $cnt ]

			set cnt [ expr $cnt - 1 ]
		}

		set cnt [ expr $cnt + 1 ]
	}

	set pop [ lindex $msg_list 0 ]

	set popstr [ lindex $pop 1 ]

	$MV.help_msg configure -text "Help:   $popstr"

	update

	#debug_proc popMsgs exit
}

proc topMsg { } \
{
	#debug_proc topMsg entry

	global msg_list

	set top [ lindex $msg_list 0 ]

	set toptype [ lindex $top 0 ]

	if { $toptype == 0 } \
		{ return [ lindex $top 1 ] } \
\
	else \
		{ return [ lindex $top 2 ] }

	#debug_proc topMsg exit
}

proc set_geometry { width height x y } \
{
	#debug_proc set_geometry entry

	global MV

	global main_height
	global main_width

	if { $width == -1 } \
		{ set width $main_width }

	if { $height == -1 } \
		{ set height $main_height }

	if { $x == -1 } \
		{ set x [ winfo rootx $MV ] }

	if { $y == -1 } \
		{ set y [ winfo rooty $MV ] }

	set geom "[ expr $width ]x[ expr $height ]+[ expr $x ]+[ expr $y ]"

	wm geometry $MV $geom

	update

	#debug_proc set_geometry exit
}

proc rootName { path } \
{
	#debug_proc rootName entry

	set tmp [ split $path . ]

	set num [ llength $tmp ]

	set root [ lindex $tmp [ expr $num - 1 ] ]

	return $root

	#debug_proc rootName exit
}

proc above { win } \
{
	#debug_proc above entry

	set y [ winfo y $win ]

	return $y

	#debug_proc above exit
}

proc below { win } \
{
	#debug_proc below entry

	set y [ expr [ winfo y $win ] + [ winfo height $win ] ]

	return $y

	#debug_proc below exit
}

proc left { win } \
{
	#debug_proc left entry

	set x [ winfo x $win ]

	return $x

	#debug_proc left exit
}

proc right { win } \
{
	#debug_proc right entry

	set x [ expr [ winfo x $win ] + [ winfo width $win ] ]

	return $x

	#debug_proc right exit
}

proc entry_setup { entry cmd } \
{
	bind $entry <Return> $cmd
}

proc butt_help { butt type msg } \
{
	#debug_proc butt_help entry

	bind $butt <Enter> [ list do_butt_help IN $butt $type $msg ]

	bind $butt <Motion> [ list do_butt_help IN $butt $type $msg ]

	if { $type == "button" } \
		{ bind $butt <ButtonPress> [ list $butt invoke ] } \
\
	elseif { $type == "close" } \
		{ bind $butt <ButtonPress> "popMsgs ; $butt invoke" }

	bind $butt <Leave> [ list do_butt_help OUT $butt $type $msg ]

	#debug_proc butt_help exit
}

proc do_butt_help { cmd butt type msg } \
{
	#debug_proc do_butt_help entry

	if { $cmd == "IN" } \
	{
		if { $type == "entry" } \
		{
			set msg [ lindex $msg 1 ]

			focus $butt
		}

		if { $type == "button" } \
			{ $butt configure -state active } \
\
		elseif { $type == "close" } \
			{ $butt configure -state active }

		if { $msg != [ topMsg ] } \
			{ setTmpMsg "Help" $msg }
	} \
\
	elseif { $cmd == "OUT" } \
	{
		popMsgs

		if { $type == "entry" } \
		{
			set cmd [ lindex $msg 0 ]

			eval $cmd
		}

		if { $type == "button" } \
			{ $butt configure -state normal } \
\
		elseif { $type == "close" } \
			{ $butt configure -state normal }
	}

	#debug_proc do_butt_help exit
}

proc set_cursor { op } \
{
	#debug_proc set_cursor entry

	global save_cursor
	global toplist

	if { $op == "busy" } \
	{
		if { $save_cursor == "none" } \
		{
			set win [ lindex $toplist 0 ]

			set save_cursor [ $win cget -cursor ]
		}

		set c [ list watch purple white ]

		foreach win $toplist \
			{ $win configure -cursor $c }
		
		update
	} \
\
	elseif { $op == "reset" } \
	{
		if { $save_cursor != "none" } \
		{
			foreach win $toplist \
				{ $win configure -cursor $save_cursor }

			set save_cursor "none"
		}
	} \
\
	else \
	{
		puts "usage: set_cursor [ busy | reset ]"
	}

	#debug_proc set_cursor exit
}

proc max_width { foo } \
{
	#debug_proc max_width entry

	set wt 0

	foreach w $foo \
	{
		set twt [ winfo width $w ]

		set wt [ max_len $wt $twt ]
	}

	return $wt

	#debug_proc max_width exit
}

proc max_len { s1 s2 } \
{
	#debug_proc max_len entry

	if { $s1 > $s2 } \
		{ return $s1 } \
\
	else \
		{ return $s2 }

	#debug_proc max_len exit
}

proc field_width { txt dcols } \
{
	#debug_proc field_width entry

	global col_width
	global main_font

	set lwt [ lindex [ get_real_text_size $txt $main_font ] 0 ]

	set cwt [ expr ( $dcols + 2 ) * $col_width ]

	return [ max_len $lwt $cwt ]

	#debug_proc field_width exit
}

proc get_base_text_size { } \
{
	#debug_proc get_base_text_size entry

	global MV

	global network_display

	global row_height
	global col_width
	global data_font

	#
	# Create Dummy Label
	#

	set DUMMY $MV.dummy.dummy

	label $DUMMY -text "xxx" -font $data_font

	restrict_bindings $DUMMY ""

	place $DUMMY -x 0 -y 0

	#
	# Run Tests
	#

	if { $network_display == "FALSE" } \
	{
		set test_list \
			"A B C D E F G H I J K L M N O P Q R S T U V W X Y Z"
	} \
\
	else \
		{ set test_list "M W X" }

	set row_height 0
	set col_width 0

	foreach c $test_list \
	{
		set wtht1 [ get_real_text_size "$c" ]

		set wtht2 [ get_real_text_size "$c$c" ]

		set wt1 [ lindex $wtht1 0 ]
		set wt2 [ lindex $wtht2 0 ]

		set wt [ expr $wt2 - $wt1 ]

		set ht [ lindex $wtht1 1 ]

		set col_width [ max_len $wt $col_width ]

		set row_height [ max_len $ht $row_height ]
	}

	destroy $DUMMY

	#debug_proc get_base_text_size exit
}

proc get_real_text_size { str } \
{
	#debug_proc get_real_text_size entry

	global MV

	set DUMMY $MV.dummy.dummy

	$DUMMY configure -text $str

	update

	set ht [ winfo height $DUMMY ]

	set wt [ winfo width $DUMMY ]

	return "$wt $ht"

	#debug_proc get_real_text_size exit
}

proc restrict_bindings { win extra } \
{
	#debug_proc restrict_bindings entry

	bindtags $win "$win $extra"

	#debug_proc restrict_bindings exit
}

proc do_xview { canvas location } \
{
	#debug_proc do_xview entry

	$canvas xview moveto 0.0

	$canvas xview scroll $location units

	#debug_proc do_xview exit
}

proc do_yview { canvas location } \
{
	#debug_proc do_yview entry

	$canvas yview moveto 0.0

	$canvas yview scroll $location units

	#debug_proc do_yview exit
}

proc debug_proc { routine inout } \
{
	global proc_debug

	if { $proc_debug == "TRUE" } \
		{ puts "(proc debug: $routine{} $inout)" }
}

