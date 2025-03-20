#!/bin/sh
# the next line restarts using wish \
# the next line restarts using wish \
    unset TCL_LIBRARY; unset TK_LIBRARY; export EPICS_CA_ADDR_LIST="129.57.255.4 129.57.255.5"; exec wish "$0" "$@"
# See fastdisplay_readme.txt for details.

# if there were a README file that would be nice
# for now settle for my annotations (frw)
# wasn't me
# didn't write this
# but I changed this code to allow us to normalize the rates by the current
# the code reads the file scalers_future.tex and extracts the lines that specify
# scaler definitions into a list  the line number in the list corresponds to the
# scaler channel -- offset by 1
# SOMEHOW the program  scaler.so  interfaces with the scaler server and gets
# scaler values in a list matching the list extracted from the tex file
# when updating these lines are matched up and the values written to screen
# the structure of the code is such that first display elements are defined and
# variables are initialized.  then the two routines DRAW and UPDATE are defined.
# several variables are shared by declaring them global
#
#
# Changes: 
# - using catch for scaler calls to improve crash resistance
# - redesign of topbar   
# - introduced blinking circle for life sign   
#       Marko, Oct 2001 
#
# - using BCM2 for clock (BCM1 dead)    
#       Marko Oct 2001
#
# - scaler map is read in once during startup now (not for every update as previously)
# - hopefully fixed the bug which made it slower and s l o w e r  and  s  l  o  w  e  r  and  s   l    o
#       Marko, Nov 2001
#
#
#
#
## Get files
regexp {(.*/)*} $argv0 string path
append string scaler.so
load $string
regexp {(.*/)*} $argv0 string path
append string epics.so
load $string

#set x [epics get hac_bcm_average]
#puts $x
#if {$argc > 0} {
#    set dname [lindex $argv 0]
#    if {$dname == "enge"} {
#	set detector 1
#    } elseif {$dname == "hks"} {
#	set detector 0
#    } else {
#	puts "unknown detector $dname\n";
#    }
#} else {
#    # Default to HKS
#    set detector 0
#}

set detector 0

# file that maps scaler channels to detectors
if {$detector == 0} {
    set mapfile BIGHAND.screen
    set layers {1 2 3 4 5 6 7 8 9}
    set layernames {"V1" "V2" "N1" "N2" "N3" "N4" "N5" "N6" "N7"}
} else {
    set mapfile ENGE.screen
    set layers {1 2 3}
    set layernames {"EH1" "EH2" "EH3"}
}

# Startup sample rate in seconds (needs to be between 2 and 60 ) was 4
set waittime 3000

# Startup normalization (0 = raw rates, 1 = normalized by beam current)
set normalize 0

# clock index, offset, gain
#BCM2 used for Q2=1
set doublecheckchan 166
set difference($doublecheckchan) 100
set oldcheck -100
set doubletol 0.05
#Oct BCM2 calibration

set clockoffset 253730.
set clockgain 0.611463e-3

#BCM1 used during Q2=0.5
##set clockchan 346
#July BCM1 calibration
##set clockoffset 249800.
##set clockgain 0.626e-3

# other variables
# Later get this from roc map
set complist {localhost hcdesk3}
set compoffset(localhost) 0
set clockchan(localhost) 834
set compoffset(hcdesk3) 1024
set clockchan(hcdesk3) 290
set ratecheck(localhost) 33
set ratecheck(hcdesk3) 15
set Stop_Go 1
set epics_bcm_name "hac_bcm_average"
set current [epics get $epics_bcm_name]
set dispcurrent 0.
set indexnumber ""
set Label_List {} 
set num 0
set totalchans 0

foreach comp $complist {
    set first_time($comp) 0
    set compcount($comp) [llength [scaler get $comp]]
    set lastchan [expr $compoffset($comp) + $compcount($comp)]
    if {$lastchan > $totalchans} {
	set totalchans $lastchan
    }
}

global Control
set Control($epics_bcm_name) [epics get $epics_bcm_name]

proc update_current {x y z} {
    global Control current epics_bcm_name
    set current $Control($epics_bcm_name)
#    puts $current
}


trace variable Control($epics_bcm_name) w "update_current"

update_current 0 0 0

# startup
proc startup {} {
    global mapfile layers scalchan_list layer_list number_list orientation_list complist Stop_Go current dispcurrent indexnumber Label_List num comp number_min number_max Control layernames

    catch {wm title . "GeN Neutron Detector Rates"}
#    wm minsize . 1250 875
    wm geometry . [expr [llength $layers]*159]x1000
    wm geometry . +0+0

    . configure -background white
    frame .topframe -background lightgrey -height 20

    label .topframe.status -text "Starting" -background green -foreground black
    label .topframe.lifesign -bg green -width 1

    button .topframe.exit -command exit -text "Exit" -relief raised -bg red -fg black

    button .topframe.stop_Go -text "Freeze" -relief raised -width 10 -bg red -fg black -command {
	if {$Stop_Go == 1} then {
	    set Stop_Go 0
	    .topframe.stop_Go configure -text "Resume" -bg green -fg black
	    .topframe.status configure -text "Frozen" -bg red -fg black
	} else {
	    set Stop_Go 1
	    .topframe.stop_Go configure -text "Freeze" -bg red -fg black
	    .topframe.status configure -text "Running: Rates are in kHz" -bg green -fg black
	    Update
	}
    }

    radiobutton .topframe.raw -variable normalize -value 0 -text "RAW" -fg black -bg white 
    radiobutton .topframe.norm -variable normalize -value 1 -text "NORM" -fg black -bg white 

#    label .topframe.time_display -text "Sampling rate (seconds):" -bg lightgrey -fg black
#    entry .topframe.time -textvariable waittime -width 4 -relief sunken -bg white -fg black

    label .topframe.ratechecklabel -text "100kHz Pulsers" -bg yellow
    foreach comp $complist {
	global ${comp}_ratecheck
	set ${comp}_ratecheck 100.0
	entry .topframe.$comp -textvariable ${comp}_ratecheck -bg yellow -relief flat -width 10
    }

    label .topframe.currentlabel -text "Beam Current (uA):" -relief flat -bg lightgrey -fg black
    entry .topframe.current -textvariable dispcurrent -relief flat -width 8 -bg lightgrey -fg black

    canvas .topframe.canvas -width 30 -height 20 -bg lightgrey

    pack .topframe -fill x 
#    pack .topframe.exit .topframe.stop_Go .topframe.raw .topframe.norm .topframe.time_display .topframe.time .topframe.canvas -side left -fill both
    pack .topframe.exit .topframe.stop_Go .topframe.raw .topframe.norm .topframe.canvas -side left -fill both

    pack .topframe.ratechecklabel -side left -fill both
    foreach comp $complist {
	pack .topframe.$comp -side left -fill both
    }

    pack  .topframe.current .topframe.currentlabel -side right -fill both

    pack  .topframe.status -side left -fill both -expand 1

    set i 0
    foreach element $layers {
	frame .$element -background purple -relief groove -borderwidth 4
#	frame .$element.titleframe -height 10
	label .$element.title -text [lindex $layernames $i] -bg lightblue
#	pack .$element.title -in .$element.titleframe
	canvas .$element.$element -width 100 -background white -height 920
#	pack .$element .$element.$element -side left -fill both -expand 1
	pack .$element.title .$element.$element -side top -fill both -expand 1 -in .$element
	pack .$element -side left -fill both -expand 1
	incr i
    }

#    set circle [.8.8 create oval 100 10 30 30 -fill black]
    set circle [.topframe.canvas create oval 7 23 23 7 -fill lightgrey -outline {}]
    
#    blink .topframe.lifesign -bg green white 500
#    blink_item .8.8 $circle -fill grey black 500
#    blink_item .topframe.canvas $circle -fill lightgrey black 500

    update idletasks


# read in configuration from map file
    set fileid [open $mapfile]
    while {[gets $fileid line] >= 0} {
	if {[regexp {^%} $line]} then {continue}
#
#  IT IS IMPORTANT THAT THE VALID ENTRIES THAT CORRESPOND TO EACH
#  DETECTOR IN THE MAP FILE MATCH THIS EXPRESSION!
#
#               1 5B 23
#               Plane CounterOrientation Scalerchan            
#
	if {[regexp {([0-9]+) *([0-9]+)([A-Z]) *([0-9]+)} $line ID layer number orientation scalchan] } then {
#	    puts "$ID|$layer|$number|$orientation|$scalchan"
	   lappend scalchan_list $scalchan
	   set layer_list($scalchan) $layer
	   set number_list($scalchan) $number
	   set orientation_list($scalchan) $orientation  
	}
    }
    close $fileid
    foreach scalchan $scalchan_list {
	set layer $layer_list($scalchan)
	set number $number_list($scalchan)
	set orientation $orientation_list($scalchan)
#	puts "$layer|$number|$orientation"
	if {[catch {set max $number_max($layer)}]} {
	    set max 0
	}
	if {$number > $max} {
	    set number_max($layer) $number
	}
	if {[catch {set min $number_min($layer)}]} {
	    set min 10000
	}
	if {$number < $min} {
	    set number_min($layer) $number
	}
    }

}

proc blink {w option val1 val2 interval} {
    $w config $option $val1
#    update idletasks
    after $interval [list blink $w $option $val2 $val1 $interval]
}

proc blink_item {w item option val1 val2 interval} {
    $w itemconfig $item $option $val1
    after $interval [list blink_item $w $item $option $val2 $val1 $interval]
}

##  Draw runs once for each entry in the source file  
proc Draw {scalchan layer number orientation error minnumber maxnumber}   {
    global ID channel chLetter chNumber comp indexnumber waittime \
	   normalize fileid difference layers num \
	   Label_List  first_time current rate

    incr num
    set vspace [expr {(28 * 32) / ($maxnumber-$minnumber+1)}]
    set X0 75
    #	set Y0 [expr {920 - 28 * $number}]
    set Y0 [expr {468 - ($number-$minnumber - ($maxnumber-$minnumber)/2)*$vspace}]
    set X1 $X0
    set Y1 $Y0
    if {$orientation == "T" || $orientation == "L"} {incr X1 -46} elseif {$orientation == "B" || $orientation == "R"} {incr X1 46}    

	if {!$error} {
#	    set Display_Number [expr round([lindex $difference($comp) [expr $scalchan - 1]] * $TimeConstant) / 1000.0]
	    set rounderr [catch {set Display_Number [expr round($rate([expr $scalchan-1]))/1000.0]}]
	    if {$rounderr} {
		set Display_Number 9999
	    }
	    if {$normalize == 1} then {
		if {$current > 0.0} then {
		    set Display_Number [expr $Display_Number / $current]
		} else {
		    set Display_Number 0.0
		}
	    } 

	    set bg_color white
	    set fg_color black
	   if {$Display_Number >= 1000.0} {
		set bg_color #ffff00000000
		set fg_color white
		set Display_Number [expr round($Display_Number)]
	    } elseif {($Display_Number >= 100.0)&&($Display_Number < 1000.0)} { 
		set bg_color #ffff7f7f7f7f
		set fg_color white
		set Display_Number [format "%5.1f" $Display_Number]
	    } elseif {($Display_Number >= 10.0)&&($Display_Number < 100.0)} { 
		set bg_color #a4a4f5f5eeee
		set fg_color black
		set Display_Number [format "%5.2f" $Display_Number]
	    } elseif {($Display_Number >= 1.0)&&($Display_Number < 10.0)} { 
		set bg_color #0000b8b8a8a8
		set fg_color black
		set Display_Number [format "%5.3f" $Display_Number]
	    } elseif {($Display_Number >= 0.1)&&($Display_Number < 1.0)} { 
		set bg_color #d3d3d3d3d3d3
		set fg_color black
		set Display_Number [expr round($Display_Number * 1000.) / 1000.]
		set Display_Number [format "%5.3f" $Display_Number]
	    } elseif {($Display_Number >= 0.01)&&($Display_Number < 0.1)} { 
		set bg_color #969696969696
		set fg_color white
		set Display_Number [expr round($Display_Number * 1000.) / 1000.]
		set Display_Number [format "%5.3f" $Display_Number]
	    } elseif {($Display_Number >= 0.001)&&($Display_Number < 0.01)} { 
	        set bg_color #646464646464
	        set fg_color white
	        set Display_Number [expr round($Display_Number * 1000.) / 1000.]
	        set Display_Number [format "%5.3f" $Display_Number]
	    } elseif {$Display_Number == 0.00} {
		set bg_color #000
		set fg_color white
		set Display_Number [format "%5.3f" $Display_Number]
	    }
	} else {
		set bg_color #000
		set fg_color white
		set Display_Number ???
	}


	if {$first_time($comp) == 0} then {
	   .$layer.$layer create text $X0 $Y0 -text $number -fill black -font -*-*-bold-*-*--12-*-*-*-*-*-*-*
	    label .$layer.$layer.canvas$num -text $Display_Number -width 6 -font -*-*-*-*-*--16-*-*-*-*-*-*-*  -background $bg_color -fg $fg_color	
	    lappend Label_List ".$layer.$layer.canvas$num"
	    .$layer.$layer create window $X1 $Y1 -window .$layer.$layer.canvas$num
	} else {
	    .$layer.$layer.canvas$num configure -text $Display_Number -background $bg_color -fg $fg_color	
	}
}
## End proc Draw

## Update rates
proc Update {} {

    global ID scalchan_list layer_list number_list orientation_list channel chLetter chNumber comp indexnumber waittime fileid difference layers newstring oldstring first_time complist scalchan num Label_List Stop_Go afterID dispcurrent current argv0 mapfile clockchan clockoffset clockgain number_min number_max epics_bcm_name ratecheck
    global doublecheckchan oldcheck doubletol rate totalchans compoffset
    

    if {$first_time($comp) != 0} then {
	after cancel $afterID
    }


    foreach comp $complist {
	global ${comp}_ratecheck
	
	set cpuoffset $compoffset($comp)
	if {$first_time($comp) == 0} then {
	    set error [catch {set oldstring($comp) [scaler get $comp]}]
	} else {
	    set oldstring($comp) $newstring($comp)
	}
	set error [catch {set newstring($comp) [scaler get $comp]}]
	set difference($comp) {}

	if {$error} {
	    .topframe.status configure -text "Server not responding, waiting..."  -bg red -fg black
	} else { 

	    .topframe.status configure -text "Running: Rates are in kHz"  -bg green -fg black

	    set length [llength $newstring($comp)] 
	    for {set listcounter 0} {$listcounter < $length} {incr listcounter} {
		set b [lindex $newstring($comp) $listcounter]
		set a [lindex $oldstring($comp) $listcounter]
		lappend difference($comp) [expr ($b - $a)]
	    }

	    if {[lindex $difference($comp) [expr $clockchan($comp) - 1]] > 0} then {
		set TimeConstant [expr {60.0 / [lindex $difference($comp) [expr $clockchan($comp) - 1]]}]
	    } else {
		set TimeConstant 0.5
#		if {$comp != "vmec12"} {
#		    .topframe.status configure -text "The clock is not running.  Cannot get rate." -background red -foreground black
#		}
	    }

	    for {set listcounter 0} {$listcounter < $length} {incr listcounter} {
		set rate([expr $listcounter+$cpuoffset]) [expr [lindex $difference($comp) $listcounter]*$TimeConstant]
		
	    }



	}

	if {!$error} {
	}

	# Copy from difference array into rate array, multiplying by the time constant

	set normrate [expr round($rate([expr $ratecheck($comp)+$cpuoffset]))/1000.0]
	set ${comp}_ratecheck $normrate

    }
    if {!$error} {
	#	    set current [expr ([lindex $difference($comp) [expr $clockchan - 1]] * $TimeConstant)]

	#	    set current [expr ($current - $clockoffset) * $clockgain]			
	#	    set current [epics get "IBCxxxxCRCUR2"]
	#	    set current 1
	set current [epics get $epics_bcm_name]
	if {$current >= 1000.0} {
	    set dispcurrent [expr round($current)]
	} elseif {$current >= 100.0} {
	    set dispcurrent [expr round($current * 10.) / 10.]
	} elseif {$current >= 10.0} {
	    set dispcurrent [expr round($current * 100.) / 100.]
	} else {
	    set dispcurrent [expr round($current * 1000.) / 1000.]
	}	    
    } else {
	#	    set current 0.
	#	    set dispcurrent 0.
    }

    foreach scalchan $scalchan_list {
	set ilayer $layer_list($scalchan)
       Draw $scalchan $layer_list($scalchan) $number_list($scalchan) $orientation_list($scalchan) $error $number_min($ilayer) $number_max($ilayer)
#       update
    }
    set num 0


    ##Start Error stuff for time input
    if {![regexp {[0-9]*} $waittime] || $waittime < 2 || $waittime > 60} then {
	set waittime 2
	.topframe.status configure -text "Error:  Rate must be between 2 and 60 seconds" -background red -foreground black
    } 

    update
    set afterID [after [expr $waittime * 1000] {if {$Stop_Go == 1} then {Update}}]
    if {!$error} {
	foreach comp $complist {
	    set first_time($comp) 1
	}
    }
}
##End proc Update

startup
Update    











