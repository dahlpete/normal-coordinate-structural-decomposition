# Normal coordinate structural decomposition
set num_modes 15
set modes_of_interest {1 2 3 8 9 14} 

proc d_oop {ref_molid obs_molid difference_vector} {
	set Fe_atom [lindex [[atomselect $ref_molid "name FE"] get {x y z}] 0]
	set ref [atomselect $ref_molid all]
	$ref moveby [vecscale -1 $Fe_atom]
	set obs [atomselect $obs_molid all]
	$obs moveby [vecscale -1 $Fe_atom]

	set vec1 [lindex [[atomselect $ref_molid "name NA"] get {x y z}] 0]
	set vec2 [lindex [[atomselect $ref_molid "name NB"] get {x y z}] 0]
	set normal [veccross $vec1 $vec2]
	set unit_normal [vecscale [expr 1/[veclength $normal]] $normal]

	# Now compute the out-of-plane distortions
	set count 0
	foreach element $difference_vector {
		set squared_oop  [expr [expr $element * [lindex $unit_normal $count]] * [expr $element * [lindex $unit_normal $count]]]
		#lappend oop $squared_oop
		lappend oop [expr $element * [lindex $unit_normal $count]]
		
		incr count
		if {$count == 3} {
			set count 0
		}
	}
	#set doop [expr sqrt([expr [vecsum $oop] / [veclength $oop]])]
	set doop [veclength $oop]

	return $doop
}
		

proc NSD {difference_vector mode_number} {

	set filename "frequency/normal_mode_${mode_number}.txt"
        set fileID [open $filename]
        set lines [split [read $fileID] "\n"]
        set cnt 0
        foreach line $lines {
		set l1 [lindex [split $line {}] 0]
		if {$cnt > 1 && [llength $l1] > 0} {
			lappend mode_vector $line
		}
		incr cnt 
	}	
	
	#set total_mass [vecsum $masses]
	#set weights [vecscale [expr 1/$total_mass] $masses]
	#puts $weights

	set count 0
	#set mass_idx 0
	#set cnt 0
	#puts "\n"
	#puts "mode ${mode_number}"
	foreach vec_val $mode_vector {
		#puts "[lindex $difference_vector $count]	$vec_val"
		#set mass [lindex $masses $mass_idx]
		#set weight [expr $mass / $total_mass]
		lappend mass_weight_vec $vec_val

		incr count
		#incr cnt
		#if {[expr {fmod($count,3)}] == 0} {
			#incr mass_idx
			#set cnt 0
		#}
	}
	
	#set mass_weight_vec [vecscale [expr 1/[veclength $mass_weight_vec]] $mass_weight_vec]
	#set mass_weighted_length [veclength $mass_weight_vec]
	#set unit_mass_weight_vec [vecscale [expr 1/$mass_weighted_length] $mass_weight_vec]

	#set proj [vecdot [vecscale [expr 1/$total_mass] $difference_vector] $mass_weight_vec]
	set proj [vecdot $difference_vector $mass_weight_vec]
	#set proj [expr abs($proj)]


	return $proj
}

set file_path "porphine_isolation/"
set heme_num [lindex $argv 0]
set snapshot [lindex $argv 1]
#set clas [lindex $argv 1]

# Load the reference structure
mol new copper_reference2.pdb
#xy_plane_transformation {0}
# Load the structure you wish to analyze
mol new "${file_path}heme${heme_num}/heme${heme_num}_snap_${snapshot}.pdb"
#mol new heme_distortions/hhcytc_porph.pdb
#mol new /Users/peterdahl/Documents/BatistaLab/omcs_structure/omcs_2chains/temperature_dependence/fully_reduced/essential_modes/heme_atoms/covarience/center_hemes/distortion_analysis/heme${heme_num}_${clas}.pdb

set nf [molinfo top get numframes]

# Read atom list
set atom_list [open "atom_names.txt"]
set lines [split [read $atom_list] "\n"]
set cnt 0
foreach line $lines {
	lappend atom_names $line
	set l1 [lindex [split $line {}] 0]
	if {$l1 == "C"} {
		lappend mass_list 12.011
	} elseif {$l1 == "N"} {
		lappend mass_list 14.007
	} elseif {$l1 == "F"} {
		lappend mass_list 55.845
	}
	incr cnt
}

set sel1 [atomselect 0 all]
if {$nf == 1} {
	set sel2 [atomselect top all]
	set transformation [measure fit $sel2 $sel1]
	$sel2 move $transformation
	#$sel2 writepdb transformed.pdb

	set atm_cnt 0
	foreach atom $atom_names {
		set l1 [lindex [split $atom {}] 0]
		
		if {[llength $l1] > 0} {
			set ref [atomselect 0 "name $atom"]
			set ref_coords [$ref get {x y z}]

			set obs [atomselect top "name $atom"]
			set obs_coords [$obs get {x y z}]
		
			set mass [lindex $mass_list $atm_cnt]
			set weight [expr sqrt($mass)]
			set coords_diff [vecscale $weight [vecsub [lindex $obs_coords 0] [lindex $ref_coords 0]]]
			puts [veclength [vecsub [lindex $obs_coords 0] [lindex $ref_coords 0]]]

			foreach coord $coords_diff {
				lappend difference_vec $coord
			}
			
			incr atm_cnt
		}		
	}
	
	for {set mode 1} {$mode <= $num_modes} {incr mode} {
		set proj [NSD $difference_vec $mode]
		lappend nsd_list $proj
	}
	set doop [d_oop 0 1 $difference_vec]
}

set doopoutname "doop_heme${heme_num}.txt"
set doopOUT [open $doopoutname a+]
set sadoutname "saddling_heme${heme_num}.txt"
set sadOUT [open $sadoutname a+]
set rufoutname "ruffling_heme${heme_num}.txt"
set rufOUT [open $rufoutname a+]
set domoutname "doming_heme${heme_num}.txt"
set domOUT [open $domoutname a+]
set wavxoutname "wavx_heme${heme_num}.txt"
set wavxOUT [open $wavxoutname a+]
set wavyoutname "wavy_heme${heme_num}.txt"
set wavyOUT [open $wavyoutname a+]
set prooutname "pro_heme${heme_num}.txt"
set proOUT [open $prooutname a+]
#
puts $doopOUT "$doop"
#
set count 0
foreach moi $modes_of_interest {
	set value [lindex $nsd_list [expr $moi-1]]
	puts $value
	if {$count == 0} {
		puts $sadOUT "$value"
	} elseif {$count == 1} {
		puts $rufOUT "$value"
	} elseif {$count == 2} {
		puts $domOUT "$value"
	} elseif {$count == 3} {
		puts $wavxOUT "$value"
	} elseif {$count == 4} {
		puts $wavyOUT "$value"
	} elseif {$count == 5} {
		puts $proOUT "$value"
	}

	lappend important_modes $value
	incr count
}

exit
