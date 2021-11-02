set atom_names "CHA CHB CHC CHD C1A C2A C3A C4A C1B C2B C3B C4B C1C C2C C3C C4C C1D C2D C3D C4D NA NB NC ND FE"

#mol new snapshots/omcs_ox_snap_0.pdb 
#for {set i 0} {$i < 200} {incr i} {
#	mol addfile snapshots/omcs_ox_snap_$i.pdb waitfor all
#}

#set nf [molinfo top get numframes]
set start [lindex $argv 0]
set nf 20
set end [expr $start + $nf]

for {set frm $start} {$frm < $end} {incr frm} {
	mol new snapshots/omcs_ox_snap_$frm.pdb

	set heme1 [atomselect top "resid 412 and chain Z and name $atom_names" frame 0]
	set heme2 [atomselect top "resid 410 and chain Z and name $atom_names" frame 0]
	set heme3 [atomselect top "resid 411 and chain Z and name $atom_names" frame 0]
	set heme4b [atomselect top "resid 409 and chain Z and name $atom_names" frame 0]
	set heme4 [atomselect top "resid 409 and chain Y and name $atom_names" frame 0]
	set heme5 [atomselect top "resid 408 and chain Y and name $atom_names" frame 0]
	set heme6 [atomselect top "resid 413 and chain Y and name $atom_names" frame 0]

	$heme1 writepdb heme1/heme1_snap_$frm.pdb
	$heme2 writepdb heme2/heme2_snap_$frm.pdb
	$heme3 writepdb heme3/heme3_snap_$frm.pdb
	$heme4b writepdb heme4b/heme4b_snap_$frm.pdb
	$heme4 writepdb heme4/heme4_snap_$frm.pdb
	$heme5 writepdb heme5/heme5_snap_$frm.pdb
	$heme6 writepdb heme6/heme6_snap_$frm.pdb

}

exit
