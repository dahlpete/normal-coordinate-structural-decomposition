mol new omcs_2chains_ox_wb_ionized.psf      
mol addfile  /home/pd455/scratch60/omcs_structure/2chains/fully_ox/cold_temp/production/omcs_fully_ox_prod.dcd first 4020 last 8000 step 20 waitfor all

set nf [molinfo top get numframes]
puts $nf

for {set i 0} {$i < $nf} {incr i} {
  #set newframenum [expr {$i + 20}]
  set newframenum $i
  set sel [atomselect top all frame $i]
  $sel writepdb omcs_2chains_ox_snap_$newframenum.pdb
}

exit
