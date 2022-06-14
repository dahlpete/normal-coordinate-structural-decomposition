#!/bin/bash

for i in {1..6};do
	for j in {0..40};do
		vmd -dispdev text -e nsd_analysis_v2.tcl -args $i $j
		
	done
done
		
