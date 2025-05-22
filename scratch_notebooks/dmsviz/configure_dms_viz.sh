#!/bin/bash

configure-dms-viz format \
   --name "H7 Functional Scores" \
   --input "293T_entry_func_effects_sequential.csv" \
   --metric "effect" \
   --metric-name "Effect" \
   --structure "4r8w_timer.pdb" \
   --heatmap-limits '-4, 0 , 1' \
   --included-chains 'A B' \
   --excluded-chains 'L H' \
   --output h7_functional_effects.json