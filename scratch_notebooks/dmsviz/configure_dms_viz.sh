#!/bin/bash

configure-dms-viz format \
   --name "H7 stability Scores 2_6" \
   --input "2-6_entry_func_effects_sequential.csv" \
   --metric "effect" \
   --metric-name "effect" \
   --structure "4r8w_trimer.pdb" \
   --heatmap-limits '-4, 0 , 1' \
   --included-chains 'A B' \
   --excluded-chains 'L H' \
   --output H7_2-6_functional.json



   # configure-dms-viz format \
   # --name "H7 stability Scores" \
   # --input "293T_entry_func_effects_sequential.csv" \
   # --metric "effect" \
   # --metric-name "Effect" \
   # --structure "4r8w_timer.pdb" \
   # --heatmap-limits '-4, 0 , 1' \
   # --included-chains 'A B' \
   # --excluded-chains 'L H' \
   # --output h7_stability_effects.json