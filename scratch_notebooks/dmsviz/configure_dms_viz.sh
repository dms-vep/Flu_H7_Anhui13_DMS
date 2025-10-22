#!/bin/bash

configure-dms-viz format \
   --name "H7 1A8 Ab selection" \
   --input "H7_ab_1A8.csv" \
   --metric "escape_mean" \
   --metric-name "escape_mean" \
   --structure "4r8w_trimer.pdb" \
   --heatmap-limits '-4, 0 , 1' \
   --included-chains 'A B' \
   --excluded-chains 'L H' \
   --output H7_1A8_ab_selection.json



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