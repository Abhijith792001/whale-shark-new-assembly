#!/bin/sh

########################################################################
#
# RepeatMasker can be slow on a whole genome, therefore it is useful
# split a genome into separate files and run RepeatMasker on each file.
# RepeatMasker outputs *.tbl files that describe repeat content.
# This is a simple script to combine RepeatMasker .tbl file output.
#
# Script by Milton Tan
#
# Usage:
### First, cat all of the input tables into a single file.
# cat *.tbl > input.tbl 
### Second, run the script
# RM_tbl_add.sh input.tbl > RM_tbl_sum.out
#
########################################################################

tbl=$1

printf "==================================================\n"
printf "file name: %13s (merged tables)\n" "$tbl"
seqs=`grep 'sequences:' "$tbl" | awk '{sum += $2} END {print sum}'`
printf "sequences: %13d \n" "$seqs"
tl=`grep 'total ' "$tbl" | awk '{sum += $3} END {print sum}'`
tle=`grep 'total ' "$tbl" | tr -d '(' | awk '{sum += $5} END {print sum}'`
printf "total length: %10d bp  (%d by excl N/X-runs)\n" "$tl" "$tle"
#gc = 
printf "GC level: \n"
bm=`grep 'bases ' "$tbl" | awk '{sum += $3} END {print sum}'`
bmp=`echo "$bm $tle" | awk '{printf "%.2f", 100*$1/$2}'`
printf "bases masked: %10d bp (%5.2f %%)\n" "$bm" "$bmp"
# 189045927/2931582331
printf "==================================================\n"
printf "               number of      length   percentage\n"
printf "               elements*    occupied  of sequence\n"
printf "%s\n" "--------------------------------------------------"
sine_l=`grep 'SINEs' "$tbl" | awk '{sum += $3} END {print sum}'`
printf "SINEs:     %12d %12d bp   %5.2f %%\n" `grep 'SINEs' "$tbl" | awk '{sum += $2} END {print sum}'` "$sine_l" `echo "$sine_l $tle" | awk '{printf "%.2f", 100*$1/$2}'`
alu_l=`grep 'ALUs' "$tbl" | awk '{sum += $3} END {print sum}'`
printf "      ALUs %12d %12d bp   %5.2f %%\n" `grep 'ALUs' "$tbl" | awk '{sum += $2} END {print sum}'` "$alu_l" `echo "$alu_l $tle" | awk '{printf "%.2f", 100*$1/$2}'`
mir_l=`grep 'MIRs' "$tbl" | awk '{sum += $3} END {print sum}'`
printf "      MIRs %12d %12d bp   %5.2f %%\n\n" `grep 'MIRs' "$tbl" | awk '{sum += $2} END {print sum}'` "$mir_l" `echo "$mir_l $tle" | awk '{printf "%.2f", 100*$1/$2}'`

line_l=`grep 'LINEs' "$tbl" | awk '{sum += $3} END {print sum}'`
printf "LINEs:     %12d %12d bp   %5.2f %%\n" `grep 'LINEs' "$tbl" | awk '{sum += $2} END {print sum}'` "$line_l" `echo "$line_l $tle" | awk '{printf "%.2f", 100*$1/$2}'`
line1_l=`grep 'LINE1' "$tbl" | awk '{sum += $3} END {print sum}'`
printf "      LINE1 %11d %12d bp   %5.2f %%\n" `grep 'LINE1' "$tbl" | awk '{sum += $2} END {print sum}'` "$line1_l" `echo "$line1_l $tle" | awk '{printf "%.2f", 100*$1/$2}'`
line2_l=`grep 'LINE2' "$tbl" | awk '{sum += $3} END {print sum}'`
printf "      LINE2 %11d %12d bp   %5.2f %%\n" `grep 'LINE2' "$tbl" | awk '{sum += $2} END {print sum}'` "$line2_l" `echo "$line2_l $tle" | awk '{printf "%.2f", 100*$1/$2}'`
l3_l=`grep 'L3' "$tbl" | awk '{sum += $3} END {print sum}'`
printf "      L3/CR1 %10d %12d bp   %5.2f %%\n\n" `grep 'L3' "$tbl" | awk '{sum += $2} END {print sum}'` "$l3_l" `echo "$l3_l $tle" | awk '{printf "%.2f", 100*$1/$2}'`

ltr_l=`grep 'LTR elements' "$tbl" | awk '{sum += $4} END {print sum}'`
printf "LTR elements: %9d %12d bp   %5.2f %%\n" `grep 'LTR elements' "$tbl" | awk '{sum += $3} END {print sum}'` "$ltr_l" `echo "$ltr_l $tle" | awk '{printf "%.2f", 100*$1/$2}'`
ervl_l=`grep 'ERVL ' "$tbl" | awk '{sum += $3} END {print sum}'`
printf "      ERVL %12d %12d bp   %5.2f %%\n" `grep 'ERVL ' "$tbl" | awk '{sum += $2} END {print sum}'` "$ervl_l" `echo "$ervl_l $tle" | awk '{printf "%.2f", 100*$1/$2}'`
ervlM_l=`grep 'ERVL-MaLRs' "$tbl" | awk '{sum += $3} END {print sum}'`
printf "      ERVL-MaLRs %6d %12d bp   %5.2f %%\n" `grep 'ERVL-MaLRs' "$tbl" | awk '{sum += $2} END {print sum}'` "$ervlM_l" `echo "$ervlM_l $tle" | awk '{printf "%.2f", 100*$1/$2}'`
erv1_l=`grep 'ERV_classI ' "$tbl" | awk '{sum += $3} END {print sum}'`
printf "      ERV_classI %6d %12d bp   %5.2f %%\n" `grep 'ERV_classI ' "$tbl" | awk '{sum += $2} END {print sum}'` "$erv1_l" `echo "$erv1_l $tle" | awk '{printf "%.2f", 100*$1/$2}'`
erv2_l=`grep 'ERV_classII' "$tbl" | awk '{sum += $3} END {print sum}'`
printf "      ERV_classII %5d %12d bp   %5.2f %%\n\n" `grep 'ERV_classII' "$tbl" | awk '{sum += $2} END {print sum}'` "$erv2_l" `echo "$erv2_l $tle" | awk '{printf "%.2f", 100*$1/$2}'`

dna_l=`grep 'DNA elements' "$tbl" | awk '{sum += $4} END {print sum}'`
printf "DNA elements: %9d %12d bp   %5.2f %%\n" `grep 'DNA elements' "$tbl" | awk '{sum += $3} END {print sum}'` "$dna_l" `echo "$dna_l $tle" | awk '{printf "%.2f", 100*$1/$2}'`
hAT_l=`grep 'hAT' "$tbl" | awk '{sum += $3} END {print sum}'`
printf "      hAT-Charlie %5d %12d bp   %5.2f %%\n" `grep 'hAT' "$tbl" | awk '{sum += $2} END {print sum}'` "$hAT_l" `echo "$hAT_l $tle" | awk '{printf "%.2f", 100*$1/$2}'`
TcM_l=`grep 'hAT' "$tbl" | awk '{sum += $3} END {print sum}'`
printf "      TcMar-Tigger %4d %12d bp   %5.2f %%\n\n" `grep 'TcM' "$tbl" | awk '{sum += $2} END {print sum}'` "$TcM_l" `echo "$TcM_l $tle" | awk '{printf "%.2f", 100*$1/$2}'`

unc_l=`grep 'Unclassified' "$tbl" | awk '{sum += $3} END {print sum}'`
printf "Unclassified: %9d %12d bp   %5.2f %%\n\n" `grep 'Unclassified' "$tbl" | awk '{sum += $2} END {print sum}'` "$unc_l" `echo "$unc_l $tle" | awk '{printf "%.2f", 100*$1/$2}'`

ir_l=`grep 'interspersed' "$tbl" | awk '{sum += $4} END {print sum}'`
printf "Total interspersed repeats: %8d bp   %5.2f %%\n\n\n" "$ir_l" `echo "$ir_l $tle" | awk '{printf "%.2f", 100*$1/$2}'`

smr_l=`grep 'Small' "$tbl" | awk '{sum += $4} END {print sum}'`
printf "Small RNA: %12d %12d bp   %5.2f %%\n" `grep 'Small' "$tbl" | awk '{sum += $3} END {print sum}'` "$smr_l" `echo "$smr_l $tle" | awk '{printf "%.2f", 100*$1/$2}'`

sat_l=`grep 'Satellites' "$tbl" | awk '{sum += $3} END {print sum}'`
printf "Satellites: %11d %12d bp   %5.2f %%\n" `grep 'Satellites' "$tbl" | awk '{sum += $2} END {print sum}'` "$sat_l" `echo "$sat_l $tle" | awk '{printf "%.2f", 100*$1/$2}'`
sim_l=`grep 'Simple' "$tbl" | awk '{sum += $4} END {print sum}'`
printf "Simple repeats: %7d %12d bp   %5.2f %%\n" `grep 'Simple' "$tbl" | awk '{sum += $3} END {print sum}'` "$sim_l" `echo "$sim_l $tle" | awk '{printf "%.2f", 100*$1/$2}'`
low_l=`grep 'Low' "$tbl" | awk '{sum += $4} END {print sum}'`
printf "Low complexity: %7d %12d bp   %5.2f %%\n" `grep 'Low' "$tbl" | awk '{sum += $3} END {print sum}'` "$low_l" `echo "$low_l $tle" | awk '{printf "%.2f", 100*$1/$2}'`


printf "==================================================\n\n"
printf "* most repeats fragmented by insertions or deletions\n"
printf "  have been counted as one element\n\n"
printf "This is a merged file. For details of each run, consult\n"
printf "individual tbl output files.\n"