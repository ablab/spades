#!/bin/sh

for f in nw_scr_*
do
    if ! diff nw_scr_orig_NA_32_32.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

for f in sg_scr_*
do
    if ! diff sg_scr_orig_NA_32_32.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

for f in sw_scr_*
do
    if ! diff sw_scr_orig_NA_32_32.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

for f in nw_stats_scr_*
do
    if ! diff nw_stats_scr_orig_NA_32_32.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

for f in sg_stats_scr_*
do
    if ! diff sg_stats_scr_orig_NA_32_32.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

for f in sw_stats_scr_*
do
    if ! diff sw_stats_scr_orig_NA_32_32.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

for f in nw_stats_mch_*
do
    if ! diff nw_stats_mch_orig_NA_32_32.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

for f in sg_stats_mch_*
do
    if ! diff sg_stats_mch_orig_NA_32_32.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

for f in sw_stats_mch_*
do
    if ! diff sw_stats_mch_orig_NA_32_32.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

for f in nw_stats_len_*
do
    if ! diff nw_stats_len_orig_NA_32_32.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

for f in sg_stats_len_*
do
    if ! diff sg_stats_len_orig_NA_32_32.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

for f in sw_stats_len_*
do
    if ! diff sw_stats_len_orig_NA_32_32.txt $f > /dev/null
    then
        echo "$f differs"
    fi
done

