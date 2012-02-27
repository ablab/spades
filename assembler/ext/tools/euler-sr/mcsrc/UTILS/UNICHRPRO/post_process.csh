#! /usr/bin/env csh


/home/mchaisso/UTILS/quality_filter.pl traces.fas traces.qual 10 15 > traces.filtered.fas

/home/mchaisso/UTILS/gfquery.csh traces.filtered.fas traces.htang.blat 9997 1200
/home/mchaisso/UTILS/gfquery.csh traces.filtered.fas traces.ucsc.blat 9998 1200

