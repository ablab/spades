############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import subprocess
import glob
   
### main function ###
def do(report_dict, total_report_horizontal, total_report_vertical, min_contig=0, output_dir=None):

    # suffixes for files with transposed and normal report tables    
    table_ext = '.txt'
    tab_ext   = '.tsv'

    print 'Summarizing...'
    print '  Creating total report...'
    total_report = total_report_horizontal + table_ext
    total_report_tab = total_report_horizontal + tab_ext
    tr_file = open(total_report, 'w')
    tab_file = open(total_report_tab, 'w')

    # calculate columns widthes
    col_widthes = [0 for i in range(len(report_dict['header']))]
    for row in report_dict.keys():            
        for id, value in enumerate(report_dict[row]):
            if len(str(value)) > col_widthes[id]:
                col_widthes[id] = len(str(value))        

    # to avoid confusions:
    if min_contig:
        tr_file.write('Only contigs of length >= ' + str(min_contig) + ' were taken into account\n\n');
    # header
    for id, value in enumerate(report_dict['header']):
        tr_file.write( ' ' + str(value).center(col_widthes[id]) + ' |')
        if id:
            tab_file.write('\t')
        tab_file.write(value)
    tr_file.write('\n')
    tab_file.write('\n')

    # metrics values
    for contig_name in sorted(report_dict.keys()):    
        if contig_name == 'header':
            continue
        for id, value in enumerate(report_dict[contig_name]):
            if id:
                tr_file.write( ' ' + str(value).rjust(col_widthes[id]) + ' |')
                tab_file.write('\t')
            else:
                tr_file.write( ' ' + str(value).ljust(col_widthes[id]) + ' |')
            tab_file.write(str(value))
        tr_file.write('\n')
        tab_file.write('\n')

    tr_file.close()
    tab_file.close()
    print '    Saved to', total_report, 'and', total_report_tab

    print '  Transposed version of total report...'   
    total_report = total_report_vertical + table_ext
    total_report_tab = total_report_vertical + tab_ext
    tr_file = open(total_report, 'w')
    tab_file = open(total_report_tab, 'w') 

    # calculate columns widthes
    col_widthes = [0 for i in range(len(report_dict.keys()))] 
    header_id = 0
    for id, col in enumerate(sorted(report_dict.keys())):
        if col == 'header':
            header_id = id
        for value in report_dict[col]:
            if len(str(value)) > col_widthes[id]:
                col_widthes[id] = len(str(value))

    # to avoid confusions:
    if min_contig:
        tr_file.write('Only contigs of length >= ' + str(min_contig) + ' were taken into account\n\n');

    # filling
    for i in range(len(report_dict['header'])):
        value = report_dict['header'][i]
        tr_file.write( ' ' + str(value).ljust(col_widthes[header_id]) + ' ')
        tab_file.write(str(value) + '\t')
        for id, contig_name in enumerate(sorted(report_dict.iterkeys())):
            if contig_name == 'header':
                continue
            value = report_dict[contig_name][i]
            tr_file.write( ' ' + str(value).ljust(col_widthes[id]) + ' ')
            tab_file.write(str(value) + '\t')
        tr_file.write('\n')
        tab_file.write('\n')
            
    tr_file.close()
    tab_file.close()
    print '    Saved to', total_report, 'and', total_report_tab

    '''
    if all_pdf != None:
        print '  Merging all pdfs...'
        pdfs = ''
        if glob.glob(output_dir + '/basic_stats/*.pdf'):
            pdfs += output_dir + '/basic_stats/*.pdf '
        if glob.glob(output_dir + '/aligned_stats/*.pdf'):
            pdfs += output_dir + '/aligned_stats/*.pdf '
        if glob.glob(output_dir + '/mauve/*.pdf'):
            pdfs += output_dir + '/mauve/*.pdf '
        if pdfs != '':
            subprocess.call(['pdftk ' + pdfs  + ' cat output ' + all_pdf], shell=True)
            print '    Saved to ' + all_pdf
        else:
            print '  There are no pdf files!'    
        
    print '  Done.'
    '''
