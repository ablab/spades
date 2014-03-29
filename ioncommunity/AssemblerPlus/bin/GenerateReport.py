import sys
import os
import json
sys.path.append(os.path.join(os.environ['DIRNAME'], 'lib'))

from ReportGenerator import *

if __name__ == '__main__':
    jsons = sys.argv[1:]
    samples = [Sample(j) for j in jsons]
    with open('startplugin.json', 'r') as f:
        config = json.load(f)
    if config['pluginconfig'].has_key('runMira'):
        mira_report = generateReport(samples, 'mira')
        with open('MIRA_report.html', 'w+') as f:
            f.write(mira_report)
    if config['pluginconfig'].has_key('runSpades'):
        spades_report = generateReport(samples, 'spades')
        with open('SPAdes_report.html', 'w+') as f:
            f.write(spades_report)
    
