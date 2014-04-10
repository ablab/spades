# creating nice HTML reports with summaries about multiple samples

import hypertext as h
from ordereddict import OrderedDict # TorrentServer VM has only python 2.6 :(
import locale # for printing numbers with thousands delimited by commas
import json
import os

import sys
# FIXME: hardcoded version
sys.path.append(os.path.join(os.environ['DIRNAME'], "bin", "quast-2.3"))
from libs.N50 import N50
from libs.fastaparser import get_lengths_from_fastafile, read_fasta

MIRA_LINK = str(h.a("MIRA", target="_blank",
                    href="http://mira-assembler.sourceforge.net"))

SPADES_LINK = str(h.a("SPAdes", target="_blank",
                    href="http://bioinf.spbau.ru/spades"))

ASSEMBLER_LINKS = { 'mira' : MIRA_LINK,
                    'spades' : SPADES_LINK }

DIRNAME = os.environ['TSP_FILEPATH_PLUGIN_DIR']

def _alertAboutAssemblerWebsite(assembler):
    link = ASSEMBLER_LINKS[assembler]
    with h.div(id="alertUser", class_="row-fluid"):
        with h.div(class_="span12"):
            with h.div(class_="alert alert-info"):
                h.button("x", class_="close", data_dismiss="alert",
                         type="button")
                h.UNESCAPED(("Assemblies were performed using %s. If you have "
                            "questions on the quality of the assembly please "
                            "refer to the %s site.") % (link, link))

class Sample(object):
    def __init__(self, info_json_filename):
        with open(info_json_filename, 'r') as info:
            self._info = json.load(info)
        self._name = self._info.get('sampleName') or '.'
        self._metrics = {}
        self._assembly_settings = {}
        self._downloads = {}

        self._has_scaffolds_cache = {}

        params = self._info['params']
        self._assembly_settings['Fraction of Reads Used'] = \
            params['fraction_of_reads']
        self._assembly_settings['RAM'] = params['RAM']
        self._assembly_settings['Barcode Minimum Read Cutoff'] = \
            params['min_reads']

        locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

        mira_info = self._info.get('mira')
        if mira_info:
            self._addAssembler('mira')
            self._loadMiraInfo(mira_info)

        spades_info = self._info.get('spades')
        if spades_info:
            self._addAssembler('spades')
            self._loadSPAdesInfo(spades_info)

    def _addAssembler(self, assembler):
        self._metrics[assembler] = { 'Contig' : OrderedDict(),
                                     'Scaffold' : OrderedDict() }
        self._assembly_settings[assembler] = OrderedDict()
        self._downloads[assembler] = OrderedDict()

    def _loadMiraInfo(self, mira_info):
        settings = self._assembly_settings['mira']
        downloads = self._downloads['mira']
        settings['Reference'] = mira_info['reference']
        settings['MIRA Version'] = mira_info['version']
        settings['Assembly Type'] = mira_info['type']
        self._loadMiraAssemblyInfo(mira_info['info'])
        downloads['Assembled Contigs (FASTA)'] = mira_info['contigs']
        downloads['Coverage Information (WIG)'] = mira_info['wig']
        downloads['Assembly Statistics (TXT)'] = mira_info['info']
        downloads['Mira Log (TXT)'] = mira_info['log']
        if mira_info.get('quastReportDir'):
            report_url = os.path.join(mira_info['quastReportDir'], 'report.html')
            downloads['QUAST report (HTML)'] = report_url

    # fills metrics from an *_info_assembly.txt file
    def _loadMiraAssemblyInfo(self, info_assembly_path):
        with open(info_assembly_path, 'r') as f:
            lines = f.readlines()
            starts = [i-1 for i, line in enumerate(lines) if line[0:4] == "===="]
            sections = [lines[start:end] for start, end \
                            in zip(starts, starts[1:] + [len(lines)])]

            def val(section_index, description):
                for line in sections[section_index]:
                    fields = [f.strip() for f in line.split(':')]
                    if fields[0] == description:
                        return fields[1]
                return None

            def ival(section_index, description):
                value = int(val(section_index, description))
                return locale.format("%d", value, grouping=True)

            def cval(description):
                return (ival(2, description), ival(3, description))

            metrics = self._metrics['mira']['Contig']
            metrics['Assembled Reads'] = (ival(0, 'Num. reads assembled'),)
            metrics['Coverage'] = (val(1, 'Avg. total coverage'),)
            metrics['Consensus Length'] = cval('Total consensus')
            metrics['Number of Contigs'] = cval('Number of contigs')
            metrics['Largest Contig'] = (ival(2, 'Largest contig'),)
            for q in (50, 90, 95):
                metrics['N%d' % q] = cval('N%d contig size' % q)

    def _loadSPAdesInfo(self, spades_info):
        settings = self._assembly_settings['spades']
        downloads = self._downloads['spades']
        settings['SPAdes Version'] = spades_info['version']
        settings['Options'] = spades_info['userOptions']
        downloads['Assembled Contigs (FASTA)'] = spades_info['contigs']
        if self.hasScaffolds('spades'):
            downloads['Assembled Scaffolds (FASTA)'] = spades_info['scaffolds']
        downloads['SPAdes Log (TXT)'] = spades_info['log']
        if spades_info.get('quastReportDir'):
            htmlpath = os.path.join(spades_info['quastReportDir'], 'report.html')
            downloads['QUAST report (HTML)'] = htmlpath

            def fillMetrics(kind, fasta_fn):
                metrics = self._metrics['spades'][kind]
                lengths_all = get_lengths_from_fastafile(fasta_fn)
                lengths_large = [l for l in lengths_all if l >= 500]

                def i(number):
                    return locale.format("%d", int(number), grouping=True)

                def cval(func):
                    return (i(func(lengths_large)), i(func(lengths_all)))

                metrics['Largest ' + kind] = (i(max(lengths_all)), )
                metrics['Total Length'] = cval(sum)
                metrics['Number of ' + kind + 's'] = cval(len)
                for q in [50, 75, 90, 95]:
                    metrics['N%s' % q] = cval(lambda x: N50(x, q))

            fillMetrics('Contig', spades_info['contigs'])
            fillMetrics('Scaffold', spades_info['scaffolds'])
            
    # unique sample ID
    def name(self):
        return self._name

    # ordered dictionary, each record is {metric name : metric value(s)},
    # either a tuple containing single value for all contigs,
    # or a tuple of two values - one for long contigs and another for all of them
    def metrics(self, assembler, kind):
        return self._metrics[assembler][kind]

    def hasScaffolds(self, assembler):
        result = self._has_scaffolds_cache.get(assembler)
        if result is not None:
            return result
        
        result = False
        if assembler == "spades":
            scaffolds_fn = self._info[assembler]['scaffolds']
            for name, seq in read_fasta(scaffolds_fn):
                if 'N' in seq:
                    result = True
                    break

        self._has_scaffolds_cache[assembler] = result
        return result

    # ordered dictionary, each record is {setting name : setting value}
    def assemblySettings(self, assembler):
        return self._assembly_settings[assembler]

    # ordered dictionary, each record is {description : path to file}
    def downloads(self, assembler):
        return self._downloads[assembler]

def _sampleSelector(samples):
    with h.div(class_="row-fluid"):
        with h.div(class_="span12", style="padding-top:20px"):
            with h.div(class_="pull-left"):
                h.strong("View Results")
                h.TEXT(":")
                with h.select(style="width:400px",
                              onchange="showSample(this.value)"):
                    for i, sample in enumerate(samples):
                        if i == 0:
                            h.option(sample.name(),
                                     value=sample.name(), selected="selected")
                        else:
                            h.option(sample.name(), value=sample.name())

def _linksToDownloads(sample, assembler):
    div_id = "collapseDownload" + sample.name()
    with h.div(class_="accordion-group"):
        with h.div(class_="accordion-heading"):
            h.a("Downloads", href="#"+div_id, class_="accordion-toggle",
                data_toggle="collapse")
        with h.div(class_="accordion-body collapse in", id=div_id):
            with h.div(class_="accordion-inner"):
                with h.div(class_="row-fluid", id="downloadPane"):
                    with h.div(class_="row-fluid"):
                        with h.div(class_="span12"):
                            h.p("Download all your assembly result files.")
                            with h.div(style="padding-top:5px"):
                                with h.p():
                                    index = 0
                                    downloads = sample.downloads(assembler)
                                    for descr, url in downloads.items():
                                        rurl = os.path.relpath(url, DIRNAME)
                                        index += 1
                                        if index % 5 != 1:
                                            h.UNESCAPED(("&nbsp;&nbsp;|"
                                                         "&nbsp;&nbsp;"))
                                        h.UNESCAPED("<a target='_blank' "
                                                    "href='%s'>"
                                                    "<i class='icon-download'>"
                                                    "</i>%s</a>" % (rurl, descr))
                                        if index % 5 == 0:
                                            h.br()

def _assemblyStats(sample, assembler):
    def _assemblySettingsTable():
        with h.table(class_="table table-condensed"):
            h.tr(h.th("Parameter"), h.th("Value"))
            for param, value in sample.assemblySettings(assembler).items():
                h.tr(h.td(param), h.td(value))

    def _assemblyMetricsTable(kind):
        with h.table(class_="table table-condensed"):
            with h.tr():
                h.th("Metric")
                with h.th():
                    h.UNESCAPED("Large " + kind + "s (&ge; 500bp)")
                h.th("All " + kind + 's')
            for metric, value in sample.metrics(assembler, kind).items():
                if len(value) == 1:
                    h.tr(h.td(metric), h.td(value[0], colspan='2'))
                else:
                    h.tr(h.td(metric), h.td(value[0]), h.td(value[1]))

    div_id = "collapseStats" + sample.name()
    with h.div(class_="accordion-group"):
        with h.div(class_="accordion-heading"):
            h.a("Assembly Statistics", href="#"+div_id,
                class_="accordion-toggle", data_toggle="collapse")
        with h.div(id=div_id, class_="accordion-body collapse in"):
            with h.div(class_="accordion-inner"):
                with h.div(class_="row-fluid", id="statPane"):
                    with h.div(class_="row-fluid"):
                        with h.div(class_="span12"):
                            if sample.name() != '.':
                                h.p("Assembly summary statistics for ",
                                    sample.name() + ".")
                            with h.div(style="padding-top:5px"):
                                _assemblySettingsTable()
                            with h.div(style="padding-top:10px"):
                                _assemblyMetricsTable("Contig")
                            if sample.hasScaffolds(assembler):
                                with h.div(style="padding-top:10px"):
                                    _assemblyMetricsTable("Scaffold")

def generateReport(samples, assembler):
    css_path = "/pluginMedia/AssemblerPlus/css/"
    css = ["kendo.common.min.css", "kendo.default.min.css", "kendo.ir.css",
           "ir.css", "app.css", "bootstrap.css", "bootstrap-custom.css",
           "bootstrap-select.css"]
    less = ["app.less"]

    js_path = "/pluginMedia/AssemblerPlus/js/"
    js = ["less-1.4.1.min.js", "jquery-1.8.2.min.js",
          "bootstrap-select.min.js", "bootstrap.min.js"]

    with h.html5() as root:
        with h.head():
            h.title("AssemblerPlus Plugin")
            h.meta(name="viewport",
                   content="width=device-width, initial-scale=1.0")
            for css_fn in css:
                h.link(href=css_path+css_fn, rel="stylesheet")
            for less_fn in less:
                h.link(href=css_path+less_fn, rel="stylesheet/less")

        with h.body():
            for js_fn in js:
                h.script('', src=js_path+js_fn)
            with h.script():
                h.UNESCAPED("function showSample(sample) { \n"
                            "  var str = sample.replace(/\./g, '\\\\.') \n"
                            "  $('div.row-fluid.sample').hide() \n"
                            "  $('div#'+str+'.row-fluid.sample').show() \n"
                            "}")
            with h.div(class_="main"):
                with h.div(class_="main-content clearfix"):
                    with h.div(class_="container-fluid"):
                        _alertAboutAssemblerWebsite(assembler)
                        if len(samples) > 1:
                            _sampleSelector(samples)
                        for i, sample in enumerate(samples):
                            display = "display:" + ("none" if i > 0 else "block")
                            with h.div(name=sample.name(), id=sample.name(),
                                       class_="row-fluid sample", style=display):
                                with h.div(class_="span12"):
                                    with h.div(style="padding:0 19px"):
                                        _linksToDownloads(sample, assembler)
                                        with h.div(style="padding-top:10px"):
                                            _assemblyStats(sample, assembler)

    return str(root)
