# Classes for parsing *_info_assembly.txt files produced by MIRA
# and creating nice HTML reports with summaries about multiple samples.

import hypertext as h
from ordereddict import OrderedDict # TorrentServer VM has only python 2.6 :(
import locale # for printing numbers with thousands delimited by commas
import json

MIRA_LINK = str(h.a("MIRA", target="_blank",
                    href="http://mira-assembler.sourceforge.net"))

class Sample(object):
    def __init__(self, name):
        self._name = name
        self._metrics = OrderedDict()
        self._assembly_settings = OrderedDict()
        self._downloads = OrderedDict() # TODO

    # fills assembly settings dictionary from 'startplugin.json'
    def loadAssemblySettings(self, startplugin_json_path):
        with open(startplugin_json_path, 'r') as f:
            config = json.load(f)
            params = config['pluginconfig']
            self._assembly_settings['Reference'] = params['agenome']
            self._assembly_settings['MIRA Version'] = params['miraversion']
            self._assembly_settings['Assembly Type'] = params['type']
            self._assembly_settings['Fraction of Reads Used'] = \
                params['fraction_of_reads']
            self._assembly_settings['RAM'] = params['RAM']
            self._assembly_settings['Barcode Minimum Read Cutoff'] = \
                params['min_reads']

    # fills metrics from an *_info_assembly.txt file
    def loadMiraAssemblyInfo(self, info_assembly_path):
        locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

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

            self._metrics['Assembled Reads'] = (ival(0, 'Num. reads assembled'),)
            self._metrics['Coverage'] = (val(1, 'Avg. total coverage'),)
            self._metrics['Consensus Length'] = cval('Total consensus')
            self._metrics['Number of Contigs'] = cval('Number of contigs')
            self._metrics['Largest Contig'] = (ival(2, 'Largest contig'),)
            for q in (50, 90, 95):
                self._metrics['N%d' % q] = cval('N%d contig size' % q)

    # unique sample ID
    def name(self):
        return self._name

    # ordered dictionary, each record is {metric name : metric value(s)},
    # either a tuple containing single value for all contigs,
    # or a tuple of two values - one for long contigs and another for all of them
    def metrics(self):
        return self._metrics

    # ordered dictionary, each record is {setting name : setting value}
    def assemblySettings(self):
        return self._assembly_settings

    # ordered dictionary, each record is {description : path to file}
    def downloads(self):
        return self._downloads

def _alertAboutAssemblerWebsites():
    with h.div(id="alertUser", class_="row-fluid"):
        with h.div(class_="span12"):
            with h.div(class_="alert alert-info"):
                h.button("x", class_="close", data_dismiss="alert",
                         type="button")
                h.UNESCAPED(("Assemblies were performed using %s. If you have "
                            "questions on the quality of the assembly please "
                            "refer to the %s site.") % (MIRA_LINK, MIRA_LINK))

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
                            h.option(value=sample.name(), selected=selected)
                        else:
                            h.option(value=sample.name())

def _linksToDownloads(sample):
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
                                    for descr, url in sample.downloads():
                                        index += 1
                                        if index % 5 == 0:
                                            h.br()
                                        elif index % 5 > 1:
                                            h.UNESCAPED(("&nbsp;&nbsp;|"
                                                         "&nbsp;&nbsp;"))
                                        h.a(h.i(class_="icon-download"),
                                            descr, target="_blank", href=url)

def _assemblyStats(sample):
    def _assemblySettingsTable():
        with h.table(class_="table table-condensed"):
            h.tr(h.th("Parameter"), h.th("Value"))
            for param, value in sample.assemblySettings().items():
                h.tr(h.td(param), h.td(value))

    def _assemblyMetricsTable():
        with h.table(class_="table table-condensed"):
            h.tr(h.th("Metric"), h.th("Large Contigs"), h.th("All Contigs"))
            for metric, value in sample.metrics().items():
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
                            h.p("Assembly summary statistics for ",
                                sample.name() + ".")
                            with h.div(style="padding-top:5px"):
                                _assemblySettingsTable()
                            with h.div(style="padding-top:10px"):
                                _assemblyMetricsTable()

# main function in this module
def generateReport(samples):
    css_path = "/pluginMedia/AssemblerPlus/css/"
    css = ["kendo.common.min.css", "kendo.default.min.css", "kendo.ir.css",
           "ir.css", "app.css", "bootstrap.css", "bootstrap-custom.css",
           "bootstrap-select.min.css"]
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
                            "  $('div.row-fluid.sample').hide() \n"
                            "  $('div#'+sample+'.row-fluid.sample').show() \n"
                            "}")
            with h.div(class_="main"):
                with h.div(class_="main-content clearfix"):
                    with h.div(class_="container-fluid"):
                        _alertAboutAssemblerWebsites()
                        if len(samples) > 1:
                            _sampleSelector(samples)
                        for i, sample in enumerate(samples):
                            display = "display:" + ("none" if i > 0 else "block")
                            with h.div(name=sample.name(), id=sample.name(),
                                       class_="row-fluid sample", style=display):
                                with h.div(class_="span12"):
                                    with h.div(style="padding:0 19px"):
                                        _linksToDownloads(sample)
                                        with h.div(style="padding-top:10px"):
                                            _assemblyStats(sample)

    return str(root)
