from flask import Flask
from flask import render_template
from flask import request
from flask import send_file

import re
from urllib import quote, unquote

from shellder import *

app = Flask(__name__)

#TODO: start session

def make_url(string):
    files = re.findall(r"(?:[\w\.]+\/)+(?:\w+\.\w+)", string)
    print files
    res = string
    for f in files:
        print f
        url = "<a href=\"get?file=" + quote(f) + "\">" + f + "</a>"
        res = res.replace(f, url)
    return res

def format_output(lines):
    return "<br/>".join(map(make_url, lines))

@app.route('/', methods=['GET'])
def index():
    pipe = Shellder(pipe_out = "/tmp/vis_out")
    return render_template("index.html", console=format_output(pipe.get_output()))

@app.route('/command', methods=['POST'])
def command():
    pipe = Shellder("/tmp/vis_in", "/tmp/vis_out")
    pipe.send(request.form['command'])
    return format_output(pipe.get_output())

@app.route('/get')
def get():
    env_path = "../../../"
    file_path = env_path + unquote(request.args.get("file", ""))
    print("Getting", file_path)
    return send_file(file_path)

if __name__ == '__main__':
    app.debug = True
    app.run()
