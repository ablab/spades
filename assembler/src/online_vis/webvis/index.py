#Flask imports
import flask
from flask import Flask
from flask import request
#System imports
from os import path
import os
import subprocess
#Aux imports
import re
from urllib import quote, unquote

from shellder import *

app = Flask(__name__)

#TODO: start session

def make_url(string):
    files = re.findall(r"(?:[\w\.]+\/)+(?:\w+\.\w+)", string)
    res = string
    for f in files:
        method = "render" if f.endswith(".dot") else "get"
        url = "<a href=\"%s?file=%s\">%s</a>" % (method, quote(f), f)
        res = res.replace(f, url)
    return res

def format_output(lines):
    return "<br/>".join(map(make_url, lines))

@app.route('/', methods=['GET'])
def index():
    pipe = Shellder(pipe_out = "/tmp/vis_out")
    return flask.render_template("index.html", console=format_output(pipe.get_output()))

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
    return flask.send_file(file_path, as_attachment=True, attachment_filename=path.basename(file_path))

@app.route('/render')
def render():
    env_path = "../../../"
    pushd = os.getcwd()
    os.chdir(env_path)
    file_path = unquote(request.args.get("file", ""))
    dirfile, _ = path.splitext(file_path)
    res_path = dirfile + ".png"
    result = open(res_path, "w")
    subprocess.call(["dot", "-Tpng", file_path], stdout=result)
    result.close()
    os.chdir(pushd)
    return flask.redirect("/get?file=" + res_path)

if __name__ == '__main__':
    app.debug = True
    app.run()
