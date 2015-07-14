#Flask imports
import flask
from flask import Flask, session, request
#System imports
from os import path
import os
import signal
import subprocess
#Aux imports
import re
from urllib import quote, unquote

from shellder import *

app = Flask(__name__)

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

#Making in-out pipes
def prepare_pipe(path):
    if not os.path.exists(path):
        os.mkfifo(path)

pipe_in = "/tmp/vis_in"
pipe_out = "/tmp/vis_out"
env_path = "../../../"

@app.route("/", methods=['GET'])
def index():
    log = ""
    if "pid" not in session:
        prepare_pipe(pipe_in)
        prepare_pipe(pipe_out)
        launcher = subprocess.Popen(["python", "launcher.py", pipe_in, pipe_out, env_path])
        print "Started launcher at", launcher.pid
        session["pid"] = launcher.pid
        log = Shellder(pipe_out = "/tmp/vis_out").get_output()
    return flask.render_template("index.html", console=log)

@app.route("/logout", methods=['GET'])
def logout():
    if "pid" in session:
        pid = session["pid"]
        session.pop("pid", None)
        print "Killing", pid
        os.kill(pid, signal.SIGTERM)
        return "You have been logged out"
    return "No opened session"

@app.route("/command", methods=['POST'])
def command():
    pipe = Shellder("/tmp/vis_in", "/tmp/vis_out")
    pipe.send(request.form["command"])
    return format_output(pipe.get_output())

@app.route("/get")
def get():
    file_path = env_path + unquote(request.args.get("file", ""))
    print("Getting", file_path)
    return flask.send_file(file_path, as_attachment=True, attachment_filename=path.basename(file_path))

@app.route("/render")
def render():
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

if __name__ == "__main__":
    app.debug = True
    app.secret_key = "somekey"
    app.run()
