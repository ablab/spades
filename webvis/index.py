from __future__ import print_function
#Flask imports
import flask
from flask import Flask, session, request
from flask.ext.session import Session
from jinja2 import Markup
#Aux imports
import re
from urllib import quote, unquote
#System imports
import os
from os import listdir
from os.path import isfile, join
#App imports
from shellder import *
from dotjson import dot_to_json
from ConfigParser import ConfigParser

# General app settings
app = Flask(__name__)
SESSION_TYPE = "filesystem"
app.config.from_object(__name__)
Session(app)

FILENAME_REGEXP = r"(?:\.{0,2}\/)[^\s]*"
SVG = ".svg"
JSON = ".json"
DOT = ".dot"

def _debug(*args):
    if app.debug:
        print(*args)

def make_url(string):
    files = re.findall(FILENAME_REGEXP, string)
    res = string
    for f in files:
        method = None
        if f.endswith(DOT):
            method = "render"
        elif f.endswith("/"):
            method = "ls"
        if method is not None:
            url = "<a href=\"/%s?path=%s\">%s</a>" % (method, quote(f), f)
            res = res.replace(f, url, 1)
    return res

def format_output(lines):
    return "".join(make_url(str(flask.escape(line))) + "<br/>" for line in lines)

env_path = ""
cache_path = "static/cache/"
shellders = dict()

@app.route("/", methods=['GET'])
def index():
    if "username" in session:
        return flask.render_template("index.html", username=session["username"])
    else:
        logged = ", ".join(shellders.keys())
        return flask.render_template("login.html", names=logged)

@app.route("/log", methods=['GET'])
def log():
    if "username" in session:
        return format_output(session["log"])
    else:
        return ""

@app.route("/login", methods=['GET'])
def login():
    session.permanent = True
    global shellders
    if "username" not in session:
        name = request.args.get("username")
        #TODO: server name validation
        if not name:
            name = "gaf"
        try:
            os.mkdir(path.join(cache_path, name))
        except OSError:
            pass
        if name not in shellders:
            try:
                launch = Shellder("/tmp/vis_in_" + name, "/tmp/vis_out_" + name, env_path)
                shellders[name] = launch
                session["log"] = launch.get_output()
            except Exception as e:
                message = "Cannot start online_vis session: " + str(e);
                return flask.render_template("logout.html", message=message)
        else:
            session["log"] = ["(the previous session log has been lost)"]
        session["username"] = name
    return flask.redirect("/")

@app.route("/logout", methods=['GET'])
def logout():
    message = "No opened session."
    if "username" in session:
        name = session.pop("username")
        shellder = shellders.pop(name, None)
        if shellder is not None:
            shellder.close()
        message = "You have been logged out."
    return flask.render_template("logout.html", message=message)

@app.route("/command", methods=['GET'])
def command():
    if session["username"] not in shellders:
        (result, complete) = (["online_vis is disconnected, please restart the session"], True)
    else:
        sh = shellders[session["username"]]
        com = request.args.get("command", "")
        if len(com):
            _debug("Sending `%s`..." % com)
            session["log"].append(">" + com)
            sh.send(com)
        (result, complete) = sh.get_output(5)
        session["log"].extend(result)
    return flask.jsonify(log=format_output(result), complete=complete)

@app.route("/get")
def get():
    full_path = augment(unquote(request.args.get("path", "")))
    try:
        return flask.send_file(full_path, as_attachment=True, attachment_filename=path.basename(full_path))
    except IOError:
        return flask.abort(404)

def make_cached(basename):
    return path.join(cache_path, session["username"], basename)

def do_render(full_path):
    file_name = path.basename(full_path)
    name_only, ext = path.splitext(file_name)
    #if not ext == DOT:
    #    raise RuntimeError()
    name = session["username"]
    res_path = make_cached(name_only + SVG)
    result = open(res_path, "w")
    subprocess.call(["dot", "-Tsvg", full_path], stdout=result)
    result.close()
    _debug("Written to", res_path)
    return res_path

@app.route("/render")
def render():
    full_path = augment(unquote(request.args.get("path", "")))
    _debug("Rendering", full_path)
    try:
        return do_render(full_path)
    except:
        flask.abort(500)

def do_graph(full_path):
    file_name = path.basename(full_path)
    name_only, ext = path.splitext(file_name)
    #if not ext == DOT:
    #    raise RuntimeError()
    input = open(full_path, "r")
    content = dot_to_json(input.read())
    input.close()
    return content

@app.route("/graph")
def graph():
    full_path = augment(unquote(request.args.get("path", "")))
    _debug("Building ", full_path)
    try:
        return do_graph(full_path)
    except:
        flask.abort(500)

@app.route("/vertex/<name>")
def vertex(name):
    if session["username"] not in shellders:
        return flask.abort(500)
    shellder = shellders[session["username"]]
    vertex_id, ext = path.splitext(name)
    res_path = next((augment(f) for f in listdir(cache_path) if f.endswith(name)), None)
    if res_path is None:
        #Render a new file
        shellder.send("draw_vertex " + vertex_id)
        out = " ".join(shellder.get_output())
        try:
            full_path = augment(re.finditer(FILENAME_REGEXP, out).next().group())
            if ext == SVG:
                return flask.redirect(do_render(full_path))
            elif ext == JSON:
                return do_graph(full_path)
            else:
                return flask.error(500)
        except:
            res_path = make_cached(vertex_id + "_err.txt")
            result = open(res_path, "w")
            result.write(out)
            result.close()
    return flask.redirect(res_path)

def augment(path):
    global env_path
    if path.startswith("/"):
        return path
    else:
        return env_path + path

@app.route("/ls")
def ls():
    path = unquote(request.args.get("path", ""))
    full_path = augment(path)
    _debug("Getting contents of", path)
    if full_path[-1] == "*":
        #For path completion
        try:
            files = subprocess.check_output("ls -pd " + full_path, shell=True).split("\n")[0:-1]
            _debug(files)
            #TODO: unify
            return " ".join(files)
        except:
            return ""
    try:
        #For accessing folder content
        content = [path + f for f in os.listdir(full_path)]
        files = [f for f in content if isfile(augment(f))]
        _debug(files)
        return format_output(files)
    except IOError as err:
        return err.strerror

if __name__ == "__main__":
    config = ConfigParser()
    config.readfp(open("webvis.cfg"))
    env_path = config.get("server", "env")
    port = config.getint("server", "port")
    app.debug = config.getboolean("server", "debug")
    app.secret_key = "somekey"
    app.run(port=port)
