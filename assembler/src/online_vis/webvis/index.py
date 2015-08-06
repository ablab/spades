#Flask imports
import flask
from flask import Flask, session, request
from flask.ext.session import Session
#Aux imports
import re
from urllib import quote, unquote
#System imports
from os import listdir
from os.path import isfile, join
#App imports
from shellder import *
from dotjson import dot_to_json

# General app settings
app = Flask(__name__)
SESSION_TYPE = "filesystem"
app.config.from_object(__name__)
Session(app)

FILENAME_REGEXP = r"(?:[\w\-\.]+\/)+(?:[\w\-]+\.\w+)"

def make_url(string):
    files = re.findall(FILENAME_REGEXP, string)
    res = string
    for f in files:
        method = "render" if f.endswith(".dot") else "get"
        url = "<a href=\"%s?file=%s\">%s</a>" % (method, quote(f), f)
        res = res.replace(f, url)
    print res
    return res

def format_output(lines):
    res = "<br/>".join(map(make_url, lines))
    print res
    return res

env_path = "../../../"
cache_path = "static/cache/"
shellder = None

@app.route("/", methods=['GET'])
def index():
    global shellder
    if shellder is None:
        shellder = Shellder("/tmp/vis_in", "/tmp/vis_out", env_path)
        session["log"] = shellder.get_output()
    return flask.render_template("index.html", console=format_output(session["log"]))

@app.route("/logout", methods=['GET'])
def logout():
    global shellder
    if shellder is not None:
        shellder.close()
        shellder = None
        return "You have been logged out"
    return "No opened session"

@app.route("/command", methods=['POST'])
def command():
    global shellder
    result = shellder.send(request.form["command"]).get_output()
    session["log"].extend(result)
    return format_output(result)

@app.route("/get")
def get():
    file_path = env_path + unquote(request.args.get("file", ""))
    return flask.send_file(file_path, as_attachment=True, attachment_filename=path.basename(file_path))

@app.route("/render")
def render():
    global cache_path
    file_path = unquote(request.args.get("file", ""))
    type = request.args.get("method", "svg")
    _, full_name = path.split(file_path)
    name_only, _ = path.splitext(full_name)
    if type == "svg":
        res_path = cache_path + name_only + ".svg"
        result = open(res_path, "w")
        subprocess.call(["dot", "-Tsvg", env_path + file_path], stdout=result)
        result.close()
        return res_path
    elif type == "json":
        input = open(env_path + file_path, "r")
        content = dot_to_json(input.read())
        if (request.args.get("result", "") == "response"):
           return content
        res_path = cache_path + name_only + ".json"
        result = open(res_path, "w")
        result.write(content)
        result.close()
        return res_path
    else:
        return "Unknown method"

@app.route("/static/cache/vertex/<vertex_id>")
def vertex(vertex_id):
    res_path = next((cache_path + f for f in listdir(cache_path) if f.endswith(vertex_id + "_.svg")), None)
    if res_path is None:
        #Render a new file
        shellder.send("draw_vertex " + vertex_id)
        file_path = re.finditer(FILENAME_REGEXP, shellder.get_output()[0]).next().group()
        _, full_name = path.split(file_path)
        name_only, _ = path.splitext(full_name)        
        res_path = cache_path + name_only + ".svg"
        result = open(res_path, "w")
        subprocess.call(["dot", "-Tsvg", env_path + file_path], stdout=result)
        result.close()
    return flask.redirect(res_path)
    
if __name__ == "__main__":
    app.debug = True
    app.secret_key = "somekey"
    app.run()
