############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

from django.shortcuts import render_to_response, get_object_or_404, redirect
from django.http import HttpResponse
from tasks.models import Task, TaskLog
from django.template import Context, loader
from django.template import RequestContext
import datetime

def index(request):
    high_latest_tasks = Task.objects.filter(priority='H', status='A').order_by('pub_date')[:20]
    low_latest_tasks = Task.objects.filter(priority='L', status='A').order_by('pub_date')[:20]
    return render_to_response('tasks/index.html', {
                                'high_latest_tasks' : high_latest_tasks,
                                'low_latest_tasks' : low_latest_tasks},
                                 context_instance=RequestContext(request))

def datails(request, task_id):
    p = get_object_or_404(Task, pk=task_id)
    task_logs = TaskLog.objects.filter(reference_task=p).order_by('-changed_date')[:20]
    return render_to_response('tasks/details.html', {
                                                'task': p,
                                                'task_logs' : task_logs},
                              context_instance=RequestContext(request))



def newtask(request):
    return render_to_response('tasks/newtask.html', {},
                              context_instance=RequestContext(request))

def addtask(request):
    if request.POST.get('cancel'):
        return index(request)
    try:
        name_error = ''
        property1_error = ''
        property2_error = ''
        property3_error = ''
        name = request.POST['name']
        if not name:
            name_error = 'Name should be defined'
        priority = request.POST['priority']
        property1 = request.POST['property1']
        if not property1:
            property1_error = 'Property1 should be defined'
        property2 = request.POST['property2']
        if not property2:
            property2_error = 'Property2 should be defined'
        property3 = request.POST['property3']        
##        if not property3:
##            property3_error = 'Property3 should be defined'
        if name_error or property1_error or property2_error or property3_error:
            return render_to_response('tasks/newtask.html', {
                                            'name_error': name_error,
                                      'property1_error' : property1_error,
                                      'property2_error' : property2_error,
                                      'property3_error' : property3_error},
                                      context_instance=RequestContext(request))    
    except (KeyError):
        return render_to_response('tasks/newtask.html', {
                                  'name_error': name_error,
                                  'property1_error' : property1_error,
                                  'property2_error' : property2_error,
                                  'property3_error' : property3_error},
                                  context_instance=RequestContext(request))
    else:
        task = Task(name = name,
                    priority = priority,
                    property1 = property1,
                    property2 = property2,
                    property3 = property3,
                    pub_date = datetime.datetime.now(),
                    status = 'A'
                    )
        task.save()
    return redirect("/tasks/" + str(task.id) + "/added")

def added(request, task_id):
    task = Task.objects.get(id = task_id)
    return render_to_response('tasks/addedpage.html', {
                              'name' : task.name},
                              context_instance=RequestContext(request))    


def remove(request, task_id):   
    if request.POST.get("remove"):        
        task = Task.objects.get(id = task_id)
        if task:
            task.status = 'R'
            task.save()
            taskLog = TaskLog(
                name = task.name,
                reference_task = task,
                changed_date = datetime.datetime.now(),
                changed_property = 'Task was deleted'
            )
            taskLog.save();
    return redirect("/tasks")

def decrease(request, task_id):   
    if request.POST.get("decrease"):        
        task = Task.objects.get(id = task_id)
        if task:
            task.priority = 'L'
            task.save()
            taskLog = TaskLog(
                  name = task.name,
                  reference_task = task,
                  changed_date = datetime.datetime.now(),
                  changed_property = 'Task priority was decreased'
                  )
            taskLog.save();
    return redirect("/tasks")

def increase(request, task_id):   
    if request.POST.get("increase"):        
        task = Task.objects.get(id = task_id)
        if task:
            task.priority = 'H'
            task.save()
            taskLog = TaskLog(
                    name = task.name,
                    reference_task = task,
                    changed_date = datetime.datetime.now(),
                    changed_property = 'Task priority was increased'
                    )
            taskLog.save();

    return redirect("/tasks")

