############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

from tasks.models import Task, TaskLog
from django.contrib import admin

class LogInline(admin.StackedInline):
    model = TaskLog

class TaskAdmin(admin.ModelAdmin):
    fieldsets = [
                 (None,                     {'fields': ['name', 
                                                        'priority',
                                                        'status',
                                                        ]}),
                 ('Properties information', {'fields': ['property1',
                                                        'property2',
                                                        'property3',
                                                        'pub_date',
                                                        ],
                                             'classes': ['collapse']})
    ]
    inlines = [LogInline]
    list_display=('name', 'priority')
    list_filter = ['pub_date', 'priority']
    search_fields = ['name']
    date_hierarchy = 'pub_date'

admin.site.register(Task, TaskAdmin)