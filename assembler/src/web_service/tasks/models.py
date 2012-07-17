from django.db import models

class Task(models.Model):
    name = models.CharField(max_length=200)
    property1 = models.CharField(max_length=30)
    property2 = models.CharField(max_length=30)
    property3 = models.CharField(max_length=30)
    pub_date = models.DateTimeField('date published')
    PRIORITY_VALUES = (
                          ('H', 'High'),
                          ('L', 'Low'),
                     )
    priority = models.CharField(max_length=1, choices=PRIORITY_VALUES)
    STATUS_VALUES = (
                           ('A', 'Active'),
                           ('R', 'Removed'),
                           )
    status = models.CharField(max_length=1, choices=STATUS_VALUES)  

    def __unicode__(self):
        return self.name


    
class TaskLog(models.Model):
    name = models.CharField(max_length=200)
    reference_task = models.ForeignKey(Task)
    changed_date = models.DateTimeField('date changed')
    changed_property = models.CharField(max_length=100)

    
    def __unicode__(self):
        return self.name


    
