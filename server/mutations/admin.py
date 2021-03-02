from django.contrib import admin

# Register your models here.
from .models import Gene, Mutation

admin.site.register(Gene)
admin.site.register(Mutation)