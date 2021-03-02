from django.urls import path

from . import views

app_name = 'mutations'
urlpatterns = [
    path('', views.MutationView.as_view(), name='mutations_table'),
    path('gene/', views.GeneView.as_view(), name='genes_table'),
    path('gene/<int:pk>/', views.GeneDetailView.as_view(), name='gene_detail')
]