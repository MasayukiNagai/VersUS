from django.views import generic
from .models import Mutation, Gene

class MutationView(generic.ListView):
    template_name= 'myapp/mutations.html'
    context_object_name = 'mutation_list'

    def get_queryset(self):
        return Mutation.objects.all()[:4]


class GeneView(generic.ListView):
    template_name = 'myapp/gene.html'
    context_object_name = 'gene_list'
    
    def get_queryset(self):
        return Gene.objects.all()[:4]


class GeneDetailView(generic.DetailView):
    model = Gene
    template_name = 'myapp/gene_detail.html'
