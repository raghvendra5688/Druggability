# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.0
#   kernelspec:
#     display_name: ProtSol
#     language: python
#     name: protsol
# ---

from DeepPurpose import DTI as models
net = models.model_pretrained(model = 'MPNN_CNN_DAVIS')

net1 = models.model_pretrained(model = 'CNN_CNN_BindingDB_IC50')

net2 = models.model_pretrained(model = 'Morgan_CNN_BindingDB_IC50')

net3 = models.model_pretrained(model = 'Morgan_AAC_BindingDB_IC50')

net4 = models.model_pretrained(model = 'Daylight_AAC_BindingDB_IC50')

net5 = models.model_pretrained(model = 'CNN_CNN_DAVIS')

net6 = models.model_pretrained(model = 'CNN_CNN_BindingDB')

net7 = models.model_pretrained(model = 'Morgan_CNN_BindingDB')

net8 = models.model_pretrained(model = 'Morgan_CNN_DAVIS')

net9 = models.model_pretrained(model = 'MPNN_CNN_BindingDB')

net10 = models.model_pretrained(model = 'MPNN_CNN_KIBA')

net11 = models.model_pretrained(model = 'Transformer_CNN_BindingDB')

net12 = models.model_pretrained(model = 'Daylight_AAC_DAVIS')

net13 = models.model_pretrained(model = 'Daylight_AAC_KIBA')

net14 = models.model_pretrained(model = 'Daylight_AAC_BindingDB')

net15 = models.model_pretrained(model = 'Morgan_AAC_BindingDB')

net16 = models.model_pretrained(model = 'Morgan_AAC_KIBA')

net17 = models.model_pretrained(model = 'Morgan_AAC_DAVIS')

