import pandas as pd
from eutils import Client
'''import gffutils

db = gffutils.create_db(
    './gencode.v19.annotation.gtf.gz',
    dbfn='gencode_v19.db',
    verbose=True,
    merge_strategy='error',
    disable_infer_transcripts=True,
    disable_infer_genes=True,
)
'''

def countgenes(dfseg,outf,species="Mus+musculus"):
	ec = Client(api_key="f7078150f88498b5a1d0ba56db5dcffc2908") # initialize NCBI tool
	targetgenes = pd.read_csv("./assets/cellcycle-reg-genes.csv",engine='python',index_col=0)
	for ind,row in targetgenes.iterrows():
		esr = ec.esearch(db='gene',term=ind+'[PREF] AND '+species+'[ORGN]')
		if(len(esr.ids) > 1):
			print("ALERT: multiple results found for",ind,"("+str(len(esr.ids))+" total)"," taking first")
		elif(len(esr.ids) == 0):
			print("NO RESUTS FOR",ind)
			continue
		#print(ind,esr.ids[0])
		egs = ec.efetch(db='gene', id=esr.ids[0]).entrezgenes[0]
		#refs = sorted([(r.acv, r.label) for r in egs.references])
		print(egs.references[0])


# '__class__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__',
# '__gt__', '__hash__', '__init__', '__le__', '__lt__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__',
# '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__',
# '_root_tag', '_xml_root', 'accession', 'acv', 'genomic_coords', 'heading', 'label', 'products', 'type', 'version'