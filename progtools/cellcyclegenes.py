import pandas as pd
from eutils import Client

def countgenes(dfseg,outf,species="Mus+musculus"):
	ec = Client(api_key="f7078150f88498b5a1d0ba56db5dcffc2908") # initialize NCBI tool
	targetgenes = pd.read_csv("./assets/cellcycle-reg-genes.csv",engine='python',index_col=0)
	for ind,row in targetgenes.iterrows():
		esr = ec.esearch(db='gene',term=ind+'[PREF] AND Homo+sapiens[ORGN]')
		if(len(esr.ids) > 1):
			print("ALERT: multiple results found for",ind,len(esr.ids)," taking first")
		elif(len(esr.ids) == 0):
			print("NO RESUTS FOR",ind)
			continue
		#print(esr.ids[0])
		egs = ec.efetch(db='gene', id=esr.ids[0])