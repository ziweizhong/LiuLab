import pickle
import sys
import pprint
import pandas as pd

df = pd.read_csv(sys.argv[1]+'.csv')

out = {}

for i in range(0,len(df)):
	out[df.at[i,'Name']] = df.at[i,'Levy_Fit']

with open('fitness-'+sys.argv[1]+'.pickle', 'wb') as handle:
		pickle.dump(out, handle)