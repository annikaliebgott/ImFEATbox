import ImFEATbox, csv
import numpy as np
with open('testimg.csv', 'r') as csvfile:
	I = np.array(list(csv.reader(csvfile, delimiter=','))).astype(np.float)
typeflag = dict()
typeflag['corr'] = False
typeflag['texture'] = False
typeflag['global'] = True
typeflag['entropy'] = False
out = ImFEATbox.GlobalFeatures.Geometrical.glcm.cFeatures(I, typeflag)
with open('python-out.csv', 'wb') as csvfile:
	writer = csv.writer(csvfile, delimiter=',', escapechar=' ', quoting=csv.QUOTE_NONE)
	writer.writerow(out)