import ImFEATbox, csv
import numpy as np
with open('testimg.csv', 'r') as csvfile:
	I = np.array(list(csv.reader(csvfile, delimiter=','))).astype(np.float)
plotflag = False
typeflag = dict()
typeflag['moments'] = True
typeflag['local'] = True
typeflag['corr'] = True
typeflag['texture'] = True
out = ImFEATbox.LocalFeatures.Line.lineprofile.cFeatures(I, typeflag, plotflag)
with open('python-out.csv', 'wb') as csvfile:
	writer = csv.writer(csvfile, delimiter=',', escapechar=' ', quoting=csv.QUOTE_NONE)
	writer.writerow(out)