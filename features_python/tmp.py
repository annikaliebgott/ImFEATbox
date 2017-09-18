import ImFEATbox, csv
import numpy as np
with open('testimg.csv', 'r') as csvfile:
	I = np.array(list(csv.reader(csvfile, delimiter=','))).astype(np.float)
gradtype = dict()
typeflag = dict()
gradtype['first'] = True
gradtype['second'] = True
typeflag['texture'] = True
typeflag['global'] = True
typeflag['gradient'] = True
out = ImFEATbox.GlobalFeatures.Intensity.gradient.cFeatures(I, typeflag, gradtype)
with open('python-out.csv', 'wb') as csvfile:
	writer = csv.writer(csvfile, delimiter=',', escapechar=' ', quoting=csv.QUOTE_NONE)
	writer.writerow(out)