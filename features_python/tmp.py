import ImFEATbox, csv
import numpy as np
with open('testimg.csv', 'r') as csvfile:
	I = np.array(list(csv.reader(csvfile, delimiter=','))).astype(np.float)
gradtype = dict()
typeflag = dict()
gradtype['first'] = False
gradtype['second'] = False
typeflag['texture'] = False
typeflag['global'] = False
typeflag['gradient'] = False
out = ImFEATbox.GlobalFeatures.Intensity.gradient.cFeatures(I, typeflag, gradtype);
with open('python-out.csv', 'wb') as csvfile:
	writer = csv.writer(csvfile, delimiter=',')
	writer.writerow(out)