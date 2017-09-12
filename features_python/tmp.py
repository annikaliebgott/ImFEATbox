import ImFEATbox, csv
import numpy as np
with open('testimg.csv', 'r') as csvfile:
	I = np.array(list(csv.reader(csvfile, delimiter=','))).astype(np.float)
typeflag = dict()
typeflag['corr'] = False
typeflag['texture'] = False
typeflag['global'] = False
typeflag['entropy'] = False
out = ImFEATbox.GlobalFeatures.Intensity.intensity.cFeatures(I, typeflag);
with open('python-out.csv', 'wb') as csvfile:
	writer = csv.writer(csvfile, delimiter=',')
	writer.writerow(out)