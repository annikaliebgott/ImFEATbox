import ImFEATbox
import csv
import numpy as np
with open('testimg.csv', 'r') as csvfile:
	I = np.array(list(csv.reader(csvfile, delimiter=','))).astype(np.float)
out = ImFEATbox.LocalFeatures.Point.law.cFeatures(I)
with open('python-out.csv', 'wb') as csvfile:
	writer = csv.writer(csvfile, delimiter=',', escapechar=' ', quoting=csv.QUOTE_NONE)
	writer.writerow(out)
