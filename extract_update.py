#identifies the most expensive parts of the miniQMC computation by analyzing the xml output info.

import xml.etree.ElementTree as ET
import sys


contents = ET.parse(sys.argv[1])
root = contents.getroot()
timing = root[1]
for timer in timing.iter('timer'):
	name = timer.find('name').text
	if name == 'Update':
		for t in timer.iter('timer'):
			name2 = t.find('name').text
			if name2 == 'Determinant::update':
				time = t.find('time_excl').text
				print(time)


