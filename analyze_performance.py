#identifies the most expensive parts of the miniQMC computation by analyzing the xml output info.

import xml.etree.ElementTree as ET
import os

xml_filenames = []
for root, dirs, files in os.walk("."):
    xml_filenames += [f for f in files if 'xml' in f[-3:]]

print(xml_filenames)

   

def get_total_time(root):
    for timer in root.iter('timer'):
        name = timer.find('name').text
        if name=='Total':
            return timer.find('time_incl').text

def get_max_times(root):
    tuples = []
    for timer in root.iter('timer'):
        name = timer.find('name').text
        time = float(timer.find('time_excl').text)
        tup = (name, time)
        tuples += [tup]
    
    tuples.sort(key=lambda x: x[1])
    
    return tuples[-3:]
            


for name in xml_filenames:
    print("Tiling: ", name[-9], " ", name[-7], " ",name[-5])
    contents = ET.parse(name)
    root = contents.getroot()
    timing = root[1]
    total_time = get_total_time(timing)
    print("total time:", total_time)

    top3 = get_max_times(root)
    print('top 3: ', top3)
    

     
