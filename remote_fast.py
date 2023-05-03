import ADC_T15
import time
import json

polys = [34, 35] #numbers of polys

#regime = 'debug'
path_res = 'data/DRS/'
regime = 'debug'

shotn = 1

try:
    with open(path_res + 'shotn_list.json', 'r') as past_file:
        list_of_data = json.load(past_file)
except FileNotFoundError:
    list_of_data = {}

polys_active = []
list_of_data[shotn] = {}
for poly in polys:
    polys_active.append(ADC_T15.DRS(port=8000+poly))
for poly in polys_active:
    print(poly.waitData())
print('all polys waits data')
for i, poly in enumerate(polys_active):
    file_name = poly.saveData()
    print(polys[i], file_name)
    list_of_data[shotn][polys[i]] = file_name

print('ALL RAW DATA WAS SAVED')

with open(path_res + 'shotn_list.json', 'w') as new_file:
    json.dump(list_of_data, new_file, indent=2)


