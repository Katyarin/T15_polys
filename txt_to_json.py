from pathlib import Path
import json
import struct

path_res = 'data/DRS/'
res_data_path = 'data/DRS/plasma/results/'
raw_data_path = 'data/DRS/plasma/raw_data/'

ch_count = 7


def to_json(shotn, save_file=False) -> list:
    concatenated_data = {}
    with open(path_res + 'shotn_list.json', 'r') as past_file:
        list_of_data = json.load(past_file)
    for poly in list_of_data[str(shotn)]:
        if len(list_of_data[str(shotn)][poly])>6:
            shotnf = list_of_data[str(shotn)][poly][6:]
        else:
            shotnf = list_of_data[str(shotn)][poly]
        print(shotnf)
        # try binary file first
        path: Path = Path('%s%s' % (raw_data_path, shotnf))
        try:
            data = __bin_to_JSON(path_in=path)
        except:
            # try text file
            path: Path = Path('%s%s' % (raw_data_path, shotnf))
            if path.is_file():
                data = __ascii_to_JSON(path_in=path)
            else:
                print('file not found: ', path)
                stop
        concatenated_data[poly] = data
    #concatenated_data = data
    if save_file:
        print('Warning! Data is ready, but saving to disk is slow')
        #path = Path('%s/' % (res_data_path, shotn))
        path = Path('%s%d/' % (res_data_path, shotn))
        if not path.is_dir():
            path.mkdir(parents=True)

        with open('%s/%d.json' % (path, shotn), 'w') as file:
            json.dump(concatenated_data, file, indent=2)
    return concatenated_data


def __bin_to_JSON(path_in: Path) -> list:
    data: list = []
    with path_in.open(mode='rb') as file:
        raw = file.read()
        count = len(raw) / (8 * 1024 * 4 + 8)
        if not int(count) == count:
            stop
        count = int(count)

        for event_ind in range(count):
            event = {
                't': 0,
                'ch': []
            }
            for ch in range(8):
                event['ch'].append(struct.unpack_from('1024f', raw, (event_ind * (8 * 1024 * 4)) + ch * 1024 * 4))
            event['t'] = struct.unpack_from('L', raw, (count * (8 * 1024 * 4)) + 8 * event_ind)[0]
            data.append(event)
    return data


def __ascii_to_JSON(path_in: Path) -> list:
    data: list = []
    with path_in.open(mode='r') as file:
        count = 0
        event = {
            't': 0,
            'ch': [[] for ch in range(ch_count + 1)]
        }
        for line in file:
            if count > 1023:
                count += 1
                if count == 1026:
                    data.append(event.copy())
                    event = {
                        't': 0,
                        'ch': [[] for ch in range(ch_count + 1)]
                    }
                    count = 0
                continue
            sp = line.split()

            for ch in range(ch_count + 1):
                event['ch'][ch].append(float(sp[1 + ch]))
            if count == 0:
                event['t'] = int(sp[-2])
            count += 1
    return data


#to_json('2023-04-26_16-44-31')
#to_json(42198, save_file=True)
#to_json('1000.bin', save_file=True)