import socket
import struct
import json


class slow_ADC():

    def __init__(self, ip_adrr='localhost', port=8002):
        ADDRESS = ip_adrr  # IP компьютера, на котором запущен клиент
        PORT = port  # порт клиента, который соединен с прибором
        self.BUFFER_SIZE = 20 + 2048 * 4 + 2048 * 2 * 8
        print('port: ', PORT)
        ENCODING = 'utf-8'
        req = 'get'
        self.REQUEST_TEXT = req.encode(ENCODING)

        self.packetSize = 20 + 2048 * 4 + 2048 * 2 * 8
        self.size_int = 4
        self.size_f = 2
        self.ch_size = 8
        self.header = 20
        self.data_size = 2048

        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.settimeout(5)

        self.sock.connect((ADDRESS, PORT))
        print('connected slow ADC, poly %i' %(port-8000))

    def getValue(self, shotn):
        try:
            if not self.sock.send(self.REQUEST_TEXT) == len(self.REQUEST_TEXT):
                print('Failed to send!')
                return -1
            data = []
            read_data = {'header': {}, 'data': []}
            while (len(data) < self.BUFFER_SIZE):
                data.extend(self.sock.recv(self.BUFFER_SIZE))
            print(len(data))
            if not (data[0] == 92 and data[1] == 35):
                print('this data strange!')
                return 1

            read_data['header']['day'] = data[2]
            read_data['header']['h'] = data[3]
            read_data['header']['m'] = data[4]
            read_data['header']['s'] = data[5]

            read_data['header']['FreqH'] = data[6]
            read_data['header']['FreqL'] = data[7]
            read_data['header']['DataWidthH'] = data[8]
            read_data['header']['DataWidthL'] = data[9]
            day = data[2]
            hour = data[3]
            minut = data[4]
            sec = data[5]
            print('time:', day, 'day %i:%i:%i' % (hour, minut, sec))

            timestamp = []
            for num in range(self.data_size):
                local_data = []
                timestamp.append(struct.unpack('i', bytearray(
                    data[self.header + num * self.size_int: self.header + (num + 1) * self.size_int]))[0])
                local_data.append(struct.unpack('i', bytearray(
                    data[self.header + num * self.size_int: self.header + (num + 1) * self.size_int]))[0])
                for ch in range(self.ch_size):
                    local_data.append(struct.unpack('H', bytearray(data[
                                                                   self.header + self.size_int * self.data_size + num * self.size_f + ch * self.data_size * self.size_f: self.header + self.size_int * self.data_size + num * self.size_f + ch * self.data_size * self.size_f + self.size_f]))[
                                          0])
                read_data['data'].append(local_data)
            with open('%s.json' % shotn, 'w') as file:
                json.dump(read_data, file)
            return 0
        except socket.timeout:
            print('Connection timeout!')
            return -3

    def close_conn(self):
        self.sock.close()


class DRS():
    def __init__(self, ip_adrr='localhost', port=8002):
        ADDRESS = ip_adrr  # IP компьютера, на котором запущен клиент
        PORT = port  # порт клиента, который соединен с прибором
        self.BUFFER_SIZE = 60
        print('port: ', PORT)
        self.ENCODING = 'utf-8'
        req = 'get fast hex'
        self.REQUEST_TEXT = req.encode(self.ENCODING)

        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.settimeout(500)

        self.sock.connect((ADDRESS, PORT))
        print('connected DRS, poly %i' % (port - 8000))

    def waitData(self):
        try:
            if not self.sock.send(self.REQUEST_TEXT) == len(self.REQUEST_TEXT):
                print('Failed to send!')
                return -1
            return 'poly wait signals'
        except socket.timeout:
            print('Connection timeout!')
            return -3

    def saveData(self):
        try:
            return self.sock.recv(self.BUFFER_SIZE).decode(self.ENCODING)
        except socket.timeout:
            print('Connection timeout!')
            return -3
