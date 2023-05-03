import ADC_T15
import time

polys = [34, 35] #numbers of polys

#regime = 'debug'
path_res = 'data/slow/'
regime = 'debug'

time_wait = 60
try:
    with open(path_res + '%s/shotn.txt' %regime, 'r') as sh_file:
        for line in sh_file:
            shotn = int(line)
            break
except FileNotFoundError:
    with open(path_res + '%s/shotn.txt' % regime, 'w') as sh_file:
        sh_file.write('1')
        shotn = 1
print(shotn)

while True:
    for poly in polys:
        poly_active = ADC_T15.slow_ADC(port=8000+poly)
        print(poly_active.getValue(regime + '/' + str(poly)+'_' + str(shotn)))
        poly_active.close_conn()
    shotn += 1
    with open('%s/shotn.txt' % regime, 'w') as sh_file:
        sh_file.write(str(shotn))
    print('I am wait')
    time.sleep(time_wait)