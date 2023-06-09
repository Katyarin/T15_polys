## Набор скриптов для считывания данных с приборов диагностики томсоновского рассеяния лазерного излучения для токамака Т15-МД

### ADC_T15.py
***
Имеет два класа - медленное АЦП (*slow ADC*) и быстрое АЦП (*DRS*)

Методы **slow ADC(ip_adrr, port)**

При инициализации нужно ввести IP адрес компьютера, на котором запущена программа client и порт этой программы (прописан в cfg-файле)
- getValue(shotn) - считывает страницу данных медленного АЦП и сохраняет ее в файле shotn.json.
- close_conn() - закрывает соединение с прибором.

Методы **DRS(ip_adrr, port)**

При инициализации нужно ввести IP адрес компьютера, на котором запущена программа client и порт этой программы (прописан в cfg-файле)

- waitData() - запускает ожидание n страниц в преднастроенном быстром АЦП. Кол-во страниц указывается в самой программе, так же, как настройки уровня тригера и т.д.
- saveData() - если все страницы с данными были записаны, сохраняет данные с быстрого АЦП в двоичный файл и с медленного АЦП в текстовый файл. Путь сохранения указывается в cfg-файле программы client. Возвращает название сохраненного файла медленного АЦП, который отличается от названия файла с быстрым АЦП лишь добавкой 'lsADC_' в начале. Чтобы файлы с разных полихроматоров не перезаписывали друг друга, можно создать отдельную директорию под данные каждого полихроматора и в cfg-файле прописать путь до нее.

### Примеры использования 
***
Скрипт **remote_slow** показывает пример использования модуля slow ADC.

Скрипт **remote_fast** показывает пример использования модуля DRS. Так же он сохраняет файл shotn_list.json, который ставит в соответсвие введенному номеру разряда название файлов данных из каждого полихроматора.

### Пример обработки данных
***

Скрипт **txt_to_json.py** переводит записанный текстовый или бинарный файл с данными быстрого АЦП в json-файл. При этом названия файлов читаются из shotn_list.json, а итоговый файл аккумулирует сырые данные со всех АЦП в один файл с названием shotn.json в папку DRS/plasma/results/shotn/ 

Скрипт **plasma_exp.py** обрабатывает сырые данные, записанные с помощью **remote_fast** со спектральной, абсолютной калибровкой и ожидаемыми сигналами, которые должны храниться в папке source.
