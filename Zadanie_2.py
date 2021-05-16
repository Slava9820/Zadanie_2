import requests
import math
import scipy.special
import numpy as np
from scipy.special import spherical_jn, spherical_yn
import matplotlib.pyplot as plt
import requests as rqst
import os


# Расчет ЭПР(RCS - Radar cross section)
def RCS(lam, r):
    summ = 0
    kr = 2 * math.pi * r / lam
    # Задаем значения функций Бесселя для n = 0 для первой итерации
    J_prev = spherical_jn(0, kr)
    Y_prev = spherical_yn(0, kr)
    H_prev = J_prev + 1j * Y_prev
    for n in range(1, 50):
        # Вычисляем значения функций Бесселя для текущей n
        J_now = spherical_jn(n, kr)
        # J_prev = spherical_jn (n - 1, kr)
        Y_now = spherical_yn(n, kr)
        # Y_prev = spherical_yn (n - 1, kr)
        H_now = J_now + 1j * Y_now
        # H_prev = J_prev + 1j * Y_prev
        # Считаем коэффициенты a и b
        a = J_now / H_now
        b = (kr * J_prev - n * J_now) / (kr * H_prev - n * H_now)
        summ += ((-1) ** n) * (n + 0.5) * (b - a)
        # Переносим значения функций Бесселя на следующий шаг
        J_prev = J_now
        Y_prev = Y_now
        H_prev = H_now
    return lam * lam * np.abs(summ) * np.abs(summ) / math.pi


# Построение графика
def graf_freq(f, p):
    plt.plot(f, p)
    plt.ylabel('RCS, [$м^2$]')
    plt.xlabel('freq, [Гц]')
    plt.grid()
    plt.show()


# Загрузить файл с фариантом задания
print('Загрузка файла с вариантом задания')

def download(url):
    r = rqst.get(url)
    return r.text


# Разобрать прочитанные данные и найти нужные данные
print('Возврат данных по номеру варианта')


def var(text, nomervar):
    t = text.splitlines()
    return t[nomervar]


if __name__ == '__main__':
    txt = download('https://jenyay.net/uploads/Student/Modelling/task_02_02.txt')
    print()
    print('Исходный файл с заданием')
    print(txt)
    line = var(txt, 10)
    print('Требуемая строка варианта 10')
    print(line)
    L = list(line.split(' '))
    for i in L[::-1]:
        if i == '':
            L.remove(i)
    D = float(L[1])
    fmin = float(L[2])
    fmax = float(L[3])
    s = []
    ls = np.linspace(fmin, fmax, 300)
    Lambda = 3e8 / ls
    for f in ls:
        p = RCS(3e8 / f, D / 2)
        s.append(p)


        # Сохранить в файл требуемого формата
        def write_point(file, ls, Lambda, s):
            file.write("\t<row>\n")
            file.write("\t\t<freq>" + str(ls) + "</freq>\n")
            file.write("\t\t<lambda>" + str(Lambda) + "</lambda>\n")
            file.write("\t\t<rcs>" + str(s) + "</rcs>\n")
            file.write("\t</row>\n")


        path = os.path.dirname(__file__)
        path = os.path.join(path, "results.xml")
        with open(path, 'w') as f:
            f.write('<?xml version="1.1" encoding="UTF-8" ?>\n')
            f.write('<data>\n')
            for i in range(len(s)):
                write_point(f, ls[i], Lambda[i], s[i])
            f.write('</data>')
    graf_freq(ls, s)
print('Запись в xml-файл')



