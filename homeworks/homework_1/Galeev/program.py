"""Домашнее задание №1
Галеев Рауль.
Взял 15 вариант, т.к. в 3-ем варианте некорректные данные для d_tub"""

# !/usr/bin/env python
# coding: utf-8

# # Расчетная часть

# In[1]:


# Standard library
import json
import math as mt

# Third Party
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate


# функция расчета плотности воды в зависимости от давления и температуры
def rho_w_kgm3(P_Mpa, T_K, ws=0):
    """
    :param rho_w_kgm3: плотность воды, кг/м^3
    :param ws: солёность, 1/ (кг/м^3)^2
    :param T_K: температура, К
    """
    rho_w_sc_kgm3 = 1000 * (1.0009 - 0.7114 * ws + 0.2605 * ws**2) ** (-1)
    return rho_w_sc_kgm3 / (1 + (T_K - 273) / 10000 * (0.269 * (T_K - 273) ** (0.637) - 0.8))


# функция расчета солености через плотность
def salinity_gg(rho_kgm3):
    """
    :param rho_kgm3: плотность воды, кг/м^3
    """
    sal = 1 / rho_kgm3 * (1.36545 * rho_kgm3 - (3838.77 * rho_kgm3 - 2.009 * rho_kgm3**2) ** 0.5)
    # если значение отрицательное, значит скорее всего плотность ниже допустимой 992 кг/м3
    if sal > 0:
        return sal
    else:
        return 0


# Расчет вязкости воды в зависимости от температуры и давления
def visc_w_cP(P_Mpa, T_K, ws=0):
    """
    :param P_Mpa: давление, Mpa
    :param ws: солёность, 1/ (кг/м^3)^2
    :param T_K: температура, К
    :param muw: Вязкость, Пас
    """
    A = 109.574 - 0.8406 * 1000 * ws + 3.1331 * 1000 * ws * ws + 8.7221 * 1000 * ws * ws * ws
    B = (
        1.1217
        - 2.6396 * ws
        + 6.7946 * ws * ws
        + 54.7119 * ws * ws * ws
        - 155.586 * ws * ws * ws * ws
    )
    muw = (
        A
        * (1.8 * T_K - 460) ** (-B)
        * (0.9994 + 0.0058 * P_Mpa + 0.6534 * (10) ** (0 - 4) * P_Mpa * P_Mpa)
    )
    return muw


# Расчет числа Рейнольдса
def Re(q_m3day, d_m, mu_mPas=0.2, rho_kgm3=1000):
    """
    :param q_m3day : дебит жидкости, м3/сут
    :param rho_kgm3 : плотность воды или жидкости, по умолчанию 1000 кг/м3, чистая вода
    :param mu_mPas  : вязкость жидкости по умолчанию 0.2 мПас
    :param d_m      : диаметр трубы, м
    """
    v_ms = q_m3day / 86400 / 3.1415 * 4 / d_m**2
    return rho_kgm3 * v_ms * d_m / mu_mPas * 1000


def friction_Jain(q_m3day, d_m=0.089, mu_mPas=0.2, rho_kgm3=1000, roughness=0.000018):
    """
    param q_m3day: дебит жидкости, м3/сут
    param rho_kgm3: плотность воды или жидкости, по умолчанию 1000 кг/м3, чистая вода
    param mu_mPas: вязкость жидкости по умолчанию 0.2 мПас
    param d_m: диаметр трубы, м
    """
    Re_val = Re(q_m3day, d_m, mu_mPas, rho_kgm3)
    if Re_val < 3000:
        return 64 / Re_val
    else:
        return 1 / (1.14 - 2 * np.log10(roughness / d_m + 21.25 / (Re_val**0.9))) ** 2


def friction_Churchill(q_m3day, d_m=0.089, mu_mPas=0.2, rho_kgm3=1000, roughness=0.000018):
    """
    param q_m3day: дебит жидкости, м3/сут
    param rho_kgm3: плотность воды или жидкости, по умолчанию 1000 кг/м3, чистая вода
    param mu_mPas: вязкость жидкости по умолчанию 0.2 мПас
    param d_m: диаметр трубы, м
    """
    Re_val = Re(q_m3day, d_m, mu_mPas, rho_kgm3)
    A = (-2.457 * np.log((7 / Re_val) ** (0.9) + 0.27 * (roughness / d_m))) ** 16
    B = (37530 / Re_val) ** 16
    return 8 * ((8 / Re_val) ** 12 + 1 / (A + B) ** 1.5) ** (1 / 12)


# переделанная функция dp/dh
def dp_dh(h, p, t, q_liq, d_m, angle, ws, roughness):
    """
    param h: глубина, м
    param p: давление, МПа
    param t: температура, К
    param q_liq: дебит жидкости, м3/сут
    param d_m: диаметр трубы, м
    param angle: угол отклонения скважины, рад
    param ws: солёность
    param roughness: шероховатость
    """

    t = func_t(h)
    global rho_global

    ksi = 1 / 10**6
    # ws = salinity_gg(rho_global)
    # расчёт солёности на каждой итерации (Если взять данный расчёт, у графика в начале будет излом)
    ws = salinity_gg(
        gamma_water * 1000
    )  # принятие солёности постоянной на всех расчётах (выбрать одно из двух)
    rho = rho_w_kgm3(p, t, ws)
    mu = visc_w_cP(p, t, ws)
    rho_global = rho
    f = friction_Jain(q_liq, d_m, mu, rho, roughness)

    dp_dl_grav = rho * 9.81 * mt.sin(angle * mt.pi / 180)
    dp_dl_fric = f * rho * (q_liq / 86400) ** 2 / d_m**5

    return ksi * (dp_dl_grav - 0.815 * dp_dl_fric)


# Third Party
from scipy.integrate import solve_ivp

# Вариант №15
# json 15

with open("../../../../../../PycharmProjects/Course_MultiPhase_Flow/15.json", "r") as file:
    data = json.load(file)

gamma_water = data[
    "gamma_water"
]  # Относительная плотность воды по пресной воде с плотностью 1000 кг/м3 в с.у., безразм.
md_vdp = data["md_vdp"]  # Измеренная глубина верхних дыр перфорации, м
d_m = data["d_tub"]  # диаметр НКТ по которой ведется закачка, м
angle = data["angle"]  # угол отклонения скважины от вертикали, град
roughness = data["roughness"]  # шероховатость
p_wh = data["p_wh"]  # давление на устье, атм (нужно перевести в МПа)
t_wh = data["t_wh"]  # температура на устье, в цельсиях
temp_grad = data["temp_grad"]  # температурный градиент град С на 100 м

d_m = d_m * 0.1  # перевод из дециметров в метры
t_wh = t_wh + 273.15  # перевод в Кельвины
p_wh = p_wh * 0.101325  # перевод из атм в МПа

ws = salinity_gg(gamma_water * 1000)

rho_global = 1000

# массив для дебита
q_m3day = []
for i in range(0, 401, 20):
    q_m3day.append(i)
q_m3day[0] = 1
print("Дебит: ", q_m3day)

# Градиент температур
N = 20
Temp = []
for i in range(21):
    Temp.append(t_wh + temp_grad * i)
print("Температуры: ", Temp)

pres = []
pres.append(p_wh)

# Сетка
t_min = 0
t_max = md_vdp
t_e = np.linspace(t_min, t_max, 21)

# Интерполяциия функции распределения температуры
x = np.linspace(0, md_vdp)
y = t_wh + temp_grad * x / 100

func_t = interpolate.interp1d(x, y, fill_value="extrapolate")


pressures = []
pres_final = []

for j in range(0, len(q_m3day)):
    result = solve_ivp(
        dp_dh,
        t_span=[0, md_vdp],
        y0=np.array([p_wh]),
        t_eval=t_e,
        args=(t_wh, q_m3day[j], d_m, angle, ws, roughness),
    )
    pres_final.append(result.y[0, -1])


# Перевод в из МегаПа в атм (pres означает p result, а не reservoir)
for i in range(len(pres_final)):
    pres_final[i] = pres_final[i] * 9.8692327

# Экспорт данных в формат json
data_export = {
    "q_liq": q_m3day,
    "p_wf": pres_final
}

with open('output.json', 'w') as json_file:
    json.dump(data_export, json_file, indent=4)

plt.plot(q_m3day, pres_final)
plt.grid()
plt.title("Зависимость давления от дебита")
plt.ylabel("Давление, атм")
plt.xlabel("Дебит, м^3/сут")
plt.show()

