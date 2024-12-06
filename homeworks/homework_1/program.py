#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import json
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# In[2]:


def calc_ws( gamma_wat: float ) -> float:
    """
    Функция для расчета солесодержания в воде
    :param gamma_wat: относительная плотность по пресной воде с плотностью 1000 кг/м3, безразм.
    :return: солесодержание в воде, г/г
    """
    ws = (1 / (gamma_wat * 1000) * (1.36545 * gamma_wat * 1000 - (3838.77 * gamma_wat * 1000 - 2.009 * (gamma_wat * 1000) ** 2) ** 0.5))
    # если значение отрицательное, значит скорее всего плотность ниже допустимой 992 кг/м3
    if ws > 0:
        return ws
    else:
        return 0


# In[3]:


def calc_rho_w(ws: float,t: float) -> float:
    """
    Функция для расчета плотности воды в зависимости от температуры и солесодержания
    :param ws: солесодержание воды, г/г
    :param t: температура, К
    :return: плотность воды, кг/м3
    """
    rho_w = 1000 * (1.0009 - 0.7114 * ws + 0.2605 * ws ** 2) ** (-1)
    return rho_w / (1 + (t - 273) * 1e-4 * (0.269 * (t - 273) ** 0.637 - 0.8))


# In[4]:


def calc_mu_w(ws: float,t: float,p: float) -> float:
    """
    Функция для расчета динамической вязкости воды по корреляции Matthews & Russel
    :param ws: солесодержание воды, г/г
    :param t: температура, К
    :param p: давление, Па
    :return: динамическая вязкость воды, сПз
    """
    a = (109.574-(0.840564 * 1000 * ws)+(3.13314 * 1000 * ws ** 2) + (8.72213 * 1000 * ws ** 3))
    b = (1.12166-2.63951 * ws+ 6.79461 * ws ** 2 + 54.7119 * ws ** 3 - 155.586 * ws ** 4)
    mu_w = (a * (1.8 * t - 460) ** (-b)* (0.9994 + 0.0058 * (p * 1e-6) + 0.6534 * 1e-4 * (p * 1e-6) ** 2))
    return mu_w


# In[5]:


def calc_n_re(rho_w: float,q_ms: float,mu_w: float,d_tub: float) -> float:
    """
    Функция для расчета числа Рейнольдса
    :param rho_w: плотность воды, кг/м3
    :param q_ms: дебит жидкости, м3/с
    :param mu_w: динамическая вязкость воды, сПз
    :param d_tub: диаметр НКТ, м
    
    :return: число Рейнольдса, безразмерн.
    """
    v = q_ms / (np.pi * d_tub ** 2 / 4)
    return rho_w * v * d_tub / mu_w * 1000


# In[6]:


def calc_ff_churchill(n_re: float,roughness: float,d_tub: float) -> float:
    """
    Функция для расчета коэффициента трения по корреляции Churchill
    :param n_re: число Рейнольдса, безразмерн.
    :param roughness: шероховатость стен трубы, м
    :param d_tub: диаметр НКТ, м

    :return: коэффициент трения, безразмерн.
    """
    a = (-2.457 * np.log((7 / n_re) ** 0.9 + 0.27 * (roughness / d_tub))) ** 16
    b = (37530 / n_re) ** 16
    ff = 8 * ((8 / n_re) ** 12 + 1 / (a + b) ** 1.5) ** (1/12)
    return ff


# In[7]:


def calc_ff_jain(n_re: float,roughness: float,d_tub: float) -> float:
    """
    Функция для расчета коэффициента трения по корреляции Jain
    :param n_re: число Рейнольдса, безразмерн.
    :param roughness: шероховатость стен трубы, м
    :param d_tub: диаметр НКТ, м

    :return: коэффициент трения, безразмерн.
    """
    if n_re < 3000:
        ff = 64 / n_re
    else:
        ff = 1 / (1.14 - 2 * np.log10(roughness / d_tub + 21.25 / (n_re**0.9))) ** 2
    return ff


# In[8]:


def calc_dp_dl_grav(rho_w: float, angle: float):
    """
    Функция для расчета градиента на гравитацию

    :param rho_w: плотность воды, кг/м3
    :param angle: угол наклона скважины к горизонтали, градусы

    :return: градиент давления на гравитацию в трубе, Па/м
    """
    dp_dl_grav = rho_w * 9.81 * np.sin(angle / 180 * np.pi)
    return dp_dl_grav


# In[9]:


def calc_dp_dl_fric(rho_w: float,mu_w: float,q_ms: float, d_tub: float,roughness: float):
    """
    Функция для расчета градиента давления на трение
    :param rho_w: плотность воды, кг/м3
    :param mu_w: динамическая вязкость воды, сПз
    :param q_ms: дебит жидкости, м3/с
    :param d_tub: диаметр НКТ, м
    :param roughness: шероховатость стен трубы, м

    :return: градиент давления в трубе, Па/м
    """
    if q_ms != 0:
        n_re = calc_n_re(rho_w, q_ms, mu_w, d_tub)
        ff = calc_ff_churchill(n_re, roughness, d_tub)
        dp_dl_fric = ff * rho_w * q_ms ** 2 / d_tub ** 5
    else:
        dp_dl_fric = 0
    return dp_dl_fric


# In[10]:


def calc_dp_dl(rho_w: float,mu_w: float,angle: float,q_ms: float, d_tub: float,roughness: float) -> float:
    """
    Функция для расчета градиента давления в трубе

    :param rho_w: плотность воды, кг/м3
    :param mu_w: динамическая вязкость воды, сПз
    :param angle: угол наклона скважины к горизонтали, градусы
    :param q_ms: дебит жидкости, м3/с
    :param d_tub: диаметр НКТ, м
    :param roughness: шероховатость стен трубы, м

    :return: градиент давления в трубе, Па/м
    """
    dp_dl_grav = calc_dp_dl_grav(rho_w, angle)
    dp_dl_fric = calc_dp_dl_fric(rho_w, mu_w, q_ms, d_tub, roughness)
    dp_dl = dp_dl_grav - 0.815 * dp_dl_fric
    return dp_dl


# In[11]:


def __integr_func(h: float, pt: tuple,temp_grad: float,gamma_wat: float,angle: float,q_ms: float,d_tub: float,roughness: float) -> tuple:
    """
    Функция для интегрирования трубы

    :param h: текущая глубина, м
    :param pt: текущее давление, Па и текущая температура, К
    :param temp_grad: геотермический градиент, К/м * (1e-2)
    :param gamma_wat: относительная плотность по пресной воде с плотностью 1000 кг/м3, безразм.
    :param angle: угол наклона скважины к горизонтали, градусы
    :param q_ms: дебит жидкости, м3/с
    :param d_tub: диаметр НКТ, м
    :param roughness: шероховатость стен трубы, м

    :return: градиенты давления, Па/м и температуры, К/м
    """
    dp_dl_cc=calc_dp_dl(calc_rho_w(calc_ws(gamma_wat),pt[1]),calc_mu_w(calc_ws(gamma_wat), pt[1], pt[0]), angle, q_ms, d_tub, roughness)
    return dp_dl_cc, temp_grad


# In[12]:


def calc_pipe( p_wh: float,t_wh: float,h0: float,md_vdp: float,temp_grad: float,gamma_wat: float,angle: float,q_ms: float, d_tub: float,roughness: float) -> tuple:
    """
    Функция для расчета давления в трубе

    :param p_wh: буферное давление, Па
    :param t_wh: температура жидкости у буферной задвижки, К
    :param h0: начальная глубина, м
    :param md_vdp: глубина верхних дыр перфорации, м
    :param temp_grad: геотермический градиент, К/м * (1e-2)
    :param gamma_wat: относительная плотность по пресной воде с плотностью 1000 кг/м3, безразм.
    :param angle: угол наклона скважины к горизонтали, градусы
    :param q_ms: дебит жидкости, м3/с
    :param d_tub: диаметр НКТ, м
    :param roughness: шероховатость стен трубы, м

    :return: давление, Па и температура, K, глубины
    """
    res = solve_ivp(fun=__integr_func, t_span=[h0, md_vdp], y0=[p_wh, t_wh],args=(temp_grad, gamma_wat, angle, q_ms, d_tub, roughness))

    return (res.y[0], res.y[1])


# In[13]:


def calc_p_wf(p_wh: float, t_wh: float,h0: float,md_vdp: float,temp_grad: float,gamma_wat: float,angle: float,q_ms: float, d_tub: float,roughness: float) -> float:
    """
    Функция для расчета давления на забое скважины
    :param p_wh: буферное давление, Па
    :param t_wh: температура жидкости у буферной задвижки, К
    :param h0: начальная глубина, м
    :param md_vdp: глубина верхних дыр перфорации, м
    :param temp_grad: геотермический градиент, К/м * (1e-2)
    :param gamma_wat: относительная плотность по пресной воде с плотностью 1000 кг/м3, безразм.
    :param angle: угол наклона скважины к горизонтали, градусы
    :param q_ms: дебит жидкости, м3/с
    :param d_tub: диаметр НКТ, м
    :param roughness: шероховатость стен трубы, м

    :return: давление на забое скважины, Па
    """
    res,_ = calc_pipe(p_wh, t_wh, h0, md_vdp, temp_grad, gamma_wat, angle, q_ms, d_tub, roughness)

    return res[-1]


# In[18]:


with open('12.json') as json_file:
    data = json.load(json_file)
Q = np.linspace(0, 400, 100)

gamma_water = data['gamma_water']
H = data['md_vdp']
d_tub = data['d_tub']
angle = data['angle']
roughness = data['roughness']
p_wh = data['p_wh'] * 101325
t_wh = data['t_wh'] + 273
temp_grad = data['temp_grad'] / 100

press = list()

for q in Q:
    p = calc_p_wf(p_wh, t_wh, 0, H, temp_grad, gamma_water, angle, q / 86400, d_tub, roughness)
    press.append(p/101325) 
res = {'q_liq': list(Q), 'p_wf': list(press)}
with open(r'output.json','w') as file:
    json.dump(res,file,indent='\t')


# In[20]:


get_ipython().system('jupyter nbconvert --to script dzpotoki.ipynb')


# In[ ]:




