# Загрузка библиотек
import matplotlib
import math as mt
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import odeint, solve_ivp
from pathlib import Path
import json


def calc_ws(gamma_wat: float) -> float:
    """
    Функция для расчета солесодержания в воде

    :param gamma_wat: относительная плотность по пресной воде с плотностью 1000 кг/м3, безразм.

    :return: солесодержание в воде, г/г
    """
    ws = (
        1
        / (gamma_wat * 1000)
        * (
            1.36545 * gamma_wat * 1000
            - (3838.77 * gamma_wat * 1000 - 2.009 * (gamma_wat * 1000) ** 2) ** 0.5
        )
    )
    # если значение отрицательное, значит скорее всего плотность ниже допустимой 992 кг/м3
    if ws > 0:
        return ws
    else:
        return 0


def calc_rho_w(ws: float, t: float) -> float:
    """
    Функция для расчета плотности воды в зависимости от температуры и солесодержания

    :param ws: солесодержание воды, г/г
    :param t: температура, К

    :return: плотность воды, кг/м3
    """
    rho_w = 1000 * (1.0009 - 0.7114 * ws + 0.2605 * ws**2) ** (-1)

    return rho_w / (1 + (t - 273) * 1e-4 * (0.269 * (t - 273) ** 0.637 - 0.8))


def calc_mu_w(ws: float, t: float, p: float) -> float:
    """
    Функция для расчета динамической вязкости воды по корреляции Matthews & Russel

    :param ws: солесодержание воды, г/г
    :param t: температура, К
    :param p: давление, Па

    :return: динамическая вязкость воды, сПз
    """
    a = (
        109.574
        - (0.840564 * 1000 * ws)
        + (3.13314 * 1000 * ws**2)
        + (8.72213 * 1000 * ws**3)
    )
    b = 1.12166 - 2.63951 * ws + 6.79461 * ws**2 + 54.7119 * ws**3 - 155.586 * ws**4

    mu_w = (
        a
        * (1.8 * t - 460) ** (-b)
        * (0.9994 + 0.0058 * (p * 1e-6) + 0.6534 * 1e-4 * (p * 1e-6) ** 2)
    )
    return mu_w


def calc_n_re(rho_w: float, q_ms: float, mu_w: float, d_tub: float) -> float:
    """
    Функция для расчета числа Рейнольдса

    :param rho_w: плотность воды, кг/м3
    :param q_ms: дебит жидкости, м3/с
    :param mu_w: динамическая вязкость воды, сПз
    :param d_tub: диаметр НКТ, м

    :return: число Рейнольдса, безразмерн.
    """
    v = q_ms / (np.pi * d_tub**2 / 4)
    return rho_w * v * d_tub / mu_w * 1000


def calc_ff_churchill(n_re: float, roughness: float, d_tub: float) -> float:
    """
    Функция для расчета коэффициента трения по корреляции Churchill

    :param n_re: число Рейнольдса, безразмерн.
    :param roughness: шероховатость стен трубы, м
    :param d_tub: диаметр НКТ, м

    :return: коэффициент трения, безразмерн.
    """
    a = (-2.457 * np.log((7 / n_re) ** 0.9 + 0.27 * (roughness / d_tub))) ** 16
    b = (37530 / n_re) ** 16

    ff = 8 * ((8 / n_re) ** 12 + 1 / (a + b) ** 1.5) ** (1 / 12)
    return ff


def calc_dp_dl_grav(rho_w: float, angle: float):
    """
    Функция для расчета градиента на гравитацию

    :param rho_w: плотность воды, кг/м3
    :param angle: угол наклона скважины к горизонтали, градусы

    :return: градиент давления на гравитацию в трубе, Па/м
    """
    dp_dl_grav = rho_w * 9.81 * np.sin(angle / 180 * np.pi)
    return dp_dl_grav


def calc_dp_dl_fric(
    rho_w: float, mu_w: float, q_ms: float, d_tub: float, roughness: float
):
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
        dp_dl_fric = ff * rho_w * q_ms**2 / d_tub**5
    else:
        dp_dl_fric = 0
    return dp_dl_fric


def calc_dp_dl(
    rho_w: float, mu_w: float, angle: float, q_ms: float, d_tub: float, roughness: float
) -> float:
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


def __integr_func(
    h: float,
    pt: tuple,
    temp_grad: float,
    gamma_wat: float,
    angle: float,
    q_ms: float,
    d_tub: float,
    roughness: float,
) -> tuple:
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

    ws = calc_ws(gamma_wat=gamma_wat)

    rho_w = calc_rho_w(ws=ws, t=pt[1])

    mu_w = calc_mu_w(ws=ws, t=pt[1], p=pt[0])

    dp_dl = calc_dp_dl(
        rho_w=rho_w, mu_w=mu_w, angle=angle, q_ms=q_ms, d_tub=d_tub, roughness=roughness
    )

    dt_dl = temp_grad / 100

    return (dp_dl, dt_dl)


def calc_pipe(
    p_wh: float,
    t_wh: float,
    h0: float,
    md_vdp: float,
    temp_grad: float,
    gamma_wat: float,
    angle: float,
    q_ms: float,
    d_tub: float,
    roughness: float,
) -> tuple:
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

    # Решение ОДУ методом Рунге-Кутты
    solution = solve_ivp(
        fun=__integr_func,  # функция
        t_span=[h0, md_vdp],  # область интегрирования
        y0=[p_wh, t_wh],  # начальные условия для задачи Коши
        args=(temp_grad, gamma_wat, angle, q_ms, d_tub, roughness),  # аргументы функции
        method="RK45",  # метод Рунге-Кутты 4 порядка
    )

    p_wh_pipe, t_wh_pipe = solution.y[0], solution.y[1]
    md_vdp_pipe = solution.t

    return (p_wh_pipe, t_wh_pipe, md_vdp_pipe)


def calc_p_wf(
    p_wh: float,
    t_wh: float,
    h0: float,
    md_vdp: float,
    temp_grad: float,
    gamma_wat: float,
    angle: float,
    q_ms: float,
    d_tub: float,
    roughness: float,
) -> float:
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

    solution = calc_pipe(
        p_wh=p_wh,
        t_wh=t_wh,
        h0=h0,
        md_vdp=md_vdp,
        temp_grad=temp_grad,
        gamma_wat=gamma_wat,
        angle=angle,
        q_ms=q_ms,
        d_tub=d_tub,
        roughness=roughness,
    )

    p_wf = solution[0][-1]

    return p_wf


# Путь к папке с проектом
project_path = Path.cwd()

# Импорт начальных данных из папки input_data
with open(project_path.parent.joinpath("input_data", "2.json"), "r") as input_file:
    input_data = json.load(input_file)

# Диапазон дебитов жидкости для генерации VLP
q_liquid = list(range(0, 400, 10))

# Расчет забойных давлений для значений дебита жидкости
pressure_wf = [
    calc_p_wf(
        p_wh=input_data["p_wh"] * 101325,
        t_wh=input_data["t_wh"] + 273,
        h0=0,
        md_vdp=input_data["md_vdp"],
        temp_grad=input_data["temp_grad"],
        gamma_wat=input_data["gamma_water"],
        angle=input_data["angle"],
        q_ms=q / 86400,
        d_tub=input_data["d_tub"],
        roughness=input_data["roughness"],
    )
    / 101325
    for q in q_liquid
]

# Данные на выходе
output_data = {"q_liq": q_liquid, "p_wf": pressure_wf}

# Запись результатов в файл
with open(project_path.joinpath("output.json"), "w") as output_file:
    json.dump(output_data, output_file, indent='\t')