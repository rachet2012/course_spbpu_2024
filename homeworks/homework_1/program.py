from json import load, dump
from pathlib import Path, WindowsPath

from scipy.integrate import solve_ivp
from numpy import arange, power, pi, log, sin, array
from numpy.typing import NDArray


def calc_ws(gamma_wat: float) -> float:
    """
    Функция для расчета солесодержания в воде.

    :param gamma_wat: Относительная плотность по пресной воде с плотностью 1000 кг/м3, безразм.

    :return: Солесодержание в воде, г/г.
    """
    ws = (1 / (gamma_wat * 1000) * (1.36545 * gamma_wat * 1000 -
                                    power(3838.77 * gamma_wat * 1000 -
                                          2.009 * power(gamma_wat * 1000, 2), 0.5)))
    # если значение отрицательное, значит скорее всего плотность ниже допустимой 992 кг/м3
    if ws > 0:
        return ws
    else:
        return 0


def calc_rho_w(ws: float, t: float) -> float:
    """
    Функция для расчета плотности воды в зависимости от температуры и солесодержания.

    :param ws: Солесодержание воды, г/г.
    :param t: Температура, К.

    :return: Плотность воды, кг/м3.
    """
    rho_w = 1000 * 1 / (1.0009 - 0.7114 * ws + 0.2605 * power(ws, 2))

    return rho_w / (1 + (t - 273) * 1e-4 * (0.269 * power(t - 273, 0.637) - 0.8))


def calc_mu_w(ws: float, t: float, p: float) -> float:
    """
    Функция для расчета динамической вязкости воды по корреляции Matthews & Russel.

    :param ws: Солесодержание воды, г/г.
    :param t: Температура, К.
    :param p: Давление, Па.

    :return: Динамическая вязкость воды, сПз.
    """
    a = (109.574 - (0.840564 * 1000 * ws) + (3.13314 * 1000 * power(ws, 2)) +
         (8.72213 * 1000 * power(ws, 3)))

    b = (1.12166 - 2.63951 * ws + 6.79461 * power(ws, 2) + 54.7119 * power(ws, 3) -
         155.586 * power(ws, 4))

    return (a / power(1.8 * t - 460, b) * (0.9994 + 0.0058 * (p * 1e-6) +
                                           0.6534 * 1e-4 * power(p * 1e-6, 2)))


def calc_n_re(rho_w: float, q_ms: float, mu_w: float, d_tub: float) -> float:
    """
    Функция для расчета числа Рейнольдса.

    :param rho_w: Плотность воды, кг/м3.
    :param q_ms: Дебит жидкости, м3/с.
    :param mu_w: Динамическая вязкость воды, сПз.
    :param d_tub: Диаметр НКТ, м.

    :return: Число Рейнольдса, безразмерн.
    """
    v = q_ms / (pi * power(d_tub, 2) / 4)

    return rho_w * v * d_tub / mu_w * 1000


def calc_ff_churchill(n_re: float, roughness: float, d_tub: float) -> float:
    """
    Функция для расчета коэффициента трения по корреляции Churchill.

    :param n_re: Число Рейнольдса, безразмерн.
    :param roughness: Щероховатость стен трубы, м.
    :param d_tub: Диаметр НКТ, м.

    :return: Коэффициент трения, безразмерн.
    """
    a = power(-2.457 * log((7 / n_re) ** 0.9 + 0.27 * (roughness / d_tub)), 16)
    b = power(37530 / n_re, 16)

    return 8 * (power(8 / n_re, 12) + 1 / (a + b) ** 1.5) ** (1 / 12)


def calc_dp_dl_grav(rho_w: float, angle: float) -> float:
    """
    Функция для расчета градиента на гравитацию.

    :param rho_w: Плотность воды, кг/м3.
    :param angle: Угол наклона скважины к горизонтали, градусы.

    :return: Градиент давления на гравитацию в трубе, Па/м.
    """
    dp_dl_grav = rho_w * 9.81 * sin(angle / 180 * pi)

    return dp_dl_grav


def calc_dp_dl_fric(rho_w: float,
                    mu_w: float,
                    q_ms: float,
                    d_tub: float,
                    roughness: float) -> float:
    """
    Функция для расчета градиента давления на трение.

    :param rho_w: Плотность воды, кг/м3.
    :param mu_w: Динамическая вязкость воды, сПз.
    :param q_ms: Дебит жидкости, м3/с.
    :param d_tub: Диаметр НКТ, м.
    :param roughness: Шероховатость стен трубы, м.

    :return: Градиент давления в трубе, Па/м.
    """

    if q_ms != 0:

        n_re = calc_n_re(rho_w, q_ms, mu_w, d_tub)
        ff = calc_ff_churchill(n_re, roughness, d_tub)

        return ff * rho_w * power(q_ms, 2) / power(d_tub, 5)
    else:
        return 0


def calc_dp_dl(rho_w: float,
               mu_w: float,
               angle: float,
               q_ms: float,
               d_tub: float,
               roughness: float) -> float:
    """
    Функция для расчета градиента давления в трубе.

    :param rho_w: Плотность воды, кг/м3.
    :param mu_w: Динамическая вязкость воды, сПз.
    :param angle: Угол наклона скважины к горизонтали, градусы.
    :param q_ms: Дебит жидкости, м3/с.
    :param d_tub: Диаметр НКТ, м.
    :param roughness: Шероховатость стен трубы, м.

    :return: Градиент давления в трубе, Па/м.
    """

    return calc_dp_dl_grav(rho_w, angle) - 0.815 * calc_dp_dl_fric(rho_w, mu_w, q_ms, d_tub, roughness)


def __integr_func(h: float,
                  pt: tuple,
                  temp_grad: float,
                  gamma_wat: float,
                  angle: float,
                  q_ms: float,
                  d_tub: float,
                  roughness: float) -> tuple:
    """
    Функция для интегрирования трубы.

    ::param h: текущая глубина, м.
    :param pt: Текущее давление, Па и текущая температура, К.
    :param temp_grad: Геотермический градиент, К/м * (1e-2).
    :param gamma_wat: Относительная плотность по пресной воде с плотностью 1000 кг/м3, безразм.
    :param angle: Угол наклона скважины к горизонтали, градусы.
    :param q_ms: Дебит жидкости, м3/с.
    :param d_tub: Диаметр НКТ, м.
    :param roughness: Шероховатость стен трубы, м.

    :return: Градиенты давления, Па/м и температуры, К/м.
    """
    ws = calc_ws(gamma_wat=gamma_wat)

    rho_w = calc_rho_w(ws=ws, t=pt[1])

    mu_w = calc_mu_w(ws=ws, t=pt[1], p=pt[0])

    dp_dl = calc_dp_dl(
        rho_w=rho_w,
        mu_w=mu_w,
        angle=angle,
        q_ms=q_ms,
        d_tub=d_tub,
        roughness=roughness
    )

    return dp_dl, temp_grad


def calc_pipe(p_wh: float,
              t_wh: float,
              h0: float,
              md_vdp: float,
              temp_grad: float,
              gamma_wat: float,
              angle: float,
              q_ms: float,
              d_tub: float,
              roughness: float) -> tuple[NDArray, NDArray, NDArray]:
    """
    Функция для расчета давления в трубе.

    :param p_wh: Буферное давление, Па.
    :param t_wh: Температура жидкости у буферной задвижки, К.
    :param h0: Начальная глубина, м.
    :param md_vdp: Глубина верхних дыр перфорации, м.
    :param temp_grad: Геотермический градиент, К/м * (1e-2).
    :param gamma_wat: Относительная плотность по пресной воде с плотностью 1000 кг/м3, безразм.
    :param angle: Угол наклона скважины к горизонтали, градусы.
    :param q_ms: Дебит жидкости, м3/с.
    :param d_tub: Диаметр НКТ, м.
    :param roughness: Шероховатость стен трубы, м.

    :return: Давление, Па и температура, K, глубины.
    """

    # для решения используется конечно-разностный метод Адамса
    result = solve_ivp(
        fun=__integr_func,
        t_span=[h0, md_vdp],
        y0=[p_wh, t_wh],
        args=(temp_grad, gamma_wat, angle, q_ms, d_tub, roughness),
        method='LSODA')

    p_wh, t_wh = result.y[0], result.y[1]
    md_vdp = result.t

    return p_wh, t_wh, md_vdp


def calc_p_wf(p_wh: float,
              t_wh: float,
              md_vdp: float,
              temp_grad: float,
              gamma_wat: float,
              angle: float,
              q_ms: float,
              d_tub: float,
              roughness: float,
              h0: float = 0) -> float:
    """
    Функция для расчета давления на забое скважины.

    :param p_wh: Буферное давление, Па.
    :param t_wh: Температура жидкости у буферной задвижки, К.
    :param h0: Начальная глубина, м.
    :param md_vdp: Глубина верхних дыр перфорации, м.
    :param temp_grad: Геотермический градиент, К/м * (1e-2).
    :param gamma_wat: Относительная плотность по пресной воде с плотностью 1000 кг/м3, безразм.
    :param angle: Угол наклона скважины к горизонтали, градусы.
    :param q_ms: Дебит жидкости, м3/с.
    :param d_tub: Диаметр НКТ, м.
    :param roughness: Шероховатость стен трубы, м.

    :return: Давление на забое скважины, Па.
    """
    p_wf = calc_pipe(p_wh=p_wh,
                     t_wh=t_wh,
                     h0=h0,
                     md_vdp=md_vdp,
                     temp_grad=temp_grad,
                     gamma_wat=gamma_wat,
                     angle=angle,
                     q_ms=q_ms,
                     d_tub=d_tub,
                     roughness=roughness)[0][-1]

    return p_wf


def read_data(variant: int, project_path: WindowsPath) -> dict:
    """
    Функция считывает данные с входного json-файла из папки data по номеру варианта.

    :param project_path: Папка с проектом.
    :param variant: Номер варианта.

    :return: Словарь со входными данными
    """

    path_to_dir = project_path / "input_data" / f"{str(variant)}.json"

    try:
        with open(path_to_dir, 'r') as file:
            return load(file)
    except FileNotFoundError:
        raise FileNotFoundError(f"Файл не найден по указанному пути: {path_to_dir}")


def get_init_data(variant: int,
                  project_path: WindowsPath,
                  q_min: float = 0,
                  q_max: float = 400) -> tuple[NDArray, float, float, float, float, float, float, float, float]:

    init_data = read_data(variant, project_path)

    q_fluid = change_m3_day_to_m3_sec(arange(q_min, q_max, 10))  # Дебит флюида, м3/сек
    t_wh = change_celsius_to_kelvin(init_data['t_wh'])  # Температура на устье скважины, К
    p_wh = change_atm_to_pa(init_data['p_wh'])  # Давление на устье, Па

    md_vdp = init_data['md_vdp']  # Глубина забоя, м
    gamma_water = init_data['gamma_water']  # Относительная плотность по пресной воде с плотностью 1000 кг/м3, безразм.
    d_tub = init_data['d_tub']  # Диаметр НКТ, м
    angle = init_data['angle']  # Угол наклона, градусы
    roughness = init_data['roughness']  # Шероховатость трубы, м
    temp_grad = change_sm_to_m(init_data['temp_grad'])  # Геотермический градиент, К/м.

    return q_fluid, t_wh, p_wh, md_vdp, gamma_water, d_tub, angle, roughness, temp_grad


def change_m3_day_to_m3_sec(q_fluid_m3_day: NDArray) -> NDArray:
    """
    Функция изменяет единицы измерения дебитов с м3/д на м3/с.

    :param q_fluid_m3_day: Массив с дебитами жидкости в м3/д.

    :return: Массив с дебитами жидкости в м3/с.
    """
    return q_fluid_m3_day / 86400


def change_m3_sec_to_m3_day(q_fluid_m3_day: NDArray) -> NDArray:
    """
    Функция изменяет единицы измерения дебитов с м3/с на м3/д.

    :param q_fluid_m3_day: Массив с дебитами жидкости в м3/с.

    :return: Массив с дебитами жидкости в м3/д.
    """
    return q_fluid_m3_day * 86400


def change_celsius_to_kelvin(t_celsius: float) -> float:
    """
    Функция изменяет единицы измерения дебитов с градусов Цельсия на Кельвины.

    :param t_celsius: Температура в градусах Цельсия.

    :return: Температура в Кельвинах.
    """
    return t_celsius + 273


def change_atm_to_pa(p_atm: float) -> float:
    """
    Функция изменяет единицы измерения давления с атмосфер на Паскали.

    :param p_atm: Давление в атмосферах.

    :return: Давление в Паскалях.
    """
    return p_atm * 101325


def change_pa_to_atm(p_pa: float) -> float:
    """
    Функция изменяет единицы измерения давления с Паскалей на атмосферы.

    :param p_pa: Давление в Паскалях.

    :return: Давление в атмосферах.
    """
    return p_pa / 101325


def change_sm_to_m(temp_grad_sm: float) -> float:
    """
    Функция изменяет единицы измерения градиента температуры с К / см на К / м.

    :param temp_grad_sm: Температурный градиент в К / см.

    :return: Давление в атмосферах.
    """
    return temp_grad_sm / 100


def save_results(q_liquid: NDArray, p_wf: NDArray, project_path: WindowsPath) -> None:
    data = {"q_liq": q_liquid.tolist(), "p_wf": p_wf.tolist()}

    with open(project_path.joinpath("output.json"), "w") as result_file:
        dump(data, result_file, indent='\t')


if __name__ == '__main__':

    # вариант расчётного задания
    variant_number = 13

    path = Path.cwd()

    (q_m3_sec,
     t_wh_kelvin,
     p_wh_pa,
     md_vdp_m,
     gamma_water_nd,
     d_tub_m,
     angle_grad,
     roughness_m,
     temp_grad_kelvin_m) = get_init_data(variant_number, path)

    p_wf_pa = array([
        calc_p_wf(p_wh=p_wh_pa,
                  t_wh=t_wh_kelvin,
                  md_vdp=md_vdp_m,
                  temp_grad=temp_grad_kelvin_m,
                  gamma_wat=gamma_water_nd,
                  angle=angle_grad,
                  q_ms=q,
                  d_tub=d_tub_m,
                  roughness=roughness_m
                  ) for q in q_m3_sec
    ])

    save_results(change_m3_sec_to_m3_day(q_m3_sec), change_pa_to_atm(p_wf_pa), path)
