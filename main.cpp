#pragma warning(disable: 4503)
#define NOMINMAX

#include <Windows.h>
#include <conio.h>
#include <tchar.h>

#include <thread>
#include <mutex>
#include <algorithm>

#include <omni\util.hpp>
#include <omni\dsp\DFT.h>
#include <omni\dsp\filter.h>
#include <gr\gr.h>

#include "Tools\sample.h"
#include "Tools\Galileo\E1_Codes.h"

/////////////////////////////////////////////////////////////////////
// режим вычислений

//#define INTEGER_MODE
#define DOUBLE_MODE

// режим вычислений
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// Параметры

// частота дискретизации входной выборки
static const double fs_MHz = 27.456;

// частота несущей входной выборки
static const double fi_MHz = 4.092;

// разрешение поиска по времени, отсчетов на чип BOC(1, 1)
static const int search_resolution = 2;

// неопределенность по частоте, кГц
// предполагается диапазон [-f, f], следует указать значение f
static const double frequency_unknown_area = 10;

// допустимый набег фазы за когерентное накопление (4 мс), радиан
static const double max_phase_change = omni::util::PI / 2;

// число некогерентных накоплений
static const int noncoherent_acc_num = 10;

// число колодцев гистограммы уровней
static const int num_of_level_bars = 100;

#ifdef INTEGER_MODE

static const double filter_multiplier = 128;
static const int filter_divider = 1;

static const double lo_multiplier = 10;
static const int lo_divider = 1;

#endif

// Параметры
/////////////////////////////////////////////////////////////////////

// коэффициенты фильтра
static const double base_filter_coeff[] =
{
    -1.36908894e-003,
    -2.54509044e-003,
    -4.49952238e-003,
    -6.56766423e-003,
    -6.75570700e-003,
    -2.18623603e-003,
    +9.84282553e-003,
    +3.05609151e-002,
    +5.87277228e-002,
    +9.04143641e-002,
    +1.19766230e-001,
    +1.40567031e-001,
    +1.48088440e-001,
    +1.40567031e-001,
    +1.19766230e-001,
    +9.04143641e-002,
    +5.87277228e-002,
    +3.05609151e-002,
    +9.84282553e-003,
    -2.18623603e-003,
    -6.75570700e-003,
    -6.56766423e-003,
    -4.49952238e-003,
    -2.54509044e-003,
    -1.36908894e-003,
};

// период кода E1 в чипах
// 1.023 Ms/s чиповая скорость первичного кода E1
// BOC(1, 1) удваивает чиповую скорость
// период первичного кода E1 4 ms
// тогда 1023 * 2 * 4 число чипов за период
static const int code_period = 1023 * 2 * 4;

// размер области неопределености по врмени
static const int time_unknown_area = code_period * search_resolution;

// чиповая скорость, Ms/s
// 1.023 Ms/s чиповая скорость первичного кода E1
// BOC(1, 1) удваивает чиповую скорость
static const double chip_rate = 2.046;

// нормализованную область неопределенности по времени определяем как блжайшую большую целую степень двойки от исходной области неопределенности
const int normalized_time_unknown_area = (int)round(pow(2, ceil(log2(time_unknown_area))));

// целевая частота дискретизации после передескритезации
static const double target_fs_MHz = chip_rate * search_resolution * normalized_time_unknown_area / (double)time_unknown_area;

// требуемая длительность выборки, мс
static const double target_duration_ms = 4 * noncoherent_acc_num;

// число гипотез по частоте
// определяем максимальный набег фазы за 4 мс и делим на допустимый набег фазы
static const int num_of_freq = ((int)ceil(2 * omni::util::PI * frequency_unknown_area /* кГц */ * 2 /* мс */ / max_phase_change)) * 2 + 1;

// число кодов Galileo
static const int num_of_codes = 32;

#ifdef INTEGER_MODE
#define DATA_T int
#endif

#ifdef DOUBLE_MODE
#define DATA_T double
#endif

struct Result_t
{
    int code;

    int f_Hz;
    int t_ns;

    DATA_T acc;

    int variant; // 0 - использовалась первая половина выборки
                 // 1 - использовалась вторая половина выборки
};

static std::mutex thread_mutex;

// функция потока поиска
// результат в формате results[code][group][frequency][time]
void search_thread_proc
(
    _In_ int begin_code,
    _In_ int end_code,
    _In_ const std::vector<std::complex<DATA_T> > & sample,
    _In_ const std::vector<std::complex<DATA_T> > ref[],
    _Out_ std::vector<struct Result_t> & results,
    _Out_ char progress[],
    _In_ int progres_bias
)
{

    // цикл по спутникам
    for (int code = begin_code; code < end_code; code++)
    {
        results[code].code = code;
        results[code].acc = 0;

        // цикл по частотным гипотезам
        for (int i = 0; i < num_of_freq; i++)
        {

            // копируем данные входной выборки для коррекции частоы
            std::vector<std::complex<DATA_T> > x(sample.begin(), sample.end());

            // текущая частотная гиротеза
            const double f_kHz = (i - (num_of_freq - 1) / 2) * frequency_unknown_area * 2 / (num_of_freq - 1);

            // коррекция частоты входной выборки в соответствии с частотной гипотезой
            for (int j = 0; j < (int)x.size(); j++)
            {
#ifdef INTEGER_MODE
                auto lo = std::polar(lo_multiplier, -omni::util::PI * 2 * f_kHz / (target_fs_MHz * 1000) * j);
                x[j] = x[j] * std::complex<DATA_T>((DATA_T)round(lo.real()), (DATA_T)round(lo.imag())) / lo_divider;
#endif
#ifdef DOUBLE_MODE
                x[j] *= std::polar(1.0, -omni::util::PI * 2 * f_kHz / (target_fs_MHz * 1000) * j);
#endif
            }

            // варианты "полувыборки" - первая или вторая половина выборки
            for (int j = 0; j < 2; j++)
            {

                // некогерентные накопления
                std::vector<DATA_T> acc(normalized_time_unknown_area);

                // цикл по некогерентным накоплениям
                for (int k = 0; k < noncoherent_acc_num; k++)
                {
                    std::vector<std::complex<DATA_T> > y(normalized_time_unknown_area);

                    std::copy(
                        &x[k * normalized_time_unknown_area + j * normalized_time_unknown_area / 2],
                        &x[k * normalized_time_unknown_area + (j + 1) * normalized_time_unknown_area / 2 - 1] + 1,
                        &y[j * normalized_time_unknown_area / 2]
                    );

                    omni::dsp::fft(y);
                    for (int l = 0; l < normalized_time_unknown_area; l++)
                    {
                        y[l] *= std::conj(ref[code][l]);
                    }
                    omni::dsp::ifft(y);

                    for (int l = 0; l < normalized_time_unknown_area; l++)
                    {
#ifdef INTEGER_MODE
                        int re = abs(y[l].real());
                        int im = abs(y[l].imag());
                        if (re > im)
                        {
                            acc[l] += re + im / 2;
                        }
                        else
                        {
                            acc[l] += im + re / 2;
                        }
#endif
#ifdef DOUBLE_MODE
                        acc[l] += std::abs(y[l]);
#endif
                    }
                }

                // обработка результатов накопления
                for (int k = 0; k < normalized_time_unknown_area; k++)
                {
                    if (results[code].acc < acc[k])
                    {
                        results[code].acc = acc[k];
                        results[code].f_Hz = (int)round(f_kHz * 1000);
                        results[code].t_ns = (int)round(k * 4e6 /* ns */ / normalized_time_unknown_area);
                        results[code].variant = j;
                    }
                }
            }
        }

        thread_mutex.lock();

        progress[progres_bias + code - begin_code] = '*';
        printf("\r%s", progress);

        thread_mutex.unlock();
    }
}

// функция потока вычисления свертки во временной области с заданным кодом и частотным сдвигом
void t_scan_proc
(
    _In_ int code,
    _In_ int f,
    _In_ int variant,
    _In_ const std::vector<std::complex<DATA_T> > & sample,
    _In_ const std::vector<std::complex<DATA_T> > ref[],
    _Out_ std::vector<DATA_T> & result
)
{

    // копируем данные входной выборки для коррекции частоы
    std::vector<std::complex<DATA_T> > x(sample.begin(), sample.end());

    // текущая частотная гипотеза
    const double f_kHz = f / 1000.0;

    // коррекция частоты входной выборки в соответствии с частотной гипотезой
    for (int i = 0; i < (int)x.size(); i++)
    {
#ifdef INTEGER_MODE
        auto lo = std::polar(lo_multiplier, -omni::util::PI * 2 * f_kHz / (target_fs_MHz * 1000) * i);
        x[i] = x[i] * std::complex<DATA_T>((DATA_T)round(lo.real()), (DATA_T)round(lo.imag())) / lo_divider;
#endif
#ifdef DOUBLE_MODE
        x[i] *= std::polar(1.0, -omni::util::PI * 2 * f_kHz / (target_fs_MHz * 1000) * i);
#endif;
    }

    // некогерентные накопления
    std::vector<DATA_T> acc(normalized_time_unknown_area);

    // цикл по некогерентным накоплениям
    for (int i = 0; i < noncoherent_acc_num; i++)
    {
        std::vector<std::complex<DATA_T> > y(normalized_time_unknown_area);

        std::copy(
            &x[i * normalized_time_unknown_area + variant * normalized_time_unknown_area / 2],
            &x[i * normalized_time_unknown_area + (variant + 1) * normalized_time_unknown_area / 2],
            &y[variant * normalized_time_unknown_area / 2]
        );

        omni::dsp::fft(y);
        for (int j = 0; j < normalized_time_unknown_area; j++)
        {
            y[j] *= std::conj(ref[code][j]);
        }
        omni::dsp::ifft(y);

        for (int j = 0; j < normalized_time_unknown_area; j++)
        {
#ifdef INTEGER_MODE
            int re = abs(y[j].real());
            int im = abs(y[j].imag());
            if (re > im)
            {
                acc[j] += re + im / 2;
            }
            else
            {
                acc[j] += im + re / 2;
            }
#endif
#ifdef DOUBLE_MODE
            acc[j] += std::norm(y[j]);
#endif
        }
    }

    std::copy(acc.begin(), acc.end(), result.begin());
}


void main()
{
    try
    {

        // коэффициенты фильтра
        DATA_T filter_coeff[_countof(base_filter_coeff)];
        for (int i = 0; i < (int)_countof(filter_coeff); i++)
        {
#ifdef INTEGER_MODE
            filter_coeff[i] = (DATA_T)round(base_filter_coeff[i] * filter_multiplier);
#endif
#ifdef DOUBLE_MODE
            filter_coeff[i] = base_filter_coeff[i];
#endif
        }

        // размер выборки для обработки
        int N = (int)ceil(target_duration_ms * (fs_MHz * 1e6) / 1000);

        printf("Loading...");

        // загрузка выборки на ПЧ
        std::vector<std::complex<DATA_T> > if_sample;
        RFD_Load(_T("D:\\PROJECTS\\ЭЛВИС\\ПО\\Samples\\28072017.rfd"), 12 * (1 << 20), N + _countof(filter_coeff), if_sample);

        printf("done\n");

        printf("Moving to baseband and filtering...");

        // фильтрация и перенос частоты на 0
        std::vector<std::complex<DATA_T> > baseband_sample(N);
        omni::dsp::FIR_Filter<std::complex<DATA_T>, DATA_T> filter(&filter_coeff[0], &filter_coeff[_countof(filter_coeff)]);
        for (int i = 0; i < (int)_countof(filter_coeff); i++)
        {
#ifdef INTEGER_MODE
            auto lo = std::polar(lo_multiplier, -2 * omni::util::PI * (fi_MHz / fs_MHz) * i);
            filter(if_sample[i] * std::complex<DATA_T>((DATA_T)round(lo.real()), (DATA_T)round(lo.imag())));
#endif
#ifdef DOUBLE_MODE
            filter(if_sample[i] * std::polar(1.0, -2 * omni::util::PI * (fi_MHz / fs_MHz) * i));
#endif
        }
        for (int i = 0; i < N; i++)
        {
            int n = i + _countof(filter_coeff);
#ifdef INTEGER_MODE
            auto lo = std::polar(lo_multiplier, -2 * omni::util::PI * (fi_MHz / fs_MHz) * n);
            baseband_sample[i] = filter(if_sample[n] * std::complex<DATA_T>((DATA_T)round(lo.real()), (DATA_T)round(lo.imag()))/ lo_divider) /
                filter_divider;
#endif
#ifdef DOUBLE_MODE
            baseband_sample[i] = filter(if_sample[n] * std::polar(1.0, -2 * omni::util::PI * (fi_MHz / fs_MHz) * n));
#endif
        }

        printf("done\n");

        printf("Resampling...");

        // передескритезация
        // реализуется методом выбора ближайшего отсчета в исходной выборке
        std::vector<std::complex<DATA_T> > sample;
        {
            double t_step = fs_MHz / target_fs_MHz;
            double t = 0;

            while (t < baseband_sample.size() - 1.0)
            {
                sample.push_back(baseband_sample[(int)round(t)]);
                t += t_step;
            }
        }

        printf("done\n");

        // DBG
        // уровни сигнала
        {
            std::map<DATA_T, int> raw_h;

            DATA_T min_value = std::numeric_limits<DATA_T>::max();
            DATA_T max_value = std::numeric_limits<DATA_T>::min();

            for (int i = 0; i < (int)sample.size(); i++)
            {
                DATA_T re = sample[i].real();
                DATA_T im = sample[i].imag();

                if (re > max_value)
                {
                    max_value = re;
                }
                if (re < min_value)
                {
                    min_value = re;
                }
                if (im > max_value)
                {
                    max_value = im;
                }
                if (im < min_value)
                {
                    min_value = im;
                }

                raw_h[re]++;
                raw_h[im]++;
            }

#ifdef INTEGER_MODE
            printf_s("Signal range [%d, %d]\n", min_value, max_value);
#endif
#ifdef DOUBLE_MODE
            printf_s("Signal range [%g, %g]\n", min_value, max_value);
#endif

            const double h_step = (raw_h.rbegin()->first - raw_h.begin()->first) / (double)num_of_level_bars;
            double h_limit = raw_h.begin()->first;
            double h_acc = 0;
            std::vector<double> y;
            std::vector<double> x;
            for (auto it = raw_h.begin(); it != raw_h.end(); ++it)
            {
                if (it->first >= h_limit)
                {
                    do
                    {
                        x.push_back(h_limit);
                        y.push_back(h_acc);
                        h_limit += h_step;
                        h_acc = 0;
                    } while (it->first > h_limit);
                }

                h_acc += it->second / (double)(sample.size() * 2);
            }
            x.push_back(h_limit + h_step);
            y.push_back(h_acc);

            gr::frame * frame = gr::create_frame(_T("Signal levels histogram"));
            std::pair<gr::frame*, gr::line*> ln = gr::create_line(frame, _T("Baseband signal"));
            gr::plot_xy(ln, y.begin(), y.end(), x.begin());
        }

        // DBG
        // спектр baseband
        {
            int probe_len = 16;
            while (probe_len < 32768)
            {
                if (probe_len * 2 >(int)sample.size())
                {
                    break;
                }
                probe_len *= 2;
            }
            std::vector<std::complex<double> > probe(probe_len);
            for (int i = 0; i < probe_len; i++)
            {
                probe[i] = std::complex<double>(sample[i].real(), sample[i].imag());
            }
            omni::dsp::fft(probe);
            omni::dsp::fft_shift(probe);
            gr::frame * frame = gr::create_frame(_T("Spectrum"));
            std::pair<gr::frame*, gr::line*> ln = gr::create_line(frame, _T("Baseband signal spectrum"));
            gr::plot_xy(ln, probe.begin(), probe.end(), gr_private::x_iterator(-0.5, 1 / (double)probe_len), gr::draw_norm);
        }

        // спектры опорных сигналов
        E1_SearchReference_Pilot<DATA_T> ref;
        double t_step = time_unknown_area / (double)normalized_time_unknown_area;
        std::vector<std::complex<DATA_T> > reference[num_of_codes];
        for (int i = 0; i < num_of_codes; i++)
        {
            
            // передискретизация опорного сигнала
            for (int j = 0; j < normalized_time_unknown_area; j++)
            {
                reference[i].push_back(ref(i, (int)round(j * t_step) / search_resolution));
                //reference[i].push_back(ref(i, (int)round(j * t_step) / search_resolution) * 4095);
            }

            // спектр
            omni::dsp::fft(reference[i]);
        }

        printf("Searching...\n");

        // результаты поиска
        std::vector<struct Result_t> results(num_of_codes);

        // нужно определить число процессоров
        SYSTEM_INFO si;
        GetSystemInfo(&si);

        // число процессоров, используемых для вычислений
        int num_of_processors_to_use = (int)si.dwNumberOfProcessors;

        // среднее число спутников, обрабатываемое одним потоком
        double average_sat_per_thread = num_of_codes / (double)num_of_processors_to_use;

        // индикатор прогресса вычислений
        char progress[num_of_codes + 2 /* [] */ + 1 /* EOL */] = { 0 };

        // очищаем и отображаем индикатор прогресса вычислений
        sprintf_s(progress, "[%*c]", num_of_codes, ' ');
        printf("\r%s", progress);

        // вычислительные потоки
        std::vector<std::thread> search_threads;

        // создание вычислительных потоков
        int code_begin = 0;
        for (int i = 0; i < num_of_processors_to_use; i++)
        {
            int last_code = (i == num_of_processors_to_use - 1) ?
                num_of_codes : (int)round((i + 1)*average_sat_per_thread);

            search_threads.push_back(std::thread(search_thread_proc, code_begin, last_code, std::cref(sample),
                reference, std::ref(results), progress, 1 + code_begin));

            code_begin = last_code;
        }

        // ожидание завершения вычислительных потоков
        for (int i = 0; i < (int)search_threads.size(); i++)
        {
            search_threads[i].join();
        }

        printf("\n");
        printf("done\n");

        printf("Post processing...");

        std::sort(results.begin(), results.end(), [](auto x, auto y) {return x.acc > y.acc; });

        printf("done\n");

        printf("+------+----------+----------+----------+\n");
        printf("| code |   f (Hz) |   t (ns) |      acc |\n");
        printf("+------+----------+----------+----------+\n");

        for (auto it = results.begin(); it != results.end(); ++it)
        {
#ifdef INTEGER_MODE
            printf("| %4d | %8d | %8d | %8d |\n", it->code, it->f_Hz, it->t_ns, it->acc);
#endif
#ifdef DOUBLE_MODE
            printf("| %4d | %8d | %8d | %8.4f |\n", it->code, it->f_Hz, it->t_ns, it->acc);
#endif
        }

        printf("+------+----------+----------+----------+\n");

        // DBG
        // свертка во временной области
        {
            printf("Graphics preparing...");

            int n = std::min(num_of_processors_to_use, 4);

            std::vector<std::vector<DATA_T> > z(n);
            for (int i = 0; i < n; i++)
            {
                z[i].resize(normalized_time_unknown_area);
            }

            // вычислительные потоки
            std::vector<std::thread> t_scan_threads;

            // создание вычислительных потоков
            for (int i = 0; i < n; i++)
            {
                t_scan_threads.push_back(std::thread(t_scan_proc, results[i].code, results[i].f_Hz, results[i].variant, std::cref(sample), reference,
                    std::ref(z[i])));
            }

            // ожидание завершения вычислительных потоков
            for (int i = 0; i < n; i++)
            {
                t_scan_threads[i].join();
            }

            gr::frame * frame = gr::create_frame(_T("Seaching results"));
            for (int i = 0; i < n; i++)
            {
                // SIR
                double i_acc = 0;
                int i_cnt = 0;
                double s = 0;
                for (int j = 0; j < normalized_time_unknown_area; j++)
                {
                    double t = (int)round(j * 4e6 /* ns */ / normalized_time_unknown_area);
                    if (fabs(t - results[i].t_ns) > 4000 /* ns */)
                    {
                        i_acc += z[i][j];
                        i_cnt++;
                    }
                    if (s < z[i][j])
                    {
                        s = z[i][j];
                    }
                }
                i_acc /= i_cnt;

                TCHAR t[80];
                _stprintf_s(t, _T("code = %d (SIR = %.1f dB)"), results[i].code, 10 * log10(s / i_acc));

                std::pair<gr::frame*, gr::line*> ln = gr::create_line(frame, t);
                gr::plot(ln, z[i].begin(), z[i].end());
            }

            printf("done\n");
        }

        _getch();

    }

    catch (Error & e)
    {
        e.message_box();
    }
}