#pragma once

#include <complex>
#include <vector>

#include "error.h"

int rfd_decode(uint8_t val);

class auto_fclose
{
public:
    auto_fclose(FILE * f) : m_f(f) { ; }
    ~auto_fclose()
    {
        fclose(m_f);
    }
protected:
    FILE * m_f;
};

template <class T>
void RFD_Load(const TCHAR * file, int64_t offset, int length, std::vector<std::complex<T> > & sample)
{    
    FILE * f;
    if (_tfopen_s(&f, file, _T("rb")))
    {
        Error(_T("Не удалось открыть файл '%s' для чтения"), file);
    }

    auto_fclose _f(f);

    _fseeki64(f, 0, SEEK_END);
    int64_t f_size = _ftelli64(f);

    if (f_size < length)
    {
        Error(_T("Размер файла '%s' менее требуемой выборки (filesize = %lli)"), file, f_size);
    }
    if ((offset % sizeof(uint32_t)) != 0)
    {
        Error(_T("Смещение в RFD-файле должно быть кратно 4 байтам"));
    }
    if (f_size < offset)
    {
        Error(_T("Указанное смещение в файле '%s' превышает размер файла (filesize = %lli)"), file, f_size);
    }
    if (f_size - offset < length)
    {
        Error(_T("После смещения от начала файла на указанную величину остаток данных в файле '%s' менее требуемой выборки (filesize = %lli)"),
            file, f_size);
    }

    _fseeki64(f, offset, SEEK_SET);

    sample.clear();

    for (;;)
    {
        uint32_t tmp;
        fread(&tmp, sizeof(tmp), 1, f);

        for (int j = 0; j < 8; j++)
        {
            sample.push_back(std::complex<T>(rfd_decode(tmp >> (j * 2)), rfd_decode(tmp >> (16 + j * 2))));
            if (sample.size() == length)
            {
                return;
            }
        }
    }
}
