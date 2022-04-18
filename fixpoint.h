#pragma once

#include <array>
#include <exception>
#include <algorithm>
#define MP_INT HIDEN_MP_TYPE
#include <gmp.h>
#undef MP_INT
#include <iostream>
#include <string>
#include <cmath>
#include <cstring>
//#include <charconv>

//	presision - total significant digits
//	decimals - digits after dot
template<size_t presision, size_t decimals, bool IsSMem = false>
class FixPoint
{
	// Pregenered numbers for multiplication and round.
	class Divider
	{
	public:
		__mpz_struct mpz_;
		static constexpr size_t cLimbs = (size_t(presision * 3.33) + 1 + 8 * sizeof(mp_limb_t) - 1) / 8 / sizeof (mp_limb_t);
		std::array<mp_limb_t, cLimbs> data_;

		void Init(size_t decimal)
		{
			std::array<mp_limb_t, cLimbs + 4> data;
			mpz_._mp_alloc = data.size();
			mpz_._mp_size = 0;
			mpz_._mp_d = &data[0];
			mpz_ui_pow_ui(&mpz_, 10, decimal);
			std::copy(data.begin(), data.begin() + data_.size(), data_.begin());
			mpz_._mp_alloc = data_.size();
			mpz_._mp_d = &data_[0];
		}

	};

	class Dividers : public std::array<Divider, presision + 1>
	{
	public:
		Dividers()
		{
			for (size_t i = 0; i < this->size(); i++)
				(*this)[i].Init(i);
		}
	};

	class AddRounder
	{
	public:
		__mpz_struct mpz_;
		static constexpr size_t cLimbs = (size_t(presision * 3.33) + 1 + 8 * sizeof(mp_limb_t) - 1) / 8 / sizeof (mp_limb_t);
		std::array<mp_limb_t, cLimbs> data_;

		void Init(size_t decimal)
		{
			std::array<mp_limb_t, cLimbs + 4> data;
			mpz_._mp_alloc = data.size();
			mpz_._mp_size = 0;
			mpz_._mp_d = &data[0];
			mpz_ui_pow_ui(&mpz_, 10, decimal);
			mpz_tdiv_q_ui(&mpz_, &mpz_, 2);
			std::copy(data.begin(), data.begin() + data_.size(), data_.begin());
			mpz_._mp_alloc = data_.size();
			mpz_._mp_d = &data_[0];
		}

	};

	class AddRounders : public std::array<AddRounder, presision + 1>
	{
	public:
		AddRounders()
		{
			for (size_t i = 0; i < this->size(); i++)
				(*this)[i].Init(i);
		}
	};

	static constexpr size_t cLimbs = (size_t(presision * 3.33) + 8 * sizeof(mp_limb_t) - 1) / 8 / sizeof (mp_limb_t);
	std::array<mp_limb_t, cLimbs + 1> data_;
	static Dividers dividers_;
	static AddRounders addRounders_;
    mutable __mpz_struct mpz_;

public:
    static FixPoint<presision, decimals, IsSMem> Get10Pow(long pow)
    {
        FixPoint<presision, decimals> ret;
        pow += long(decimals);
        if (pow >= long(presision))
                pow = long(presision) - 1;
        if (pow < 0 )
                pow = 0;
        std::copy(dividers_[pow].data_.begin(), dividers_[pow].data_.begin() + ret.data_.size(), ret.data_.begin());
        ret.mpz_._mp_size = dividers_[pow].mpz_._mp_size;
        return  ret;
    }

    static FixPoint<presision, decimals, IsSMem> zero()
    {
        static FixPoint<presision, decimals, IsSMem> ret("0");
        return ret;
    }

	FixPoint()
	{
		mpz_._mp_alloc = data_.size();
		mpz_._mp_size = 0;
		mpz_._mp_d = &data_[0];
		//mpz_init(&mpz_);
	}

	FixPoint(const char* strVal)
	  :FixPoint()
	{
		*this = strVal;
	}

	FixPoint(const std::string_view& strVal)
	  :FixPoint()
	{
		*this = strVal;
	}

	FixPoint(double val)
	  :FixPoint()
	{
		*this = val;
	}

    FixPoint(const FixPoint<presision, decimals, IsSMem>& val)
	  :FixPoint()
	{
		*this = val;
	}
/*
	FixPoint(long val)
		:FixPoint()
	{
		*this = val;
	}

	FixPoint(unsigned long val)
		:FixPoint()
	{
		*this = val;
	}
*/
	~FixPoint()
	{
	}

    FixPoint<presision, decimals, IsSMem>& operator= (const FixPoint<presision, decimals, IsSMem>& from)
	{
		mpz_._mp_size = from.mpz_._mp_size;
        memcpy(&data_[0], &from.data_[0], data_.size() * sizeof(data_[0]));
        //std::copy(from.data_.begin(), from.data_.begin() + data_.size(), data_.begin());
		return *this;
	}
/*
    FixPoint<presision, decimals, IsSMem>& operator= (long from)
	{
		mpz_set_si(&mpz_, from);
		if (abs(mpz_._mp_size) > cLimbs)
			throw std::runtime_error("FixPoint overflow");
		return *this;
	}

    FixPoint<presision, decimals, IsSMem>& operator= (unsigned long from)
	{
		mpz_set_ui(&mpz_, from);
		if (abs(mpz_._mp_size) > cLimbs)
			throw std::runtime_error("FixPoint overflow");
		return *this;
	}
*/
    FixPoint<presision, decimals, IsSMem>& operator= (double from)
	{
        if (IsSMem)
            mpz_._mp_d = &data_[0];
        mpz_set_d(&mpz_, from*pow(double(10.), double(decimals)));
		if (abs(mpz_._mp_size) > cLimbs)
			throw std::runtime_error("FixPoint overflow");
		return *this;
	}

    FixPoint<presision, decimals, IsSMem>& operator= (const std::string& from)
	{
		*this = from.c_str();
		return *this;
	}

    FixPoint<presision, decimals, IsSMem>& operator= (const std::string_view& from)
	{
        char* tm = (char*)alloca(from.size() + 1);
        memcpy(tm, from.data(), from.size());
        tm[from.size()] = 0;
        *this = tm;
		return *this;
	}

    FixPoint<presision, decimals, IsSMem>& operator= (const char* from)
	{
        if (IsSMem)
            mpz_._mp_d = &data_[0];
        if (from != nullptr)
		{
			long decimal = decimals;
			const char* pointPlace = strchr(from, '.');
			const char* expPlace = strchr(from, 'E');
			if (expPlace == nullptr)
				expPlace = strchr(from, 'e');
			char buff[presision + 20];
			char* curWrite = buff;
			const char* curRead = from;
			if (pointPlace)
			{
				memcpy(curWrite, curRead, pointPlace - curRead);
				curWrite += pointPlace - curRead;
				curRead = pointPlace + 1;
				decimal -= expPlace ? (expPlace - pointPlace - 1) : (strlen(pointPlace + 1));
			}
			size_t len = (expPlace ? expPlace - curRead : strlen(curRead));
			memcpy(curWrite, curRead, len);
			curWrite += len;
			curRead += len;
			*curWrite = 0;

			if (expPlace)
				decimal += atol(expPlace + 1);

			mpz_set_str(&mpz_, buff, 10);
			std::array<mp_limb_t, cLimbs + 4> data;
			__mpz_struct mpz;
			mpz._mp_alloc = data.size();
			mpz._mp_size = 0;
			mpz._mp_d = &data[0];
			if (decimal >= 0)
			{
				mpz_mul(&mpz, &mpz_, &dividers_[decimal].mpz_);
			}
			else if (decimal < 0)
			{
				mpz_tdiv_q(&mpz, &mpz_, &dividers_[-decimal].mpz_);
			}
			std::copy(data.begin(), data.begin() + data_.size(), data_.begin());
			mpz_._mp_size = mpz._mp_size;
			if (abs(mpz_._mp_size) > cLimbs)
				throw std::runtime_error("FixPoint overflow");

		}
		else
		{
			mpz_._mp_size = 0;
		}
		return *this;
	}

	operator std::string () const
	{
		char buff[presision + 5];
		to_string(buff);
		return std::string(buff);
	}

	std::string to_string() const
	{
		return (std::string)*this;
	}

	char* to_string(char* outBuff) const
	{
        if (IsSMem)
            mpz_._mp_d = const_cast<mp_limb_t*>(&data_[0]);
        *outBuff = 0;
		std::array<mp_limb_t, cLimbs + 1> dataQ;
		__mpz_struct mpzQ;
		mpzQ._mp_alloc = dataQ.size();
		mpzQ._mp_size = 0;
		mpzQ._mp_d = &dataQ[0];
		std::array<mp_limb_t, cLimbs + 1> dataR;
		__mpz_struct mpzR;
		mpzR._mp_alloc = dataR.size();
		mpzR._mp_size = 0;
		mpzR._mp_d = &dataR[0];
		mpz_tdiv_qr(&mpzQ, &mpzR, &mpz_, &(dividers_[decimals].mpz_));
		//std::cout << dividers_[decimals].data_[0] << " " << mpzQ._mp_d[0] << std::endl;
		mpz_abs(&mpzR, &mpzR);
        if (mpz_sgn(&mpz_) < 0)
            strcpy(outBuff, "-");
        mpz_abs(&mpzQ, &mpzQ);
        mpz_get_str(outBuff + strlen(outBuff), 10, &mpzQ);
		if (mpzR._mp_size > 0)
		{
			char* curWrite = outBuff + strlen(outBuff);
			*curWrite++ = '.';
			char buff[decimals + 5];
			mpz_get_str(buff, 10, &mpzR);
			size_t rsize = strlen(buff);
			for (size_t i = 0; i < decimals - rsize; i++)
				*curWrite++ = '0';
			while(*(buff + rsize - 1) == '0')
				rsize--;
			memcpy(curWrite, buff, rsize);
			curWrite += rsize;
			*curWrite-- = 0;
		}
		return outBuff;
	}

    FixPoint<presision, decimals, IsSMem> operator+ (const FixPoint<presision, decimals, IsSMem>& other) const
	{
        if (IsSMem)
        {
            mpz_._mp_d = &data_[0];
            other.mpz_._mp_d = &other.data_[0];
        }
        FixPoint<presision, decimals, IsSMem> ret;
		mpz_add(&ret.mpz_, &mpz_, &other.mpz_);
		if (abs(ret.mpz_._mp_size) > ret.cLimbs)
			throw std::runtime_error("FixPoint overflow");
		return ret;
	}

    FixPoint<presision, decimals, IsSMem>& operator+= (const FixPoint<presision, decimals, IsSMem>& other)
	{
        if (IsSMem)
        {
            mpz_._mp_d = &data_[0];
            other.mpz_._mp_d = &other.data_[0];
        }
        mpz_add(&mpz_, &mpz_, &other.mpz_);
		if (abs(mpz_._mp_size) > cLimbs)
			throw std::runtime_error("FixPoint overflow");
		return *this;
	}

    FixPoint<presision, decimals, IsSMem> operator- (const FixPoint<presision, decimals, IsSMem>& other) const
	{
        if (IsSMem)
        {
            mpz_._mp_d = &data_[0];
            other.mpz_._mp_d = &other.data_[0];
        }
        FixPoint<presision, decimals, IsSMem> ret;
		mpz_sub(&ret.mpz_, &mpz_, &other.mpz_);
		if (abs(ret.mpz_._mp_size) > ret.cLimbs)
			throw std::runtime_error("FixPoint overflow");
		return ret;
	}

    FixPoint<presision, decimals, IsSMem>& operator-= (const FixPoint<presision, decimals, IsSMem>& other)
	{
        if (IsSMem)
        {
            mpz_._mp_d = &data_[0];
            other.mpz_._mp_d = &other.data_[0];
        }
        mpz_sub(&mpz_, &mpz_, &other.mpz_);
		if (abs(mpz_._mp_size) > cLimbs)
			throw std::runtime_error("FixPoint overflow");
		return *this;
	}

    FixPoint<presision, decimals, IsSMem> operator* (const FixPoint<presision, decimals, IsSMem>& other) const
	{
        if (IsSMem)
        {
            mpz_._mp_d = &data_[0];
            other.mpz_._mp_d = &other.data_[0];
        }
        FixPoint<presision, decimals, IsSMem> ret;

        std::array<mp_limb_t, 2 * cLimbs + 4> data;
		__mpz_struct mpz;
		mpz._mp_alloc = data.size();
		mpz._mp_size = 0;
		mpz._mp_d = &data[0];
		mpz_mul(&mpz, &mpz_, &other.mpz_);
		mpz_tdiv_q(&mpz, &mpz, &other.dividers_[decimals].mpz_);
		std::copy(data.begin(), data.begin() + ret.data_.size(), ret.data_.begin());
		ret.mpz_._mp_size = mpz._mp_size;

		//mpz_mul(&ret.mpz_, &mpz_, &other.mpz_);
		//mpz_tdiv_q(&ret.mpz_, &ret.mpz_, &other.dividers_[decimals].mpz_);
		if (abs(ret.mpz_._mp_size) > ret.cLimbs)
			throw std::runtime_error("FixPoint overflow");
		return ret;
	}

    FixPoint<presision, decimals, IsSMem>& operator*= (const FixPoint<presision, decimals, IsSMem>& other)
	{
        if (IsSMem)
        {
            mpz_._mp_d = &data_[0];
            other.mpz_._mp_d = &other.data_[0];
        }
        std::array<mp_limb_t, 2 * cLimbs + 4> data;
		__mpz_struct mpz;
		mpz._mp_alloc = data.size();
		mpz._mp_size = 0;
		mpz._mp_d = &data[0];
		mpz_mul(&mpz, &mpz_, &other.mpz_);
		mpz_tdiv_q(&mpz, &mpz, &other.dividers_[decimals].mpz_);
		std::copy(data.begin(), data.begin() + data_.size(), data_.begin());
		mpz_._mp_size = mpz._mp_size;
		if (abs(mpz_._mp_size) > cLimbs)
			throw std::runtime_error("FixPoint overflow");
		return *this;
	}

    FixPoint<presision, decimals, IsSMem> operator/ (const FixPoint<presision, decimals, IsSMem>& other) const
	{
        if (IsSMem)
        {
            mpz_._mp_d = &data_[0];
            other.mpz_._mp_d = &other.data_[0];
        }
        FixPoint<presision, decimals, IsSMem> ret;
        std::array<mp_limb_t, 2 * cLimbs + 4> data;
		__mpz_struct mpz;
		mpz._mp_alloc = data.size();
		mpz._mp_size = 0;
		mpz._mp_d = &data[0];
		mpz_mul(&mpz, &mpz_, &other.dividers_[decimals].mpz_);
        mpz_fdiv_q(&mpz, &mpz, &other.mpz_);
		std::copy(data.begin(), data.begin() + ret.data_.size(), ret.data_.begin());
		ret.mpz_._mp_size = mpz._mp_size;
		if (abs(ret.mpz_._mp_size) > ret.cLimbs)
			throw std::runtime_error("FixPoint overflow");
		return ret;
	}

    FixPoint<presision, decimals, IsSMem>& operator/= (const FixPoint<presision, decimals, IsSMem>& other)
	{
        if (IsSMem)
        {
            mpz_._mp_d = &data_[0];
            other.mpz_._mp_d = &other.data_[0];
        }
        std::array<mp_limb_t, 2 * cLimbs + 4> data;
		__mpz_struct mpz;
		mpz._mp_alloc = data.size();
		mpz._mp_size = 0;
		mpz._mp_d = &data[0];
        mpz_mul(&mpz, &mpz_, &other.dividers_[decimals].mpz_);
        mpz_fdiv_q(&mpz_, &mpz, &other.mpz_);
		if (abs(mpz_._mp_size) > cLimbs)
			throw std::runtime_error("FixPoint overflow");
		return *this;
	}

    FixPoint<presision, decimals, IsSMem> operator- () const
	{
        if (IsSMem)
            mpz_._mp_d = &data_[0];
        FixPoint<presision, decimals, IsSMem> ret;
		mpz_neg(&ret.mpz_, &mpz_);
		return ret;
	}

    bool operator< (const FixPoint<presision, decimals, IsSMem>& other) const
	{
        if (IsSMem)
        {
            mpz_._mp_d = &data_[0];
            other.mpz_._mp_d = &other.data_[0];
        }
        return mpz_cmp(&mpz_ , &other.mpz_) < 0;
	}

    bool operator> (const FixPoint<presision, decimals, IsSMem>& other) const
	{
        if (IsSMem)
        {
            mpz_._mp_d = &data_[0];
            other.mpz_._mp_d = &other.data_[0];
        }
        return mpz_cmp(&mpz_ , &other.mpz_) > 0;
	}

    bool operator== (const FixPoint<presision, decimals, IsSMem>& other) const
	{
        if (IsSMem)
        {
            mpz_._mp_d = const_cast<mp_limb_t*>(&data_[0]);
            other.mpz_._mp_d = const_cast<mp_limb_t*>(&other.data_[0]);
        }
        return mpz_cmp(&mpz_ , &other.mpz_) == 0;
	}

	bool operator< (double other) const
	{
        FixPoint<presision, decimals, IsSMem> tm(other);
        return *this < tm;
//		return mpz_cmp_si(&mpz_ , other * pow(10, decimals)) < 0;
	}

	bool operator> (double other) const
	{
        FixPoint<presision, decimals, IsSMem> tm(other);
        return *this > tm;
//		return mpz_cmp_si(&mpz_ , other * pow(10, decimals)) > 0;
	}

	bool operator== (double other) const
	{
        FixPoint<presision, decimals, IsSMem> tm(other);
        return *this == tm;
//		return mpz_cmp_si(&mpz_ , other * pow(10, decimals)) == 0;
	}

	template<typename T>
	bool operator!= (const T& other) const
	{
		return !operator==(other);
	}

	template<typename T>
	bool operator<= (const T& other) const
	{
		return !operator>(other);
	}

	template<typename T>
	bool operator>= (const T& other) const
	{
		return !operator<(other);
	}

	// round to lower value
    FixPoint<presision, decimals, IsSMem>& round(long dec)
    {
        if (IsSMem)
            mpz_._mp_d = &data_[0];
        if (dec > long(decimals))
                dec = long(decimals);
        if (dec <=  long(decimals) - long(presision))
                dec = long(decimals) - long(presision) + 1;
        std::array<mp_limb_t, cLimbs + 1> dataQ;
		__mpz_struct mpzQ;
		mpzQ._mp_alloc = dataQ.size();
		mpzQ._mp_size = 0;
		mpzQ._mp_d = &dataQ[0];
		std::array<mp_limb_t, cLimbs + 4> dataR;
		__mpz_struct mpzR;
		mpzR._mp_alloc = dataR.size();
		mpzR._mp_size = 0;
		mpzR._mp_d = &dataR[0];
		mpz_tdiv_qr(&mpzQ, &mpzR, &mpz_, &dividers_[decimals - dec].mpz_);
		mpz_sub(&mpz_, &mpz_, &mpzR);
		return *this;
	}

    // round to lower value
    FixPoint<presision, decimals, IsSMem>& round(const FixPoint<presision, decimals, IsSMem>& quant)
    {
        if (IsSMem)
            mpz_._mp_d = &data_[0];
        if (quant == zero())
            // nothing to round
            return *this;

        std::array<mp_limb_t, cLimbs + 1> dataQ;
        __mpz_struct mpzQ;
        mpzQ._mp_alloc = dataQ.size();
        mpzQ._mp_size = 0;
        mpzQ._mp_d = &dataQ[0];
        std::array<mp_limb_t, cLimbs + 4> dataR;
        __mpz_struct mpzR;
        mpzR._mp_alloc = dataR.size();
        mpzR._mp_size = 0;
        mpzR._mp_d = &dataR[0];
        mpz_tdiv_qr(&mpzQ, &mpzR, &mpz_, &quant.mpz_);
        mpz_sub(&mpz_, &mpz_, &mpzR);
        return *this;
    }

    // mathemetical round
    FixPoint<presision, decimals, IsSMem>& mround(long dec)
	{
        if (IsSMem)
            mpz_._mp_d = &data_[0];
        if (dec > long(decimals))
                dec = long(decimals);
        if (dec <=  long(decimals) - long(presision))
                dec = long(decimals) - long(presision) + 1;
        std::array<mp_limb_t, cLimbs + 1> dataQ;
		__mpz_struct mpzQ;
		mpzQ._mp_alloc = dataQ.size();
		mpzQ._mp_size = 0;
		mpzQ._mp_d = &dataQ[0];
		std::array<mp_limb_t, cLimbs + 4> dataR;
		__mpz_struct mpzR;
		mpzR._mp_alloc = dataR.size();
		mpzR._mp_size = 0;
		mpzR._mp_d = &dataR[0];
		mpz_add(&mpz_, &mpz_, &addRounders_[decimals - dec].mpz_);
		mpz_tdiv_qr(&mpzQ, &mpzR, &mpz_, &dividers_[decimals - dec].mpz_);
		mpz_sub(&mpz_, &mpz_, &mpzR);
		return *this;
	}

    // mathemetical round
    FixPoint<presision, decimals, IsSMem>& mround(const FixPoint<presision, decimals, IsSMem>& quant)
    {
        if (IsSMem)
            mpz_._mp_d = &data_[0];
        if (quant == zero())
            // nothing to round
            return *this;

        std::array<mp_limb_t, cLimbs + 1> dataQ;
        __mpz_struct mpzQ;
        mpzQ._mp_alloc = dataQ.size();
        mpzQ._mp_size = 0;
        mpzQ._mp_d = &dataQ[0];
        std::array<mp_limb_t, cLimbs + 4> dataR;
        __mpz_struct mpzR;
        mpzR._mp_alloc = dataR.size();
        mpzR._mp_size = 0;
        mpzR._mp_d = &dataR[0];
        //  TODO : needs optimization.
        FixPoint<presision, decimals, IsSMem> halfQuant = quant / "2";
        mpz_add(&mpz_, &mpz_, &halfQuant.mpz_);
        mpz_tdiv_qr(&mpzQ, &mpzR, &mpz_, &quant.mpz_);
        mpz_sub(&mpz_, &mpz_, &mpzR);
        return *this;
    }

    double get_d()
	{
        if (IsSMem)
            mpz_._mp_d = &data_[0];
        return mpz_get_d(&mpz_) / pow(10., decimals);
	}

};

template<size_t presision, size_t decimals, bool IsSMem>
typename FixPoint<presision, decimals, IsSMem>::Dividers FixPoint<presision, decimals, IsSMem>::dividers_;
template<size_t presision, size_t decimals, bool IsSMem>
typename FixPoint<presision, decimals, IsSMem>::AddRounders FixPoint<presision, decimals, IsSMem>::addRounders_;

template<size_t presision, size_t decimals, bool IsSMem = false>
inline std::basic_ostream<char>& operator<<(std::basic_ostream<char>& writer, const FixPoint<presision, decimals, IsSMem>& val)
{
	writer << (std::string)val;
	return writer;
}

template<size_t presision, size_t decimals, bool IsSMem, typename T1, typename = std::enable_if_t<std::is_convertible<T1, FixPoint<presision, decimals, IsSMem>>::value>>
FixPoint<presision, decimals, IsSMem> operator+ (const FixPoint<presision, decimals, IsSMem>& dec1, const T1& dec2)
{
    return dec1 + FixPoint<presision, decimals, IsSMem>(dec2);
}

template<size_t presision, size_t decimals, bool IsSMem, typename T1, typename = std::enable_if_t<std::is_convertible<T1, FixPoint<presision, decimals, IsSMem>>::value>>
FixPoint<presision, decimals, IsSMem> operator+ (const T1& dec1, const FixPoint<presision, decimals, IsSMem>& dec2)
{
    return FixPoint<presision, decimals, IsSMem>(dec1) + dec2;
}

template<size_t presision, size_t decimals, bool IsSMem, typename T1, typename = std::enable_if_t<std::is_convertible<T1, FixPoint<presision, decimals, IsSMem>>::value>>
FixPoint<presision, decimals, IsSMem> operator- (const FixPoint<presision, decimals, IsSMem>& dec1, const T1& dec2)
{
    return dec1 - FixPoint<presision, decimals, IsSMem>(dec2);
}

template<size_t presision, size_t decimals, bool IsSMem, typename T1, typename = std::enable_if_t<std::is_convertible<T1, FixPoint<presision, decimals, IsSMem>>::value>>
FixPoint<presision, decimals, IsSMem> operator- (const T1& dec1, const FixPoint<presision, decimals, IsSMem>& dec2)
{
    return FixPoint<presision, decimals, IsSMem>(dec1) - dec2;
}

template<size_t presision, size_t decimals, bool IsSMem, typename T1, typename = std::enable_if_t<std::is_convertible<T1, FixPoint<presision, decimals, IsSMem>>::value>>
FixPoint<presision, decimals, IsSMem> operator* (const FixPoint<presision, IsSMem, decimals>& dec1, const T1& dec2)
{
    return dec1 * FixPoint<presision, decimals, IsSMem>(dec2);
}

template<size_t presision, size_t decimals, bool IsSMem, typename T1, typename = std::enable_if_t<std::is_convertible<T1, FixPoint<presision, decimals, IsSMem>>::value>>
FixPoint<presision, decimals, IsSMem> operator* (const T1& dec1, const FixPoint<presision, decimals, IsSMem>& dec2)
{
    return FixPoint<presision, decimals, IsSMem>(dec1) * dec2;
}

template<size_t presision, size_t decimals, bool IsSMem, typename T1, typename = std::enable_if_t<std::is_convertible<T1, FixPoint<presision, decimals, IsSMem>>::value>>
FixPoint<presision, decimals, IsSMem> operator/ (const FixPoint<presision, decimals, IsSMem>& dec1, const T1& dec2)
{
    return dec1 / FixPoint<presision, decimals, IsSMem>(dec2);
}

template<size_t presision, size_t decimals, bool IsSMem, typename T1, typename = std::enable_if_t<std::is_convertible<T1, FixPoint<presision, decimals, IsSMem>>::value>>
FixPoint<presision, decimals, IsSMem> operator/ (const T1& dec1, const FixPoint<presision, decimals, IsSMem>& dec2)
{
    return FixPoint<presision, decimals, IsSMem>(dec1) / dec2;
}



