#pragma once

#include <array>
#include <exception>
#include <algorithm>
#include <gmp.h>
#include <iostream>

//	presision - total significant digits
//	decimals - digits after dot
template<size_t presision, size_t decimals>
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
			mpz_sub_ui(&mpz_, &mpz_, 1);
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
	__mpz_struct mpz_;

public:
	FixPoint(const char* strVal = nullptr)
	{
		mpz_._mp_alloc = data_.size();
		mpz_._mp_size = 0;
		mpz_._mp_d = &data_[0];
		//mpz_init(&mpz_);
		*this = strVal;
	}

	~FixPoint()
	{
	}

	FixPoint<presision, decimals>& operator= (const FixPoint<presision, decimals>& from)
	{
		mpz_._mp_size = from.mpz_._mp_size;
		data_ = from.data_;
		return *this;
	}

	FixPoint<presision, decimals>& operator= (unsigned long from)
	{
		mpz_set_ui(&mpz_, from);
		if (mpz_._mp_size > cLimbs)
			throw std::runtime_error("FixPoint overflow");
		return *this;
	}

	FixPoint<presision, decimals>& operator= (const char* from)
	{
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
			if (mpz_._mp_size > cLimbs)
				throw std::runtime_error("FixPoint overflow");

		}
		else
		{
			mpz_._mp_size = 0;
		}
		return *this;
	}

	operator std::string ()
	{
		char buff[presision + 5];
		to_string(buff);
		return std::string(buff);
	}

	char* to_string(char* outBuff)
	{
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
		mpz_tdiv_qr(&mpzQ, &mpzR, &mpz_, &dividers_[decimals].mpz_);
		mpz_get_str(outBuff, 10, &mpzQ);
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

	FixPoint<presision, decimals> operator+ (const FixPoint<presision, decimals>& other)
	{
		FixPoint<presision, decimals> ret;
		mpz_add(&ret.mpz_, &mpz_, &other.mpz_);
		if (ret.mpz_._mp_size > ret.cLimbs)
			throw std::runtime_error("FixPoint overflow");
		return ret;
	}

	FixPoint<presision, decimals>& operator+= (const FixPoint<presision, decimals>& other)
	{
		mpz_add(&mpz_, &mpz_, &other.mpz_);
		if (mpz_._mp_size > cLimbs)
			throw std::runtime_error("FixPoint overflow");
		return *this;
	}

	FixPoint<presision, decimals> operator- (const FixPoint<presision, decimals>& other)
	{
		FixPoint<presision, decimals> ret;
		mpz_sub(&ret->mpz_, &mpz_, &other.mpz_);
		if (ret.mpz_._mp_size < -ret.cLimbs)
			throw std::runtime_error("FixPoint overflow");
		return ret;
	}

	FixPoint<presision, decimals>& operator-= (const FixPoint<presision, decimals>& other)
	{
		mpz_sub(&mpz_, &mpz_, &other.mpz_);
		if (mpz_._mp_size < -cLimbs)
			throw std::runtime_error("FixPoint overflow");
		return *this;
	}

	FixPoint<presision, decimals> operator* (const FixPoint<presision, decimals>& other)
	{
		FixPoint<presision, decimals> ret;

		std::array<mp_limb_t, cLimbs + other.cLimbs + 4> data;
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
		if (ret.mpz_._mp_size > ret.cLimbs)
			throw std::runtime_error("FixPoint overflow");
		return ret;
	}

	FixPoint<presision, decimals>& operator*= (const FixPoint<presision, decimals>& other)
	{
		std::array<mp_limb_t, cLimbs + other.cLimbs + 4> data;
		__mpz_struct mpz;
		mpz._mp_alloc = data.size();
		mpz._mp_size = 0;
		mpz._mp_d = &data[0];
		mpz_mul(&mpz, &mpz_, &other.mpz_);
		mpz_tdiv_q(&mpz, &mpz, &other.dividers_[decimals].mpz_);
		std::copy(data.begin(), data.begin() + data_.size(), data_.begin());
		mpz_._mp_size = mpz._mp_size;
		if (mpz_._mp_size > cLimbs)
			throw std::runtime_error("FixPoint overflow");
		return *this;
	}

	FixPoint<presision, decimals> operator/ (const FixPoint<presision, decimals>& other)
	{
		FixPoint<presision, decimals> ret;
		std::array<mp_limb_t, cLimbs + other.cLimbs + 4> data;
		__mpz_struct mpz;
		mpz._mp_alloc = data.size();
		mpz._mp_size = 0;
		mpz._mp_d = &data[0];
		mpz_mul(&mpz, &mpz_, &other.dividers_[decimals].mpz_);
		mpz_tdiv_q(&mpz, &mpz, &other.mpz_);
		std::copy(data.begin(), data.begin() + ret.data_.size(), ret.data_.begin());
		ret.mpz_._mp_size = mpz._mp_size;
		if (ret.mpz_._mp_size > ret.cLimbs)
			throw std::runtime_error("FixPoint overflow");
		return ret;
	}

	FixPoint<presision, decimals>& operator/= (const FixPoint<presision, decimals>& other)
	{
		std::array<mp_limb_t, cLimbs + other.cLimbs + 4> data;
		__mpz_struct mpz;
		mpz._mp_alloc = data.size();
		mpz._mp_size = 0;
		mpz._mp_d = &data[0];
		mpz_mul(&mpz, &mpz_, &other.divider_.mpz_);
		mpz_tdiv_q(&mpz, &mpz, &other.mpz_);
		std::copy(data.begin(), data.begin() + data_.size(), data_.begin());
		mpz_._mp_size = mpz._mp_size;
		if (mpz_._mp_size > cLimbs)
			throw std::runtime_error("FixPoint overflow");
		return *this;
	}

	bool operator< (const FixPoint<presision, decimals>& other)
	{
		return mpz_cmp(&mpz_ , &other.mpz_) < 0;
	}

	bool operator> (const FixPoint<presision, decimals>& other)
	{
		return mpz_cmp(&mpz_ , &other.mpz_) > 0;
	}

	bool operator== (const FixPoint<presision, decimals>& other)
	{
		return mpz_cmp(&mpz_ , &other.mpz_) == 0;
	}

	// round to lower value
	void round(long dec)
	{
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
	}

	// mathemetical round
	void mround(long dec)
	{
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
	}

};

template<size_t presision, size_t decimals>
typename FixPoint<presision, decimals>::Dividers FixPoint<presision, decimals>::dividers_;
template<size_t presision, size_t decimals>
typename FixPoint<presision, decimals>::AddRounders FixPoint<presision, decimals>::addRounders_;

//#define ADDITIONAL_TESTS

#define MP_INT HIDEN_MP_TYPE
#include <gmpxx.h>
#undef MP_INT
#include <cassert>
#ifdef ADDITIONAL_TESTS
#include "../gmpwrap.h"
extern "C" {
#include <mpdecimal.h>
}
#endif	//#ifdef ADDITIONAL_TESTS


inline void FixPointPerfTest()
{
	//cout << sizeof(unsigned long int) << endl;
	timespec t1;
	timespec t2;
	long long t;

	using decimal = FixPoint<37, 19>;
	decimal tm1;
	tm1 = 10100000000000000000ul;
	decimal tm2;
	tm2 = 10200000000000000000ul;

	//	from string tests
	decimal tm3;
	tm3 = "1.02";
	assert(tm2 == tm3);
	tm3 = "102e-2";
	assert(tm2 == tm3);
	tm3 = "102E-2";
	assert(tm2 == tm3);
	tm3 = "10.2e-1";
	assert(tm2 == tm3);
	tm3 = "0.0102e2";
	assert(tm2 == tm3);
	tm3 = "1";
	decimal tm4;
	tm4 = 10000000000000000000ul;
	assert(tm4 == tm3);
	tm3 = "1.0302";
	assert(tm1 * tm2 == tm3);
	assert(tm2 / tm1 == "1.00990099009900990095");
	tm3 = "1.0099009900990099009";
	tm3.round(5);
	assert(tm3 == "1.0099");
	tm3 = "1.0099009900990099009";
	tm3.round(-1);
	assert(tm3 == "0");
	tm3 = "1.0099009900990099009";
	tm3.round(1);
	assert(tm3 == "1");
	tm3 = "1.0099009900990099009";
	tm3.round(0);
	assert(tm3 == "1");
	tm3 = "1.0099009900990099009";
	tm3.mround(3);
	assert(tm3 == "1.01");
	tm3 = "1.0099009900990099009";
	tm3.mround(-1);
	assert(tm3 == "0");
	tm3 = "1.0099009900990099009";
	tm3.mround(0);
	assert(tm3 == "1");
	tm3 = "1.0099004900990099009";
	tm3.mround(6);
	assert(tm3 == "1.0099");


	// Performance tests.
	std::string tmpstring;
	////////////
	/// FixPoint
	clock_gettime(CLOCK_REALTIME, &t1);
	for (int i = 0; i < 1000000; ++i)
	{
		decimal tm3 = tm1 + tm2;
	}
	clock_gettime(CLOCK_REALTIME, &t2);
	t = t2.tv_sec * 1000000000 + t2.tv_nsec - t1.tv_sec * 1000000000 - t1.tv_nsec;
	std::cout << "FixPoint sum time: " << (t/1000000) << "ms." << ((t/1000) % 1000) << "."  << (t % 1000) << std::endl;

	clock_gettime(CLOCK_REALTIME, &t1);
	for (int i = 0; i < 1000000; ++i)
	{
		decimal tm3 = tm1 * tm2;
	}
	clock_gettime(CLOCK_REALTIME, &t2);
	t = t2.tv_sec * 1000000000 + t2.tv_nsec - t1.tv_sec * 1000000000 - t1.tv_nsec;
	std::cout << "FixPoint mul time: " << (t/1000000) << "ms." << ((t/1000) % 1000) << "."  << (t % 1000) << std::endl;

	clock_gettime(CLOCK_REALTIME, &t1);
	for (int i = 0; i < 1000000; ++i)
	{
		decimal tm3 = tm1 / tm2;
	}
	clock_gettime(CLOCK_REALTIME, &t2);
	t = t2.tv_sec * 1000000000 + t2.tv_nsec - t1.tv_sec * 1000000000 - t1.tv_nsec;
	std::cout << "FixPoint div time: " << (t/1000000) << "ms." << ((t/1000) % 1000) << "."  << (t % 1000) << std::endl;

	clock_gettime(CLOCK_REALTIME, &t1);
	tm1 = 10100000000000000000ul;
	for (int i = 0; i < 1000000; ++i)
	{
		decimal tm2;
		tm2 = 10200000000000000000ul;
		decimal tm3 = tm1 / tm2;
	}
	clock_gettime(CLOCK_REALTIME, &t2);

	t = t2.tv_sec * 1000000000 + t2.tv_nsec - t1.tv_sec * 1000000000 - t1.tv_nsec;
	std::cout << "FixPoint div with construct time: " << (t/1000000) << "ms." << ((t/1000) % 1000) << "."  << (t % 1000) << std::endl;

	clock_gettime(CLOCK_REALTIME, &t1);
	tm1 = 10100000000000000000ul;
	for (int i = 0; i < 1000000; ++i)
	{
		tmpstring = tm1;
	}
	clock_gettime(CLOCK_REALTIME, &t2);

	t = t2.tv_sec * 1000000000 + t2.tv_nsec - t1.tv_sec * 1000000000 - t1.tv_nsec;
	std::cout << "FixPoint to string time: " << (t/1000000) << "ms." << ((t/1000) % 1000) << "."  << (t % 1000) << std::endl;

	clock_gettime(CLOCK_REALTIME, &t1);
	tm1 = 10100000000000000000ul;
	for (int i = 0; i < 1000000; ++i)
	{
		tm1 = "1.0018";
	}
	clock_gettime(CLOCK_REALTIME, &t2);

	t = t2.tv_sec * 1000000000 + t2.tv_nsec - t1.tv_sec * 1000000000 - t1.tv_nsec;
	std::cout << "FixPoint from string time: " << (t/1000000) << "ms." << ((t/1000) % 1000) << "."  << (t % 1000) << std::endl;

	//////////////////////
	/// GMP MPQ
	mpq_class tz1(101, 100), tz2(102, 100);

	clock_gettime(CLOCK_REALTIME, &t1);
	for (int i = 0; i < 1000000; ++i)
	{
		mpq_class tz3 = tz1 + tz2;
	}
	clock_gettime(CLOCK_REALTIME, &t2);

	t = t2.tv_sec * 1000000000 + t2.tv_nsec - t1.tv_sec * 1000000000 - t1.tv_nsec;
	std::cout << "MPQ sum time: " << (t/1000000) << "ms." << ((t/1000) % 1000) << "."  << (t % 1000) << std::endl;

	clock_gettime(CLOCK_REALTIME, &t1);
	for (int i = 0; i < 1000000; ++i)
	{
		mpq_class tz3 = tz1 * tz2;
	}
	clock_gettime(CLOCK_REALTIME, &t2);

	t = t2.tv_sec * 1000000000 + t2.tv_nsec - t1.tv_sec * 1000000000 - t1.tv_nsec;
	std::cout << "MPQ mul time: " << (t/1000000) << "ms." << ((t/1000) % 1000) << "."  << (t % 1000) << std::endl;

	clock_gettime(CLOCK_REALTIME, &t1);
	for (int i = 0; i < 1000000; ++i)
	{
		mpq_class tz3 = tz1 / tz2;
	}
	clock_gettime(CLOCK_REALTIME, &t2);

	t = t2.tv_sec * 1000000000 + t2.tv_nsec - t1.tv_sec * 1000000000 - t1.tv_nsec;
	std::cout << "MPQ div time: " << (t/1000000) << "ms." << ((t/1000) % 1000) << "."  << (t % 1000) << std::endl;

	clock_gettime(CLOCK_REALTIME, &t1);
	tz1 = mpq_class (101, 100);
	for (int i = 0; i < 1000000; ++i)
	{
		mpq_class tz2(102, 100);
		mpq_class tz3 = tz1 / tz2;
	}
	clock_gettime(CLOCK_REALTIME, &t2);

	t = t2.tv_sec * 1000000000 + t2.tv_nsec - t1.tv_sec * 1000000000 - t1.tv_nsec;
	std::cout << "MPQ div with construct time: " << (t/1000000) << "ms." << ((t/1000) % 1000) << "."  << (t % 1000) << std::endl;

#ifdef ADDITIONAL_TESTS
	clock_gettime(CLOCK_REALTIME, &t1);
	for (int i = 0; i < 1000000; ++i)
	{
		tmpstring = to_string(tz1);
	}
	clock_gettime(CLOCK_REALTIME, &t2);

	t = t2.tv_sec * 1000000000 + t2.tv_nsec - t1.tv_sec * 1000000000 - t1.tv_nsec;
	std::cout << "MPQ to string time: " << (t/1000000) << "ms." << ((t/1000) % 1000) << "."  << (t % 1000) << std::endl;

	clock_gettime(CLOCK_REALTIME, &t1);
	for (int i = 0; i < 1000000; ++i)
	{
		tz1 = to_mpq("1.0018");
	}
	clock_gettime(CLOCK_REALTIME, &t2);

	t = t2.tv_sec * 1000000000 + t2.tv_nsec - t1.tv_sec * 1000000000 - t1.tv_nsec;
	std::cout << "MPQ from string time: " << (t/1000000) << "ms." << ((t/1000) % 1000) << "."  << (t % 1000) << std::endl;

	/////////////
	/// MpDecimal compare
	mpd_context_t* context = new mpd_context_t();
	mpd_init(context, 128);
	context->traps = 0;
	mpd_t *td1;
	mpd_t *td2;
	mpd_t *td3;
	uint32_t status{0};
	td1 = mpd_qnew();
	mpd_qset_string(td1, "1.01", context, &status);
	td2 = mpd_qnew();
	mpd_qset_string(td2, "1.02", context, &status);
	td3 = mpd_qnew();
	mpd_qmul(td3, td1, td2, context, &status);
	std::cout << mpd_qformat(td3, "f", context, &status) << std::endl;

	clock_gettime(CLOCK_REALTIME, &t1);
	for (int i = 0; i < 1000000; ++i)
	{
		mpd_del(td3);
		td3 = mpd_qnew();
		mpd_qadd(td3, td1, td2, context, &status);
	}
	clock_gettime(CLOCK_REALTIME, &t2);
	mpd_del(td1);
	mpd_del(td2);
	mpd_del(td3);

	t = t2.tv_sec * 1000000000 + t2.tv_nsec - t1.tv_sec * 1000000000 - t1.tv_nsec;
	std::cout << "MPD sum time: " << (t/1000000) << "ms." << ((t/1000) % 1000) << "."  << (t % 1000) << std::endl;

#endif	//#ifdef ADDITIONAL_TESTS



}

