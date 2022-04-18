#pragma once

#include "fixpoint.h"

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
	tm1 = 1.01;
	decimal tm2;
	tm2 = 1.02;

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
	tm4 = 1.;
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
	assert(tm2 > tm1);
	assert(tm1 < tm2);
	assert(tm2 >= tm1);
	assert(tm1 <= tm2);
	//assert(decimal(1900.) == decimal("1900"));


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
	tm1 = "1.01";
	for (int i = 0; i < 1000000; ++i)
	{
		decimal tm22;
		tm22 = tm2;
		decimal tm3 = tm1 / tm22;
	}
	clock_gettime(CLOCK_REALTIME, &t2);

	t = t2.tv_sec * 1000000000 + t2.tv_nsec - t1.tv_sec * 1000000000 - t1.tv_nsec;
	std::cout << "FixPoint div with construct time: " << (t/1000000) << "ms." << ((t/1000) % 1000) << "."  << (t % 1000) << std::endl;

	clock_gettime(CLOCK_REALTIME, &t1);
	tm1 = "1.01";
	for (int i = 0; i < 1000000; ++i)
	{
		tmpstring = tm1;
	}
	clock_gettime(CLOCK_REALTIME, &t2);

	t = t2.tv_sec * 1000000000 + t2.tv_nsec - t1.tv_sec * 1000000000 - t1.tv_nsec;
	std::cout << "FixPoint to string time: " << (t/1000000) << "ms." << ((t/1000) % 1000) << "."  << (t % 1000) << std::endl;

	clock_gettime(CLOCK_REALTIME, &t1);
	tm1 = "1.01";
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

