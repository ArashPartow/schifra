/*
(**************************************************************************)
(*                                                                        *)
(*                                Schifra                                 *)
(*                Reed-Solomon Error Correcting Code Library              *)
(*                                                                        *)
(* Release Version 0.0.1                                                  *)
(* http://www.schifra.com                                                 *)
(* Copyright (c) 2000-2010 Arash Partow, All Rights Reserved.             *)
(*                                                                        *)
(* The Schifra Reed-Solomon error correcting code library and all its     *)
(* components are supplied under the terms of the General Schifra License *)
(* agreement. The contents of the Schifra Reed-Solomon error correcting   *)
(* code library and all its components may not be copied or disclosed     *)
(* except in accordance with the terms of that agreement.                 *)
(*                                                                        *)
(* URL: http://www.schifra.com/license.html                               *)
(*                                                                        *)
(**************************************************************************)
*/


#ifndef INCLUDE_SCHIFRA_UTILITES_HPP
#define INCLUDE_SCHIFRA_UTILITES_HPP


#include <cstddef>

#ifdef WIN32
   #include <windows.h>
#else
  #include <sys/time.h>
  #include <sys/types.h>
#endif

namespace schifra
{

   namespace utils
   {

      const std::size_t high_bits_in_char[256] = {
                                                   0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
                                                   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
                                                   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
                                                   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                                                   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
                                                   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                                                   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                                                   3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
                                                   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
                                                   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                                                   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                                                   3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
                                                   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                                                   3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
                                                   3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
                                                   4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
                                                 };

      template<typename T>
      inline std::size_t hamming_distance_element(const T v1, const T v2)
      {
         std::size_t distance = 0;
         const unsigned char* it1 = reinterpret_cast<const unsigned char*>(&v1);
         const unsigned char* it2 = reinterpret_cast<const unsigned char*>(&v2);
         for (std::size_t i = 0; i < sizeof(T); ++i, ++it1, ++it2)
         {
            distance += high_bits_in_char[((*it1) ^ (*it2)) & 0xFF];
         }
         return distance;
      }

      inline std::size_t hamming_distance(const unsigned char data1[], const unsigned char data2[], const std::size_t length)
      {
         std::size_t distance = 0;
         const unsigned char* it1 = data1;
         const unsigned char* it2 = data2;
         for (std::size_t i = 0; i < length; ++i, ++it1, ++it2)
         {
            distance += high_bits_in_char[((*it1) ^ (*it2)) & 0xFF];
         }
         return distance;
      }

      template<typename ForwardIterator>
      inline std::size_t hamming_distance(ForwardIterator it1_begin, ForwardIterator it2_begin, ForwardIterator it1_end)
      {
         std::size_t distance = 0;
         ForwardIterator it1 = it1_begin;
         ForwardIterator it2 = it2_begin;
         for (; it1 != it1_end; ++it1, ++it2)
         {
            distance += hamming_distance_element(*it1,*it2);
         }
         return distance;
      }

      #ifdef WIN32

      class timer
      {
      public:

         timer()      { QueryPerformanceFrequency(&clock_frequency); }
         void start() { QueryPerformanceCounter(&start_time);        }
         void stop()  { QueryPerformanceCounter(&stop_time);         }

         double time()
         {
            return (1.0 *(stop_time.QuadPart - start_time.QuadPart)) / (1.0 * clock_frequency.QuadPart);
         }

       private:

        LARGE_INTEGER start_time;
        LARGE_INTEGER stop_time;
        LARGE_INTEGER clock_frequency;
      };

      #else

      class timer
      {
      public:

         void start() { gettimeofday(&start_time, 0); }
         void stop()  { gettimeofday(&stop_time,  0); }

         double time()
         {
            double diff = (stop_time.tv_sec - start_time.tv_sec) * 1000000.0;
            if (stop_time.tv_usec > start_time.tv_usec)
               diff += (1.0 * (stop_time.tv_usec - start_time.tv_usec));
            else if (stop_time.tv_usec < start_time.tv_usec)
               diff -= (1.0 * (start_time.tv_usec - stop_time.tv_usec));

            return (diff / 1000000.0);
         }

       private:

        struct timeval start_time;
        struct timeval stop_time;
        struct timeval clock_frequency;
      };

      #endif

   } // namespace utils

} // namespace schifra


#endif
