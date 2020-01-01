/*
(**************************************************************************)
(*                                                                        *)
(*                                Schifra                                 *)
(*                Reed-Solomon Error Correcting Code Library              *)
(*                                                                        *)
(* Release Version 0.0.1                                                  *)
(* http://www.schifra.com                                                 *)
(* Copyright (c) 2000-2020 Arash Partow, All Rights Reserved.             *)
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

#if defined(_WIN32) || defined(__WIN32__) || defined(WIN32)
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

      template <typename T>
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

      template <typename ForwardIterator>
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

      class timer
      {
      public:

         #if defined(_WIN32) || defined(__WIN32__) || defined(WIN32)
            timer()
            : in_use_(false)
            {
               QueryPerformanceFrequency(&clock_frequency_);
            }

            inline void start()
            {
               in_use_ = true;
               QueryPerformanceCounter(&start_time_);
            }

            inline void stop()
            {
               QueryPerformanceCounter(&stop_time_);
               in_use_ = false;
            }

            inline double time() const
            {
               return (1.0 * (stop_time_.QuadPart - start_time_.QuadPart)) / (1.0 * clock_frequency_.QuadPart);
            }

         #else

            timer()
            : in_use_(false)
            {
               start_time_.tv_sec  = 0;
               start_time_.tv_usec = 0;
               stop_time_.tv_sec   = 0;
               stop_time_.tv_usec  = 0;
            }

            inline void start()
            {
               in_use_ = true;
               gettimeofday(&start_time_,0);
            }

            inline void stop()
            {
               gettimeofday(&stop_time_, 0);
               in_use_ = false;
            }

            inline unsigned long long int usec_time() const
            {
               if (!in_use_)
               {
                  if (stop_time_.tv_sec >= start_time_.tv_sec)
                  {
                     return 1000000 * (stop_time_.tv_sec  - start_time_.tv_sec ) +
                                      (stop_time_.tv_usec - start_time_.tv_usec);
                  }
                  else
                     return std::numeric_limits<unsigned long long int>::max();
               }
               else
                  return std::numeric_limits<unsigned long long int>::max();
            }

            inline double time() const
            {
               return usec_time() * 0.000001;
            }

         #endif

            inline bool in_use() const
            {
               return in_use_;
            }

      private:

            bool in_use_;

         #if defined(_WIN32) || defined(__WIN32__) || defined(WIN32)
            LARGE_INTEGER start_time_;
            LARGE_INTEGER stop_time_;
            LARGE_INTEGER clock_frequency_;
         #else
            struct timeval start_time_;
            struct timeval stop_time_;
         #endif
      };

   } // namespace utils

} // namespace schifra


#endif
