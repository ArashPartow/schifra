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


#ifndef INCLUDE_SCHIFRA_REED_SOLOMON_SPPED_EVALUATOR_HPP
#define INCLUDE_SCHIFRA_REED_SOLOMON_SPPED_EVALUATOR_HPP


#include <cstddef>
#include <cstdio>
#include <iostream>
#include <string>

#include "schifra_galois_field.hpp"
#include "schifra_sequential_root_generator_polynomial_creator.hpp"
#include "schifra_reed_solomon_block.hpp"
#include "schifra_reed_solomon_encoder.hpp"
#include "schifra_reed_solomon_decoder.hpp"
#include "schifra_reed_solomon_file_encoder.hpp"
#include "schifra_reed_solomon_file_decoder.hpp"
#include "schifra_error_processes.hpp"
#include "schifra_utilities.hpp"


namespace schifra
{

   namespace reed_solomon
   {

      template <std::size_t code_length, std::size_t fec_length>
      void create_messages(const encoder<code_length,fec_length>& rs_encoder,
                           std::vector< block<code_length,fec_length> >& original_block_list,
                           const bool full_test_set = false)
      {
         const std::size_t data_length  = code_length - fec_length;
         std::vector<std::string> message_list;
         if (full_test_set)
         {
            for (unsigned int i = 0; i < 256; ++i)
            {
               message_list.push_back(std::string(data_length,static_cast<unsigned char>(i)));
            }
         }
         else
         {
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0x00)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0xAA)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0xA5)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0xAC)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0xCA)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0x5A)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0xCC)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0xF0)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0x0F)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0xFF)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0x92)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0x6D)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0x77)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0x7A)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0xA7)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0xE5)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0xEB)));
         }

         std::string tmp_str = std::string(data_length,static_cast<unsigned char>(0x00));

         for (std::size_t i = 0; i < data_length; ++i)
         {
            tmp_str[i] = static_cast<unsigned char>(i);
         }

         message_list.push_back(tmp_str);

         for (int i = data_length - 1; i >= 0; --i)
         {
            tmp_str[i] = static_cast<unsigned char>(i);
         }

         message_list.push_back(tmp_str);

         for (std::size_t i = 0; i < data_length; ++i)
         {
            tmp_str[i] = (((i & 0x01) == 1) ? static_cast<unsigned char>(i) : 0x00);
         }

         message_list.push_back(tmp_str);

         for (std::size_t i = 0; i < data_length; ++i)
         {
            tmp_str[i] = (((i & 0x01) == 0) ? static_cast<unsigned char>(i) : 0x00);
         }

         message_list.push_back(tmp_str);

         for (int i = data_length - 1; i >= 0; --i)
         {
            tmp_str[i] = (((i & 0x01) == 1) ? static_cast<unsigned char>(i) : 0x00);
         }

         message_list.push_back(tmp_str);

         for (int i = data_length - 1; i >= 0; --i)
         {
            tmp_str[i] = (((i & 0x01) == 0) ? static_cast<unsigned char>(i) : 0x00);
         }

         message_list.push_back(tmp_str);

         tmp_str = std::string(data_length,static_cast<unsigned char>(0x00));

         for (std::size_t i = 0; i < (data_length >> 1); ++i)
         {
            tmp_str[i] = static_cast<unsigned char>(0xFF);
         }

         message_list.push_back(tmp_str);

         tmp_str = std::string(data_length,static_cast<unsigned char>(0xFF)) ;

         for (std::size_t i = 0; i < (data_length >> 1); ++i)
         {
            tmp_str[i] = static_cast<unsigned char>(0x00);
         }

         message_list.push_back(tmp_str);

         for (std::size_t i = 0; i < message_list.size(); ++i)
         {
            block<code_length,fec_length> current_block;
            rs_encoder.encode(message_list[i],current_block);
            original_block_list.push_back(current_block);
         }
      }

      template <std::size_t field_descriptor,
                std::size_t gen_poly_index,
                std::size_t code_length,
                std::size_t fec_length,
                typename RSEncoder = encoder<code_length,fec_length>,
                typename RSDecoder = decoder<code_length,fec_length>,
                std::size_t data_length = code_length - fec_length>
      struct all_errors_decoder_speed_test
      {
      public:

         all_errors_decoder_speed_test(const std::size_t prim_poly_size, const unsigned int prim_poly[])
         {
            galois::field field(field_descriptor,prim_poly_size,prim_poly);
            galois::field_polynomial generator_polynomial(field);

            if (
                 !make_sequential_root_generator_polynomial(field,
                                                            gen_poly_index,
                                                            fec_length,
                                                            generator_polynomial)
               )
            {
               return;
            }

            RSEncoder rs_encoder(field,generator_polynomial);
            RSDecoder rs_decoder(field,gen_poly_index);

            std::vector< block<code_length,fec_length> > original_block;

            create_messages<code_length,fec_length>(rs_encoder,original_block);

            std::vector<block<code_length,fec_length> > rs_block;
            std::vector<std::size_t> block_index_list;

            for (std::size_t block_index = 0; block_index < original_block.size(); ++block_index)
            {
               for (std::size_t error_count = 1; error_count <= (fec_length >> 1); ++error_count)
               {
                  for (std::size_t start_position = 0; start_position < code_length; ++start_position)
                  {
                     block<code_length,fec_length> block = original_block[block_index];
                     corrupt_message_all_errors(block,error_count,start_position,1);
                     rs_block.push_back(block);
                     block_index_list.push_back(block_index);
                  }
               }
            }

            const std::size_t max_iterations = 100;
            std::size_t blocks_decoded       = 0;
            std::size_t block_failures       = 0;

            schifra::utils::timer timer;
            timer.start();

            for (std::size_t j = 0; j < max_iterations; ++j)
            {
               for (std::size_t i = 0; i < rs_block.size(); ++i)
               {
                  if (!rs_decoder.decode(rs_block[i]))
                  {
                     std::cout << "Decoding Failure!" << std::endl;
                     block_failures++;
                  }
                  else if (!are_blocks_equivelent(rs_block[i],original_block[block_index_list[i]]))
                  {
                     std::cout << "Error Correcting Failure!" << std::endl;
                     block_failures++;
                  }
                  else
                     blocks_decoded++;
               }
            }

            timer.stop();

            double time = timer.time();
            double mbps = ((max_iterations * rs_block.size() * data_length) * 8.0) / (1048576.0 * time);

            print_codec_properties();

            if (block_failures == 0)
               printf("Blocks decoded: %8d  Time:%8.3fsec  Rate:%8.3fMbps\n",
                      static_cast<int>(blocks_decoded),
                      time,
                      mbps);
            else
               std::cout << "Blocks decoded: " << blocks_decoded << "\tDecode Failures: " << block_failures <<"\tTime: " << time <<"sec\tRate: " << mbps << "Mbps" << std::endl;
         }

         void print_codec_properties()
         {
            printf("[All Errors Test] Codec: RS(%03d,%03d,%03d) ",
                   static_cast<int>(code_length),
                   static_cast<int>(data_length),
                   static_cast<int>(fec_length));
         }
      };

      template <std::size_t field_descriptor,
                std::size_t gen_poly_index,
                std::size_t code_length,
                std::size_t fec_length,
                typename RSEncoder = encoder<code_length,fec_length>,
                typename RSDecoder = decoder<code_length,fec_length>,
                std::size_t data_length = code_length - fec_length>
      struct all_erasures_decoder_speed_test
      {
      public:

         all_erasures_decoder_speed_test(const std::size_t prim_poly_size, const unsigned int prim_poly[])
         {
            galois::field field(field_descriptor,prim_poly_size,prim_poly);
            galois::field_polynomial generator_polynomial(field);

            if (
                 !make_sequential_root_generator_polynomial(field,
                                                            gen_poly_index,
                                                            fec_length,
                                                            generator_polynomial)
               )
            {
               return;
            }

            RSEncoder rs_encoder(field,generator_polynomial);
            RSDecoder rs_decoder(field,gen_poly_index);

            std::vector< block<code_length,fec_length> > original_block;

            create_messages<code_length,fec_length>(rs_encoder,original_block);

            std::vector<block<code_length,fec_length> > rs_block;
            std::vector<erasure_locations_t> erasure_list;
            std::vector<std::size_t> block_index_list;

            for (std::size_t block_index = 0; block_index < original_block.size(); ++block_index)
            {
               for (std::size_t erasure_count = 1; erasure_count <= fec_length; ++erasure_count)
               {
                  for (std::size_t start_position = 0; start_position < code_length; ++start_position)
                  {
                     block<code_length,fec_length> block = original_block[block_index];
                     erasure_locations_t erasures;
                     corrupt_message_all_erasures(block,erasures,erasure_count,start_position,1);

                     if (erasure_count != erasures.size())
                     {
                        std::cout << "all_erasures_decoder_speed_test() - Failed to properly generate erasures list.  Details:";
                        std::cout << "(" << block_index << "," << erasure_count << "," << start_position << ")" << std::endl;
                     }

                     rs_block.push_back(block);
                     erasure_list.push_back(erasures);
                     block_index_list.push_back(block_index);
                  }
               }
            }

            const std::size_t max_iterations = 100;
            std::size_t blocks_decoded       =   0;
            std::size_t block_failures       =   0;

            schifra::utils::timer timer;
            timer.start();

            for (std::size_t j = 0; j < max_iterations; ++j)
            {
               for (std::size_t i = 0; i < rs_block.size(); ++i)
               {
                  if (!rs_decoder.decode(rs_block[i],erasure_list[i]))
                  {
                     std::cout << "Decoding Failure!" << std::endl;
                     block_failures++;
                  }
                  else if (!are_blocks_equivelent(rs_block[i],original_block[block_index_list[i]]))
                  {
                     std::cout << "Error Correcting Failure!" << std::endl;
                     block_failures++;
                  }
                  else
                     blocks_decoded++;
               }
            }

            timer.stop();

            double time = timer.time();
            double mbps = ((max_iterations * rs_block.size() * data_length) * 8.0) / (1048576.0 * time);

            print_codec_properties();

            if (block_failures == 0)
               printf("Blocks decoded: %8d  Time:%8.3fsec  Rate:%8.3fMbps\n",
                      static_cast<int>(blocks_decoded),
                      time,
                      mbps);
            else
               std::cout << "Blocks decoded: " << blocks_decoded << "\tDecode Failures: " << block_failures <<"\tTime: " << time <<"sec\tRate: " << mbps << "Mbps" << std::endl;
         }

         void print_codec_properties()
         {
            printf("[All Erasures Test] Codec: RS(%03d,%03d,%03d) ",
                   static_cast<int>(code_length),
                   static_cast<int>(data_length),
                   static_cast<int>(fec_length));
         }

      };

      void speed_test_00()
      {
         all_errors_decoder_speed_test<8,120,255,  2>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_errors_decoder_speed_test<8,120,255,  4>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_errors_decoder_speed_test<8,120,255,  6>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_errors_decoder_speed_test<8,120,255,  8>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_errors_decoder_speed_test<8,120,255, 10>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_errors_decoder_speed_test<8,120,255, 12>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_errors_decoder_speed_test<8,120,255, 14>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_errors_decoder_speed_test<8,120,255, 16>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_errors_decoder_speed_test<8,120,255, 18>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_errors_decoder_speed_test<8,120,255, 20>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_errors_decoder_speed_test<8,120,255, 32>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_errors_decoder_speed_test<8,120,255, 48>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_errors_decoder_speed_test<8,120,255, 64>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_errors_decoder_speed_test<8,120,255, 80>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_errors_decoder_speed_test<8,120,255, 96>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_errors_decoder_speed_test<8,120,255,128>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
      }

      void speed_test_01()
      {
         all_erasures_decoder_speed_test<8,120,255,  2>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_erasures_decoder_speed_test<8,120,255,  4>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_erasures_decoder_speed_test<8,120,255,  6>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_erasures_decoder_speed_test<8,120,255,  8>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_erasures_decoder_speed_test<8,120,255, 10>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_erasures_decoder_speed_test<8,120,255, 12>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_erasures_decoder_speed_test<8,120,255, 14>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_erasures_decoder_speed_test<8,120,255, 16>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_erasures_decoder_speed_test<8,120,255, 18>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_erasures_decoder_speed_test<8,120,255, 20>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_erasures_decoder_speed_test<8,120,255, 32>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_erasures_decoder_speed_test<8,120,255, 48>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_erasures_decoder_speed_test<8,120,255, 64>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_erasures_decoder_speed_test<8,120,255, 80>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_erasures_decoder_speed_test<8,120,255, 96>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
         all_erasures_decoder_speed_test<8,120,255,128>(galois::primitive_polynomial_size06,galois::primitive_polynomial06);
      }

   } // namespace reed_solomon

} // namespace schifra

#endif
