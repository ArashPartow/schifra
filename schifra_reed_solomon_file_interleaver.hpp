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


#ifndef INCLUDE_SCHIFRA_REED_SOLOMON_FILE_INTERLEAVER_HPP
#define INCLUDE_SCHIFRA_REED_SOLOMON_FILE_INTERLEAVER_HPP


#include <iostream>
#include <string>

#include "schifra_reed_solomon_interleaving.hpp"
#include "schifra_fileio.hpp"


namespace schifra
{

   namespace reed_solomon
   {

      template <std::size_t block_length, std::size_t stack_size>
      class file_interleaver
      {
      public:

         file_interleaver(const std::string& input_file_name,
                          const std::string& output_file_name)
         {
            std::size_t remaining_bytes = schifra::fileio::file_size(input_file_name);

            if (0 == remaining_bytes)
            {
               std::cout << "reed_solomon::file_interleaver() - Error: input file has ZERO size." << std::endl;
               return;
            }

            std::ifstream in_stream(input_file_name.c_str(),std::ios::binary);

            if (!in_stream)
            {
               std::cout << "reed_solomon::file_interleaver() - Error: input file could not be opened." << std::endl;
               return;
            }

            std::ofstream out_stream(output_file_name.c_str(),std::ios::binary);

            if (!out_stream)
            {
               std::cout << "reed_solomon::file_interleaver() - Error: output file could not be created." << std::endl;
               return;
            }

            while (remaining_bytes >= (block_length * stack_size))
            {
               process_block(in_stream,out_stream);
               remaining_bytes -= (block_length * stack_size);
            }

            if (remaining_bytes > 0)
            {
               process_incomplete_block(in_stream,out_stream,remaining_bytes);
            }

            in_stream.close();
            out_stream.close();
         }

      private:

         inline void process_block(std::ifstream& in_stream,
                                   std::ofstream& out_stream)
         {
            for (std::size_t i = 0; i < stack_size; ++i)
            {
               in_stream.read(&block_stack_[i][0],static_cast<std::streamsize>(block_length));
            }

            interleave<char,block_length,stack_size>(block_stack_);

            for (std::size_t i = 0; i < stack_size; ++i)
            {
               out_stream.write(&block_stack_[i][0],static_cast<std::streamsize>(block_length));
            }
         }

         inline void process_incomplete_block(std::ifstream& in_stream,
                                              std::ofstream& out_stream,
                                              const std::size_t amount)
         {
            std::size_t complete_row_count = amount / block_length;
            std::size_t remainder = amount % block_length;

            for (std::size_t i = 0; i < complete_row_count; ++i)
            {
               in_stream.read(&block_stack_[i][0],static_cast<std::streamsize>(block_length));
            }

            if (remainder != 0)
            {
               in_stream.read(&block_stack_[complete_row_count][0],static_cast<std::streamsize>(remainder));
            }

            if (remainder == 0)
               interleave<char,block_length,stack_size>(block_stack_,complete_row_count);
            else
               interleave<char,block_length>(block_stack_,complete_row_count + 1,remainder);

            for (std::size_t i = 0; i < complete_row_count; ++i)
            {
               out_stream.write(&block_stack_[i][0],static_cast<std::streamsize>(block_length));
            }

            if (remainder != 0)
            {
               out_stream.write(&block_stack_[complete_row_count][0],static_cast<std::streamsize>(remainder));
            }
         }

         data_block<char,block_length> block_stack_[stack_size];

      };

      template <std::size_t block_length, std::size_t stack_size>
      class file_deinterleaver
      {
      public:

         file_deinterleaver(const std::string& input_file_name,
                            const std::string& output_file_name)
         {
            std::size_t input_file_size = schifra::fileio::file_size(input_file_name);

            if (input_file_size == 0)
            {
               std::cout << "reed_solomon::file_deinterleaver() - Error: input file has ZERO size." << std::endl;
               return;
            }

            std::ifstream in_stream(input_file_name.c_str(),std::ios::binary);

            if (!in_stream)
            {
               std::cout << "reed_solomon::file_deinterleaver() - Error: input file could not be opened." << std::endl;
               return;
            }

            std::ofstream out_stream(output_file_name.c_str(),std::ios::binary);

            if (!out_stream)
            {
               std::cout << "reed_solomon::file_deinterleaver() - Error: output file could not be created." << std::endl;
               return;
            }

            for (std::size_t i = 0; i < (input_file_size / (block_length * stack_size)); ++i)
            {
               process_block(in_stream,out_stream);
            }

            if ((input_file_size % (block_length * stack_size)) != 0)
            {
               process_incomplete_block(in_stream,out_stream,(input_file_size % (block_length * stack_size)));
            }

            in_stream.close();
            out_stream.close();
         }

      private:

         inline void process_block(std::ifstream& in_stream,
                                   std::ofstream& out_stream)
         {
            for (std::size_t i = 0; i < stack_size; ++i)
            {
               in_stream.read(&block_stack_[i][0],static_cast<std::streamsize>(block_length));
            }

            deinterleave<char,block_length,stack_size>(block_stack_);

            for (std::size_t i = 0; i < stack_size; ++i)
            {
               out_stream.write(&block_stack_[i][0],static_cast<std::streamsize>(block_length));
            }
         }

         inline void process_incomplete_block(std::ifstream& in_stream,
                                              std::ofstream& out_stream,
                                              const std::size_t amount)
         {
            std::size_t complete_row_count = amount / block_length;
            std::size_t remainder = amount % block_length;

            for (std::size_t i = 0; i < complete_row_count; ++i)
            {
               in_stream.read(&block_stack_[i][0],static_cast<std::streamsize>(block_length));
            }

            if (remainder != 0)
            {
               in_stream.read(&block_stack_[complete_row_count][0],static_cast<std::streamsize>(remainder));
            }

            if (remainder == 0)
               deinterleave<char,block_length>(block_stack_,complete_row_count);
            else
               deinterleave<char,block_length>(block_stack_,complete_row_count + 1,remainder);

            for (std::size_t i = 0; i < complete_row_count; ++i)
            {
               out_stream.write(&block_stack_[i][0],static_cast<std::streamsize>(block_length));
            }

            if (remainder != 0)
            {
               out_stream.write(&block_stack_[complete_row_count][0],static_cast<std::streamsize>(remainder));
            }
         }

         data_block<char,block_length> block_stack_[stack_size];

      };

   } // namespace reed_solomon

} // namespace schifra

#endif
