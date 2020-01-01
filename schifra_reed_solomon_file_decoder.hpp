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


#ifndef INCLUDE_SCHIFRA_REED_SOLOMON_FILE_DECODER_HPP
#define INCLUDE_SCHIFRA_REED_SOLOMON_FILE_DECODER_HPP


#include <iostream>
#include <fstream>

#include "schifra_reed_solomon_block.hpp"
#include "schifra_reed_solomon_decoder.hpp"
#include "schifra_fileio.hpp"


namespace schifra
{

   namespace reed_solomon
   {

      template <std::size_t code_length, std::size_t fec_length, std::size_t data_length = code_length - fec_length>
      class file_decoder
      {
      public:

         typedef decoder<code_length,fec_length> decoder_type;
         typedef typename decoder_type::block_type block_type;

         file_decoder(const decoder_type& decoder,
                      const std::string& input_file_name,
                      const std::string& output_file_name)
         : current_block_index_(0)
         {
            std::size_t remaining_bytes = schifra::fileio::file_size(input_file_name);

            if (remaining_bytes == 0)
            {
               std::cout << "reed_solomon::file_decoder() - Error: input file has ZERO size." << std::endl;
               return;
            }

            std::ifstream in_stream(input_file_name.c_str(),std::ios::binary);
            if (!in_stream)
            {
               std::cout << "reed_solomon::file_decoder() - Error: input file could not be opened." << std::endl;
               return;
            }

            std::ofstream out_stream(output_file_name.c_str(),std::ios::binary);
            if (!out_stream)
            {
               std::cout << "reed_solomon::file_decoder() - Error: output file could not be created." << std::endl;
               return;
            }

            current_block_index_ = 0;

            while (remaining_bytes >= code_length)
            {
               process_complete_block(decoder,in_stream,out_stream);
               remaining_bytes -= code_length;
               current_block_index_++;
            }

            if (remaining_bytes > 0)
            {
               process_partial_block(decoder,in_stream,out_stream,remaining_bytes);
            }

            in_stream.close();
            out_stream.close();
         }

      private:

         inline void process_complete_block(const decoder_type& decoder,
                                            std::ifstream& in_stream,
                                            std::ofstream& out_stream)
         {
            in_stream.read(&buffer_[0],static_cast<std::streamsize>(code_length));
            copy<char,code_length,fec_length>(buffer_,code_length,block_);

            if (!decoder.decode(block_))
            {
               std::cout << "reed_solomon::file_decoder.process_complete_block() - Error during decoding of block " << current_block_index_ << "!" << std::endl;
               return;
            }

            for (std::size_t i = 0; i < data_length; ++i)
            {
               buffer_[i] = static_cast<char>(block_[i]);
            }

            out_stream.write(&buffer_[0],static_cast<std::streamsize>(data_length));
         }

         inline void process_partial_block(const decoder_type& decoder,
                                           std::ifstream& in_stream,
                                           std::ofstream& out_stream,
                                           const std::size_t& read_amount)
         {
            if (read_amount <= fec_length)
            {
               std::cout << "reed_solomon::file_decoder.process_partial_block() - Error during decoding of block " << current_block_index_ << "!" << std::endl;
               return;
            }

            in_stream.read(&buffer_[0],static_cast<std::streamsize>(read_amount));

            for (std::size_t i = 0; i < (read_amount - fec_length); ++i)
            {
               block_.data[i] = static_cast<typename block_type::symbol_type>(buffer_[i]);
            }

            if ((read_amount - fec_length) < data_length)
            {
               for (std::size_t i = (read_amount - fec_length); i < data_length; ++i)
               {
                  block_.data[i] = 0;
               }
            }

            for (std::size_t i = 0; i < fec_length; ++i)
            {
               block_.fec(i) = static_cast<typename block_type::symbol_type>(buffer_[(read_amount - fec_length) + i]);
            }

            if (!decoder.decode(block_))
            {
               std::cout << "reed_solomon::file_decoder.process_partial_block() - Error during decoding of block " << current_block_index_ << "!" << std::endl;
               return;
            }

            for (std::size_t i = 0; i < (read_amount - fec_length); ++i)
            {
               buffer_[i] = static_cast<char>(block_.data[i]);
            }

            out_stream.write(&buffer_[0],static_cast<std::streamsize>(read_amount - fec_length));
         }

         block_type block_;
         std::size_t current_block_index_;
         char buffer_[code_length];
      };

   } // namespace reed_solomon

} // namespace schifra

#endif
