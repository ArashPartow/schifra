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


#ifndef INCLUDE_SCHIFRA_REED_SOLOMON_FILE_ENCODER_HPP
#define INCLUDE_SCHIFRA_REED_SOLOMON_FILE_ENCODER_HPP


#include <cstring>
#include <iostream>
#include <fstream>

#include "schifra_reed_solomon_block.hpp"
#include "schifra_reed_solomon_encoder.hpp"
#include "schifra_fileio.hpp"


namespace schifra
{

   namespace reed_solomon
   {

      template <std::size_t code_length, std::size_t fec_length, std::size_t data_length = code_length - fec_length>
      class file_encoder
      {
      public:

         typedef encoder<code_length,fec_length> encoder_type;
         typedef typename encoder_type::block_type block_type;

         file_encoder(const encoder_type& encoder,
                      const std::string& input_file_name,
                      const std::string& output_file_name)
         {
            std::size_t remaining_bytes = schifra::fileio::file_size(input_file_name);
            if (remaining_bytes == 0)
            {
               std::cout << "reed_solomon::file_encoder() - Error: input file has ZERO size." << std::endl;
               return;
            }

            std::ifstream in_stream(input_file_name.c_str(),std::ios::binary);
            if (!in_stream)
            {
               std::cout << "reed_solomon::file_encoder() - Error: input file could not be opened." << std::endl;
               return;
            }

            std::ofstream out_stream(output_file_name.c_str(),std::ios::binary);
            if (!out_stream)
            {
               std::cout << "reed_solomon::file_encoder() - Error: output file could not be created." << std::endl;
               return;
            }

            std::memset(data_buffer_,0,sizeof(data_buffer_));
            std::memset(fec_buffer_ ,0,sizeof(fec_buffer_ ));

            while (remaining_bytes >= data_length)
            {
               process_block(encoder,in_stream,out_stream,data_length);
               remaining_bytes -= data_length;
            }

            if (remaining_bytes > 0)
            {
               process_block(encoder,in_stream,out_stream,remaining_bytes);
            }

            in_stream.close();
            out_stream.close();
         }

      private:

         inline void process_block(const encoder_type& encoder,
                                   std::ifstream& in_stream,
                                   std::ofstream& out_stream,
                                   const std::size_t& read_amount)
         {
            in_stream.read(&data_buffer_[0],static_cast<std::streamsize>(read_amount));
            for (std::size_t i = 0; i < read_amount; ++i)
            {
               block_.data[i] = (data_buffer_[i] & 0xFF);
            }

            if (read_amount < data_length)
            {
               for (std::size_t i = read_amount; i < data_length; ++i)
               {
                  block_.data[i] = 0x00;
               }
            }

            if (!encoder.encode(block_))
            {
               std::cout << "reed_solomon::file_encoder.process_block() - Error during encoding of block!" << std::endl;
               return;
            }

            for (std::size_t i = 0; i < fec_length; ++i)
            {
               fec_buffer_[i] = static_cast<char>(block_.fec(i) & 0xFF);
            }

            out_stream.write(&data_buffer_[0],static_cast<std::streamsize>(read_amount));
            out_stream.write(&fec_buffer_[0],fec_length);
         }

         block_type block_;
         char data_buffer_[data_length];
         char fec_buffer_[fec_length];
      };

   } // namespace reed_solomon

} // namespace schifra

#endif
