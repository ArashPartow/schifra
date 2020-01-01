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


#ifndef INCLUDE_SCHIFRA_REED_SOLOMON_BLOCK_HPP
#define INCLUDE_SCHIFRA_REED_SOLOMON_BLOCK_HPP


#include <iostream>
#include <string>

#include "schifra_galois_field.hpp"
#include "schifra_ecc_traits.hpp"


namespace schifra
{

   namespace reed_solomon
   {

      template <std::size_t code_length, std::size_t fec_length, std::size_t data_length = code_length - fec_length>
      struct block
      {
      public:

         typedef galois::field_symbol symbol_type;
         typedef traits::reed_solomon_triat<code_length,fec_length,data_length> trait;
         typedef traits::symbol<code_length> symbol;
         typedef block<code_length,fec_length,data_length> block_t;

         enum error_t
         {
            e_no_error       = 0,
            e_encoder_error0 = 1,
            e_encoder_error1 = 2,
            e_decoder_error0 = 3,
            e_decoder_error1 = 4,
            e_decoder_error2 = 5,
            e_decoder_error3 = 6,
            e_decoder_error4 = 7
         };

         block()
         : errors_detected (0),
           errors_corrected(0),
           zero_numerators (0),
           unrecoverable(false),
           error(e_no_error)
         {
            traits::validate_reed_solomon_block_parameters<code_length,fec_length,data_length>();
         }

         block(const std::string& _data, const std::string& _fec)
         : errors_detected (0),
           errors_corrected(0),
           zero_numerators (0),
           unrecoverable(false),
           error(e_no_error)
         {
            traits::validate_reed_solomon_block_parameters<code_length,fec_length,data_length>();

            for (std::size_t i = 0; i < data_length; ++i)
            {
               data[i] = static_cast<galois::field_symbol>(_data[i]);
            }

            for (std::size_t i = 0; i < fec_length; ++i)
            {
               data[i + data_length] = static_cast<galois::field_symbol>(_fec[i]);
            }
         }

         galois::field_symbol& operator[](const std::size_t& index)
         {
            return data[index];
         }

         const galois::field_symbol& operator[](const std::size_t& index) const
         {
            return data[index];
         }

         galois::field_symbol& operator()(const std::size_t& index)
         {
            return operator[](index);
         }

         galois::field_symbol& fec(const std::size_t& index)
         {
            return data[data_length + index];
         }

         bool data_to_string(std::string& data_str) const
         {
            if (data_str.length() != data_length)
            {
               return false;
            }

            for (std::size_t i = 0; i < data_length; ++i)
            {
               data_str[i] = static_cast<char>(data[i]);
            }

            return true;
         }

         bool fec_to_string(std::string& fec_str) const
         {
            if (fec_str.length() != fec_length)
            {
               return false;
            }

            for (std::size_t i = 0; i < fec_length; ++i)
            {
               fec_str[i] = static_cast<char>(data[data_length + i]);
            }

            return true;
         }

         std::string fec_to_string() const
         {
            std::string fec_str(fec_length,0x00);
            fec_to_string(fec_str);
            return fec_str;
         }

         void clear(galois::field_symbol value = 0)
         {
            for (std::size_t i = 0; i < code_length; ++i)
            {
               data[i] = value;
            }
         }

         void clear_data(galois::field_symbol value = 0)
         {
            for (std::size_t i = 0; i < data_length; ++i)
            {
               data[i] = value;
            }
         }

         void clear_fec(galois::field_symbol value = 0)
         {
            for (std::size_t i = 0; i < fec_length; ++i)
            {
               data[data_length + i] = value;
            }
         }

         void reset(galois::field_symbol value = 0)
         {
            clear(value);
            errors_detected  = 0;
            errors_corrected = 0;
            zero_numerators  = 0;
            unrecoverable    = false;
            error            = e_no_error;
         }

         template <typename BlockType>
         void copy_state(const BlockType& b)
         {
            errors_detected  = b.errors_detected;
            errors_corrected = b.errors_corrected;
            zero_numerators  = b.zero_numerators;
            unrecoverable    = b.unrecoverable;
            error            = static_cast<error_t>(b.error);
         }

         inline std::string error_as_string() const
         {
            switch (error)
            {
               case e_no_error       : return "No Error";
               case e_encoder_error0 : return "Invalid Encoder";
               case e_encoder_error1 : return "Incompatible Generator Polynomial";
               case e_decoder_error0 : return "Invalid Decoder";
               case e_decoder_error1 : return "Decoder Failure - Non-zero Syndrome";
               case e_decoder_error2 : return "Decoder Failure - Too Many Errors/Erasures";
               case e_decoder_error3 : return "Decoder Failure - Invalid Symbol Correction";
               case e_decoder_error4 : return "Decoder Failure - Invalid Codeword Correction";
               default               : return "Invalid Error Code";
            }
         }

         std::size_t  errors_detected;
         std::size_t errors_corrected;
         std::size_t  zero_numerators;
         bool           unrecoverable;
         error_t                error;
         galois::field_symbol data[code_length];
      };

      template <std::size_t code_length, std::size_t fec_length>
      inline void copy(const block<code_length,fec_length>& src_block, block<code_length,fec_length>& dest_block)
      {
         for (std::size_t index = 0; index < code_length; ++index)
         {
            dest_block.data[index] = src_block.data[index];
         }
      }

      template <typename T, std::size_t code_length, std::size_t fec_length>
      inline void copy(const T src_data[], block<code_length,fec_length>& dest_block)
      {
         for (std::size_t index = 0; index < (code_length - fec_length); ++index, ++src_data)
         {
            dest_block.data[index] = static_cast<typename block<code_length,fec_length>::symbol_type>(*src_data);
         }
      }

      template <typename T, std::size_t code_length, std::size_t fec_length>
      inline void copy(const T src_data[],
                       const std::size_t& src_length,
                       block<code_length,fec_length>& dest_block)
      {
         for (std::size_t index = 0; index < src_length; ++index, ++src_data)
         {
            dest_block.data[index] = static_cast<typename block<code_length,fec_length>::symbol_type>(*src_data);
         }
      }

      template <std::size_t code_length, std::size_t fec_length, std::size_t stack_size>
      inline void copy(const block<code_length,fec_length>  src_block_stack[stack_size],
                             block<code_length,fec_length> dest_block_stack[stack_size])
      {
         for (std::size_t row = 0; row < stack_size; ++row)
         {
            copy(src_block_stack[row], dest_block_stack[row]);
         }
      }

      template <typename T, std::size_t code_length, std::size_t fec_length, std::size_t stack_size>
      inline bool copy(const T src_data[],
                       const std::size_t src_length,
                       block<code_length,fec_length> dest_block_stack[stack_size])
      {
         const std::size_t data_length = code_length - fec_length;

         if (src_length > (stack_size * data_length))
         {
            return false;
         }

         const std::size_t row_count =  src_length / data_length;

         for (std::size_t row = 0; row < row_count; ++row, src_data += data_length)
         {
            copy(src_data, dest_block_stack[row]);
         }

         if ((src_length % data_length) != 0)
         {
            copy(src_data, src_length % data_length, dest_block_stack[row_count]);
         }

         return true;
      }

      template <typename T, std::size_t code_length, std::size_t fec_length>
      inline void full_copy(const block<code_length,fec_length>& src_block,
                            T dest_data[])
      {
         for (std::size_t i = 0; i < code_length; ++i, ++dest_data)
         {
            (*dest_data) = static_cast<T>(src_block[i]);
         }
      }

      template <typename T, std::size_t code_length, std::size_t fec_length, std::size_t stack_size>
      inline void copy(const block<code_length,fec_length> src_block_stack[stack_size],
                       T dest_data[])
      {
         const std::size_t data_length = code_length - fec_length;

         for (std::size_t i = 0; i < stack_size; ++i)
         {
            for (std::size_t j = 0; j < data_length; ++j, ++dest_data)
            {
               (*dest_data) = static_cast<T>(src_block_stack[i][j]);
            }
         }
      }

      template <std::size_t code_length, std::size_t fec_length>
      inline std::ostream& operator<<(std::ostream& os, const block<code_length,fec_length>& rs_block)
      {
         for (std::size_t i = 0; i < code_length; ++i)
         {
            os << static_cast<char>(rs_block[i]);
         }

         return os;
      }

      template <typename T, std::size_t block_length>
      struct data_block
      {
      public:

         typedef T value_type;

               T& operator[](const std::size_t index)        { return data[index]; }
         const T& operator[](const std::size_t index) const  { return data[index]; }

               T* begin()        { return data; }
         const T* begin() const  { return data; }

               T* end()          { return data + block_length; }
         const T* end()   const  { return data + block_length; }

         void clear(T value = 0)
         {
            for (std::size_t i = 0; i < block_length; ++i)
            {
               data[i] = value;
            }
         }

      private:

         T data[block_length];
      };

      template <typename T, std::size_t block_length>
      inline void copy(const data_block<T,block_length>& src_block, data_block<T,block_length>& dest_block)
      {
         for (std::size_t index = 0; index < block_length; ++index)
         {
            dest_block[index] = src_block[index];
         }
      }

      template <typename T, std::size_t block_length, std::size_t stack_size>
      inline void copy(const data_block<T,block_length>  src_block_stack[stack_size],
                             data_block<T,block_length> dest_block_stack[stack_size])
      {
         for (std::size_t row = 0; row < stack_size; ++row)
         {
            copy(src_block_stack[row], dest_block_stack[row]);
         }
      }

      template <typename T, std::size_t block_length>
      inline void full_copy(const data_block<T,block_length>& src_block, T dest_data[])
      {
         for (std::size_t i = 0; i < block_length; ++i, ++dest_data)
         {
            (*dest_data) = static_cast<T>(src_block[i]);
         }
      }

      typedef std::vector<std::size_t> erasure_locations_t;

   } // namespace reed_solomon

} // namepsace schifra

#endif
