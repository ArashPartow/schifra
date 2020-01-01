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


#ifndef INCLUDE_SCHIFRA_ECC_TRAITS_HPP
#define INCLUDE_SCHIFRA_ECC_TRAITS_HPP


namespace schifra
{
   namespace traits
   {

      template <std::size_t code_length> struct symbol;
                                        /* bits per symbol */
      template <> struct symbol<    3>  { enum {size =  2}; };
      template <> struct symbol<    7>  { enum {size =  3}; };
      template <> struct symbol<   15>  { enum {size =  4}; };
      template <> struct symbol<   31>  { enum {size =  5}; };
      template <> struct symbol<   63>  { enum {size =  6}; };
      template <> struct symbol<  127>  { enum {size =  7}; };
      template <> struct symbol<  255>  { enum {size =  8}; };
      template <> struct symbol<  511>  { enum {size =  9}; };
      template <> struct symbol< 1023>  { enum {size = 10}; };
      template <> struct symbol< 2047>  { enum {size = 11}; };
      template <> struct symbol< 4195>  { enum {size = 12}; };
      template <> struct symbol< 8191>  { enum {size = 13}; };
      template <> struct symbol<16383>  { enum {size = 14}; };
      template <> struct symbol<32768>  { enum {size = 15}; };
      template <> struct symbol<65535>  { enum {size = 16}; };

      /* Credits: Modern C++ Design - Andrei Alexandrescu */
      template <bool> class __static_assert__
      {
      public:

          __static_assert__(...) {}
      };

      template <> class __static_assert__<true> {};
      template <> class __static_assert__<false>;

      template <std::size_t code_length, std::size_t fec_length, std::size_t data_length>
      struct validate_reed_solomon_code_parameters
      {
      private:

         __static_assert__<(code_length > 0)> assertion1;
         __static_assert__<(code_length > fec_length)> assertion2;
         __static_assert__<(code_length > data_length)> assertion3;
         __static_assert__<(code_length == fec_length + data_length)> assertion4;
      };

      template <std::size_t code_length, std::size_t fec_length, std::size_t data_length>
      struct validate_reed_solomon_block_parameters
      {
      private:

         __static_assert__<(code_length > 0)> assertion1;
         __static_assert__<(code_length > fec_length)> assertion2;
         __static_assert__<(code_length > data_length)> assertion3;
         __static_assert__<(code_length == fec_length + data_length)> assertion4;
      };

      template <typename Encoder, typename Decoder>
      struct equivalent_encoder_decoder
      {
      private:

         __static_assert__<(Encoder::trait::code_length == Decoder::trait::code_length)> assertion1;
         __static_assert__<(Encoder::trait::fec_length  == Decoder::trait::fec_length) > assertion2;
         __static_assert__<(Encoder::trait::data_length == Decoder::trait::data_length)> assertion3;
      };

      template <std::size_t code_length_, std::size_t fec_length_, std::size_t data_length_ = code_length_ - fec_length_>
      class reed_solomon_triat
      {
      public:

         typedef validate_reed_solomon_code_parameters<code_length_,fec_length_,data_length_> vrscp;

         enum { code_length = code_length_ };
         enum { fec_length  = fec_length_  };
         enum { data_length = data_length_ };
      };

   }

} // namespace schifra

#endif
