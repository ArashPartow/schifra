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


#ifndef INCLUDE_SCHIFRA_REED_SOLOMON_BITIO_HPP
#define INCLUDE_SCHIFRA_REED_SOLOMON_BITIO_HPP


#include <iostream>


namespace schifra
{

   namespace reed_solomon
   {

      namespace bitio
      {

         template <std::size_t symbol_bit_count> class convert_data_to_symbol;

         template <>
         class convert_data_to_symbol<2>
         {
         public:

            template <typename BitBlock>
            convert_data_to_symbol(const BitBlock data[], const std::size_t data_length, int symbol[])
            {
               const BitBlock* d_it = &  data[0];
               int*            s_it = &symbol[0];

               for (std::size_t i = 0; i < data_length; ++i, ++d_it, s_it+=4)
               {
                  (* s_it     ) =  (*d_it)       & 0x03;
                  (*(s_it + 1)) = ((*d_it) >> 2) & 0x03;
                  (*(s_it + 2)) = ((*d_it) >> 4) & 0x03;
                  (*(s_it + 3)) = ((*d_it) >> 6) & 0x03;
               }
            }
         };

         template <>
         class convert_data_to_symbol<4>
         {
         public:

            template <typename BitBlock>
            convert_data_to_symbol(const BitBlock data[], const std::size_t data_length, int symbol[])
            {
               const BitBlock* d_it = &  data[0];
               int*            s_it = &symbol[0];

               for (std::size_t i = 0; i < data_length; ++i, ++d_it, s_it+=2)
               {
                  (* s_it     ) =  (*d_it)       & 0x0F;
                  (*(s_it + 1)) = ((*d_it) >> 4) & 0x0F;
               }
            }
         };

         template <>
         class convert_data_to_symbol<8>
         {
         public:

            template <typename BitBlock>
            convert_data_to_symbol(const BitBlock data[], const std::size_t data_length, int symbol[])
            {
               const BitBlock* d_it = &  data[0];
               int*            s_it = &symbol[0];

               for (std::size_t i = 0; i < data_length; ++i, ++d_it, ++s_it)
               {
                  (*s_it) =  (*d_it) & 0xFF;
               }
            }
         };

         template <>
         class convert_data_to_symbol<16>
         {
         public:

            template <typename BitBlock>
            convert_data_to_symbol(const BitBlock data[], const std::size_t data_length, int symbol[])
            {
               const BitBlock* d_it = &  data[0];
               int*            s_it = &symbol[0];

               for (std::size_t i = 0; i < data_length; i+=2, d_it+=2, ++s_it)
               {
                  (*s_it)  = (*d_it) & 0x000000FF;
                  (*s_it) |= (static_cast<int>((*(d_it + 1))) << 8) & 0x0000FF00;
               }
            }
         };

         template <>
         class convert_data_to_symbol<24>
         {
         public:

            template <typename BitBlock>
            convert_data_to_symbol(const BitBlock data[], const std::size_t data_length, int symbol[])
            {
               BitBlock* d_it = &  data[0];
               int*      s_it = &symbol[0];

               for (std::size_t i = 0; i < data_length; i+=3, d_it+=3, ++s_it)
               {
                  (*s_it) |= (*d_it) & 0x000000FF;
                  (*s_it) |= (static_cast<int>((*(d_it + 1))) <<  8) & 0x0000FF00;
                  (*s_it) |= (static_cast<int>((*(d_it + 2))) << 16) & 0x00FF0000;
               }
            }
         };

         template <std::size_t symbol_bit_count> class convert_symbol_to_data;

         template <>
         class convert_symbol_to_data<4>
         {
         public:

            template <typename BitBlock>
            convert_symbol_to_data(const int symbol[], BitBlock data[], const std::size_t data_length)
            {
               BitBlock*  d_it = &  data[0];
               const int* s_it = &symbol[0];

               for (std::size_t i = 0; i < data_length; ++i, ++d_it, ++s_it)
               {
                  (*d_it)  =        (*s_it) & 0x0000000F;
                  (*d_it) |= ((*(s_it + 1)) & 0x0000000F) << 4;
               }
            }
         };

         template <>
         class convert_symbol_to_data<8>
         {
         public:
            template <typename BitBlock>
            convert_symbol_to_data(const int symbol[], BitBlock data[], const std::size_t data_length)
            {
               BitBlock*  d_it = &  data[0];
               const int* s_it = &symbol[0];

               for (std::size_t i = 0; i < data_length; ++i, ++d_it, ++s_it)
               {
                  (*d_it) = static_cast<BitBlock>((*s_it) & 0xFF);
               }
            }
         };

         template <>
         class convert_symbol_to_data<16>
         {
         public:

            template <typename BitBlock>
            convert_symbol_to_data(const int symbol[], BitBlock data[], const std::size_t data_length)
            {
               BitBlock*  d_it = &  data[0];
               const int* s_it = &symbol[0];

               for (std::size_t i = 0; i < data_length; ++i, ++d_it, ++s_it)
               {
                  (*d_it) = (*s_it) & 0xFFFF;
               }
            }
         };

      } // namespace bitio

   } // namespace reed_solomon

} // namespace schifra


#endif
