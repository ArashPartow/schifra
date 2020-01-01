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


/*
   Description: This example will demonstrate how convert an array
                of chars to a 16-bit per symbol reed solomon block.
*/


#include <cstddef>
#include <iostream>
#include <string>

#include "schifra_reed_solomon_block.hpp"
#include "schifra_reed_solomon_bitio.hpp"


int main()
{
   /* Finite Field Parameter */
   const std::size_t field_descriptor = 16;

   /* Reed Solomon Code Parameters */
   const std::size_t code_length = 65535; //(2^16 - 1)
   const std::size_t fec_length  = 100;
   const std::size_t data_length = code_length - fec_length;

   std::string message(data_length << 1,static_cast<unsigned char>(0xAA));
   schifra::reed_solomon::block<code_length,fec_length> rs_block;
   schifra::reed_solomon::bitio::convert_data_to_symbol<field_descriptor>(message.c_str(),
                                                                          message.size(),
                                                                          rs_block.data);

   return 0;
}
