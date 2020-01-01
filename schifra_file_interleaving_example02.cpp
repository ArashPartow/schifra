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


#include <cstddef>
#include <iostream>
#include <string>

#include "schifra_reed_solomon_file_interleaver.hpp"


int main()
{
   const std::size_t code_length = 255;
   const std::size_t stack_size  = 255;

   const std::string input_file_name  = "output.inter";
   const std::string output_file_name = "output.deinter";

   schifra::reed_solomon::file_deinterleaver<code_length,stack_size>(input_file_name,output_file_name);

   return 0;
}
