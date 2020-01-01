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
   const std::size_t code_length =  10;
   const std::size_t stack_size  = 255;

   const std::string input_file_name   = "input.dat";
   const std::string output_file_name  = "output.inter";
   const std::string output_file_name2 = "output.deinter";

   schifra::reed_solomon::file_interleaver<code_length,stack_size>(input_file_name,output_file_name);
   schifra::reed_solomon::file_deinterleaver<code_length,stack_size>(output_file_name,output_file_name2);

   if (!schifra::fileio::files_identical(input_file_name,output_file_name2))
   {
      std::cout << "ERROR - Input file and deinterleaved files are not equivalent!" << std::endl;
      return 1;
   }

   return 0;
}
