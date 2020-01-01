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


#include <iostream>

#include "schifra_reed_solomon_codec_validator.hpp"


int main()
{
   bool codec_validation_result = schifra::reed_solomon::codec_validation_test00() &&
                                  schifra::reed_solomon::codec_validation_test01() ;

   if (codec_validation_result)
   {
      std::cout << "Schifra Reed-Solomon Codec Successfully Validated!" << std::endl;
      return 0;
   }
   else
   {
      std::cout << "Schifra Reed-Solomon Codec Validation Failure!" << std::endl;
      return 1;
   }
}
