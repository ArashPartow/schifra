#### Description
Schifra is a very robust, highly optimized and extremely configurable
Reed-Solomon error correcting code library for software applications
implemented in C++. Schifra supports standard, shortened and punctured
Reed-Solomon codes. It also has support for stacked product codes and
interleaving.

Schifra provides a concise, predictable and deterministic interface
which lends itself to easy and seamless integration into the
development of complex data communication projects requiring Reed-
Solomon error correcting code capabilities.

----

#### Download
http://www.schifra.com

----

#### Features
* Errors and Erasures
* Supported Symbol Sizes - 2 to 32 bits
* Variable Code Block Length
* User defined primitive polynomial and finite field
* Accurate and Validated Reed-Solomon Codecs - Complete combinatorial errors and erasures unit testing
* Supported Architectures For Optimizations - x86-32, x86-64, PowerPC, m68k, XScale
* Supported Reed-Solomon Codes - Intelsat 1-4, DVB(S and T), MPEG-2 TSP, VDL Mode 2-4, MIL-STD-188-165a, ITU-T G.709, CCSDS (Basis transform), CIRC, ETS 300-421, ETS 300-429, xDSL, PostBar, MaxiCode, DataMatrix and many more...
* Shortened, Punctured and Concatenated Reed-Solomon Codes - WiMAX IEEE 802.16d standard
* Product Codes
* Standard and Algebraic Interleavers
* Special Optimized Decoder - For cases of 2t = 2, 4, 6, 16 and 32
* DO-178B Level A Certified Reed-Solomon Codec - RTCA DO-224A for VDL mode 2 and 3, RTCA DO-242A (ADS-B), Certified for installation on-board airborne systems
* Reed-Solomon Based channel code for Erasure Channels

----

#### Compatible C++ Compilers
+ GNU Compiler Collection (3.1+)
+ Intel® C++ Compiler (8.x+)
+ Clang/LLVM (1.1+)
+ PGI C++ (10.x+)
+ Microsoft Visual Studio C++ Compiler (7.1+)

----

#### C++ Encoder/Decoder Example
This example will demonstrate using C++ how to instantiate a Reed-
Solomon encoder and decoder, add the full amount of possible errors,
correct the errors, and output the various pieces of relevant
information. The Reed-Solomon code's properties are as follows:


+ Symbol size: 8-bits
+ Codeword length: 255
+ Number of data symbols: 223
+ Number of FEC symbols: 32
+ Finite Field: GF(2<sup>8</sup>)
+ Finite Field polynomial: 1x<sup>8</sup> + 1x<sup>7</sup> + 0x<sup>6</sup> + 0x<sup>5</sup> + 0x<sup>4</sup> + 0x<sup>3</sup> + 1x<sup>2</sup> + 1x<sup>1</sup> + 1x<sup>0</sup>
+ Generator polynomial roots: 32
+ Generator polynomial field index: 120th element (32 consecutive roots)

```c++
#include <cstddef>
#include <iostream>
#include <string>

#include "schifra_galois_field.hpp"
#include "schifra_galois_field_polynomial.hpp"
#include "schifra_sequential_root_generator_polynomial_creator.hpp"
#include "schifra_reed_solomon_encoder.hpp"
#include "schifra_reed_solomon_decoder.hpp"
#include "schifra_reed_solomon_block.hpp"
#include "schifra_error_processes.hpp"

int main()
{
   /* Finite Field Parameters */
   const std::size_t field_descriptor                =   8;
   const std::size_t generator_polynomial_index      = 120;
   const std::size_t generator_polynomial_root_count =  32;

   /* Reed Solomon Code Parameters */
   const std::size_t code_length = 255;
   const std::size_t fec_length  =  32;
   const std::size_t data_length = code_length - fec_length;

   /* Instantiate Finite Field and Generator Polynomials */
   schifra::galois::field field(field_descriptor,
                                schifra::galois::primitive_polynomial_size06,
                                schifra::galois::primitive_polynomial06);

   schifra::galois::field_polynomial generator_polynomial(field);

   if (
        !schifra::make_sequential_root_generator_polynomial(
                     field,
                     generator_polynomial_index,
                     generator_polynomial_root_count,
                     generator_polynomial)
      )
   {
      std::cout << "Error - Failed to create sequential root generator!" << std::endl;
      return 1;
   }

   /* Instantiate Encoder and Decoder (Codec) */
   typedef schifra::reed_solomon::encoder<code_length,fec_length,data_length> encoder_t;
   typedef schifra::reed_solomon::decoder<code_length,fec_length,data_length> decoder_t;

   encoder_t encoder(field,generator_polynomial);
   decoder_t decoder(field,generator_polynomial_index);

   std::string message = "A professional is a person who knows more and more about less"
                         "and less until they know everything about nothing";

   /* Pad message with nulls up until the code-word length */
   message.resize(code_length,0x00);

   /* Instantiate RS Block For Codec */
   schifra::reed_solomon::block<code_length,fec_length> block;

   /* Transform message into Reed-Solomon encoded codeword */
   if (!encoder.encode(message,block))
   {
      std::cout << "Error - Critical encoding failure!" << std::endl;
      return 1;
   }

   /* Add errors at every 3rd location starting at position zero */
   schifra::corrupt_message_all_errors00(block,0,3);

   if (!decoder.decode(block))
   {
      std::cout << "Error - Critical decoding failure!" << std::endl;
      return 1;
   }
   else if (!schifra::is_block_equivelent(block,message))
   {
      std::cout << "Error - Error correction failed!" << std::endl;
      return 1;
   }

   block.data_to_string(message);

   return 0;
}
```
