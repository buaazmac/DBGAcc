/// @file include-private/xed3-dynamic-part1-capture.h

// This file was automatically generated.
// Do not edit this file.

#if !defined(INCLUDE_PRIVATE_XED3_DYNAMIC_PART1_CAPTURE_H)
# define INCLUDE_PRIVATE_XED3_DYNAMIC_PART1_CAPTURE_H
/*BEGIN_LEGAL

Copyright (c) 2018 Intel Corporation

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
  
END_LEGAL */
#include "xed-internal-header.h"
#include "xed3-nt-capture.h"
static XED_INLINE xed_error_enum_t xed3_dynamic_decode_part1(xed_decoded_inst_t* d);

static XED_INLINE xed_error_enum_t xed3_dynamic_decode_part1(xed_decoded_inst_t* d)
{
xed3_capture_nt_OSZ_NONTERM(d);
if (xed3_operand_get_error(d)) {
return xed3_operand_get_error(d);
}
xed3_capture_nt_ASZ_NONTERM(d);
if (xed3_operand_get_error(d)) {
return xed3_operand_get_error(d);
}
return XED_ERROR_NONE;
}
#endif
