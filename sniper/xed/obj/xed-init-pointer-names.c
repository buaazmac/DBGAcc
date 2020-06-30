/// @file xed-init-pointer-names.c

// This file was automatically generated.
// Do not edit this file.

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
#include "xed-init-pointer-names.h"
#include <string.h>
const char* xed_pointer_name[XED_MAX_POINTER_NAMES];
const char* xed_pointer_name_suffix[XED_MAX_POINTER_NAMES];
void xed_init_pointer_names(void)
{
   memset((void*)xed_pointer_name,0,sizeof(const char*)*XED_MAX_POINTER_NAMES);
   xed_pointer_name[1] = "byte ";
   xed_pointer_name[2] = "word ";
   xed_pointer_name[4] = "dword ";
   xed_pointer_name[8] = "qword ";
   xed_pointer_name[16] = "xmmword ";
   xed_pointer_name[32] = "ymmword ";
   xed_pointer_name[64] = "zmmword ";
   memset((void*)xed_pointer_name_suffix,0,sizeof(const char*)*XED_MAX_POINTER_NAMES);
   xed_pointer_name_suffix[1] = "b ";
   xed_pointer_name_suffix[2] = "w ";
   xed_pointer_name_suffix[4] = "l ";
   xed_pointer_name_suffix[8] = "q ";
   xed_pointer_name_suffix[16] = "x ";
   xed_pointer_name_suffix[32] = "y ";
   xed_pointer_name_suffix[64] = "z ";
}
