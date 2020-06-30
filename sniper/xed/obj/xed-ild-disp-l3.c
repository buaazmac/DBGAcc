/// @file xed-ild-disp-l3.c

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
/*Array declaration*/
xed_bits_t xed_lookup_BRDISPz_BRDISP_WIDTH[4];

/*Array initialization*/
void xed_lookup_function_init_BRDISPz_BRDISP_WIDTH(void)
{
xed_lookup_BRDISPz_BRDISP_WIDTH[1]=0x10;
xed_lookup_BRDISPz_BRDISP_WIDTH[2]=0x20;
xed_lookup_BRDISPz_BRDISP_WIDTH[3]=0x20;
}
/*Array declaration*/
xed_bits_t xed_lookup_MEMDISPv_DISP_WIDTH[4];

/*Array initialization*/
void xed_lookup_function_init_MEMDISPv_DISP_WIDTH(void)
{
xed_lookup_MEMDISPv_DISP_WIDTH[1]=0x10;
xed_lookup_MEMDISPv_DISP_WIDTH[2]=0x20;
xed_lookup_MEMDISPv_DISP_WIDTH[3]=0x40;
}
void xed_ild_disp_l3_init(void)
{
xed_lookup_function_init_BRDISPz_BRDISP_WIDTH();
xed_lookup_function_init_MEMDISPv_DISP_WIDTH();
}
