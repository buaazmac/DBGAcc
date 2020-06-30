/// @file xed-operand-ctype-map.c

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
#include "xed-operand-ctype-map.h"
static xed_operand_ctype_enum_t xed_operand_ctype[XED_OPERAND_LAST];
static unsigned int  xed_operand_bits[XED_OPERAND_LAST];
void xed_init_operand_ctypes(void)
{
   xed_operand_ctype[XED_OPERAND_AGEN]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_LAST_F2F3]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_DISP_WIDTH]=XED_OPERAND_CTYPE_XED_UINT8_T;
   xed_operand_ctype[XED_OPERAND_USING_DEFAULT_SEGMENT0]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_USING_DEFAULT_SEGMENT1]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_HINT]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_SAE]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_MODE_FIRST_PREFIX]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_UIMM1]=XED_OPERAND_CTYPE_XED_UINT8_T;
   xed_operand_ctype[XED_OPERAND_UIMM0]=XED_OPERAND_CTYPE_XED_UINT64_T;
   xed_operand_ctype[XED_OPERAND_SMODE]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_RM]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_REP]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_AMD3DNOW]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_MAP]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_VEXVALID]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_SIBINDEX]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_SIB]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_NOMINAL_OPCODE]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_SEG1]=XED_OPERAND_CTYPE_XED_REG_ENUM_T;
   xed_operand_ctype[XED_OPERAND_SEG0]=XED_OPERAND_CTYPE_XED_REG_ENUM_T;
   xed_operand_ctype[XED_OPERAND_PTR]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_INDEX]=XED_OPERAND_CTYPE_XED_REG_ENUM_T;
   xed_operand_ctype[XED_OPERAND_ILD_F2]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_SCALE]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_ESRC]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_NREXES]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_NEED_MEMDISP]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_POS_SIB]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_NPREFIXES]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_HAS_SIB]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_EOSZ]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_ELEMENT_SIZE]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_POS_DISP]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_UBIT]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_VEXDEST210]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_VEXDEST3]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_CET]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_CLDEMOTE]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_P4]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_MODEP55C]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_ICLASS]=XED_OPERAND_CTYPE_XED_ICLASS_ENUM_T;
   xed_operand_ctype[XED_OPERAND_IMM_WIDTH]=XED_OPERAND_CTYPE_XED_UINT8_T;
   xed_operand_ctype[XED_OPERAND_BCRC]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_ERROR]=XED_OPERAND_CTYPE_XED_ERROR_ENUM_T;
   xed_operand_ctype[XED_OPERAND_NELEM]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_IMM0SIGNED]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_REG8]=XED_OPERAND_CTYPE_XED_REG_ENUM_T;
   xed_operand_ctype[XED_OPERAND_REG6]=XED_OPERAND_CTYPE_XED_REG_ENUM_T;
   xed_operand_ctype[XED_OPERAND_REG7]=XED_OPERAND_CTYPE_XED_REG_ENUM_T;
   xed_operand_ctype[XED_OPERAND_REG4]=XED_OPERAND_CTYPE_XED_REG_ENUM_T;
   xed_operand_ctype[XED_OPERAND_REG5]=XED_OPERAND_CTYPE_XED_REG_ENUM_T;
   xed_operand_ctype[XED_OPERAND_REG2]=XED_OPERAND_CTYPE_XED_REG_ENUM_T;
   xed_operand_ctype[XED_OPERAND_REG3]=XED_OPERAND_CTYPE_XED_REG_ENUM_T;
   xed_operand_ctype[XED_OPERAND_REG0]=XED_OPERAND_CTYPE_XED_REG_ENUM_T;
   xed_operand_ctype[XED_OPERAND_REG1]=XED_OPERAND_CTYPE_XED_REG_ENUM_T;
   xed_operand_ctype[XED_OPERAND_BASE0]=XED_OPERAND_CTYPE_XED_REG_ENUM_T;
   xed_operand_ctype[XED_OPERAND_BASE1]=XED_OPERAND_CTYPE_XED_REG_ENUM_T;
   xed_operand_ctype[XED_OPERAND_MOD]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_SEG_OVD]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_REX]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_REXB]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_POS_IMM]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_CHIP]=XED_OPERAND_CTYPE_XED_CHIP_ENUM_T;
   xed_operand_ctype[XED_OPERAND_REXW]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_ILD_F3]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_REXR]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_ENCODER_PREFERRED]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_REG]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_REXX]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_PREFIX66]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_REXRR]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_ASZ]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_MASK]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_MEM1]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_EASZ]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_POS_IMM1]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_MEM_WIDTH]=XED_OPERAND_CTYPE_XED_UINT16_T;
   xed_operand_ctype[XED_OPERAND_LZCNT]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_MEM0]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_BRDISP_WIDTH]=XED_OPERAND_CTYPE_XED_UINT8_T;
   xed_operand_ctype[XED_OPERAND_IMM1_BYTES]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_TZCNT]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_DF64]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_LOCK]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_HAS_MODRM]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_ZEROING]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_SRM]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_LLRC]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_NEEDREX]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_SKIP_OSZ]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_OUTREG]=XED_OPERAND_CTYPE_XED_REG_ENUM_T;
   xed_operand_ctype[XED_OPERAND_DEFAULT_SEG]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_NSEG_PREFIXES]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_VEX_C4]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_POS_MODRM]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_BCAST]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_TYPE]=XED_OPERAND_CTYPE_XED_OPERAND_ELEMENT_TYPE_ENUM_T;
   xed_operand_ctype[XED_OPERAND_DISP]=XED_OPERAND_CTYPE_XED_INT64_T;
   xed_operand_ctype[XED_OPERAND_VEX_PREFIX]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_DUMMY]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_NOREX]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_ROUNDC]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_SIBBASE]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_POS_NOMINAL_OPCODE]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_OUT_OF_BYTES]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_IMM1]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_IMM0]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_NO_SCALE_DISP8]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_RELBR]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_ILD_SEG]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_DF32]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_REALMODE]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_MODRM_BYTE]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_MODE]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_MPXMODE]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_SIBSCALE]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_VL]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_OSZ]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_WBNOINVD]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_MAX_BYTES]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_FIRST_F2F3]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_MODEP5]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_ctype[XED_OPERAND_VEXDEST4]=XED_OPERAND_CTYPE_XED_BITS_T;
   xed_operand_bits[XED_OPERAND_AGEN]=1;
   xed_operand_bits[XED_OPERAND_LAST_F2F3]=2;
   xed_operand_bits[XED_OPERAND_DISP_WIDTH]=8;
   xed_operand_bits[XED_OPERAND_USING_DEFAULT_SEGMENT0]=1;
   xed_operand_bits[XED_OPERAND_USING_DEFAULT_SEGMENT1]=1;
   xed_operand_bits[XED_OPERAND_HINT]=3;
   xed_operand_bits[XED_OPERAND_SAE]=1;
   xed_operand_bits[XED_OPERAND_MODE_FIRST_PREFIX]=1;
   xed_operand_bits[XED_OPERAND_UIMM1]=8;
   xed_operand_bits[XED_OPERAND_UIMM0]=64;
   xed_operand_bits[XED_OPERAND_SMODE]=2;
   xed_operand_bits[XED_OPERAND_RM]=3;
   xed_operand_bits[XED_OPERAND_REP]=2;
   xed_operand_bits[XED_OPERAND_AMD3DNOW]=1;
   xed_operand_bits[XED_OPERAND_MAP]=4;
   xed_operand_bits[XED_OPERAND_VEXVALID]=3;
   xed_operand_bits[XED_OPERAND_SIBINDEX]=3;
   xed_operand_bits[XED_OPERAND_SIB]=1;
   xed_operand_bits[XED_OPERAND_NOMINAL_OPCODE]=8;
   xed_operand_bits[XED_OPERAND_SEG1]=16;
   xed_operand_bits[XED_OPERAND_SEG0]=16;
   xed_operand_bits[XED_OPERAND_PTR]=1;
   xed_operand_bits[XED_OPERAND_INDEX]=16;
   xed_operand_bits[XED_OPERAND_ILD_F2]=1;
   xed_operand_bits[XED_OPERAND_SCALE]=4;
   xed_operand_bits[XED_OPERAND_ESRC]=4;
   xed_operand_bits[XED_OPERAND_NREXES]=8;
   xed_operand_bits[XED_OPERAND_NEED_MEMDISP]=6;
   xed_operand_bits[XED_OPERAND_POS_SIB]=8;
   xed_operand_bits[XED_OPERAND_NPREFIXES]=8;
   xed_operand_bits[XED_OPERAND_HAS_SIB]=1;
   xed_operand_bits[XED_OPERAND_EOSZ]=2;
   xed_operand_bits[XED_OPERAND_ELEMENT_SIZE]=9;
   xed_operand_bits[XED_OPERAND_POS_DISP]=8;
   xed_operand_bits[XED_OPERAND_UBIT]=1;
   xed_operand_bits[XED_OPERAND_VEXDEST210]=3;
   xed_operand_bits[XED_OPERAND_VEXDEST3]=1;
   xed_operand_bits[XED_OPERAND_CET]=1;
   xed_operand_bits[XED_OPERAND_CLDEMOTE]=1;
   xed_operand_bits[XED_OPERAND_P4]=1;
   xed_operand_bits[XED_OPERAND_MODEP55C]=1;
   xed_operand_bits[XED_OPERAND_ICLASS]=16;
   xed_operand_bits[XED_OPERAND_IMM_WIDTH]=8;
   xed_operand_bits[XED_OPERAND_BCRC]=1;
   xed_operand_bits[XED_OPERAND_ERROR]=8;
   xed_operand_bits[XED_OPERAND_NELEM]=4;
   xed_operand_bits[XED_OPERAND_IMM0SIGNED]=1;
   xed_operand_bits[XED_OPERAND_REG8]=16;
   xed_operand_bits[XED_OPERAND_REG6]=16;
   xed_operand_bits[XED_OPERAND_REG7]=16;
   xed_operand_bits[XED_OPERAND_REG4]=16;
   xed_operand_bits[XED_OPERAND_REG5]=16;
   xed_operand_bits[XED_OPERAND_REG2]=16;
   xed_operand_bits[XED_OPERAND_REG3]=16;
   xed_operand_bits[XED_OPERAND_REG0]=16;
   xed_operand_bits[XED_OPERAND_REG1]=16;
   xed_operand_bits[XED_OPERAND_BASE0]=16;
   xed_operand_bits[XED_OPERAND_BASE1]=16;
   xed_operand_bits[XED_OPERAND_MOD]=2;
   xed_operand_bits[XED_OPERAND_SEG_OVD]=3;
   xed_operand_bits[XED_OPERAND_REX]=1;
   xed_operand_bits[XED_OPERAND_REXB]=1;
   xed_operand_bits[XED_OPERAND_POS_IMM]=8;
   xed_operand_bits[XED_OPERAND_CHIP]=16;
   xed_operand_bits[XED_OPERAND_REXW]=1;
   xed_operand_bits[XED_OPERAND_ILD_F3]=1;
   xed_operand_bits[XED_OPERAND_REXR]=1;
   xed_operand_bits[XED_OPERAND_ENCODER_PREFERRED]=1;
   xed_operand_bits[XED_OPERAND_REG]=3;
   xed_operand_bits[XED_OPERAND_REXX]=1;
   xed_operand_bits[XED_OPERAND_PREFIX66]=1;
   xed_operand_bits[XED_OPERAND_REXRR]=1;
   xed_operand_bits[XED_OPERAND_ASZ]=1;
   xed_operand_bits[XED_OPERAND_MASK]=3;
   xed_operand_bits[XED_OPERAND_MEM1]=1;
   xed_operand_bits[XED_OPERAND_EASZ]=2;
   xed_operand_bits[XED_OPERAND_POS_IMM1]=8;
   xed_operand_bits[XED_OPERAND_MEM_WIDTH]=16;
   xed_operand_bits[XED_OPERAND_LZCNT]=1;
   xed_operand_bits[XED_OPERAND_MEM0]=1;
   xed_operand_bits[XED_OPERAND_BRDISP_WIDTH]=8;
   xed_operand_bits[XED_OPERAND_IMM1_BYTES]=8;
   xed_operand_bits[XED_OPERAND_TZCNT]=1;
   xed_operand_bits[XED_OPERAND_DF64]=1;
   xed_operand_bits[XED_OPERAND_LOCK]=1;
   xed_operand_bits[XED_OPERAND_HAS_MODRM]=2;
   xed_operand_bits[XED_OPERAND_ZEROING]=1;
   xed_operand_bits[XED_OPERAND_SRM]=3;
   xed_operand_bits[XED_OPERAND_LLRC]=2;
   xed_operand_bits[XED_OPERAND_NEEDREX]=1;
   xed_operand_bits[XED_OPERAND_SKIP_OSZ]=1;
   xed_operand_bits[XED_OPERAND_OUTREG]=16;
   xed_operand_bits[XED_OPERAND_DEFAULT_SEG]=2;
   xed_operand_bits[XED_OPERAND_NSEG_PREFIXES]=8;
   xed_operand_bits[XED_OPERAND_VEX_C4]=1;
   xed_operand_bits[XED_OPERAND_POS_MODRM]=8;
   xed_operand_bits[XED_OPERAND_BCAST]=5;
   xed_operand_bits[XED_OPERAND_TYPE]=4;
   xed_operand_bits[XED_OPERAND_DISP]=64;
   xed_operand_bits[XED_OPERAND_VEX_PREFIX]=2;
   xed_operand_bits[XED_OPERAND_DUMMY]=1;
   xed_operand_bits[XED_OPERAND_NOREX]=1;
   xed_operand_bits[XED_OPERAND_ROUNDC]=3;
   xed_operand_bits[XED_OPERAND_SIBBASE]=3;
   xed_operand_bits[XED_OPERAND_POS_NOMINAL_OPCODE]=8;
   xed_operand_bits[XED_OPERAND_OUT_OF_BYTES]=1;
   xed_operand_bits[XED_OPERAND_IMM1]=1;
   xed_operand_bits[XED_OPERAND_IMM0]=1;
   xed_operand_bits[XED_OPERAND_NO_SCALE_DISP8]=1;
   xed_operand_bits[XED_OPERAND_RELBR]=1;
   xed_operand_bits[XED_OPERAND_ILD_SEG]=8;
   xed_operand_bits[XED_OPERAND_DF32]=1;
   xed_operand_bits[XED_OPERAND_REALMODE]=1;
   xed_operand_bits[XED_OPERAND_MODRM_BYTE]=8;
   xed_operand_bits[XED_OPERAND_MODE]=2;
   xed_operand_bits[XED_OPERAND_MPXMODE]=1;
   xed_operand_bits[XED_OPERAND_SIBSCALE]=2;
   xed_operand_bits[XED_OPERAND_VL]=2;
   xed_operand_bits[XED_OPERAND_OSZ]=1;
   xed_operand_bits[XED_OPERAND_WBNOINVD]=1;
   xed_operand_bits[XED_OPERAND_MAX_BYTES]=8;
   xed_operand_bits[XED_OPERAND_FIRST_F2F3]=2;
   xed_operand_bits[XED_OPERAND_MODEP5]=1;
   xed_operand_bits[XED_OPERAND_VEXDEST4]=1;
}
xed_operand_ctype_enum_t xed_operand_get_ctype(xed_operand_enum_t opname)
{
    xed_assert(opname <XED_OPERAND_LAST);
    return xed_operand_ctype[opname];
}
unsigned int xed_operand_decider_get_width(xed_operand_enum_t opname)
{
    xed_assert(opname <XED_OPERAND_LAST);
    return xed_operand_bits[opname];
}
