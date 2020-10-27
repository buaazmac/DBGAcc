#ifndef __SIM_API
#define __SIM_API


#define SIM_CMD_ROI_TOGGLE      0  // Deprecated, for compatibility with programs compiled long ago
#define SIM_CMD_ROI_START       1
#define SIM_CMD_ROI_END         2
#define SIM_CMD_MHZ_SET         3
#define SIM_CMD_MARKER          4
#define SIM_CMD_USER            5
#define SIM_CMD_INSTRUMENT_MODE 6
#define SIM_CMD_MHZ_GET         7
#define SIM_CMD_IN_SIMULATOR    8
#define SIM_CMD_PROC_ID         9
#define SIM_CMD_THREAD_ID       10
#define SIM_CMD_NUM_PROCS       11
#define SIM_CMD_NUM_THREADS     12
#define SIM_CMD_NAMED_MARKER    13
#define SIM_CMD_SET_THREAD_NAME 14

#define SIM_PIM_OFFLOAD_START 15 // [ZMAC] magic number for starting PIM offloading phase
#define SIM_PIM_OFFLOAD_END 16 // [ZMAC] magic number for ending PIM offloading phase

#define SIM_NDP_DATA_RANGE_GLOBAL 17 // [ZMAC] magic number of indicating global NDP data range
#define SIM_NDP_DATA_RANGE_LOCAL 18 // [ZMAC] magic number of indicating local NDP data range

#define SIM_OPT_INSTRUMENT_DETAILED    0
#define SIM_OPT_INSTRUMENT_WARMUP      1
#define SIM_OPT_INSTRUMENT_FASTFORWARD 2

#if defined(ARM_64)

#define SimMagic0(cmd) ({                       \
   unsigned long _cmd = (cmd), _res;            \
   asm volatile (           \
   "mov x1, %[x]\n"         \
   "\tbfm x0, x0, 0, 0\n"   \
   : [ret]"=r"(_res)        \
   : [x]"r"(_cmd)           \
   );                       \
})

#define SimMagic1(cmd, arg0) ({              \
   unsigned long _cmd = (cmd), _arg0 = (arg0), _res; \
   asm volatile (           \
   "mov x1, %[x]\n"         \
   "\tmov x2, %[y]\n"       \
   "\tbfm x0, x0, 0, 0\n"   \
   : [ret]"=r"(_res)        \
   : [x]"r"(_cmd),          \
     [y]"r"(_arg0)          \
   : "x2", "x1"                   \
   );                       \
})

#define SimMagic2(cmd, arg0, arg1) ({        \
   unsigned long _cmd = (cmd), _arg0 = (arg0), _arg1 = (arg1), _res; \
   asm volatile (           \
   "mov x1, %[x]\n"         \
   "\tmov x2, %[y]\n"       \
   "\tmov x3, %[z]\n"       \
   "\tbfm x0, x0, 0, 0\n"   \
   : [ret]"=r"(_res)        \
   : [x]"r"(_cmd),          \
     [y]"r"(_arg0),          \
     [z]"r"(_arg1)          \
   : "x1", "x2", "x3"                   \
   );                       \
})

#else  // end ARM_64

#if defined(__i386)
   #define MAGIC_REG_A "eax"
   #define MAGIC_REG_B "edx" // Required for -fPIC support
   #define MAGIC_REG_C "ecx"
#else
   #define MAGIC_REG_A "rax"
   #define MAGIC_REG_B "rbx"
   #define MAGIC_REG_C "rcx"
#endif

#define SimMagic0(cmd) ({                    \
   unsigned long _cmd = (cmd), _res;         \
   __asm__ __volatile__ (                    \
   "mov %1, %%" MAGIC_REG_A "\n"             \
   "\txchg %%bx, %%bx\n"                     \
   : "=a" (_res)           /* output    */   \
   : "g"(_cmd)             /* input     */   \
      );                   /* clobbered */   \
   _res;                                     \
})

#define SimMagic1(cmd, arg0) ({              \
   unsigned long _cmd = (cmd), _arg0 = (arg0), _res; \
   __asm__ __volatile__ (                    \
   "mov %1, %%" MAGIC_REG_A "\n"             \
   "\tmov %2, %%" MAGIC_REG_B "\n"           \
   "\txchg %%bx, %%bx\n"                     \
   : "=a" (_res)           /* output    */   \
   : "g"(_cmd),                              \
     "g"(_arg0)            /* input     */   \
   : "%" MAGIC_REG_B );    /* clobbered */   \
   _res;                                     \
})

#define SimMagic2(cmd, arg0, arg1) ({        \
   unsigned long _cmd = (cmd), _arg0 = (arg0), _arg1 = (arg1), _res; \
   __asm__ __volatile__ (                    \
   "mov %1, %%" MAGIC_REG_A "\n"             \
   "\tmov %2, %%" MAGIC_REG_B "\n"           \
   "\tmov %3, %%" MAGIC_REG_C "\n"           \
   "\txchg %%bx, %%bx\n"                     \
   : "=a" (_res)           /* output    */   \
   : "g"(_cmd),                              \
     "g"(_arg0),                             \
     "g"(_arg1)            /* input     */   \
   : "%" MAGIC_REG_B, "%" MAGIC_REG_C ); /* clobbered */ \
   _res;                                     \
})

#endif

#define SimRoiStart()             SimMagic0(SIM_CMD_ROI_START)
#define SimRoiEnd()               SimMagic0(SIM_CMD_ROI_END)
#define SimGetProcId()            SimMagic0(SIM_CMD_PROC_ID)
#define SimGetThreadId()          SimMagic0(SIM_CMD_THREAD_ID)
#define SimSetThreadName(name)    SimMagic1(SIM_CMD_SET_THREAD_NAME, (unsigned long)(name))
#define SimGetNumProcs()          SimMagic0(SIM_CMD_NUM_PROCS)
#define SimGetNumThreads()        SimMagic0(SIM_CMD_NUM_THREADS)
#define SimSetFreqMHz(proc, mhz)  SimMagic2(SIM_CMD_MHZ_SET, proc, mhz)
#define SimSetOwnFreqMHz(mhz)     SimSetFreqMHz(SimGetProcId(), mhz)
#define SimGetFreqMHz(proc)       SimMagic1(SIM_CMD_MHZ_GET, proc)
#define SimGetOwnFreqMHz()        SimGetFreqMHz(SimGetProcId())
#define SimMarker(arg0, arg1)     SimMagic2(SIM_CMD_MARKER, arg0, arg1)
#define SimNamedMarker(arg0, str) SimMagic2(SIM_CMD_NAMED_MARKER, arg0, (unsigned long)(str))
#define SimUser(cmd, arg)         SimMagic2(SIM_CMD_USER, cmd, arg)
#define SimSetInstrumentMode(opt) SimMagic1(SIM_CMD_INSTRUMENT_MODE, opt)
#define SimInSimulator()          (SimMagic0(SIM_CMD_IN_SIMULATOR)!=SIM_CMD_IN_SIMULATOR)

#define SimPimOffloadStart()      SimMagic0(SIM_PIM_OFFLOAD_START) // [ZMAC] interface for starting PIM offloading
#define SimPimOffloadEnd()        SimMagic0(SIM_PIM_OFFLOAD_END) // [ZMAC] interface for ending PIM offloading

#define SimNdpRangeGlobal(arg0, arg1)   SimMagic2(SIM_NDP_DATA_RANGE_GLOBAL, arg0, arg1) // [NDP] indicate a data range as NDP data that distributed globally 
#define SimNdpRangeLocal(arg0, arg1)    SimMagic2(SIM_NDP_DATA_RANGE_LOCAL, arg0, arg1) // [NDP] indicate a data range as NDP data to a specific core

#endif /* __SIM_API */
