#include "error.h"

const int mainHasArgs = 0;
const int mainPreserveDelimiter = 0;
const int warnUnstable = 0;
void CreateConfigVarTable(void) {
  initConfigVarTable();
  installConfigVar("printModuleInitOrder", "bool", "Built-in", /* private = */ 0, /* deprecated = */ 0, "printModuleInitOrder is deprecated"
                   , /* unstable = */ 1, "The variable 'printModuleInitOrder' is unstable and its interface is subject to change in the future"
                   );
  installConfigVar("dataParTasksPerLocale", "int(64)", "Built-in", /* private = */ 0, /* deprecated = */ 0, "dataParTasksPerLocale is deprecated"
                   , /* unstable = */ 1, "The variable 'dataParTasksPerLocale' is unstable and its interface is subject to change in the future"
                   );
  installConfigVar("dataParIgnoreRunningTasks", "bool", "Built-in", /* private = */ 0, /* deprecated = */ 0, "dataParIgnoreRunningTasks is deprecated"
                   , /* unstable = */ 1, "The variable 'dataParIgnoreRunningTasks' is unstable and its interface is subject to change in the future"
                   );
  installConfigVar("dataParMinGranularity", "int(64)", "Built-in", /* private = */ 0, /* deprecated = */ 0, "dataParMinGranularity is deprecated"
                   , /* unstable = */ 1, "The variable 'dataParMinGranularity' is unstable and its interface is subject to change in the future"
                   );
  installConfigVar("memTrack", "bool", "Built-in", /* private = */ 0, /* deprecated = */ 0, "memTrack is deprecated"
                   , /* unstable = */ 1, "The variable 'memTrack' is unstable and its interface is subject to change in the future"
                   );
  installConfigVar("memStats", "bool", "Built-in", /* private = */ 0, /* deprecated = */ 0, "memStats is deprecated"
                   , /* unstable = */ 1, "The variable 'memStats' is unstable and its interface is subject to change in the future"
                   );
  installConfigVar("memLeaksByType", "bool", "Built-in", /* private = */ 0, /* deprecated = */ 0, "memLeaksByType is deprecated"
                   , /* unstable = */ 1, "The variable 'memLeaksByType' is unstable and its interface is subject to change in the future"
                   );
  installConfigVar("memLeaks", "bool", "Built-in", /* private = */ 0, /* deprecated = */ 0, "memLeaks is deprecated"
                   , /* unstable = */ 1, "The variable 'memLeaks' is unstable and its interface is subject to change in the future"
                   );
  installConfigVar("memMax", "uint(64)", "Built-in", /* private = */ 0, /* deprecated = */ 0, "memMax is deprecated"
                   , /* unstable = */ 1, "The variable 'memMax' is unstable and its interface is subject to change in the future"
                   );
  installConfigVar("memThreshold", "uint(64)", "Built-in", /* private = */ 0, /* deprecated = */ 0, "memThreshold is deprecated"
                   , /* unstable = */ 1, "The variable 'memThreshold' is unstable and its interface is subject to change in the future"
                   );
  installConfigVar("memLog", "string", "Built-in", /* private = */ 0, /* deprecated = */ 0, "memLog is deprecated"
                   , /* unstable = */ 1, "The variable 'memLog' is unstable and its interface is subject to change in the future"
                   );
  installConfigVar("memLeaksLog", "string", "Built-in", /* private = */ 0, /* deprecated = */ 0, "memLeaksLog is deprecated"
                   , /* unstable = */ 1, "The variable 'memLeaksLog' is unstable and its interface is subject to change in the future"
                   );
  installConfigVar("memLeaksByDesc", "string", "Built-in", /* private = */ 0, /* deprecated = */ 0, "memLeaksByDesc is deprecated"
                   , /* unstable = */ 1, "The variable 'memLeaksByDesc' is unstable and its interface is subject to change in the future"
                   );
  installConfigVar("debugGpu", "bool", "Built-in", /* private = */ 0, /* deprecated = */ 0, "debugGpu is deprecated"
                   , /* unstable = */ 1, "The variable 'debugGpu' is unstable and its interface is subject to change in the future"
                   );
  installConfigVar("gpuNoCpuModeWarning", "bool", "Built-in", /* private = */ 0, /* deprecated = */ 0, "gpuNoCpuModeWarning is deprecated"
                   , /* unstable = */ 1, "The variable 'gpuNoCpuModeWarning' is unstable and its interface is subject to change in the future"
                   );
  installConfigVar("gpuSyncWithHostAfterGpuOp", "bool", "Built-in", /* private = */ 0, /* deprecated = */ 0, "gpuSyncWithHostAfterGpuOp is deprecated"
                   , /* unstable = */ 1, "The variable 'gpuSyncWithHostAfterGpuOp' is unstable and its interface is subject to change in the future"
                   );
  installConfigVar("gpuUseStreamPerTask", "bool", "Built-in", /* private = */ 0, /* deprecated = */ 0, "gpuUseStreamPerTask is deprecated"
                   , /* unstable = */ 1, "The variable 'gpuUseStreamPerTask' is unstable and its interface is subject to change in the future"
                   );
  installConfigVar("enableGpuP2P", "bool", "Built-in", /* private = */ 0, /* deprecated = */ 0, "enableGpuP2P is deprecated"
                   , /* unstable = */ 1, "The variable 'enableGpuP2P' is unstable and its interface is subject to change in the future"
                   );
  installConfigVar("defaultHashTableResizeThreshold", "real(64)", "Built-in", /* private = */ 0, /* deprecated = */ 0, "defaultHashTableResizeThreshold is deprecated"
                   , /* unstable = */ 0, "defaultHashTableResizeThreshold is unstable"
                   );
  installConfigVar("numLocales", "int(64)", "Built-in", /* private = */ 0, /* deprecated = */ 0, "numLocales is deprecated"
                   , /* unstable = */ 0, "numLocales is unstable"
                   );
  installConfigVar("np", "int(64)", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "np is deprecated"
                   , /* unstable = */ 0, "np is unstable"
                   );
  installConfigVar("nworms", "int(64)", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "nworms is deprecated"
                   , /* unstable = */ 0, "nworms is unstable"
                   );
  installConfigVar("nsteps", "int(64)", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "nsteps is deprecated"
                   , /* unstable = */ 0, "nsteps is unstable"
                   );
  installConfigVar("fdogic", "real(64)", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "fdogic is deprecated"
                   , /* unstable = */ 0, "fdogic is unstable"
                   );
  installConfigVar("walldrive", "bool", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "walldrive is deprecated"
                   , /* unstable = */ 0, "walldrive is unstable"
                   );
  installConfigVar("fdogicwall", "real(64)", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "fdogicwall is deprecated"
                   , /* unstable = */ 0, "fdogicwall is unstable"
                   );
  installConfigVar("fdep", "real(64)", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "fdep is deprecated"
                   , /* unstable = */ 0, "fdep is unstable"
                   );
  installConfigVar("dogic_fdep", "real(64)", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "dogic_fdep is deprecated"
                   , /* unstable = */ 0, "dogic_fdep is unstable"
                   );
  installConfigVar("fdepwall", "real(64)", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "fdepwall is deprecated"
                   , /* unstable = */ 0, "fdepwall is unstable"
                   );
  installConfigVar("rwall", "int(64)", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "rwall is deprecated"
                   , /* unstable = */ 0, "rwall is unstable"
                   );
  installConfigVar("dt", "real(64)", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "dt is deprecated"
                   , /* unstable = */ 0, "dt is unstable"
                   );
  installConfigVar("kspring", "real(64)", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "kspring is deprecated"
                   , /* unstable = */ 0, "kspring is unstable"
                   );
  installConfigVar("k2spring", "real(64)", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "k2spring is deprecated"
                   , /* unstable = */ 0, "k2spring is unstable"
                   );
  installConfigVar("kbend", "real(64)", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "kbend is deprecated"
                   , /* unstable = */ 0, "kbend is unstable"
                   );
  installConfigVar("length0", "real(64)", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "length0 is deprecated"
                   , /* unstable = */ 0, "length0 is unstable"
                   );
  installConfigVar("rcut", "real(64)", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "rcut is deprecated"
                   , /* unstable = */ 0, "rcut is unstable"
                   );
  installConfigVar("save_interval", "int(64)", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "save_interval is deprecated"
                   , /* unstable = */ 0, "save_interval is unstable"
                   );
  installConfigVar("boundary", "int(64)", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "boundary is deprecated"
                   , /* unstable = */ 0, "boundary is unstable"
                   );
  installConfigVar("fluid_cpl", "bool", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "fluid_cpl is deprecated"
                   , /* unstable = */ 0, "fluid_cpl is unstable"
                   );
  installConfigVar("debug", "bool", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "debug is deprecated"
                   , /* unstable = */ 0, "debug is unstable"
                   );
  installConfigVar("thermo", "bool", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "thermo is deprecated"
                   , /* unstable = */ 0, "thermo is unstable"
                   );
  installConfigVar("thermow", "bool", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "thermow is deprecated"
                   , /* unstable = */ 0, "thermow is unstable"
                   );
  installConfigVar("kbt", "real(64)", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "kbt is deprecated"
                   , /* unstable = */ 0, "kbt is unstable"
                   );
  installConfigVar("fluid_rho", "real(64)", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "fluid_rho is deprecated"
                   , /* unstable = */ 0, "fluid_rho is unstable"
                   );
  installConfigVar("sigma", "real(64)", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "sigma is deprecated"
                   , /* unstable = */ 0, "sigma is unstable"
                   );
  installConfigVar("worm_particle_mass", "real(64)", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "worm_particle_mass is deprecated"
                   , /* unstable = */ 0, "worm_particle_mass is unstable"
                   );
  installConfigVar("L", "real(64)", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "L is deprecated"
                   , /* unstable = */ 0, "L is unstable"
                   );
  installConfigVar("restart_filename", "string", "amatter3d", /* private = */ 0, /* deprecated = */ 0, "restart_filename is deprecated"
                   , /* unstable = */ 0, "restart_filename is unstable"
                   );
}


