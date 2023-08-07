#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _asin_reg(void);
extern void _cldif_CA3_NKCC1_HCO3_reg(void);
extern void _gabaA_Cl_HCO3_reg(void);
extern void _tonic_reg(void);
extern void _vecevent_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"asin.mod\"");
    fprintf(stderr, " \"cldif_CA3_NKCC1_HCO3.mod\"");
    fprintf(stderr, " \"gabaA_Cl_HCO3.mod\"");
    fprintf(stderr, " \"tonic.mod\"");
    fprintf(stderr, " \"vecevent.mod\"");
    fprintf(stderr, "\n");
  }
  _asin_reg();
  _cldif_CA3_NKCC1_HCO3_reg();
  _gabaA_Cl_HCO3_reg();
  _tonic_reg();
  _vecevent_reg();
}

#if defined(__cplusplus)
}
#endif
