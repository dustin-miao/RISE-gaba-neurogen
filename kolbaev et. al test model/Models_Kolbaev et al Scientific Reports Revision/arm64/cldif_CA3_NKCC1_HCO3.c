/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 static void _difusfunc(ldifusfunc2_t, NrnThread*);
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__cldif_CA3_NKCC1_HCO3
#define _nrn_initial _nrn_initial__cldif_CA3_NKCC1_HCO3
#define nrn_cur _nrn_cur__cldif_CA3_NKCC1_HCO3
#define _nrn_current _nrn_current__cldif_CA3_NKCC1_HCO3
#define nrn_jacob _nrn_jacob__cldif_CA3_NKCC1_HCO3
#define nrn_state _nrn_state__cldif_CA3_NKCC1_HCO3
#define _net_receive _net_receive__cldif_CA3_NKCC1_HCO3 
#define factors factors__cldif_CA3_NKCC1_HCO3 
#define state state__cldif_CA3_NKCC1_HCO3 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define cli0 _p[0]
#define cli0_columnindex 0
#define clo0 _p[1]
#define clo0_columnindex 1
#define hco3i0 _p[2]
#define hco3i0_columnindex 2
#define hco3o0 _p[3]
#define hco3o0_columnindex 3
#define ehco3_help _p[4]
#define ehco3_help_columnindex 4
#define ecl_help _p[5]
#define ecl_help_columnindex 5
#define cl (_p + 6)
#define cl_columnindex 6
#define hco3 (_p + 10)
#define hco3_columnindex 10
#define icl _p[14]
#define icl_columnindex 14
#define ihco3 _p[15]
#define ihco3_columnindex 15
#define cli _p[16]
#define cli_columnindex 16
#define hco3i _p[17]
#define hco3i_columnindex 17
#define hco3o _p[18]
#define hco3o_columnindex 18
#define ActPump _p[19]
#define ActPump_columnindex 19
#define Dcl (_p + 20)
#define Dcl_columnindex 20
#define Dhco3 (_p + 24)
#define Dhco3_columnindex 24
#define _g _p[28]
#define _g_columnindex 28
#define _ion_icl	*_ppvar[0]._pval
#define _ion_cli	*_ppvar[1]._pval
#define _style_cl	*((int*)_ppvar[2]._pvoid)
#define _ion_dicldv	*_ppvar[3]._pval
#define _ion_ihco3	*_ppvar[4]._pval
#define _ion_hco3i	*_ppvar[5]._pval
#define _style_hco3	*((int*)_ppvar[6]._pvoid)
#define _ion_dihco3dv	*_ppvar[7]._pval
#define diam	*_ppvar[8]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_factors(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_cldif_CA3_NKCC1_HCO3", _hoc_setdata,
 "factors_cldif_CA3_NKCC1_HCO3", _hoc_factors,
 0, 0
};
 /* declare global and static user variables */
#define DCl DCl_cldif_CA3_NKCC1_HCO3
 double DCl = 2;
#define cli_Start cli_Start_cldif_CA3_NKCC1_HCO3
 double cli_Start = 10;
#define hco3i_Start hco3i_Start_cldif_CA3_NKCC1_HCO3
 double hco3i_Start = 16;
#define tau_hco3 tau_hco3_cldif_CA3_NKCC1_HCO3
 double tau_hco3 = 1000;
#define tau_passive tau_passive_cldif_CA3_NKCC1_HCO3
 double tau_passive = 321000;
#define tau_NKCC1 tau_NKCC1_cldif_CA3_NKCC1_HCO3
 double tau_NKCC1 = 174000;
#define vrat vrat_cldif_CA3_NKCC1_HCO3
 double vrat[4];
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "DCl_cldif_CA3_NKCC1_HCO3", "um2/ms",
 "tau_NKCC1_cldif_CA3_NKCC1_HCO3", "ms",
 "tau_passive_cldif_CA3_NKCC1_HCO3", "ms",
 "tau_hco3_cldif_CA3_NKCC1_HCO3", "ms",
 "cli_Start_cldif_CA3_NKCC1_HCO3", "mM",
 "hco3i_Start_cldif_CA3_NKCC1_HCO3", "mM",
 "cli0_cldif_CA3_NKCC1_HCO3", "mM",
 "clo0_cldif_CA3_NKCC1_HCO3", "mM",
 "hco3i0_cldif_CA3_NKCC1_HCO3", "mM",
 "hco3o0_cldif_CA3_NKCC1_HCO3", "mM",
 "cl_cldif_CA3_NKCC1_HCO3", "mM",
 "hco3_cldif_CA3_NKCC1_HCO3", "mM",
 "ehco3_help_cldif_CA3_NKCC1_HCO3", "mV",
 "ecl_help_cldif_CA3_NKCC1_HCO3", "mV",
 0,0
};
 static double cl0 = 0;
 static double delta_t = 0.01;
 static double hco30 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "DCl_cldif_CA3_NKCC1_HCO3", &DCl_cldif_CA3_NKCC1_HCO3,
 "tau_NKCC1_cldif_CA3_NKCC1_HCO3", &tau_NKCC1_cldif_CA3_NKCC1_HCO3,
 "tau_passive_cldif_CA3_NKCC1_HCO3", &tau_passive_cldif_CA3_NKCC1_HCO3,
 "tau_hco3_cldif_CA3_NKCC1_HCO3", &tau_hco3_cldif_CA3_NKCC1_HCO3,
 "cli_Start_cldif_CA3_NKCC1_HCO3", &cli_Start_cldif_CA3_NKCC1_HCO3,
 "hco3i_Start_cldif_CA3_NKCC1_HCO3", &hco3i_Start_cldif_CA3_NKCC1_HCO3,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 "vrat_cldif_CA3_NKCC1_HCO3", vrat_cldif_CA3_NKCC1_HCO3, 4,
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[9]._i
 static void _ode_synonym(int, double**, Datum**);
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"cldif_CA3_NKCC1_HCO3",
 "cli0_cldif_CA3_NKCC1_HCO3",
 "clo0_cldif_CA3_NKCC1_HCO3",
 "hco3i0_cldif_CA3_NKCC1_HCO3",
 "hco3o0_cldif_CA3_NKCC1_HCO3",
 0,
 "ehco3_help_cldif_CA3_NKCC1_HCO3",
 "ecl_help_cldif_CA3_NKCC1_HCO3",
 0,
 "cl_cldif_CA3_NKCC1_HCO3[4]",
 "hco3_cldif_CA3_NKCC1_HCO3[4]",
 0,
 0};
 static Symbol* _morphology_sym;
 static Symbol* _cl_sym;
 static int _type_icl;
 static Symbol* _hco3_sym;
 static int _type_ihco3;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 29, _prop);
 	/*initialize range parameters*/
 	cli0 = 50;
 	clo0 = 133.5;
 	hco3i0 = 16;
 	hco3o0 = 26;
 	_prop->param = _p;
 	_prop->param_size = 29;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 10, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_morphology_sym);
 	_ppvar[8]._pval = &prop_ion->param[0]; /* diam */
 prop_ion = need_memb(_cl_sym);
  _type_icl = prop_ion->_type;
 nrn_check_conc_write(_prop, prop_ion, 1);
 nrn_promote(prop_ion, 3, 0);
 	_ppvar[0]._pval = &prop_ion->param[3]; /* icl */
 	_ppvar[1]._pval = &prop_ion->param[1]; /* cli */
 	_ppvar[2]._pvoid = (void*)(&(prop_ion->dparam[0]._i)); /* iontype for cl */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_dicldv */
 prop_ion = need_memb(_hco3_sym);
  _type_ihco3 = prop_ion->_type;
 nrn_check_conc_write(_prop, prop_ion, 1);
 nrn_promote(prop_ion, 3, 0);
 	_ppvar[4]._pval = &prop_ion->param[3]; /* ihco3 */
 	_ppvar[5]._pval = &prop_ion->param[1]; /* hco3i */
 	_ppvar[6]._pvoid = (void*)(&(prop_ion->dparam[0]._i)); /* iontype for hco3 */
 	_ppvar[7]._pval = &prop_ion->param[4]; /* _ion_dihco3dv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 "cl_cldif_CA3_NKCC1_HCO3", 1e-10,
 "hco3_cldif_CA3_NKCC1_HCO3", 1e-10,
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _cldif_CA3_NKCC1_HCO3_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("cl", -1.0);
 	ion_reg("hco3", -1.0);
 	_morphology_sym = hoc_lookup("morphology");
 	_cl_sym = hoc_lookup("cl_ion");
 	_hco3_sym = hoc_lookup("hco3_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 29, 10);
  hoc_register_dparam_semantics(_mechtype, 0, "cl_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "cl_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "#cl_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cl_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "hco3_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "hco3_ion");
  hoc_register_dparam_semantics(_mechtype, 6, "#hco3_ion");
  hoc_register_dparam_semantics(_mechtype, 7, "hco3_ion");
  hoc_register_dparam_semantics(_mechtype, 9, "cvodeieq");
  hoc_register_dparam_semantics(_mechtype, 8, "diam");
 	nrn_writes_conc(_mechtype, 0);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_synonym(_mechtype, _ode_synonym);
 	hoc_register_ldifus1(_difusfunc);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 cldif_CA3_NKCC1_HCO3 /Users/dustinmiao/Documents/school/rise/RISE-gaba-neurogen/kolbaev et. al test model/Models_Kolbaev et al Scientific Reports Revision/cldif_CA3_NKCC1_HCO3.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 
#define FARADAY _nrnunit_FARADAY[_nrnunit_use_legacy_]
static double _nrnunit_FARADAY[2] = {0x1.34c0c8b92a9b7p+3, 9.64853}; /* 9.64853321233100125 */
 
#define PI _nrnunit_PI[_nrnunit_use_legacy_]
static double _nrnunit_PI[2] = {0x1.921fb54442d18p+1, 3.14159}; /* 3.14159265358979312 */
 
#define F _nrnunit_F[_nrnunit_use_legacy_]
static double _nrnunit_F[2] = {0x1.78e555060882cp+16, 96485.3}; /* 96485.3321233100141 */
 
#define R _nrnunit_R[_nrnunit_use_legacy_]
static double _nrnunit_R[2] = {0x1.0a1013e8990bep+3, 8.3145}; /* 8.3144626181532395 */
 static double _zfactors_done ;
 static double _zfrat [ 4 ] ;
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int factors();
 extern double *_getelm();
 
#define _MATELM1(_row,_col)	*(_getelm(_row + 1, _col + 1))
 
#define _RHS1(_arg) _coef1[_arg + 1]
 static double *_coef1;
 
#define _linmat1  1
 static void* _sparseobj1;
 static void* _cvsparseobj1;
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[8], _dlist1[8]; static double *_temp1;
 static int state();
 
static int  factors (  ) {
   double _lr , _ldr2 ;
 _lr = 1.0 / 2.0 ;
   _ldr2 = _lr / ( 4.0 - 1.0 ) / 2.0 ;
   vrat [ 0 ] = 0.0 ;
   _zfrat [ 0 ] = 2.0 * _lr ;
   {int  _li ;for ( _li = 0 ; _li <= 4 - 2 ; _li ++ ) {
     vrat [ _li ] = vrat [ _li ] + PI * ( _lr - _ldr2 / 2.0 ) * 2.0 * _ldr2 ;
     _lr = _lr - _ldr2 ;
     _zfrat [ _li + 1 ] = 2.0 * PI * _lr / ( 2.0 * _ldr2 ) ;
     _lr = _lr - _ldr2 ;
     vrat [ _li + 1 ] = PI * ( _lr + _ldr2 / 2.0 ) * 2.0 * _ldr2 ;
     } }
    return 0; }
 
static void _hoc_factors(void) {
  double _r;
   _r = 1.;
 factors (  );
 hoc_retpushx(_r);
}
 
static int state ()
 {_reset=0;
 {
   double b_flux, f_flux, _term; int _i;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<8;_i++){
  	_RHS1(_i) = -_dt1*(_p[_slist1[_i]] - _p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} 
for (_i=0; _i < 4; _i++) {
  	_RHS1(_i + 0) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 0, _i + 0) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 4; _i++) {
  	_RHS1(_i + 4) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 4, _i + 4) *= ( diam * diam * vrat [ ((int) _i ) ]);  } }
 if ( cli0 >= cl [ 0 ] ) {
     ActPump = 1.0 ;
     }
   else {
     ActPump = 0.0 ;
     }
   /* COMPARTMENT _li , diam * diam * vrat [ ((int) _i ) ] {
     cl }
   */
 /* LONGITUDINAL_DIFFUSION _li , DCl * diam * diam * vrat [ ((int) _i ) ] {
     cl }
   */
 /* ~ cl [ 0 ] < < ( ( icl * PI * diam / FARADAY ) + ActPump * ( diam * diam * vrat [ 0 ] * ( cli0 - cl [ 0 ] ) / tau_NKCC1 ) + ( diam * diam * vrat [ 0 ] * ( cli0 - cl [ 0 ] ) / tau_passive ) )*/
 f_flux = b_flux = 0.;
 _RHS1( 0 +  0) += (b_flux =   ( ( icl * PI * diam / FARADAY ) + ActPump * ( diam * diam * vrat [ 0 ] * ( cli0 - cl [ 0 ] ) / tau_NKCC1 ) + ( diam * diam * vrat [ 0 ] * ( cli0 - cl [ 0 ] ) / tau_passive ) ) );
 /*FLUX*/
  {int  _li ;for ( _li = 0 ; _li <= 4 - 2 ; _li ++ ) {
     /* ~ cl [ _li ] <-> cl [ _li + 1 ] ( DCl * _zfrat [ _li + 1 ] , DCl * _zfrat [ _li + 1 ] )*/
 f_flux =  DCl * _zfrat [ _li + 1 ] * cl [ _li] ;
 b_flux =  DCl * _zfrat [ _li + 1 ] * cl [ _li + 1] ;
 _RHS1( 0 +  _li) -= (f_flux - b_flux);
 _RHS1( 0 +  _li + 1) += (f_flux - b_flux);
 
 _term =  DCl * _zfrat [ _li + 1 ] ;
 _MATELM1( 0 +  _li ,0 +  _li)  += _term;
 _MATELM1( 0 +  _li + 1 ,0 +  _li)  -= _term;
 _term =  DCl * _zfrat [ _li + 1 ] ;
 _MATELM1( 0 +  _li ,0 +  _li + 1)  -= _term;
 _MATELM1( 0 +  _li + 1 ,0 +  _li + 1)  += _term;
 /*REACTION*/
  } }
   cli = cl [ 0 ] ;
   /* COMPARTMENT _li , diam * diam * vrat [ ((int) _i ) ] {
     hco3 }
   */
 /* LONGITUDINAL_DIFFUSION _li , DCl * diam * diam * vrat [ ((int) _i ) ] {
     hco3 }
   */
 /* ~ hco3 [ 0 ] < < ( ( ihco3 * PI * diam / FARADAY ) + ( diam * diam * vrat [ 0 ] * ( hco3i0 - hco3 [ 0 ] ) / tau_hco3 ) )*/
 f_flux = b_flux = 0.;
 _RHS1( 4 +  0) += (b_flux =   ( ( ihco3 * PI * diam / FARADAY ) + ( diam * diam * vrat [ 0 ] * ( hco3i0 - hco3 [ 0 ] ) / tau_hco3 ) ) );
 /*FLUX*/
  {int  _li ;for ( _li = 0 ; _li <= 4 - 2 ; _li ++ ) {
     /* ~ hco3 [ _li ] <-> hco3 [ _li + 1 ] ( DCl * _zfrat [ _li + 1 ] , DCl * _zfrat [ _li + 1 ] )*/
 f_flux =  DCl * _zfrat [ _li + 1 ] * hco3 [ _li] ;
 b_flux =  DCl * _zfrat [ _li + 1 ] * hco3 [ _li + 1] ;
 _RHS1( 4 +  _li) -= (f_flux - b_flux);
 _RHS1( 4 +  _li + 1) += (f_flux - b_flux);
 
 _term =  DCl * _zfrat [ _li + 1 ] ;
 _MATELM1( 4 +  _li ,4 +  _li)  += _term;
 _MATELM1( 4 +  _li + 1 ,4 +  _li)  -= _term;
 _term =  DCl * _zfrat [ _li + 1 ] ;
 _MATELM1( 4 +  _li ,4 +  _li + 1)  -= _term;
 _MATELM1( 4 +  _li + 1 ,4 +  _li + 1)  += _term;
 /*REACTION*/
  } }
   hco3i = hco3 [ 0 ] ;
     } return _reset;
 }
 
/*CVODE ode begin*/
 static int _ode_spec1() {_reset=0;{
 double b_flux, f_flux, _term; int _i;
 {int _i; for(_i=0;_i<8;_i++) _p[_dlist1[_i]] = 0.0;}
 if ( cli0 >= cl [ 0 ] ) {
   ActPump = 1.0 ;
   }
 else {
   ActPump = 0.0 ;
   }
 /* COMPARTMENT _li , diam * diam * vrat [ ((int) _i ) ] {
   cl }
 */
 /* LONGITUDINAL_DIFFUSION _li , DCl * diam * diam * vrat [ ((int) _i ) ] {
   cl }
 */
 /* ~ cl [ 0 ] < < ( ( icl * PI * diam / FARADAY ) + ActPump * ( diam * diam * vrat [ 0 ] * ( cli0 - cl [ 0 ] ) / tau_NKCC1 ) + ( diam * diam * vrat [ 0 ] * ( cli0 - cl [ 0 ] ) / tau_passive ) )*/
 f_flux = b_flux = 0.;
 Dcl [ 0] += (b_flux =   ( ( icl * PI * diam / FARADAY ) + ActPump * ( diam * diam * vrat [ 0 ] * ( cli0 - cl [ 0 ] ) / tau_NKCC1 ) + ( diam * diam * vrat [ 0 ] * ( cli0 - cl [ 0 ] ) / tau_passive ) ) );
 /*FLUX*/
  {int  _li ;for ( _li = 0 ; _li <= 4 - 2 ; _li ++ ) {
   /* ~ cl [ _li ] <-> cl [ _li + 1 ] ( DCl * _zfrat [ _li + 1 ] , DCl * _zfrat [ _li + 1 ] )*/
 f_flux =  DCl * _zfrat [ _li + 1 ] * cl [ _li] ;
 b_flux =  DCl * _zfrat [ _li + 1 ] * cl [ _li + 1] ;
 Dcl [ _li] -= (f_flux - b_flux);
 Dcl [ _li + 1] += (f_flux - b_flux);
 
 /*REACTION*/
  } }
 cli = cl [ 0 ] ;
 /* COMPARTMENT _li , diam * diam * vrat [ ((int) _i ) ] {
   hco3 }
 */
 /* LONGITUDINAL_DIFFUSION _li , DCl * diam * diam * vrat [ ((int) _i ) ] {
   hco3 }
 */
 /* ~ hco3 [ 0 ] < < ( ( ihco3 * PI * diam / FARADAY ) + ( diam * diam * vrat [ 0 ] * ( hco3i0 - hco3 [ 0 ] ) / tau_hco3 ) )*/
 f_flux = b_flux = 0.;
 Dhco3 [ 0] += (b_flux =   ( ( ihco3 * PI * diam / FARADAY ) + ( diam * diam * vrat [ 0 ] * ( hco3i0 - hco3 [ 0 ] ) / tau_hco3 ) ) );
 /*FLUX*/
  {int  _li ;for ( _li = 0 ; _li <= 4 - 2 ; _li ++ ) {
   /* ~ hco3 [ _li ] <-> hco3 [ _li + 1 ] ( DCl * _zfrat [ _li + 1 ] , DCl * _zfrat [ _li + 1 ] )*/
 f_flux =  DCl * _zfrat [ _li + 1 ] * hco3 [ _li] ;
 b_flux =  DCl * _zfrat [ _li + 1 ] * hco3 [ _li + 1] ;
 Dhco3 [ _li] -= (f_flux - b_flux);
 Dhco3 [ _li + 1] += (f_flux - b_flux);
 
 /*REACTION*/
  } }
 hco3i = hco3 [ 0 ] ;
 for (_i=0; _i < 4; _i++) { _p[_dlist1[_i + 0]] /= ( diam * diam * vrat [ ((int) _i ) ]);}
 for (_i=0; _i < 4; _i++) { _p[_dlist1[_i + 4]] /= ( diam * diam * vrat [ ((int) _i ) ]);}
   } return _reset;
 }
 
/*CVODE matsol*/
 static int _ode_matsol1() {_reset=0;{
 double b_flux, f_flux, _term; int _i;
   b_flux = f_flux = 0.;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<8;_i++){
  	_RHS1(_i) = _dt1*(_p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} 
for (_i=0; _i < 4; _i++) {
  	_RHS1(_i + 0) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 0, _i + 0) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 4; _i++) {
  	_RHS1(_i + 4) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 4, _i + 4) *= ( diam * diam * vrat [ ((int) _i ) ]);  } }
 if ( cli0 >= cl [ 0 ] ) {
 ActPump = 1.0 ;
 }
 else {
 ActPump = 0.0 ;
 }
 /* COMPARTMENT _li , diam * diam * vrat [ ((int) _i ) ] {
 cl }
 */
 /* LONGITUDINAL_DIFFUSION _li , DCl * diam * diam * vrat [ ((int) _i ) ] {
 cl }
 */
 /* ~ cl [ 0 ] < < ( ( icl * PI * diam / FARADAY ) + ActPump * ( diam * diam * vrat [ 0 ] * ( cli0 - cl [ 0 ] ) / tau_NKCC1 ) + ( diam * diam * vrat [ 0 ] * ( cli0 - cl [ 0 ] ) / tau_passive ) )*/
 /*FLUX*/
  {int  _li ;for ( _li = 0 ; _li <= 4 - 2 ; _li ++ ) {
 /* ~ cl [ _li ] <-> cl [ _li + 1 ] ( DCl * _zfrat [ _li + 1 ] , DCl * _zfrat [ _li + 1 ] )*/
 _term =  DCl * _zfrat [ _li + 1 ] ;
 _MATELM1( 0 +  _li ,0 +  _li)  += _term;
 _MATELM1( 0 +  _li + 1 ,0 +  _li)  -= _term;
 _term =  DCl * _zfrat [ _li + 1 ] ;
 _MATELM1( 0 +  _li ,0 +  _li + 1)  -= _term;
 _MATELM1( 0 +  _li + 1 ,0 +  _li + 1)  += _term;
 /*REACTION*/
  } }
 cli = cl [ 0 ] ;
 /* COMPARTMENT _li , diam * diam * vrat [ ((int) _i ) ] {
 hco3 }
 */
 /* LONGITUDINAL_DIFFUSION _li , DCl * diam * diam * vrat [ ((int) _i ) ] {
 hco3 }
 */
 /* ~ hco3 [ 0 ] < < ( ( ihco3 * PI * diam / FARADAY ) + ( diam * diam * vrat [ 0 ] * ( hco3i0 - hco3 [ 0 ] ) / tau_hco3 ) )*/
 /*FLUX*/
  {int  _li ;for ( _li = 0 ; _li <= 4 - 2 ; _li ++ ) {
 /* ~ hco3 [ _li ] <-> hco3 [ _li + 1 ] ( DCl * _zfrat [ _li + 1 ] , DCl * _zfrat [ _li + 1 ] )*/
 _term =  DCl * _zfrat [ _li + 1 ] ;
 _MATELM1( 4 +  _li ,4 +  _li)  += _term;
 _MATELM1( 4 +  _li + 1 ,4 +  _li)  -= _term;
 _term =  DCl * _zfrat [ _li + 1 ] ;
 _MATELM1( 4 +  _li ,4 +  _li + 1)  -= _term;
 _MATELM1( 4 +  _li + 1 ,4 +  _li + 1)  += _term;
 /*REACTION*/
  } }
 hco3i = hco3 [ 0 ] ;
   } return _reset;
 }
 
/*CVODE end*/
 
static int _ode_count(int _type){ return 8;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  icl = _ion_icl;
  cli = _ion_cli;
  ihco3 = _ion_ihco3;
  hco3i = _ion_hco3i;
     _ode_spec1 ();
  _ion_cli = cli;
  _ion_hco3i = hco3i;
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 8; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 static void _ode_synonym(int _cnt, double** _pp, Datum** _ppd) { 
 	int _i; 
	for (_i=0; _i < _cnt; ++_i) {_p = _pp[_i]; _ppvar = _ppd[_i];
 _ion_cli =  cl [ 0 ] ;
 _ion_hco3i =  hco3 [ 0 ] ;
 }}
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _cvode_sparse(&_cvsparseobj1, 8, _dlist1, _p, _ode_matsol1, &_coef1);
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  icl = _ion_icl;
  cli = _ion_cli;
  ihco3 = _ion_ihco3;
  hco3i = _ion_hco3i;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_cl_sym, _ppvar, 0, 3);
   nrn_update_ion_pointer(_cl_sym, _ppvar, 1, 1);
   nrn_update_ion_pointer(_cl_sym, _ppvar, 3, 4);
   nrn_update_ion_pointer(_hco3_sym, _ppvar, 4, 3);
   nrn_update_ion_pointer(_hco3_sym, _ppvar, 5, 1);
   nrn_update_ion_pointer(_hco3_sym, _ppvar, 7, 4);
 }
 static void* _difspace1;
extern double nrn_nernst_coef();
static double _difcoef1(int _i, double* _p, Datum* _ppvar, double* _pdvol, double* _pdfcdc, Datum* _thread, NrnThread* _nt) {
   *_pdvol =  diam * diam * vrat [ ((int) _i ) ] ;
 if (_i ==  0) {
  *_pdfcdc = nrn_nernst_coef(_type_icl)*( ( ( _ion_dicldv  * PI * diam / FARADAY ) + ActPump * ( diam * diam * vrat [ 0 ] * ( cli0 - cl [ 0 ] ) / tau_NKCC1 ) + ( diam * diam * vrat [ 0 ] * ( cli0 - cl [ 0 ] ) / tau_passive ) ));
 }else{ *_pdfcdc=0.;}
   return DCl * diam * diam * vrat [ ((int) _i ) ] ;
}
 static void* _difspace2;
extern double nrn_nernst_coef();
static double _difcoef2(int _i, double* _p, Datum* _ppvar, double* _pdvol, double* _pdfcdc, Datum* _thread, NrnThread* _nt) {
   *_pdvol =  diam * diam * vrat [ ((int) _i ) ] ;
 if (_i ==  0) {
  *_pdfcdc = nrn_nernst_coef(_type_ihco3)*( ( ( _ion_dihco3dv  * PI * diam / FARADAY ) + ( diam * diam * vrat [ 0 ] * ( hco3i0 - hco3 [ 0 ] ) / tau_hco3 ) ));
 }else{ *_pdfcdc=0.;}
   return DCl * diam * diam * vrat [ ((int) _i ) ] ;
}
 static void _difusfunc(ldifusfunc2_t _f, NrnThread* _nt) {int _i;
  for (_i=0; _i < 4; ++_i) (*_f)(_mechtype, _difcoef1, &_difspace1, _i,  6, 20 , _nt);
  for (_i=0; _i < 4; ++_i) (*_f)(_mechtype, _difcoef2, &_difspace2, _i,  10, 24 , _nt);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
 for (_i=0; _i<4; _i++) cl[_i] = cl0;
 for (_i=0; _i<4; _i++) hco3[_i] = hco30;
 {
   if ( _zfactors_done  == 0.0 ) {
     _zfactors_done = 1.0 ;
     factors ( _threadargs_ ) ;
     }
   cli = cli_Start ;
   hco3i = hco3i0 ;
   hco3o = hco3o0 ;
   {int  _li ;for ( _li = 0 ; _li <= 4 - 1 ; _li ++ ) {
     cl [ _li ] = cli ;
     } }
   {int  _li ;for ( _li = 0 ; _li <= 4 - 1 ; _li ++ ) {
     hco3 [ _li ] = hco3i ;
     } }
   ehco3_help = log ( hco3i / hco3o ) * ( 1000.0 ) * ( celsius + 273.15 ) * R / F ;
   ecl_help = log ( cli / clo0 ) * ( 1000.0 ) * ( celsius + 273.15 ) * R / F ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  icl = _ion_icl;
  cli = _ion_cli;
  ihco3 = _ion_ihco3;
  hco3i = _ion_hco3i;
 initmodel();
  _ion_cli = cli;
  nrn_wrote_conc(_cl_sym, (&(_ion_cli)) - 1, _style_cl);
  _ion_hco3i = hco3i;
  nrn_wrote_conc(_hco3_sym, (&(_ion_hco3i)) - 1, _style_hco3);
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{
} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 
}}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
double _dtsav = dt;
if (secondorder) { dt *= 0.5; }
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  icl = _ion_icl;
  cli = _ion_cli;
  ihco3 = _ion_ihco3;
  hco3i = _ion_hco3i;
 { error = sparse(&_sparseobj1, 8, _slist1, _dlist1, _p, &t, dt, state,&_coef1, _linmat1);
 if(error){fprintf(stderr,"at line 75 in file cldif_CA3_NKCC1_HCO3.mod:\n		SOLVE state METHOD sparse\n"); nrn_complain(_p); abort_run(error);}
    if (secondorder) {
    int _i;
    for (_i = 0; _i < 8; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 } {
   ecl_help = log ( cli / clo0 ) * ( 1000.0 ) * ( celsius + 273.15 ) * R / F ;
   ehco3_help = log ( hco3i / hco3o0 ) * ( 1000.0 ) * ( celsius + 273.15 ) * R / F ;
   }
  _ion_cli = cli;
  _ion_hco3i = hco3i;
}}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 for(_i=0;_i<4;_i++){_slist1[0+_i] = cl_columnindex + _i;  _dlist1[0+_i] = Dcl_columnindex + _i;}
 for(_i=0;_i<4;_i++){_slist1[4+_i] = hco3_columnindex + _i;  _dlist1[4+_i] = Dhco3_columnindex + _i;}
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/dustinmiao/Documents/school/rise/RISE-gaba-neurogen/kolbaev et. al test model/Models_Kolbaev et al Scientific Reports Revision/cldif_CA3_NKCC1_HCO3.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "\n"
  "Chloride accumulation and diffusion with decay (time constant tau) to resting level cli0.\n"
  "The decay approximates a reversible chloride pump with first order kinetics.\n"
  "To eliminate the chloride pump, just use this hoc statement\n"
  "To make the time constant effectively \"infinite\".\n"
  "tau and the resting level are both RANGE variables\n"
  "\n"
  "Diffusion model is modified from Ca diffusion model in Hines & Carnevale:\n"
  "Expanding NEURON with NMODL, Neural Computation 12: 839-851, 2000 (Example 8)\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX cldif_CA3_NKCC1_HCO3\n"
  "	USEION cl READ icl WRITE cli VALENCE -1 : Ion cl, use cl current to calculate cl internal concentration\n"
  "	USEION hco3 READ ihco3 WRITE hco3i VALENCE -1: Ion HCO3, use HCO3 internal concentration to calculate the external concentration\n"
  "	GLOBAL vrat		:vrat must be GLOBAL, so it does not change with position. vrat = volumes of concentric shells\n"
  "	RANGE tau, cli0, clo0, hco3i0, hco3o0, egaba, delta_egaba, init_egaba, ehco3_help, ecl_help : all of these change with position\n"
  "}\n"
  "\n"
  "DEFINE Nannuli 4\n"
  "\n"
  "UNITS {\n"
  "	(molar) = (1/liter)\n"
  "	(mM) = (millimolar)\n"
  "	(um) = (micron)\n"
  "	(mA) = (milliamp)\n"
  "	(mV)    = (millivolt)\n"
  "	FARADAY = (faraday) (10000 coulomb)\n"
  "	PI = (pi) (1)\n"
  "	F = (faraday) (coulombs)\n"
  "	R = (k-mole)  (joule/degC)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	DCl = 2 (um2/ms) : Kuner & Augustine, Neuron 27: 447 : diffusion coefficient of cl\n"
  "	tau_NKCC1 = 174000 (ms)   : 174 s From Kolbaev, Lombardi kilb (in Prep) - kinetics after Cl decline\n"
  "	tau_passive = 321000 (ms) : 321 s From Kolbaev, Lombardi Kilb (in prep) - kinetics after bumetanid washin\n"
  "        tau_hco3 = 1000 (ms) : tau for Bicarbonate, just an arbitrary value\n"
  "	cli0 = 50 (mM) : basal Cl internal concentration\n"
  "	cli_Start = 10 (mM) :Cl- concentration at start\n"
  "	clo0 = 133.5 (mM) : basal Cl external concentration\n"
  "	hco3i0 = 16	(mM) : basal HCO3 internal concentration\n"
  "	hco3o0 = 26	(mM) : basal HCO3 external concentration\n"
  "	hco3i_Start = 16 (mM) : Cl- concentration at start\n"
  "	celsius = 31    (degC)\n"
  "\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	diam 	(um)\n"
  "	icl 	(mA/cm2) : Cl current\n"
  "        ihco3 	(mA/cm2) : HCO3- current current\n"
  "	cli 	(mM) : Cl internal concentration\n"
  "	hco3i	(mM) : HCO3 internal concentration\n"
  "	hco3o	(mM) : HCO3 external concentration\n"
  "	vrat[Nannuli]	: numeric value of vrat[i] equals the volume\n"
  "			: of annulus i of a 1um diameter cylinder\n"
  "			: multiply by diam^2 to get volume per um length\n"
  "	ehco3_help 	(mV)\n"
  "	ecl_help	(mV)\n"
  "	ActPump   :Binary value that defines if active inward pumping of passive outward diffusion\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	: cl[0] is equivalent to cli\n"
  "	: cl[] are very small, so specify absolute tolerance\n"
  "	cl[Nannuli]	(mM) <1e-10>\n"
  "        hco3[Nannuli]	(mM) <1e-10>\n"
  "}\n"
  "\n"
  "\n"
  "BREAKPOINT {\n"
  "		SOLVE state METHOD sparse\n"
  "		ecl_help = log(cli/clo0)*(1000)*(celsius + 273.15)*R/F\n"
  "                ehco3_help = log(hco3i/hco3o0)*(1000)*(celsius + 273.15)*R/F\n"
  "}\n"
  "\n"
  "LOCAL factors_done\n"
  "\n"
  "INITIAL {\n"
  "	if (factors_done == 0) {  	: flag becomes 1 in the first segment\n"
  "		factors_done = 1	: all subsequent segments will have\n"
  "		factors()		: vrat = 0 unless vrat is GLOBAL. We make sure that vrat is applied to the shell volumes\n"
  "	}\n"
  "	cli = cli_Start\n"
  "	hco3i = hco3i0\n"
  "	hco3o = hco3o0\n"
  "	FROM i=0 TO Nannuli-1 { : So that at the begining the Cl [] is the same in all shells ( steady state)\n"
  "		cl[i] = cli\n"
  "	}\n"
  "        FROM i=0 TO Nannuli-1 { : So that at the begining the HCO3 [] is the same in all shells ( steady state)\n"
  "		hco3[i] = hco3i\n"
  "	}\n"
  "	ehco3_help = log(hco3i/hco3o)*(1000)*(celsius + 273.15)*R/F : Nerst eq for HCO3 at time 0\n"
  "	ecl_help = log(cli/clo0)*(1000)*(celsius + 273.15)*R/F\n"
  "}\n"
  "\n"
  "LOCAL frat[Nannuli]	: scales the rate constants for model geometry\n"
  "\n"
  "PROCEDURE factors() {\n"
  "	LOCAL r, dr2\n"
  "	r = 1/2			: starts at edge (half diam), diam = 1, length = 1\n"
  "	dr2 = r/(Nannuli-1)/2	: full thickness of outermost annulus,\n"
  "				: half thickness of all other annuli\n"
  "	vrat[0] = 0\n"
  "	frat[0] = 2*r		: = diam\n"
  "	FROM i=0 TO Nannuli-2 {\n"
  "		vrat[i] = vrat[i] + PI*(r-dr2/2)*2*dr2	: interior half\n"
  "		r = r - dr2\n"
  "		frat[i+1] = 2*PI*r/(2*dr2)	: outer radius of annulus Ai+1/delta_r=2PI*r*1/delta_r\n"
  "						: div by distance between centers\n"
  "		r = r - dr2\n"
  "		vrat[i+1] = PI*(r+dr2/2)*2*dr2	: outer half of annulus\n"
  "	}\n"
  "}\n"
  "\n"
  "KINETIC state {\n"
  "    if (cli0 >= cl[0]) { : Under this condition the NKCC1 mediates active Cl- uptake ( positive inward flux)\n"
  "		  ActPump = 1\n"
  "		}\n"
  "		else {     : Under this condition NKCC1 should be not functional ( negative inward flux)\n"
  "		  ActPump = 0\n"
  "		}\n"
  "\n"
  "  	COMPARTMENT i, diam*diam*vrat[i] {cl}\n"
  "		LONGITUDINAL_DIFFUSION i, DCl*diam*diam*vrat[i] {cl}\n"
  "				~ cl[0] << ((icl*PI*diam/FARADAY) + ActPump*(diam*diam*vrat[0]*(cli0 - cl[0])/tau_NKCC1) + (diam*diam*vrat[0]*(cli0 - cl[0])/tau_passive)) : icl is Cl- influx\n"
  "	 	FROM i=0 TO Nannuli-2 {\n"
  "		~ cl[i] <-> cl[i+1]	(DCl*frat[i+1], DCl*frat[i+1])\n"
  "                }\n"
  "	        cli = cl[0]\n"
  "\n"
  "        COMPARTMENT i, diam*diam*vrat[i] {hco3}\n"
  "		LONGITUDINAL_DIFFUSION i, DCl*diam*diam*vrat[i] {hco3}\n"
  "				~ hco3[0] << ((ihco3*PI*diam/FARADAY)  + (diam*diam*vrat[0]*(hco3i0 - hco3[0])/tau_hco3)) : ihco3 is HCO3- influx\n"
  "	 	FROM i=0 TO Nannuli-2 {\n"
  "		~ hco3[i] <-> hco3[i+1]	(DCl*frat[i+1], DCl*frat[i+1])\n"
  "                }\n"
  "	        hco3i = hco3[0]\n"
  "}\n"
  ;
#endif
