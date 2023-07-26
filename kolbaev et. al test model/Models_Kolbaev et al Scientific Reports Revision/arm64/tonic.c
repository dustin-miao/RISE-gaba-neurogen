/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
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
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__tonic
#define _nrn_initial _nrn_initial__tonic
#define nrn_cur _nrn_cur__tonic
#define _nrn_current _nrn_current__tonic
#define nrn_jacob _nrn_jacob__tonic
#define nrn_state _nrn_state__tonic
#define _net_receive _net_receive__tonic 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define leak _p[0]
#define leak_columnindex 0
#define P _p[1]
#define P_columnindex 1
#define HCO3e _p[2]
#define HCO3e_columnindex 2
#define HCO3i _p[3]
#define HCO3i_columnindex 3
#define icl _p[4]
#define icl_columnindex 4
#define i _p[5]
#define i_columnindex 5
#define ihco3 _p[6]
#define ihco3_columnindex 6
#define e _p[7]
#define e_columnindex 7
#define ehco3 _p[8]
#define ehco3_columnindex 8
#define ecl _p[9]
#define ecl_columnindex 9
#define v _p[10]
#define v_columnindex 10
#define _g _p[11]
#define _g_columnindex 11
#define _ion_ecl	*_ppvar[0]._pval
#define _ion_icl	*_ppvar[1]._pval
#define _ion_dicldv	*_ppvar[2]._pval
 
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
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
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
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_tonic", _hoc_setdata,
 0, 0
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "leak_tonic", 0, 1e+09,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "leak_tonic", "siemens/cm2",
 "HCO3e_tonic", "mM",
 "HCO3i_tonic", "mM",
 "icl_tonic", "milliamp/cm2",
 "i_tonic", "milliamp/cm2",
 "ihco3_tonic", "milliamp/cm2",
 "e_tonic", "mV",
 "ehco3_tonic", "mV",
 0,0
};
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"tonic",
 "leak_tonic",
 "P_tonic",
 "HCO3e_tonic",
 "HCO3i_tonic",
 0,
 "icl_tonic",
 "i_tonic",
 "ihco3_tonic",
 "e_tonic",
 "ehco3_tonic",
 0,
 0,
 0};
 static Symbol* _cl_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 12, _prop);
 	/*initialize range parameters*/
 	leak = 0.01;
 	P = 0.18;
 	HCO3e = 26;
 	HCO3i = 16;
 	_prop->param = _p;
 	_prop->param_size = 12;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_cl_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ecl */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* icl */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dicldv */
 
}
 static void _initlists();
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _tonic_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("cl", -1.0);
 	_cl_sym = hoc_lookup("cl_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 12, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "cl_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "cl_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "cl_ion");
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 tonic /Users/dustinmiao/Documents/school/rise/RISE-gaba-neurogen/kolbaev et. al test model/Models_Kolbaev et al Scientific Reports Revision/tonic.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 
#define F _nrnunit_F[_nrnunit_use_legacy_]
static double _nrnunit_F[2] = {0x1.78e555060882cp+16, 96485.3}; /* 96485.3321233100141 */
 
#define R _nrnunit_R[_nrnunit_use_legacy_]
static double _nrnunit_R[2] = {0x1.0a1013e8990bep+3, 8.3145}; /* 8.3144626181532395 */
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_cl_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_cl_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_cl_sym, _ppvar, 2, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
 {
   ehco3 = log ( HCO3i / HCO3e ) * ( 1000.0 ) * ( celsius + 273.15 ) * R / F ;
   e = P * ehco3 + ( 1.0 - P ) * ecl ;
   }

}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
  ecl = _ion_ecl;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   icl = ( 1.0 - P ) * leak * ( v - ecl ) ;
   ihco3 = P * leak * ( v - ehco3 ) ;
   i = icl + ihco3 ;
   e = P * ehco3 + ( 1.0 - P ) * ecl ;
   }
 _current += icl;
 _current += ihco3;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
  ecl = _ion_ecl;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dicl;
  _dicl = icl;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dicldv += (_dicl - icl)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_icl += icl ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 
}
 
}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type) {

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/dustinmiao/Documents/school/rise/RISE-gaba-neurogen/kolbaev et. al test model/Models_Kolbaev et al Scientific Reports Revision/tonic.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "Modified Proddutur A, Yu J, Elgammal FS, Santhakumar V (2013)\n"
  "Tonic Inhibition with changing Cl- concentration\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON{\n"
  "	SUFFIX tonic\n"
  "	USEION cl READ ecl WRITE icl VALENCE -1\n"
  "	NONSPECIFIC_CURRENT ihco3\n"
  "	RANGE  icl, leak, P, e,ihco3,ehco3 \n"
  "	RANGE P, HCO3e, HCO3i, i\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(uS)  = (micromho)\n"
  "	(nA)  = (nanoamp)\n"
  "	(mV)  = (millivolt)\n"
  "	(mM)    = (milli/liter)\n"
  "	F 	  = (faraday) (coulombs)\n"
  "	R     = (k-mole)  (joule/degC)\n"
  "}\n"
  "\n"
  "\n"
  "PARAMETER {\n"
  "	leak    = 0.01 (siemens/cm2) <0, 1e9>\n"
  "	P    = 0.18		: HCO3/Cl relative permeability\n"
  "	celsius = 31    (degC)\n"
  "	HCO3e   = 26	(mM)	: extracellular HCO3- concentration\n"
  "	HCO3i   = 16	(mM)	: intracellular HCO3- concentration\n"
  "}\n"
  "\n"
  "ASSIGNED{\n"
  "	icl (milliamp/cm2)\n"
  "	i (milliamp/cm2)\n"
  "	ihco3 (milliamp/cm2)\n"
  "	v	(mV)		: postsynaptic voltage \n"
  "	e	(mV)		: reversal potential for GABAR	\n"
  "	ecl	(mV)		: equilibrium potential for Cl-\n"
  "    ehco3	(mV)		: equilibrium potential for HCO3-	\n"
  "}\n"
  "\n"
  "INITIAL { \n"
  "\n"
  "	ehco3 = log(HCO3i/HCO3e)*(1000)*(celsius + 273.15)*R/F\n"
  "	e = P*ehco3 + (1-P)*ecl\n"
  "}\n"
  "\n"
  "\n"
  "BREAKPOINT {\n"
  "	icl = (1-P)*leak*(v-ecl)\n"
  "	ihco3 = P*leak*(v-ehco3)\n"
  "	i = icl + ihco3\n"
  "	e = P*ehco3 + (1-P)*ecl\n"
  "}\n"
  ;
#endif
