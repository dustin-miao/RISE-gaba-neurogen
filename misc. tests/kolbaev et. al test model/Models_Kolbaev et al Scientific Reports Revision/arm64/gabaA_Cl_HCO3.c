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
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__gaba
#define _nrn_initial _nrn_initial__gaba
#define nrn_cur _nrn_cur__gaba
#define _nrn_current _nrn_current__gaba
#define nrn_jacob _nrn_jacob__gaba
#define nrn_state _nrn_state__gaba
#define _net_receive _net_receive__gaba 
#define state state__gaba 
 
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
#define tau1 _p[0]
#define tau1_columnindex 0
#define tau2 _p[1]
#define tau2_columnindex 1
#define P _p[2]
#define P_columnindex 2
#define icl _p[3]
#define icl_columnindex 3
#define ihco3 _p[4]
#define ihco3_columnindex 4
#define i _p[5]
#define i_columnindex 5
#define g _p[6]
#define g_columnindex 6
#define e _p[7]
#define e_columnindex 7
#define A _p[8]
#define A_columnindex 8
#define B _p[9]
#define B_columnindex 9
#define factor _p[10]
#define factor_columnindex 10
#define ecl _p[11]
#define ecl_columnindex 11
#define ehco3 _p[12]
#define ehco3_columnindex 12
#define DA _p[13]
#define DA_columnindex 13
#define DB _p[14]
#define DB_columnindex 14
#define _g _p[15]
#define _g_columnindex 15
#define _tsav _p[16]
#define _tsav_columnindex 16
#define _nd_area  *_ppvar[0]._pval
#define _ion_ecl	*_ppvar[2]._pval
#define _ion_icl	*_ppvar[3]._pval
#define _ion_dicldv	*_ppvar[4]._pval
#define _ion_ehco3	*_ppvar[5]._pval
#define _ion_ihco3	*_ppvar[6]._pval
#define _ion_dihco3dv	*_ppvar[7]._pval
 
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

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(Object* _ho) { void* create_point_process(int, Object*);
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt(void*);
 static double _hoc_loc_pnt(void* _vptr) {double loc_point_process(int, void*);
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(void* _vptr) {double has_loc_point(void*);
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(void* _vptr) {
 double get_loc_point_process(void*); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 0, 0
};
 /* declare global and static user variables */
#define total total_gaba
 double total = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "tau2", 1e-09, 1e+09,
 "tau1", 1e-09, 1e+09,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "total_gaba", "uS",
 "tau1", "ms",
 "tau2", "ms",
 "A", "uS",
 "B", "uS",
 "icl", "nA",
 "ihco3", "nA",
 "i", "nA",
 "g", "uS",
 "e", "mV",
 0,0
};
 static double A0 = 0;
 static double B0 = 0;
 static double delta_t = 0.01;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "total_gaba", &total_gaba,
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
 static void _hoc_destroy_pnt(void* _vptr) {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[8]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"gaba",
 "tau1",
 "tau2",
 "P",
 0,
 "icl",
 "ihco3",
 "i",
 "g",
 "e",
 0,
 "A",
 "B",
 0,
 0};
 static Symbol* _cl_sym;
 static Symbol* _hco3_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 17, _prop);
 	/*initialize range parameters*/
 	tau1 = 0.1;
 	tau2 = 80;
 	P = 0.18;
  }
 	_prop->param = _p;
 	_prop->param_size = 17;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 9, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_cl_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[2]._pval = &prop_ion->param[0]; /* ecl */
 	_ppvar[3]._pval = &prop_ion->param[3]; /* icl */
 	_ppvar[4]._pval = &prop_ion->param[4]; /* _ion_dicldv */
 prop_ion = need_memb(_hco3_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[5]._pval = &prop_ion->param[0]; /* ehco3 */
 	_ppvar[6]._pval = &prop_ion->param[3]; /* ihco3 */
 	_ppvar[7]._pval = &prop_ion->param[4]; /* _ion_dihco3dv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _net_receive(Point_process*, double*, double);
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _gabaA_Cl_HCO3_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("cl", -1.0);
 	ion_reg("hco3", -1.0);
 	_cl_sym = hoc_lookup("cl_ion");
 	_hco3_sym = hoc_lookup("hco3_ion");
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 0,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 17, 9);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "cl_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cl_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "cl_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "hco3_ion");
  hoc_register_dparam_semantics(_mechtype, 6, "hco3_ion");
  hoc_register_dparam_semantics(_mechtype, 7, "hco3_ion");
  hoc_register_dparam_semantics(_mechtype, 8, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 gaba /Users/dustinmiao/Documents/school/rise/RISE-gaba-neurogen/kolbaev et. al test model/Models_Kolbaev et al Scientific Reports Revision/gabaA_Cl_HCO3.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 
#define F _nrnunit_F[_nrnunit_use_legacy_]
static double _nrnunit_F[2] = {0x1.78e555060882cp+16, 96485.3}; /* 96485.3321233100141 */
 
#define R _nrnunit_R[_nrnunit_use_legacy_]
static double _nrnunit_R[2] = {0x1.0a1013e8990bep+3, 8.3145}; /* 8.3144626181532395 */
static int _reset;
static char *modelname = "GABAergic conductance with changing Cl- concentration";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int state(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   DA = - A / tau1 ;
   DB = - B / tau2 ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 DA = DA  / (1. - dt*( ( - 1.0 ) / tau1 )) ;
 DB = DB  / (1. - dt*( ( - 1.0 ) / tau2 )) ;
  return 0;
}
 /*END CVODE*/
 static int state () {_reset=0;
 {
    A = A + (1. - exp(dt*(( - 1.0 ) / tau1)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau1 ) - A) ;
    B = B + (1. - exp(dt*(( - 1.0 ) / tau2)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau2 ) - B) ;
   }
  return 0;
}
 
static void _net_receive (Point_process* _pnt, double* _args, double _lflag) 
{    _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t; {
     if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = A;
    double __primary = (A + _args[0] * factor) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau1 ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau1 ) - __primary );
    A += __primary;
  } else {
 A = A + _args[0] * factor ;
     }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = B;
    double __primary = (B + _args[0] * factor) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau2 ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau2 ) - __primary );
    B += __primary;
  } else {
 B = B + _args[0] * factor ;
     }
 total = total + _args[0] ;
   } }
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ecl = _ion_ecl;
  ehco3 = _ion_ehco3;
     _ode_spec1 ();
   }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
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
  ecl = _ion_ecl;
  ehco3 = _ion_ehco3;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_cl_sym, _ppvar, 2, 0);
   nrn_update_ion_pointer(_cl_sym, _ppvar, 3, 3);
   nrn_update_ion_pointer(_cl_sym, _ppvar, 4, 4);
   nrn_update_ion_pointer(_hco3_sym, _ppvar, 5, 0);
   nrn_update_ion_pointer(_hco3_sym, _ppvar, 6, 3);
   nrn_update_ion_pointer(_hco3_sym, _ppvar, 7, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  A = A0;
  B = B0;
 {
   double _ltp ;
 total = 0.0 ;
   if ( tau1 / tau2 > .9999 ) {
     tau1 = .9999 * tau2 ;
     }
   A = 0.0 ;
   B = 0.0 ;
   _ltp = ( tau1 * tau2 ) / ( tau2 - tau1 ) * log ( tau2 / tau1 ) ;
   factor = - exp ( - _ltp / tau1 ) + exp ( - _ltp / tau2 ) ;
   factor = 1.0 / factor ;
   e = P / ( 1.0 + P ) * ehco3 + 1.0 / ( 1.0 + P ) * ecl ;
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
 _tsav = -1e20;
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
  ehco3 = _ion_ehco3;
 initmodel();
  }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   g = B - A ;
   icl = 1.0 / ( 1.0 + P ) * g * ( v - ecl ) ;
   ihco3 = P / ( 1.0 + P ) * g * ( v - ehco3 ) ;
   i = icl + ihco3 ;
   e = P / ( 1.0 + P ) * ehco3 + P / ( 1.0 + P ) * ecl ;
   }
 _current += icl;
 _current += ihco3;

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
  ecl = _ion_ecl;
  ehco3 = _ion_ehco3;
 _g = _nrn_current(_v + .001);
 	{ double _dihco3;
 double _dicl;
  _dicl = icl;
  _dihco3 = ihco3;
 _rhs = _nrn_current(_v);
  _ion_dicldv += (_dicl - icl)/.001 * 1.e2/ (_nd_area);
  _ion_dihco3dv += (_dihco3 - ihco3)/.001 * 1.e2/ (_nd_area);
 	}
 _g = (_g - _rhs)/.001;
  _ion_icl += icl * 1.e2/ (_nd_area);
  _ion_ihco3 += ihco3 * 1.e2/ (_nd_area);
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
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
  ecl = _ion_ecl;
  ehco3 = _ion_ehco3;
 { error =  state();
 if(error){fprintf(stderr,"at line 104 in file gabaA_Cl_HCO3.mod:\n\n"); nrn_complain(_p); abort_run(error);}
 }  }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = A_columnindex;  _dlist1[0] = DA_columnindex;
 _slist1[1] = B_columnindex;  _dlist1[1] = DB_columnindex;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/dustinmiao/Documents/school/rise/RISE-gaba-neurogen/kolbaev et. al test model/Models_Kolbaev et al Scientific Reports Revision/gabaA_Cl_HCO3.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "\n"
  "Synaptic GABAergic mechanism\n"
  "\n"
  "Reversal potential Egaba is changing according to [Cl-]i change (due to Cl- influx, which we hypothesize to be significant). Bicarbonate (HCO3) flows through the GABAR too, and therefore Egaba is also [HCO3]i/[HCO3]o -dependent\n"
  "igaba = icl + ihco3 (we assume icl and ihco3 to be mutually independent)\n"
  "\n"
  "Two state kinetic scheme synapse described by rise time tau1,\n"
  "and decay time constant tau2. The normalized peak condunductance is 1.\n"
  "Decay time MUST be greater than rise time.\n"
  "\n"
  "The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is\n"
  " A = a*exp(-t/tau1) and\n"
  " G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))\n"
  "	where tau1 < tau2\n"
  "\n"
  "If tau2-tau1 -> 0 then we have a alphasynapse.\n"
  "and if tau1 -> 0 then we have just single exponential decay.\n"
  "\n"
  "The factor is evaluated in the\n"
  "initial block such that an event of weight 1 generates a\n"
  "peak conductance of 1.\n"
  "\n"
  "Because the solution is a sum of exponentials, the\n"
  "coupled equations can be solved as a pair of independent equations\n"
  "by the more efficient cnexp method.\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "\n"
  "TITLE GABAergic conductance with changing Cl- concentration\n"
  "\n"
  "NEURON {\n"
  "	POINT_PROCESS gaba\n"
  "\n"
  "	USEION cl READ ecl WRITE icl VALENCE -1\n"
  "        USEION hco3 READ ehco3 WRITE ihco3 VALENCE -1\n"
  "\n"
  "	RANGE tau1, tau2, g\n"
  "	RANGE P, i\n"
  "\n"
  "	RANGE icl, ihco3, ehco3, e\n"
  "	GLOBAL total\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(mA)    = (milliamp)\n"
  "	(nA)    = (nanoamp)\n"
  "	(mV)    = (millivolt)\n"
  "	(uS)  = (micromho)\n"
  "	(mM)    = (milli/liter)\n"
  "	F = (faraday) (coulombs)\n"
  "	R = (k-mole)  (joule/degC)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	tau1	=.1	(ms)	<1e-9,1e9>\n"
  "	tau2	= 80	(ms)	<1e-9,1e9>\n"
  "\n"
  "	P       = 0.18		: HCO3/Cl relative permeability\n"
  "\n"
  "	celsius = 31    (degC)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v	(mV)		: postsynaptic voltage\n"
  "\n"
  "	icl	(nA)		: chloride current = 1/(1+P)*g*(v - ecl)\n"
  "	ihco3	(nA)		: bicarb current = P/(1+P)*g*(v - ehco3)\n"
  "	i	(nA)		: total current generated by this mechanism\n"
  "				: = icl + ihco3\n"
  "	g 	(uS)		: total conductance, split between bicarb (P/(1+P)*g)\n"
  "				: and chloride (1/(1+P)*g)\n"
  "	factor\n"
  "	total	(uS)\n"
  "\n"
  "	ecl	(mV)		: equilibrium potential for Cl-\n"
  "	ehco3	(mV)		: equilibrium potential for HCO3-\n"
  "\n"
  "	e	(mV)		: reversal potential for GABAR\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	A (uS)\n"
  "	B (uS)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	LOCAL tp\n"
  "	total = 0\n"
  "	if (tau1/tau2 > .9999) {\n"
  "		tau1 = .9999*tau2\n"
  "	}\n"
  "	A = 0\n"
  "	B = 0\n"
  "	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)\n"
  "	factor = -exp(-tp/tau1) + exp(-tp/tau2)\n"
  "	factor = 1/factor\n"
  "	e = P/(1+P)*ehco3 + 1/(1+P)*ecl\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE state METHOD cnexp\n"
  "\n"
  "	g = B - A\n"
  "\n"
  "	icl = 1/(1+P)*g*(v-ecl)\n"
  "\n"
  "	ihco3 = P/(1+P)*g*(v-ehco3)\n"
  "	i = icl + ihco3\n"
  "	e = P/(1+P)*ehco3 + P/(1+P)*ecl\n"
  "\n"
  "}\n"
  "\n"
  "DERIVATIVE state {\n"
  "	A' = -A/tau1\n"
  "	B' = -B/tau2\n"
  "}\n"
  "\n"
  "NET_RECEIVE(weight (uS)) {\n"
  "	A = A + weight*factor\n"
  "	B = B + weight*factor\n"
  "	total = total+weight\n"
  "}\n"
  ;
#endif
