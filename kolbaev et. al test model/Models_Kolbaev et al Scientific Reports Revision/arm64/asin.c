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
 
#define nrn_init _nrn_init__asin
#define _nrn_initial _nrn_initial__asin
#define nrn_cur _nrn_cur__asin
#define _nrn_current _nrn_current__asin
#define nrn_jacob _nrn_jacob__asin
#define nrn_state _nrn_state__asin
#define _net_receive _net_receive__asin 
 
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
 /* declaration of user functions */
 static void _hoc_exp_i(void);
 static void _hoc_mytime(void);
 static void _hoc_my_sin(void);
 static void _hoc_my_asin(void);
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
 "setdata_asin", _hoc_setdata,
 "exp_i_asin", _hoc_exp_i,
 "mytime_asin", _hoc_mytime,
 "my_sin_asin", _hoc_my_sin,
 "my_asin_asin", _hoc_my_asin,
 0, 0
};
#define exp_i exp_i_asin
#define mytime mytime_asin
#define my_sin my_sin_asin
#define my_asin my_asin_asin
 extern double exp_i( double );
 extern double mytime( );
 extern double my_sin( double );
 extern double my_asin( double );
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 0,0
};
 static double v = 0;
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
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"asin",
 0,
 0,
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 0, _prop);
 	/*initialize range parameters*/
 	_prop->param = _p;
 	_prop->param_size = 0;
 
}
 static void _initlists();
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _asin_reg() {
	int _vectorized = 0;
  _initlists();
 	register_mech(_mechanism, nrn_alloc,(void*)0, (void*)0, (void*)0, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 asin /Users/dustinmiao/Documents/school/rise/RISE-gaba-neurogen/kolbaev et. al test model/Models_Kolbaev et al Scientific Reports Revision/asin.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
/*VERBATIM*/
#include <math.h>
#include <time.h>
 
double my_asin (  double _larg ) {
   double _lmy_asin;
 
/*VERBATIM*/
	double ret;
	ret=asin(*getarg(1));
 _lmy_asin = ret ;
   
return _lmy_asin;
 }
 
static void _hoc_my_asin(void) {
  double _r;
   _r =  my_asin (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double my_sin (  double _larg ) {
   double _lmy_sin;
 
/*VERBATIM*/
	double ret;
	ret=sin(*getarg(1));
 _lmy_sin = ret ;
   
return _lmy_sin;
 }
 
static void _hoc_my_sin(void) {
  double _r;
   _r =  my_sin (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double exp_i (  double _larg ) {
   double _lexp_i;
 
/*VERBATIM*/
	double ret=1.0;
	double euler=exp(1.0);
	int n=0;
	for (n=0;n<(*getarg(1));n++) {
		ret*=euler;
	}
 _lexp_i = ret ;
   
return _lexp_i;
 }
 
static void _hoc_exp_i(void) {
  double _r;
   _r =  exp_i (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double mytime (  ) {
   double _lmytime;
 
/*VERBATIM*/
	double ret = 0.0;
	ret = (double)time(0)/3600.0;
 _lmytime = ret ;
   
return _lmytime;
 }
 
static void _hoc_mytime(void) {
  double _r;
   _r =  mytime (  );
 hoc_retpushx(_r);
}

static void initmodel() {
  int _i; double _save;_ninits++;
{

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
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{
} return _current;
}

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
}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/dustinmiao/Documents/school/rise/RISE-gaba-neurogen/kolbaev et. al test model/Models_Kolbaev et al Scientific Reports Revision/asin.mod";
static const char* nmodl_file_text = 
  "VERBATIM\n"
  "#include <math.h>\n"
  "#include <time.h>\n"
  "ENDVERBATIM\n"
  "\n"
  "FUNCTION my_asin(arg) {\n"
  "	VERBATIM\n"
  "	double ret;\n"
  "	ret=asin(*getarg(1));\n"
  "	ENDVERBATIM\n"
  "	my_asin=ret\n"
  "}\n"
  "\n"
  "FUNCTION my_sin(arg) {\n"
  "	VERBATIM\n"
  "	double ret;\n"
  "	ret=sin(*getarg(1));\n"
  "	ENDVERBATIM\n"
  "	my_sin=ret\n"
  "}\n"
  "\n"
  "FUNCTION exp_i(arg) {\n"
  "	VERBATIM\n"
  "	double ret=1.0;\n"
  "	double euler=exp(1.0);\n"
  "	int n=0;\n"
  "	for (n=0;n<(*getarg(1));n++) {\n"
  "		ret*=euler;\n"
  "	}\n"
  "	ENDVERBATIM\n"
  "	exp_i=ret\n"
  "}\n"
  "\n"
  "FUNCTION mytime() {\n"
  "	VERBATIM\n"
  "	double ret = 0.0;\n"
  "	ret = (double)time(0)/3600.0;\n"
  "	ENDVERBATIM\n"
  "	mytime = ret\n"
  "}\n"
  "\n"
  ;
#endif
