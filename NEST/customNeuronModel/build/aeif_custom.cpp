/*
*  aeif_custom.cpp
*
*  This file is part of NEST.
*
*  Copyright (C) 2004 The NEST Initiative
*
*  NEST is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 2 of the License, or
*  (at your option) any later version.
*
*  NEST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
*
*  2018-06-20 14:53:30.749853
*/

// C++ includes:
#include <limits>

// Includes from libnestutil:
#include "numerics.h"

// Includes from nestkernel:
#include "exceptions.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"
#include "lockptrdatum.h"

#include "aeif_custom.h"


/* ----------------------------------------------------------------
* Recordables map
* ---------------------------------------------------------------- */
nest::RecordablesMap<aeif_custom> aeif_custom::recordablesMap_;

namespace nest
{
  // Override the create() method with one call to RecordablesMap::insert_()
  // for each quantity to be recorded.
  template <> void RecordablesMap<aeif_custom>::create(){
  // use standard names whereever you can for consistency!  

  insert_("V_m", &aeif_custom::get_V_m);

  insert_("w", &aeif_custom::get_w);
  }
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * Note: the implementation is empty. The initialization is of variables
 * is a part of the aeif_custom's constructor.
 * ---------------------------------------------------------------- */
aeif_custom::Parameters_::Parameters_(){}

aeif_custom::State_::State_(){}

/* ----------------------------------------------------------------
* Parameter and state extractions and manipulation functions
* ---------------------------------------------------------------- */

aeif_custom::Buffers_::Buffers_(aeif_custom &n):
  logger_(n), __s( 0 ), __c( 0 ), __e( 0 ){
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

aeif_custom::Buffers_::Buffers_(const Buffers_ &, aeif_custom &n):
  logger_(n), __s( 0 ), __c( 0 ), __e( 0 ){
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */
aeif_custom::aeif_custom():Archiving_Node(), P_(), S_(), B_(*this)
{
  recordablesMap_.create();
  // use a default `good` enough value for the absolute error.
  // it cab be adjusted via `SetStatus`
  P_.__gsl_error_tol = 1e-3;
  
  P_.C_m = 50.0*1.0; // as pF
  
  P_.t_ref = 2.0*1.0; // as ms
  
  P_.V_reset = 15.0*1.0; // as mV
  
  P_.g_L = 2.5*1.0; // as nS
  
  P_.E_L = 0.0*1.0; // as mV
  
  P_.I_e = 0*1.0; // as pA
  
  P_.b = 1.0*1.0; // as pA
  
  P_.tau_w = 1000.0*1.0; // as ms
  
  P_.V_th = 20.0*1.0; // as mV
  
  P_.V_peak = 20.0*1.0; // as mV
  
  P_.I_stim = 0*1.0; // as pA
  
  S_.ode_state[State_::V_m] = 15*1.0; // as mV
  
  S_.ode_state[State_::w] = 0*1.0; // as pA
}

aeif_custom::aeif_custom(const aeif_custom& __n):
  Archiving_Node(), P_(__n.P_), S_(__n.S_), B_(__n.B_, *this){
  P_.C_m = __n.P_.C_m;
  P_.t_ref = __n.P_.t_ref;
  P_.V_reset = __n.P_.V_reset;
  P_.g_L = __n.P_.g_L;
  P_.E_L = __n.P_.E_L;
  P_.I_e = __n.P_.I_e;
  P_.b = __n.P_.b;
  P_.tau_w = __n.P_.tau_w;
  P_.V_th = __n.P_.V_th;
  P_.V_peak = __n.P_.V_peak;
  P_.I_stim = __n.P_.I_stim;
  
  
  S_.ode_state[State_::V_m] = __n.S_.ode_state[State_::V_m];
  S_.ode_state[State_::w] = __n.S_.ode_state[State_::w];
  
  V_.RefractoryCounts = __n.V_.RefractoryCounts;
  V_.r = __n.V_.r;
  
}

aeif_custom::~aeif_custom(){ 
  // GSL structs may not have been allocated, so we need to protect destruction
  if (B_.__s)
    gsl_odeiv_step_free( B_.__s );
  if (B_.__c)
    gsl_odeiv_control_free( B_.__c );
  if (B_.__e)
    gsl_odeiv_evolve_free( B_.__e );
}

/* ----------------------------------------------------------------
* Node initialization functions
* ---------------------------------------------------------------- */

void aeif_custom::init_state_(const Node& proto){
  const aeif_custom& pr = downcast<aeif_custom>(proto);
  S_ = pr.S_;
}



extern "C" inline int aeif_custom_dynamics(double, const double ode_state[], double f[], void* pnode){
  typedef aeif_custom::State_ State_;
  // get access to node so we can almost work as in a member function
  assert( pnode );
  const aeif_custom& node = *( reinterpret_cast< aeif_custom* >( pnode ) );

  // ode_state[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.ode_state[].
  
  f[State_::V_m] = (-(node.get_g_L()) * (ode_state[State_::V_m] - node.get_E_L()) - ode_state[State_::w] + node.get_I_e() + node.get_I_stim() + node.B_.spikesExc_grid_sum_ + node.B_.spikesInh_grid_sum_) / node.get_C_m();
  f[State_::w] = -(ode_state[State_::w]) / node.get_tau_w();
  return GSL_SUCCESS;
}



void aeif_custom::init_buffers_(){
  get_spikesInh().clear(); //includes resize
  get_spikesExc().clear(); //includes resize
  get_currents().clear(); //includes resize
  
  B_.logger_.reset(); // includes resize
  Archiving_Node::clear_history();
  
  if ( B_.__s == 0 ){
    B_.__s = gsl_odeiv_step_alloc( gsl_odeiv_step_rkf45, 2 );
  } else {
    gsl_odeiv_step_reset( B_.__s );
  }

  if ( B_.__c == 0 ){
    B_.__c = gsl_odeiv_control_y_new( P_.__gsl_error_tol, 0.0 );
  } else {
    gsl_odeiv_control_init( B_.__c, P_.__gsl_error_tol, 0.0, 1.0, 0.0 );
  }

  if ( B_.__e == 0 ){
    B_.__e = gsl_odeiv_evolve_alloc( 2 );
  } else {
    gsl_odeiv_evolve_reset( B_.__e );
  }

  B_.__sys.function = aeif_custom_dynamics;
  B_.__sys.jacobian = NULL;
  B_.__sys.dimension = 2;
  B_.__sys.params = reinterpret_cast< void* >( this );
  B_.__step = nest::Time::get_resolution().get_ms();
  B_.__integration_step = nest::Time::get_resolution().get_ms();
}

void aeif_custom::calibrate(){
  B_.logger_.init();
  
  
  V_.RefractoryCounts =nest::Time(nest::Time::ms((double) P_.t_ref)).get_steps();
  
  
  V_.r =0;
}

/* ----------------------------------------------------------------
* Update and spike handling functions
* ---------------------------------------------------------------- */

/*
 *
 */
void aeif_custom::update(nest::Time const & origin,const long from, const long to){
  double __t = 0;

  for ( long lag = from ; lag < to ; ++lag ) {
    B_.spikesInh_grid_sum_ = get_spikesInh().get_value(lag);
    B_.spikesExc_grid_sum_ = get_spikesExc().get_value(lag);
    B_.currents_grid_sum_ = get_currents().get_value(lag);
      
    

    __t = 0;
    // numerical integration with adaptive step size control:
    // ------------------------------------------------------
    // gsl_odeiv_evolve_apply performs only a single numerical
    // integration step, starting from t and bounded by step;
    // the while-loop ensures integration over the whole simulation
    // step (0, step] if more than one integration step is needed due
    // to a small integration step size;
    // note that (t+IntegrationStep > step) leads to integration over
    // (t, step] and afterwards setting t to step, but it does not
    // enforce setting IntegrationStep to step-t; this is of advantage
    // for a consistent and efficient integration across subsequent
    // simulation intervals
    while ( __t < B_.__step )
    {
      const int status = gsl_odeiv_evolve_apply(B_.__e,
                                                B_.__c,
                                                B_.__s,
                                                &B_.__sys,              // system of ODE
                                                &__t,                   // from t
                                                B_.__step,              // to t <= step
                                                &B_.__integration_step, // integration step size
                                                S_.ode_state);          // neuronal state

      if ( status != GSL_SUCCESS ) {
        throw nest::GSLSolverFailure( get_name(), status );
      }
    }
    



    if (V_.r>0) {
      

      V_.r = V_.r-1;
      

      S_.ode_state[State_::V_m] = P_.V_reset;
    }else if(S_.ode_state[State_::V_m]>=P_.V_peak) {
      

      V_.r = V_.RefractoryCounts;
      

      S_.ode_state[State_::V_m] = P_.V_reset;
      

      S_.ode_state[State_::w] += P_.b;
      
      set_spiketime(nest::Time::step(origin.get_steps()+lag+1));
      nest::SpikeEvent se;
      nest::kernel().event_delivery_manager.send(*this, se, lag);
    } /* if end */


    // voltage logging
    B_.logger_.record_data(origin.get_steps()+lag);
  }

}

// Do not move this function as inline to h-file. It depends on
// universal_data_logger_impl.h being included here.
void aeif_custom::handle(nest::DataLoggingRequest& e){
  B_.logger_.handle(e);
}


void aeif_custom::handle(nest::SpikeEvent &e){
  assert(e.get_delay() > 0);
  
  const double weight = e.get_weight();
  const double multiplicity = e.get_multiplicity();
  
  if ( weight < 0.0 ){ // inhibitory
    get_spikesInh().
        add_value(e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
                       // ensure conductance is positive 
                       -1 *  weight * multiplicity );
  }
  if ( weight >= 0.0 ){ // excitatory
    get_spikesExc().
        add_value(e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
                       weight * multiplicity );
  }
}

void aeif_custom::handle(nest::CurrentEvent& e){
  assert(e.get_delay() > 0);

  const double current=e.get_current();
  const double weight=e.get_weight();

  // add weighted current; HEP 2002-10-04
  get_currents().add_value(
               e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
               weight * current );
  
}
