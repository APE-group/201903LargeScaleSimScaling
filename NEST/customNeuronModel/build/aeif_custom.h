
/*
*  aeif_custom.h
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
*  2018-06-20 14:53:30.658439
*/
#ifndef AEIF_CUSTOM
#define AEIF_CUSTOM

#include "config.h"


#ifdef HAVE_GSL

// External includes:
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// forwards the declaration of the function
/**
 * Function computing right-hand side of ODE for GSL solver.
 * @note Must be declared here so we can befriend it in class.
 * @note Must have C-linkage for passing to GSL. Internally, it is
 *       a first-class C++ function, but cannot be a member function
 *       because of the C-linkage.
 * @note No point in declaring it inline, since it is called
 *       through a function pointer.
 * @param void* Pointer to model neuron instance.
 */
extern "C" inline int aeif_custom_dynamics( double, const double y[], double f[], void* pnode );


// Includes from nestkernel:
#include "archiving_node.h"
#include "connection.h"
#include "event.h"
#include "nest_types.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"


// Includes from sli:
#include "dictdatum.h"

/* BeginDocumentation
  Name: aeif_custom.

  Description:  
    
  Name: aeif_cond_alpha -  Conductance based exponential integrate-and-fire neuron
                           model according to Brette and Gerstner (2005).

  Description:
  aeif_cond_alpha is the adaptive exponential integrate and fire neuron according
  to Brette and Gerstner (2005).
  Synaptic conductances are modelled as alpha-functions.

  This implementation uses the embedded 4th order Runge-Kutta-Fehlberg solver with
  adaptive step size to integrate the differential equation.

  The membrane potential is given by the following differential equation:
  C dV/dt= -g_L(V-E_L)+g_L*Delta_T*exp((V-V_T)/Delta_T)-g_e(t)(V-E_e)
                                                       -g_i(t)(V-E_i)-w +I_e

  and

  tau_w * dw/dt= a(V-E_L) -W

  Author: Marc-Oliver Gewaltig

  Sends: SpikeEvent

  Receives: SpikeEvent, CurrentEvent, DataLoggingRequest

  References: Brette R and Gerstner W (2005) Adaptive Exponential
              Integrate-and-Fire Model as an Effective Description of Neuronal
              Activity. J Neurophysiol 94:3637-3642

  SeeAlso: iaf_cond_alpha, aeif_cond_exp



  Parameters:
  The following parameters can be set in the status dictionary.
  C_m [pF]  membrane parameters
281.0pF        Membrane Capacitance in pF
  t_ref [ms] 0.0ms          Refractory period in ms
  V_reset [mV] -60.0mV      Reset Potential in mV
  g_L [nS] 30.0nS           Leak Conductance in nS
  E_L [mV] -70.6mV          Leak reversal Potential (aka resting potential) in mV
  I_e [pA]  Constant Current in pA
  b [pA]  spike adaptation parameters
a nS = 4nS                Subthreshold adaptation
0.5pA             pike-triggered adaptation
Delta_T mV = 2.0mV        Slope factor
  tau_w [ms] Delta_T mV = 2.0mV        Slope factor
144.0ms        Adaptation time constant
  V_th [mV] -50.4mV         Threshold Potential in mV
  V_peak [mV] 0mV           Spike detection threshold
  I_stim [pA]  Input current injected by CurrentEvent.
 This variable is used to transport the current applied into the
 _dynamics function computing the derivative of the state vector.
  

  Dynamic state variables:
  

  Initial values:
  V_m [mV]  Membrane potential
  w [pA]  Spike-adaptation current
  

  References: Empty

  Sends: nest::SpikeEvent

  Receives: Spike, Current, DataLoggingRequest
*/
class aeif_custom : public nest::Archiving_Node{
public:
  /**
  * The constructor is only used to create the model prototype in the model manager.
  */
  aeif_custom();

  /**
  * The copy constructor is used to create model copies and instances of the model.
  * @node The copy constructor needs to initialize the parameters and the state.
  *       Initialization of buffers and interal variables is deferred to
  *       @c init_buffers_() and @c calibrate().
  */
  aeif_custom(const aeif_custom &);

  /**
  * Releases resources.
  */
  ~aeif_custom();

  /**
   * Import sets of overloaded virtual functions.
   * @see Technical Issues / Virtual Functions: Overriding, Overloading, and
   * Hiding
   */
  using nest::Node::handles_test_event;
  using nest::Node::handle;

  /**
  * Used to validate that we can send nest::SpikeEvent to desired target:port.
  */
  nest::port send_test_event(nest::Node& target, nest::rport receptor_type, nest::synindex, bool);

  /**
  * @defgroup mynest_handle Functions handling incoming events.
  * We tell nest that we can handle incoming events of various types by
  * defining @c handle() and @c connect_sender() for the given event.
  * @{
  */
  void handle(nest::SpikeEvent &);        //! accept spikes
  void handle(nest::CurrentEvent &);      //! accept input current
  void handle(nest::DataLoggingRequest &);//! allow recording with multimeter

  nest::port handles_test_event(nest::SpikeEvent&, nest::port);
  nest::port handles_test_event(nest::CurrentEvent&, nest::port);
  nest::port handles_test_event(nest::DataLoggingRequest&, nest::port);
  /** @} */

  // SLI communication functions:
  void get_status(DictionaryDatum &) const;
  void set_status(const DictionaryDatum &);

private:
  //! Reset parameters and state of neuron.

  //! Reset state of neuron.
  void init_state_(const Node& proto);

  //! Reset internal buffers of neuron.
  void init_buffers_();

  //! Initialize auxiliary quantities, leave parameters and state untouched.
  void calibrate();

  //! Take neuron through given time interval
  void update(nest::Time const &, const long, const long);

  // The next two classes need to be friends to access the State_ class/member
  friend class nest::RecordablesMap<aeif_custom>;
  friend class nest::UniversalDataLogger<aeif_custom>;

  /**
  * Free parameters of the neuron.
  *
  *
  *
  * These are the parameters that can be set by the user through @c SetStatus.
  * They are initialized from the model prototype when the node is created.
  * Parameters do not change during calls to @c update() and are not reset by
  * @c ResetNetwork.
  *
  * @note Parameters_ need neither copy constructor nor @c operator=(), since
  *       all its members are copied properly by the default copy constructor
  *       and assignment operator. Important:
  *       - If Parameters_ contained @c Time members, you need to define the
  *         assignment operator to recalibrate all members of type @c Time . You
  *         may also want to define the assignment operator.
  *       - If Parameters_ contained members that cannot copy themselves, such
  *         as C-style arrays, you need to define the copy constructor and
  *         assignment operator to copy those members.
  */
  struct Parameters_{
        
        

    //!  membrane parameters
    //! 281.0pF        Membrane Capacitance in pF
    double C_m;

    //! 0.0ms          Refractory period in ms
    double t_ref;

    //! -60.0mV      Reset Potential in mV
    double V_reset;

    //! 30.0nS           Leak Conductance in nS
    double g_L;

    //! -70.6mV          Leak reversal Potential (aka resting potential) in mV
    double E_L;

    //!  Constant Current in pA
    double I_e;

    //!  spike adaptation parameters
    //! a nS = 4nS                Subthreshold adaptation
    //! 0.5pA             pike-triggered adaptation
    //! Delta_T mV = 2.0mV        Slope factor
    double b;

    //! Delta_T mV = 2.0mV        Slope factor
    //! 144.0ms        Adaptation time constant
    double tau_w;

    //! -50.4mV         Threshold Potential in mV
    double V_th;

    //! 0mV           Spike detection threshold
    double V_peak;

    //!  Input current injected by CurrentEvent.
    //!  This variable is used to transport the current applied into the
    //!  _dynamics function computing the derivative of the state vector.
    double I_stim;

    double __gsl_error_tol;
    /** Initialize parameters to their default values. */
    Parameters_();
  };

  /**
  * Dynamic state of the neuron.
  *
  *
  *
  * These are the state variables that are advanced in time by calls to
  * @c update(). In many models, some or all of them can be set by the user
  * through @c SetStatus. The state variables are initialized from the model
  * prototype when the node is created. State variables are reset by @c ResetNetwork.
  *
  * @note State_ need neither copy constructor nor @c operator=(), since
  *       all its members are copied properly by the default copy constructor
  *       and assignment operator. Important:
  *       - If State_ contained @c Time members, you need to define the
  *         assignment operator to recalibrate all members of type @c Time . You
  *         may also want to define the assignment operator.
  *       - If State_ contained members that cannot copy themselves, such
  *         as C-style arrays, you need to define the copy constructor and
  *         assignment operator to copy those members.
  */
  struct State_{
    //! Symbolic indices to the elements of the state vector y
    enum StateVecElems{
    //  Membrane potential
      V_m,
      //  Spike-adaptation current
      w,
      STATE_VEC_SIZE
    };
    //! state vector, must be C-array for GSL solver
    double ode_state[STATE_VEC_SIZE];    

    State_();
  };

  /**
  * Internal variables of the neuron.
  *
  *
  *
  * These variables must be initialized by @c calibrate, which is called before
  * the first call to @c update() upon each call to @c Simulate.
  * @node Variables_ needs neither constructor, copy constructor or assignment operator,
  *       since it is initialized by @c calibrate(). If Variables_ has members that
  *       cannot destroy themselves, Variables_ will need a destructor.
  */
  struct Variables_ {    

    //!  refractory time in steps
    //!  counts number of tick during the refractory period
    long RefractoryCounts;
        

    //!  counts number of tick during the refractory period
    long r;
    
  };

  /**
    * Buffers of the neuron.
    * Ususally buffers for incoming spikes and data logged for analog recorders.
    * Buffers must be initialized by @c init_buffers_(), which is called before
    * @c calibrate() on the first call to @c Simulate after the start of NEST,
    * ResetKernel or ResetNetwork.
    * @node Buffers_ needs neither constructor, copy constructor or assignment operator,
    *       since it is initialized by @c init_nodes_(). If Buffers_ has members that
    *       cannot destroy themselves, Buffers_ will need a destructor.
    */
  struct Buffers_ {
    Buffers_(aeif_custom &);
    Buffers_(const Buffers_ &, aeif_custom &);

    /** Logger for all analog data */
    nest::UniversalDataLogger<aeif_custom> logger_;
    
    inline nest::RingBuffer& get_spikesInh() {return spikesInh;}
    //!< Buffer incoming nSs through delay, as sum
    nest::RingBuffer spikesInh;
    double spikesInh_grid_sum_;
    
    inline nest::RingBuffer& get_spikesExc() {return spikesExc;}
    //!< Buffer incoming nSs through delay, as sum
    nest::RingBuffer spikesExc;
    double spikesExc_grid_sum_;
    
    //!< Buffer incoming pAs through delay, as sum
    nest::RingBuffer currents;
    inline nest::RingBuffer& get_currents() {return currents;}
    double currents_grid_sum_;
    /** GSL ODE stuff */
    gsl_odeiv_step* __s;    //!< stepping function
    gsl_odeiv_control* __c; //!< adaptive stepsize control function
    gsl_odeiv_evolve* __e;  //!< evolution function
    gsl_odeiv_system __sys; //!< struct describing system

    // IntergrationStep_ should be reset with the neuron on ResetNetwork,
    // but remain unchanged during calibration. Since it is initialized with
    // step_, and the resolution cannot change after nodes have been created,
    // it is safe to place both here.
    double __step;             //!< step size in ms
    double __integration_step; //!< current integration time step, updated by GSL
    };
  inline double get_V_m() const {
    return S_.ode_state[State_::V_m];
  }
  inline void set_V_m(const double __v) {
    S_.ode_state[State_::V_m] = __v;
  }

  inline double get_w() const {
    return S_.ode_state[State_::w];
  }
  inline void set_w(const double __v) {
    S_.ode_state[State_::w] = __v;
  }

  inline double get_C_m() const {
    return P_.C_m;
  }
  inline void set_C_m(const double __v) {
    P_.C_m = __v;
  }

  inline double get_t_ref() const {
    return P_.t_ref;
  }
  inline void set_t_ref(const double __v) {
    P_.t_ref = __v;
  }

  inline double get_V_reset() const {
    return P_.V_reset;
  }
  inline void set_V_reset(const double __v) {
    P_.V_reset = __v;
  }

  inline double get_g_L() const {
    return P_.g_L;
  }
  inline void set_g_L(const double __v) {
    P_.g_L = __v;
  }

  inline double get_E_L() const {
    return P_.E_L;
  }
  inline void set_E_L(const double __v) {
    P_.E_L = __v;
  }

  inline double get_I_e() const {
    return P_.I_e;
  }
  inline void set_I_e(const double __v) {
    P_.I_e = __v;
  }

  inline double get_b() const {
    return P_.b;
  }
  inline void set_b(const double __v) {
    P_.b = __v;
  }

  inline double get_tau_w() const {
    return P_.tau_w;
  }
  inline void set_tau_w(const double __v) {
    P_.tau_w = __v;
  }

  inline double get_V_th() const {
    return P_.V_th;
  }
  inline void set_V_th(const double __v) {
    P_.V_th = __v;
  }

  inline double get_V_peak() const {
    return P_.V_peak;
  }
  inline void set_V_peak(const double __v) {
    P_.V_peak = __v;
  }

  inline double get_I_stim() const {
    return P_.I_stim;
  }
  inline void set_I_stim(const double __v) {
    P_.I_stim = __v;
  }

  inline long get_RefractoryCounts() const {
    return V_.RefractoryCounts;
  }
  inline void set_RefractoryCounts(const long __v) {
    V_.RefractoryCounts = __v;
  }

  inline long get_r() const {
    return V_.r;
  }
  inline void set_r(const long __v) {
    V_.r = __v;
  }


  
  inline nest::RingBuffer& get_spikesInh() {return B_.get_spikesInh();};
  
  inline nest::RingBuffer& get_spikesExc() {return B_.get_spikesExc();};
  
  inline nest::RingBuffer& get_currents() {return B_.get_currents();};
  

  // Generate function header
  
  /**
  * @defgroup pif_members Member variables of neuron model.
  * Each model neuron should have precisely the following four data members,
  * which are one instance each of the parameters, state, buffers and variables
  * structures. Experience indicates that the state and variables member should
  * be next to each other to achieve good efficiency (caching).
  * @note Devices require one additional data member, an instance of the @c Device
  *       child class they belong to.
  * @{
  */
  Parameters_ P_;  //!< Free parameters.
  State_      S_;  //!< Dynamic state.
  Variables_  V_;  //!< Internal Variables
  Buffers_    B_;  //!< Buffers.

  //! Mapping of recordables names to access functions
  static nest::RecordablesMap<aeif_custom> recordablesMap_;

  friend int aeif_custom_dynamics( double, const double y[], double f[], void* pnode );
  
/** @} */
}; /* neuron aeif_custom */

inline nest::port aeif_custom::send_test_event(
    nest::Node& target, nest::rport receptor_type, nest::synindex, bool){
  // You should usually not change the code in this function.
  // It confirms that the target of connection @c c accepts @c nest::SpikeEvent on
  // the given @c receptor_type.
  nest::SpikeEvent e;
  e.set_sender(*this);
  return target.handles_test_event(e, receptor_type);
}

inline nest::port aeif_custom::handles_test_event(nest::SpikeEvent&, nest::port receptor_type){
  
    // You should usually not change the code in this function.
    // It confirms to the connection management system that we are able
    // to handle @c SpikeEvent on port 0. You need to extend the function
    // if you want to differentiate between input ports.
    if (receptor_type != 0)
      throw nest::UnknownReceptorType(receptor_type, get_name());
    return 0;
}



inline nest::port aeif_custom::handles_test_event(
    nest::CurrentEvent&, nest::port receptor_type){
  // You should usually not change the code in this function.
  // It confirms to the connection management system that we are able
  // to handle @c CurrentEvent on port 0. You need to extend the function
  // if you want to differentiate between input ports.
  if (receptor_type != 0)
  throw nest::UnknownReceptorType(receptor_type, get_name());
  return 0;
}

inline nest::port aeif_custom::handles_test_event(
    nest::DataLoggingRequest& dlr, nest::port receptor_type){
  // You should usually not change the code in this function.
  // It confirms to the connection management system that we are able
  // to handle @c DataLoggingRequest on port 0.
  // The function also tells the built-in UniversalDataLogger that this node
  // is recorded from and that it thus needs to collect data during simulation.
  if (receptor_type != 0)
  throw nest::UnknownReceptorType(receptor_type, get_name());

  return B_.logger_.connect_logging_device(dlr, recordablesMap_);
}

// TODO call get_status on used or internal components
inline void aeif_custom::get_status(DictionaryDatum &__d) const{  
  def<double>(__d, "C_m", get_C_m());
      
  def<double>(__d, "t_ref", get_t_ref());
      
  def<double>(__d, "V_reset", get_V_reset());
      
  def<double>(__d, "g_L", get_g_L());
      
  def<double>(__d, "E_L", get_E_L());
      
  def<double>(__d, "I_e", get_I_e());
      
  def<double>(__d, "b", get_b());
      
  def<double>(__d, "tau_w", get_tau_w());
      
  def<double>(__d, "V_th", get_V_th());
      
  def<double>(__d, "V_peak", get_V_peak());
      
  def<double>(__d, "I_stim", get_I_stim());
      
  def<double>(__d, "V_m", get_V_m());
      
  def<double>(__d, "w", get_w());
    

  (*__d)[nest::names::recordables] = recordablesMap_.get_list();
  
  def< double >(__d, nest::names::gsl_error_tol, P_.__gsl_error_tol);
  if ( P_.__gsl_error_tol <= 0. ){
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }
  

}

inline void aeif_custom::set_status(const DictionaryDatum &__d){

  double tmp_C_m = get_C_m();
  updateValue<double>(__d, "C_m", tmp_C_m);


  double tmp_t_ref = get_t_ref();
  updateValue<double>(__d, "t_ref", tmp_t_ref);


  double tmp_V_reset = get_V_reset();
  updateValue<double>(__d, "V_reset", tmp_V_reset);


  double tmp_g_L = get_g_L();
  updateValue<double>(__d, "g_L", tmp_g_L);


  double tmp_E_L = get_E_L();
  updateValue<double>(__d, "E_L", tmp_E_L);


  double tmp_I_e = get_I_e();
  updateValue<double>(__d, "I_e", tmp_I_e);


  double tmp_b = get_b();
  updateValue<double>(__d, "b", tmp_b);


  double tmp_tau_w = get_tau_w();
  updateValue<double>(__d, "tau_w", tmp_tau_w);


  double tmp_V_th = get_V_th();
  updateValue<double>(__d, "V_th", tmp_V_th);


  double tmp_V_peak = get_V_peak();
  updateValue<double>(__d, "V_peak", tmp_V_peak);


  double tmp_I_stim = get_I_stim();
  updateValue<double>(__d, "I_stim", tmp_I_stim);


  double tmp_V_m = get_V_m();
  updateValue<double>(__d, "V_m", tmp_V_m);

  

  double tmp_w = get_w();
  updateValue<double>(__d, "w", tmp_w);

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  Archiving_Node::set_status(__d);

  // if we get here, temporaries contain consistent set of properties


  set_C_m(tmp_C_m);



  set_t_ref(tmp_t_ref);



  set_V_reset(tmp_V_reset);



  set_g_L(tmp_g_L);



  set_E_L(tmp_E_L);



  set_I_e(tmp_I_e);



  set_b(tmp_b);



  set_tau_w(tmp_tau_w);



  set_V_th(tmp_V_th);



  set_V_peak(tmp_V_peak);



  set_I_stim(tmp_I_stim);



  set_V_m(tmp_V_m);



  set_w(tmp_w);


  
  updateValue< double >(__d, nest::names::gsl_error_tol, P_.__gsl_error_tol);
  if ( P_.__gsl_error_tol <= 0. ){
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }
  
};

#endif /* #ifndef AEIF_CUSTOM */
#endif /* HAVE GSL */