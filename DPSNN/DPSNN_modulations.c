// DPSNN_modulations.c
// Distributed Plastic Spiking Neural Network, Simulation Engine
// DPSNN_*.*
// AUTHOR: Pier Stanislao Paolucci (Roma, Italy, 2011),
// AUTHOR: Elena Pastorelli (2013-...)
// AUTHOR: Francesco Simula (2018-...)
// AUTHOR: ...
// AUTHOR: plus other members of INFN Lab, Roma, Italy

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/sysinfo.h>
#include <unistd.h>

#include "DPSNN_debug.h"
#include "DPSNN_random.h"
#include "DPSNN_neuron.h"
#include "DPSNN_messagePassing.h"
#include "DPSNN_spike.h"
#include "DPSNN_localNet.h"

void localNetClass::modulations_ms() {
  lnp_dynPar.fatigueDynamic.additionalFatigue=0.0;
  lnp_dynPar.inputCurrentDynamic.multiplierInhCurrent=1.0;
  lnp_dynPar.inputCurrentDynamic.multiplierExcCurrent=1.0;
 
  if((thisSimTimeStep_ms >= lnp_par.startModulations_ms) &&
     (thisSimTimeStep_ms <= lnp_par.stopModulations_ms))
    {
        fatigueModulation_ms();
	inhInputCurrentModulation_ms();
	excInputCurrentModulation_ms();
    };
};

void localNetClass::inhInputCurrentModulation_ms() {
  lnp_dynPar.inputCurrentDynamic.multiplierInhCurrent=1.0;

  if((thisSimTimeStep_ms >= lnp_par.startModulations_ms) &&
     (thisSimTimeStep_ms <= lnp_par.stopModulations_ms))
    {
      double effectiveTimeFromModulationStart_ms;
      double totalModulationPeriod_ms;
      double modulusOfEffectiveTime_ms;
      double transitionDuration_ms;
      double stationaryDuration_ms;
      double cosTime_ms;
      transitionDuration_ms = lnp_par.modulations.multiplierInhCurrentTransition_ms;
      stationaryDuration_ms = lnp_par.modulations.multiplierInhCurrentStationary_ms;
      
      totalModulationPeriod_ms =
	2 * transitionDuration_ms +
	2 * stationaryDuration_ms;
      
      effectiveTimeFromModulationStart_ms =
	(double)thisSimTimeStep_ms -
	lnp_par.startModulations_ms;

      modulusOfEffectiveTime_ms =
	fmod(effectiveTimeFromModulationStart_ms, totalModulationPeriod_ms);
      
      //first transition period
      if(modulusOfEffectiveTime_ms < transitionDuration_ms)
        cosTime_ms = modulusOfEffectiveTime_ms;
      else
        //first stationary period, at other fatigue value
	if (modulusOfEffectiveTime_ms <
	     transitionDuration_ms + stationaryDuration_ms)
	     cosTime_ms = transitionDuration_ms;
	else
	  //transition back period
	  if(modulusOfEffectiveTime_ms <
	     2 * transitionDuration_ms + stationaryDuration_ms)
	    cosTime_ms = modulusOfEffectiveTime_ms - stationaryDuration_ms;
	  else
	    //second stationary period, on original value
	    if(modulusOfEffectiveTime_ms <
	       2 * transitionDuration_ms + 2 *stationaryDuration_ms)
	      cosTime_ms = 2 * transitionDuration_ms;
	      else {
		printf(
		  "ERROR: modulations() should never assume modulusOfEffectiveTime_ms =%f\n",
		  modulusOfEffectiveTime_ms);
		fflush(stdout); exit(0);
	      };
     
       	  lnp_dynPar.inputCurrentDynamic.multiplierInhCurrent = 1.0 +
	    (lnp_par.modulations.multiplierInhCurrentExtreme - 1.0) * 0.5*
	     ( 1.0 - cos( (2 * 3.1415926536 * cosTime_ms)
			  / (2 * transitionDuration_ms)));
    };

 DPSNNverboseStart(true,1,0);
    if(lnp_par.loc_h==0&& (thisSimTimeStep_ms%100==0)) {
           printf("localNet_sim() -  modulations_ms() at %d ms: multiplierInhInputCurrent=%f, multiplierExtreme=%f, Transition_ms=%f Stationary_ms=%f\n",
		  thisSimTimeStep_ms, 
		  lnp_dynPar.inputCurrentDynamic.multiplierInhCurrent,
		  lnp_par.modulations.multiplierInhCurrentExtreme,
		  lnp_par.modulations.multiplierInhCurrentTransition_ms,
                  lnp_par.modulations.multiplierInhCurrentStationary_ms );
           fflush(stdout);
        }
   DPSNNverboseEnd();     
};

void localNetClass::excInputCurrentModulation_ms() {
  lnp_dynPar.inputCurrentDynamic.multiplierExcCurrent=1.0;

  if((thisSimTimeStep_ms >= lnp_par.startModulations_ms) &&
     (thisSimTimeStep_ms <= lnp_par.stopModulations_ms))
    {
      double effectiveTimeFromModulationStart_ms;
      double totalModulationPeriod_ms;
      double modulusOfEffectiveTime_ms;
      double transitionDuration_ms;
      double stationaryDuration_ms;
      double cosTime_ms;
      transitionDuration_ms = lnp_par.modulations.multiplierExcCurrentTransition_ms;
      stationaryDuration_ms = lnp_par.modulations.multiplierExcCurrentStationary_ms;
      
      totalModulationPeriod_ms =
	2 * transitionDuration_ms +
	2 * stationaryDuration_ms;
      
      effectiveTimeFromModulationStart_ms =
	(double)thisSimTimeStep_ms -
	lnp_par.startModulations_ms;

      modulusOfEffectiveTime_ms =
	fmod(effectiveTimeFromModulationStart_ms, totalModulationPeriod_ms);
      
      //first transition period
      if(modulusOfEffectiveTime_ms < transitionDuration_ms)
        cosTime_ms = modulusOfEffectiveTime_ms;
      else
        //first stationary period, at other fatigue value
	if (modulusOfEffectiveTime_ms <
	     transitionDuration_ms + stationaryDuration_ms)
	     cosTime_ms = transitionDuration_ms;
	else
	  //transition back period
	  if(modulusOfEffectiveTime_ms <
	     2 * transitionDuration_ms + stationaryDuration_ms)
	    cosTime_ms = modulusOfEffectiveTime_ms - stationaryDuration_ms;
	  else
	    //second stationary period, on original value
	    if(modulusOfEffectiveTime_ms <
	       2 * transitionDuration_ms + 2 *stationaryDuration_ms)
	      cosTime_ms = 2 * transitionDuration_ms;
	      else {
		printf(
		  "ERROR: modulations() should never assume modulusOfEffectiveTime_ms =%f\n",
		  modulusOfEffectiveTime_ms);
		fflush(stdout); exit(0);
	      };
     
       	  lnp_dynPar.inputCurrentDynamic.multiplierExcCurrent = 1.0 +
	    (lnp_par.modulations.multiplierExcCurrentExtreme - 1.0) * 0.5*
	     ( 1.0 - cos( (2 * 3.1415926536 * cosTime_ms)
			  / (2 * transitionDuration_ms)));
    };

 DPSNNverboseStart(true,1,0);
    if(lnp_par.loc_h==0&& (thisSimTimeStep_ms%100==0)) {
           printf("localNet_sim() -  modulations_ms() at %d ms: multiplierExcInputCurrent=%f, multiplierExtreme=%f, Transition_ms=%f Stationary_ms=%f\n",
		  thisSimTimeStep_ms, 
		  lnp_dynPar.inputCurrentDynamic.multiplierExcCurrent,
		  lnp_par.modulations.multiplierExcCurrentExtreme,
		  lnp_par.modulations.multiplierExcCurrentTransition_ms,
                  lnp_par.modulations.multiplierExcCurrentStationary_ms );
           fflush(stdout);
        }
   DPSNNverboseEnd();
};

void localNetClass::fatigueModulation_ms() {
  lnp_dynPar.fatigueDynamic.additionalFatigue=0.0;

  if((thisSimTimeStep_ms >= lnp_par.startModulations_ms) &&
     (thisSimTimeStep_ms <= lnp_par.stopModulations_ms))
    {
      double effectiveTimeFromModulationStart_ms;
      double totalModulationPeriod_ms;
      double modulusOfEffectiveTime_ms;
      double transitionDuration_ms;
      double stationaryDuration_ms;
      double cosTime_ms;
      transitionDuration_ms = lnp_par.modulations.additionalFatigueTransition_ms;
      stationaryDuration_ms = lnp_par.modulations.additionalFatigueStationary_ms;
      
      totalModulationPeriod_ms =
	2 * transitionDuration_ms +
	2 * stationaryDuration_ms;
      
      effectiveTimeFromModulationStart_ms =
	(double)thisSimTimeStep_ms -
	lnp_par.startModulations_ms;

      modulusOfEffectiveTime_ms =
	fmod(effectiveTimeFromModulationStart_ms, totalModulationPeriod_ms);
      
      //first transition period
      if(modulusOfEffectiveTime_ms < transitionDuration_ms)
        cosTime_ms = modulusOfEffectiveTime_ms;
      else
        //first stationary period, at other fatigue value
	if (modulusOfEffectiveTime_ms <
	     transitionDuration_ms + stationaryDuration_ms)
	     cosTime_ms = transitionDuration_ms;
	else
	  //transition back period
	  if(modulusOfEffectiveTime_ms <
	     2 * transitionDuration_ms + stationaryDuration_ms)
	    cosTime_ms = modulusOfEffectiveTime_ms - stationaryDuration_ms;
	  else
	    //second stationary period, on original value
	    if(modulusOfEffectiveTime_ms <
	       2 * transitionDuration_ms + 2 *stationaryDuration_ms)
	      cosTime_ms = 2 * transitionDuration_ms;
	      else {
		printf(
		  "ERROR: modulations() should never assume modulusOfEffectiveTime_ms =%f\n",
		  modulusOfEffectiveTime_ms);
		fflush(stdout); exit(0);
	      };
     
       	  lnp_dynPar.fatigueDynamic.additionalFatigue =
             lnp_par.modulations.additionalFatigueAmplitude * 0.5*
	     ( 1.0 - cos( (2 * 3.1415926536 * cosTime_ms)
			  / (2 * transitionDuration_ms)));
    }

 DPSNNverboseStart(true,1,0);
    if(lnp_par.loc_h==0&& (thisSimTimeStep_ms%100==0)) {
           printf("localNet_sim() -  modulations_ms() at %d ms: additionalFatigue=%f, fatigueAmplitude=%f, Transition_ms=%f Stationary_ms=%f\n",
		  thisSimTimeStep_ms, 
		  lnp_dynPar.fatigueDynamic.additionalFatigue,
		  lnp_par.modulations.additionalFatigueAmplitude,
		  lnp_par.modulations.additionalFatigueTransition_ms,
                  lnp_par.modulations.additionalFatigueStationary_ms );
           fflush(stdout);
        }
   DPSNNverboseEnd();     
};
