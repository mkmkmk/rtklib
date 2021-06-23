/*
 * NavObservable.h
 *
 *  Created on: Apr 7, 2020
 *      Author: M. Krej
 */

#ifndef GNSS_LIB_NAVOBSERVABLE_H_
#define GNSS_LIB_NAVOBSERVABLE_H_


//---------------------
/*
 * Struktura na pomiary
 */
typedef struct
{
    double  rx_tow;
    double  tow;
    double  carrier_Doppler_hz;
    double  carrier_phase_cyc;
    double  pseudorange_m;
    int prn;
    int valid;
} NavObservable;


#endif /* GNSS_LIB_NAVOBSERVABLE_H_ */
