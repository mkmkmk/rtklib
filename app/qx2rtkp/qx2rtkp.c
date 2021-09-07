#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <stddef.h>
#include <stdarg.h>
#include "rtklib.h"

#include <fcntl.h>
#include <errno.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>

#include "NavObservable.h"

#define SAVE_TRACE (1)

// temp powrót do TOW w Q.30
#define OLD_OBSO_Q30_TOW (0)


#define MAX_CHAN (256)

#define Q_32_dbl ((double)(1LL << 32))
#define Q_10_dbl ((double)(1LL << 10))
#define Q_30_dbl ((double)(1LL << 30))


static obs_t obss={0};          /* observation data */
static nav_t navs={0};          /* navigation data */
static sta_t stas[MAXRCV];      /* station infomation */
static int nepoch=0;            /* number of observation epochs */
static char proc_rov [64]="";   /* rover for current processing */
static char proc_base[64]="";   /* base station for current processing */


// #define REF_POS {4849202.3940,-360328.9929,4114913.1862}

// Kobyłka ref v2 2019-06-03
// #define REF_POS {3642332.408, 1411096.212, 5025380.626}

// Kobyłka ref v3 2021-05-01
#define REF_POS {3642325.000, 1411093.342, 5025370.336}

#define ELEVATION_MASK (6)
#define NUM_BANDS (1)

const prcopt_t prcopt_qx =    { /* defaults processing options */
    //PMODE_STATIC,0,NUM_BANDS,SYS_GPS,   /* mode,soltype,nf,navsys */
    PMODE_PPP_STATIC,0,NUM_BANDS,SYS_GPS,   /* mode,soltype,nf,navsys */
    ELEVATION_MASK*D2R,{{0,0}},            /* elmin,snrmask */
    0,1,1,1,                    /* sateph,modear,glomodear,bdsmodear */
    5,0,10,                     /* glomodear,maxout,minlock,minfix */
    IONOOPT_BRDC,TROPOPT_SAAS,0,0,     /* estion,esttrop,dynamics,tidecorr */
    //IONOOPT_IFLC,TROPOPT_SAAS,0,0,
    1,0,1,0,0,                  /* niter,codesmooth,intpref,sbascorr,sbassatsel */
    0,0,                        /* rovpos,refpos */
    {100.0,100.0},              /* eratio[] */
    {100.0,0.003,0.003,0.0,1.0}, /* err[] */
    {30.0,0.03,0.3},            /* std[] */
    {1E-4,1E-3,1E-4,1,1E-1}, /* prn[] co to??*/
    5E-12,                      /* sclkstab */
    {5.0,0.9999,0.20},          /* thresar */
    0.0,0.0,0.05,               /* elmaskar,almaskhold,thresslip */
    30.0,30.0,30.0,             /* maxtdif,maxinno,maxgdop */
    {0},{6378137},REF_POS,                /* baseline,ru,rb */
    {"*","*"},                    /* anttype */
    {{0}},{{0}},{0}             /* antdel,pcv,exsats */
};




extern int showmsg(char *format, ...)
{
    va_list arg;
    va_start(arg,format); vfprintf(stderr,format,arg); va_end(arg);
    fprintf(stderr,"\r");
    return 0;
}
extern void settspan(gtime_t ts, gtime_t te) {}
extern void settime(gtime_t time) {}

/* show message and check break ----------------------------------------------*/
static int checkbrk(const char *format, ...)
{
    va_list arg;
    char buff[1024],*p=buff;
    if (!*format) return showmsg("");
    va_start(arg,format);
    p+=vsprintf(p,format,arg);
    va_end(arg);
    if (*proc_rov&&*proc_base) sprintf(p," (%s-%s)",proc_rov,proc_base);
    else if (*proc_rov ) sprintf(p," (%s)",proc_rov );
    else if (*proc_base) sprintf(p," (%s)",proc_base);
    return showmsg(buff);
}
/* read obs and nav data -----------------------------------------------------*/
static int readobsnav(gtime_t ts, gtime_t te, double ti, char **infile,
                      const int *index, int n, const prcopt_t *prcopt,
                      obs_t *obs, nav_t *nav, sta_t *sta)
{
    int i,j,ind=0,nobs=0,rcv=1;

    trace(3,"readobsnav: ts=%s n=%d\n",time_str(ts,0),n);

    obs->data=NULL; obs->n =obs->nmax =0;
    nav->eph =NULL; nav->n =nav->nmax =0;
    nav->geph=NULL; nav->ng=nav->ngmax=0;
    nav->seph=NULL; nav->ns=nav->nsmax=0;
    nepoch=0;

    for (i=0;i<n;i++) {
        if (checkbrk("")) return 0;

        if (index[i]!=ind) {
            if (obs->n>nobs) rcv++;
            ind=index[i]; nobs=obs->n;
        }
        /* read rinex obs and nav file */
        if (readrnxt(infile[i],rcv,ts,te,ti,prcopt->rnxopt[rcv<=1?0:1],obs,nav,
                     rcv<=2?sta+rcv-1:NULL)<0) {
            checkbrk("error : insufficient memory");
            trace(1,"insufficient memory\n");
            return 0;
        }
    }
    /*if (obs->n<=0) {
        checkbrk("error : no obs data");
        trace(1,"no obs data\n");
        return 0;
    }*/
    if (nav->n<=0&&nav->ng<=0&&nav->ns<=0) {
        checkbrk("error : no nav data");
        trace(1,"no nav data\n");
        return 0;
    }
    /* sort observation data */
    //nepoch=sortobs(obs);

    /* delete duplicated ephemeris */
    uniqnav(nav);

    /* set time span for progress display */
    if (ts.time==0||te.time==0) {
        for (i=0;   i<obs->n;i++) if (obs->data[i].rcv==1) break;
        for (j=obs->n-1;j>=0;j--) if (obs->data[j].rcv==1) break;
        if (i<j) {
            if (ts.time==0) ts=obs->data[i].time;
            if (te.time==0) te=obs->data[j].time;
            settspan(ts,te);
        }
    }
    return 1;
}
/* free obs and nav data -----------------------------------------------------*/
static void freeobsnav(obs_t *obs, nav_t *nav)
{
    trace(3,"freeobsnav:\n");

    free(obs->data); obs->data=NULL; obs->n =obs->nmax =0;
    free(nav->eph ); nav->eph =NULL; nav->n =nav->nmax =0;
    free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;
    free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
}


int read_binary_obs_qx(FILE *fobs, NavObservable *o, int *chan)
{
    int64_t tmp;
    //int32_t tmp32;
    uint8_t frame_id;
    uint8_t ichan;
    char iprn;

    o->valid = 0;

    if (fread(&frame_id, 1, 1, fobs) != 1)
        return 0;

    if (frame_id != 0x10 && frame_id != 0x11)
    {
        printf("*** not implemented frame ID=%d\n", frame_id);
        return 0;
    }

    int isDualSec = frame_id == 0x11;
    
    if (fread(&ichan, 1, 1, fobs) != 1)
        return 0;

    if (chan)
        *chan = ichan;

    if (fread(&tmp, sizeof(tmp), 1, fobs) != 1)
        return 0;

#if OLD_OBSO_Q30_TOW
#warning OLD_OBSO_Q30_TOW
        o->rx_tow = tmp / Q_30_dbl;
#else
        o->rx_tow = tmp / Q_32_dbl;
#endif

    o->tow = 0;
    o->carrier_Doppler_hz = 0;

    if (fread(&tmp, sizeof(tmp), 1, fobs) != 1)
        return 0;

    //o->carrier_phase_cyc = tmp / Q_32_dbl;
    o->carrier_phase_cyc = tmp / Q_10_dbl;

    //#warning TEMP DBG !!!
    //o->carrier_phase_cyc = tmp / Q_10_dbl / 255.75 ;/// 4092.0;
    //o->carrier_phase_cyc = (int)(-tmp / Q_10_dbl) / 409200;

    int64_t phase_quant_err = (int64_t)(o->carrier_phase_cyc * Q_10_dbl) - tmp;
    if (phase_quant_err)
        printf("\n*** phase__quant_err = %ld !!! \n\n"  , phase_quant_err);


    if (fread(&tmp, sizeof(tmp), 1, fobs) != 1)
        return 0;

    //o->pseudorange_m = tmp / Q_32_dbl;
    o->pseudorange_m = tmp / Q_10_dbl;

    if (fread(&iprn, 1, 1, fobs) != 1)
        return 0;

    if (iprn < 0)
    {
        printf("*** not implemented PRN = %d\n", iprn);
        return 0;
    }

    // w PRN b.8 jest zakodowane info czy dual secondary
    o->prn = iprn + (!isDualSec ? 0 : 0x80);

    o->valid = 1;
    //o->valid = !isDualSec;

    return 1;
}

void printdiff(FILE *diffCsv, rtk_t *rtk)
{
    sol_t sol={{0}};
    double rb[3]={0};

    int i;
    int week;
    double gpst;
    double err3d, dd, err2d, drr[3] = {0};
    /*double pos[3] = {0},*/
    double refpos[3] = {0}, enu[3] = {0};

    sol = rtk->sol;
    for (i = 0; i < 3; i++)
        rb[i] = rtk->rb[i];

    /*ecef2pos(sol.rr, pos);*/
    ecef2pos(rb, refpos);

    err3d = 0;
    for (i = 0; i < 3; i++)
    {
        dd = sol.rr[i] - rb[i];
        drr[i] = dd;
        err3d += dd * dd;
    }
    err3d = sqrt(err3d);

    ecef2enu(refpos, drr, enu);
    err2d = sqrt(enu[0]*enu[0] + enu[1]*enu[1]);

    gpst = time2gpst(rtk->sol.time, &week);

    /*fprintf(diffCsv, "%.12g; %g; %g; %g; %g; %g\n", gpst, err3d, error_LLH_m, d_ls_pvt->get_gdop(), max_abs_resp, max_abs_resc);*/
    printf("%.12g; %g; %g\n", gpst, err3d, err2d);
    fprintf(diffCsv, "%.12g; %g; %g\n", gpst, err3d, err2d);

}


// [[Rcpp
int main(int argc, char **argv)
{

    if (SAVE_TRACE)
    {
        traceopen("qx2rtkp.trace");
        tracelevel(5);
    }

    /*
        brdc1490.21n observables.qx 8
    */

    if(argc < 4)
    {
        printf("args: brdc-path qx-path chan-num\n");
        return 1;
    }
    char *brdcpath = argv[1];
    char *qxpath = argv[2];

    printf("READING BRDC :\n  %s\n", brdcpath);
    printf("READING .qx  :\n  %s\n", qxpath);

    int chan_num = atoi(argv[3]);

    printf("chan_num = %d\n", chan_num);

    if (chan_num < 0 || chan_num > MAX_CHAN)
    {
        printf("chan_num: out of range\n");
        return 1;
    }

    int isGalileo = 0;
    if (argc >= 5)
        isGalileo = atoi(argv[4]) != 0;
    printf("isGalileo = %d\n", isGalileo);

    if (isGalileo)
        printf("GAL_AS_GPS_L1L2 = %d\n", GAL_AS_GPS_L1L2);

    FILE *fqxobs = fopen(qxpath, "rb");
    if (fqxobs == NULL)
    {
        printf(".qx FILE OPEN FAILED \n");
        return 1;
    }

    FILE *diffCsv;
    char diffPath[2000];
    sprintf(diffPath, "%s%s", qxpath, ".diff.csv");
    diffCsv = fopen(diffPath, "w");
    fprintf(diffCsv, "time; diff_3D [m]; diff_2D [m]\n");

    gtime_t ts = {0}, te = {0};
    double ti = 0.0;
    int index[100] = {0,1,2,3,4,5,6,7,8};

    prcopt_t opt = prcopt_qx;
    prcopt_t *popt = &opt;

    if (isGalileo)
        popt->navsys = SYS_GAL;

    char *infile[1] = { brdcpath };

    printf("reading files ...");
    if (!readobsnav(ts,te,ti,infile,index,1,popt,&obss,&navs,stas))
    {
        printf("read errors \n");
        return 0;
    }

    rtk_t rtk;
    rtkinit(&rtk, popt);

    NavObservable obs_row[MAX_CHAN];

    for (int i = 0; i < MAX_CHAN; ++i)
        obs_row[i].valid = 0;

    NavObservable obs;

    int ch;
    int num_obs = 0;

    double prev_rx_tow =0;

    while (read_binary_obs_qx(fqxobs, &obs, &ch))
    {

        if (ch >= chan_num)
            printf("WARNING chan >= chan_num, chan=%d, chan_num=%d \n", ch, chan_num);

        int tm_chng = fabs(obs.rx_tow - prev_rx_tow) > 1e-12;

        if (tm_chng || obs_row[ch].valid)
        {
            if (!tm_chng)
                printf("WARNGING: no time change && new obs for the same PRN, PRN=%d, rxTOW=%.20g\n", obs.prn, obs.rx_tow);

            if (0)
            {
                printf("%06d | ", (int)prev_rx_tow);
                for (int i = 0; i < chan_num; ++i)
                {
                    int prn = obs_row[i].prn;
                    int sec = 0;
                    if (prn & 0x80)
                    {
                        sec = 1;
                        prn -= 0x80;
                    }

                    if (obs_row[i].valid)
                        printf("%02d%s ", prn, sec ? "*" : " ");
                    else
                        printf("--  ");
                }
                printf("\n");
            }

            obsd_t obs[MAXOBS * 2];
            int nobs = 0;

            for (int i = 0; i < chan_num; ++i)
            {
                if (!obs_row[i].valid)
                    continue;

                int sat = obs_row[i].prn;
                int sec = 0;
                if (sat & 0x80)
                {
                    sec = 1;
                    sat -= 0x80;
                }
                if (isGalileo)
                    sat += NSATGPS + NSATGLO;

                if (0 && sec)
                    continue;

                int band = sec ? 2 : 0;
                if (GAL_AS_GPS_L1L2)
                    band = sec;

                int xobs = nobs;

                // łaczenie obserwant tej samej PRN
                for (int ii = 0; ii < nobs; ++ii)
                {
                    if (obs[ii].sat != sat)
                        continue;
                    xobs = ii;
                    break;
                }

                // tj. w rtklib_conversions.cc
                obs[xobs].D[band] = obs_row[i].carrier_Doppler_hz;
                obs[xobs].P[band] = obs_row[i].pseudorange_m;
                obs[xobs].L[band] = obs_row[i].carrier_phase_cyc;

                switch (band)
                {
                    case 0:
                        obs[xobs].code[band] = (unsigned char) CODE_L1C;
                        break;
                    case 1:
                        obs[xobs].code[band] = (unsigned char) CODE_L2S;
                        break;
                    case 2:
                        obs[xobs].code[band] = (unsigned char) CODE_L5X;
                        break;
                }

                double CN0_dB_Hz_est = 55; // todo
                if (CN0_dB_Hz_est > 63.75)
                        CN0_dB_Hz_est = 63.75;
                if (CN0_dB_Hz_est < 0.0)
                        CN0_dB_Hz_est = 0.0;
                int CN0_dB_Hz = .5 + CN0_dB_Hz_est / 0.25;
                obs[xobs].SNR[band] = CN0_dB_Hz;

                // todo LLI
                if (CN0_dB_Hz_est <= -0.001)
                    obs[xobs].LLI[band] = 1;

                obs[xobs].sat = sat;

                //todo week!
                int week = navs.eph[0].week;
                obs[xobs].time = gpst2time(week, obs_row[i].rx_tow);
                if (obs_row[i].rx_tow < 18.0)
                    obs[xobs].time = timeadd(obs[xobs].time, 604800);

                obs[xobs].rcv = 1;
                if (xobs == nobs)
                    nobs++;
            }

            if (nobs > 0)
            {
                if (rtkpos(&rtk, obs, nobs, &navs))
                {
                    printdiff(diffCsv, &rtk);
                }
            }

            for (int i = 0; i < chan_num; ++i)
                obs_row[i].valid = 0;
        }
        prev_rx_tow = obs.rx_tow;
        obs_row[ch] = obs;
        num_obs++;

    }

    freeobsnav(&obss,&navs);
    rtkfree(&rtk);

    printf("EOF \n");

    printf("num_obs = %d\n", num_obs);

    fclose(diffCsv);
    diffCsv = 0;

    fclose(fqxobs);
    fqxobs = 0;

    printf("DONE\n");

    if (SAVE_TRACE)
        traceclose();

    return 0;

}

