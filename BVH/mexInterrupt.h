#ifndef MEX_INTERRUPT_H
#define MEX_INTERRUPT_H
/* ------------------------------------------------------------------------- *
 *  mexInterrupt.h  --  Ctrl-C support for long-running MEX query loops.
 *
 *  MATLAB services Ctrl-C only between interpreter statements, so a MEX stuck
 *  in a tight compute loop cannot normally be cancelled.  libut exposes the
 *  (undocumented but stable across releases) entry point utIsInterruptPending(),
 *  which reports whether a Ctrl-C is queued.  Link the MEX with  -lut .
 *
 *  IMPORTANT -- we only POLL, we never CONSUME the interrupt: utIsInterruptPending()
 *  leaves the Ctrl-C queued in MATLAB.  On detection we simply stop the loop and
 *  let mexFunction return NORMALLY.  MATLAB then services the still-pending Ctrl-C
 *  and aborts the calling .m wrapper with the usual "Operation terminated by user".
 *  We deliberately do NOT call mexErrMsgIdAndTxt: in a C++ MEX that longjmps and
 *  skips C++ destructors (node pools, std::vector, ...) -> leaks.  Returning
 *  normally runs every destructor, so this path is leak-free.
 *
 *  Threading -- utIsInterruptPending() touches MATLAB internals and must be
 *  called only from the MATLAB main thread.  Inside an OpenMP parallel region the
 *  master (omp_get_thread_num()==0) IS the thread that entered the region, i.e.
 *  the MEX main thread.  So ONLY the master polls libut; every thread reads the
 *  shared monotonic (false->true) abort flag and bails.  The flag is a per-TU
 *  static (each MEX is one translation unit), so it must be reset with
 *  mexClearInterrupt() at gateway entry (a MEX stays loaded across calls).
 *
 *  Usage:
 *    at gateway entry:            mexClearInterrupt();
 *    top of each loop iteration:  if( mexInterrupted( i ) ) break;   // i = loop counter
 * ------------------------------------------------------------------------- */

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __cplusplus
extern "C" bool utIsInterruptPending( void );
#else
extern bool utIsInterruptPending( void );
#endif

/* shared monotonic abort flag -- one instance per translation unit (per MEX) */
static volatile bool g_mexAbort = false;

static inline void mexClearInterrupt( void ) { g_mexAbort = false; }
static inline bool mexWasInterrupted( void ) { return g_mexAbort;  }

/* poll granularity: consult libut once every 1<<MEX_INT_LOG2 iterations
 * (only on the master thread), so the overhead per point is a single flag read */
#ifndef MEX_INT_LOG2
#define MEX_INT_LOG2 12                       /* every 4096 iterations */
#endif

/* returns true once a Ctrl-C has been seen; the caller should then stop the loop */
static inline bool mexInterrupted( long long i )
{
  /* per-iteration cost is exactly this: a flag read + a masked compare.       */
  if( g_mexAbort ) return true;                      /* every thread sees it   */
  if( ( i & ( ( 1LL << MEX_INT_LOG2 ) - 1 ) ) != 0 ) return false;  /* throttle */
  /* the rest runs only once per 1<<MEX_INT_LOG2 iterations:                    */
#ifdef _OPENMP
  if( omp_get_thread_num() != 0 ) return false;      /* only the master polls   */
#endif
  if( utIsInterruptPending() ) { g_mexAbort = true; return true; }
  return false;
}

#endif /* MEX_INTERRUPT_H */
