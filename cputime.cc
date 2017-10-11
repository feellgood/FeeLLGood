#include <time.h>
#include <sys/times.h>
#include <unistd.h>

double cputime(void) {
   	static struct tms tmsbuf;
	static double ticks = -1.0;
	         
	if (ticks < 0.0) {
	   ticks = (double)sysconf(_SC_CLK_TCK);
			 }
	   times(&tmsbuf);
	   return tmsbuf.tms_utime/ticks;
}
