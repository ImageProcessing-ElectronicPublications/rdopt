
/*
** Copyright 1995,1996 by Viresh Ratnakar, Miron Livny
**
** Permission to use and copy this software and its documentation
** for any non-commercial purpose and without fee is hereby granted,
** provided that the above copyright notice appear in all copies and that
** both that copyright notice and this permission notice appear in
** supporting documentation.
**
**
** The University of Wisconsin and the copyright holders make no
** representations about the suitability of this software for any
** purpose.  It is provided "as is" without express or implied warranty.
**
**
** THE UNIVERSITY OF WISCONSIN AND THE COPYRIGHT HOLDERS DISCLAIM ALL
** WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES
** OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE UNIVERSITY OF
** WISCONSIN OR THE COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
** OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
** OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
** OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
** OR PERFORMANCE OF THIS SOFTWARE.
**
** Author:  Viresh Ratnakar
** Email:   ratnakar@cs.wisc.edu
**
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <limits.h>
#include <math.h>
#include <sys/wait.h>
#include <unistd.h>
#include <signal.h>
#include <sys/time.h>


#include "rdopt.h"

/***
**** WAIT proc used by parents to wait for child to terminate
***/


static int ProcWait(char *proc, pid_t pid)
{
    int procwpid, procstatus;

    if ( ((procwpid = waitpid(pid, &procstatus, 0)) == -1) ||
            (!WIFEXITED(procstatus)) ||
            (WEXITSTATUS(procstatus) != 0) ) return(0);
    return(1);
}

extern boolean GetActualBpp(OptimJob *job, char *qfile, FFLOAT *ans)
{
    int piped[2], pid;
    int num_bytes;
    char *args[30], args_store[30][100];
    int last_arg, i;
    int c;

#ifdef DONT_HAVE_CJPEG
    return(0);
#endif

    if (job->TheImage.ImKind == IM_RAW) return(0);
    else if (!strcmp(job->TheImage.ImFileName,"")) return(0);

    if (pipe(piped) == -1) return(0);

    if ((pid=fork()) == -1)
    {
        fprintf(stdout,"Too many processes: could not fork\n");
        return(0);
    }

    if (pid==0)
    {
        /* child */
        close(piped[0]);
        if (dup2(piped[1],1) == -1)
        {
            fprintf(stdout,"Stdvqe could not dup2\n");
            exit(1);
        }

        last_arg = 0;

        strcpy(args_store[last_arg],"cjpeg");
        last_arg++;

        strcpy(args_store[last_arg],"-optimize");
        last_arg++;

        strcpy(args_store[last_arg],"-dct");
        last_arg++;
        strcpy(args_store[last_arg],"float");
        last_arg++;

        strcpy(args_store[last_arg],"-qtables");
        last_arg++;
        strcpy(args_store[last_arg],qfile);
        last_arg++;

        if (job->ThreshSpan > 0)
        {
            strcpy(args_store[last_arg],"-thresh");
            last_arg++;
            strcpy(args_store[last_arg],qfile);
            last_arg++;
        }

        if (job->TheImage.NumComponents == 1)
        {
            strcpy(args_store[last_arg],"-grayscale");
            last_arg++;
        }

        if (job->NumTables == 3)
        {
            strcpy(args_store[last_arg],"-qslots");
            last_arg++;
            strcpy(args_store[last_arg],"0,1,2");
            last_arg++;
        }

        strcpy(args_store[last_arg],job->TheImage.ImFileName);
        last_arg++;

        for (i=0; i<last_arg; i++) args[i] = args_store[i];
        args[last_arg] = (char *) 0;
        execvp("cjpeg",args);
    }
    else
    {
        /* parent */
        close(piped[1]);
        /* read bytes from pipe */
        num_bytes = 0;
        while (read(piped[0],(char *) &c, 1) == 1) num_bytes++;
        close(piped[0]);

        *ans = ((FFLOAT) (num_bytes*8))/job->NumPixelsForBpp;

        if (!ProcWait("cjpeg",pid)) return(0);

        return(1);
    }

}

#define BUFF_SIZE 1024

extern boolean Compress(OptimJob *job, char *qfile, char *cfile,
                        int *nbytes)
{
    int piped[2], pid;
    int num_bytes;
    char *args[30], args_store[30][100];
    int last_arg, i, fd;
    char cbuff[BUFF_SIZE];

#ifdef DONT_HAVE_CJPEG
    return(0);
#endif

    if (job->TheImage.ImKind == IM_RAW) return(0);
    else if (!strcmp(job->TheImage.ImFileName,"")) return(0);

    if (pipe(piped) == -1) return(0);

    if ((pid=fork()) == -1)
    {
        fprintf(stdout,"Too many processes: could not fork\n");
        return(0);
    }

    if (pid==0)
    {
        /* child */
        close(piped[0]);
        if (dup2(piped[1],1) == -1)
        {
            fprintf(stdout,"Stdvqe could not dup2\n");
            exit(1);
        }

        last_arg = 0;

        strcpy(args_store[last_arg],"cjpeg");
        last_arg++;

        strcpy(args_store[last_arg],"-optimize");
        last_arg++;

        strcpy(args_store[last_arg],"-dct");
        last_arg++;
        strcpy(args_store[last_arg],"float");
        last_arg++;

        strcpy(args_store[last_arg],"-qtables");
        last_arg++;
        strcpy(args_store[last_arg],qfile);
        last_arg++;

        if (job->ThreshSpan > 0)
        {
            strcpy(args_store[last_arg],"-thresh");
            last_arg++;
            strcpy(args_store[last_arg],qfile);
            last_arg++;
        }

        if (job->TheImage.NumComponents == 1)
        {
            strcpy(args_store[last_arg],"-grayscale");
            last_arg++;
        }

        if (job->NumTables == 3)
        {
            strcpy(args_store[last_arg],"-qslots");
            last_arg++;
            strcpy(args_store[last_arg],"0,1,2");
            last_arg++;
        }

        strcpy(args_store[last_arg],job->TheImage.ImFileName);
        last_arg++;

        for (i=0; i<last_arg; i++) args[i] = args_store[i];
        args[last_arg] = (char *) 0;
        execvp("cjpeg",args);
    }
    else
    {
        /* parent */
        close(piped[1]);

        /* open file for writing */
        if ((fd = open(cfile, O_WRONLY | O_TRUNC | O_CREAT, S_IRUSR | S_IWUSR)) < 0)
        {
            ProcWait("cjpeg",pid);
            return(FALSE);
        }

        /* read bytes from pipe */
        num_bytes = 0;
        while ((i = read(piped[0],cbuff, BUFF_SIZE)) > 0)
        {
            num_bytes += i;
            write(fd,cbuff,i);
        }

        close(fd);
        close(piped[0]);

        *nbytes = num_bytes;

        if (!ProcWait("cjpeg",pid)) return(0);

        return(TRUE);
    }

}

