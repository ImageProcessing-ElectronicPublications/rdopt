/*
** Copyright 1995 by Viresh Ratnakar, Miron Livny
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

#include "precision.h"
#include "qentry.h"


extern int MapQentry(int n, int q)
{
    /*** dependence on n later **/
#if (QTABBITS==8)
    return(q);
#elif (QTABBITS==12)
    if (q <= 16) return(q);
    else if (q <= 32) return(2*q);
    else if (q <= 64) return(4*q);
    else if (q <= 128) return(8*q);
    else return(16*q);
#elif (QTABBITS==16)
    if (q <= 4) return(q);
    else if (q <= 16) return(2*q);
    else if (q <= 32) return(4*q);
    else if (q <= 64) return(16*q);
    else if (q <= 128) return(64*q);
    else return(256*q);
#endif
}

extern int UnMapQentry(int n, int q)
{

#if (QTABBITS==8)
    return(q);
#elif (QTABBITS==12)
    if (q <= 16) return(q);
    else if (q <= 32*2)
    {
        if (q < 17*2) return(17);
        else return(q/2);
    }
    else if (q <= 64*4)
    {
        if (q < 33*4) return(33);
        else return(q/4);
    }
    else if (q <= 128*8)
    {
        if (q < 65*8) return(65);
        else return(q/8);
    }
    else
    {
        if (q < 129*16) return(129);
        else return(q/16);
    }
#elif (QTABBITS==16)
    if (q <= 4) return(q);
    else if (q <= 16*2)
    {
        if (q < 5*2) return(5);
        else return(q/2);
    }
    else if (q <= 32*4)
    {
        if (q < 17*4) return(17);
        else return(q/4);
    }
    else if (q <= 64*16)
    {
        if (q < 33*16) return(33);
        else return(q/16);
    }
    else if (q <= 128*64)
    {
        if (q < 65*64) return(65);
        else return(q/64);
    }
    else
    {
        if (q < 129*256) return(129);
        else return(q/256);
    }
#endif
}
