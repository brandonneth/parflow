/*BHEADER**********************************************************************
 * (c) 1995   The Regents of the University of California
 *
 * See the file COPYRIGHT_and_DISCLAIMER for a complete copyright
 * notice, contact person, and disclaimer.
 *
 * $Revision: 1.1.1.1 $
 *********************************************************************EHEADER*/
#include "amps.h"

char *amps_recv(src)
int src;
{
   void *ret;
   int check_again;
   int not_found;
   AMPS_ShMemBuffer *buf;
   
   not_found = 1;

   /* first look on the local queue */
   buf = amps_shmem_info -> buffers_local[amps_rank];
   while(buf)
   {
      if( *(int *)(buf -> data) == src)
      {
	 /* remove the found buffer from the list */
	 if(buf -> prev)
	    buf -> prev -> next = buf -> next;
	 else
	    /* start of list */
	    amps_shmem_info -> buffers_local[amps_rank]
	       = buf -> next;
	 
	 if (buf -> next )
	    buf -> next -> prev = buf -> prev;
	 else
	    amps_shmem_info -> buffers_local_end[amps_rank] = 
	       buf -> prev;

	 /* return the buffer that was found */
	 not_found = 0;
	 break;
      }
      buf = buf -> next;
   }
   

   check_again = 1; 
   while(not_found & check_again)
   {
      /* try and get semaphore */
      if(uscpsema(amps_shmem_info -> sema[amps_rank]))
      {
	 /* there was a waiting message */
	 ussetlock(amps_shmem_info -> locks[amps_rank]);
      
	 buf = amps_shmem_info -> buffers[amps_rank];
	 
	 /* remove the first buffer from the list */
	 if(buf -> prev)
	    buf -> prev -> next = buf -> next;
	 else
	    /* start of list */
	    amps_shmem_info -> buffers[amps_rank] = buf -> next;
	 
	 if (buf -> next )
	    buf -> next -> prev = buf -> prev;
	 else
	    amps_shmem_info -> buffers_end[amps_rank] = buf -> prev;
	    
	 usunsetlock(amps_shmem_info -> locks[amps_rank]);

	 /* if this is the one we are looking for then we are done */
	 if( *(int *)(buf -> data) == src)
	 {
	    not_found = 0;
	 }
	 else
	 {
	    /* not the buffer we are looking for add to the local buffer
	       list */
	    buf -> prev = amps_shmem_info -> buffers_local_end[amps_rank];
	    buf -> next = NULL;
	    
	    if( amps_shmem_info -> buffers_local_end[amps_rank] )
	       amps_shmem_info -> buffers_local_end[amps_rank] -> next = buf;
	    
	    amps_shmem_info -> buffers_local_end[amps_rank] = buf;
	    
	    if ( !amps_shmem_info -> buffers_local[amps_rank] )
	       amps_shmem_info -> buffers_local[amps_rank] = buf;
	    
	 }
      }
      else
	 check_again = 0;
   }

   if(!not_found)
   {
      ret = (buf -> data) + sizeof(double);
      usfree(buf, amps_arena);
   }
   else
      ret = NULL;

   return ret;
}



amps_Handle amps_IRecv(comm, source, invoice)
amps_Comm comm;
int source;
amps_Invoice invoice;
{
   char *buffer;

   if((buffer = amps_recv(source)))
   {
      /* we have recvd it */
      amps_ClearInvoice(invoice);
      amps_unpack(comm, invoice, buffer);
      AMPS_PACK_FREE_LETTER(comm, invoice, buffer);
      return NULL;
   }
   else
      /* did not recv it */
      return amps_NewHandle(comm, source, invoice, NULL);
}